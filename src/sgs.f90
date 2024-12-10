! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
#define _CHANNEL
module mod_sgs
  use mpi
  use mod_precision, only: rp,sp,dp,i8,MPI_REAL_RP
  use mod_common_mpi, only: ierr
  use mod_params, only: c_smag,big
  use mod_typedef, only: bound
  use mod_bound, only: boundp,bounduvw
  implicit none
  private
  public cmpt_sgs
  contains
  !
  subroutine cmpt_sgs(sgstype,n,ng,lo,hi,cbcvel,cbcsgs,bcs,nb,is_bound,lwm,l,dl,dli,zc,zf,dzc,dzf, &
                      dzci,dzfi,visc,hwm,u,v,w,bcuf,bcvf,bcwf,bcu_mag,bcv_mag,bcw_mag,visct)
    !
    ! compute subgrid viscosity at cell centers
    ! the LES with the dynamic model is ~2 times the cost of the LES with the static one.
    ! Acceleration can be further achieved by calling only one boundp to do the data
    ! exchange of multiple variables. Note that it does not save time to simply mask
    ! wall bc's in boundp. The dynamic version yields quite good results for Re_tau=
    ! 395,550,1000, with <=5% errors in the friction coefficient. Clipping is necessary
    ! to avoid negative values of averaged eddy viscosity (duct flow).
    !
    ! 2D filter is used in the first off-wall layer (Bae, Orlandi, LESGO and Balaras, 1995).
    ! It is difficult to do 3D filtering of Sij in the first layer, though feasible for
    ! velocity.
    !
    implicit none
    character(len=*), intent(in) :: sgstype
    integer , intent(in ), dimension(3) :: n,ng,lo,hi
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcsgs
    type(Bound), intent(in   ) :: bcs
    type(Bound), intent(inout) :: bcuf,bcvf,bcwf
    type(Bound), intent(in   ) :: bcu_mag,bcv_mag,bcw_mag
    integer , intent(in ), dimension(0:1,3)      :: nb,lwm
    logical , intent(in ), dimension(0:1,3)      :: is_bound
    real(rp), intent(in ), dimension(3)          :: l,dl,dli
    real(rp), intent(in ), dimension(0:)         :: zc,zf,dzc,dzf,dzci,dzfi
    real(rp), intent(in )                        :: visc,hwm
    real(rp), intent(in ), dimension(0:,0:,0:)   :: u,v,w
    real(rp), intent(out), dimension(0:,0:,0:)   :: visct
    real(rp), allocatable, dimension(:,:,:)  , save :: s0,uc,vc,wc,uf,vf,wf,alph2
    real(rp), target,  allocatable, dimension(:,:,:,:), save :: wk,sij,mij
    real(rp), pointer, contiguous , dimension(:,:,:,:), save :: lij
    real(rp) :: dxi,dyi,visci,tauw(2),dw(6),dw_plus, &
                fd,del,dw_min,tauw_s,mij_s(6),lij_s(6)
    real(rp), save :: is_wall(6)
    logical, save :: is_first = .true.
    integer :: i,j,k,m,ind_min
    real(rp), parameter :: one_third = 1._rp/3._rp
    !
    select case(trim(sgstype))
    case('none')
      if(is_first) then
        is_first = .false.
        !$acc kernels default(present) async(1)
        visct(:,:,:) = 0._rp
        !$acc end kernels
      end if
    case('smag')
      if(is_first) then
        is_first = .false.
        allocate(s0(0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk(0:n(1)+1,0:n(2)+1,0:n(3)+1,3))
        !$acc enter data create(s0,wk) async(1)
        is_wall = 0._rp
        if(is_bound(0,1).and.cbcvel(0,1,1)=='D') is_wall(1) = 1._rp
        if(is_bound(1,1).and.cbcvel(1,1,1)=='D') is_wall(2) = 1._rp
        if(is_bound(0,2).and.cbcvel(0,2,2)=='D') is_wall(3) = 1._rp
        if(is_bound(1,2).and.cbcvel(1,2,2)=='D') is_wall(4) = 1._rp
        if(is_bound(0,3).and.cbcvel(0,3,3)=='D') is_wall(5) = 1._rp
        if(is_bound(1,3).and.cbcvel(1,3,3)=='D') is_wall(6) = 1._rp
        !$acc enter data copyin(is_wall) async(1)
      end if
      !$acc kernels default(present) async(1)
      wk(:,:,:,1) = u(:,:,:)
      wk(:,:,:,2) = v(:,:,:)
      wk(:,:,:,3) = w(:,:,:)
      !$acc end kernels
      call extrapolate(n,is_bound,dzci,wk(:,:,:,1),iface=1,lwm=lwm)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,2),iface=2,lwm=lwm)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,3),iface=3,lwm=lwm)
      call strain_rate(n,dli,dzci,dzfi,wk(:,:,:,1),wk(:,:,:,2),wk(:,:,:,3),s0)
      !
      dxi = dli(1)
      dyi = dli(2)
      visci = 1._rp/visc
      !
      !$acc parallel loop collapse(3) default(present) async(1) &
      !$acc private(dw,tauw,ind_min,fd,del,tauw_s,dw_min,dw_plus)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            if(sum(is_wall(1:6))==0) then
              ! triperiodic
              fd = 1._rp
            else
              ! van Driest damping function
              dw(1) = dl(1)*(i-0.5)
              dw(2) = dl(1)*(n(1)-i+0.5)
              dw(3) = dl(2)*(j-0.5)
              dw(4) = dl(2)*(n(2)-j+0.5)
              dw(5) = zc(k)
              dw(6) = l(3)-zc(k)
              dw = dw*is_wall + big*(1._rp-is_wall)
              ind_min = minloc(dw,1)
              dw_min = dw(ind_min)
              !
              if(ind_min==1) then
                tauw(1) = v(1,j,k)-v(0,j,k)+v(1,j-1,k)-v(0,j-1,k)
                tauw(2) = w(1,j,k)-w(0,j,k)+w(1,j,k-1)-w(0,j,k-1)
                tauw_s = sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dxi
              else if(ind_min==2) then
                tauw(1) = v(n(1),j,k)-v(n(1)+1,j,k)+v(n(1),j-1,k)-v(n(1)+1,j-1,k)
                tauw(2) = w(n(1),j,k)-w(n(1)+1,j,k)+w(n(1),j,k-1)-w(n(1)+1,j,k-1)
                tauw_s = sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dxi
              else if(ind_min==3) then
                tauw(1) = u(i,1,k)-u(i,0,k)+u(i-1,1,k)-u(i-1,0,k)
                tauw(2) = w(i,1,k)-w(i,0,k)+w(i,1,k-1)-w(i,0,k-1)
                tauw_s  = sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dyi
              else if(ind_min==4) then
                tauw(1) = u(i,n(2),k)-u(i,n(2)+1,k)+u(i-1,n(2),k)-u(i-1,n(2)+1,k)
                tauw(2) = w(i,n(2),k)-w(i,n(2)+1,k)+w(i,n(2),k-1)-w(i,n(2)+1,k-1)
                tauw_s  = sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dyi
              else if(ind_min==5) then
                tauw(1) = u(i,j,1)-u(i,j,0)+u(i-1,j,1)-u(i-1,j,0)
                tauw(2) = v(i,j,1)-v(i,j,0)+v(i,j-1,1)-v(i,j-1,0)
                tauw_s  = sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dzci(0)
              else if(ind_min==6) then
                tauw(1) = u(i,j,n(3))-u(i,j,n(3)+1)+u(i-1,j,n(3))-u(i-1,j,n(3)+1)
                tauw(2) = v(i,j,n(3))-v(i,j,n(3)+1)+v(i,j-1,n(3))-v(i,j-1,n(3)+1)
                tauw_s  = sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dzci(n(3))
              end if
              tauw_s = 0.5_rp*visc*tauw_s
              dw_plus = dw_min*sqrt(tauw_s)*visci
              fd = 1._rp-exp(-dw_plus/25._rp)
            end if
            !
            del = (dl(1)*dl(2)*dzf(k))**(one_third)
            visct(i,j,k) = (c_smag*del*fd)**2*s0(i,j,k)
          end do
        end do
      end do
    case('dsmag')
      if(is_first) then
        is_first = .false.
        allocate(uc (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 vc (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wc (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 uf (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 vf (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wf (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 s0 (0:n(1)+1,0:n(2)+1,0:n(3)+1  ), &
                 wk (0:n(1)+1,0:n(2)+1,0:n(3)+1,6), &
                 sij(0:n(1)+1,0:n(2)+1,0:n(3)+1,6), &
                 mij(0:n(1)+1,0:n(2)+1,0:n(3)+1,6))
        allocate(alph2(0:n(1)+1,0:n(2)+1,0:n(3)+1))
        !$acc enter data create(uc,vc,wc,uf,vf,wf) async(1)
        !$acc enter data create(s0,wk,sij,mij) async(1)
        !$acc enter data create(alph2) async(1)
        call cmpt_alph2(n,is_bound,cbcvel,alph2)
      end if
      !
      !$acc kernels default(present) async(1)
      wk(:,:,:,1) = u(:,:,:)
      wk(:,:,:,2) = v(:,:,:)
      wk(:,:,:,3) = w(:,:,:)
      !$acc end kernels
      call extrapolate(n,is_bound,dzci,wk(:,:,:,1),iface=1,lwm=lwm)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,2),iface=2,lwm=lwm)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,3),iface=3,lwm=lwm)
      call strain_rate(n,dli,dzci,dzfi,wk(:,:,:,1),wk(:,:,:,2),wk(:,:,:,3),s0,sij)
      !
      !$acc kernels default(present) async(1)
      visct(:,:,:) = s0(:,:,:)
      !$acc end kernels
      !
      ! Mij
      !
      ! periodic/patched bc's are updated, wall bc ghost points are set
      ! by extrapolation from the interior.
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,s0)
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,sij(:,:,:,1))
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,sij(:,:,:,2))
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,sij(:,:,:,3))
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,sij(:,:,:,4))
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,sij(:,:,:,5))
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,sij(:,:,:,6))
      !$acc parallel loop collapse(3) default(present) async(1)
      do k = 0,n(3)+1
        do j = 0,n(2)+1
          do i = 0,n(1)+1
            wk(i,j,k,1) = s0(i,j,k)*sij(i,j,k,1)
            wk(i,j,k,2) = s0(i,j,k)*sij(i,j,k,2)
            wk(i,j,k,3) = s0(i,j,k)*sij(i,j,k,3)
            wk(i,j,k,4) = s0(i,j,k)*sij(i,j,k,4)
            wk(i,j,k,5) = s0(i,j,k)*sij(i,j,k,5)
            wk(i,j,k,6) = s0(i,j,k)*sij(i,j,k,6)
          end do
        end do
      end do
#if !defined(_FILTER_2D)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,1),iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,2),iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,3),iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,4),iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,5),iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,6),iface=0,cbc=cbcvel)
      call filter3d(n,wk(:,:,:,1),mij(:,:,:,1))
      call filter3d(n,wk(:,:,:,2),mij(:,:,:,2))
      call filter3d(n,wk(:,:,:,3),mij(:,:,:,3))
      call filter3d(n,wk(:,:,:,4),mij(:,:,:,4))
      call filter3d(n,wk(:,:,:,5),mij(:,:,:,5))
      call filter3d(n,wk(:,:,:,6),mij(:,:,:,6))
      !
      !$acc kernels default(present) async(1)
      wk(:,:,:,1) = u(:,:,:)
      wk(:,:,:,2) = v(:,:,:)
      wk(:,:,:,3) = w(:,:,:)
      !$acc end kernels
      call extrapolate(n,is_bound,dzci,wk(:,:,:,1),iface=1,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,2),iface=2,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,3),iface=3,cbc=cbcvel)
      call filter3d(n,wk(:,:,:,1),uf)
      call filter3d(n,wk(:,:,:,2),vf)
      call filter3d(n,wk(:,:,:,3),wf)
#else
      call filter2d(n,wk(:,:,:,1),mij(:,:,:,1))
      call filter2d(n,wk(:,:,:,2),mij(:,:,:,2))
      call filter2d(n,wk(:,:,:,3),mij(:,:,:,3))
      call filter2d(n,wk(:,:,:,4),mij(:,:,:,4))
      call filter2d(n,wk(:,:,:,5),mij(:,:,:,5))
      call filter2d(n,wk(:,:,:,6),mij(:,:,:,6))
      !
      call filter2d(n,u,uf)
      call filter2d(n,v,vf)
      call filter2d(n,w,wf)
#endif
      ! when not using a wall model, all bc's are used for computing strain rate, and
      ! bcv/v/wf are equal to bcu/v/w=0 for no-slip walls. When using a wall model,
      ! wall bc's for wall-parallel velocity are not used due to extrapolation, but
      ! the wall bc's for the wall-normal velicity are used to impose no-penetration
      ! on the filtered velocity. Also, periodic/patched bc's are used. Hence,
      ! there is no need to update wall model bc's in bounduvw. Filtered velocity
      ! should satisfy the wall bc's, evidenced by the fact that each mode satisfies
      ! the no-slip/no-penetration bc's if a spectral filter is applied.
      call bounduvw(cbcvel,n,bcuf,bcvf,bcwf,bcu_mag,bcv_mag,bcw_mag,nb,is_bound,lwm, &
                    l,dl,zc,zf,dzc,dzf,visc,hwm,.false.,.false.,uf,vf,wf)
      call extrapolate(n,is_bound,dzci,uf,iface=1,lwm=lwm)
      call extrapolate(n,is_bound,dzci,vf,iface=2,lwm=lwm)
      call extrapolate(n,is_bound,dzci,wf,iface=3,lwm=lwm)
      call strain_rate(n,dli,dzci,dzfi,uf,vf,wf,s0,sij)
      !$acc parallel loop collapse(3) default(present) async(1)
      do k = 1,n(3)
        do j = 1,n(2)
          do i = 1,n(1)
            !$acc loop seq
            do m = 1,6
              mij(i,j,k,m) = 2._rp*(mij(i,j,k,m)-alph2(i,j,k)*s0(i,j,k)*sij(i,j,k,m))
            end do
          end do
        end do
      end do
      !
      ! Lij, stored in sij
      !
      lij => sij
      call interpolate(n,u,v,w,uc,vc,wc)
      ! periodic/patched bc's are updated, wall bc ghost points are set
      ! by extrapolation from the interior.
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,uc)
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,vc)
      call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,wc)
      !$acc parallel loop collapse(3) default(present) async(1)
      do k = 0,n(3)+1
        do j = 0,n(2)+1
          do i = 0,n(1)+1
            wk(i,j,k,1) = uc(i,j,k)*uc(i,j,k)
            wk(i,j,k,2) = vc(i,j,k)*vc(i,j,k)
            wk(i,j,k,3) = wc(i,j,k)*wc(i,j,k)
            wk(i,j,k,4) = uc(i,j,k)*vc(i,j,k)
            wk(i,j,k,5) = uc(i,j,k)*wc(i,j,k)
            wk(i,j,k,6) = vc(i,j,k)*wc(i,j,k)
          end do
        end do
      end do
#if !defined(_FILTER_2D)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,1),iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,2),iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,3),iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,4),iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,5),iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wk(:,:,:,6),iface=0,cbc=cbcvel)
      call filter3d(n,wk(:,:,:,1),lij(:,:,:,1))
      call filter3d(n,wk(:,:,:,2),lij(:,:,:,2))
      call filter3d(n,wk(:,:,:,3),lij(:,:,:,3))
      call filter3d(n,wk(:,:,:,4),lij(:,:,:,4))
      call filter3d(n,wk(:,:,:,5),lij(:,:,:,5))
      call filter3d(n,wk(:,:,:,6),lij(:,:,:,6))
      !
      call extrapolate(n,is_bound,dzci,uc,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,vc,iface=0,cbc=cbcvel)
      call extrapolate(n,is_bound,dzci,wc,iface=0,cbc=cbcvel)
      call filter3d(n,uc,uf)
      call filter3d(n,vc,vf)
      call filter3d(n,wc,wf)
#else
      call filter2d(n,wk(:,:,:,1),lij(:,:,:,1))
      call filter2d(n,wk(:,:,:,2),lij(:,:,:,2))
      call filter2d(n,wk(:,:,:,3),lij(:,:,:,3))
      call filter2d(n,wk(:,:,:,4),lij(:,:,:,4))
      call filter2d(n,wk(:,:,:,5),lij(:,:,:,5))
      call filter2d(n,wk(:,:,:,6),lij(:,:,:,6))
      !
      call filter2d(n,uc,uf)
      call filter2d(n,vc,vf)
      call filter2d(n,wc,wf)
#endif
      !$acc parallel loop collapse(3) default(present) async(1) &
      !$acc private(mij_s,lij_s)
      do k = 1,n(3)
        do j = 1,n(2)
          do i = 1,n(1)
            !
            mij_s = mij(i,j,k,1:6)
            lij_s = lij(i,j,k,1:6)
            !
            lij_s(1) = lij_s(1) - uf(i,j,k)*uf(i,j,k)
            lij_s(2) = lij_s(2) - vf(i,j,k)*vf(i,j,k)
            lij_s(3) = lij_s(3) - wf(i,j,k)*wf(i,j,k)
            lij_s(4) = lij_s(4) - uf(i,j,k)*vf(i,j,k)
            lij_s(5) = lij_s(5) - uf(i,j,k)*wf(i,j,k)
            lij_s(6) = lij_s(6) - vf(i,j,k)*wf(i,j,k)
            !
            wk(i,j,k,1) = mij_s(1)*lij_s(1) + &
                          mij_s(2)*lij_s(2) + &
                          mij_s(3)*lij_s(3) + &
                         (mij_s(4)*lij_s(4) + &
                          mij_s(5)*lij_s(5) + &
                          mij_s(6)*lij_s(6))*2._rp
            wk(i,j,k,2) = mij_s(1)*mij_s(1) + &
                          mij_s(2)*mij_s(2) + &
                          mij_s(3)*mij_s(3) + &
                         (mij_s(4)*mij_s(4) + &
                          mij_s(5)*mij_s(5) + &
                          mij_s(6)*mij_s(6))*2._rp
          end do
        end do
      end do
#if defined(_DIT)
      call ave0d_dit(ng,lo,hi,l,dl,dzf,wk(:,:,:,1))
      call ave0d_dit(ng,lo,hi,l,dl,dzf,wk(:,:,:,2))
#elif defined(_CHANNEL)
      call ave1d_channel(ng,lo,hi,3,l,dl,dzf,wk(:,:,:,1))
      call ave1d_channel(ng,lo,hi,3,l,dl,dzf,wk(:,:,:,2))
#elif defined(_DUCT)
      call ave2d_duct(ng,lo,hi,1,l,dl,dzf,wk(:,:,:,1))
      call ave2d_duct(ng,lo,hi,1,l,dl,dzf,wk(:,:,:,2))
#elif defined(_CAVITY)
      !
#endif
      ! cs = (c_smag*del)^2
      !$acc parallel loop collapse(3) default(present) async(1)
      do k = 1,n(3)
        do j = 1,n(2)
          do i = 1,n(1)
            visct(i,j,k) = visct(i,j,k)*wk(i,j,k,1)/wk(i,j,k,2)
            visct(i,j,k) = max(visct(i,j,k),0._rp)
          end do
        end do
      end do
    case('amd')
      print*, 'ERROR: AMD model not yet implemented'
    case default
      print*, 'ERROR: unknown SGS model'
    end select
  end subroutine cmpt_sgs
  !
  subroutine ave0d_dit(ng,lo,hi,l,dl,dz,p)
    !
    ! average a variable over three directions
    !
    ! ng    -> global domain sizes
    ! lo,hi -> upper and lower extents of the input array
    ! dl,l  -> uniform grid spacing and length arrays
    ! dz    -> local z grid spacing array (should work also with the global one)
    ! p     -> 3D scalar field
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,lo,hi
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(lo(3)-1:) :: dz
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(dp) :: p0d
    real(dp) :: grid_area_ratio
    integer :: i,j,k
    !
    grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
    p0d = 0._rp
    !$acc data copy(p0d) async(1)
    !$acc parallel loop collapse(3) default(present) reduction(+:p0d) async(1)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          p0d = p0d + p(i,j,k)*grid_area_ratio*dz(k)/l(3)
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,p0d,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    !$acc data copyin(p0d) async(1)
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          p(i,j,k) = p0d
        end do
      end do
    end do
    !$acc end data
  end subroutine ave0d_dit
  !
  subroutine ave1d_channel(ng,lo,hi,idir,l,dl,dz,p)
    !
    ! average a variable over two directions
    ! adapted from out1d
    !
    ! ng    -> global domain sizes
    ! lo,hi -> upper and lower extents of the input array
    ! idir  -> direction of the profile
    ! dl,l  -> uniform grid spacing and length arrays
    ! dz    -> local z grid spacing array (should work also with the global one)
    ! p     -> 3D scalar field
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(lo(3)-1:) :: dz
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(dp) :: p1d(ng(idir))
    real(dp) :: grid_area_ratio,p1d_s
    integer :: i,j,k
    !
    select case(idir)
    case(3)
      grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
      !$acc data copyout(p1d) async(1)
      !$acc kernels default(present) async(1)
      p1d(:) = 0._rp
      !$acc end kernels
      !$acc parallel loop gang default(present) private(p1d_s) async(1)
      do k=lo(3),hi(3)
        p1d_s = 0._rp
        !$acc loop vector collapse(2) reduction(+:p1d_s)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            p1d_s = p1d_s + p(i,j,k)
          end do
        end do
        p1d(k) = p1d_s*grid_area_ratio
      end do
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      !$acc data copyin(p1d) async(1)
      !$acc kernels default(present) async(1)
      do k=lo(3),hi(3)
        p(:,:,k) = p1d(k)
      end do
      !$acc end kernels
      !$acc end data
    case(2)
      grid_area_ratio = dl(1)/(l(1)*l(3))
      !$acc data copyout(p1d) async(1)
      !$acc kernels default(present) async(1)
      p1d(:) = 0._rp
      !$acc end kernels
      !$acc parallel loop gang default(present) private(p1d_s) async(1)
      do j=lo(2),hi(2)
        p1d_s = 0._rp
        !$acc loop vector collapse(2) reduction(+:p1d_s)
        do k=lo(3),hi(3)
          do i=lo(1),hi(1)
            p1d_s = p1d_s + p(i,j,k)*dz(k)
          end do
        end do
        p1d(j) = p1d_s*grid_area_ratio
      end do
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(2),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      !$acc data copyin(p1d) async(1)
      !$acc kernels default(present) async(1)
      do j=lo(2),hi(2)
        p(:,j,:) = p1d(j)
      end do
      !$acc end kernels
      !$acc end data
    case(1)
      grid_area_ratio = dl(2)/(l(2)*l(3))
      !$acc data copyout(p1d) async(1)
      !$acc kernels default(present) async(1)
      p1d(:) = 0._rp
      !$acc end kernels
      !$acc parallel loop gang default(present) private(p1d_s) async(1)
      do i=lo(1),hi(1)
        p1d_s = 0._rp
        !$acc loop vector collapse(2) reduction(+:p1d_s)
        do k=lo(3),hi(3)
          do j=lo(2),hi(2)
            p1d_s = p1d_s + p(i,j,k)*dz(k)
          end do
        end do
        p1d(i) = p1d_s*grid_area_ratio
      end do
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(1),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      !$acc data copyin(p1d) async(1)
      !$acc kernels default(present) async(1)
      do i=lo(1),hi(1)
        p(i,:,:) = p1d(i)
      end do
      !$acc end kernels
      !$acc end data
    end select
  end subroutine ave1d_channel
  !
  subroutine ave2d_duct(ng,lo,hi,idir,l,dl,dz,p)
    !
    ! average a variable over one direction,
    ! adapted from out2d_duct
    !
    implicit none
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(lo(3)-1:) :: dz
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(dp), allocatable, dimension(:,:) :: p2d
    real(dp) :: grid_area_ratio,p2d_s
    integer :: i,j,k
    !
    select case(idir) ! streamwise direction
    case(3)
    case(2)
      grid_area_ratio = dl(2)/l(2)
      allocate(p2d(ng(1),ng(3)))
      !$acc data copyout(p2d) async(1)
      !$acc kernels default(present) async(1)
      p2d(:,:) = 0._rp
      !$acc end kernels
      !$acc parallel loop gang collapse(2) default(present) private(p2d_s) async(1)
      do k=lo(3),hi(3)
        do i=lo(1),hi(1)
          p2d_s = 0._rp
          !$acc loop vector reduction(+:p2d_s)
          do j=lo(2),hi(2)
            p2d_s = p2d_s + p(i,j,k)
          end do
          p2d(i,k) = p2d_s
        end do
      end do
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p2d(1,1),ng(1)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      !$acc data copyin(p2d) async(1)
      !$acc kernels default(present) async(1)
      do j=lo(2),hi(2)
        p(lo(1):hi(1),j,lo(3):hi(3)) = p2d(lo(1):hi(1),lo(3):hi(3))*grid_area_ratio
      end do
      !$acc end kernels
      !$acc end data
    case(1)
      grid_area_ratio = dl(1)/l(1)
      allocate(p2d(ng(2),ng(3)))
      !$acc data copyout(p2d) async(1)
      !$acc kernels default(present) async(1)
      p2d(:,:) = 0._rp
      !$acc end kernels
      !$acc parallel loop gang collapse(2) default(present) private(p2d_s) async(1)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          p2d_s = 0._rp
          !$acc loop vector reduction(+:p2d_s)
          do i=lo(1),hi(1)
            p2d_s = p2d_s + p(i,j,k)
          end do
          p2d(j,k) = p2d_s
        end do
      end do
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p2d(1,1),ng(2)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      !$acc data copyin(p2d) async(1)
      !$acc kernels default(present) async(1)
      do i=lo(1),hi(1)
        p(i,lo(2):hi(2),lo(3):hi(3)) = p2d(lo(2):hi(2),lo(3):hi(3))*grid_area_ratio
      end do
      !$acc end kernels
      !$acc end data
    end select
  end subroutine ave2d_duct
  !
  subroutine filter3d(n,p,pf)
    !
    ! 3D top-hat filter, second-order trapezoidal rule
    ! it is useless to define temporary variables to touch array elements, since
    ! each element is used only once. The current implementation is ~50% the cost of
    ! computing eight cubes (Bae's code and Davidson's book).
    !
    implicit none
    integer , intent(in ), dimension(3) :: n
    real(rp), intent(in ), dimension(0:,0:,0:) :: p
    real(rp), intent(out), dimension(0:,0:,0:) :: pf
    real(rp) :: p_mmm,p_cmm,p_pmm,p_mcm,p_ccm,p_pcm,p_mpm,p_cpm,p_ppm, &
                p_mmc,p_cmc,p_pmc,p_mcc,p_ccc,p_pcc,p_mpc,p_cpc,p_ppc, &
                p_mmp,p_cmp,p_pmp,p_mcp,p_ccp,p_pcp,p_mpp,p_cpp,p_ppp
    integer :: i,j,k,m
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 1,n(3)
      do j = 1,n(2)
        do i = 1,n(1)
          !
          p_mmm = p(i-1,j-1,k-1)
          p_cmm = p(i  ,j-1,k-1)
          p_pmm = p(i+1,j-1,k-1)
          p_mcm = p(i-1,j  ,k-1)
          p_ccm = p(i  ,j  ,k-1)
          p_pcm = p(i+1,j  ,k-1)
          p_mpm = p(i-1,j+1,k-1)
          p_cpm = p(i  ,j+1,k-1)
          p_ppm = p(i+1,j+1,k-1)
          
          p_mmc = p(i-1,j-1,k  )
          p_cmc = p(i  ,j-1,k  )
          p_pmc = p(i+1,j-1,k  )
          p_mcc = p(i-1,j  ,k  )
          p_ccc = p(i  ,j  ,k  )
          p_pcc = p(i+1,j  ,k  )
          p_mpc = p(i-1,j+1,k  )
          p_cpc = p(i  ,j+1,k  )
          p_ppc = p(i+1,j+1,k  )
          
          p_mmp = p(i-1,j-1,k+1)
          p_cmp = p(i  ,j-1,k+1)
          p_pmp = p(i+1,j-1,k+1)
          p_mcp = p(i-1,j  ,k+1)
          p_ccp = p(i  ,j  ,k+1)
          p_pcp = p(i+1,j  ,k+1)
          p_mpp = p(i-1,j+1,k+1)
          p_cpp = p(i  ,j+1,k+1)
          p_ppp = p(i+1,j+1,k+1)
          !
          pf(i,j,k) = (8._rp*(p_ccc) + &
                       4._rp*(p_mcc + p_cmc + p_ccm + &
                              p_pcc + p_cpc + p_ccp) + &
                       2._rp*(p_cmm + p_mcm + p_mmc + &
                              p_cpm + p_pcm + p_pmc + &
                              p_cmp + p_mcp + p_mpc + &
                              p_cpp + p_pcp + p_ppc) + &
                       1._rp*(p_mmm + p_pmm + p_mpm + p_ppm +&
                              p_mmp + p_pmp + p_mpp + p_ppp) &
                       )/64._rp
        end do
      end do
    end do
  end subroutine filter3d
  !
  subroutine extrapolate(n,is_bound,dzci,p,iface,cbc,lwm)
    !
    ! linear extrapolation of wall-parallel velocity
    !
    ! called before filter/strain_rate for 2D filter/one-sided difference. The
    ! extrapolation requires >= 2 off-wall layers. Wall-normal velocity is
    ! extrapolated when stored at cell centers (iface=0), but is not
    ! when stored at cell faces (iface=1,2,3). i.e., no-penetration bc.
    !
    implicit none
    integer , intent(in)   , dimension(3)        :: n
    logical , intent(in)   , dimension(0:1,3)    :: is_bound
    real(rp), intent(in)   , dimension(0:)       :: dzci
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    integer , intent(in) :: iface
    character(len=1), intent(in), dimension(0:1,3,3), optional :: cbc
    integer, intent(in), dimension(0:1,3), optional :: lwm
    logical, dimension(0:1,3) :: is_done
    real(rp) :: dzc(0:n(3)+1),factor0,factor1
    integer :: i,j,k
    !
    dzc = 1._rp/dzci
    !
    if(present(cbc)) then
      factor0 = 1._rp
      factor1 = 1._rp
      is_done(:,1) = (is_bound(:,1).and.cbc(:,1,1)=='D'.and.iface/=1)
      is_done(:,2) = (is_bound(:,2).and.cbc(:,2,2)=='D'.and.iface/=2)
      is_done(:,3) = (is_bound(:,3).and.cbc(:,3,3)=='D'.and.iface/=3)
    elseif(present(lwm)) then
      factor0 = dzc(0)*dzci(1)
      factor1 = dzc(n(3))*dzci(n(3)-1)
      is_done(:,1) = (is_bound(:,1).and.lwm(:,1)/=0.and.iface/=1)
      is_done(:,2) = (is_bound(:,2).and.lwm(:,2)/=0.and.iface/=2)
      is_done(:,3) = (is_bound(:,3).and.lwm(:,3)/=0.and.iface/=3)
    end if
    !
    if(is_done(0,1)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do k = 0,n(3)+1
        do j = 0,n(2)+1
          p(0,j,k) = 2._rp*p(1,j,k) - p(2,j,k)
        end do
      end do
    end if
    if(is_done(1,1)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do k = 0,n(3)+1
        do j = 0,n(2)+1
          p(n(1)+1,j,k) = 2._rp*p(n(1),j,k) - p(n(1)-1,j,k)
        end do
      end do
    end if
    if(is_done(0,2)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do k = 0,n(3)+1
        do i = 0,n(1)+1
          p(i,0,k) = 2._rp*p(i,1,k) - p(i,2,k)
        end do
      end do
    end if
    if(is_done(1,2)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do k = 0,n(3)+1
        do i = 0,n(1)+1
          p(i,n(2)+1,k) = 2._rp*p(i,n(2),k) - p(i,n(2)-1,k)
        end do
      end do
    end if
    if(is_done(0,3)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do j = 0,n(2)+1
        do i = 0,n(1)+1
          p(i,j,0) = (1._rp+factor0)*p(i,j,1) - factor0*p(i,j,2)
        end do
      end do
    end if
    if(is_done(1,3)) then
      !$acc parallel loop collapse(2) default(present) async(1)
      do j = 0,n(2)+1
        do i = 0,n(1)+1
          p(i,j,n(3)+1) = (1._rp+factor1)*p(i,j,n(3)) - factor1*p(i,j,n(3)-1)
        end do
      end do
    end if
  end subroutine extrapolate
  !
  subroutine cmpt_alph2(n,is_bound,cbc,alph2)
    !
    ! compute filter ratio used in the dynamic Smagorinsky model
    ! set alph2=2.52 in the first off-wall layer to yield more accurate results than 4.00
    ! in practical simulations. Specifically, the near-wall velocity profile is more
    ! accurate. The effect also depends on the grid aspect ratio, with more obvious
    ! effects at AR=1 than AR=2. The choice has negligible influence on WRLES.
    !
    implicit none
    integer, intent(in), dimension(3)        :: n
    logical, intent(in), dimension(0:1,3)    :: is_bound
    character(len=1), intent(in), dimension(0:1,3,3), optional :: cbc
    real(rp), intent(out), dimension(0:,0:,0:) :: alph2
    !
#if !defined(_FILTER_2D)
    !$acc kernels default(present) async(1)
    alph2 = 4.00_rp
    !$acc end kernels
    if(is_bound(0,1).and.cbc(0,1,1)=='D') then
      !$acc kernels default(present) async(1)
      alph2(1,:,:) = 2.52_rp
      !$acc end kernels
    end if
    if(is_bound(1,1).and.cbc(1,1,1)=='D') then
      !$acc kernels default(present) async(1)
      alph2(n(1),:,:) = 2.52_rp
      !$acc end kernels
    end if
    if(is_bound(0,2).and.cbc(0,2,2)=='D') then
      !$acc kernels default(present) async(1)
      alph2(:,1,:) = 2.52_rp
      !$acc end kernels
    end if
    if(is_bound(1,2).and.cbc(1,2,2)=='D') then
      !$acc kernels default(present) async(1)
      alph2(:,n(2),:) = 2.52_rp
      !$acc end kernels
    end if
    if(is_bound(0,3).and.cbc(0,3,3)=='D') then
      !$acc kernels default(present) async(1)
      alph2(:,:,1) = 2.52_rp
      !$acc end kernels
    end if
    if(is_bound(1,3).and.cbc(1,3,3)=='D') then
      !$acc kernels default(present) async(1)
      alph2(:,:,n(3)) = 2.52_rp
      !$acc end kernels
    end if
#else
    !$acc kernels default(present) async(1)
    alph2 = 2.52_rp
    !$acc end kernels
#endif
  end subroutine cmpt_alph2
  !
  subroutine filter2d(n,p,pf)
    !
    ! 2D top-hat filter, second-order trapezoidal rule
    ! filtering along only the homogeneous directions might be more suitable,
    ! so cs and del can be taken out of the second filtering operation,
    ! as indicated by R. Agrawal (2022)
    ! this subroutine only applies to channel flow at the moment
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(0:,0:,0:) :: p
    real(rp), intent(out), dimension(0:,0:,0:) :: pf
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k = 1,n(3)
      do j = 1,n(2)
        do i = 1,n(1)
          pf(i,j,k) = (4._rp*(p(i,j,k)) + &
                       2._rp*(p(i-1,j,k) + p(i,j-1,k) + p(i+1,j,k) + p(i,j+1,k)) + &
                       1._rp*(p(i-1,j-1,k) + p(i+1,j-1,k) + p(i-1,j+1,k) + p(i+1,j+1,k)))/16._rp
        end do
      end do
    end do
  end subroutine filter2d
  !
  subroutine interpolate(n,u,v,w,uc,vc,wc)
    !
    ! interpolate velocity to cell centers, equivalent to reconstruction (FV)
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(0:,0:,0:) :: uc,vc,wc
    integer :: i,j,k
    !
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          uc(i,j,k) = 0.5_rp*(u(i,j,k)+u(i-1,j,k))
          vc(i,j,k) = 0.5_rp*(v(i,j,k)+v(i,j-1,k))
          wc(i,j,k) = 0.5_rp*(w(i,j,k)+w(i,j,k-1))
        end do
      end do
    end do
  end subroutine interpolate
  !
  subroutine cmpt_dw_plus(cbc,n,is_bound,l,dl,dli,zc,dzc,dzci,visc,u,v,w,dw_plus)
    !
    ! inner-scaled distance to the nearest wall. We assume that a wall only
    ! affects its neighboring block, which requires that block to have enough
    ! off-wall height. Perfect partitioning has <= 2 blocks between two
    ! opposite walls. dw_plus is calculated based on minimum distance dw,
    ! instead of dw_plus, so the implementation ensures the same dw_plus
    ! under different partitionings.
    !
    ! distance to the nearest wall. Identification of walls is based on
    ! non-penetration boundary conditions, which applied to both the no-slip
    ! and free-slip (wall model) cases.
    !
    ! it is unacceptable to assume zero velocity at the wall. For no-slip walls,
    ! the velocity at the wall is zero. When a wall model is applied, tauw must 
    ! be computed using the first off-wall and ghost cells. It is incorrect to
    ! assume non-slip wall, which can lead to large errors.
    ! 
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in ), dimension(3)        :: n
    logical , intent(in ), dimension(0:1,3)    :: is_bound
    real(rp), intent(in ), dimension(3)        :: l,dl,dli
    real(rp), intent(in ), dimension(0:)       :: zc,dzc,dzci
    real(rp), intent(in )                      :: visc
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(0:,0:,0:) :: dw_plus
    real(rp), allocatable, dimension( :, :, :), save :: dw
    logical, save :: is_first = .true.
    real(rp) :: tauw(2),tauw_tot,this_dw,visci
    integer :: i,j,k
    !
    visci = 1._rp/visc
    !
    if(is_first) then
      is_first = .false.
      allocate(dw(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      !$acc enter data create(dw) async(1)
    end if
    !$acc kernels default(present) async(1)
    dw(:,:,:) = big
    !$acc end kernels
    !
    if(is_bound(0,1).and.cbc(0,1,1)=='D') then
      !$acc parallel loop collapse(3) default(present) private(tauw,tauw_tot,this_dw) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            this_dw = dl(1)*(i-0.5)
            if(this_dw < dw(i,j,k)) then
              dw(i,j,k) = this_dw
              tauw(1) = v(1,j,k)-v(0,j,k)+v(1,j-1,k)-v(0,j-1,k)
              tauw(2) = w(1,j,k)-w(0,j,k)+w(1,j,k-1)-w(0,j,k-1)
              tauw_tot= 0.5_rp*visc*sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dli(1)
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    if(is_bound(1,1).and.cbc(1,1,1)=='D') then
      !$acc parallel loop collapse(3) default(present) private(tauw,tauw_tot,this_dw) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            this_dw = dl(1)*(n(1)-i+0.5)
            if(this_dw < dw(i,j,k)) then
              dw(i,j,k) = this_dw
              tauw(1) = v(n(1),j,k)-v(n(1)+1,j,k)+v(n(1),j-1,k)-v(n(1)+1,j-1,k)
              tauw(2) = w(n(1),j,k)-w(n(1)+1,j,k)+w(n(1),j,k-1)-w(n(1)+1,j,k-1)
              tauw_tot= 0.5_rp*visc*sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))*dli(1)
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    if(is_bound(0,2).and.cbc(0,2,2)=='D') then
      !$acc parallel loop collapse(3) default(present) private(tauw,tauw_tot,this_dw) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            this_dw = dl(2)*(j-0.5)
            if(this_dw < dw(i,j,k)) then
              dw(i,j,k) = this_dw
              tauw(1) = u(i,1,k)-u(i,0,k)+u(i-1,1,k)-u(i-1,0,k)
              tauw(2) = w(i,1,k)-w(i,0,k)+w(i,1,k-1)-w(i,0,k-1)
              tauw_tot= 0.5_rp*visc*sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))*dli(2)
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    if(is_bound(1,2).and.cbc(1,2,2)=='D') then
      !$acc parallel loop collapse(3) default(present) private(tauw,tauw_tot,this_dw) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            this_dw = dl(2)*(n(2)-j+0.5)
            if(this_dw < dw(i,j,k)) then
              dw(i,j,k) = this_dw
              tauw(1) = u(i,n(2),k)-u(i,n(2)+1,k)+u(i-1,n(2),k)-u(i-1,n(2)+1,k)
              tauw(2) = w(i,n(2),k)-w(i,n(2)+1,k)+w(i,n(2),k-1)-w(i,n(2)+1,k-1)
              tauw_tot= 0.5_rp*visc*sqrt(tauw(1)*tauw(1) + tauw(2)*tauw(2))*dli(2)
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    if(is_bound(0,3).and.cbc(0,3,3)=='D') then
      !$acc parallel loop collapse(3) default(present) private(tauw,tauw_tot,this_dw) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            this_dw = zc(k)
            if(this_dw < dw(i,j,k)) then
              dw(i,j,k) = this_dw
              tauw(1) = u(i,j,1)-u(i,j,0)+u(i-1,j,1)-u(i-1,j,0)
              tauw(2) = v(i,j,1)-v(i,j,0)+v(i,j-1,1)-v(i,j-1,0)
              tauw_tot= 0.5_rp*visc*sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dzci(0)
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
    if(is_bound(1,3).and.cbc(1,3,3)=='D') then
      !$acc parallel loop collapse(3) default(present) private(tauw,tauw_tot,this_dw) async(1)
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            this_dw = l(3)-zc(k)
            if(this_dw < dw(i,j,k)) then
              dw(i,j,k) = this_dw
              tauw(1) = u(i,j,n(3))-u(i,j,n(3)+1)+u(i-1,j,n(3))-u(i-1,j,n(3)+1)
              tauw(2) = v(i,j,n(3))-v(i,j,n(3)+1)+v(i,j-1,n(3))-v(i,j-1,n(3)+1)
              tauw_tot= 0.5_rp*visc*sqrt(tauw(1)*tauw(1)+tauw(2)*tauw(2))*dzci(n(3))
              dw_plus(i,j,k) = this_dw*sqrt(tauw_tot)*visci
            end if
          end do
        end do
      end do
    end if
  end subroutine cmpt_dw_plus
  !
  subroutine strain_rate(n,dli,dzci,dzfi,u,v,w,s0,sij)
    !
    ! compute the strain rate field
    !
    ! Sij is first computed at (or averaged to) cell center, then s0=sqrt(2SijSij)
    ! at cell center (Bae and Orlandi's codes). The implementation is efficient, since
    ! it avoids repetitive computation of derivatives. Costa first averages SijSij to
    ! cell center, then computes s0. This leads to larger s0, especially when Sij
    ! at the cell edges have opposite signs. The second and third loops cannot combine.
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(3)        :: dli
    real(rp), intent(in ), dimension(0:)       :: dzci,dzfi
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(0:,0:,0:) :: s0
    real(rp), intent(out), dimension(0:,0:,0:,1:), optional :: sij
    real(rp) :: s11,s22,s33,s12,s13,s23
    real(rp) :: dxi,dyi
    integer :: i,j,k
    real(rp) :: u_mcm,u_ccm,u_mmc,u_cmc,u_mcc,u_ccc,u_mpc,u_cpc,u_mcp,u_ccp, &
                v_cmm,v_ccm,v_mmc,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_cmp,v_ccp, &
                w_cmm,w_mcm,w_ccm,w_pcm,w_cpm,w_cmc,w_mcc,w_ccc,w_pcc,w_cpc
    !
    dxi = dli(1)
    dyi = dli(2)
    !
    !$acc parallel loop collapse(3) default(present) async(1) &
    !$acc private(u_mcm,u_ccm,u_mmc,u_cmc,u_mcc,u_ccc,u_mpc,u_cpc,u_mcp,u_ccp) &
    !$acc private(v_cmm,v_ccm,v_mmc,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_cmp,v_ccp) &
    !$acc private(w_cmm,w_mcm,w_ccm,w_pcm,w_cpm,w_cmc,w_mcc,w_ccc,w_pcc,w_cpc) &
    !$acc private(s11,s22,s33,s12,s13,s23)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !
          u_mcm = u(i-1,j  ,k-1)
          u_ccm = u(i  ,j  ,k-1)
          u_mmc = u(i-1,j-1,k  )
          u_cmc = u(i  ,j-1,k  )
          u_mcc = u(i-1,j  ,k  )
          u_ccc = u(i  ,j  ,k  )
          u_mpc = u(i-1,j+1,k  )
          u_cpc = u(i  ,j+1,k  )
          u_mcp = u(i-1,j  ,k+1)
          u_ccp = u(i  ,j  ,k+1)
          !
          v_cmm = v(i  ,j-1,k-1)
          v_ccm = v(i  ,j  ,k-1)
          v_mmc = v(i-1,j-1,k  )
          v_cmc = v(i  ,j-1,k  )
          v_pmc = v(i+1,j-1,k  )
          v_mcc = v(i-1,j  ,k  )
          v_ccc = v(i  ,j  ,k  )
          v_pcc = v(i+1,j  ,k  )
          v_cmp = v(i  ,j-1,k+1)
          v_ccp = v(i  ,j  ,k+1)
          !
          w_cmm = w(i  ,j-1,k-1)
          w_mcm = w(i-1,j  ,k-1)
          w_ccm = w(i  ,j  ,k-1)
          w_pcm = w(i+1,j  ,k-1)
          w_cpm = w(i  ,j+1,k-1)
          w_cmc = w(i  ,j-1,k  )
          w_mcc = w(i-1,j  ,k  )
          w_ccc = w(i  ,j  ,k  )
          w_pcc = w(i+1,j  ,k  )
          w_cpc = w(i  ,j+1,k  )
          !
          s11 = (u_ccc-u_mcc)*dxi
          s22 = (v_ccc-v_cmc)*dyi
          s33 = (w_ccc-w_ccm)*dzfi(k)
          s12 = .125_rp*((u_cpc-u_ccc)*dyi + (v_pcc-v_ccc)*dxi + &
                         (u_ccc-u_cmc)*dyi + (v_pmc-v_cmc)*dxi + &
                         (u_mpc-u_mcc)*dyi + (v_ccc-v_mcc)*dxi + &
                         (u_mcc-u_mmc)*dyi + (v_cmc-v_mmc)*dxi)
          s13 = .125_rp*((u_ccp-u_ccc)*dzci(k  ) + (w_pcc-w_ccc)*dxi + &
                         (u_ccc-u_ccm)*dzci(k-1) + (w_pcm-w_ccm)*dxi + &
                         (u_mcp-u_mcc)*dzci(k  ) + (w_ccc-w_mcc)*dxi + &
                         (u_mcc-u_mcm)*dzci(k-1) + (w_ccm-w_mcm)*dxi)
          s23 = .125_rp*((v_ccp-v_ccc)*dzci(k  ) + (w_cpc-w_ccc)*dyi + &
                         (v_ccc-v_ccm)*dzci(k-1) + (w_cpm-w_ccm)*dyi + &
                         (v_cmp-v_cmc)*dzci(k  ) + (w_ccc-w_cmc)*dyi + &
                         (v_cmc-v_cmm)*dzci(k-1) + (w_ccm-w_cmm)*dyi)
          s0(i,j,k) = sqrt(2._rp*(s11**2+s22**2+s33**2+2._rp*(s12**2+s13**2+s23**2)))
          if(present(sij)) then
            sij(i,j,k,1:6) = (/s11,s22,s33,s12,s13,s23/)
          end if
        end do
      end do
    end do
  end subroutine strain_rate
end module mod_sgs