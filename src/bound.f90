! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_bound
  use mpi
  use mod_common_mpi, only: ierr,halo,ipencil_axis
  use mod_precision, only: rp,sp,dp,i8,MPI_REAL_RP
  use mod_typedef   , only: bound
  use mod_wallmodel    , only: updt_wallmodelbc
  implicit none
  private
  public boundp,bounduvw,cmpt_rhs_b,updt_rhs_b,initbc
  contains
  subroutine bounduvw(cbc,n,bcu,bcv,bcw,bcu_mag,bcv_mag,bcw_mag,nb,is_bound,lwm,l,dl,zc,zf,dzc,dzf, &
                      visc,h,index_wm,is_updt_wm,is_correc,u,v,w)
    !
    ! imposes velocity boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer         , intent(in), dimension(3) :: n
    type(bound)     , intent(inout) :: bcu,bcv,bcw
    type(bound)     , intent(in) :: bcu_mag,bcv_mag,bcw_mag
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    integer , intent(in), dimension(0:1,3) :: lwm,index_wm
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc,zf,dzc,dzf
    real(rp), intent(in) :: visc,h
    logical , intent(in) :: is_updt_wm,is_correc
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    character(len=1), dimension(0:1,3,3) :: cbc_w
    logical :: impose_norm_bc
    integer :: i,j,k,idir,nh
    !
    nh = 1
    !
#if !defined(_OPENACC)
    do idir = 1,3
      call updthalo(nh,halo(idir),nb(:,idir),idir,u)
      call updthalo(nh,halo(idir),nb(:,idir),idir,v)
      call updthalo(nh,halo(idir),nb(:,idir),idir,w)
    end do
#else
    call updthalo_gpu(nh,cbc(0,:,1)//cbc(1,:,1)==['PP','PP','PP'],u)
    call updthalo_gpu(nh,cbc(0,:,2)//cbc(1,:,2)==['PP','PP','PP'],v)
    call updthalo_gpu(nh,cbc(0,:,3)//cbc(1,:,3)==['PP','PP','PP'],w)
#endif
    !
    ! impose_norm_bc=0, the correction does not change wall-normal velocity on walls,
    ! this also holds when a wall model is used due to non-penetrating bc
    impose_norm_bc = (.not.is_correc).or.(cbc(0,1,1)//cbc(1,1,1) == 'PP')
    if(is_bound(0,1)) then
      if(impose_norm_bc)  call set_bc(cbc(0,1,1),0,1,nh,.false.,bcu%x,dl(1),u)
      if(lwm(0,1)==0) then
                          call set_bc(cbc(0,1,2),0,1,nh,.true. ,bcv%x,dl(1),v)
                          call set_bc(cbc(0,1,3),0,1,nh,.true. ,bcw%x,dl(1),w)
      end if
    end if
    if(is_bound(1,1)) then
      if(impose_norm_bc)  call set_bc(cbc(1,1,1),1,1,nh,.false.,bcu%x,dl(1),u)
      if(lwm(1,1)==0) then    
                          call set_bc(cbc(1,1,2),1,1,nh,.true. ,bcv%x,dl(1),v)
                          call set_bc(cbc(1,1,3),1,1,nh,.true. ,bcw%x,dl(1),w)
      end if
    end if
    impose_norm_bc = (.not.is_correc).or.(cbc(0,2,2)//cbc(1,2,2) == 'PP')
    if(is_bound(0,2)) then
      if(impose_norm_bc)  call set_bc(cbc(0,2,2),0,2,nh,.false.,bcv%y,dl(2),v)
      if(lwm(0,2)==0) then    
                          call set_bc(cbc(0,2,1),0,2,nh,.true. ,bcu%y,dl(2),u)
                          call set_bc(cbc(0,2,3),0,2,nh,.true. ,bcw%y,dl(2),w)
      end if
    end if
    if(is_bound(1,2)) then
      if(impose_norm_bc)  call set_bc(cbc(1,2,2),1,2,nh,.false.,bcv%y,dl(2),v)
      if(lwm(1,2)==0) then    
                          call set_bc(cbc(1,2,1),1,2,nh,.true. ,bcu%y,dl(2),u)
                          call set_bc(cbc(1,2,3),1,2,nh,.true. ,bcw%y,dl(2),w)
      end if
    end if
    impose_norm_bc = (.not.is_correc).or.(cbc(0,3,3)//cbc(1,3,3) == 'PP')
    if(is_bound(0,3)) then
      if(impose_norm_bc)  call set_bc(cbc(0,3,3),0,3,nh,.false.,bcw%z,dzf(0)   ,w)
      if(lwm(0,3)==0) then     
                          call set_bc(cbc(0,3,1),0,3,nh,.true. ,bcu%z,dzc(0)   ,u)
                          call set_bc(cbc(0,3,2),0,3,nh,.true. ,bcv%z,dzc(0)   ,v)
      end if
    end if
    if(is_bound(1,3)) then
      if(impose_norm_bc)  call set_bc(cbc(1,3,3),1,3,nh,.false.,bcw%z,dzf(n(3)),w)
      if(lwm(1,3)==0) then    
                          call set_bc(cbc(1,3,1),1,3,nh,.true. ,bcu%z,dzc(n(3)),u)
                          call set_bc(cbc(1,3,2),1,3,nh,.true. ,bcv%z,dzc(n(3)),v)
      end if
    end if
    !
    ! Neumann bc must be zero at some locations to let the ghost points have the same
    ! values from the bc implementations of different walls. For example, w(1,5,4)
    ! must be set as zero by both the upper wall (non-penetrating) and the front wall
    ! (Neumann). Wall models can guarantee this at those locations, since they do
    ! produce a zero gradient in the direction where the velocity component is zero.
    ! 
    ! bcu/v/w have correct values for all the loop locations in updt_wallmodelbc.
    ! The only difference is that for square duct/six-wall cases, uh,vh and wh may
    ! have incorrect values at the end locations. However, it is only the zero-value
    ! component (wall-normal) that is essentially used, which guarantees
    ! zero-gradient bc there, leading to correct values at the end ghost points.
    ! For two-wall channels, uh,vh and wh do have correct values at the end locations.
    ! Therefore, it makes no difference whether the wall-parallel velocity at the wall
    ! is set to zero or not, as it is never utilized.
    !
    ! updt_wallmodelbc, ~25% of bounduvw time is saved, equivalent to ~1% of the total time,
    ! The computational cost of the log-law wall model is negligible.
    !
    if(is_updt_wm) then
      call updt_wallmodelbc(n,is_bound,lwm,l,dl,zc,zf,dzc,dzf,visc,h,index_wm,u,v,w, &
                            bcu,bcv,bcw,bcu_mag,bcv_mag,bcw_mag)
    end if
    !
    if(is_bound(0,1).and.lwm(0,1)/=0) then
      call set_bc(cbc(0,1,2),0,1,nh,.true. ,bcv%x,dl(1),v)
      call set_bc(cbc(0,1,3),0,1,nh,.true. ,bcw%x,dl(1),w)
    end if
    if(is_bound(1,1).and.lwm(1,1)/=0) then
      call set_bc(cbc(1,1,2),1,1,nh,.true. ,bcv%x,dl(1),v)
      call set_bc(cbc(1,1,3),1,1,nh,.true. ,bcw%x,dl(1),w)
    end if
    if(is_bound(0,2).and.lwm(0,2)/=0) then
      call set_bc(cbc(0,2,1),0,2,nh,.true. ,bcu%y,dl(2),u)
      call set_bc(cbc(0,2,3),0,2,nh,.true. ,bcw%y,dl(2),w)
    end if
    if(is_bound(1,2).and.lwm(1,2)/=0) then
      call set_bc(cbc(1,2,1),1,2,nh,.true. ,bcu%y,dl(2),u)
      call set_bc(cbc(1,2,3),1,2,nh,.true. ,bcw%y,dl(2),w)
    end if
    if(is_bound(0,3).and.lwm(0,3)/=0) then
      call set_bc(cbc(0,3,1),0,3,nh,.true. ,bcu%z,dzc(0)   ,u)
      call set_bc(cbc(0,3,2),0,3,nh,.true. ,bcv%z,dzc(0)   ,v)
    end if
    if(is_bound(1,3).and.lwm(1,3)/=0) then
      call set_bc(cbc(1,3,1),1,3,nh,.true. ,bcu%z,dzc(n(3)),u)
      call set_bc(cbc(1,3,2),1,3,nh,.true. ,bcv%z,dzc(n(3)),v)
    end if
    !
    ! for square duct/six-wall cases, further add a loop for the end locations
    ! if all values of uh,vh and wh need to be correct at the end locations.
    ! This need should never arise.
    !
  end subroutine bounduvw
  !
  subroutine boundp(cbc,n,bcp,nb,is_bound,dl,dzc,p)
    !
    ! imposes boundary conditions for cell-centered variables (p, pp,visct)
    ! ghost cells at corners also updated, e.g. p(0,0,1), p(0,0,0)
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer         , intent(in), dimension(3) :: n
    type(bound)     , intent(in) :: bcp
    integer , intent(in), dimension(0:1,3) :: nb
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(3 ) :: dl
    real(rp), intent(in), dimension(0:) :: dzc
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    integer :: idir,nh
    !
    nh = 1
    !
#if !defined(_OPENACC)
    do idir = 1,3
      call updthalo(nh,halo(idir),nb(:,idir),idir,p)
    end do
#else
    call updthalo_gpu(nh,cbc(0,:)//cbc(1,:)==['PP','PP','PP'],p)
#endif
    !
    if(is_bound(0,1)) then
      call set_bc(cbc(0,1),0,1,nh,.true.,bcp%x,dl(1),p)
    end if
    if(is_bound(1,1)) then
      call set_bc(cbc(1,1),1,1,nh,.true.,bcp%x,dl(1),p)
    end if
    if(is_bound(0,2)) then
      call set_bc(cbc(0,2),0,2,nh,.true.,bcp%y,dl(2),p)
     end if
    if(is_bound(1,2)) then
      call set_bc(cbc(1,2),1,2,nh,.true.,bcp%y,dl(2),p)
    end if
    if(is_bound(0,3)) then
      call set_bc(cbc(0,3),0,3,nh,.true.,bcp%z,dzc(0)   ,p)
    end if
    if(is_bound(1,3)) then
      call set_bc(cbc(1,3),1,3,nh,.true.,bcp%z,dzc(n(3)),p)
    end if
  end subroutine boundp
  !
  subroutine set_bc(ctype,ibound,idir,nh,centered,bc,dr,p)
    implicit none
    character(len=1), intent(in) :: ctype
    integer , intent(in) :: ibound,idir,nh
    logical , intent(in) :: centered
    real(rp), intent(in), dimension(1-nh:,1-nh:,0:) :: bc ! bc, two faces
    real(rp), intent(in) :: dr
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    real(rp) :: sgn
    integer  :: n,dh,n1,n2
    !
    n = size(p,idir)-2*nh
    sgn = 1._rp
    if(ctype == 'D'.and.centered) then
      sgn = -1._rp
    end if
    do dh=0,nh-1
      select case(ctype)
      case('P')
        !
        ! n.b.: this periodic BC imposition assumes that the subroutine is only called for the
        !       non-decomposed directions, for which n is the domain length in index space;
        !       note that the is_bound(:,:) mask above (set under initmpi.f90) is only true along
        !       the (undecomposed) pencil direction;
        !       along decomposed directions, periodicity is naturally set via the halo exchange
        !
        ! wall normal velocity at k=n+1 is never used, it can be set to any values. It must be 
        ! problematic if it is used, since we do not have a corresponding value at k=-1. This 
        ! keeps true even when there is SGS/wall model/additional viscous term.
        !
        select case(idir)
        case(1)
          !$acc kernels default(present) async(1)
          p(  0-dh,:,:) = p(n-dh,:,:)
          p(n+1+dh,:,:) = p(1+dh,:,:)
          !$acc end kernels
        case(2)
          !$acc kernels default(present) async(1)
          p(:,  0-dh,:) = p(:,n-dh,:)
          p(:,n+1+dh,:) = p(:,1+dh,:)
          !$acc end kernels
        case(3)
          !$acc kernels default(present) async(1)
          p(:,:,  0-dh) = p(:,:,n-dh)
          p(:,:,n+1+dh) = p(:,:,1+dh)
          !$acc end kernels
        end select
      case('D')
        if(centered) then
          select case(idir)
          case(1)
            if     (ibound == 0) then
              !$acc kernels default(present) async(1)
              p(  0-dh,:,:) = 2._rp*bc(:,:,ibound)+sgn*p(1+dh,:,:)
              !$acc end kernels
            else if(ibound == 1) then
              !$acc kernels default(present) async(1)
              p(n+1+dh,:,:) = 2._rp*bc(:,:,ibound)+sgn*p(n-dh,:,:)
              !$acc end kernels
            end if
          case(2)
            if     (ibound == 0) then
              !$acc kernels default(present) async(1)
              p(:,  0-dh,:) = 2._rp*bc(:,:,ibound)+sgn*p(:,1+dh,:)
              !$acc end kernels
            else if(ibound == 1) then
              !$acc kernels default(present) async(1)
              p(:,n+1+dh,:) = 2._rp*bc(:,:,ibound)+sgn*p(:,n-dh,:)
              !$acc end kernels
            end if
          case(3)
            if     (ibound == 0) then
              !$acc kernels default(present) async(1)
              p(:,:,  0-dh) = 2._rp*bc(:,:,ibound)+sgn*p(:,:,1+dh)
              !$acc end kernels
            else if(ibound == 1) then
              !$acc kernels default(present) async(1)
              p(:,:,n+1+dh) = 2._rp*bc(:,:,ibound)+sgn*p(:,:,n-dh)
              !$acc end kernels
            end if
          end select
        else if(.not.centered) then
          select case(idir)
          case(1)
            if     (ibound == 0) then
              !$acc kernels default(present) async(1)
              p(0-dh,:,:) = bc(:,:,ibound)
              !$acc end kernels
            else if(ibound == 1) then
              !$acc kernels default(present) async(1)
              p(n+1 ,:,:) = p(n-1,:,:) ! unused
              p(n+dh,:,:) = bc(:,:,ibound)
              !$acc end kernels
            end if
          case(2)
            if     (ibound == 0) then
              !$acc kernels default(present) async(1)
              p(:,0-dh,:) = bc(:,:,ibound)
              !$acc end kernels
            else if(ibound == 1) then
              !$acc kernels default(present) async(1)
              p(:,n+1 ,:) = p(:,n-1,:) ! unused
              p(:,n+dh,:) = bc(:,:,ibound)
              !$acc end kernels
            end if
          case(3)
            if     (ibound == 0) then
              !$acc kernels default(present) async(1)
              p(:,:,0-dh) = bc(:,:,ibound)
              !$acc end kernels
            else if(ibound == 1) then
              !$acc kernels default(present) async(1)
              p(:,:,n+1 ) = p(:,:,n-1) ! unused
              p(:,:,n+dh) = bc(:,:,ibound)
              !$acc end kernels
            end if
          end select
        end if
      case('N')
        if(centered) then
          select case(idir)
          case(1)
            if     (ibound == 0) then
              !$acc kernels default(present) async(1)
              p(  0-dh,:,:) = -dr*bc(:,:,ibound)+sgn*p(1+dh,:,:)
              !$acc end kernels
            else if(ibound == 1) then
              !$acc kernels default(present) async(1)
              p(n+1+dh,:,:) = dr*bc(:,:,ibound)+sgn*p(n-dh,:,:)
              !$acc end kernels
            end if
          case(2)
            if     (ibound == 0) then
              !$acc kernels default(present) async(1)
              p(:,  0-dh,:) = -dr*bc(:,:,ibound)+sgn*p(:,1+dh,:)
              !$acc end kernels
            else if(ibound == 1) then
              !$acc kernels default(present) async(1)
              p(:,n+1+dh,:) = dr*bc(:,:,ibound)+sgn*p(:,n-dh,:)
              !$acc end kernels
            end if
          case(3)
            if     (ibound == 0) then
              !$acc kernels default(present) async(1)
              p(:,:,  0-dh) = -dr*bc(:,:,ibound)+sgn*p(:,:,1+dh)
              !$acc end kernels
            else if(ibound == 1) then
              !$acc kernels default(present) async(1)
              p(:,:,n+1+dh) = dr*bc(:,:,ibound)+sgn*p(:,:,n-dh)
              !$acc end kernels
            end if
          end select
        else if(.not.centered) then
          select case(idir)
          case(1)
            if     (ibound == 0) then
              !$acc kernels default(present) async(1)
              !p(0,:,:) = 1./3.*(-2.*factor+4.*p(1  ,:,:)-p(2  ,:,:)) ! second-order approximation of the first derivative at the boundary
              p(0-dh,:,:) = -dr*bc(:,:,ibound) + p(  1+dh,:,:)
              !$acc end kernels
            else if(ibound == 1) then
              !$acc kernels default(present) async(1)
              !p(n,:,:) = 1./3.*(-2.*factor+4.*p(n-1,:,:)-p(n-2,:,:))
              p(n+1 ,:,:) = p(n,:,:) ! unused
              p(n+dh,:,:) = dr*bc(:,:,ibound) + p(n-1-dh,:,:)
              !$acc end kernels
            end if
          case(2)
            if     (ibound == 0) then
              !$acc kernels default(present) async(1)
              !p(:,0  ,:) = 1./3.*(-2.*factor+4.*p(:,1,:)-p(:,2  ,:))
              p(:,0-dh,:) = -dr*bc(:,:,ibound) + p(:,  1+dh,:)
              !$acc end kernels
            else if(ibound == 1) then
              !$acc kernels default(present) async(1)
              !p(:,n,:) = 1./3.*(-2.*factor+4.*p(:,n-1,:)-p(:,n-2,:))
              p(:,n+1 ,:) = p(:,n,:) ! unused
              p(:,n+dh,:) = dr*bc(:,:,ibound) + p(:,n-1-dh,:)
              !$acc end kernels
            end if
          case(3)
            if     (ibound == 0) then
              !$acc kernels default(present) async(1)
              !p(:,:,0) = 1./3.*(-2.*factor+4.*p(:,:,1  )-p(:,:,2  ))
              p(:,:,0-dh) = -dr*bc(:,:,ibound) + p(:,:,  1+dh)
              !$acc end kernels
            else if(ibound == 1) then
              !$acc kernels default(present) async(1)
              !p(:,:,n) = 1./3.*(-2.*factor+4.*p(:,:,n-1)-p(:,:,n-2))
              p(:,:,n+1 ) = p(:,:,n) ! unused
              p(:,:,n+dh) = dr*bc(:,:,ibound) + p(:,:,n-1-dh)
              !$acc end kernels
            end if
          end select
        end if
      end select
    end do
  end subroutine set_bc
  !
  subroutine inflow(idir,is_bound,vel2d,u,v,w)
    implicit none
    integer , intent(in   )  :: idir
    logical , intent(in   ), dimension(0:1,3) :: is_bound
    real(rp), intent(in   ), dimension(0:,0:   ) :: vel2d
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    integer :: i,j,k
    integer, dimension(3) :: n
    !
    select case(idir)
      case(1) ! x direction
        if(is_bound(0,1)) then
          n(:) = shape(u) - 2*1
          i = 0
          !$acc parallel loop collapse(2) default(present) async(1)
          do k=1,n(3)
            do j=1,n(2)
              u(i,j,k) = vel2d(j,k)
            end do
          end do
        end if
      case(2) ! y direction
        if(is_bound(0,2)) then
          n(:) = shape(v) - 2*1
          j = 0
          !$acc parallel loop collapse(2) default(present) async(1)
          do k=1,n(3)
            do i=1,n(1)
              v(i,j,k) = vel2d(i,k)
            end do
          end do
        end if
      case(3) ! z direction
        if(is_bound(0,3)) then
          n(:) = shape(w) - 2*1
          k = 0
          !$acc parallel loop collapse(2) default(present) async(1)
          do j=1,n(2)
            do i=1,n(1)
              w(i,j,k) = vel2d(i,j)
            end do
          end do
        end if
    end select
  end subroutine inflow
  !
  subroutine cmpt_rhs_b(ng,dl,dzc,dzf,cbc,bc,c_or_f,rhsbx,rhsby,rhsbz)
    !
    ! compute values added to the right hand side
    !
    implicit none
    integer , intent(in), dimension(3) :: ng
    real(rp), intent(in), dimension(3 ) :: dl
    real(rp), intent(in), dimension(0:) :: dzc,dzf
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    type(bound)     , intent(in) :: bc
    character(len=1), intent(in), dimension(3) :: c_or_f
    real(rp), intent(out), dimension(:,:,0:), optional :: rhsbx
    real(rp), intent(out), dimension(:,:,0:), optional :: rhsby
    real(rp), intent(out), dimension(:,:,0:), optional :: rhsbz
    real(rp), allocatable, dimension(:)     , save     :: dxc01,dxf01,dyc01,dyf01, &
                                                          dzc01_c,dzf01_c,dzc01_f,dzf01_f
    logical, save :: is_first = .true.
    !
    if(is_first) then
      is_first = .false.
      allocate(dxc01  (0:1), &
               dxf01  (0:1), &
               dyc01  (0:1), &
               dyf01  (0:1), &
               dzc01_c(0:1), &
               dzf01_c(0:1), &
               dzc01_f(0:1), &
               dzf01_f(0:1))
      dxc01 = [dl(1) ,dl(1)]
      dxf01 = [dl(1) ,dl(1)]
      dyc01 = [dl(2) ,dl(2)]
      dyf01 = [dl(2) ,dl(2)]
      dzc01_c = [dzc(0),dzc(ng(3)  )]
      dzf01_c = [dzf(1),dzf(ng(3)  )]
      dzc01_f = [dzc(1),dzc(ng(3)-1)]
      dzf01_f = [dzf(1),dzf(ng(3)  )]
      !$acc enter data copyin(dxc01,dxf01,dyc01,dyf01,dzc01_c,dzf01_c,dzc01_f,dzf01_f) async(1)
    end if
    !
    if(present(rhsbx)) then
      call bc_rhs(cbc(:,1),bc%x,dxc01  ,dxf01,c_or_f(1),rhsbx) ! x-direction
    end if
    if(present(rhsby)) then
      call bc_rhs(cbc(:,2),bc%y,dyc01  ,dyf01,c_or_f(2),rhsby) ! y-direction
    end if
    if(present(rhsbz)) then
      if(c_or_f(3) == 'c') then
        call bc_rhs(cbc(:,3),bc%z,dzc01_c,dzf01_c,c_or_f(3),rhsbz) ! z-direction
      else if(c_or_f(3) == 'f') then
        call bc_rhs(cbc(:,3),bc%z,dzc01_f,dzf01_f,c_or_f(3),rhsbz) ! z-direction
      end if
    end if
  end subroutine cmpt_rhs_b
  !
  subroutine bc_rhs(cbc,bc,dlc,dlf,c_or_f,rhs)
    implicit none
    character(len=1), intent(in), dimension(0:1) :: cbc
    real(rp), intent(in), dimension(0:,0:,0:) :: bc ! bc, two faces in a direction
    real(rp), intent(in), dimension(0:1) :: dlc,dlf
    real(rp), intent(out), dimension(:,:,0:) :: rhs
    character(len=1), intent(in) :: c_or_f ! c -> cell-centered; f -> face-centered
    real(rp) :: sgn
    integer :: ibound,n1,n2
    !
    n1 = size(bc,1) - 2
    n2 = size(bc,2) - 2
    !
    select case(c_or_f)
    case('c')
      do ibound = 0,1
        select case(cbc(ibound))
        case('P')
          !$acc kernels default(present) async(1)
          rhs(:,:,ibound) =  0._rp
          !$acc end kernels
        case('D')
          !$acc kernels default(present) async(1)
          rhs(:,:,ibound) = -2._rp*bc(1:n1,1:n2,ibound)/dlc(ibound)/dlf(ibound)
          !$acc end kernels
        case('N')
          if(ibound == 0) then
            sgn =  1._rp
          else
            sgn = -1._rp
          end if
          !$acc kernels default(present) async(1)
          rhs(:,:,ibound) = sgn*bc(1:n1,1:n2,ibound)/dlf(ibound)
          !$acc end kernels
        end select
      end do
    case('f')
      do ibound = 0,1
        select case(cbc(ibound))
        case('P')
          !$acc kernels default(present) async(1)
          rhs(:,:,ibound) =  0._rp
          !$acc end kernels
        case('D')
          !$acc kernels default(present) async(1)
          rhs(:,:,ibound) = -bc(1:n1,1:n2,ibound)/dlc(ibound)/dlf(ibound)
          !$acc end kernels
        case('N')
          if(ibound == 0) then
            sgn =  1._rp
          else
            sgn = -1._rp
          end if
          !$acc kernels default(present) async(1)
          rhs(:,:,ibound) = sgn*bc(1:n1,1:n2,ibound)/dlc(ibound)
          !$acc end kernels
        end select
      end do
    end select
  end subroutine bc_rhs
  !
  subroutine updt_rhs_b(c_or_f,cbc,n,is_bound,rhsbx,rhsby,rhsbz,p)
    implicit none
    character(len=1), intent(in), dimension(3    ) :: c_or_f
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(:,:,0:), optional :: rhsbx,rhsby,rhsbz
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    integer , dimension(3) :: q
    integer :: idir
    integer :: nn
    q(:) = 0
    do idir = 1,3
      if(c_or_f(idir) == 'f'.and.cbc(1,idir) == 'D') q(idir) = 1
    end do
    !
    if(present(rhsbx)) then
      if(is_bound(0,1)) then
        !$acc kernels default(present) async(1)
        p(1 ,1:n(2),1:n(3)) = p(1 ,1:n(2),1:n(3)) + rhsbx(:,:,0)
        !$acc end kernels
      end if
      if(is_bound(1,1)) then
        nn = n(1)-q(1)
        !$acc kernels default(present) async(1)
        p(nn,1:n(2),1:n(3)) = p(nn,1:n(2),1:n(3)) + rhsbx(:,:,1)
        !$acc end kernels
      end if
    end if
    if(present(rhsby)) then
      if(is_bound(0,2)) then
        !$acc kernels default(present) async(1)
        p(1:n(1),1 ,1:n(3)) = p(1:n(1),1 ,1:n(3)) + rhsby(:,:,0)
        !$acc end kernels
      end if
      if(is_bound(1,2)) then
        nn = n(2)-q(2)
        !$acc kernels default(present) async(1)
        p(1:n(1),nn,1:n(3)) = p(1:n(1),nn,1:n(3)) + rhsby(:,:,1)
        !$acc end kernels
      end if
    end if
    if(present(rhsbz)) then
      if(is_bound(0,3)) then
        !$acc kernels default(present) async(1)
        p(1:n(1),1:n(2),1 ) = p(1:n(1),1:n(2),1 ) + rhsbz(:,:,0)
        !$acc end kernels
      end if
      if(is_bound(1,3)) then
        nn = n(3)-q(3)
        !$acc kernels default(present) async(1)
        p(1:n(1),1:n(2),nn) = p(1:n(1),1:n(2),nn) + rhsbz(:,:,1)
        !$acc end kernels
      end if
    end if
  end subroutine updt_rhs_b
  !
  subroutine updthalo(nh,halo,nb,idir,p)
    implicit none
    integer , intent(in) :: nh ! number of ghost points
    integer , intent(in) :: halo
    integer , intent(in), dimension(0:1) :: nb
    integer , intent(in) :: idir
    real(rp), dimension(1-nh:,1-nh:,1-nh:), intent(inout) :: p
    integer , dimension(3) :: lo,hi
#if defined(_ASYNC_HALO)
    integer :: requests(4)
#endif
    !
    !  this subroutine updates the halo that store info
    !  from the neighboring computational sub-domain
    !
    if(idir == ipencil_axis) return
    lo(:) = lbound(p)+nh
    hi(:) = ubound(p)-nh
    select case(idir)
    case(1) ! x direction
#if !defined(_ASYNC_HALO)
      call MPI_SENDRECV(p(lo(1)     ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                        p(hi(1)+1   ,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call MPI_SENDRECV(p(hi(1)-nh+1,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                        p(lo(1)-nh  ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
#else
      call MPI_IRECV( p(hi(1)+1  ,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                      MPI_COMM_WORLD,requests(1),ierr)
      call MPI_IRECV( p(lo(1)-nh ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),1, &
                      MPI_COMM_WORLD,requests(2),ierr)
      call MPI_ISEND(p(lo(1)     ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                      MPI_COMM_WORLD,requests(3),ierr)
      call MPI_ISEND(p(hi(1)-nh+1,lo(2)-nh,lo(3)-nh),1,halo,nb(1),1, &
                      MPI_COMM_WORLD,requests(4),ierr)
      call MPI_WAITALL(4,requests,MPI_STATUSES_IGNORE,ierr)
#endif
    case(2) ! y direction
#if !defined(_ASYNC_HALO)
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)     ,lo(3)-nh),1,halo,nb(0),0, &
                        p(lo(1)-nh,hi(2)+1   ,lo(3)-nh),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call MPI_SENDRECV(p(lo(1)-nh,hi(2)-nh+1,lo(3)-nh),1,halo,nb(1),0, &
                        p(lo(1)-nh,lo(2)-nh  ,lo(3)-nh),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
#else
      call MPI_IRECV(p(lo(1)-nh,hi(2)+1   ,lo(3)-nh),1,halo,nb(1),0, &
                      MPI_COMM_WORLD,requests(1),ierr)
      call MPI_IRECV(p(lo(1)-nh,lo(2)-nh  ,lo(3)-nh),1,halo,nb(0),1, &
                      MPI_COMM_WORLD,requests(2),ierr)
      call MPI_ISEND(p(lo(1)-nh,lo(2)     ,lo(3)-nh),1,halo,nb(0),0, &
                      MPI_COMM_WORLD,requests(3),ierr)
      call MPI_ISEND(p(lo(1)-nh,hi(2)-nh+1,lo(3)-nh),1,halo,nb(1),1, &
                      MPI_COMM_WORLD,requests(4),ierr)
      call MPI_WAITALL(4,requests,MPI_STATUSES_IGNORE,ierr)
#endif
    case(3) ! z direction
#if !defined(_ASYNC_HALO)
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)-nh,lo(3)     ),1,halo,nb(0),0, &
                        p(lo(1)-nh,lo(2)-nh,hi(3)+1   ),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)-nh,hi(3)-nh+1),1,halo,nb(1),0, &
                        p(lo(1)-nh,lo(2)-nh,lo(3)-nh  ),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
#else
      call MPI_IRECV(p(lo(1)-nh,lo(2)-nh,hi(3)+1   ),1,halo,nb(1),0, &
                      MPI_COMM_WORLD,requests(1),ierr)
      call MPI_IRECV(p(lo(1)-nh,lo(2)-nh,lo(3)-nh  ),1,halo,nb(0),1, &
                      MPI_COMM_WORLD,requests(2),ierr)
      call MPI_ISEND(p(lo(1)-nh,lo(2)-nh,lo(3)     ),1,halo,nb(0),0, &
                      MPI_COMM_WORLD,requests(3),ierr)
      call MPI_ISEND(p(lo(1)-nh,lo(2)-nh,hi(3)-nh+1),1,halo,nb(1),1, &
                      MPI_COMM_WORLD,requests(4),ierr)
      call MPI_WAITALL(4,requests,MPI_STATUSES_IGNORE,ierr)
#endif
    end select
  end subroutine updthalo
#if defined(_OPENACC)
  subroutine updthalo_gpu(nh,periods,p)
    use mod_precision, only: rp,sp,dp,i8,MPI_REAL_RP
    use cudecomp
    use mod_common_cudecomp, only: work => work_halo, &
                                   ch => handle,gd => gd_halo, &
                                   dtype => cudecomp_real_rp, &
                                   istream => istream_acc_queue_1
    implicit none
    integer , intent(in) :: nh
    logical , intent(in) :: periods(3)
    real(rp), intent(inout), dimension(1-nh:,1-nh:,1-nh:) :: p
    integer :: istat
    !$acc host_data use_device(p,work)
    select case(ipencil_axis)
    case(1)
      istat = cudecompUpdateHalosX(ch,gd,p,work,dtype,[nh,nh,nh],periods,2,stream=istream)
      istat = cudecompUpdateHalosX(ch,gd,p,work,dtype,[nh,nh,nh],periods,3,stream=istream)
    case(2)
      istat = cudecompUpdateHalosY(ch,gd,p,work,dtype,[nh,nh,nh],periods,1,stream=istream)
      istat = cudecompUpdateHalosY(ch,gd,p,work,dtype,[nh,nh,nh],periods,3,stream=istream)
    case(3)
      istat = cudecompUpdateHalosZ(ch,gd,p,work,dtype,[nh,nh,nh],periods,1,stream=istream)
      istat = cudecompUpdateHalosZ(ch,gd,p,work,dtype,[nh,nh,nh],periods,2,stream=istream)
    end select
    !$acc end host_data
  end subroutine updthalo_gpu
#endif
  !
  subroutine initbc(sgstype,cbcvel,bcvel,bcpre,bcsgs,bcu,bcv,bcw,bcp,bcs,bcu_mag,bcv_mag,bcw_mag, &
                    bcuf,bcvf,bcwf,n,is_bound,lwm,l,zc,dl,dzc,h,index_wm)
    !
    ! initialize bcu,bcv,bcw,bcp,bcs,bcu_mag,bcv_mag,bcw_mag,bcuf,bcvf,bcwf
    !
    implicit none
    character(len=*), intent(in) :: sgstype
    character(len=1), intent(inout), dimension(0:1,3,3) :: cbcvel
    real(rp)   , intent(in), dimension(0:1,3,3) :: bcvel
    real(rp)   , intent(in), dimension(0:1,3) :: bcpre,bcsgs
    type(bound), intent(inout) :: bcu,bcv,bcw,bcp,bcs,bcu_mag,bcv_mag,bcw_mag,bcuf,bcvf,bcwf
    integer , intent(in), dimension(3) :: n
    logical , intent(in), dimension(0:1,3) :: is_bound
    integer , intent(in), dimension(0:1,3) :: lwm
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc,dzc
    real(rp), intent(in) :: h
    integer,  intent(out), dimension(0:1,3) :: index_wm
    integer :: i,j,k,i1,i2,j1,j2,k1,k2,ivel,idir
    !
    do idir = 1,3
      do i = 0,1
        if(lwm(i,idir)/=0) then
          do ivel = 1,3
            if(ivel==idir) then
              cbcvel(i,idir,ivel) = 'D'
            else
              cbcvel(i,idir,ivel) = 'N'
            end if
          end do
        end if
      end do
    end do
    !
    ! unnecessary to consider wall model here, because
    ! bcu,bcv,bcw are overwritten by updt_wallmodelbc before set_bc, and that
    ! imp3d never used for WMLES. Only imp1d and explicit can be applied.
    !
    bcu%x(:,:,0) = bcvel(0,1,1)
    bcv%x(:,:,0) = bcvel(0,1,2)
    bcw%x(:,:,0) = bcvel(0,1,3)
    bcu%x(:,:,1) = bcvel(1,1,1)
    bcv%x(:,:,1) = bcvel(1,1,2)
    bcw%x(:,:,1) = bcvel(1,1,3)
    bcu%y(:,:,0) = bcvel(0,2,1)
    bcv%y(:,:,0) = bcvel(0,2,2)
    bcw%y(:,:,0) = bcvel(0,2,3)
    bcu%y(:,:,1) = bcvel(1,2,1)
    bcv%y(:,:,1) = bcvel(1,2,2)
    bcw%y(:,:,1) = bcvel(1,2,3)
    bcu%z(:,:,0) = bcvel(0,3,1)
    bcv%z(:,:,0) = bcvel(0,3,2)
    bcw%z(:,:,0) = bcvel(0,3,3)
    bcu%z(:,:,1) = bcvel(1,3,1)
    bcv%z(:,:,1) = bcvel(1,3,2)
    bcw%z(:,:,1) = bcvel(1,3,3)
    !
    bcp%x(:,:,0) = bcpre(0,1)
    bcp%x(:,:,1) = bcpre(1,1)
    bcp%y(:,:,0) = bcpre(0,2)
    bcp%y(:,:,1) = bcpre(1,2)
    bcp%z(:,:,0) = bcpre(0,3)
    bcp%z(:,:,1) = bcpre(1,3)
    !
    bcs%x(:,:,0) = bcsgs(0,1)
    bcs%x(:,:,1) = bcsgs(1,1)
    bcs%y(:,:,0) = bcsgs(0,2)
    bcs%y(:,:,1) = bcsgs(1,2)
    bcs%z(:,:,0) = bcsgs(0,3)
    bcs%z(:,:,1) = bcsgs(1,3)
    !
    ! magnitude (+/-) of velocity used for wall model
    !
    bcu_mag = bcu
    bcv_mag = bcv
    bcw_mag = bcw
    !
    bcuf = bcu
    bcvf = bcv
    bcwf = bcw
    !
    ! find the index required for interpolation to the wall model height.
    ! The stored index corresponds to the cells far from a wall, i.e., i2,j2,k2.
    ! Remmeber to set h strightly higher than the first cell center, and lower 
    ! than the last cell center (h=h-eps)
    !
    if(is_bound(0,1).and.lwm(0,1)/=0) then
      i = 1
      do while((i-0.5)*dl(1) < h)
        i = i + 1
      end do
      i2 = i
      i1 = i - 1
      index_wm(0,1) = i2
    end if
    if(is_bound(1,1).and.lwm(1,1)/=0) then
      i = n(1)
      do while((n(1)-i+0.5)*dl(1) < h)
        i = i - 1
      end do
      i2 = i
      i1 = i + 1
      index_wm(1,1) = i2
    end if
    if(is_bound(0,2).and.lwm(0,2)/=0) then
      j = 1
      do while((j-0.5)*dl(2) < h)
        j = j + 1
      end do
      j2 = j
      j1 = j - 1
      index_wm(0,2) = j2
    end if
    if(is_bound(1,2).and.lwm(1,2)/=0) then
      j = n(2)
      do while((n(2)-j+0.5)*dl(2) < h)
        j = j - 1
      end do
      j2 = j
      j1 = j + 1
      index_wm(1,2) = j2
    end if
    if(is_bound(0,3).and.lwm(0,3)/=0) then
      k = 1
      do while(zc(k) < h)
        k = k + 1
      end do
      k2 = k
      k1 = k - 1
      index_wm(0,3) = k2
    end if
    if(is_bound(1,3).and.lwm(1,3)/=0) then
      k = n(3)
      do while(l(3)-zc(k) < h)
        k = k - 1
      end do
      k2 = k
      k1 = k + 1
      index_wm(1,3) = k2
    end if
    !
  end subroutine initbc
end module mod_bound
