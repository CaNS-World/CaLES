! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_mom
  use mpi
  use mod_common_mpi, only: ierr
  use mod_precision, only: rp,sp,dp,i8,MPI_REAL_RP
  implicit none
  private
  public cmpt_wallshear,bulk_forcing,mom_xyz_ad
  contains
  !
  subroutine mom_xyz_ad(nx,ny,nz,dxi,dyi,dzci,dzfi,visc,u,v,w,visct,dudt,dvdt,dwdt,dudtd,dvdtd,dwdtd)
    !
    ! lump all r.h.s. of momentum terms (excluding pressure) into a single fast kernel
    ! interpolation of eddy viscosity, verified by visct = x+y+z
    ! extra cross-derivatives, verified by u,v,w = x*y*z
    ! The calculation of the viscous flux for the first off-wall layer of cells uses
    ! u(k=1) and u(k=2). It does not involve the wall, which is the reason that no special
    ! treatment is needed for the viscous flux at the first cells. The incorrect layer
    ! locates only between the wall and u(k=1).
    !
    integer , intent(in   ) :: nx,ny,nz
    real(rp), intent(in   ) :: dxi,dyi
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi
    real(rp), intent(in   ) :: visc
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u,v,w,visct
    real(rp), intent(inout), dimension( :, :, :) :: dudt,dvdt,dwdt
    real(rp), intent(inout), dimension( :, :, :), optional :: dudtd,dvdtd,dwdtd
    integer  :: i,j,k
    real(rp) :: u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp, &
                v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp, &
                w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp, &
                s_ccm,s_pcm,s_cpm,s_cmc,s_pmc,s_mcc,s_ccc,s_pcc,s_mpc,s_cpc,s_cmp,s_mcp,s_ccp, &
                s_ppc,s_pcp,s_cpp
    real(rp) :: uu_ip,uu_im,vu_jp,vu_jm,wu_kp,wu_km, &
                uv_ip,uv_im,vv_jp,vv_jm,wv_kp,wv_km, &
                uw_ip,uw_im,vw_jp,vw_jm,ww_kp,ww_km
    real(rp) :: dudx_ip,dudx_im,dudy_jp,dudy_jm,dudz_kp,dudz_km, &
                dvdx_ip,dvdx_im,dvdy_jp,dvdy_jm,dvdz_kp,dvdz_km, &
                dwdx_ip,dwdx_im,dwdy_jp,dwdy_jm,dwdz_kp,dwdz_km
    real(rp) ::                 dvdx_jp,dvdx_jm,dwdx_kp,dwdx_km, &
                dudy_ip,dudy_im,                dwdy_kp,dwdy_km, &
                dudz_ip,dudz_im,dvdz_jp,dvdz_jm
    real(rp) :: visc_ip,visc_im,visc_jp,visc_jm,visc_kp,visc_km
    real(rp) :: dudt_s ,dvdt_s ,dwdt_s , &
                dudtd_s,dvdtd_s,dwdtd_s, &
                dudtd_xy_s,dudtd_z_s   , &
                dvdtd_xy_s,dvdtd_z_s   , &
                dwdtd_xy_s,dwdtd_z_s
    !
    !$acc parallel loop collapse(3) default(present) async(1) &
    !$acc private(u_ccm,u_pcm,u_cpm,u_cmc,u_pmc,u_mcc,u_ccc,u_pcc,u_mpc,u_cpc,u_cmp,u_mcp,u_ccp) &
    !$acc private(v_ccm,v_pcm,v_cpm,v_cmc,v_pmc,v_mcc,v_ccc,v_pcc,v_mpc,v_cpc,v_cmp,v_mcp,v_ccp) &
    !$acc private(w_ccm,w_pcm,w_cpm,w_cmc,w_pmc,w_mcc,w_ccc,w_pcc,w_mpc,w_cpc,w_cmp,w_mcp,w_ccp) &
    !$acc private(s_ccm,s_pcm,s_cpm,s_cmc,s_pmc,s_mcc,s_ccc,s_pcc,s_mpc,s_cpc,s_cmp,s_mcp,s_ccp) &
    !$acc private(s_ppc,s_pcp,s_cpp) &
    !$acc private(uu_ip,uu_im,vu_jp,vu_jm,wu_kp,wu_km) &
    !$acc private(uv_ip,uv_im,vv_jp,vv_jm,wv_kp,wv_km) &
    !$acc private(uw_ip,uw_im,vw_jp,vw_jm,ww_kp,ww_km) &
    !$acc private(dudx_ip,dudx_im,dudy_jp,dudy_jm,dudz_kp,dudz_km) &
    !$acc private(dvdx_ip,dvdx_im,dvdy_jp,dvdy_jm,dvdz_kp,dvdz_km) &
    !$acc private(dwdx_ip,dwdx_im,dwdy_jp,dwdy_jm,dwdz_kp,dwdz_km) &
    !$acc private(                dvdx_jp,dvdx_jm,dwdx_kp,dwdx_km) &
    !$acc private(dudy_ip,dudy_im,                dwdy_kp,dwdy_km) &
    !$acc private(dudz_ip,dudz_im,dvdz_jp,dvdz_jm                ) &
    !$acc private(visc_ip,visc_im,visc_jp,visc_jm,visc_kp,visc_km) &
    !$acc private(dudt_s ,dvdt_s ,dwdt_s ) &
    !$acc private(dudtd_s,dvdtd_s,dwdtd_s) &
    !$acc private(dudtd_xy_s,dudtd_z_s) &
    !$acc private(dvdtd_xy_s,dvdtd_z_s) &
    !$acc private(dwdtd_xy_s,dwdtd_z_s)
    do k=1,nz
      do j=1,ny
        do i=1,nx
          !
          ! touch u,v,w sequentially
          !
          u_ccm = u(i  ,j  ,k-1)
          u_pcm = u(i+1,j  ,k-1)
          u_cpm = u(i  ,j+1,k-1)
          u_cmc = u(i  ,j-1,k  )
          u_pmc = u(i+1,j-1,k  )
          u_mcc = u(i-1,j  ,k  )
          u_ccc = u(i  ,j  ,k  )
          u_pcc = u(i+1,j  ,k  )
          u_mpc = u(i-1,j+1,k  )
          u_cpc = u(i  ,j+1,k  )
          u_cmp = u(i  ,j-1,k+1)
          u_mcp = u(i-1,j  ,k+1)
          u_ccp = u(i  ,j  ,k+1)
          !
          v_ccm = v(i  ,j  ,k-1)
          v_pcm = v(i+1,j  ,k-1)
          v_cpm = v(i  ,j+1,k-1)
          v_cmc = v(i  ,j-1,k  )
          v_pmc = v(i+1,j-1,k  )
          v_mcc = v(i-1,j  ,k  )
          v_ccc = v(i  ,j  ,k  )
          v_pcc = v(i+1,j  ,k  )
          v_mpc = v(i-1,j+1,k  )
          v_cpc = v(i  ,j+1,k  )
          v_cmp = v(i  ,j-1,k+1)
          v_mcp = v(i-1,j  ,k+1)
          v_ccp = v(i  ,j  ,k+1)
          !
          w_ccm = w(i  ,j  ,k-1)
          w_pcm = w(i+1,j  ,k-1)
          w_cpm = w(i  ,j+1,k-1)
          w_cmc = w(i  ,j-1,k  )
          w_pmc = w(i+1,j-1,k  )
          w_mcc = w(i-1,j  ,k  )
          w_ccc = w(i  ,j  ,k  )
          w_pcc = w(i+1,j  ,k  )
          w_mpc = w(i-1,j+1,k  )
          w_cpc = w(i  ,j+1,k  )
          w_cmp = w(i  ,j-1,k+1)
          w_mcp = w(i-1,j  ,k+1)
          w_ccp = w(i  ,j  ,k+1)
          !
          s_ccm = visct(i  ,j  ,k-1)
          s_pcm = visct(i+1,j  ,k-1)
          s_cpm = visct(i  ,j+1,k-1)
          s_cmc = visct(i  ,j-1,k  )
          s_pmc = visct(i+1,j-1,k  )
          s_mcc = visct(i-1,j  ,k  )
          s_ccc = visct(i  ,j  ,k  )
          s_pcc = visct(i+1,j  ,k  )
          s_mpc = visct(i-1,j+1,k  )
          s_cpc = visct(i  ,j+1,k  )
          s_cmp = visct(i  ,j-1,k+1)
          s_mcp = visct(i-1,j  ,k+1)
          s_ccp = visct(i  ,j  ,k+1)
          !
          s_ppc = visct(i+1,j+1,k  )
          s_pcp = visct(i+1,j  ,k+1)
          s_cpp = visct(i  ,j+1,k+1)
          !
          ! x diffusion
          !
          visc_ip = s_pcc
          visc_im = s_ccc
          visc_jp = 0.25_rp*(s_ccc+s_pcc+s_cpc+s_ppc)
          visc_jm = 0.25_rp*(s_ccc+s_pcc+s_cmc+s_pmc)
          visc_kp = 0.25_rp*(s_ccc+s_pcc+s_ccp+s_pcp)
          visc_km = 0.25_rp*(s_ccc+s_pcc+s_ccm+s_pcm)
          !
          dudx_ip = (u_pcc-u_ccc)*dxi
          dudx_im = (u_ccc-u_mcc)*dxi
          dudy_jp = (u_cpc-u_ccc)*dyi
          dudy_jm = (u_ccc-u_cmc)*dyi
          dudz_kp = (u_ccp-u_ccc)*dzci(k  )
          dudz_km = (u_ccc-u_ccm)*dzci(k-1)
          !
          dudx_ip = dudx_ip
          dudx_im = dudx_im
          dvdx_jp = (v_pcc-v_ccc)*dxi
          dvdx_jm = (v_pmc-v_cmc)*dxi
          dwdx_kp = (w_pcc-w_ccc)*dxi
          dwdx_km = (w_pcm-w_ccm)*dxi
          !
          ! x advection
          !
          uu_ip  = 0.25_rp*(u_pcc+u_ccc)*(u_ccc+u_pcc)
          uu_im  = 0.25_rp*(u_mcc+u_ccc)*(u_ccc+u_mcc)
          vu_jp  = 0.25_rp*(v_pcc+v_ccc)*(u_ccc+u_cpc)
          vu_jm  = 0.25_rp*(v_pmc+v_cmc)*(u_ccc+u_cmc)
          wu_kp  = 0.25_rp*(w_pcc+w_ccc)*(u_ccc+u_ccp)
          wu_km  = 0.25_rp*(w_pcm+w_ccm)*(u_ccc+u_ccm)
          !
          dudtd_xy_s = &
                        visc*(dudx_ip-dudx_im)*dxi + & ! d(dudx)/dx
                        visc*(dudy_jp-dudy_jm)*dyi     ! d(dudy)/dy
          dudtd_z_s  = &
                        visc*(dudz_kp-dudz_km)*dzfi(k) ! d(dudz)/dz
          dudt_s     = &
                        -(uu_ip-uu_im)*dxi - &
                         (vu_jp-vu_jm)*dyi - &
                         (wu_kp-wu_km)*dzfi(k) &
                        +(visc_ip*(dudx_ip+dudx_ip)-visc_im*(dudx_im+dudx_im))*dxi + & ! d(dudx+dudx)/dx
                         (visc_jp*(dudy_jp+dvdx_jp)-visc_jm*(dudy_jm+dvdx_jm))*dyi + & ! d(dudy+dvdx)/dy
                         (visc_kp*(dudz_kp+dwdx_kp)-visc_km*(dudz_km+dwdx_km))*dzfi(k) ! d(dudz+dwdx)/dz
          !
          ! y diffusion
          !
          visc_ip = 0.25_rp*(s_ccc+s_cpc+s_pcc+s_ppc)
          visc_im = 0.25_rp*(s_ccc+s_cpc+s_mcc+s_mpc)
          visc_jp = s_cpc
          visc_jm = s_ccc
          visc_kp = 0.25_rp*(s_ccc+s_cpc+s_ccp+s_cpp)
          visc_km = 0.25_rp*(s_ccc+s_cpc+s_ccm+s_cpm)
          !
          dvdx_ip = (v_pcc-v_ccc)*dxi
          dvdx_im = (v_ccc-v_mcc)*dxi
          dvdy_jp = (v_cpc-v_ccc)*dyi
          dvdy_jm = (v_ccc-v_cmc)*dyi
          dvdz_kp = (v_ccp-v_ccc)*dzci(k  )
          dvdz_km = (v_ccc-v_ccm)*dzci(k-1)
          !
          dudy_ip = (u_cpc-u_ccc)*dyi
          dudy_im = (u_mpc-u_mcc)*dyi
          dvdy_jp = dvdy_jp
          dvdy_jm = dvdy_jm
          dwdy_kp = (w_cpc-w_ccc)*dyi
          dwdy_km = (w_cpm-w_ccm)*dyi
          !
          ! y advection
          !
          uv_ip  = 0.25_rp*(u_ccc+u_cpc)*(v_ccc+v_pcc)
          uv_im  = 0.25_rp*(u_mcc+u_mpc)*(v_ccc+v_mcc)
          vv_jp  = 0.25_rp*(v_ccc+v_cpc)*(v_ccc+v_cpc)
          vv_jm  = 0.25_rp*(v_ccc+v_cmc)*(v_ccc+v_cmc)
          wv_kp  = 0.25_rp*(w_ccc+w_cpc)*(v_ccc+v_ccp)
          wv_km  = 0.25_rp*(w_ccm+w_cpm)*(v_ccc+v_ccm)
          !
          dvdtd_xy_s = &
                        visc*(dvdx_ip-dvdx_im)*dxi + & ! d(dvdx)/dx
                        visc*(dvdy_jp-dvdy_jm)*dyi     ! d(dvdy)/dy
          dvdtd_z_s  = &
                        visc*(dvdz_kp-dvdz_km)*dzfi(k) ! d(dvdz)/dz
          dvdt_s     = &
                        -(uv_ip-uv_im)*dxi - &
                         (vv_jp-vv_jm)*dyi - &
                         (wv_kp-wv_km)*dzfi(k) &
                        +(visc_ip*(dvdx_ip+dudy_ip)-visc_im*(dvdx_im+dudy_im))*dxi + & ! d(dvdx+dudy)/dx
                         (visc_jp*(dvdy_jp+dvdy_jp)-visc_jm*(dvdy_jm+dvdy_jm))*dyi + & ! d(dvdy+dvdy)/dy
                         (visc_kp*(dvdz_kp+dwdy_kp)-visc_km*(dvdz_km+dwdy_km))*dzfi(k) ! d(dvdz+dwdy)/dz
          !
          ! z diffusion
          !
          visc_ip = 0.25_rp*(s_ccc+s_ccp+s_pcc+s_pcp)
          visc_im = 0.25_rp*(s_ccc+s_ccp+s_mcc+s_mcp)
          visc_jp = 0.25_rp*(s_ccc+s_ccp+s_cpc+s_cpp)
          visc_jm = 0.25_rp*(s_ccc+s_ccp+s_cmc+s_cmp)
          visc_kp = s_ccp
          visc_km = s_ccc
          !
          dwdx_ip = (w_pcc-w_ccc)*dxi
          dwdx_im = (w_ccc-w_mcc)*dxi
          dwdy_jp = (w_cpc-w_ccc)*dyi
          dwdy_jm = (w_ccc-w_cmc)*dyi
          dwdz_kp = (w_ccp-w_ccc)*dzfi(k+1)
          dwdz_km = (w_ccc-w_ccm)*dzfi(k  )
          !
          dudz_ip = (u_ccp-u_ccc)*dzci(k  )
          dudz_im = (u_mcp-u_mcc)*dzci(k  )
          dvdz_jp = (v_ccp-v_ccc)*dzci(k  )
          dvdz_jm = (v_cmp-v_cmc)*dzci(k  )
          dwdz_kp = dwdz_kp
          dwdz_km = dwdz_km
          !
          ! z advection
          !
          uw_ip  = 0.25_rp*(u_ccc+u_ccp)*(w_ccc+w_pcc)
          uw_im  = 0.25_rp*(u_mcc+u_mcp)*(w_ccc+w_mcc)
          vw_jp  = 0.25_rp*(v_ccc+v_ccp)*(w_ccc+w_cpc)
          vw_jm  = 0.25_rp*(v_cmc+v_cmp)*(w_ccc+w_cmc)
          ww_kp  = 0.25_rp*(w_ccc+w_ccp)*(w_ccc+w_ccp)
          ww_km  = 0.25_rp*(w_ccc+w_ccm)*(w_ccc+w_ccm)
          !
          dwdtd_xy_s =  &
                          visc*(dwdx_ip-dwdx_im)*dxi + & ! d(dwdx)/dx
                          visc*(dwdy_jp-dwdy_jm)*dyi     ! d(dwdy)/dy
          dwdtd_z_s =   &
                          visc*(dwdz_kp-dwdz_km)*dzci(k) ! d(dwdz)/dz
          dwdt_s     =  &
                         -(uw_ip-uw_im)*dxi - &
                          (vw_jp-vw_jm)*dyi - &
                          (ww_kp-ww_km)*dzci(k) &
                         +(visc_ip*(dwdx_ip+dudz_ip)-visc_im*(dwdx_im+dudz_im))*dxi + & ! d(dwdx+dudz)/dx
                          (visc_jp*(dwdy_jp+dvdz_jp)-visc_jm*(dwdy_jm+dvdz_jm))*dyi + & ! d(dwdy+dvdz)/dy
                          (visc_kp*(dwdz_kp+dwdz_kp)-visc_km*(dwdz_km+dwdz_km))*dzci(k) ! d(dwdz+dwdz)/dz
#if defined(_IMPDIFF)
#if defined(_IMPDIFF_1D)
          dudt_s = dudt_s + dudtd_xy_s
          dvdt_s = dvdt_s + dvdtd_xy_s
          dwdt_s = dwdt_s + dwdtd_xy_s
          dudtd_s = dudtd_z_s
          dvdtd_s = dvdtd_z_s
          dwdtd_s = dwdtd_z_s
#else
          dudtd_s = dudtd_xy_s + dudtd_z_s
          dvdtd_s = dvdtd_xy_s + dvdtd_z_s
          dwdtd_s = dwdtd_xy_s + dwdtd_z_s
#endif
          dudt( i,j,k) = dudt_s
          dvdt( i,j,k) = dvdt_s
          dwdt( i,j,k) = dwdt_s
          dudtd(i,j,k) = dudtd_s
          dvdtd(i,j,k) = dvdtd_s
          dwdtd(i,j,k) = dwdtd_s
#else
          dudt_s = dudt_s + dudtd_xy_s + dudtd_z_s
          dvdt_s = dvdt_s + dvdtd_xy_s + dvdtd_z_s
          dwdt_s = dwdt_s + dwdtd_xy_s + dwdtd_z_s
          dudt(i,j,k) = dudt_s
          dvdt(i,j,k) = dvdt_s
          dwdt(i,j,k) = dwdt_s
#endif
        end do
      end do
    end do
    ! dudt,  explicit terms
    ! dudtd, implicit terms
  end subroutine mom_xyz_ad
  !
  subroutine bulk_forcing(n,is_forced,f,u,v,w)
    integer , intent(in   ), dimension(3) :: n
    logical , intent(in   ), dimension(3) :: is_forced
    real(rp), intent(in   ), dimension(3) :: f
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp) :: ff
    if(is_forced(1)) then
      ff = f(1)
      !$acc kernels default(present) async(1)
      u(1:n(1),1:n(2),1:n(3)) = u(1:n(1),1:n(2),1:n(3)) + ff
      !$acc end kernels
    end if
    if(is_forced(2)) then
      ff = f(2)
      !$acc kernels default(present) async(1)
      v(1:n(1),1:n(2),1:n(3)) = v(1:n(1),1:n(2),1:n(3)) + ff
      !$acc end kernels
    end if
    if(is_forced(3)) then
      ff = f(3)
      !$acc kernels default(present) async(1)
      w(1:n(1),1:n(2),1:n(3)) = w(1:n(1),1:n(2),1:n(3)) + ff
      !$acc end kernels
    end if
  end subroutine bulk_forcing
  !
  subroutine cmpt_wallshear(n,is_cmpt,is_bound,l,dli,dzci,dzfi,visc,u,v,w,taux,tauy,tauz)
    use mod_params, only: cbcpre
    implicit none
    integer , intent(in ), dimension(3) :: n
    logical , intent(in ), dimension(    3) :: is_cmpt
    logical , intent(in ), dimension(0:1,3) :: is_bound
    real(rp), intent(in ), dimension(3)     :: l,dli
    real(rp), intent(in ), dimension(0:)    :: dzci,dzfi
    real(rp), intent(in )                   :: visc
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out), dimension(3) :: taux,tauy,tauz
    real(rp) :: dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dvdzp,dvdzm, &
                dwdxp,dwdxm,dwdyp,dwdym
    !
    ! n.b.: replace scalars with reduction of tau(1:3,1:3) once the
    !       nvfortran bug for array reductions on Pascal architectures
    !       is solved; this subroutine is not used in production anyway
    !
    real(rp) :: tau21,tau31,tau12,tau32,tau13,tau23
    integer :: i,j,k,nx,ny,nz
    real(rp) :: dxi,dyi,lx,ly,lz
    real(rp) :: tau(3,3)
    !
    nx = n(1); ny = n(2); nz = n(3)
    dxi = dli(1); dyi = dli(2)
    lx = l(1); ly = l(2); lz = l(3)
    tau21 = 0._rp
    tau31 = 0._rp
    if(is_cmpt(1)) then
      !$acc data copy(tau21) async(1)
      if(is_bound(0,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
        !$acc parallel loop collapse(2) default(present) private(dudyp) reduction(+:tau21) async(1)
        do k=1,nz
          do i=1,nx
            dudyp = (u(i,1 ,k)-u(i,0   ,k))*dyi*visc
            tau21 = tau21 + dudyp/(dxi*dzfi(k)*lx*lz)
          end do
        end do
      end if
      if(is_bound(1,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
        !$acc parallel loop collapse(2) default(present) private(dudym) reduction(+:tau21) async(1)
        do k=1,nz
          do i=1,nx
            dudym = (u(i,ny,k)-u(i,ny+1,k))*dyi*visc
            tau21 = tau21 + dudym/(dxi*dzfi(k)*lx*lz)
          end do
        end do
      end if
      !$acc end data
      !$acc data copy(tau31) async(1)
      if(is_bound(0,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
        !$acc parallel loop collapse(2) default(present) private(dudzp) reduction(+:tau31) async(1)
        do j=1,ny
          do i=1,nx
            dudzp = (u(i,j,1 )-u(i,j,0   ))*dzci(0)*visc
            tau31 = tau31 + dudzp/(dxi*dyi*lx*ly)
          end do
        end do
      end if
      if(is_bound(1,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
        !$acc parallel loop collapse(2) default(present) private(dudzm) reduction(+:tau31) async(1)
        do j=1,ny
          do i=1,nx
            dudzm = (u(i,j,nz)-u(i,j,nz+1))*dzci(nz)*visc
            tau31 = tau31 + dudzm/(dxi*dyi*lx*ly)
          end do
        end do
      end if
      !$acc end data
    end if
    !
    tau12 = 0._rp
    tau32 = 0._rp
    if(is_cmpt(2)) then
      !$acc data copy(tau12) async(1)
      if(is_bound(0,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
        !$acc parallel loop collapse(2) default(present) private(dvdxp) reduction(+:tau12) async(1)
        do k=1,nz
          do j=1,ny
            dvdxp = (v(1  ,j,k)-v(0  ,j,k))*dxi*visc
            tau12 = tau12 + dvdxp/(dyi*dzfi(k)*ly*lz)
          end do
        end do
      end if
      if(is_bound(1,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
        !$acc parallel loop collapse(2) default(present) private(dvdxm) reduction(+:tau12) async(1)
        do k=1,nz
          do j=1,ny
            dvdxm = (v(nx,j,k)-v(nx+1,j,k))*dxi*visc
            tau12 = tau12 + dvdxm/(dyi*dzfi(k)*ly*lz)
          end do
        end do
      end if
      !$acc end data
      !$acc data copy(tau32) async(1)
      if(is_bound(0,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
        !$acc parallel loop collapse(2) default(present) private(dvdzp) reduction(+:tau32) async(1)
        do j=1,ny
          do i=1,nx
            dvdzp = (v(i,j,1 )-v(i,j,0   ))*dzci(0)*visc
            tau32 = tau32 + dvdzp/(dxi*dyi*lx*ly)
          end do
        end do
      end if
      if(is_bound(1,3).and.cbcpre(0,3)//cbcpre(1,3) /= 'PP') then
        !$acc parallel loop collapse(2) default(present) private(dvdzm) reduction(+:tau32) async(1)
        do j=1,ny
          do i=1,nx
            dvdzm = (v(i,j,nz)-v(i,j,nz+1))*dzci(nz)*visc
            tau32 = tau32 + dvdzm/(dxi*dyi*lx*ly)
          end do
        end do
      end if
      !$acc end data
    end if
    !
    tau13 = 0._rp
    tau23 = 0._rp
    if(is_cmpt(3)) then
      !$acc data copy(tau13) async(1)
      if(is_bound(0,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
        !$acc parallel loop collapse(2) default(present) private(dwdxp) reduction(+:tau13) async(1)
        do k=1,nz
          do j=1,ny
            dwdxp = (w(1 ,j,k)-w(0   ,j,k))*dxi*visc
            tau13 = tau13 + dwdxp/(dyi*dzfi(k)*ly*lz)
          end do
        end do
      end if
      if(is_bound(1,1).and.cbcpre(0,1)//cbcpre(1,1) /= 'PP') then
        !$acc parallel loop collapse(2) default(present) private(dwdxm) reduction(+:tau13) async(1)
        do k=1,nz
          do j=1,ny
            dwdxm = (w(nx,j,k)-w(nx+1,j,k))*dxi*visc
            tau13 = tau13 + dwdxm/(dyi*dzfi(k)*ly*lz)
          end do
        end do
      end if
      !$acc end data
      !$acc data copy(tau23) async(1)
      if(is_bound(0,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
        !$acc parallel loop collapse(2) default(present) private(dwdyp) reduction(+:tau23) async(1)
        do k=1,nz
          do i=1,nx
            dwdyp = (w(i,1,k )-w(i,0   ,k))*dyi*visc
            tau23 = tau23 + dwdyp/(dxi*dzfi(k)*lx*lz)
          end do
        end do
      end if
      if(is_bound(1,2).and.cbcpre(0,2)//cbcpre(1,2) /= 'PP') then
        !$acc parallel loop collapse(2) default(present) private(dwdym) reduction(+:tau23) async(1)
        do k=1,nz
          do i=1,nx
            dwdym = (w(i,ny,k)-w(i,ny+1,k))*dyi*visc
            tau23 = tau23 + dwdym/(dxi*dzfi(k)*lx*lz)
          end do
        end do
      end if
      !$acc end data
    end if
    !$acc wait(1)
    tau(:,:) = 0._rp
    tau(2,1) = tau21
    tau(3,1) = tau31
    tau(1,2) = tau12
    tau(3,2) = tau32
    tau(1,3) = tau13
    tau(2,3) = tau23
    call MPI_ALLREDUCE(MPI_IN_PLACE,tau(1,1),9,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    taux(:) = tau(:,1)
    tauy(:) = tau(:,2)
    tauz(:) = tau(:,3)
  end subroutine cmpt_wallshear
end module mod_mom
