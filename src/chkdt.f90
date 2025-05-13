! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_chkdt
  use mpi
  use mod_common_mpi, only:ierr
  use mod_precision, only: rp,sp,dp,i8,MPI_REAL_RP
  use mod_params, only: eps
  implicit none
  private
  public chkdt
  contains
  subroutine chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
    !
    ! compute maximum allowed time step, refer to Pieter Wesseling (P200)
    ! for the stability conditions of the advective and diffusion terms
    !
    ! the eddy viscosity term is taken into account in the calculation of dt.
    ! It is acceptable not to consider it, since it is larger than both
    ! the viscous and advective terms. In WRLES, the viscous term is commonly
    ! implicitly treated, so the eddy viscosity term does not influence dtmax,
    ! even when the grid is very fine
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzci,dzfi
    real(rp), intent(in) :: visc
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(out) :: dtmax
    real(rp) :: dxi,dyi,dzi
    real(rp) :: uc,vc,wc
    real(rp) :: this_dti,this_dtid,dti,dtid
    integer :: i,j,k
    !
    dxi = 1._rp/dl(1)
    dyi = 1._rp/dl(2)
    dzi = 1._rp/dl(3)
    !
    dti  = 0._rp
    dtid = 0._rp
    !$acc data copy(dti,dtid) async(1)
    !$acc parallel loop collapse(3) default(present) reduction(max:dti,dtid) async(1) &
    !$acc private(uc,vc,wc,this_dti,this_dtid)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          uc = 0.5_rp*abs(u(i-1,j,k)+u(i,j,k))
          vc = 0.5_rp*abs(v(i,j-1,k)+v(i,j,k))
          wc = 0.5_rp*abs(w(i,j,k-1)+w(i,j,k))
          this_dti = uc*dxi+vc*dyi+wc*dzfi(k)
          dti = max(dti,this_dti)
          !
#if defined(_IMPDIFF) && !defined(_IMPDIFF_1D)
          this_dtid = 0._rp
#else
          this_dtid = visc*(dxi*dxi+dyi*dyi)
#if !defined(_IMPDIFF_1D)
          this_dtid = this_dtid + visc*(dzfi(k)*dzfi(k))
#endif
#endif
          dtid  = max(dtid,this_dtid)
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    if(dti  <= eps) dti  = 1._rp
    if(dtid <= eps) dtid = eps
    dtmax = min(0.4125_rp/dtid,1.732_rp/dti) ! viscous CFL could be 1.5
    call MPI_ALLREDUCE(MPI_IN_PLACE,dtmax,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
  end subroutine chkdt
end module mod_chkdt