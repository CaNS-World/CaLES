! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_initgrid
  use mod_param, only:pi
  use mod_precision, only: rp,sp,dp,i8,MPI_REAL_RP
  implicit none
  private
  public initgrid
  contains
  subroutine initgrid(gtype,n,gr,lz,dzc,dzf,zc,zf)
    !
    ! initializes the non-uniform grid along z
    !
    implicit none
    integer, parameter :: CLUSTER_TWO_END              = 1, &
                          CLUSTER_ONE_END              = 2, &
                          CLUSTER_ONE_END_R            = 3, &
                          CLUSTER_MIDDLE               = 4, &
                          CLUSTER_NATURAL              = 5, &
                          CLUSTER_WALL_MODEL           = 6
    integer , intent(in ) :: gtype,n
    real(rp), intent(in ) :: gr,lz
    real(rp), intent(out), dimension(0:n+1) :: dzc,dzf,zc,zf
    real(rp) :: z0
    integer :: k
    procedure (), pointer :: gridpoint => null()
    select case(gtype)
    case(CLUSTER_TWO_END)
      gridpoint => gridpoint_cluster_two_end
    case(CLUSTER_ONE_END)
      gridpoint => gridpoint_cluster_one_end
    case(CLUSTER_ONE_END_R)
      gridpoint => gridpoint_cluster_one_end_r
    case(CLUSTER_MIDDLE)
      gridpoint => gridpoint_cluster_middle
    case(CLUSTER_NATURAL)
      gridpoint => gridpoint_cluster_natural
    case(CLUSTER_WALL_MODEL)
      gridpoint => gridpoint_cluster_wall_model
    case default
      gridpoint => gridpoint_cluster_two_end
    end select
    !
    ! step 1) determine coordinates of cell faces zf
    !
    zf(0) = 0.
    do k=1,n
      z0  = (k-0.)/(1.*n)
      call gridpoint(k,n,gr,z0,zf(k))
      zf(k) = zf(k)*lz
    end do
    !
    ! step 2) determine grid spacing between faces dzf
    !
    do k=1,n
      dzf(k) = zf(k)-zf(k-1)
    end do
    dzf(0  ) = dzf(1) !halo cell same as the first cell
    dzf(n+1) = dzf(n) !halo cell same as the last cell
    !
    ! step 3) determine grid spacing between centers dzc
    !
    do k=0,n
      dzc(k) = .5*(dzf(k)+dzf(k+1))
    end do
    dzc(n+1) = dzc(n)
    !
    ! step 4) compute coordinates of cell centers zc and faces zf
    !
    zc(0)    = -dzc(0)/2.
    zf(0)    = 0.
    do k=1,n+1
      zc(k) = zc(k-1) + dzc(k-1)
      zf(k) = zf(k-1) + dzf(k)
    end do
  end subroutine initgrid
  !
  ! grid stretching functions
  ! see e.g., Fluid Flow Phenomena -- A Numerical Toolkit, by P. Orlandi
  !           Pirozzoli et al. JFM 788, 614–639 (commented)
  !
  subroutine gridpoint_cluster_two_end(kg,nzg,alpha,z0,z)
    !
    ! clustered at the two sides
    !
    implicit none
    integer, intent(in) :: kg,nzg
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    if(alpha /= 0.) then
      z = 0.5*(1.+tanh((z0-0.5)*alpha)/tanh(alpha/2.))
      !z = 0.5*(1.+erf( (z0-0.5)*alpha)/erf( alpha/2.))
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_two_end
  subroutine gridpoint_cluster_one_end(kg,nzg,alpha,z0,z)
    !
    ! clustered at the lower side
    !
    implicit none
    integer, intent(in) :: kg,nzg
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    if(alpha /= 0.) then
      z = 1.0*(1.+tanh((z0-1.0)*alpha)/tanh(alpha/1.))
      !z = 1.0*(1.+erf( (z0-1.0)*alpha)/erf( alpha/1.))
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_one_end
  subroutine gridpoint_cluster_one_end_r(kg,nzg,alpha,r0,r)
    !
    ! clustered at the upper side
    !
    implicit none
    integer, intent(in) :: kg,nzg
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    if(alpha /= 0._rp) then
      r = 1._rp-1.0_rp*(1._rp+tanh((1._rp-r0-1.0_rp)*alpha)/tanh(alpha/1._rp))
      !r = 1._rp-1.0_rp*(1._rp+erf( (1._rp-r0-1.0_rp)*alpha)/erf( alpha/1._rp))
    else
      r = r0
    end if
  end subroutine gridpoint_cluster_one_end_r
  subroutine gridpoint_cluster_middle(kg,nzg,alpha,z0,z)
    !
    ! clustered in the middle
    !
    implicit none
    integer, intent(in) :: kg,nzg
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    if(alpha /= 0.) then
      if(     z0 <= 0.5) then
        z = 0.5*(1.-1.+tanh(2.*alpha*(z0-0.))/tanh(alpha))
        !z = 0.5*(1.-1.+erf( 2.*alpha*(z0-0.))/erf( alpha))
      else if(z0  > 0.5) then
        z = 0.5*(1.+1.+tanh(2.*alpha*(z0-1.))/tanh(alpha))
        !z = 0.5*(1.+1.+erf( 2.*alpha*(z0-1.))/erf( alpha))
      end if
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_middle
  subroutine gridpoint_cluster_wall_model(kg,nzg,alpha,z0,z)
    !
    ! clustered using Larsson's formula
    !
    implicit none
    integer, intent(in) :: kg,nzg
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    real(rp) :: dzc
    !
    dzc = 0.1*32./nzg ! dzc = 0.1 for nzg = 32
    z = z0 - (dzc*nzg/2.-1.)/(2.*pi)*sin(2.*pi*z0)
  end subroutine gridpoint_cluster_wall_model
  subroutine gridpoint_cluster_natural(kg,nzg,dummy,z0,z)
    !
    ! a physics-based, 'natural' grid stretching function for wall-bounded turbulence
    ! see Pirozzoli & Orlandi, JCP 439 - 110408 (2021)
    !
    ! clustered at the two sides
    !
    implicit none
    real(rp), parameter :: kb_p     = 32._rp,    &
                           alpha_p  = pi/1.5_rp, &
                           c_eta_p  = 0.8_rp,    &
                           dyp_p    = 0.05_rp
    integer , intent(in ) :: kg,nzg
    real(rp), intent(in)  :: dummy,z0
    real(rp), intent(out) :: z
    real(rp) :: kb,alpha,c_eta,dyp
    real(rp) :: retau,n,k
    !
    ! handle input parameters
    kb    = kb_p
    alpha = alpha_p
    c_eta = c_eta_p
    dyp   = dyp_p
    ! determine retau
    n = nzg/2._rp
    retau = 1._rp/(1._rp+(n/kb)**2)*(dyp*n+(3._rp/4._rp*alpha*c_eta*n)**(4._rp/3._rp)*(n/kb)**2)
    if(kg==1) print*,'Grid targeting Retau = ',retau
    k = 1._rp*min(kg,(nzg-kg))
    ! dermine z/(2h)
    z = 1._rp/(1._rp+(k/kb)**2)*(dyp*k+(3._rp/4._rp*alpha*c_eta*k)**(4._rp/3._rp)*(k/kb)**2)/(2._rp*retau)
    if( kg > nzg-kg ) z = 1._rp-z
  end subroutine gridpoint_cluster_natural
end module mod_initgrid
