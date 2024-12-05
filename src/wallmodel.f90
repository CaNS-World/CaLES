! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_wallmodel
  use, intrinsic :: ieee_arithmetic, only: is_nan => ieee_is_nan, is_finite => ieee_is_finite
  use mpi
  use mod_precision, only: rp, sp, dp, i8, MPI_REAL_RP
  use mod_typedef, only: bound
  use mod_param, only: kap_log, b_log, eps
  use mod_smartredis, only: InitSmartRedis, ExchangeDataSmartRedis
  implicit none
  private
  public :: init_wallmodels, updt_wallmodelbcs

  integer, parameter :: WM_LOG = 1
  integer, parameter :: WM_LAM = 2
  integer, parameter :: WM_DRL = 3
  integer, parameter :: MAX_WALL_MODELS = 10

  type :: WallModelInput
    real(rp), allocatable :: vel1(:,:)
    real(rp), allocatable :: vel2(:,:)
  end type WallModelInput

  type :: WallModelOutput
    real(rp), allocatable :: tauw1(:,:)
    real(rp), allocatable :: tauw2(:,:)
  end type WallModelOutput

  abstract interface
    subroutine wallmodel_interface(visc, hwm, inputs, outputs)
      import :: WallModelInput, WallModelOutput
      import :: rp
      real(rp), intent(in) :: visc, hwm
      type(WallModelInput), intent(in) :: inputs
      type(WallModelOutput), intent(out) :: outputs
    end subroutine wallmodel_interface
  end interface

  type :: WallModelProcedure
    procedure(wallmodel_interface), pointer, nopass :: ptr => null()
  end type WallModelProcedure

  type(WallModelProcedure) :: wallmodel(MAX_WALL_MODELS)

  contains

  subroutine updt_wallmodelbcs(n, is_bound, lwm, l, dl, zc, zf, dzc, dzf, visc, hwm, u, v, w, &
                               bcu, bcv, bcw, bcu_mag, bcv_mag, bcw_mag)
    implicit none
    integer, intent(in), dimension(3) :: n
    logical, intent(in), dimension(0:1,3) :: is_bound
    integer, intent(in), dimension(0:1,3) :: lwm
    real(rp), intent(in), dimension(3) :: l, dl
    real(rp), intent(in), dimension(0:) :: zc, zf, dzc, dzf
    real(rp), intent(in) :: visc, hwm
    real(rp), intent(in), dimension(0:,0:,0:) :: u, v, w
    type(Bound), intent(inout) :: bcu, bcv, bcw
    type(Bound), intent(in) :: bcu_mag, bcv_mag, bcw_mag
    type(WallModelInput) :: inputs
    type(WallModelOutput) :: outputs
    logical, save :: is_first = .true.
    integer, dimension(0:1,3), save :: hwm_index
    integer :: nh, mtype, idir, ibound, cell_index

    if (is_first) then
      is_first = .false.
      call init_wallmodels(n, is_bound, lwm, l, dl, zc, hwm, hwm_index)
    end if

    nh = 1

    do ibound = 0, 1
      do idir = 1, 3
        if (is_bound(ibound, idir) .and. lwm(ibound, idir) /= 0) then
          mtype = lwm(ibound, idir)
          cell_index = hwm_index(ibound, idir)
          call cmpt_wallmodelbc(n, ibound, idir, nh, mtype, l, dl, zc, zf, dzc, dzf, &
                                visc, hwm, cell_index, u, v, w, bcu, bcv, bcw, bcu_mag, bcv_mag, bcw_mag)
        end if
      end do
    end do
  end subroutine updt_wallmodelbcs

  ! find the cell_index required for interpolation to the wall model height.
  ! The stored cell_index corresponds to the cells far from a wall, i.e., i2,j2,k2.
  ! Remmeber to set hwm strightly higher than the first cell center, and lower
  ! than the last cell center (hwm=hwm-eps)

  subroutine init_wallmodels(n, is_bound, lwm, l, dl, zc, hwm, hwm_index)
    implicit none
    integer, intent(in) :: n(3)
    logical, intent(in) :: is_bound(0:1,3)
    integer, intent(in) :: lwm(0:1,3)
    real(rp), intent(in) :: l(3), dl(3), zc(0:), hwm
    integer, intent(out) :: hwm_index(0:1,3)
    integer :: i, j, k, i1, i2, j1, j2, k1, k2

    wallmodel(WM_LOG)%ptr => wallmodel_loglaw
    wallmodel(WM_LAM)%ptr => wallmodel_laminar
    wallmodel(WM_DRL)%ptr => wallmodel_DRL

    if(is_bound(0,1).and.lwm(0,1)/=0) then ! to remove if statement
      i = 1
      do while((i-0.5)*dl(1) < hwm)
        i = i + 1
      end do
      i2 = i
      i1 = i - 1
      hwm_index(0,1) = i2
    end if
    if(is_bound(1,1).and.lwm(1,1)/=0) then
      i = n(1)
      do while((n(1)-i+0.5)*dl(1) < hwm)
        i = i - 1
      end do
      i2 = i
      i1 = i + 1
      hwm_index(1,1) = i2
    end if
    if(is_bound(0,2).and.lwm(0,2)/=0) then
      j = 1
      do while((j-0.5)*dl(2) < hwm)
        j = j + 1
      end do
      j2 = j
      j1 = j - 1
      hwm_index(0,2) = j2
    end if
    if(is_bound(1,2).and.lwm(1,2)/=0) then
      j = n(2)
      do while((n(2)-j+0.5)*dl(2) < hwm)
        j = j - 1
      end do
      j2 = j
      j1 = j + 1
      hwm_index(1,2) = j2
    end if
    if(is_bound(0,3).and.lwm(0,3)/=0) then
      k = 1
      do while(zc(k) < hwm)
        k = k + 1
      end do
      k2 = k
      k1 = k - 1
      hwm_index(0,3) = k2
    end if
    if(is_bound(1,3).and.lwm(1,3)/=0) then
      k = n(3)
      do while(l(3)-zc(k) < hwm)
        k = k - 1
      end do
      k2 = k
      k1 = k + 1
      hwm_index(1,3) = k2
    end if
  end subroutine init_wallmodels

  subroutine cmpt_wallmodelbc(n, ibound, idir, nh, mtype, l, dl, zc, zf, dzc, dzf, visc, hwm, cell_index, &
                              u, v, w, bcu, bcv, bcw, bcu_mag, bcv_mag, bcw_mag)
    implicit none
    integer, intent(in) :: n(3), ibound, idir, nh, mtype
    real(rp), intent(in) :: l(3), dl(3)
    real(rp), intent(in), dimension(0:) :: zc, zf, dzc, dzf
    real(rp), intent(in) :: visc, hwm
    integer, intent(in) :: cell_index
    real(rp), intent(in), dimension(0:,0:,0:) :: u, v, w
    type(Bound), intent(inout) :: bcu, bcv, bcw
    type(Bound), intent(in) :: bcu_mag, bcv_mag, bcw_mag
    type(WallModelInput) :: inputs
    type(WallModelOutput) :: outputs

    call prepare_wallmodel_input(n, idir, ibound, cell_index, 1, l, dl, zc, zf, dzc, dzf, &
                                 visc, hwm, u, v, w, bcu_mag, bcv_mag, bcw_mag, inputs)
    call wallmodel(mtype)%ptr(visc, hwm, inputs, outputs)
    call assign_wallmodelbc(idir, ibound, 1, visc, outputs, bcu, bcv, bcw)

    call prepare_wallmodel_input(n, idir, ibound, cell_index, 2, l, dl, zc, zf, dzc, dzf, &
                                 visc, hwm, u, v, w, bcu_mag, bcv_mag, bcw_mag, inputs)
    call wallmodel(mtype)%ptr(visc, hwm, inputs, outputs)
    call assign_wallmodelbc(idir, ibound, 2, visc, outputs, bcu, bcv, bcw)

  end subroutine cmpt_wallmodelbc

  subroutine prepare_wallmodel_input(n, idir, ibound, cell_index, ivel, l, dl, zc, zf, dzc, dzf, &
                     visc, hwm, u, v, w, bcu_mag, bcv_mag, bcw_mag, inputs)
    implicit none
    integer, intent(in) :: n(3), idir, ibound, cell_index, ivel
    real(rp), intent(in) :: l(3), dl(3)
    real(rp), intent(in), dimension(0:) :: zc, zf, dzc, dzf
    real(rp), intent(in) :: visc, hwm
    real(rp), intent(in), dimension(0:,0:,0:) :: u, v, w
    type(Bound), intent(in) :: bcu_mag, bcv_mag, bcw_mag
    type(WallModelInput), intent(out) :: inputs
    integer :: i, j, k
    integer :: i1, i2, j1, j2, k1, k2
    real(rp) :: coef, wei
    real(rp) :: u1, u2, v1, v2, w1, w2
    real(rp) :: u_mag, v_mag, w_mag

    if (allocated(inputs%vel1)) deallocate(inputs%vel1)
    if (allocated(inputs%vel2)) deallocate(inputs%vel2)

    select case(idir)
    case(1)  ! x-direction walls: vel1 = v, vel2 = w
      allocate(inputs%vel1(0:n(2)+1, 0:n(3)+1))
      allocate(inputs%vel2(0:n(2)+1, 0:n(3)+1))
      inputs%vel1 = 0._rp  ! unused elements have zero value
      inputs%vel2 = 0._rp
      if (ibound == 0) then
        i2 = cell_index
        i1 = cell_index - 1
        coef = (hwm - (i1 - 0.5_rp) * dl(1)) / dl(1)
      else
        i2 = cell_index
        i1 = cell_index + 1
        coef = (hwm - (n(1) - i1 + 0.5_rp) * dl(1)) / dl(1)
      end if
      if (ivel == 1) then
        do k = 1, n(3)
          do j = 0, n(2)
            v1 = v(i1, j, k)
            v2 = v(i2, j, k)
            w1 = 0.25_rp * (w(i1, j, k) + w(i1, j+1, k) + w(i1, j, k-1) + w(i1, j+1, k-1))
            w2 = 0.25_rp * (w(i2, j, k) + w(i2, j+1, k) + w(i2, j, k-1) + w(i2, j+1, k-1))
            v_mag = bcv_mag%x(j, k, ibound)
            w_mag = 0.25_rp * (bcw_mag%x(j, k  , ibound) + bcw_mag%x(j+1, k  , ibound) + &
                               bcw_mag%x(j, k-1, ibound) + bcw_mag%x(j+1, k-1, ibound))
            inputs%vel1(j, k) = vel_relative(v1, v2, coef, v_mag)
            inputs%vel2(j, k) = vel_relative(w1, w2, coef, w_mag)
          end do
        end do
      else
        do k = 0, n(3)
          do j = 1, n(2)
            wei = (zf(k) - zc(k)) / dzc(k)
            v1 = 0.5_rp * ((1._rp - wei) * (v(i1, j-1, k  ) + v(i1, j, k  )) + &
                                    wei  * (v(i1, j-1, k+1) + v(i1, j, k+1)))
            v2 = 0.5_rp * ((1._rp - wei) * (v(i2, j-1, k  ) + v(i2, j, k  )) + &
                                    wei  * (v(i2, j-1, k+1) + v(i2, j, k+1)))
            w1 = w(i1, j, k)
            w2 = w(i2, j, k)
            v_mag = 0.5_rp * ((1._rp - wei) * (bcv_mag%x(j-1, k  , ibound) + bcv_mag%x(j, k  , ibound)) + &
                                       wei  * (bcv_mag%x(j-1, k+1, ibound) + bcv_mag%x(j, k+1, ibound)))
            w_mag = bcw_mag%x(j, k, ibound)
            inputs%vel1(j, k) = vel_relative(v1, v2, coef, v_mag)
            inputs%vel2(j, k) = vel_relative(w1, w2, coef, w_mag)
          end do
        end do
      end if
    case(2)  ! y-direction walls: vel1 = u, vel2 = w
      allocate(inputs%vel1(0:n(1)+1, 0:n(3)+1))
      allocate(inputs%vel2(0:n(1)+1, 0:n(3)+1))
      inputs%vel1 = 0._rp
      inputs%vel2 = 0._rp
      if (ibound == 0) then
        j2 = cell_index
        j1 = cell_index - 1
        coef = (hwm - (j1 - 0.5_rp) * dl(2)) / dl(2)
      else
        j2 = cell_index
        j1 = cell_index + 1
        coef = (hwm - (n(2) - j1 + 0.5_rp) * dl(2)) / dl(2)
      end if
      if (ivel == 1) then
        do k = 1, n(3)
          do i = 0, n(1)
            u1 = u(i, j1, k)
            u2 = u(i, j2, k)
            w1 = 0.25_rp * (w(i, j1, k) + w(i+1, j1, k) + w(i, j1, k-1) + w(i+1, j1, k-1))
            w2 = 0.25_rp * (w(i, j2, k) + w(i+1, j2, k) + w(i, j2, k-1) + w(i+1, j2, k-1))
            u_mag = bcu_mag%y(i, k, ibound)
            w_mag = 0.25_rp * (bcw_mag%y(i, k  , ibound) + bcw_mag%y(i+1, k  , ibound) + &
                               bcw_mag%y(i, k-1, ibound) + bcw_mag%y(i+1, k-1, ibound))
            inputs%vel1(i, k) = vel_relative(u1, u2, coef, u_mag)
            inputs%vel2(i, k) = vel_relative(w1, w2, coef, w_mag)
          end do
        end do
      else
        do k = 0, n(3)
          do i = 1, n(1)
            wei = (zf(k) - zc(k)) / dzc(k)
            u1 = 0.5_rp * ((1._rp - wei) * (u(i-1, j1, k  ) + u(i, j1, k  )) + &
                                    wei  * (u(i-1, j1, k+1) + u(i, j1, k+1)))
            u2 = 0.5_rp * ((1._rp - wei) * (u(i-1, j2, k  ) + u(i, j2, k  )) + &
                                    wei  * (u(i-1, j2, k+1) + u(i, j2, k+1)))
            w1 = w(i, j1, k)
            w2 = w(i, j2, k)
            u_mag = 0.5_rp * ((1._rp - wei) * (bcu_mag%y(i-1, k  , ibound) + bcu_mag%y(i, k  , ibound)) + &
                                       wei  * (bcu_mag%y(i-1, k+1, ibound) + bcu_mag%y(i, k+1, ibound)))
            w_mag = bcw_mag%y(i, k, ibound)
            inputs%vel1(i, k) = vel_relative(u1, u2, coef, u_mag)
            inputs%vel2(i, k) = vel_relative(w1, w2, coef, w_mag)
          end do
        end do
      end if
    case(3)  ! z-direction walls: vel1 = u, vel2 = v
      allocate(inputs%vel1(0:n(1)+1, 0:n(2)+1))
      allocate(inputs%vel2(0:n(1)+1, 0:n(2)+1))
      inputs%vel1 = 0._rp
      inputs%vel2 = 0._rp
      if (ibound == 0) then
        k2 = cell_index
        k1 = cell_index - 1
        coef = (hwm - zc(k1)) / dzc(k1)
      else
        k2 = cell_index
        k1 = cell_index + 1
        coef = (hwm - (l(3) - zc(k1))) / dzc(k2)
      end if
      if (ivel == 1) then
        do j = 1, n(2)
          do i = 0, n(1)
            u1 = u(i, j, k1)
            u2 = u(i, j, k2)
            v1 = 0.25_rp * (v(i, j, k1) + v(i+1, j, k1) + v(i, j-1, k1) + v(i+1, j-1, k1))
            v2 = 0.25_rp * (v(i, j, k2) + v(i+1, j, k2) + v(i, j-1, k2) + v(i+1, j-1, k2))
            u_mag = bcu_mag%z(i, j, ibound)
            v_mag = 0.25_rp * (bcv_mag%z(i, j  , ibound) + bcv_mag%z(i+1, j  , ibound) + &
                               bcv_mag%z(i, j-1, ibound) + bcv_mag%z(i+1, j-1, ibound))
            inputs%vel1(i, j) = vel_relative(u1, u2, coef, u_mag)
            inputs%vel2(i, j) = vel_relative(v1, v2, coef, v_mag)
          end do
        end do
      else
        do j = 0, n(2)
          do i = 1, n(1)
            u1 = 0.25_rp * (u(i-1, j, k1) + u(i, j, k1) + u(i-1, j+1, k1) + u(i, j+1, k1))
            u2 = 0.25_rp * (u(i-1, j, k2) + u(i, j, k2) + u(i-1, j+1, k2) + u(i, j+1, k2))
            v1 = v(i, j, k1)
            v2 = v(i, j, k2)
            u_mag = 0.25_rp * (bcu_mag%z(i-1, j  , ibound) + bcu_mag%z(i, j  , ibound) + &
                               bcu_mag%z(i-1, j+1, ibound) + bcu_mag%z(i, j+1, ibound))
            v_mag = bcv_mag%z(i, j, ibound)
            inputs%vel1(i, j) = vel_relative(u1, u2, coef, u_mag)
            inputs%vel2(i, j) = vel_relative(v1, v2, coef, v_mag)
          end do
        end do
      end if
    end select
    ! Compute gradients if required by certain wall models
    ! Placeholder: Implement gradient computations as needed
    ! For example:
    ! if (some_condition_for_gradients) then
    !   allocate(inputs%grad_u(...))
    !   ... compute gradients ...
    ! end if
  end subroutine prepare_wallmodel_input
  
  subroutine assign_wallmodelbc(idir, ibound, ivel, visc, outputs, bcu, bcv, bcw)
    implicit none
    integer, intent(in) :: idir, ibound, ivel
    real(rp), intent(in) :: visc
    type(WallModelOutput), intent(in) :: outputs
    type(Bound), intent(inout) :: bcu, bcv, bcw
    real(rp) :: visci, sgn
    
    visci = 1._rp/visc

    if (ibound == 0) then
      sgn = 1._rp
    else
      sgn = -1._rp
    end if
    
    select case(idir)
    case(1)
      if (ivel == 1) then
        bcv%x(:,:,ibound) = sgn*visci*outputs%tauw1
      else
        bcw%x(:,:,ibound) = sgn*visci*outputs%tauw2
      end if
    case(2)
      if (ivel == 1) then
        bcu%y(:,:,ibound) = sgn*visci*outputs%tauw1
      else
        bcw%y(:,:,ibound) = sgn*visci*outputs%tauw2
      end if
    case(3)
      if (ivel == 1) then
        bcu%z(:,:,ibound) = sgn*visci*outputs%tauw1
      else
        bcv%z(:,:,ibound) = sgn*visci*outputs%tauw2
      end if
    end select
  end subroutine assign_wallmodelbc

  subroutine wallmodel_DRL(visc, hwm, inputs, outputs)
    real(rp), intent(in) :: visc, hwm
    type(WallModelInput), intent(in) :: inputs
    type(WallModelOutput), intent(out) :: outputs
    integer :: n1, n2, i, j
    real(rp) :: upar, utau, conv, utau_old, f, fp, tauw_tot
    logical, save :: is_first = .true.

    n1 = size(inputs%vel1, 1) - 2
    n2 = size(inputs%vel1, 2) - 2

    allocate(outputs%tauw1(0:n1+1, 0:n2+1))
    allocate(outputs%tauw2(0:n1+1, 0:n2+1))
    outputs%tauw1 = 0._rp
    outputs%tauw2 = 0._rp

    if (is_first) then
      is_first = .false.
      call InitSmartRedis()
    end if
    CALL ExchangeDataSmartRedis(inputs%vel1, inputs%vel2, outputs%tauw1, outputs%tauw2)

  end subroutine wallmodel_DRL
  
  subroutine wallmodel_loglaw(visc, hwm, inputs, outputs)
    real(rp), intent(in) :: visc, hwm
    type(WallModelInput), intent(in) :: inputs
    type(WallModelOutput), intent(out) :: outputs
    integer :: n1, n2, i, j
    real(rp) :: upar, utau, conv, utau_old, f, fp, tauw_tot

    n1 = size(inputs%vel1, 1) - 2
    n2 = size(inputs%vel1, 2) - 2

    allocate(outputs%tauw1(0:n1+1, 0:n2+1))
    allocate(outputs%tauw2(0:n1+1, 0:n2+1))
    outputs%tauw1 = 0._rp
    outputs%tauw2 = 0._rp

    do j = 0, n2+1
      do i = 0, n1+1
        upar = sqrt(inputs%vel1(i,j)**2 + inputs%vel2(i,j)**2)
        utau = max(sqrt(upar/hwm*visc), visc/hwm*exp(-kap_log*b_log))
        conv = 1._rp
        do while(conv > 0.5e-4_rp)
          utau_old = utau
          f = upar/utau - 1._rp/kap_log*log(hwm*utau/visc) - b_log
          fp = -1._rp/utau*(upar/utau + 1._rp/kap_log)
          utau = abs(utau - f / fp)
          conv = abs(utau / utau_old - 1._rp)
        end do
        tauw_tot = utau**2
        outputs%tauw1(i,j) = tauw_tot * inputs%vel1(i,j) / (upar + eps)
        outputs%tauw2(i,j) = tauw_tot * inputs%vel2(i,j) / (upar + eps)
      end do
    end do
  end subroutine wallmodel_loglaw

  subroutine wallmodel_laminar(visc, hwm, inputs, outputs)
    real(rp), intent(in) :: visc, hwm
    type(WallModelInput), intent(in) :: inputs
    type(WallModelOutput), intent(out) :: outputs
    integer :: i, j
    real(rp) :: upar, umax, del, tauw_tot
    integer :: n1, n2

    n1 = size(inputs%vel1, 1) - 2
    n2 = size(inputs%vel1, 2) - 2

    allocate(outputs%tauw1(0:n1+1, 0:n2+1))
    allocate(outputs%tauw2(0:n1+1, 0:n2+1))
    outputs%tauw1 = 0._rp
    outputs%tauw2 = 0._rp

    del = 1._rp

    do j = 0, n2+1
      do i = 0, n1+1
        upar = sqrt(inputs%vel1(i,j)**2 + inputs%vel2(i,j)**2)
        umax = upar / (hwm / del * (2._rp - hwm / del))
        tauw_tot = 2._rp / del * umax * visc
        outputs%tauw1(i,j) = tauw_tot * inputs%vel1(i,j) / (upar + eps)
        outputs%tauw2(i,j) = tauw_tot * inputs%vel2(i,j) / (upar + eps)
      end do
    end do
  end subroutine wallmodel_laminar

  function vel_relative(v1,v2,coef,bcv_mag)
    !
    ! compute relative velocity to a wall
    !
    implicit none
    real(rp), intent(in) :: v1,v2,coef,bcv_mag
    real(rp) :: vel_relative
    !
    !$acc routine seq
    vel_relative = (1._rp-coef)*v1 + coef*v2
    vel_relative = vel_relative - bcv_mag
  end function vel_relative

end module mod_wallmodel