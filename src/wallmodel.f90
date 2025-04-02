! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_wallmodel
  use, intrinsic :: ieee_arithmetic, only: is_nan => ieee_is_nan, is_finite => ieee_is_finite
  use, intrinsic :: ieee_exceptions, only: ieee_get_flag, ieee_set_flag, ieee_underflow
  use mpi
  use smartredis_mpi, only: init_smartredis_mpi, put_step_type, put_state, &
                            put_reward, get_action
  use mod_precision, only: rp
  use mod_typedef,   only: Bound, BoundProfile, BoundInteger
  use mod_params,    only: kap_log, b_log, eps, tag, db_clustered, &
                           total_time_steps, agent_interval, action_interval
  use mod_bound,     only: boundp
  implicit none
  private
  public :: init_wallmodel_heights, compute_and_apply_wall_stress

  integer, parameter :: WM_LOG = 1  ! Log law
  integer, parameter :: WM_LAM = 2  ! Laminar
  integer, parameter :: WM_DRL = 3  ! DRL
  
  type :: WallState
    type(Bound)        :: vel1, vel2, vel, hwm, visc
    type(BoundInteger) :: hwm_idx
  end type WallState

  type :: WallStress
    type(Bound) :: tauw1, tauw2, tauw
  end type WallStress

  type :: PerformanceMetric
    type(Bound)        :: vel1, vel1_profile_err
    type(BoundProfile) :: vel1_profile
  end type PerformanceMetric

  type :: FlattenedState
    real(rp), allocatable, dimension(:) :: vel1, vel2, vel, hwm, visc
    integer,  allocatable, dimension(:) :: hwm_idx
  end type FlattenedState

  type :: FlattenedStress
    real(rp), allocatable, dimension(:) :: tauw1, tauw2, tauw
  end type FlattenedStress

  type :: FlattenedMetric
    real(rp), allocatable, dimension(:)   :: vel1, vel1_profile_err
    real(rp), allocatable, dimension(:,:) :: vel1_profile
  end type FlattenedMetric

  abstract interface
    subroutine wallmodel_interface(visc, hwm, state, stress, metric)
      import :: rp, FlattenedState, FlattenedStress, FlattenedMetric
      real(rp),               intent(in)    :: visc, hwm
      type(FlattenedState),   intent(in)    :: state
      type(FlattenedStress),  intent(inout) :: stress
      type(FlattenedMetric),  intent(in), optional :: metric
    end subroutine wallmodel_interface
  end interface
  
  type :: WallModelProcedure
    procedure(wallmodel_interface), pointer, nopass :: ptr => null()
  end type WallModelProcedure

  type(WallModelProcedure) :: wallmodel_dispatch_table(3)

  contains

  subroutine init_bound(var, n, initial_val)
    type(Bound), intent(inout) :: var
    integer, intent(in), dimension(3) :: n
    real(rp), intent(in) :: initial_val
    
    allocate(var%x(0:n(2)+1, 0:n(3)+1, 0:1)); var%x = initial_val
    allocate(var%y(0:n(1)+1, 0:n(3)+1, 0:1)); var%y = initial_val
    allocate(var%z(0:n(1)+1, 0:n(2)+1, 0:1)); var%z = initial_val
  end subroutine init_bound

  subroutine init_bound_integer(var, n, initial_val)
    type(BoundInteger), intent(inout) :: var
    integer, intent(in), dimension(3) :: n
    integer, intent(in) :: initial_val
    
    allocate(var%x(0:n(2)+1, 0:n(3)+1, 0:1)); var%x = initial_val
    allocate(var%y(0:n(1)+1, 0:n(3)+1, 0:1)); var%y = initial_val
    allocate(var%z(0:n(1)+1, 0:n(2)+1, 0:1)); var%z = initial_val
  end subroutine init_bound_integer

  subroutine init_bound_profile(var, n, initial_val)
    type(BoundProfile), intent(inout) :: var
    integer, intent(in), dimension(3) :: n
    real(rp), intent(in) :: initial_val
    
    allocate(var%x(1:n(1), 0:n(2)+1, 0:n(3)+1, 0:1)); var%x = initial_val
    allocate(var%y(1:n(2), 0:n(1)+1, 0:n(3)+1, 0:1)); var%y = initial_val
    allocate(var%z(1:n(3), 0:n(1)+1, 0:n(2)+1, 0:1)); var%z = initial_val
  end subroutine init_bound_profile

  ! Find the cell_index required for interpolation to the wall model height.
  ! The stored cell_index corresponds to the cells far from a wall, i.e., i2, j2, k2.
  ! Remmeber to set hwm strightly higher than the first cell center, and lower
  ! than the last cell center (hwm=hwm-eps)
  subroutine init_wallmodel_heights(n, is_bound, lwm, l, dl, zc, hwm, hwm_idx, wall_state)
    implicit none
    integer, intent(in), dimension(3)      :: n
    logical, intent(in), dimension(0:1, 3) :: is_bound
    integer, intent(in), dimension(0:1, 3) :: lwm
    real(rp), intent(in), dimension(3)     :: l, dl
    real(rp), intent(in), dimension(0:)    :: zc
    real(rp), intent(in)                   :: hwm
    integer, intent(out), dimension(0:1, 3) :: hwm_idx
    type(WallState), intent(inout)         :: wall_state
    integer, allocatable, dimension(:)     :: seed
    integer :: i, j, k, i1, i2, j1, j2, k1, k2, ncells

    if (is_bound(0, 1) .and. lwm(0, 1) /= 0) then ! to remove statement
      i = 1
      do while ((i - 0.5) * dl(1) < hwm)
        i = i + 1
      end do
      i2 = i
      i1 = i - 1
      hwm_idx(0, 1) = i2
    end if

    if (is_bound(1, 1) .and. lwm(1, 1) /= 0) then
      i = n(1)
      do while ((n(1) - i + 0.5) * dl(1) < hwm)
        i = i - 1
      end do
      i2 = i
      i1 = i + 1
      hwm_idx(1, 1) = i2
    end if

    if (is_bound(0, 2) .and. lwm(0, 2) /= 0) then
      j = 1
      do while ((j - 0.5) * dl(2) < hwm)
        j = j + 1
      end do
      j2 = j
      j1 = j - 1
      hwm_idx(0, 2) = j2
    end if

    if (is_bound(1, 2) .and. lwm(1, 2) /= 0) then
      j = n(2)
      do while ((n(2) - j + 0.5) * dl(2) < hwm)
        j = j - 1
      end do
      j2 = j
      j1 = j + 1
      hwm_idx(1, 2) = j2
    end if

    if (is_bound(0, 3) .and. lwm(0, 3) /= 0) then
      ncells = n(1) * n(2)
      call random_seed(size = ncells)
      allocate(seed(ncells))
      seed = 123
      call random_seed(put = seed)
      call random_number(wall_state%hwm%z(1:n(1), 1:n(2), 0))
      wall_state%hwm%z(1:n(1), 1:n(2), 0) = 0.075_rp + (0.150_rp - 0.075_rp) * wall_state%hwm%z(1:n(1), 1:n(2), 0)
      
      do j = 1, n(2)
        do i = 1, n(1)
          k = 1
          do while (zc(k) < wall_state%hwm%z(i, j, 0))
            k = k + 1
          end do
          k2 = k
          k1 = k - 1
          wall_state%hwm_idx%z(i, j, 0) = k2
        end do
      end do
    end if

    if (is_bound(1, 3) .and. lwm(1, 3) /= 0) then
      ncells = n(1) * n(2)
      call random_seed(size = ncells)
      if (allocated(seed)) deallocate(seed)
      allocate(seed(ncells))
      seed = 123
      call random_seed(put = seed)
      call random_number(wall_state%hwm%z(1:n(1), 1:n(2), 1))
      wall_state%hwm%z(1:n(1), 1:n(2), 1) = 0.075_rp + (0.150_rp - 0.075_rp) * wall_state%hwm%z(1:n(1), 1:n(2), 1)
      
      do j = 1, n(2)
        do i = 1, n(1)
          k = n(3)
          do while (l(3) - zc(k) < wall_state%hwm%z(i, j, 1))
            k = k - 1
          end do
          k2 = k
          k1 = k + 1
          wall_state%hwm_idx%z(i, j, 1) = k2
        end do
      end do
    end if
  end subroutine init_wallmodel_heights

  subroutine compute_and_apply_wall_stress(n, nb, is_bound, lwm, l, dl, zc, zf, dzc, dzf, visc, hwm, u, v, w, &
                                           cbcsgs, bcu, bcv, bcw, bcs, bcu_mag, bcv_mag, bcw_mag)
    implicit none
    integer, intent(in), dimension(3)      :: n
    integer, intent(in), dimension(0:1, 3) :: nb
    logical, intent(in), dimension(0:1, 3) :: is_bound
    integer, intent(in), dimension(0:1, 3) :: lwm
    real(rp), intent(in), dimension(3)     :: l, dl
    real(rp), intent(in), dimension(0:)    :: zc, zf, dzc, dzf
    real(rp), intent(in)                   :: visc, hwm
    real(rp), intent(in), dimension(0:, 0:, 0:) :: u, v, w
    character(len=1), intent(in), dimension(0:1, 3) :: cbcsgs
    type(Bound), intent(inout) :: bcu, bcv, bcw
    type(Bound), intent(in) :: bcs
    type(Bound), intent(in) :: bcu_mag, bcv_mag, bcw_mag
    type(WallState), save :: wall_state
    type(WallStress), save :: wall_stress
    type(PerformanceMetric), save :: performance_metric
    type(FlattenedState), save :: flattened_state
    type(FlattenedStress), save :: flattened_stress
    type(FlattenedMetric), save :: flattened_metric
    real(rp), allocatable, dimension(:,:,:,:), save:: stress_field

    logical, save :: is_first = .true.
    integer, save :: istep = 0
    
    integer, save :: n_points, n_points_x, n_points_y, n_points_z
    integer, save :: interval(3)
    integer, dimension(0:1, 3), save :: hwm_idx
    
    integer :: mtype, idir, ibound, cell_index
    integer :: i, j, k, i0, i1, j0, j1, ierr, i_point, i_var

    real(rp) :: tmp, tmp2

    if (is_first) then
      is_first = .false.
      
      wallmodel_dispatch_table(WM_LOG)%ptr => wallmodel_loglaw
      wallmodel_dispatch_table(WM_LAM)%ptr => wallmodel_laminar
      wallmodel_dispatch_table(WM_DRL)%ptr => wallmodel_DRL
      
      call init_bound(wall_state%vel1, n, 0._rp)
      call init_bound(wall_state%vel2, n, 0._rp)
      call init_bound(wall_state%vel, n, 0._rp)
      call init_bound(wall_state%hwm, n, 0._rp)
      call init_bound(wall_state%visc, n, 0._rp)
      call init_bound_integer(wall_state%hwm_idx, n, 0)
  
      call init_bound(wall_stress%tauw1, n, 0._rp)
      call init_bound(wall_stress%tauw2, n, 0._rp)
      call init_bound(wall_stress%tauw, n, 0._rp)
      
      call init_bound(performance_metric%vel1, n, 0._rp)
      call init_bound(performance_metric%vel1_profile_err, n, 0._rp)
      call init_bound_profile(performance_metric%vel1_profile, n, 0._rp)

      call init_wallmodel_heights(n, is_bound, lwm, l, dl, zc, hwm, hwm_idx, wall_state)

      n_points = 0
      interval(1:3) = agent_interval  ! agent_interval is a single value
      
      n_points_x = (n(2)/interval(2)) * (n(3)/interval(3))
      n_points_y = (n(1)/interval(1)) * (n(3)/interval(3))
      n_points_z = (n(1)/interval(1)) * (n(2)/interval(2))
      
      if (is_bound(0, 1) .and. lwm(0, 1) /= 0) n_points = n_points + n_points_x
      if (is_bound(1, 1) .and. lwm(1, 1) /= 0) n_points = n_points + n_points_x
      if (is_bound(0, 2) .and. lwm(0, 2) /= 0) n_points = n_points + n_points_y
      if (is_bound(1, 2) .and. lwm(1, 2) /= 0) n_points = n_points + n_points_y
      if (is_bound(0, 3) .and. lwm(0, 3) /= 0) n_points = n_points + n_points_z
      if (is_bound(1, 3) .and. lwm(1, 3) /= 0) n_points = n_points + n_points_z

      allocate(flattened_state%vel1(n_points));    flattened_state%vel1    = 0._rp
      allocate(flattened_state%vel2(n_points));    flattened_state%vel2    = 0._rp
      allocate(flattened_state%vel(n_points));     flattened_state%vel     = 0._rp
      allocate(flattened_state%hwm(n_points));     flattened_state%hwm     = 0._rp
      allocate(flattened_state%visc(n_points));    flattened_state%visc    = 0._rp
      allocate(flattened_state%hwm_idx(n_points)); flattened_state%hwm_idx = 0

      allocate(flattened_stress%tauw1(n_points)); flattened_stress%tauw1 = 0._rp
      allocate(flattened_stress%tauw2(n_points)); flattened_stress%tauw2 = 0._rp
      allocate(flattened_stress%tauw(n_points));  flattened_stress%tauw  = 0._rp
      
      allocate(flattened_metric%vel1(n_points));               flattened_metric%vel1             = 0._rp
      allocate(flattened_metric%vel1_profile_err(n_points));   flattened_metric%vel1_profile_err = 0._rp
      allocate(flattened_metric%vel1_profile(n_points, n(3))); flattened_metric%vel1_profile     = 0._rp

      allocate(stress_field(0:n(1)+1, 0:n(2)+1, 0:n(3)+1, 3)); 
      stress_field = 0._rp
      
      istep = 0
    else
      istep = istep + 1
    end if

    call compute_wall_data(n, is_bound, lwm, l, dl, zc, zf, dzc, dzf, visc, hwm, hwm_idx, &
                           u, v, w, bcu_mag, bcv_mag, bcw_mag, wall_state, performance_metric)

    if (mod(istep, action_interval) == 0) then

      i_point = 1
      if (is_bound(0, 3) .and. lwm(0, 3) /= 0) then
        call coarsen_and_flatten_wall_data(wall_state, performance_metric, flattened_state, &
                                           flattened_metric, n, interval, 0, n_points_z, &
                                           i_point)
        i_point = i_point + n_points_z
      end if

      if (is_bound(1, 3) .and. lwm(1, 3) /= 0) then
        call coarsen_and_flatten_wall_data(wall_state, performance_metric, flattened_state, &
                                           flattened_metric, n, interval, 1, n_points_z, &
                                           i_point)
        i_point = i_point + n_points_z
      end if
      
      if (i_point /= n_points + 1) then
        print*, 'ERROR: i_point /= n_points + 1.'
      end if

      mtype = maxval(lwm(0:1, 1:3))
      call wallmodel_dispatch_table(mtype)%ptr(visc, hwm, flattened_state, flattened_stress, flattened_metric)

      i_point = 1
      if (is_bound(0, 3) .and. lwm(0, 3) /= 0) then
        flattened_metric%vel1_profile_err(i_point:i_point+n_points_z-1) = flattened_stress%tauw1(i_point:i_point+n_points_z-1)
        i_point = i_point + n_points_z
      end if

      if (is_bound(1, 3) .and. lwm(1, 3) /= 0) then
        flattened_metric%vel1_profile_err(i_point:i_point+n_points_z-1) = flattened_stress%tauw1(i_point:i_point+n_points_z-1)
        i_point = i_point + n_points_z
      end if

      i_point = 1
      if (is_bound(0, 3) .and. lwm(0, 3) /= 0) then
        call map_stress_to_sparse_grid(flattened_stress, stress_field, n, interval, 0, n_points_z, i_point)
        i_point = i_point + n_points_z
      end if

      if (is_bound(1, 3) .and. lwm(1, 3) /= 0) then
        call map_stress_to_sparse_grid(flattened_stress, stress_field, n, interval, 1, n_points_z, i_point)
        i_point = i_point + n_points_z
      end if

      ! n + 1 filled, but 0 not filled
      ! cbcsgs and bcs must be assigned. bcs is zero on the walls, so the cells have opposite values
      ! on both sides of the wall, which brings zero wall shear stress
      ! square duct has not been tested
      do i_var = 1, 3
        call boundp(cbcsgs, n, bcs, nb, is_bound, dl, dzc, stress_field(:,:,:,i_var))
      end do

      if (is_bound(0, 3) .and. lwm(0, 3) /= 0) then
        call interpolate_stress_field(stress_field, n, interval, 1)
      end if

      if (is_bound(1, 3) .and. lwm(1, 3) /= 0) then
        call interpolate_stress_field(stress_field, n, interval, n(3))
      end if

      ! 0 and other empty ghost cells filled
      do i_var = 1, 3
        call boundp(cbcsgs, n, bcs, nb, is_bound, dl, dzc, stress_field(:,:,:,i_var))
      end do

      if (is_bound(0, 3) .and. lwm(0, 3) /= 0) then
        wall_stress%tauw1%z(:,:,0) = stress_field(:,:,1   ,1)
        wall_stress%tauw2%z(:,:,0) = stress_field(:,:,1   ,2)
      end if

      if (is_bound(1, 3) .and. lwm(1, 3) /= 0) then
        wall_stress%tauw1%z(:,:,1) = stress_field(:,:,n(3),1)
        wall_stress%tauw2%z(:,:,1) = stress_field(:,:,n(3),2)
      end if

      do idir = 1, 3
        do ibound = 0, 1
          if (is_bound(ibound, idir) .and. lwm(ibound, idir) /= 0) then
            call apply_wall_stress_bc(idir, ibound, visc, wall_stress, bcu, bcv, bcw)
          end if
        end do
      end do
    end if
  end subroutine compute_and_apply_wall_stress

  subroutine compute_wall_data(n, is_bound, lwm, l, dl, zc, zf, dzc, dzf, visc, hwm, hwm_idx, u, v, w, &
                               bcu_mag, bcv_mag, bcw_mag, wall_state, performance_metric)
    implicit none
    integer,  intent(in), dimension(3)      :: n
    logical,  intent(in), dimension(0:1, 3) :: is_bound
    integer,  intent(in), dimension(0:1, 3) :: lwm
    real(rp), intent(in), dimension(3)      :: l, dl
    real(rp), intent(in), dimension(0:)     :: zc, zf, dzc, dzf
    real(rp), intent(in)                    :: visc, hwm
    integer,  intent(in), dimension(0:1, 3) :: hwm_idx
    real(rp), intent(in), dimension(0:,0:,0:) :: u, v, w
    type(Bound), intent(in) :: bcu_mag, bcv_mag, bcw_mag
    type(WallState),     intent(inout) :: wall_state
    type(PerformanceMetric), intent(inout) :: performance_metric
    
    real(rp) :: coef, wei, u1, u2, v1, v2, w1, w2, u_mag, v_mag, w_mag, uh, vh, wh, this_hwm
    integer  :: i1, i2, j1, j2, k1, k2, i, j, k, ibound, idir, cell_index
    logical, save  :: is_first = .true.
    integer, save  :: istep, n_samples
    real(rp), allocatable, save :: u_ref(:), u_ref_0(:), u_ref_1(:), u_profile(:)
    real(rp), allocatable, save :: u_profile_ave(:,:,:,:)
    real(rp) :: dummy
    integer  :: ierr
    
    if (is_first) then
      is_first   = .false.
      istep      = 0
      n_samples  = 1
      performance_metric%vel1%z = 0._rp
      allocate(u_ref_0(n(3)))
      allocate(u_ref_1(n(3)))
      allocate(u_profile_ave(n(3), n(1), n(2), 0:1))
      u_profile_ave = 0._rp
      allocate(u_ref(n(3)))
      allocate(u_profile(n(3)))
    else
      istep       = istep + 1
      n_samples = n_samples + 1
    end if
    if (mod(istep, action_interval) == 1) then
      n_samples = 1
      performance_metric%vel1%z = 0._rp
      u_profile_ave = 0._rp
    end if

    do idir = 1, 3
      do ibound = 0, 1
        if (is_bound(ibound, idir) .and. lwm(ibound, idir) /= 0) then

          select case(idir)
          case(3)
            do j = 1, n(2)
              do i = 1, n(1)
                cell_index = wall_state%hwm_idx%z(i, j, ibound)
                this_hwm = wall_state%hwm%z(i, j, ibound)
                if (ibound == 0) then
                  k2 = cell_index
                  k1 = cell_index - 1
                  coef = (this_hwm - zc(k1)) / dzc(k1)
                  u_ref = u_ref_0
                else
                  k2 = cell_index
                  k1 = cell_index + 1
                  coef = (this_hwm - (l(3) - zc(k1))) / dzc(k2)
                  u_ref = u_ref_1
                end if
                u1 = 0.5_rp * (u(i - 1, j, k1) + u(i, j, k1))
                u2 = 0.5_rp * (u(i - 1, j, k2) + u(i, j, k2))
                v1 = 0.5_rp * (v(i, j - 1, k1) + v(i, j, k1))
                v2 = 0.5_rp * (v(i, j - 1, k2) + v(i, j, k2))
                u_mag = 0.5_rp * (bcu_mag%z(i - 1, j, ibound) + bcu_mag%z(i, j, ibound))
                v_mag = 0.5_rp * (bcv_mag%z(i, j - 1, ibound) + bcv_mag%z(i, j, ibound))
                uh = vel_relative(u1, u2, coef, u_mag)
                vh = vel_relative(v1, v2, coef, v_mag)
                wall_state%vel1%z(i, j, ibound) = uh
                wall_state%vel2%z(i, j, ibound) = vh
                wall_state% vel%z(i, j, ibound) = sqrt(uh**2 + vh**2)
                wall_state%visc%z(i, j, ibound) = visc * 20000._rp
                performance_metric%vel1%z(i, j, ibound) = ((n_samples - 1) / float(n_samples)) * performance_metric%vel1%z(i, j, ibound) + &
                                                          (               1  / float(n_samples)) * uh
                u_profile = 0.5_rp * (u(i - 1, j, 1:n(3)) + u(i, j, 1:n(3)))
                u_profile_ave(:, i, j, ibound) = ((n_samples - 1) / float(n_samples)) * u_profile_ave(:, i, j, ibound) + &
                                                 (1                 / float(n_samples)) * u_profile
                performance_metric%vel1_profile_err%z(i, j, ibound) = sum(dzf(1:n(3)) * (u_profile_ave(:, i, j, ibound) - u_ref)**2)
                if (ibound == 0) then
                  performance_metric%vel1_profile%z(:, i, j, ibound) = u_profile_ave(1:n(3)   , i, j, ibound)
                else
                  performance_metric%vel1_profile%z(:, i, j, ibound) = u_profile_ave(n(3):1:-1, i, j, ibound)
                end if
              end do
            end do
          end select
        end if
      end do
    end do
  end subroutine compute_wall_data
  
  subroutine apply_wall_stress_bc(idir, ibound, visc, wall_stress, bcu, bcv, bcw)
    implicit none
    integer, intent(in) :: idir, ibound
    real(rp), intent(in) :: visc
    type(WallStress), intent(in) :: wall_stress
    type(Bound), intent(inout) :: bcu, bcv, bcw
    real(rp) :: visci, sgn
    integer :: nx, ny, ierr
    
    visci = 1._rp / visc
    if (ibound == 0) then
      sgn =  1._rp
    else
      sgn = -1._rp
    end if

    select case(idir)
    case(3)
      nx = size(bcu%z, 1) - 2
      ny = size(bcu%z, 2) - 2
      bcu%z(0:nx, 1:ny, ibound) = sgn * visci * 0.5_rp * (wall_stress%tauw1%z(0:nx  , 1:ny  , ibound) + &
                                                          wall_stress%tauw1%z(1:nx+1, 1:ny  , ibound))
      bcv%z(1:nx, 0:ny, ibound) = sgn * visci * 0.5_rp * (wall_stress%tauw2%z(1:nx  , 0:ny  , ibound) + &
                                                          wall_stress%tauw2%z(1:nx  , 1:ny+1, ibound))
    end select
  end subroutine apply_wall_stress_bc

  subroutine wallmodel_loglaw(visc, hwm, state, stress, metric)
    implicit none
    real(rp), intent(in) :: visc, hwm
    type(FlattenedState), intent(in) :: state
    type(FlattenedStress), intent(inout) :: stress
    type(FlattenedMetric), intent(in), optional :: metric
    real(rp) :: u1, u2, upar, utau, conv, utau_old, f, fp, tauw_tot, tauw1, tauw2, this_hwm
    integer :: n_points, i, i_stag

    n_points = size(state%vel1)

    do i = 1, n_points
      this_hwm = state%hwm(i)
      u1 = state%vel1(i)
      u2 = state%vel2(i)
      upar = sqrt(u1**2 + u2**2)
      utau = max(sqrt(upar / this_hwm * visc), visc / this_hwm * exp(-kap_log * b_log))
      conv = 1._rp
      do while (conv > 0.5e-4_rp)
        utau_old = utau
        f = upar / utau - 1._rp / kap_log * log(this_hwm * utau / visc) - b_log
        fp = -1._rp / utau * (upar / utau + 1._rp / kap_log)
        utau = abs(utau - f / fp)
        conv = abs(utau / utau_old - 1._rp)
      end do
      tauw_tot = utau**2
      tauw1 = tauw_tot * u1 / (upar + eps)
      tauw2 = tauw_tot * u2 / (upar + eps)
      stress%tauw1(i) = tauw1
      stress%tauw2(i) = tauw2
    end do
  end subroutine wallmodel_loglaw

  subroutine wallmodel_laminar(visc, hwm, state, stress, metric)
    implicit none
    real(rp), intent(in) :: visc, hwm
    type(FlattenedState), intent(in) :: state
    type(FlattenedStress), intent(inout) :: stress
    type(FlattenedMetric), intent(in), optional :: metric
    real(rp) :: u1, u2, upar, umax, del, tauw_tot, tauw1, tauw2, this_hwm
    integer :: n_points, i, ierr

    n_points = size(state%vel1)
    del = 1._rp

    do i = 1, n_points
      this_hwm = state%hwm(i)
      u1 = state%vel1(i)
      u2 = state%vel2(i)
      upar = sqrt(u1**2 + u2**2)
      umax = upar / (this_hwm / del * (2._rp - this_hwm / del))
      tauw_tot = 2._rp / del * umax * visc
      tauw1 = tauw_tot * u1 / (upar + eps)
      tauw2 = tauw_tot * u2 / (upar + eps)
      stress%tauw1(i) = tauw1
      stress%tauw2(i) = tauw2
    end do
  end subroutine wallmodel_laminar

  subroutine wallmodel_DRL(visc, hwm, state, stress, metric)
    implicit none
    real(rp), intent(in) :: visc, hwm
    type(FlattenedState), intent(in) :: state
    type(FlattenedStress), intent(inout) :: stress
    type(FlattenedMetric), intent(in), optional :: metric
    real(rp) :: u1, u2, upar, tauw_tot, tauw1, tauw2
    logical, save :: is_first = .true.
    integer, save :: istep, n_points, n_vars_state, n_vars_action, n_vars_reward
    integer :: i, i_var

    real(rp), allocatable, dimension(:,:), save :: state_array
    real(rp), allocatable, dimension(:,:), save :: action
    real(rp), allocatable, dimension(:,:), save :: reward

    if (is_first) then
      is_first = .false.
      istep = 0
      call init_smartredis_mpi(db_clustered, MPI_COMM_WORLD)
      n_points      = size(state%vel1)
      n_vars_state  = 5
      n_vars_action = 3
      n_vars_reward = 2
      allocate(state_array(n_vars_state, n_points))
      allocate(action(n_vars_action, n_points))
      allocate(reward(n_vars_reward, n_points))
      print*, "n_points = ", n_points
    else
      istep = istep + action_interval
    end if
    
    state_array(1, :) = state%vel1
    state_array(2, :) = state%vel2
    state_array(3, :) = state%vel
    state_array(4, :) = state%hwm
    state_array(5, :) = state%visc
    
    if (present(metric)) then
      reward(1, :) = metric%vel1
      reward(2, :) = metric%vel1_profile_err
    end if
    
    print*, "istep = ", istep, "total_time_steps = ", total_time_steps
    call put_state(trim(adjustl(tag))//".state", shape(state_array(3:5,:)), state_array(3:5,:))
    if (istep /= 0) then
      call put_reward(trim(adjustl(tag))//".reward", shape(reward(1:n_vars_reward,:)), reward(1:n_vars_reward,:))
    end if
    call get_action(trim(adjustl(tag))//".action", shape(action(3,:)), action(3,:))
    
    stress%tauw = action(3, :)
    
    do i = 1, n_points
      u1 = state%vel1(i)
      u2 = state%vel2(i)
      upar = sqrt(u1**2 + u2**2)
      tauw_tot = stress%tauw(i)
      tauw1 = tauw_tot * u1 / (upar + eps)
      tauw2 = tauw_tot * u2 / (upar + eps)
      stress%tauw1(i) = tauw1
      stress%tauw2(i) = tauw2
    end do
  end subroutine wallmodel_DRL

  function vel_relative(v1, v2, coef, bcv_mag)
    implicit none
    real(rp), intent(in) :: v1, v2, coef, bcv_mag
    real(rp) :: vel_relative
    !$acc routine seq
    vel_relative = (1._rp - coef) * v1 + coef * v2
    vel_relative = vel_relative - bcv_mag
  end function vel_relative

  subroutine coarsen_and_flatten_wall_data(wall_state, performance_metric, flattened_state, &
                                           flattened_metric, n, interval, ibound, n_points_z, i_point)
    implicit none
    type(WallState), intent(in) :: wall_state
    type(PerformanceMetric), intent(in) :: performance_metric
    type(FlattenedState), intent(inout) :: flattened_state
    type(FlattenedMetric), intent(inout) :: flattened_metric
    integer, intent(in), dimension(3) :: n, interval
    integer, intent(in) :: ibound, n_points_z, i_point
    
    flattened_state%vel1(i_point:i_point+n_points_z-1) = reshape(wall_state%vel1%z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_state%vel2(i_point:i_point+n_points_z-1) = reshape(wall_state%vel2%z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_state%vel (i_point:i_point+n_points_z-1) = reshape(wall_state%vel %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_state%hwm (i_point:i_point+n_points_z-1) = reshape(wall_state%hwm %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_state%visc(i_point:i_point+n_points_z-1) = reshape(wall_state%visc%z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    
    flattened_metric%vel1(i_point:i_point+n_points_z-1) = reshape(performance_metric%vel1%z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_metric%vel1_profile_err(i_point:i_point+n_points_z-1) = reshape(performance_metric%vel1_profile_err%z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    
    flattened_metric%vel1_profile(i_point:i_point+n_points_z-1, 1:n(3)) = reshape(performance_metric%vel1_profile%z(1:n(3), 1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z, n(3)/))
  end subroutine coarsen_and_flatten_wall_data

  subroutine map_stress_to_sparse_grid(flattened_stress, stress_field, n, interval, &
                                       ibound, n_points_z, i_point)
    implicit none
    type(FlattenedStress), intent(in) :: flattened_stress
    real(rp), intent(inout), dimension(0:, 0:, 0:, 1:) :: stress_field
    integer, intent(in), dimension(3) :: n, interval
    integer, intent(in) :: ibound, n_points_z, i_point
    integer :: k_pos
    
    k_pos = merge(1, n(3), ibound == 0)
    stress_field(1:n(1):interval(1), 1:n(2):interval(2), k_pos, 1) = reshape(flattened_stress%tauw1(i_point:i_point+n_points_z-1), (/n(1)/interval(1), n(2)/interval(2)/))
    stress_field(1:n(1):interval(1), 1:n(2):interval(2), k_pos, 2) = reshape(flattened_stress%tauw2(i_point:i_point+n_points_z-1), (/n(1)/interval(1), n(2)/interval(2)/))
    stress_field(1:n(1):interval(1), 1:n(2):interval(2), k_pos, 3) = reshape(flattened_stress%tauw (i_point:i_point+n_points_z-1), (/n(1)/interval(1), n(2)/interval(2)/))
  end subroutine map_stress_to_sparse_grid

  subroutine interpolate_stress_field(stress_field, n, interval, k_pos)
    implicit none
    real(rp), intent(inout) :: stress_field(0:, 0:, 0:, 1:)
    integer, intent(in), dimension(3) :: n, interval
    integer, intent(in) :: k_pos
    integer :: i, j, i0, i1, j0, j1
    
    do j = 1, n(2)
      do i = 1, n(1)
        i0 = ((i - 1) / interval(1)) * interval(1) + 1
        i1 = i0 + interval(1)
        j0 = ((j - 1) / interval(2)) * interval(2) + 1
        j1 = j0 + interval(2)
        stress_field(i, j, k_pos, :) = (stress_field(i0, j0, k_pos, :) * (i1 - i) * (j1 - j) + &
                                        stress_field(i1, j0, k_pos, :) * (i - i0) * (j1 - j) + &
                                        stress_field(i0, j1, k_pos, :) * (i1 - i) * (j - j0) + &
                                        stress_field(i1, j1, k_pos, :) * (i - i0) * (j - j0)) / &
                                        ((i1 - i0) * (j1 - j0))
      end do
    end do
  end subroutine interpolate_stress_field

end module mod_wallmodel