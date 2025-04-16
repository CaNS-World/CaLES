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
  use mod_common_mpi, only: ierr,myid
  use smartredis_mpi, only: init_smartredis_mpi, put_step_type, put_state, &
                            put_reward, get_action
  use mod_precision, only: rp
  use mod_typedef,   only: Bound, BoundProfile, BoundInteger
  use mod_params,    only: kap_log, b_log, eps, tag, db_clustered, &
                           agent_interval, action_interval, &
                           hwm_min, hwm_max, tauw_ref_min, tauw_ref_max, &
                           cfd_seed
  use mod_bound,     only: boundp
  implicit none
  private
  public :: init_wallmodel_heights, compute_and_apply_wall_stress

  integer, parameter :: WM_LOG = 1  ! Log law
  integer, parameter :: WM_LAM = 2  ! Laminar
  integer, parameter :: WM_DRL = 3  ! DRL
  
  type :: WallState
    type(Bound)        :: vel1, vel2, vel, hwm, visc, &
                          s1, s2, s1_old
    type(BoundInteger) :: hwm_idx
  end type WallState

  type :: WallStress
    type(Bound) :: tauw1, tauw2, tauw, & 
                   tauw1_prev, tauw2_prev, tauw_prev
  end type WallStress

  type :: PerformanceMetric
    type(Bound)        :: tauw1, tauw2, tauw, &
                          tauw1_prev, tauw2_prev, tauw_prev, & 
                          vel1, vel1_profile_err
    type(BoundProfile) :: vel1_profile
  end type PerformanceMetric

  type :: FlattenedState
    real(rp), allocatable, dimension(:) :: vel1, vel2, vel, hwm, visc, &
                                           s1, s2, s1_old
    integer,  allocatable, dimension(:) :: hwm_idx
  end type FlattenedState

  type :: FlattenedStress
    real(rp), allocatable, dimension(:) :: tauw1, tauw2, tauw, &
                                           tauw1_prev, tauw2_prev, tauw_prev
  end type FlattenedStress

  type :: FlattenedMetric
    real(rp), allocatable, dimension(:)   :: tauw1, tauw2, tauw, &
                                             tauw1_prev, tauw2_prev, tauw_prev, &
                                             vel1, vel1_profile_err
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
  subroutine init_wallmodel_heights(n, is_bound, lwm, l, dl, zc, hwm, hwm_idx, state)
    implicit none
    integer, intent(in), dimension(3)      :: n
    logical, intent(in), dimension(0:1, 3) :: is_bound
    integer, intent(in), dimension(0:1, 3) :: lwm
    real(rp), intent(in), dimension(3)     :: l, dl
    real(rp), intent(in), dimension(0:)    :: zc
    real(rp), intent(in)                   :: hwm
    integer, intent(out), dimension(0:1, 3) :: hwm_idx
    type(WallState), intent(inout)         :: state
    integer, allocatable, dimension(:)     :: seed
    integer :: i, j, k, i1, i2, j1, j2, k1, k2, seed_size
    real(rp), allocatable, dimension(:,:) :: random_values

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

    call random_seed(size = seed_size)
    allocate(seed(seed_size))
    do k = 1, seed_size
      seed(k) = cfd_seed + k + myid * 100
    end do
    call random_seed(put = seed)

    allocate(random_values(n(1), n(2)))

    if (is_bound(0, 3) .and. lwm(0, 3) /= 0) then
      call random_number(random_values)
      state%hwm%z(1:n(1), 1:n(2), 0) = hwm_min + random_values * (hwm_max - hwm_min)
      
      do j = 1, n(2)
        do i = 1, n(1)
          k = 1
          do while (zc(k) < state%hwm%z(i, j, 0))
            k = k + 1
          end do
          k2 = k
          k1 = k - 1
          state%hwm_idx%z(i, j, 0) = k2
        end do
      end do
    end if

    if (is_bound(1, 3) .and. lwm(1, 3) /= 0) then
      call random_number(random_values)
      state%hwm%z(1:n(1), 1:n(2), 1) = hwm_min + random_values * (hwm_max - hwm_min)
      
      do j = 1, n(2)
        do i = 1, n(1)
          k = n(3)
          do while (l(3) - zc(k) < state%hwm%z(i, j, 1))
            k = k - 1
          end do
          k2 = k
          k1 = k + 1
          state%hwm_idx%z(i, j, 1) = k2
        end do
      end do
    end if

    deallocate(random_values)
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
    real(rp), allocatable, dimension(:,:,:,:), save :: stress_field
    integer , allocatable, dimension(:)   :: seed
    real(rp), allocatable, dimension(:,:) :: random_values

    logical, save :: is_first = .true.
    integer, save :: istep = 0
    
    integer, save :: n_points, n_points_x, n_points_y, n_points_z
    integer, save :: interval(3)
    integer, dimension(0:1, 3), save :: hwm_idx
    
    integer :: mtype, idir, ibound, cell_index
    integer :: i, j, k, i0, i1, j0, j1, i_point, i_var
    integer :: seed_size

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
      call init_bound(wall_state%s1, n, 0._rp)
      call init_bound(wall_state%s2, n, 0._rp)
      call init_bound(wall_state%s1_old, n, 0._rp)
      call init_bound_integer(wall_state%hwm_idx, n, 0)
  
      call init_bound(wall_stress%tauw1, n, 0._rp)
      call init_bound(wall_stress%tauw2, n, 0._rp)
      call init_bound(wall_stress%tauw, n, 0._rp)
      call init_bound(wall_stress%tauw1_prev, n, 0._rp)
      call init_bound(wall_stress%tauw2_prev, n, 0._rp)
      call init_bound(wall_stress%tauw_prev, n, 0._rp)
      
      call init_bound(performance_metric%tauw1, n, 0._rp)
      call init_bound(performance_metric%tauw2, n, 0._rp)
      call init_bound(performance_metric%tauw, n, 0._rp)
      call init_bound(performance_metric%tauw1_prev, n, 0._rp)
      call init_bound(performance_metric%tauw2_prev, n, 0._rp)
      call init_bound(performance_metric%tauw_prev, n, 0._rp)
      call init_bound(performance_metric%vel1, n, 0._rp)
      call init_bound(performance_metric%vel1_profile_err, n, 0._rp)
      call init_bound_profile(performance_metric%vel1_profile, n, 0._rp)

      call init_wallmodel_heights(n, is_bound, lwm, l, dl, zc, hwm, hwm_idx, wall_state)
      !
      ! Initialize wall_stress required for the first step to compute the state
      !
      call random_seed(size = seed_size)
      allocate(seed(seed_size))
      do k = 1, seed_size
        seed(k) = cfd_seed + k + myid * 100
      end do
      call random_seed(put = seed)

      allocate(random_values(0:n(1)+1, 0:n(2)+1))
      call random_number(random_values(0:n(1)+1, 0:n(2)+1))
      random_values = tauw_ref_min + random_values * (tauw_ref_max - tauw_ref_min)
      wall_stress%tauw1%z(0:n(1)+1, 0:n(2)+1, 0) = random_values
      wall_stress%tauw2%z(0:n(1)+1, 0:n(2)+1, 0) = 0._rp
      wall_stress%tauw %z(0:n(1)+1, 0:n(2)+1, 0) = random_values
      call random_number(random_values(0:n(1)+1, 0:n(2)+1))
      random_values = tauw_ref_min + random_values * (tauw_ref_max - tauw_ref_min)
      wall_stress%tauw1%z(0:n(1)+1, 0:n(2)+1, 1) = random_values
      wall_stress%tauw2%z(0:n(1)+1, 0:n(2)+1, 1) = 0._rp
      wall_stress%tauw %z(0:n(1)+1, 0:n(2)+1, 1) = random_values
      deallocate(random_values)
      
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
      allocate(flattened_state%s1(n_points));      flattened_state%s1      = 0._rp
      allocate(flattened_state%s2(n_points));      flattened_state%s2      = 0._rp
      allocate(flattened_state%s1_old(n_points));  flattened_state%s1_old  = 0._rp
      allocate(flattened_state%hwm_idx(n_points)); flattened_state%hwm_idx = 0

      allocate(flattened_stress%tauw1(n_points)); flattened_stress%tauw1 = 0._rp
      allocate(flattened_stress%tauw2(n_points)); flattened_stress%tauw2 = 0._rp
      allocate(flattened_stress%tauw(n_points));  flattened_stress%tauw  = 0._rp
      allocate(flattened_stress%tauw1_prev(n_points)); flattened_stress%tauw1_prev = 0._rp
      allocate(flattened_stress%tauw2_prev(n_points)); flattened_stress%tauw2_prev = 0._rp
      allocate(flattened_stress%tauw_prev(n_points));  flattened_stress%tauw_prev  = 0._rp
      
      allocate(flattened_metric%tauw1(n_points)); flattened_metric%tauw1 = 0._rp
      allocate(flattened_metric%tauw2(n_points)); flattened_metric%tauw2 = 0._rp
      allocate(flattened_metric%tauw(n_points));  flattened_metric%tauw  = 0._rp
      allocate(flattened_metric%tauw1_prev(n_points)); flattened_metric%tauw1_prev = 0._rp
      allocate(flattened_metric%tauw2_prev(n_points)); flattened_metric%tauw2_prev = 0._rp
      allocate(flattened_metric%tauw_prev(n_points));  flattened_metric%tauw_prev  = 0._rp
      allocate(flattened_metric%vel1(n_points));               flattened_metric%vel1             = 0._rp
      allocate(flattened_metric%vel1_profile_err(n_points));   flattened_metric%vel1_profile_err = 0._rp
      allocate(flattened_metric%vel1_profile(n_points, n(3))); flattened_metric%vel1_profile     = 0._rp

      allocate(stress_field(0:n(1)+1, 0:n(2)+1, 0:n(3)+1, 3)); stress_field = 0._rp

      istep = 0
    else
      istep = istep + 1
    end if

    call compute_wall_data(n, is_bound, lwm, l, dl, zc, zf, dzc, dzf, visc, hwm, hwm_idx, &
                           u, v, w, bcu_mag, bcv_mag, bcw_mag, wall_state, wall_stress, &
                           performance_metric)

    if (mod(istep, action_interval) == 0) then

      i_point = 1
      if (is_bound(0, 3) .and. lwm(0, 3) /= 0) then
        call coarsen_and_flatten_wall_data(wall_state, wall_stress, performance_metric, flattened_state, &
                                           flattened_stress, flattened_metric, n, interval, 0, n_points_z, &
                                           i_point)
        i_point = i_point + n_points_z
      end if

      if (is_bound(1, 3) .and. lwm(1, 3) /= 0) then
        call coarsen_and_flatten_wall_data(wall_state, wall_stress, performance_metric, flattened_state, &
                                           flattened_stress, flattened_metric, n, interval, 1, n_points_z, &
                                           i_point)
        i_point = i_point + n_points_z
      end if

      if (i_point /= n_points + 1) then
        print*, 'ERROR: i_point /= n_points + 1.'
      end if

      ! State of wall shear stress, s_n
      if (is_bound(0, 3) .and. lwm(0, 3) /= 0) then
        wall_stress%tauw1_prev%z(:,:,0) = wall_stress%tauw1%z(:,:,0)
        wall_stress%tauw2_prev%z(:,:,0) = wall_stress%tauw2%z(:,:,0)
        wall_stress%tauw_prev %z(:,:,0) = wall_stress%tauw %z(:,:,0)
      end if

      if (is_bound(1, 3) .and. lwm(1, 3) /= 0) then
        wall_stress%tauw1_prev%z(:,:,1) = wall_stress%tauw1%z(:,:,1)
        wall_stress%tauw2_prev%z(:,:,1) = wall_stress%tauw2%z(:,:,1)
        wall_stress%tauw_prev %z(:,:,1) = wall_stress%tauw %z(:,:,1)
      end if

      mtype = maxval(lwm(0:1, 1:3))
      call wallmodel_dispatch_table(mtype)%ptr(visc, hwm, flattened_state, flattened_stress, flattened_metric)

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

      ! State of wall shear stress, s_n+1
      if (is_bound(0, 3) .and. lwm(0, 3) /= 0) then
        wall_stress%tauw1%z(:,:,0) = stress_field(:,:,1   ,1)
        wall_stress%tauw2%z(:,:,0) = stress_field(:,:,1   ,2)
        wall_stress%tauw %z(:,:,0) = stress_field(:,:,1   ,3)
      end if

      if (is_bound(1, 3) .and. lwm(1, 3) /= 0) then
        wall_stress%tauw1%z(:,:,1) = stress_field(:,:,n(3),1)
        wall_stress%tauw2%z(:,:,1) = stress_field(:,:,n(3),2)
        wall_stress%tauw %z(:,:,1) = stress_field(:,:,n(3),3)
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
                               bcu_mag, bcv_mag, bcw_mag, state, stress, metric)
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
    type(WallState),     intent(inout) :: state
    type(WallStress),    intent(in) :: stress
    type(PerformanceMetric), intent(inout) :: metric
    
    real(rp) :: coef, sgn, u1, u2, v1, v2, w1, w2, u_mag, v_mag, w_mag, uh, vh, wh, this_hwm
    real(rp) :: vel_1, vel_2, vel_h, tauw1, tauw2, tauw, tauw1_prev, tauw2_prev, tauw_prev
    real(rp) :: del_v, dveldz, this_hwm_plus, vel_h_plus, dveldz_plus, kap, b, utau
    real(rp) :: s1, s2, s1_old
    integer  :: i1, i2, j1, j2, k1, k2, i, j, k, ibound, idir, cell_index
    logical, save  :: is_first = .true.
    integer, save  :: istep, n_samples
    real(rp), allocatable, save :: u_ref(:), u_ref_0(:), u_ref_1(:), u_profile(:)
    real(rp), allocatable, save :: u_profile_ave(:,:,:,:)
    
    if (is_first) then
      is_first   = .false.
      istep      = 0
      n_samples  = 1
      metric%vel1%z = 0._rp
      allocate(u_ref_0(n(3)))
      allocate(u_ref_1(n(3)))
      allocate(u_profile_ave(n(3), n(1), n(2), 0:1))
      u_profile_ave = 0._rp
      allocate(u_ref(n(3)))
      allocate(u_profile(n(3)))
    else
      istep       = istep + 1
      n_samples   = n_samples + 1
    end if
    if (mod(istep, action_interval) == 1) then
      n_samples = 1
      metric%vel1%z = 0._rp
      u_profile_ave = 0._rp
    end if

    do idir = 1, 3
      do ibound = 0, 1
        if (is_bound(ibound, idir) .and. lwm(ibound, idir) /= 0) then

          select case(idir)
          case(3)
            do j = 1, n(2)
              do i = 1, n(1)
                cell_index = state%hwm_idx%z(i, j, ibound)
                this_hwm = state%hwm%z(i, j, ibound)
                if (ibound == 0) then
                  k2 = cell_index
                  k1 = cell_index - 1
                  coef = (this_hwm - zc(k1)) / dzc(k1)
                  u_ref = u_ref_0
                  sgn =  1._rp
                else
                  k2 = cell_index
                  k1 = cell_index + 1
                  coef = (this_hwm - (l(3) - zc(k1))) / dzc(k2)
                  u_ref = u_ref_1
                  sgn = -1._rp
                end if
                u1 = 0.5_rp * (u(i - 1, j, k1) + u(i, j, k1))
                v1 = 0.5_rp * (v(i, j - 1, k1) + v(i, j, k1))
                u2 = 0.5_rp * (u(i - 1, j, k2) + u(i, j, k2))
                v2 = 0.5_rp * (v(i, j - 1, k2) + v(i, j, k2))
                uh = vel_relative(u1, u2, coef, 0._rp)
                vh = vel_relative(v1, v2, coef, 0._rp)
                vel_1 = sqrt(u1**2 + v1**2)
                vel_2 = sqrt(u2**2 + v2**2)
                vel_h = sqrt(uh**2 + vh**2)
                !
                ! Local spatial average should benefit the performance
                !
                tauw1_prev = stress%tauw1_prev%z(i, j, ibound) ! s_n
                tauw2_prev = stress%tauw2_prev%z(i, j, ibound) ! s_n
                tauw_prev  = stress%tauw_prev %z(i, j, ibound) ! s_n
                tauw1 = stress%tauw1%z(i, j, ibound) ! s_n+1
                tauw2 = stress%tauw2%z(i, j, ibound) ! s_n+1
                tauw  = stress%tauw %z(i, j, ibound) ! s_n+1
                !
                ! Wall units based on tauw (not tauw1), assuming that
                ! vel_1, vel_2 and tauw are along the same direction,
                ! and that the velocity profile is monotonically increasing.
                ! This is a common assumption in equilibrium wall models.
                ! How about computing dveldz from dudz and dvdz?
                ! Local spatial average should benefit
                ! 
                utau = sqrt(tauw)
                del_v = visc/utau
                dveldz = sgn * (vel_2 - vel_1) / (zc(k2) - zc(k1))
                this_hwm_plus = this_hwm / del_v
                vel_h_plus = vel_h / utau
                dveldz_plus = dveldz * del_v / utau
                kap = 1._rp / (this_hwm_plus * dveldz_plus)
                b = vel_h_plus - 1._rp / kap * log(this_hwm_plus)
                s1_old = 1._rp / kap
                s1 = (1._rp / kap - 1._rp / kap_log) * log(this_hwm_plus)
                s2 = b
                !
                state%vel1  %z(i, j, ibound) = uh
                state%vel2  %z(i, j, ibound) = vh
                state%vel   %z(i, j, ibound) = vel_h
                state%visc  %z(i, j, ibound) = visc
                state%s1_old%z(i, j, ibound) = s1_old
                state%s1    %z(i, j, ibound) = s1
                state%s2    %z(i, j, ibound) = s2
                !
                ! Reward based on tauw1 = tauw_ref and tauw2 = 0
                ! Reward considers the wall shear stress at s_n and s_n+1
                !
                metric%tauw1_prev%z(i, j, ibound) = tauw1_prev ! s_n
                metric%tauw2_prev%z(i, j, ibound) = tauw2_prev ! s_n
                metric%tauw_prev %z(i, j, ibound) = tauw_prev ! s_n
                metric%tauw1     %z(i, j, ibound) = tauw1 ! s_n+1
                metric%tauw2     %z(i, j, ibound) = tauw2 ! s_n+1
                metric%tauw      %z(i, j, ibound) = tauw ! s_n+1

                ! performance_metric%vel1%z(i, j, ibound) = ((n_samples - 1) / float(n_samples)) * performance_metric%vel1%z(i, j, ibound) + &
                !                                           (             1  / float(n_samples)) * uh
                ! u_profile = 0.5_rp * (u(i - 1, j, 1:n(3)) + u(i, j, 1:n(3)))
                ! u_profile_ave(:, i, j, ibound) = ((n_samples - 1) / float(n_samples)) * u_profile_ave(:, i, j, ibound) + &
                !                                  (1               / float(n_samples)) * u_profile
                ! performance_metric%vel1_profile_err%z(i, j, ibound) = sum(dzf(1:n(3)) * (u_profile_ave(:, i, j, ibound) - u_ref)**2)
                ! if (ibound == 0) then
                !   performance_metric%vel1_profile%z(:, i, j, ibound) = u_profile_ave(1:n(3)   , i, j, ibound)
                ! else
                !   performance_metric%vel1_profile%z(:, i, j, ibound) = u_profile_ave(n(3):1:-1, i, j, ibound)
                ! end if
              end do
            end do
          end select
        end if
      end do
    end do
    ! 
  end subroutine compute_wall_data

  function vel_relative(v1, v2, coef, bcv_mag)
    implicit none
    real(rp), intent(in) :: v1, v2, coef, bcv_mag
    real(rp) :: vel_relative
    !$acc routine seq
    vel_relative = (1._rp - coef) * v1 + coef * v2
    vel_relative = vel_relative - bcv_mag
  end function vel_relative

  subroutine coarsen_and_flatten_wall_data(state, stress, metric, flattened_state, flattened_stress, &
                                           flattened_metric, n, interval, ibound, n_points_z, i_point)
    implicit none
    type(WallState), intent(in) :: state
    type(WallStress), intent(in) :: stress
    type(PerformanceMetric), intent(in) :: metric
    type(FlattenedState), intent(inout) :: flattened_state
    type(FlattenedStress), intent(inout) :: flattened_stress
    type(FlattenedMetric), intent(inout) :: flattened_metric
    integer, intent(in), dimension(3) :: n, interval
    integer, intent(in) :: ibound, n_points_z, i_point
    
    flattened_state%vel1   (i_point:i_point+n_points_z-1) = reshape(state%vel1   %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_state%vel2   (i_point:i_point+n_points_z-1) = reshape(state%vel2   %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_state%vel    (i_point:i_point+n_points_z-1) = reshape(state%vel    %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_state%hwm    (i_point:i_point+n_points_z-1) = reshape(state%hwm    %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_state%visc   (i_point:i_point+n_points_z-1) = reshape(state%visc   %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_state%hwm_idx(i_point:i_point+n_points_z-1) = reshape(state%hwm_idx%z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_state%s1     (i_point:i_point+n_points_z-1) = reshape(state%s1     %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_state%s2     (i_point:i_point+n_points_z-1) = reshape(state%s2     %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_state%s1_old (i_point:i_point+n_points_z-1) = reshape(state%s1_old %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))

    flattened_stress%tauw1(i_point:i_point+n_points_z-1) = reshape(stress%tauw1%z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_stress%tauw2(i_point:i_point+n_points_z-1) = reshape(stress%tauw2%z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_stress%tauw (i_point:i_point+n_points_z-1) = reshape(stress%tauw %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    
    flattened_metric%tauw1     (i_point:i_point+n_points_z-1) = reshape(metric%tauw1     %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_metric%tauw2     (i_point:i_point+n_points_z-1) = reshape(metric%tauw2     %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_metric%tauw      (i_point:i_point+n_points_z-1) = reshape(metric%tauw      %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_metric%tauw1_prev(i_point:i_point+n_points_z-1) = reshape(metric%tauw1_prev%z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_metric%tauw2_prev(i_point:i_point+n_points_z-1) = reshape(metric%tauw2_prev%z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    flattened_metric%tauw_prev (i_point:i_point+n_points_z-1) = reshape(metric%tauw_prev %z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))

    ! flattened_metric%vel1(i_point:i_point+n_points_z-1) = reshape(metric%vel1%z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    ! flattened_metric%vel1_profile_err(i_point:i_point+n_points_z-1) = reshape(metric%vel1_profile_err%z(1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z/))
    ! flattened_metric%vel1_profile(i_point:i_point+n_points_z-1, 1:n(3)) = reshape(metric%vel1_profile%z(1:n(3), 1:n(1):interval(1), 1:n(2):interval(2), ibound), (/n_points_z, n(3)/))
  end subroutine coarsen_and_flatten_wall_data

  subroutine wallmodel_loglaw(visc, hwm, flattened_state, flattened_stress, flattened_metric)
    implicit none
    real(rp), intent(in) :: visc, hwm
    type(FlattenedState), intent(in) :: flattened_state
    type(FlattenedStress), intent(inout) :: flattened_stress
    type(FlattenedMetric), intent(in), optional :: flattened_metric
    real(rp) :: u1, u2, upar, utau, conv, utau_old, f, fp, tauw_tot, tauw1, tauw2, this_hwm
    integer :: n_points, i

    n_points = size(flattened_state%vel1)

    do i = 1, n_points
      this_hwm = flattened_state%hwm(i)
      u1 = flattened_state%vel1(i)
      u2 = flattened_state%vel2(i)
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
      flattened_stress%tauw1(i) = tauw1
      flattened_stress%tauw2(i) = tauw2
    end do
  end subroutine wallmodel_loglaw

  subroutine wallmodel_laminar(visc, hwm, flattened_state, flattened_stress, flattened_metric)
    implicit none
    real(rp), intent(in) :: visc, hwm
    type(FlattenedState), intent(in) :: flattened_state
    type(FlattenedStress), intent(inout) :: flattened_stress
    type(FlattenedMetric), intent(in), optional :: flattened_metric
    real(rp) :: u1, u2, upar, umax, del, tauw_tot, tauw1, tauw2, this_hwm
    integer :: n_points, i

    n_points = size(flattened_state%vel1)
    del = 1._rp

    do i = 1, n_points
      this_hwm = flattened_state%hwm(i)
      u1 = flattened_state%vel1(i)
      u2 = flattened_state%vel2(i)
      upar = sqrt(u1**2 + u2**2)
      umax = upar / (this_hwm / del * (2._rp - this_hwm / del))
      tauw_tot = 2._rp / del * umax * visc
      tauw1 = tauw_tot * u1 / (upar + eps)
      tauw2 = tauw_tot * u2 / (upar + eps)
      flattened_stress%tauw1(i) = tauw1
      flattened_stress%tauw2(i) = tauw2
    end do
  end subroutine wallmodel_laminar

  subroutine wallmodel_DRL(visc, hwm, flattened_state, flattened_stress, flattened_metric)
    implicit none
    real(rp), intent(in) :: visc, hwm
    type(FlattenedState), intent(in) :: flattened_state
    type(FlattenedStress), intent(inout) :: flattened_stress
    type(FlattenedMetric), intent(in), optional :: flattened_metric
    real(rp) :: u1, u2, upar, tauw_tot, tauw1, tauw2, factor
    integer :: n_points, i
    logical, save :: is_first = .true.
    integer, save :: istep

    real(rp), allocatable, dimension(:,:), save :: drl_state
    real(rp), allocatable, dimension(:,:), save :: drl_action
    real(rp), allocatable, dimension(:,:), save :: drl_reward

    n_points = size(flattened_state%vel1)

    if (is_first) then
      is_first = .false.
      istep = 0
      call init_smartredis_mpi(db_clustered, MPI_COMM_WORLD)
      allocate(drl_state (2, n_points))
      allocate(drl_reward(2, n_points))
      allocate(drl_action(1, n_points))
    else
      istep = istep + action_interval
    end if

    if (myid == 0) then
      print*, "istep = ", istep
    end if
    
    drl_state(1, :) = flattened_state%s1
    drl_state(2, :) = flattened_state%s2
    call put_state(trim(adjustl(tag))//".state", shape(drl_state), drl_state)
    !
    if (istep /= 0) then
      drl_reward(1, :) = flattened_metric%tauw1
      drl_reward(2, :) = flattened_metric%tauw1_prev
      call put_reward(trim(adjustl(tag))//".reward", shape(drl_reward), drl_reward)
    end if
    !
    call get_action(trim(adjustl(tag))//".action", shape(drl_action), drl_action)
    
    do i = 1, n_points
      factor = drl_action(1, i)
      tauw_tot = factor * flattened_stress%tauw(i) ! tauw_n (s_n) -> tauw_n+1 (s_n+1)
      u1 = flattened_state%vel1(i)
      u2 = flattened_state%vel2(i)
      upar = sqrt(u1**2 + u2**2)
      tauw1 = tauw_tot * u1 / (upar + eps)
      tauw2 = tauw_tot * u2 / (upar + eps)
      flattened_stress%tauw1(i) = tauw1 ! tauw_n+1 (s_n+1)
      flattened_stress%tauw2(i) = tauw2 ! tauw_n+1 (s_n+1)
      flattened_stress%tauw (i) = tauw_tot ! tauw_n+1 (s_n+1)
    end do
  end subroutine wallmodel_DRL

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
  
  subroutine apply_wall_stress_bc(idir, ibound, visc, stress, bcu, bcv, bcw)
    implicit none
    integer, intent(in) :: idir, ibound
    real(rp), intent(in) :: visc
    type(WallStress), intent(in) :: stress
    type(Bound), intent(inout) :: bcu, bcv, bcw
    real(rp) :: visci, sgn
    integer :: nx, ny
    
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
      bcu%z(0:nx, 1:ny, ibound) = sgn * visci * 0.5_rp * (stress%tauw1%z(0:nx  , 1:ny  , ibound) + &
                                                          stress%tauw1%z(1:nx+1, 1:ny  , ibound))
      bcv%z(1:nx, 0:ny, ibound) = sgn * visci * 0.5_rp * (stress%tauw2%z(1:nx  , 0:ny  , ibound) + &
                                                          stress%tauw2%z(1:nx  , 1:ny+1, ibound))
    end select
  end subroutine apply_wall_stress_bc

end module mod_wallmodel