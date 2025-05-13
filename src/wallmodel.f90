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
  use mod_precision,  only: rp
  use mod_typedef,    only: Bound, BoundInteger, BoundProfile
  use mod_params,     only: kap_log, b_log, eps, tag, db_clustered, &
                            agent_interval, action_start_time, time_duration_per_action, &
                            hwm_min, hwm_max, tauw_ref_min, tauw_ref_max, &
                            cfd_seed
  use mod_bound,      only: boundp
  implicit none
  private
  public :: compute_and_assign_wall_stress

  integer, parameter :: WM_LOG = 1  ! Log law
  integer, parameter :: WM_LAM = 2  ! Laminar
  integer, parameter :: WM_DRL = 3  ! DRL
  
  type :: WallState
    type(Bound)        :: vel_x1, vel_x2, hwm, &
                          hwm_plus, velh_plus, dveldz_plus
    type(BoundInteger) :: hwm_idx
  end type WallState

  type :: WallStress
    type(Bound) :: tauw1, tauw2, tauw, &
                   tauw1_prev, tauw2_prev, tauw_prev
  end type WallStress

  type :: PerformanceMetric
    type(Bound)        :: tauw1, tauw2, tauw, & 
                          tauw1_prev, tauw2_prev, tauw_prev
  end type PerformanceMetric

  type :: FlattenedState
    real(rp), allocatable, dimension(:) :: vel_x1, vel_x2, hwm, &
                                           hwm_plus, velh_plus, dveldz_plus
    integer,  allocatable, dimension(:) :: hwm_idx
  end type FlattenedState

  type :: FlattenedStress
    real(rp), allocatable, dimension(:) :: tauw1, tauw2, tauw, &
                                           tauw1_prev, tauw2_prev, tauw_prev
    real(rp), allocatable, dimension(:,:) :: action
  end type FlattenedStress

  type :: FlattenedMetric
    real(rp), allocatable, dimension(:)   :: tauw1, tauw2, tauw, &
                                             tauw1_prev, tauw2_prev, tauw_prev
  end type FlattenedMetric

  abstract interface
    subroutine wallmodel_interface(visc, state, stress, metric)
      import :: rp, FlattenedState, FlattenedStress, FlattenedMetric
      real(rp),               intent(in)    :: visc
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

  subroutine compute_and_assign_wall_stress(n, nb, is_bound, lwm, l, dl, zc, zf, dzc, dzf, visc, u, v, w, &
                                           cbcsgs, bcu, bcv, bcw, bcs, bcu_mag, bcv_mag, bcw_mag, time)
    implicit none
    integer, intent(in), dimension(3)      :: n
    integer, intent(in), dimension(0:1, 3) :: nb
    logical, intent(in), dimension(0:1, 3) :: is_bound
    integer, intent(in), dimension(0:1, 3) :: lwm
    real(rp), intent(in), dimension(3)     :: l, dl
    real(rp), intent(in), dimension(0:)    :: zc, zf, dzc, dzf
    real(rp), intent(in)                   :: visc
    real(rp), intent(in), dimension(0:, 0:, 0:) :: u, v, w
    character(len=1), intent(in), dimension(0:1, 3) :: cbcsgs
    type(Bound), intent(inout) :: bcu, bcv, bcw
    type(Bound), intent(in) :: bcs
    type(Bound), intent(in) :: bcu_mag, bcv_mag, bcw_mag
    real(rp), intent(in) :: time
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
    integer, save :: n_points, n_points_x, n_points_y, n_points_z
    integer, save :: interval(3)
    integer, dimension(0:1, 3), save :: hwm_idx
    
    real(rp), save :: action_next_time
    integer :: mtype, idir, ibound, cell_index
    integer :: i, j, k, i0, i1, j0, j1, i_point, i_var
    integer :: seed_size

    if (is_first) then
      is_first = .false.

      action_next_time = max(time, action_start_time)
      wallmodel_dispatch_table(WM_LOG)%ptr => wallmodel_loglaw
      wallmodel_dispatch_table(WM_LAM)%ptr => wallmodel_laminar
      wallmodel_dispatch_table(WM_DRL)%ptr => wallmodel_DRL

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

      allocate(wall_state%vel_x1     %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               wall_state%vel_x1     %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               wall_state%vel_x1     %z(0:n(1)+1, 0:n(2)+1, 0:1), &
               wall_state%vel_x2     %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               wall_state%vel_x2     %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               wall_state%vel_x2     %z(0:n(1)+1, 0:n(2)+1, 0:1), &
               wall_state%hwm        %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               wall_state%hwm        %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               wall_state%hwm        %z(0:n(1)+1, 0:n(2)+1, 0:1), &
               wall_state%hwm_plus   %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               wall_state%hwm_plus   %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               wall_state%hwm_plus   %z(0:n(1)+1, 0:n(2)+1, 0:1), &
               wall_state%velh_plus  %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               wall_state%velh_plus  %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               wall_state%velh_plus  %z(0:n(1)+1, 0:n(2)+1, 0:1), &
               wall_state%dveldz_plus%x(0:n(2)+1, 0:n(3)+1, 0:1), &
               wall_state%dveldz_plus%y(0:n(1)+1, 0:n(3)+1, 0:1), &
               wall_state%dveldz_plus%z(0:n(1)+1, 0:n(2)+1, 0:1), &
               wall_state%hwm_idx    %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               wall_state%hwm_idx    %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               wall_state%hwm_idx    %z(0:n(1)+1, 0:n(2)+1, 0:1))

      allocate(wall_stress%tauw1     %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               wall_stress%tauw1     %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               wall_stress%tauw1     %z(0:n(1)+1, 0:n(2)+1, 0:1), &
               wall_stress%tauw2     %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               wall_stress%tauw2     %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               wall_stress%tauw2     %z(0:n(1)+1, 0:n(2)+1, 0:1), &
               wall_stress%tauw      %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               wall_stress%tauw      %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               wall_stress%tauw      %z(0:n(1)+1, 0:n(2)+1, 0:1), &
               wall_stress%tauw1_prev%x(0:n(2)+1, 0:n(3)+1, 0:1), &
               wall_stress%tauw1_prev%y(0:n(1)+1, 0:n(3)+1, 0:1), &
               wall_stress%tauw1_prev%z(0:n(1)+1, 0:n(2)+1, 0:1), &
               wall_stress%tauw2_prev%x(0:n(2)+1, 0:n(3)+1, 0:1), &
               wall_stress%tauw2_prev%y(0:n(1)+1, 0:n(3)+1, 0:1), &
               wall_stress%tauw2_prev%z(0:n(1)+1, 0:n(2)+1, 0:1), &
               wall_stress%tauw_prev %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               wall_stress%tauw_prev %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               wall_stress%tauw_prev %z(0:n(1)+1, 0:n(2)+1, 0:1))

      allocate(performance_metric%tauw1     %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               performance_metric%tauw1     %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               performance_metric%tauw1     %z(0:n(1)+1, 0:n(2)+1, 0:1), &
               performance_metric%tauw2     %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               performance_metric%tauw2     %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               performance_metric%tauw2     %z(0:n(1)+1, 0:n(2)+1, 0:1), &
               performance_metric%tauw      %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               performance_metric%tauw      %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               performance_metric%tauw      %z(0:n(1)+1, 0:n(2)+1, 0:1), &
               performance_metric%tauw1_prev%x(0:n(2)+1, 0:n(3)+1, 0:1), &
               performance_metric%tauw1_prev%y(0:n(1)+1, 0:n(3)+1, 0:1), &
               performance_metric%tauw1_prev%z(0:n(1)+1, 0:n(2)+1, 0:1), &
               performance_metric%tauw2_prev%x(0:n(2)+1, 0:n(3)+1, 0:1), &
               performance_metric%tauw2_prev%y(0:n(1)+1, 0:n(3)+1, 0:1), &
               performance_metric%tauw2_prev%z(0:n(1)+1, 0:n(2)+1, 0:1), &
               performance_metric%tauw_prev %x(0:n(2)+1, 0:n(3)+1, 0:1), &
               performance_metric%tauw_prev %y(0:n(1)+1, 0:n(3)+1, 0:1), &
               performance_metric%tauw_prev %z(0:n(1)+1, 0:n(2)+1, 0:1))

      allocate(flattened_state%vel_x1     (n_points), &
               flattened_state%vel_x2     (n_points), &
               flattened_state%hwm        (n_points), &
               flattened_state%hwm_plus   (n_points), &
               flattened_state%velh_plus  (n_points), &
               flattened_state%dveldz_plus(n_points), &
               flattened_state%hwm_idx    (n_points))

      allocate(flattened_stress%tauw1     (n_points), &
               flattened_stress%tauw2     (n_points), &
               flattened_stress%tauw      (n_points), &
               flattened_stress%tauw1_prev(n_points), &
               flattened_stress%tauw2_prev(n_points), &
               flattened_stress%tauw_prev (n_points), &
               flattened_stress%action    (n_points,2))
      
      allocate(flattened_metric%tauw1     (n_points), &
               flattened_metric%tauw2     (n_points), &
               flattened_metric%tauw      (n_points), &
               flattened_metric%tauw1_prev(n_points), &
               flattened_metric%tauw2_prev(n_points), &
               flattened_metric%tauw_prev (n_points))

      allocate(stress_field(0:n(1)+1, 0:n(2)+1, 0:n(3)+1, 3))

      ! Find the cell_index required for interpolation to the wall model height.
      ! The stored cell_index corresponds to the cells far from a wall, i.e., i2, j2, k2.
      ! Remmeber to set hwm strightly higher than the first cell center, and lower
      ! than the last cell center (hwm=hwm-eps)
      !
      call random_seed(size = seed_size)
      allocate(seed(seed_size))
      do k = 1, seed_size
        seed(k) = cfd_seed + k + myid * 100
      end do
      call random_seed(put = seed)
      
      allocate(random_values(0:n(1)+1, 0:n(2)+1))

      if (is_bound(0, 3) .and. lwm(0, 3) /= 0) then
        call random_number(random_values)
        wall_state%hwm%z(0:n(1)+1, 0:n(2)+1, 0) = hwm_min + random_values * (hwm_max - hwm_min)
        do j = 1, n(2)
          do i = 1, n(1)
            k = 1
            do while (zc(k) < wall_state%hwm%z(i, j, 0))
              k = k + 1
            end do
            wall_state%hwm_idx%z(i, j, 0) = k
          end do
        end do
      end if

      if (is_bound(1, 3) .and. lwm(1, 3) /= 0) then
        call random_number(random_values)
        wall_state%hwm%z(0:n(1)+1, 0:n(2)+1, 1) = hwm_min + random_values * (hwm_max - hwm_min)
        do j = 1, n(2)
          do i = 1, n(1)
            k = n(3)
            do while (l(3) - zc(k) < wall_state%hwm%z(i, j, 1))
              k = k - 1
            end do
            wall_state%hwm_idx%z(i, j, 1) = k
          end do
        end do
      end if
      deallocate(random_values)
      !
      ! Initialize wall_stress to compute the initial state and reward (unused)
      !
      wall_stress%tauw1_prev%z(:,:,0) = 0._rp
      wall_stress%tauw2_prev%z(:,:,0) = 0._rp
      wall_stress%tauw_prev %z(:,:,0) = 0._rp
      wall_stress%tauw1_prev%z(:,:,1) = 0._rp
      wall_stress%tauw2_prev%z(:,:,1) = 0._rp
      wall_stress%tauw_prev %z(:,:,1) = 0._rp
      !
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

      !$acc enter data copyin(interval) async(1)

      !$acc enter data async(1) &
      !$acc copyin(wall_state,               &
      !$acc        wall_state%vel_x1     %x, &
      !$acc        wall_state%vel_x1     %y, &
      !$acc        wall_state%vel_x1     %z, &
      !$acc        wall_state%vel_x2     %x, &
      !$acc        wall_state%vel_x2     %y, &
      !$acc        wall_state%vel_x2     %z, &
      !$acc        wall_state%hwm        %x, & 
      !$acc        wall_state%hwm        %y, & 
      !$acc        wall_state%hwm        %z, &
      !$acc        wall_state%hwm_plus   %x, &
      !$acc        wall_state%hwm_plus   %y, &
      !$acc        wall_state%hwm_plus   %z, &
      !$acc        wall_state%velh_plus  %x, &
      !$acc        wall_state%velh_plus  %y, &
      !$acc        wall_state%velh_plus  %z, &
      !$acc        wall_state%dveldz_plus%x, &
      !$acc        wall_state%dveldz_plus%y, &
      !$acc        wall_state%dveldz_plus%z, &
      !$acc        wall_state%hwm_idx    %x, &
      !$acc        wall_state%hwm_idx    %y, &
      !$acc        wall_state%hwm_idx    %z) 

      !$acc enter data async(1) &
      !$acc copyin(wall_stress,              &
      !$acc        wall_stress%tauw1     %x, &
      !$acc        wall_stress%tauw1     %y, &
      !$acc        wall_stress%tauw1     %z, &
      !$acc        wall_stress%tauw2     %x, &
      !$acc        wall_stress%tauw2     %y, &
      !$acc        wall_stress%tauw2     %z, &
      !$acc        wall_stress%tauw      %x, &
      !$acc        wall_stress%tauw      %y, &
      !$acc        wall_stress%tauw      %z, &
      !$acc        wall_stress%tauw1_prev%x, &
      !$acc        wall_stress%tauw1_prev%y, &
      !$acc        wall_stress%tauw1_prev%z, &
      !$acc        wall_stress%tauw2_prev%x, &
      !$acc        wall_stress%tauw2_prev%y, &
      !$acc        wall_stress%tauw2_prev%z, &
      !$acc        wall_stress%tauw_prev %x, &
      !$acc        wall_stress%tauw_prev %y, &
      !$acc        wall_stress%tauw_prev %z)

      !$acc enter data async(1) &
      !$acc create(performance_metric,              &
      !$acc        performance_metric%tauw1     %x, &
      !$acc        performance_metric%tauw1     %y, &
      !$acc        performance_metric%tauw1     %z, &
      !$acc        performance_metric%tauw2     %x, &
      !$acc        performance_metric%tauw2     %y, &
      !$acc        performance_metric%tauw2     %z, &
      !$acc        performance_metric%tauw      %x, &
      !$acc        performance_metric%tauw      %y, &
      !$acc        performance_metric%tauw      %z, &
      !$acc        performance_metric%tauw1_prev%x, &
      !$acc        performance_metric%tauw1_prev%y, &
      !$acc        performance_metric%tauw1_prev%z, &
      !$acc        performance_metric%tauw2_prev%x, &
      !$acc        performance_metric%tauw2_prev%y, &
      !$acc        performance_metric%tauw2_prev%z, &
      !$acc        performance_metric%tauw_prev %x, &
      !$acc        performance_metric%tauw_prev %y, &
      !$acc        performance_metric%tauw_prev %z)

      !$acc enter data async(1) &
      !$acc create(flattened_state            , &
      !$acc        flattened_state%vel_x1     , &
      !$acc        flattened_state%vel_x2     , &
      !$acc        flattened_state%hwm        , &
      !$acc        flattened_state%hwm_plus   , &
      !$acc        flattened_state%velh_plus  , &
      !$acc        flattened_state%dveldz_plus, &
      !$acc        flattened_state%hwm_idx    )

      !$acc enter data async(1) &
      !$acc create(flattened_stress           , &
      !$acc        flattened_stress%tauw1     , &
      !$acc        flattened_stress%tauw2     , &
      !$acc        flattened_stress%tauw      , &
      !$acc        flattened_stress%tauw1_prev, &
      !$acc        flattened_stress%tauw2_prev, &
      !$acc        flattened_stress%tauw_prev , &
      !$acc        flattened_stress%action    )
      !$acc enter data async(1) &
      
      !$acc create(flattened_metric           , &
      !$acc        flattened_metric%tauw1     , &
      !$acc        flattened_metric%tauw2     , &
      !$acc        flattened_metric%tauw      , &
      !$acc        flattened_metric%tauw1_prev, &
      !$acc        flattened_metric%tauw2_prev, &
      !$acc        flattened_metric%tauw_prev )

      !$acc enter data create(stress_field) async(1)
      
    end if

    ! Every cfd step, compute the wall state and performance_metric,
    ! regardless of time_duration_per_action, so time average can be
    ! conducted in compute_wall_data.
    call compute_wall_data(n, is_bound, lwm, l, dl, zc, zf, dzc, dzf, visc, &
                           u, v, w, bcu_mag, bcv_mag, bcw_mag, wall_state, wall_stress, &
                           performance_metric)

    if (time + eps >= action_next_time) then

      action_next_time = action_next_time + time_duration_per_action

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

      ! State of wall shear stress, s_n
      if (is_bound(0, 3) .and. lwm(0, 3) /= 0) then
        !$acc kernels default(present) async(1)
        wall_stress%tauw1_prev%z(:,:,0) = wall_stress%tauw1%z(:,:,0)
        wall_stress%tauw2_prev%z(:,:,0) = wall_stress%tauw2%z(:,:,0)
        wall_stress%tauw_prev %z(:,:,0) = wall_stress%tauw %z(:,:,0)
        !$acc end kernels
      end if

      if (is_bound(1, 3) .and. lwm(1, 3) /= 0) then
        !$acc kernels default(present) async(1)
        wall_stress%tauw1_prev%z(:,:,1) = wall_stress%tauw1%z(:,:,1)
        wall_stress%tauw2_prev%z(:,:,1) = wall_stress%tauw2%z(:,:,1)
        wall_stress%tauw_prev %z(:,:,1) = wall_stress%tauw %z(:,:,1)
        !$acc end kernels
      end if

      mtype = maxval(lwm(0:1, 1:3))
      call wallmodel_dispatch_table(mtype)%ptr(visc, flattened_state, flattened_stress, flattened_metric)

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
        call interpolate_stress_field(stress_field, n, interval, 1   )
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
        !$acc kernels default(present) async(1)
        wall_stress%tauw1%z(:,:,0) = stress_field(:,:,1   ,1)
        wall_stress%tauw2%z(:,:,0) = stress_field(:,:,1   ,2)
        wall_stress%tauw %z(:,:,0) = stress_field(:,:,1   ,3)
        !$acc end kernels
      end if

      if (is_bound(1, 3) .and. lwm(1, 3) /= 0) then
        !$acc kernels default(present) async(1)
        wall_stress%tauw1%z(:,:,1) = stress_field(:,:,n(3),1)
        wall_stress%tauw2%z(:,:,1) = stress_field(:,:,n(3),2)
        wall_stress%tauw %z(:,:,1) = stress_field(:,:,n(3),3)
        !$acc end kernels
      end if
    end if
    ! Every cfd step, assign_wall_stress_bc, regardless of time_duration_per_action.
    ! The initial stress is used before action_start_time, and the stress computed by 
    ! a wall model is used from action_start_time.
    do idir = 1, 3
      do ibound = 0, 1
        if (is_bound(ibound, idir) .and. lwm(ibound, idir) /= 0) then
          call assign_wall_stress_bc(idir, ibound, visc, wall_stress, bcu, bcv, bcw)
        end if
      end do
    end do
  end subroutine compute_and_assign_wall_stress

  subroutine compute_wall_data(n, is_bound, lwm, l, dl, zc, zf, dzc, dzf, visc, u, v, w, &
                               bcu_mag, bcv_mag, bcw_mag, state, stress, metric)
    implicit none
    integer,  intent(in), dimension(3)      :: n
    logical,  intent(in), dimension(0:1, 3) :: is_bound
    integer,  intent(in), dimension(0:1, 3) :: lwm
    real(rp), intent(in), dimension(3)      :: l, dl
    real(rp), intent(in), dimension(0:)     :: zc, zf, dzc, dzf
    real(rp), intent(in)                    :: visc
    real(rp), intent(in), dimension(0:,0:,0:) :: u, v, w
    type(Bound), intent(in) :: bcu_mag, bcv_mag, bcw_mag
    type(WallState),     intent(inout) :: state
    type(WallStress),    intent(in) :: stress
    type(PerformanceMetric), intent(inout) :: metric
    
    real(rp) :: coef, sgn, u1, u2, v1, v2, w1, w2, u_mag, v_mag, w_mag, uh, vh, wh, this_hwm
    real(rp) :: vel1, vel2, velh, tauw1, tauw2, tauw, tauw1_prev, tauw2_prev, tauw_prev
    real(rp) :: del_v, dveldz, this_hwm_plus, velh_plus, dveldz_plus, utau
    integer  :: i1, i2, j1, j2, k1, k2, i, j, k, ibound, idir, cell_index

    do idir = 1, 3
      do ibound = 0, 1
        if (is_bound(ibound, idir) .and. lwm(ibound, idir) /= 0) then

          select case(idir)
          case(3)
            !$acc parallel loop collapse(2) default(present) async(1) &
            !$acc private(u1,u2,v1,v2,uh,vh,vel1,vel2,velh, &
            !$acc         tauw1,tauw2,tauw,tauw1_prev,tauw2_prev,tauw_prev, &
            !$acc         del_v,dveldz,this_hwm_plus,velh_plus,dveldz_plus,utau, &
            !$acc         cell_index,this_hwm,coef,sgn,k1,k2)
            do j = 1, n(2)
              do i = 1, n(1)
                cell_index = state%hwm_idx%z(i, j, ibound)
                this_hwm = state%hwm%z(i, j, ibound)
                ! The if statement is efficient, as it is based on loop-invariant ibound,
                ! rather than loop indices i and j
                if (ibound == 0) then
                  k2 = cell_index
                  k1 = cell_index - 1
                  coef = (this_hwm - zc(k1)) / dzc(k1)
                  sgn =  1._rp
                else
                  k2 = cell_index
                  k1 = cell_index + 1
                  coef = (this_hwm - (l(3) - zc(k1))) / dzc(k2)
                  sgn = -1._rp
                end if
                u1 = 0.5_rp * (u(i - 1, j, k1) + u(i, j, k1))
                v1 = 0.5_rp * (v(i, j - 1, k1) + v(i, j, k1))
                u2 = 0.5_rp * (u(i - 1, j, k2) + u(i, j, k2))
                v2 = 0.5_rp * (v(i, j - 1, k2) + v(i, j, k2))
                uh = vel_relative(u1, u2, coef, 0._rp)
                vh = vel_relative(v1, v2, coef, 0._rp)
                vel1 = sqrt(u1**2 + v1**2)
                vel2 = sqrt(u2**2 + v2**2)
                velh = sqrt(uh**2 + vh**2)
                !
                ! Local spatial average should benefit the performance
                !
                tauw1_prev = stress%tauw1_prev%z(i, j, ibound) ! s_n
                tauw2_prev = stress%tauw2_prev%z(i, j, ibound) ! s_n
                tauw_prev  = stress%tauw_prev %z(i, j, ibound) ! s_n
                tauw1      = stress%tauw1     %z(i, j, ibound) ! s_n+1
                tauw2      = stress%tauw2     %z(i, j, ibound) ! s_n+1
                tauw       = stress%tauw      %z(i, j, ibound) ! s_n+1
                !
                ! Wall units based on tauw (not tauw1), assuming that
                ! vel1, vel2 and tauw are along the same direction,
                ! and that the velocity profile is monotonically increasing.
                ! This is a common assumption in equilibrium wall models.
                ! The current implementation should be more reasonable than only using
                ! the x-direction velocity and shear stress. If using only the x-direction
                ! info, we have to set the y-direction velocity and shear stress to zero.
                ! How about computing dveldz from dudz and dvdz?
                ! 
                utau = sqrt(tauw)
                del_v = visc/utau
                dveldz = sgn * (vel2 - vel1) / (zc(k2) - zc(k1))
                this_hwm_plus = this_hwm / del_v
                velh_plus = velh / utau
                dveldz_plus = dveldz * del_v / utau
                !
                state%vel_x1     %z(i, j, ibound) = uh
                state%vel_x2     %z(i, j, ibound) = vh
                state%hwm_plus   %z(i, j, ibound) = this_hwm_plus
                state%velh_plus  %z(i, j, ibound) = velh_plus
                state%dveldz_plus%z(i, j, ibound) = dveldz_plus
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
    integer :: i, j, idx, ix, jy
    integer :: nx_coarse, ny_coarse

    nx_coarse = n(1)/interval(1)
    ny_coarse = n(2)/interval(2)

    !$acc parallel loop collapse(2) default(present) private(idx,ix,jy) async(1)
    do j = 1, ny_coarse
      do i = 1, nx_coarse
        ix = (i-1)*interval(1) + 1
        jy = (j-1)*interval(2) + 1
        idx = i_point + (i-1) + (j-1)*nx_coarse
        flattened_state%vel_x1     (idx) = state%vel_x1     %z(ix, jy, ibound)
        flattened_state%vel_x2     (idx) = state%vel_x2     %z(ix, jy, ibound)
        flattened_state%hwm        (idx) = state%hwm        %z(ix, jy, ibound)
        flattened_state%hwm_plus   (idx) = state%hwm_plus   %z(ix, jy, ibound)
        flattened_state%velh_plus  (idx) = state%velh_plus  %z(ix, jy, ibound)
        flattened_state%dveldz_plus(idx) = state%dveldz_plus%z(ix, jy, ibound)
        flattened_state%hwm_idx    (idx) = state%hwm_idx    %z(ix, jy, ibound)
        
        flattened_stress%tauw1     (idx) = stress%tauw1     %z(ix, jy, ibound)
        flattened_stress%tauw2     (idx) = stress%tauw2     %z(ix, jy, ibound)
        flattened_stress%tauw      (idx) = stress%tauw      %z(ix, jy, ibound)
        
        flattened_metric%tauw1     (idx) = metric%tauw1     %z(ix, jy, ibound)
        flattened_metric%tauw2     (idx) = metric%tauw2     %z(ix, jy, ibound)
        flattened_metric%tauw      (idx) = metric%tauw      %z(ix, jy, ibound)
        flattened_metric%tauw1_prev(idx) = metric%tauw1_prev%z(ix, jy, ibound)
        flattened_metric%tauw2_prev(idx) = metric%tauw2_prev%z(ix, jy, ibound)
        flattened_metric%tauw_prev (idx) = metric%tauw_prev %z(ix, jy, ibound)
      end do
    end do

  end subroutine coarsen_and_flatten_wall_data

  subroutine wallmodel_loglaw(visc, flattened_state, flattened_stress, flattened_metric)
    implicit none
    real(rp), intent(in) :: visc
    type(FlattenedState), intent(in) :: flattened_state
    type(FlattenedStress), intent(inout) :: flattened_stress
    type(FlattenedMetric), intent(in), optional :: flattened_metric
    real(rp) :: u1, u2, upar, utau, conv, utau_old, f, fp, tauw_tot, tauw1, tauw2, this_hwm
    integer :: n_points, i

    n_points = size(flattened_state%vel_x1)
    !$acc parallel loop collapse(1) default(present) async(1) &
    !$acc private(u1,u2,upar,utau,conv,utau_old,f,fp,tauw_tot,tauw1,tauw2,this_hwm)
    do i = 1, n_points
      this_hwm = flattened_state%hwm(i)
      u1 = flattened_state%vel_x1(i)
      u2 = flattened_state%vel_x2(i)
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

  subroutine wallmodel_laminar(visc, flattened_state, flattened_stress, flattened_metric)
    implicit none
    real(rp), intent(in) :: visc
    type(FlattenedState), intent(in) :: flattened_state
    type(FlattenedStress), intent(inout) :: flattened_stress
    type(FlattenedMetric), intent(in), optional :: flattened_metric
    real(rp) :: u1, u2, upar, umax, del, tauw_tot, tauw1, tauw2, this_hwm
    integer :: n_points, i

    n_points = size(flattened_state%vel_x1)
    del = 1._rp
    !$acc parallel loop collapse(1) default(present) async(1) &
    !$acc private(u1,u2,upar,umax,tauw_tot,tauw1,tauw2,this_hwm)
    do i = 1, n_points
      this_hwm = flattened_state%hwm(i)
      u1 = flattened_state%vel_x1(i)
      u2 = flattened_state%vel_x2(i)
      upar = sqrt(u1**2 + u2**2)
      umax = upar / (this_hwm / del * (2._rp - this_hwm / del))
      tauw_tot = 2._rp / del * umax * visc
      tauw1 = tauw_tot * u1 / (upar + eps)
      tauw2 = tauw_tot * u2 / (upar + eps)
      flattened_stress%tauw1(i) = tauw1
      flattened_stress%tauw2(i) = tauw2
    end do
  end subroutine wallmodel_laminar

  subroutine wallmodel_DRL(visc, flattened_state, flattened_stress, flattened_metric)
    implicit none
    real(rp), intent(in) :: visc
    type(FlattenedState), intent(in) :: flattened_state
    type(FlattenedStress), intent(inout) :: flattened_stress
    type(FlattenedMetric), intent(in), optional :: flattened_metric
    real(rp) :: u1, u2, upar, tauw_tot, tauw1, tauw2, this_action
    integer :: n_points, i
    logical, save :: is_first = .true.

    real(rp), allocatable, dimension(:,:), save :: drl_state
    real(rp), allocatable, dimension(:,:), save :: drl_action
    real(rp), allocatable, dimension(:,:), save :: drl_reward

    n_points = size(flattened_state%vel_x1)

    if (is_first) then
      call init_smartredis_mpi(db_clustered, MPI_COMM_WORLD)
      allocate(drl_state (3, n_points))
      allocate(drl_reward(2, n_points))
      allocate(drl_action(1, n_points))
    end if
    !$acc update async(1) &
    !$acc self(flattened_state %hwm_plus   , &
    !$acc      flattened_state %velh_plus  , &
    !$acc      flattened_state %dveldz_plus, &
    !$acc      flattened_metric%tauw1      , &
    !$acc      flattened_metric%tauw1_prev )
    !$acc wait(1)
    drl_state(1, :) = flattened_state%hwm_plus
    drl_state(2, :) = flattened_state%velh_plus
    drl_state(3, :) = flattened_state%dveldz_plus
    call put_state(trim(adjustl(tag))//".state", shape(drl_state), drl_state)
    !
    if (is_first) then
      is_first = .false.
    else
      drl_reward(1, :) = flattened_metric%tauw1
      drl_reward(2, :) = flattened_metric%tauw1_prev
      call put_reward(trim(adjustl(tag))//".reward", shape(drl_reward), drl_reward)
    end if
    !
    call get_action(trim(adjustl(tag))//".action", shape(drl_action), drl_action)

    flattened_stress%action(:, 1) = drl_action(1, :)
    !$acc update device(flattened_stress%action) async(1)
    !$acc wait(1)
    !
    !$acc parallel loop collapse(1) default(present) async(1) &
    !$acc private(this_action,tauw_tot,u1,u2,upar,tauw1,tauw2)
    do i = 1, n_points
      this_action = flattened_stress%action(i, 1)
      tauw_tot = this_action * flattened_stress%tauw(i) ! tauw_n (s_n) -> tauw_n+1 (s_n+1)
      u1 = flattened_state%vel_x1(i)
      u2 = flattened_state%vel_x2(i)
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
    integer :: k_pos, i, j, idx, ix, jy
    integer :: nx_coarse, ny_coarse
    
    k_pos = merge(1, n(3), ibound == 0)
    nx_coarse = n(1)/interval(1)
    ny_coarse = n(2)/interval(2)
    
    !$acc parallel loop collapse(2) default(present) private(idx,ix,jy) async(1)
    do j = 1, ny_coarse
      do i = 1, nx_coarse
        ix = (i-1)*interval(1) + 1
        jy = (j-1)*interval(2) + 1
        idx = i_point + (i-1) + (j-1)*nx_coarse
        stress_field(ix, jy, k_pos, 1) = flattened_stress%tauw1(idx)
        stress_field(ix, jy, k_pos, 2) = flattened_stress%tauw2(idx)
        stress_field(ix, jy, k_pos, 3) = flattened_stress%tauw (idx)
      end do
    end do
  end subroutine map_stress_to_sparse_grid

  subroutine interpolate_stress_field(stress_field, n, interval, k_pos)
    implicit none
    real(rp), intent(inout) :: stress_field(0:, 0:, 0:, 1:)
    integer, intent(in), dimension(3) :: n, interval
    integer, intent(in) :: k_pos
    integer :: i, j, i0, i1, j0, j1

    !$acc parallel loop collapse(2) default(present) async(1) &
    !$acc private(i0,i1,j0,j1)
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
  
  subroutine assign_wall_stress_bc(idir, ibound, visc, stress, bcu, bcv, bcw)
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
      !$acc kernels default(present) async(1)
      bcu%z(0:nx, 1:ny, ibound) = sgn * visci * 0.5_rp * (stress%tauw1%z(0:nx  , 1:ny  , ibound) + &
                                                          stress%tauw1%z(1:nx+1, 1:ny  , ibound))
      bcv%z(1:nx, 0:ny, ibound) = sgn * visci * 0.5_rp * (stress%tauw2%z(1:nx  , 0:ny  , ibound) + &
                                                          stress%tauw2%z(1:nx  , 1:ny+1, ibound))
      !$acc end kernels
    end select
  end subroutine assign_wall_stress_bc

end module mod_wallmodel
