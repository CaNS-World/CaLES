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
  ! block
  !   logical :: underflowed
  !   call ieee_get_flag(ieee_underflow, underflowed)
  !   if (underflowed) then
  !     print *, "Underflow detected at", i, j, "value = ", outputs_array_3d(i,j,1   ,:)
  !     call ieee_set_flag(ieee_underflow, .false.)
  !   end if
  ! end block
  use mpi
  use mod_precision, only: rp
  use mod_typedef, only: Bound,BoundProfile,BoundInteger
  use mod_smartredis, only: init_smartredis,put_step_type,put_state,put_reward,get_action
  use mod_params, only: kap_log,b_log,eps,tag,db_clustered, &
                        total_time_steps,agent_interval,action_interval
  use mod_bound, only: boundp
  implicit none
  private
  public :: init_wallmodel, wallmodel

  integer, parameter :: WM_LOG = 1
  integer, parameter :: WM_LAM = 2
  integer, parameter :: WM_DRL = 3
  integer, parameter :: MAX_WALL_MODELS = 10
  
  ! the whole wall model can be a class in future refactoring
  ! for DRL models, convert inputs and outputs to arrays, an additional subroutine may be desired
  ! but for other models, just keep them
  ! better to rename input and output as state and action
  type :: WallModelInput
    type(Bound) :: vel1,vel2,vel,hwm,visc
    type(BoundInteger) :: hwm_index
  end type WallModelInput

  type :: WallModelOutput
    type(Bound) :: tauw1,tauw2,tauw
  end type WallModelOutput

  type :: WallModelReward
    type(Bound) :: vel1,vel1_profile_err
    type(BoundProfile) :: vel1_profile
  end type WallModelReward

  abstract interface
    subroutine wallmodel_interface(visc,hwm,inputs_array,outputs_array,rewards_array)
      import :: rp
      real(rp), intent(in) :: visc,hwm
      real(rp), intent(in), dimension(:,:) :: inputs_array
      real(rp), intent(out), dimension(:,:) :: outputs_array
      real(rp), intent(in), dimension(:,:), optional :: rewards_array
    end subroutine wallmodel_interface
  end interface

  type :: WallModelProcedure
    procedure(wallmodel_interface), pointer, nopass :: ptr => null()
  end type WallModelProcedure

  type(WallModelProcedure) :: wallmodel_proc(MAX_WALL_MODELS)

  contains

  subroutine wallmodel(n,nb,is_bound,lwm,l,dl,zc,zf,dzc,dzf,visc,hwm,u,v,w, &
                       cbcsgs,bcu,bcv,bcw,bcs,bcu_mag,bcv_mag,bcw_mag)
    implicit none
    integer, intent(in), dimension(3) :: n
    integer, intent(in), dimension(0:1,3) :: nb
    logical, intent(in), dimension(0:1,3) :: is_bound
    integer, intent(in), dimension(0:1,3) :: lwm
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc,zf,dzc,dzf
    real(rp), intent(in) :: visc,hwm
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    character(len=1), intent(in), dimension(0:1,3) :: cbcsgs
    type(Bound), intent(inout) :: bcu,bcv,bcw
    type(Bound), intent(in) :: bcs
    type(Bound), intent(in) :: bcu_mag,bcv_mag,bcw_mag
    type(WallModelInput), save :: inputs
    type(WallModelOutput), save :: outputs
    type(WallModelReward), save :: rewards
    real(rp), allocatable, dimension(:,:), save :: inputs_array
    real(rp), allocatable, dimension(:,:), save :: rewards_array
    real(rp), allocatable, dimension(:,:), save :: outputs_array
    real(rp), allocatable, dimension(:,:,:,:), save:: outputs_array_3d

    logical, save :: is_first = .true.
    integer, save :: n_points,n_points_x,n_points_y,n_points_z,interval(3)
    integer, dimension(0:1,3), save :: hwm_index
    integer :: mtype,idir,ibound,cell_index
    integer :: i,j,i0,i1,j0,j1,ierr,i_points,i_var,k
    integer, save :: istep

    if(is_first) then
      is_first = .false.
      wallmodel_proc(WM_LOG)%ptr => wallmodel_loglaw
      wallmodel_proc(WM_LAM)%ptr => wallmodel_laminar
      wallmodel_proc(WM_DRL)%ptr => wallmodel_DRL
      allocate(inputs%vel1%x(0:n(2)+1,0:n(3)+1,0:1)); inputs%vel1%x = 0._rp
      allocate(inputs%vel2%x(0:n(2)+1,0:n(3)+1,0:1)); inputs%vel2%x = 0._rp
      allocate(inputs% vel%x(0:n(2)+1,0:n(3)+1,0:1)); inputs% vel%x = 0._rp
      allocate(inputs% hwm%x(0:n(2)+1,0:n(3)+1,0:1)); inputs% hwm%x = 0._rp
      allocate(inputs%visc%x(0:n(2)+1,0:n(3)+1,0:1)); inputs%visc%x = 0._rp
      allocate(inputs%vel1%y(0:n(1)+1,0:n(3)+1,0:1)); inputs%vel1%y = 0._rp
      allocate(inputs%vel2%y(0:n(1)+1,0:n(3)+1,0:1)); inputs%vel2%y = 0._rp
      allocate(inputs% vel%y(0:n(1)+1,0:n(3)+1,0:1)); inputs% vel%y = 0._rp
      allocate(inputs% hwm%y(0:n(1)+1,0:n(3)+1,0:1)); inputs% hwm%y = 0._rp
      allocate(inputs%visc%y(0:n(1)+1,0:n(3)+1,0:1)); inputs%visc%y = 0._rp
      allocate(inputs%vel1%z(0:n(1)+1,0:n(2)+1,0:1)); inputs%vel1%z = 0._rp
      allocate(inputs%vel2%z(0:n(1)+1,0:n(2)+1,0:1)); inputs%vel2%z = 0._rp
      allocate(inputs% vel%z(0:n(1)+1,0:n(2)+1,0:1)); inputs% vel%z = 0._rp
      allocate(inputs% hwm%z(0:n(1)+1,0:n(2)+1,0:1)); inputs% hwm%z = 0._rp
      allocate(inputs%visc%z(0:n(1)+1,0:n(2)+1,0:1)); inputs%visc%z = 0._rp

      allocate(inputs%hwm_index%x(0:n(2)+1,0:n(3)+1,0:1)); inputs%hwm_index%x = 0
      allocate(inputs%hwm_index%y(0:n(1)+1,0:n(3)+1,0:1)); inputs%hwm_index%y = 0
      allocate(inputs%hwm_index%z(0:n(1)+1,0:n(2)+1,0:1)); inputs%hwm_index%z = 0

      allocate(rewards%vel1%x(0:n(2)+1,0:n(3)+1,0:1)); rewards%vel1%x = 0._rp
      allocate(rewards%vel1%y(0:n(1)+1,0:n(3)+1,0:1)); rewards%vel1%y = 0._rp
      allocate(rewards%vel1%z(0:n(1)+1,0:n(2)+1,0:1)); rewards%vel1%z = 0._rp

      allocate(rewards%vel1_profile_err%x(0:n(2)+1,0:n(3)+1,0:1)); rewards%vel1_profile_err%x = 0._rp
      allocate(rewards%vel1_profile_err%y(0:n(1)+1,0:n(3)+1,0:1)); rewards%vel1_profile_err%y = 0._rp
      allocate(rewards%vel1_profile_err%z(0:n(1)+1,0:n(2)+1,0:1)); rewards%vel1_profile_err%z = 0._rp

      allocate(rewards%vel1_profile%x(1:n(1),0:n(2)+1,0:n(3)+1,0:1)); rewards%vel1_profile%x = 0._rp
      allocate(rewards%vel1_profile%y(1:n(2),0:n(1)+1,0:n(3)+1,0:1)); rewards%vel1_profile%y = 0._rp
      allocate(rewards%vel1_profile%z(1:n(3),0:n(1)+1,0:n(2)+1,0:1)); rewards%vel1_profile%z = 0._rp

      allocate(outputs%tauw1%x(0:n(2)+1,0:n(3)+1,0:1)); outputs%tauw1%x = 0._rp
      allocate(outputs%tauw2%x(0:n(2)+1,0:n(3)+1,0:1)); outputs%tauw2%x = 0._rp
      allocate(outputs% tauw%x(0:n(2)+1,0:n(3)+1,0:1)); outputs% tauw%x = 0._rp
      allocate(outputs%tauw1%y(0:n(1)+1,0:n(3)+1,0:1)); outputs%tauw1%y = 0._rp
      allocate(outputs%tauw2%y(0:n(1)+1,0:n(3)+1,0:1)); outputs%tauw2%y = 0._rp
      allocate(outputs% tauw%y(0:n(1)+1,0:n(3)+1,0:1)); outputs% tauw%y = 0._rp
      allocate(outputs%tauw1%z(0:n(1)+1,0:n(2)+1,0:1)); outputs%tauw1%z = 0._rp
      allocate(outputs%tauw2%z(0:n(1)+1,0:n(2)+1,0:1)); outputs%tauw2%z = 0._rp
      allocate(outputs% tauw%z(0:n(1)+1,0:n(2)+1,0:1)); outputs% tauw%z = 0._rp

      n_points = 0
      interval(1:3) = agent_interval ! agent_interval is a single value
      n_points_x = (n(2)/interval(2))*(n(3)/interval(3)) ! only internal points for inputs_array and outputs_array
      n_points_y = (n(1)/interval(1))*(n(3)/interval(3))
      n_points_z = (n(1)/interval(1))*(n(2)/interval(2))
      if(is_bound(0,1).and.lwm(0,1) /= 0) n_points = n_points + n_points_x
      if(is_bound(1,1).and.lwm(1,1) /= 0) n_points = n_points + n_points_x
      if(is_bound(0,2).and.lwm(0,2) /= 0) n_points = n_points + n_points_y
      if(is_bound(1,2).and.lwm(1,2) /= 0) n_points = n_points + n_points_y
      if(is_bound(0,3).and.lwm(0,3) /= 0) n_points = n_points + n_points_z
      if(is_bound(1,3).and.lwm(1,3) /= 0) n_points = n_points + n_points_z
      allocate( inputs_array(n_points,5     ));  inputs_array = 0._rp
      allocate(rewards_array(n_points,3+n(3))); rewards_array = 0._rp
      allocate(outputs_array(n_points,3     )); outputs_array = 0._rp
      allocate(outputs_array_3d(0:n(1)+1,0:n(2)+1,0:n(3)+1,3)); outputs_array_3d = 0._rp
      istep = 0
      call init_wallmodel(n,is_bound,lwm,l,dl,zc,hwm,hwm_index,inputs)
    else 
      istep = istep + 1
    end if

    ! in the future, a separate subroutine for handling the bc's of 2D arrays may be desired
    ! compute at cell centers u, v, ...
    ! send to drl to compute tauw, stored in 2D arrays
    ! store 2D arrays in 3D arrays
    ! implement boundary conditions for tauw, only periodic bc's are used (repitition can be avoided)
    ! store 3D arrays to 2D arrays
    ! impose u and v bc's
    call wallmodel_input(n,is_bound,lwm,l,dl,zc,zf,dzc,dzf,visc,hwm,hwm_index,u,v,w, &
                         bcu_mag,bcv_mag,bcw_mag,inputs,rewards)

    if(mod(istep,action_interval) == 0) then  ! move to main file, but it has to be here for now

      i_points = 1
      ! if(is_bound(0,1).and.lwm(0,1) /= 0) then
      !   inputs_array (i_points:i_points+n_points_x-1,1,1) = reshape(inputs%vel1_p1%x(1:n(2),1:n(3),0),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,1,2) = reshape(inputs%vel2_p1%x(1:n(2),1:n(3),0),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,1,3) = reshape(inputs% vel_p1%x(1:n(2),1:n(3),0),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,1,4) = reshape(inputs% hwm_p1%x(1:n(2),1:n(3),0),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,1,5) = reshape(inputs%visc_p1%x(1:n(2),1:n(3),0),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,2,1) = reshape(inputs%vel1_p2%x(1:n(2),1:n(3),0),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,2,2) = reshape(inputs%vel2_p2%x(1:n(2),1:n(3),0),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,2,3) = reshape(inputs% vel_p2%x(1:n(2),1:n(3),0),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,2,4) = reshape(inputs% hwm_p2%x(1:n(2),1:n(3),0),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,2,5) = reshape(inputs%visc_p2%x(1:n(2),1:n(3),0),(/n_points_x/))
      !   rewards_array(i_points:i_points+n_points_x-1,1,1) = reshape(rewards%vel_p1%x(1:n(2),1:n(3),0),(/n_points_x/))
      !   rewards_array(i_points:i_points+n_points_x-1,2,1) = reshape(rewards%vel_p2%x(1:n(2),1:n(3),0),(/n_points_x/))
      !   i_points = i_points + n_points_x
      ! end if
      ! if(is_bound(1,1).and.lwm(1,1) /= 0) then
      !   inputs_array (i_points:i_points+n_points_x-1,1,1) = reshape(inputs%vel1_p1%x(1:n(2),1:n(3),1),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,1,2) = reshape(inputs%vel2_p1%x(1:n(2),1:n(3),1),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,1,3) = reshape(inputs% vel_p1%x(1:n(2),1:n(3),1),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,1,4) = reshape(inputs% hwm_p1%x(1:n(2),1:n(3),1),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,1,5) = reshape(inputs%visc_p1%x(1:n(2),1:n(3),1),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,2,1) = reshape(inputs%vel1_p2%x(1:n(2),1:n(3),1),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,2,2) = reshape(inputs%vel2_p2%x(1:n(2),1:n(3),1),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,2,3) = reshape(inputs% vel_p2%x(1:n(2),1:n(3),1),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,2,4) = reshape(inputs% hwm_p2%x(1:n(2),1:n(3),1),(/n_points_x/))
      !   inputs_array (i_points:i_points+n_points_x-1,2,5) = reshape(inputs%visc_p2%x(1:n(2),1:n(3),1),(/n_points_x/))
      !   rewards_array(i_points:i_points+n_points_x-1,1,1) = reshape(rewards%vel_p1%x(1:n(2),1:n(3),1),(/n_points_x/))
      !   rewards_array(i_points:i_points+n_points_x-1,2,1) = reshape(rewards%vel_p2%x(1:n(2),1:n(3),1),(/n_points_x/))
      !   i_points = i_points + n_points_x
      ! end if
      ! if(is_bound(0,2).and.lwm(0,2) /= 0) then
      !   inputs_array (i_points:i_points+n_points_y-1,1,1) = reshape(inputs%vel1_p1%y(1:n(1),1:n(3),0),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,1,2) = reshape(inputs%vel2_p1%y(1:n(1),1:n(3),0),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,1,3) = reshape(inputs% vel_p1%y(1:n(1),1:n(3),0),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,1,4) = reshape(inputs% hwm_p1%y(1:n(1),1:n(3),0),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,1,5) = reshape(inputs%visc_p1%y(1:n(1),1:n(3),0),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,2,1) = reshape(inputs%vel1_p2%y(1:n(1),1:n(3),0),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,2,2) = reshape(inputs%vel2_p2%y(1:n(1),1:n(3),0),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,2,3) = reshape(inputs% vel_p2%y(1:n(1),1:n(3),0),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,2,4) = reshape(inputs% hwm_p2%y(1:n(1),1:n(3),0),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,2,5) = reshape(inputs%visc_p2%y(1:n(1),1:n(3),0),(/n_points_y/))
      !   rewards_array(i_points:i_points+n_points_y-1,1,1) = reshape(rewards%vel_p1%y(1:n(1),1:n(3),0),(/n_points_y/))
      !   rewards_array(i_points:i_points+n_points_y-1,2,1) = reshape(rewards%vel_p2%y(1:n(1),1:n(3),0),(/n_points_y/))
      !   i_points = i_points + n_points_y
      ! end if
      ! if(is_bound(1,2).and.lwm(1,2) /= 0) then
      !   inputs_array (i_points:i_points+n_points_y-1,1,1) = reshape(inputs%vel1_p1%y(1:n(1),1:n(3),1),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,1,2) = reshape(inputs%vel2_p1%y(1:n(1),1:n(3),1),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,1,3) = reshape(inputs% vel_p1%y(1:n(1),1:n(3),1),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,1,4) = reshape(inputs% hwm_p1%y(1:n(1),1:n(3),1),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,1,5) = reshape(inputs%visc_p1%y(1:n(1),1:n(3),1),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,2,1) = reshape(inputs%vel1_p2%y(1:n(1),1:n(3),1),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,2,2) = reshape(inputs%vel2_p2%y(1:n(1),1:n(3),1),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,2,3) = reshape(inputs% vel_p2%y(1:n(1),1:n(3),1),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,2,4) = reshape(inputs% hwm_p2%y(1:n(1),1:n(3),1),(/n_points_y/))
      !   inputs_array (i_points:i_points+n_points_y-1,2,5) = reshape(inputs%visc_p2%y(1:n(1),1:n(3),1),(/n_points_y/))
      !   rewards_array(i_points:i_points+n_points_y-1,1,1) = reshape(rewards%vel_p1%y(1:n(1),1:n(3),1),(/n_points_y/))
      !   rewards_array(i_points:i_points+n_points_y-1,2,1) = reshape(rewards%vel_p2%y(1:n(1),1:n(3),1),(/n_points_y/))
      !   i_points = i_points + n_points_y
      ! end if
      if(is_bound(0,3).and.lwm(0,3) /= 0) then
        inputs_array (i_points:i_points+n_points_z-1,1) = reshape(inputs %vel1%z(1:n(1):interval(1),1:n(2):interval(2),0),(/n_points_z/))
        inputs_array (i_points:i_points+n_points_z-1,2) = reshape(inputs %vel2%z(1:n(1):interval(1),1:n(2):interval(2),0),(/n_points_z/))
        inputs_array (i_points:i_points+n_points_z-1,3) = reshape(inputs % vel%z(1:n(1):interval(1),1:n(2):interval(2),0),(/n_points_z/))
        inputs_array (i_points:i_points+n_points_z-1,4) = reshape(inputs % hwm%z(1:n(1):interval(1),1:n(2):interval(2),0),(/n_points_z/))
        inputs_array (i_points:i_points+n_points_z-1,5) = reshape(inputs %visc%z(1:n(1):interval(1),1:n(2):interval(2),0),(/n_points_z/))
        rewards_array(i_points:i_points+n_points_z-1,1) = reshape(rewards%vel1%z(1:n(1):interval(1),1:n(2):interval(2),0),(/n_points_z/)) ! u
        rewards_array(i_points:i_points+n_points_z-1,2) = outputs_array(i_points:i_points+n_points_z-1,1) ! taux
        rewards_array(i_points:i_points+n_points_z-1,3) = reshape(rewards%vel1_profile_err%z(1:n(1):interval(1),1:n(2):interval(2),0),(/n_points_z/)) ! u_profile_err
        do k = 1,n(3)
          rewards_array(i_points:i_points+n_points_z-1,3+k) = reshape(rewards%vel1_profile%z(k,1:n(1):interval(1),1:n(2):interval(2),0),(/n_points_z/)) ! u_profile
        end do
        i_points = i_points + n_points_z
      end if
      if(is_bound(1,3).and.lwm(1,3) /= 0) then
        inputs_array (i_points:i_points+n_points_z-1,1) = reshape(inputs %vel1%z(1:n(1):interval(1),1:n(2):interval(2),1),(/n_points_z/))
        inputs_array (i_points:i_points+n_points_z-1,2) = reshape(inputs %vel2%z(1:n(1):interval(1),1:n(2):interval(2),1),(/n_points_z/))
        inputs_array (i_points:i_points+n_points_z-1,3) = reshape(inputs % vel%z(1:n(1):interval(1),1:n(2):interval(2),1),(/n_points_z/))
        inputs_array (i_points:i_points+n_points_z-1,4) = reshape(inputs % hwm%z(1:n(1):interval(1),1:n(2):interval(2),1),(/n_points_z/))
        inputs_array (i_points:i_points+n_points_z-1,5) = reshape(inputs %visc%z(1:n(1):interval(1),1:n(2):interval(2),1),(/n_points_z/))
        rewards_array(i_points:i_points+n_points_z-1,1) = reshape(rewards%vel1%z(1:n(1):interval(1),1:n(2):interval(2),1),(/n_points_z/))
        rewards_array(i_points:i_points+n_points_z-1,2) = outputs_array(i_points:i_points+n_points_z-1,1)
        rewards_array(i_points:i_points+n_points_z-1,3) = reshape(rewards%vel1_profile_err%z(1:n(1):interval(1),1:n(2):interval(2),1),(/n_points_z/))
        do k = 1,n(3)
          rewards_array(i_points:i_points+n_points_z-1,3+k) = reshape(rewards%vel1_profile%z(k,1:n(1):interval(1),1:n(2):interval(2),1),(/n_points_z/))
        end do
        i_points = i_points + n_points_z
      end if

      if(i_points /= n_points + 1) then
        print *, 'ERROR: i_points /= n_points + 1.'
      end if

      mtype = maxval(lwm(0:1,1:3)) ! only allow a single type of wall model for now
      call wallmodel_proc(mtype)%ptr(visc,hwm,inputs_array,outputs_array,rewards_array)

      i_points = 1
      ! if(is_bound(0,1).and.lwm(0,1) /= 0) then
      !   outputs%tauw1_p1%x(1:n(2),1:n(3),0) = reshape(outputs_array(i_points:i_points+n_points_x-1,1,1),(/n(2),n(3)/))
      !   outputs%tauw2_p1%x(1:n(2),1:n(3),0) = reshape(outputs_array(i_points:i_points+n_points_x-1,1,2),(/n(2),n(3)/))
      !   outputs%tauw1_p2%x(1:n(2),1:n(3),0) = reshape(outputs_array(i_points:i_points+n_points_x-1,2,1),(/n(2),n(3)/))
      !   outputs%tauw2_p2%x(1:n(2),1:n(3),0) = reshape(outputs_array(i_points:i_points+n_points_x-1,2,2),(/n(2),n(3)/))
      !   i_points = i_points + n_points_x
      ! end if
      ! if(is_bound(1,1).and.lwm(1,1) /= 0) then
      !   outputs%tauw1_p1%x(1:n(2),1:n(3),1) = reshape(outputs_array(i_points:i_points+n_points_x-1,1,1),(/n(2),n(3)/))
      !   outputs%tauw2_p1%x(1:n(2),1:n(3),1) = reshape(outputs_array(i_points:i_points+n_points_x-1,1,2),(/n(2),n(3)/))
      !   outputs%tauw1_p2%x(1:n(2),1:n(3),1) = reshape(outputs_array(i_points:i_points+n_points_x-1,2,1),(/n(2),n(3)/))
      !   outputs%tauw2_p2%x(1:n(2),1:n(3),1) = reshape(outputs_array(i_points:i_points+n_points_x-1,2,2),(/n(2),n(3)/))
      !   i_points = i_points + n_points_x
      ! end if
      ! if(is_bound(0,2).and.lwm(0,2) /= 0) then
      !   outputs%tauw1_p1%y(1:n(1),1:n(3),0) = reshape(outputs_array(i_points:i_points+n_points_y-1,1,1),(/n(1),n(3)/))
      !   outputs%tauw2_p1%y(1:n(1),1:n(3),0) = reshape(outputs_array(i_points:i_points+n_points_y-1,1,2),(/n(1),n(3)/))
      !   outputs%tauw1_p2%y(1:n(1),1:n(3),0) = reshape(outputs_array(i_points:i_points+n_points_y-1,2,1),(/n(1),n(3)/))
      !   outputs%tauw2_p2%y(1:n(1),1:n(3),0) = reshape(outputs_array(i_points:i_points+n_points_y-1,2,2),(/n(1),n(3)/))
      !   i_points = i_points + n_points_y
      ! end if
      ! if(is_bound(1,2).and.lwm(1,2) /= 0) then
      !   outputs%tauw1_p1%y(1:n(1),1:n(3),1) = reshape(outputs_array(i_points:i_points+n_points_y-1,1,1),(/n(1),n(3)/))
      !   outputs%tauw2_p1%y(1:n(1),1:n(3),1) = reshape(outputs_array(i_points:i_points+n_points_y-1,1,2),(/n(1),n(3)/))
      !   outputs%tauw1_p2%y(1:n(1),1:n(3),1) = reshape(outputs_array(i_points:i_points+n_points_y-1,2,1),(/n(1),n(3)/))
      !   outputs%tauw2_p2%y(1:n(1),1:n(3),1) = reshape(outputs_array(i_points:i_points+n_points_y-1,2,2),(/n(1),n(3)/))
      !   i_points = i_points + n_points_y
      ! end if
      if(is_bound(0,3).and.lwm(0,3) /= 0) then
        outputs_array_3d(1:n(1):interval(1),1:n(2):interval(2),1   ,:) = reshape(outputs_array(i_points:i_points+n_points_z-1,:),(/n(1)/interval(1),n(2)/interval(2),3/))
        i_points = i_points + n_points_z
      end if
      if(is_bound(1,3).and.lwm(1,3) /= 0) then
        outputs_array_3d(1:n(1):interval(1),1:n(2):interval(2),n(3),:) = reshape(outputs_array(i_points:i_points+n_points_z-1,:),(/n(1)/interval(1),n(2)/interval(2),3/))
        i_points = i_points + n_points_z
      end if
      ! n + 1 filled, but 0 not filled
      ! cbcsgs and bcs must be assigned. bcs is zero on the walls, so the cells have opposite values
      ! on both sides of the wall, which brings zero wall shear stress
      ! square duct has not been tested
      do i_var = 1,3
        call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,outputs_array_3d(:,:,:,i_var))
      end do

      if(is_bound(0,3).and.lwm(0,3) /= 0) then
        ! bilinear interpolation, uniform grid
        ! the interpolated shear stress is smoother than the original one
        do j = 1,n(2)
          do i = 1,n(1)
            i0 = ((i-1)/interval(1))*interval(1) + 1
            i1 = i0 + interval(1)
            j0 = ((j-1)/interval(2))*interval(2) + 1
            j1 = j0 + interval(2)
            outputs_array_3d(i,j,1   ,:) = (outputs_array_3d(i0,j0,1   ,:)*(i1-i)*(j1-j) + &
                                            outputs_array_3d(i1,j0,1   ,:)*(i-i0)*(j1-j) + &
                                            outputs_array_3d(i0,j1,1   ,:)*(i1-i)*(j-j0) + &
                                            outputs_array_3d(i1,j1,1   ,:)*(i-i0)*(j-j0))/((i1-i0)*(j1-j0))
          end do
        end do
      end if

      if(is_bound(1,3).and.lwm(1,3) /= 0) then
        ! bilinear interpolation, uniform grid
        do j = 1,n(2)
          do i = 1,n(1)
            i0 = ((i-1)/interval(1))*interval(1) + 1
            i1 = i0 + interval(1)
            j0 = ((j-1)/interval(2))*interval(2) + 1
            j1 = j0 + interval(2)
            outputs_array_3d(i,j,n(3),:) = (outputs_array_3d(i0,j0,n(3),:)*(i1-i)*(j1-j) + &
                                            outputs_array_3d(i1,j0,n(3),:)*(i-i0)*(j1-j) + &
                                            outputs_array_3d(i0,j1,n(3),:)*(i1-i)*(j-j0) + &
                                            outputs_array_3d(i1,j1,n(3),:)*(i-i0)*(j-j0))/((i1-i0)*(j1-j0))
          end do
        end do
      end if

      ! 0 and other empty ghost cells filled
      do i_var = 1,3
        call boundp(cbcsgs,n,bcs,nb,is_bound,dl,dzc,outputs_array_3d(:,:,:,i_var))
      end do
      
      if(is_bound(0,3).and.lwm(0,3) /= 0) then
        outputs%tauw1%z(:,:,0) = outputs_array_3d(:,:,1   ,1)
        outputs%tauw2%z(:,:,0) = outputs_array_3d(:,:,1   ,2)
      end if

      if(is_bound(1,3).and.lwm(1,3) /= 0) then
        outputs%tauw1%z(:,:,1) = outputs_array_3d(:,:,n(3),1)
        outputs%tauw2%z(:,:,1) = outputs_array_3d(:,:,n(3),2)
      end if

      do idir = 1,3
        do ibound = 0,1
          if(is_bound(ibound,idir).and.lwm(ibound,idir) /= 0) then
            call assign_wallmodelbc(idir,ibound,visc,outputs,bcu,bcv,bcw)
          end if
        end do
      end do
    end if
  end subroutine wallmodel

  ! find the cell_index required for interpolation to the wall model height.
  ! The stored cell_index corresponds to the cells far from a wall, i.e., i2,j2,k2.
  ! Remmeber to set hwm strightly higher than the first cell center, and lower
  ! than the last cell center (hwm=hwm-eps)

  subroutine init_wallmodel(n,is_bound,lwm,l,dl,zc,hwm,hwm_index,inputs)
    implicit none
    integer, intent(in) :: n(3)
    logical, intent(in) :: is_bound(0:1,3)
    integer, intent(in) :: lwm(0:1,3)
    real(rp), intent(in) :: l(3),dl(3),zc(0:),hwm
    integer, intent(out) :: hwm_index(0:1,3)
    type(WallModelInput), intent(inout) :: inputs
    integer, allocatable, dimension(:) :: seed
    integer :: i,j,k,i1,i2,j1,j2,k1,k2,ncells

    if(is_bound(0,1).and.lwm(0,1) /= 0) then ! to remove statement
      i = 1
      do while((i-0.5)*dl(1) < hwm)
        i = i + 1
      end do
      i2 = i
      i1 = i - 1
      hwm_index(0,1) = i2
    end if
    if(is_bound(1,1).and.lwm(1,1) /= 0) then
      i = n(1)
      do while((n(1)-i+0.5)*dl(1) < hwm)
        i = i - 1
      end do
      i2 = i
      i1 = i + 1
      hwm_index(1,1) = i2
    end if
    if(is_bound(0,2).and.lwm(0,2) /= 0) then
      j = 1
      do while((j-0.5)*dl(2) < hwm)
        j = j + 1
      end do
      j2 = j
      j1 = j - 1
      hwm_index(0,2) = j2
    end if
    if(is_bound(1,2).and.lwm(1,2) /= 0) then
      j = n(2)
      do while((n(2)-j+0.5)*dl(2) < hwm)
        j = j - 1
      end do
      j2 = j
      j1 = j + 1
      hwm_index(1,2) = j2
    end if
    if(is_bound(0,3).and.lwm(0,3) /= 0) then
      ncells = n(1)*n(2)
      call random_seed(size = ncells)
      allocate(seed(ncells))
      seed = 123
      call random_seed(put = seed)
      call random_number(inputs%hwm%z(1:n(1),1:n(2),0))
      inputs%hwm%z(1:n(1),1:n(2),0) = 0.075_rp + (0.150_rp-0.075_rp)*inputs%hwm%z(1:n(1),1:n(2),0)
      !
      do j = 1,n(2)
        do i = 1,n(1)
          k = 1
          do while(zc(k) < inputs%hwm%z(i,j,0))
            k = k + 1
          end do
          k2 = k
          k1 = k - 1
          inputs%hwm_index%z(i,j,0) = k2
        end do
      end do
    end if
    if(is_bound(1,3).and.lwm(1,3) /= 0) then
      ncells = n(1)*n(2)
      call random_seed(size = ncells)
      if(allocated(seed)) deallocate(seed)
      allocate(seed(ncells))
      seed = 123 ! bottom and top wall model heights are the same for now
      call random_seed(put = seed)
      call random_number(inputs%hwm%z(1:n(1),1:n(2),1))
      inputs%hwm%z(1:n(1),1:n(2),1) = 0.075_rp + (0.150_rp-0.075_rp)*inputs%hwm%z(1:n(1),1:n(2),1)
      !
      do j = 1,n(2)
        do i = 1,n(1)
          k = n(3)
          do while(l(3)-zc(k) < inputs%hwm%z(i,j,1))
            k = k - 1
          end do
          k2 = k
          k1 = k + 1
          inputs%hwm_index%z(i,j,1) = k2
        end do
      end do
    end if
  end subroutine init_wallmodel

  subroutine wallmodel_input(n,is_bound,lwm,l,dl,zc,zf,dzc,dzf,visc,hwm,hwm_index,u,v,w, &
                             bcu_mag,bcv_mag,bcw_mag,inputs,rewards)
    implicit none
    integer, intent(in) :: n(3)
    logical, intent(in) :: is_bound(0:1,3)
    integer, intent(in) :: lwm(0:1,3)
    real(rp), intent(in) :: l(3),dl(3)
    real(rp), intent(in), dimension(0:) :: zc,zf,dzc,dzf
    real(rp), intent(in) :: visc,hwm
    integer, intent(in) :: hwm_index(0:1,3)
    real(rp), intent(in), dimension(0:,0:,0:) :: u,v,w
    type(Bound), intent(in) :: bcu_mag,bcv_mag,bcw_mag
    type(WallModelInput), intent(inout) :: inputs
    type(WallModelReward), intent(inout) :: rewards
    real(rp) :: coef,wei,u1,u2,v1,v2,w1,w2,u_mag,v_mag,w_mag,uh,vh,wh,this_hwm
    integer :: i1,i2,j1,j2,k1,k2,i,j,k,ibound,idir,cell_index
    logical, save :: is_first = .true.
    integer, save :: istep,num_samples
    real(rp), allocatable, dimension(:), save :: u_ref,u_ref_0,u_ref_1,u_profile
    real(rp), allocatable, dimension(:,:,:,:), save :: u_profile_ave
    real(rp) :: dummy
    integer :: ierr
    
    if(is_first) then
      is_first = .false.
      istep = 0
      num_samples = 1
      rewards%vel1%z = 0._rp
      allocate(u_ref_0(n(3)))
      allocate(u_ref_1(n(3)))
      ! open(123, file='stats-single-point-chan-05200.out', status='old', action='read')
      ! do i = 1, n(3) ! we assume the z direction is split into 2 parts, only one part is loaded
      !   read(123,*) dummy, dummy, u_ref_0(i)
      ! end do
      ! u_ref_1 = u_ref_0(n(3):1:-1)
      ! close(123)
      allocate(u_profile_ave(n(3),n(1),n(2),0:1)) ! only support case 3
      u_profile_ave = 0._rp
      allocate(u_ref(n(3)))
      allocate(u_profile(n(3)))
      !
    else
      istep = istep + 1
      num_samples = num_samples + 1
    end if
    if(mod(istep,action_interval) == 1) then
      num_samples = 1
      rewards%vel1%z = 0._rp
      u_profile_ave = 0._rp
    end if

    do idir = 1,3
      do ibound = 0,1
        if(is_bound(ibound,idir).and.lwm(ibound,idir) /= 0) then

          select case(idir)
          ! case(1)  ! x-direction walls: vel1 = v, vel2 = w
          !   if(ibound == 0) then
          !     i2 = cell_index
          !     i1 = cell_index - 1
          !     coef = (hwm-(i1-0.5_rp)*dl(1))/dl(1)
          !   else
          !     i2 = cell_index
          !     i1 = cell_index + 1
          !     coef = (hwm-(n(1)-i1+0.5_rp)*dl(1))/dl(1)
          !   end if
          !   do k = 1,n(3)
          !     ! do j = 0,n(2)
          !     do j = 1,n(2)
          !       v1 = v(i1,j,k)
          !       v2 = v(i2,j,k)
          !       w1 = 0.25_rp*(w(i1,j,k) + w(i1,j+1,k) + w(i1,j,k-1) + w(i1,j+1,k-1))
          !       w2 = 0.25_rp*(w(i2,j,k) + w(i2,j+1,k) + w(i2,j,k-1) + w(i2,j+1,k-1))
          !       v_mag = bcv_mag%x(j,k,ibound)
          !       w_mag = 0.25_rp*(bcw_mag%x(j,k  ,ibound) + bcw_mag%x(j+1,k  ,ibound) + &
          !                        bcw_mag%x(j,k-1,ibound) + bcw_mag%x(j+1,k-1,ibound))
          !       vh = vel_relative(v1,v2,coef,v_mag)
          !       wh = vel_relative(w1,w2,coef,w_mag)
          !       inputs%vel1_p1%x(j,k,ibound) = vh
          !       inputs%vel2_p1%x(j,k,ibound) = wh
          !       inputs% vel_p1%x(j,k,ibound) = sqrt(vh**2 + wh**2)
          !       inputs% hwm_p1%x(j,k,ibound) = hwm
          !       inputs%visc_p1%x(j,k,ibound) = visc
          !       rewards%vel_p1%x(j,k,ibound) = -abs(sqrt(vh**2 + wh**2) - 0.853_rp)
          !     end do
          !   end do
          !   ! do k = 0,n(3)
          !   do k = 1,n(3)
          !     do j = 1,n(2)
          !       wei = (zf(k)-zc(k))/dzc(k)
          !       v1 = 0.5_rp*((1._rp-wei)*(v(i1,j-1,k  ) + v(i1,j,k  )) + &
          !                           wei *(v(i1,j-1,k+1) + v(i1,j,k+1)))
          !       v2 = 0.5_rp*((1._rp-wei)*(v(i2,j-1,k  ) + v(i2,j,k  )) + &
          !                           wei *(v(i2,j-1,k+1) + v(i2,j,k+1)))
          !       w1 = w(i1,j,k)
          !       w2 = w(i2,j,k)
          !       v_mag = 0.5_rp*((1._rp-wei)*(bcv_mag%x(j-1,k  ,ibound) + bcv_mag%x(j,k  ,ibound)) + &
          !                              wei *(bcv_mag%x(j-1,k+1,ibound) + bcv_mag%x(j,k+1,ibound)))
          !       w_mag = bcw_mag%x(j,k,ibound)
          !       vh = vel_relative(v1,v2,coef,v_mag)
          !       wh = vel_relative(w1,w2,coef,w_mag)
          !       inputs%vel1_p2%x(j,k,ibound) = vh
          !       inputs%vel2_p2%x(j,k,ibound) = wh
          !       inputs% vel_p2%x(j,k,ibound) = sqrt(vh**2 + wh**2)
          !       inputs% hwm_p2%x(j,k,ibound) = hwm
          !       inputs%visc_p2%x(j,k,ibound) = visc
          !       rewards%vel_p2%x(j,k,ibound) = -abs(sqrt(vh**2 + wh**2) - 0.853_rp)
          !     end do
          !   end do
          ! case(2)  ! y-direction walls: vel1 = u, vel2 = w
          !   if(ibound == 0) then
          !     j2 = cell_index
          !     j1 = cell_index - 1
          !     coef = (hwm-(j1-0.5_rp)*dl(2))/dl(2)
          !   else
          !     j2 = cell_index
          !     j1 = cell_index + 1
          !     coef = (hwm-(n(2)-j1+0.5_rp)*dl(2))/dl(2)
          !   end if
          !   do k = 1,n(3)
          !     ! do i = 0,n(1)
          !     do i = 1,n(1)
          !       u1 = u(i,j1,k)
          !       u2 = u(i,j2,k)
          !       w1 = 0.25_rp*(w(i,j1,k) + w(i+1,j1,k) + w(i,j1,k-1) + w(i+1,j1,k-1))
          !       w2 = 0.25_rp*(w(i,j2,k) + w(i+1,j2,k) + w(i,j2,k-1) + w(i+1,j2,k-1))
          !       u_mag = bcu_mag%y(i,k,ibound)
          !       w_mag = 0.25_rp*(bcw_mag%y(i,k  ,ibound) + bcw_mag%y(i+1,k  ,ibound) + &
          !                        bcw_mag%y(i,k-1,ibound) + bcw_mag%y(i+1,k-1,ibound))
          !       uh = vel_relative(u1,u2,coef,u_mag)
          !       wh = vel_relative(w1,w2,coef,w_mag)
          !       inputs%vel1_p1%y(i,k,ibound) = uh
          !       inputs%vel2_p1%y(i,k,ibound) = wh
          !       inputs% vel_p1%y(i,k,ibound) = sqrt(uh**2 + wh**2)
          !       inputs% hwm_p1%y(i,k,ibound) = hwm
          !       inputs%visc_p1%y(i,k,ibound) = visc
          !       rewards%vel_p1%y(i,k,ibound) = -abs(sqrt(uh**2 + wh**2) - 0.853_rp)
          !     end do
          !   end do
          !   ! do k = 0,n(3)
          !   do k = 1,n(3)
          !     do i = 1,n(1)
          !       wei = (zf(k)-zc(k))/dzc(k)
          !       u1 = 0.5_rp*((1._rp-wei)*(u(i-1,j1,k  ) + u(i,j1,k  )) + &
          !                           wei *(u(i-1,j1,k+1) + u(i,j1,k+1)))
          !       u2 = 0.5_rp*((1._rp-wei)*(u(i-1,j2,k  ) + u(i,j2,k  )) + &
          !                           wei *(u(i-1,j2,k+1) + u(i,j2,k+1)))
          !       w1 = w(i,j1,k)
          !       w2 = w(i,j2,k)
          !       u_mag = 0.5_rp*((1._rp-wei)*(bcu_mag%y(i-1,k  ,ibound) + bcu_mag%y(i,k  ,ibound)) + &
          !                              wei *(bcu_mag%y(i-1,k+1,ibound) + bcu_mag%y(i,k+1,ibound)))
          !       w_mag = bcw_mag%y(i,k,ibound)
          !       uh = vel_relative(u1,u2,coef,u_mag)
          !       wh = vel_relative(w1,w2,coef,w_mag)
          !       inputs%vel1_p2%y(i,k,ibound) = uh
          !       inputs%vel2_p2%y(i,k,ibound) = wh
          !       inputs% vel_p2%y(i,k,ibound) = sqrt(uh**2 + wh**2)
          !       inputs% hwm_p2%y(i,k,ibound) = hwm
          !       inputs%visc_p2%y(i,k,ibound) = visc
          !       rewards%vel_p2%y(i,k,ibound) = -abs(sqrt(uh**2 + wh**2) - 0.853_rp)
          !     end do
          !   end do
          ! case(3)  ! z-direction walls: vel1 = u, vel2 = v
          !   if(ibound == 0) then
          !     k2 = cell_index
          !     k1 = cell_index - 1
          !     coef = (hwm-zc(k1))/dzc(k1)
          !   else
          !     k2 = cell_index
          !     k1 = cell_index + 1
          !     coef = (hwm-(l(3)-zc(k1)))/dzc(k2)
          !   end if
          !   do j = 1,n(2)
          !     ! do i = 0,n(1)
          !     do i = 1,n(1)
          !       u1 = u(i,j,k1)
          !       u2 = u(i,j,k2)
          !       v1 = 0.25_rp*(v(i,j,k1) + v(i+1,j,k1) + v(i,j-1,k1) + v(i+1,j-1,k1))
          !       v2 = 0.25_rp*(v(i,j,k2) + v(i+1,j,k2) + v(i,j-1,k2) + v(i+1,j-1,k2))
          !       u_mag = bcu_mag%z(i,j,ibound)
          !       v_mag = 0.25_rp*(bcv_mag%z(i,j  ,ibound) + bcv_mag%z(i+1,j  ,ibound) + &
          !                        bcv_mag%z(i,j-1,ibound) + bcv_mag%z(i+1,j-1,ibound))
          !       uh = vel_relative(u1,u2,coef,u_mag)
          !       vh = vel_relative(v1,v2,coef,v_mag)
          !       inputs%vel1_p1%z(i,j,ibound) = uh
          !       inputs%vel2_p1%z(i,j,ibound) = vh
          !       inputs% vel_p1%z(i,j,ibound) = sqrt(uh**2 + vh**2)
          !       inputs% hwm_p1%z(i,j,ibound) = hwm
          !       inputs%visc_p1%z(i,j,ibound) = visc
          !       rewards%vel_p1%z(i,j,ibound) = -abs(sqrt(uh**2 + vh**2) - 0.853_rp)
          !     end do
          !   end do
          !   ! do j = 0,n(2)
          !   do j = 1,n(2)
          !     do i = 1,n(1)
          !       u1 = 0.25_rp*(u(i-1,j,k1) + u(i,j,k1) + u(i-1,j+1,k1) + u(i,j+1,k1))
          !       u2 = 0.25_rp*(u(i-1,j,k2) + u(i,j,k2) + u(i-1,j+1,k2) + u(i,j+1,k2))
          !       v1 = v(i,j,k1)
          !       v2 = v(i,j,k2)
          !       u_mag = 0.25_rp*(bcu_mag%z(i-1,j  ,ibound) + bcu_mag%z(i,j  ,ibound) + &
          !                        bcu_mag%z(i-1,j+1,ibound) + bcu_mag%z(i,j+1,ibound))
          !       v_mag = bcv_mag%z(i,j,ibound)
          !       uh = vel_relative(u1,u2,coef,u_mag)
          !       vh = vel_relative(v1,v2,coef,v_mag)
          !       inputs%vel1_p2%z(i,j,ibound) = uh
          !       inputs%vel2_p2%z(i,j,ibound) = vh
          !       inputs% vel_p2%z(i,j,ibound) = sqrt(uh**2 + vh**2)
          !       inputs% hwm_p2%z(i,j,ibound) = hwm
          !       inputs%visc_p2%z(i,j,ibound) = visc
          !       rewards%vel_p2%z(i,j,ibound) = -abs(sqrt(uh**2 + vh**2) - 0.853_rp)
          case(3)  ! z-direction walls: vel1 = u, vel2 = v

            do j = 1,n(2) ! cell centers, should consider interval here
              do i = 1,n(1)
                !
                ! cell_index = hwm_index(ibound,idir)
                cell_index = inputs%hwm_index%z(i,j,ibound)
                this_hwm = inputs%hwm%z(i,j,ibound)
                ! write(55,'(3i10,f10.3,i10)') i, j, ibound, this_hwm, cell_index
                if(ibound == 0) then
                  k2 = cell_index
                  k1 = cell_index - 1
                  coef = (this_hwm-zc(k1))/dzc(k1)
                  u_ref = u_ref_0
                else
                  k2 = cell_index
                  k1 = cell_index + 1
                  coef = (this_hwm-(l(3)-zc(k1)))/dzc(k2)
                  u_ref = u_ref_1
                end if
                !
                u1 = 0.5_rp*(u(i-1,j,k1) + u(i,j,k1))
                u2 = 0.5_rp*(u(i-1,j,k2) + u(i,j,k2))
                v1 = 0.5_rp*(v(i,j-1,k1) + v(i,j,k1))
                v2 = 0.5_rp*(v(i,j-1,k2) + v(i,j,k2))
                u_mag = 0.5_rp*(bcu_mag%z(i-1,j,ibound) + bcu_mag%z(i,j,ibound))
                v_mag = 0.5_rp*(bcv_mag%z(i,j-1,ibound) + bcv_mag%z(i,j,ibound))
                uh = vel_relative(u1,u2,coef,u_mag)
                vh = vel_relative(v1,v2,coef,v_mag)
                inputs%vel1%z(i,j,ibound) = uh
                inputs%vel2%z(i,j,ibound) = vh
                inputs% vel%z(i,j,ibound) = sqrt(uh**2 + vh**2)
                inputs%visc%z(i,j,ibound) = visc*20000._rp   ! 
                ! rewards%vel%z(i,j,ibound) = -abs(sqrt(uh**2 + vh**2) - 0.8045_rp)
                rewards%vel1%z(i,j,ibound) = ((num_samples-1)/float(num_samples))*rewards%vel1%z(i,j,ibound) + & 
                                             (1              /float(num_samples))*uh ! move this out?
                ! reward contains the time-averaged velocity
                ! now, the reference velocity is the WMLES velocity, so we do not have to exclude the wall-modeled layer
                u_profile = 0.5_rp*(u(i-1,j,1:n(3)) + u(i,j,1:n(3)))
                u_profile_ave(:,i,j,ibound) = ((num_samples-1)/float(num_samples))*u_profile_ave(:,i,j,ibound) + &
                                              (1              /float(num_samples))*u_profile
                rewards%vel1_profile_err%z(i,j,ibound) = sum(dzf(1:n(3))*(u_profile_ave(:,i,j,ibound) - u_ref)**2) ! weighted L2 norm
                ! time-averaged profile, always from wall to center
                if(ibound == 0) then
                  rewards%vel1_profile%z(:,i,j,ibound) = u_profile_ave(1:n(3)   ,i,j,ibound)
                else
                  rewards%vel1_profile%z(:,i,j,ibound) = u_profile_ave(n(3):1:-1,i,j,ibound)
                end if
              end do
            end do
            ! if(ibound==1) then
            !   call MPI_FINALIZE(ierr)
            !   stop
            ! end if
            ! do j = 0,n(2) ! edge midpoints
            !   do i = 0,n(1)
            !     u1 = 0.5_rp*(u(i,j,k1) + u(i,j+1,k1))
            !     u2 = 0.5_rp*(u(i,j,k2) + u(i,j+1,k2))
            !     v1 = 0.5_rp*(v(i,j,k1) + v(i+1,j,k1))
            !     v2 = 0.5_rp*(v(i,j,k2) + v(i+1,j,k2))
            !     u_mag = 0.5_rp*(bcu_mag%z(i,j,ibound) + bcu_mag%z(i,j+1,ibound))
            !     v_mag = 0.5_rp*(bcv_mag%z(i,j,ibound) + bcv_mag%z(i+1,j,ibound))
            !     uh = vel_relative(u1,u2,coef,u_mag)
            !     vh = vel_relative(v1,v2,coef,v_mag)
            !     inputs%vel1%z(i,j,ibound) = uh
            !     inputs%vel2%z(i,j,ibound) = vh
            !     inputs% vel%z(i,j,ibound) = sqrt(uh**2 + vh**2)
            !     inputs% hwm%z(i,j,ibound) = hwm
            !     inputs%visc%z(i,j,ibound) = visc
            !     ! rewards%vel%z(i,j,ibound) = -abs(sqrt(uh**2 + vh**2) - 0.8045_rp)
            !     rewards%vel%z(i,j,ibound) = ((num_samples-1)/float(num_samples))*rewards%vel%z(i,j,ibound) + & 
            !                                 (1              /float(num_samples))*sqrt(uh**2 + vh**2)
            !   end do
            ! end do
          end select
          ! if(mod(istep,action_interval) == action_interval-1) then
          !   do k = 1,n(3)
          !     write(55,'(*(e20.10))') zc(k), u_profile_ave(k,1,1,ibound), u_ref(k)
          !   end do
          !   call MPI_FINALIZE(ierr)
          !   stop
          ! end if
        end if
      end do
    end do
  end subroutine wallmodel_input
  
  subroutine assign_wallmodelbc(idir,ibound,visc,outputs,bcu,bcv,bcw)
    implicit none
    integer, intent(in) :: idir,ibound
    real(rp), intent(in) :: visc
    type(WallModelOutput), intent(in) :: outputs
    type(Bound), intent(inout) :: bcu,bcv,bcw
    real(rp) :: visci,sgn
    integer :: nx,ny
    
    visci = 1._rp/visc
    if(ibound == 0) then
      sgn =  1._rp
    else
      sgn = -1._rp
    end if

    select case(idir)
    ! case(1)
    !   bcv%x(:,:,ibound) = sgn*visci*outputs%tauw1_p1%x(:,:,ibound)
    !   bcw%x(:,:,ibound) = sgn*visci*outputs%tauw2_p2%x(:,:,ibound)
    ! case(2)
    !   bcu%y(:,:,ibound) = sgn*visci*outputs%tauw1_p1%y(:,:,ibound)
    !   bcw%y(:,:,ibound) = sgn*visci*outputs%tauw2_p2%y(:,:,ibound)
    case(3)
      ! bcu%z(:,:,ibound) = sgn*visci*outputs%tauw1_p1%z(:,:,ibound)
      ! bcv%z(:,:,ibound) = sgn*visci*outputs%tauw2_p2%z(:,:,ibound)
      nx = size(bcu%z,1) - 2
      ny = size(bcu%z,2) - 2
      bcu%z(0:nx,1:ny,ibound) = sgn*visci*0.5_rp*(outputs%tauw1%z(0:nx  ,1:ny  ,ibound) + &
                                                  outputs%tauw1%z(1:nx+1,1:ny  ,ibound))
      bcv%z(1:nx,0:ny,ibound) = sgn*visci*0.5_rp*(outputs%tauw2%z(1:nx  ,0:ny  ,ibound) + &
                                                  outputs%tauw2%z(1:nx  ,1:ny+1,ibound))
    end select
  end subroutine assign_wallmodelbc

  subroutine wallmodel_loglaw(visc,hwm,inputs_array,outputs_array,rewards_array)
    implicit none
    real(rp), intent(in) :: visc,hwm
    real(rp), intent(in), dimension(:,:) :: inputs_array
    real(rp), intent(out), dimension(:,:) :: outputs_array
    real(rp), intent(in), dimension(:,:), optional :: rewards_array
    real(rp) :: u1,u2,upar,utau,conv,utau_old,f,fp,tauw_tot,tauw1,tauw2,this_hwm
    integer :: n_points,i,i_stag

    n_points = size(inputs_array,1)

    do i = 1,n_points
      this_hwm = inputs_array(i,4)
      u1 = inputs_array(i,1)
      u2 = inputs_array(i,2)
      upar = sqrt(u1**2 + u2**2)
      utau = max(sqrt(upar/this_hwm*visc),visc/this_hwm*exp(-kap_log*b_log))
      conv = 1._rp
      do while(conv > 0.5e-4_rp)
        utau_old = utau
        f = upar/utau - 1._rp/kap_log*log(this_hwm*utau/visc) - b_log
        fp = -1._rp/utau*(upar/utau + 1._rp/kap_log)
        utau = abs(utau - f/fp)
        conv = abs(utau/utau_old - 1._rp)
      end do
      tauw_tot = utau**2
      tauw1 = tauw_tot*u1/(upar+eps)
      tauw2 = tauw_tot*u2/(upar+eps)
      outputs_array(i,1) = tauw1
      outputs_array(i,2) = tauw2
    end do
  end subroutine wallmodel_loglaw

  subroutine wallmodel_laminar(visc,hwm,inputs_array,outputs_array,rewards_array)
    implicit none
    real(rp), intent(in) :: visc,hwm
    real(rp), intent(in), dimension(:,:) :: inputs_array
    real(rp), intent(out), dimension(:,:) :: outputs_array
    real(rp), intent(in), dimension(:,:), optional :: rewards_array
    real(rp) :: u1,u2,upar,umax,del,tauw_tot,tauw1,tauw2,this_hwm
    integer :: n_points,i,i_stag

    n_points = size(inputs_array,1)
    del = 1._rp

    do i = 1,n_points
      this_hwm = inputs_array(i,4)
      u1 = inputs_array(i,1)
      u2 = inputs_array(i,2)
      upar = sqrt(u1**2 + u2**2)
      umax = upar/(this_hwm/del*(2._rp-this_hwm/del))
      tauw_tot = 2._rp/del*umax*visc
      tauw1 = tauw_tot*u1/(upar+eps)
      tauw2 = tauw_tot*u2/(upar+eps)
      outputs_array(i,1) = tauw1
      outputs_array(i,2) = tauw2
    end do
  end subroutine wallmodel_laminar

  subroutine wallmodel_DRL(visc,hwm,inputs_array,outputs_array,rewards_array)
    implicit none
    real(rp), intent(in) :: visc,hwm
    real(rp), intent(in), dimension(:,:) :: inputs_array
    real(rp), intent(out), dimension(:,:) :: outputs_array
    real(rp), intent(in), dimension(:,:), optional :: rewards_array
    real(rp) :: u1,u2,upar,tauw_tot,tauw1,tauw2
    logical, save :: is_first = .true.
    integer, save :: istep,n_points,n_vars_state,n_vars_action,n_vars_reward
    integer :: i,i_var

    real(rp), allocatable, dimension(:,:), save :: state
    real(rp), allocatable, dimension(:,:), save :: action
    real(rp), allocatable, dimension(:,:), save :: reward

    if(is_first) then
      is_first = .false.
      istep = 0 ! starting from 0 even when the simulation is restarted
      call init_smartredis(db_clustered)
      n_points      = size(inputs_array ,1)
      n_vars_state  = size(inputs_array ,2)
      n_vars_action = size(outputs_array,2)
      n_vars_reward = size(rewards_array,2)
      allocate(state (n_vars_state ,n_points))
      allocate(action(n_vars_action,n_points))
      allocate(reward(n_vars_reward,n_points))
      print *, "n_points = ", n_points
    else
      istep = istep + action_interval
    end if
    do i_var = 1,n_vars_state
      state (i_var,:) = inputs_array (:,i_var)
    end do
    do i_var = 1,n_vars_reward
      reward(i_var,:) = rewards_array(:,i_var)
    end do
    print *, "istep = ", istep, "total_time_steps = ", total_time_steps
    ! s_t, r_t
    call put_state( "ensemble_"//trim(adjustl(tag))//".state" ,shape(state (3:5            ,:)),state (3:5            ,:))
    if(istep /= 0) then
      call put_reward("ensemble_"//trim(adjustl(tag))//".reward",shape(reward(1:n_vars_reward,:)),reward(1:n_vars_reward,:))
    end if
    ! a_t, no action for the ending state
    call get_action("ensemble_"//trim(adjustl(tag))//".action",shape(action(3              ,:)),action(3              ,:))
    do i_var = 1,n_vars_action
      outputs_array(:,i_var) = action(i_var,:)
    end do
    !
    do i = 1,n_points
      u1 = inputs_array(i,1)
      u2 = inputs_array(i,2)
      upar = sqrt(u1**2 + u2**2)
      tauw_tot = outputs_array(i,3)
      tauw1 = tauw_tot*u1/(upar+eps)
      tauw2 = tauw_tot*u2/(upar+eps)
      outputs_array(i,1) = tauw1
      outputs_array(i,2) = tauw2
    end do

  end subroutine wallmodel_DRL

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