module mod_wallmodel
  use, intrinsic :: ieee_arithmetic, only: is_nan => ieee_is_nan, is_finite => ieee_is_finite
  use mpi
  use mod_precision
  use mod_typedef, only: bound
  use mod_param, only: kap_log, b_log, eps
  implicit none
  private
  public :: updt_wallmodelbc, init_wallmodels

  ! Define wall model types (extend as needed)
  integer, parameter :: WM_LAM = 1, WM_LOG = 2
  integer, parameter :: max_wall_models = 10

  ! Abstract interface for wall models
  abstract interface
    subroutine wallmodel_interface(inputs, outputs)
      import :: WallModelInput, WallModelOutput
      type(WallModelInput), intent(in) :: inputs
      type(WallModelOutput), intent(out) :: outputs
    end subroutine wallmodel_interface
  end interface

  ! Procedure pointer array for wall models
  procedure(wallmodel_interface), pointer :: wallmodel_ptr(max_wall_models)

  ! Derived type for wall model inputs
  type :: WallModelInput
    integer :: idir            ! Direction index (1=x, 2=y, 3=z)
    integer :: ibound          ! Boundary index (0 or 1)
    integer :: index           ! Wall index
    real(rp) :: h              ! Height from wall
    real(rp) :: visc           ! Viscosity
    real(rp), allocatable :: vel1(:,:)  ! Relevant velocity component 1
    real(rp), allocatable :: vel2(:,:)  ! Relevant velocity component 2
    ! Add additional input variables as needed for new wall models
  end type WallModelInput

  ! Derived type for wall model outputs
  type :: WallModelOutput
    real(rp), allocatable :: tauw1(:,:) ! Wall shear stress for vel1
    real(rp), allocatable :: tauw2(:,:) ! Wall shear stress for vel2
    ! Add additional output variables as needed for new wall models
  end type WallModelOutput

contains

  ! Initialization routine to assign wall model subroutines
  subroutine init_wallmodels()
    integer :: m
    ! Initialize all procedure pointers to null
    do m = 1, max_wall_models
      nullify(wallmodel_ptr(m))
    end do
    ! Assign specific wall models to their corresponding mtype
    wallmodel_ptr(WM_LAM) => wallmodel_lam
    wallmodel_ptr(WM_LOG) => wallmodel_log
    ! Assign additional wall models here, e.g.,
    ! wallmodel_ptr(WM_NEW) => wallmodel_new
  end subroutine init_wallmodels

  ! Dispatcher to call the appropriate wall model based on mtype
  subroutine wallmodel_dispatcher(mtype, inputs, outputs)
    integer, intent(in) :: mtype
    type(WallModelInput), intent(in) :: inputs
    type(WallModelOutput), intent(out) :: outputs

    if (mtype >= 1 .and. mtype <= max_wall_models) then
      if (associated(wallmodel_ptr(mtype))) then
        call wallmodel_ptr(mtype)(inputs, outputs)
      else
        print *, "Error: Wall model type ", mtype, " is not implemented."
        stop
      end if
    else
      print *, "Error: Invalid wall model type ", mtype
      stop
    end if
  end subroutine wallmodel_dispatcher

  ! Specific Wall Model Implementations

  ! Logarithmic Wall Model
  subroutine wallmodel_log(inputs, outputs)
    type(WallModelInput), intent(in) :: inputs
    type(WallModelOutput), intent(out) :: outputs
    integer :: i, j
    real(rp) :: upar, utau, conv, utau_old, f, fp, tauw_tot

    ! Allocate output arrays
    allocate(outputs%tauw1(size(inputs%vel1,1), size(inputs%vel1,2)))
    allocate(outputs%tauw2(size(inputs%vel2,1), size(inputs%vel2,2)))
    ! Initialize if necessary
    outputs%tauw1 = 0.0_rp
    outputs%tauw2 = 0.0_rp

    do j = 1, size(inputs%vel1,2)
      do i = 1, size(inputs%vel1,1)
        upar = sqrt(inputs%vel1(i,j)**2 + inputs%vel2(i,j)**2)
        utau = max(sqrt(upar/inputs%h*inputs%visc), inputs%visc/inputs%h*exp(-kap_log*b_log))
        conv = 1._rp
        do while(conv > 0.5e-4_rp)
          utau_old = utau
          f = upar/utau - 1._rp/kap_log*log(inputs%h*utau/inputs%visc) - b_log
          fp = -1._rp/utau*(upar/utau + 1._rp/kap_log)
          utau = abs(utau - f / fp)
          conv = abs(utau / utau_old - 1._rp)
        end do
        tauw_tot = utau**2
        outputs%tauw1(i,j) = tauw_tot * inputs%vel1(i,j) / (upar + eps)
        outputs%tauw2(i,j) = tauw_tot * inputs%vel2(i,j) / (upar + eps)
      end do
    end do
  end subroutine wallmodel_log

  ! Laminar Wall Model
  subroutine wallmodel_lam(inputs, outputs)
    type(WallModelInput), intent(in) :: inputs
    type(WallModelOutput), intent(out) :: outputs
    integer :: i, j
    real(rp) :: upar, umax, del, tauw_tot

    ! Allocate output arrays
    allocate(outputs%tauw1(size(inputs%vel1,1), size(inputs%vel1,2)))
    allocate(outputs%tauw2(size(inputs%vel2,1), size(inputs%vel2,2)))
    ! Initialize if necessary
    outputs%tauw1 = 0.0_rp
    outputs%tauw2 = 0.0_rp

    del = 1.0_rp  ! 1.0 is channel half height; adjust as needed
    do j = 1, size(inputs%vel1,2)
      do i = 1, size(inputs%vel1,1)
        upar = sqrt(inputs%vel1(i,j)**2 + inputs%vel2(i,j)**2)
        umax = upar / (inputs%h/del * (2._rp - inputs%h/del))
        tauw_tot = 2._rp / del * umax * inputs%visc
        outputs%tauw1(i,j) = tauw_tot * inputs%vel1(i,j) / (upar + eps)
        outputs%tauw2(i,j) = tauw_tot * inputs%vel2(i,j) / (upar + eps)
      end do
    end do
  end subroutine wallmodel_lam

  ! Add additional wall model subroutines here following the same pattern
  ! Example:
  ! subroutine wallmodel_new(inputs, outputs)
  !   type(WallModelInput), intent(in) :: inputs
  !   type(WallModelOutput), intent(out) :: outputs
  !   ! Implement new wall model logic
  ! end subroutine wallmodel_new

  ! Main Boundary Condition Update Subroutine
  subroutine updt_wallmodelbc(n, is_bound, lwm, l, dl, zc, zf, dzc, dzf, visc, h, index_wm, u, v, w, &
                              bcu, bcv, bcw, bcu_mag, bcv_mag, bcw_mag)
    implicit none
    integer, intent(in), dimension(3) :: n
    logical, intent(in), dimension(0:1,3) :: is_bound
    integer, intent(in), dimension(0:1,3) :: lwm, index_wm
    real(rp), intent(in), dimension(3) :: l, dl
    real(rp), intent(in), dimension(:) :: zc, zf, dzc, dzf
    real(rp), intent(in) :: visc, h
    real(rp), intent(in), dimension(:,:,:), intent(in) :: u, v, w
    type(bound), intent(inout) :: bcu, bcv, bcw
    type(bound), intent(in) :: bcu_mag, bcv_mag, bcw_mag

    type(WallModelInput) :: inputs
    type(WallModelOutput) :: outputs
    integer :: nh, mtype, idir, ibound, index

    nh = 1  ! Assuming nh=1 as per original code

    ! Loop over each possible wall boundary
    do ibound = 0, 1
      do idir = 1, 3
        if (is_bound(ibound, idir) .and. lwm(ibound, idir) /= 0) then
          mtype = lwm(ibound, idir)
          index = index_wm(ibound, idir)
          ! Compute and assign vel1 and vel2 based on direction and boundary
          call cmpt_wallmodelbc(n, ibound, idir, nh, mtype, l, dl, zc, zf, dzc, dzf, &
                                visc, h, index, u, v, w, bcu, bcv, bcw, bcu_mag, bcv_mag, bcw_mag)
        end if
      end do
    end do
  end subroutine updt_wallmodelbc

  ! Subroutine to compute wall model boundary conditions
  subroutine cmpt_wallmodelbc(n, ibound, idir, nh, mtype, l, dl, zc, zf, dzc, dzf, visc, h, index, &
                               u, v, w, bcvel1, bcvel2, bcvel3, bcu_mag, bcv_mag, bcw_mag)
    implicit none
    integer, intent(in) :: n(3), ibound, idir, nh, mtype
    real(rp), intent(in) :: l(3), dl(3)
    real(rp), intent(in) :: zc(:), zf(:), dzc(:), dzf(:)
    real(rp), intent(in) :: visc, h
    integer, intent(in) :: index
    real(rp), intent(in), dimension(:,:,:) :: u, v, w
    type(bound), intent(inout) :: bcvel1, bcvel2, bcvel3
    type(bound), intent(in) :: bcu_mag, bcv_mag, bcw_mag

    type(WallModelInput) :: inputs
    type(WallModelOutput) :: outputs

    ! Prepare WallModelInput
    call prepare_wallmodel_input(n, idir, ibound, index, l, dl, zc, zf, dzc, dzf, &
                                 visc, h, u, v, w, inputs)

    ! Call the dispatcher to compute wall shear stresses
    call wallmodel_dispatcher(mtype, inputs, outputs)

    ! Update boundary condition velocities based on outputs
    call update_boundary_conditions(idir, ibound, index, outputs, bcvel1, bcvel2, bcvel3)

    ! Deallocate dynamically allocated input arrays to prevent memory leaks
    if (allocated(inputs%vel1)) deallocate(inputs%vel1)
    if (allocated(inputs%vel2)) deallocate(inputs%vel2)
    ! Add deallocation for additional allocatable components if used
  end subroutine cmpt_wallmodelbc

  ! Subroutine to prepare WallModelInput based on wall direction and boundary
  subroutine prepare_wallmodel_input(n, idir, ibound, index, l, dl, zc, zf, dzc, dzf, &
                                     visc, h, u, v, w, inputs)
    implicit none
    integer, intent(in) :: n(3), idir, ibound, index
    real(rp), intent(in) :: l(3), dl(3)
    real(rp), intent(in) :: zc(:), zf(:), dzc(:), dzf(:)
    real(rp), intent(in) :: visc, h
    real(rp), intent(in), dimension(:,:,:) :: u, v, w
    type(WallModelInput), intent(out) :: inputs
    integer :: i, j, k

    ! Initialize inputs
    inputs%idir = idir
    inputs%ibound = ibound
    inputs%index = index
    inputs%h = h
    inputs%visc = visc
    inputs%vel1 = []
    inputs%vel2 = []

    ! Determine the velocity components based on direction
    select case(idir)
      case(1)  ! x-direction walls: vel1 = v, vel2 = w
        allocate(inputs%vel1(n(2), n(3)))
        allocate(inputs%vel2(n(2), n(3)))
        do j = 1, n(2)
          do k = 1, n(3)
            inputs%vel1(j,k) = v(index,j,k)
            inputs%vel2(j,k) = w(index,j,k)
          end do
        end do
      case(2)  ! y-direction walls: vel1 = u, vel2 = w
        allocate(inputs%vel1(n(1), n(3)))
        allocate(inputs%vel2(n(1), n(3)))
        do i = 1, n(1)
          do k = 1, n(3)
            inputs%vel1(i,k) = u(i,index,k)
            inputs%vel2(i,k) = w(i,index,k)
          end do
        end do
      case(3)  ! z-direction walls: vel1 = u, vel2 = v
        allocate(inputs%vel1(n(1), n(2)))
        allocate(inputs%vel2(n(1), n(2)))
        do i = 1, n(1)
          do j = 1, n(2)
            inputs%vel1(i,j) = u(i,j,index)
            inputs%vel2(i,j) = v(i,j,index)
          end do
        end do
      case default
        print *, "Error: Invalid wall direction ", idir
        stop
    end select

    ! Compute gradients if required by certain wall models
    ! Placeholder: Implement gradient computations as needed
    ! For example:
    ! if (some_condition_for_gradients) then
    !   allocate(inputs%grad_u(...))
    !   ... compute gradients ...
    ! end if
  end subroutine prepare_wallmodel_input

  ! Subroutine to update boundary condition velocities based on wall model outputs
  subroutine update_boundary_conditions(idir, ibound, index, outputs, bcu, bcv, bcw)
    implicit none
    integer, intent(in) :: idir, ibound, index
    type(WallModelOutput), intent(in) :: outputs
    type(bound), intent(inout) :: bcu, bcv, bcw

    select case(idir)
      case(1)  ! x-direction walls: update bcv and bcw
        if (ibound == 0) then
          bcv%x(index, :, :) = outputs%tauw1
          bcw%x(index, :, :) = outputs%tauw2
        else
          bcv%x(index, :, :) = outputs%tauw1
          bcw%x(index, :, :) = outputs%tauw2
        end if
      case(2)  ! y-direction walls: update bcu and bcw
        if (ibound == 0) then
          bcu%y(index, :, :) = outputs%tauw1
          bcw%y(index, :, :) = outputs%tauw2
        else
          bcu%y(index, :, :) = outputs%tauw1
          bcw%y(index, :, :) = outputs%tauw2
        end if
      case(3)  ! z-direction walls: update bcu and bcv
        if (ibound == 0) then
          bcu%z(index, :, :) = outputs%tauw1
          bcv%z(index, :, :) = outputs%tauw2
        else
          bcu%z(index, :, :) = outputs%tauw1
          bcv%z(index, :, :) = outputs%tauw2
        end if
      case default
        print *, "Error: Invalid wall direction ", idir
        stop
    end select
  end subroutine update_boundary_conditions

end module mod_wallmodel
