!=======================================================================
! comm_channel.f90
! Contains communication routines specific to the CHANNEL case
!=======================================================================
module Comm_Channel
  use mpi
  use Comm_Interface
  implicit none
  private
  public :: exchange_data_channel

contains

  subroutine exchange_data_wmodel(client, U, firstTimeStep, lastTimeStep, comm)
    !---------------------------------------------------------------------
    ! MODULES
    use mpi
    use Comm_Interface
    implicit none
    !---------------------------------------------------------------------
    ! INPUT/OUTPUT VARIABLES
    type(SmartRedisClient), intent(inout) :: client
    real, intent(inout) :: U(:,:,:,:,:)
    logical, intent(in) :: firstTimeStep, lastTimeStep
    integer, intent(in) :: comm
    !---------------------------------------------------------------------
    ! LOCAL VARIABLES
    integer :: ierr, rank, size, i
    character(len=255) :: key
    integer :: local_quantity_size, total_quantity_size
    integer, allocatable :: recv_counts(:), displs(:)
    real, allocatable :: local_quantity(:), global_quantity(:)
    real :: action_scalar
    !---------------------------------------------------------------------

    ! Get MPI rank and size
    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, size, ierr)

    ! Step 1: Compute quantities of interest (e.g., wall shear stress, velocity profiles)
    ! Placeholder computation: Replace with actual computations relevant to CHANNEL
    local_quantity_size = 1  ! Adjust based on actual quantities
    allocate(local_quantity(local_quantity_size))
    local_quantity(1) = compute_channel_quantity(U)

    ! Step 2: Gather quantities to root
    if (rank == 0) then
      allocate(recv_counts(size))
      allocate(displs(size))
    end if

    call MPI_Gather(local_quantity_size, 1, MPI_INTEGER, &
                    recv_counts, 1, MPI_INTEGER, 0, comm, ierr)

    ! Calculate displacements and total size on root
    if (rank == 0) then
      displs(1) = 0
      do i = 2, size
        displs(i) = displs(i-1) + recv_counts(i-1)
      end do
      total_quantity_size = displs(size) + recv_counts(size)
      allocate(global_quantity(total_quantity_size))
    end if

    ! Gather all local quantities to global_quantity on root
    call MPI_Gatherv(local_quantity, local_quantity_size, MPI_REAL, &
                     global_quantity, recv_counts, displs, MPI_REAL, &
                     0, comm, ierr)

    ! Step 3: On root, send the global quantities to the ML component
    if (rank == 0) then
      key = 'channel_state_tensor'
      call send_tensor(client, global_quantity, key, comm)
      deallocate(global_quantity)
    end if

    deallocate(local_quantity)
    if (rank == 0) then
      deallocate(recv_counts, displs)
    end if

    ! Step 4: On root, receive actions from the ML component
    if (rank == 0) then
      key = 'channel_action_tensor'
      call receive_tensor(client, action_scalar, key, comm)
    end if

    ! Step 5: Broadcast the action to all processes
    call MPI_Bcast(action_scalar, 1, MPI_REAL, 0, comm, ierr)

    ! Step 6: Apply the action to the simulation
    call apply_action_channel(U, action_scalar, firstTimeStep, lastTimeStep)

  end subroutine exchange_data_channel

  !-----------------------------------------------------------------------
  !> Computes a specific quantity relevant to the CHANNEL simulation
  real function compute_channel_quantity(U)
    real, intent(in) :: U(:,:,:,:,:)
    ! Implement the actual computation for CHANNEL here
    compute_channel_quantity = sum(U)  ! Replace with actual computation
  end function compute_channel_quantity

  !-----------------------------------------------------------------------
  !> Applies the received action to the simulation for the CHANNEL case
  subroutine apply_action_channel(U, action, firstTimeStep, lastTimeStep)
    real, intent(inout) :: U(:,:,:,:,:)
    real, intent(in) :: action
    logical, intent(in) :: firstTimeStep, lastTimeStep
    ! Implement CHANNEL-specific action application logic here
    ! Example: Modify U based on the received action
    U = U * action  ! Replace with actual application logic
  end subroutine apply_action_channel

end module Comm_Channel