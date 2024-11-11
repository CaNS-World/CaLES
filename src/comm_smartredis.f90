!=======================================================================
! comm_smartredis.f90
! Implements core communication functions using SmartRedis
!=======================================================================
module Comm_SmartRedis
  use mpi
  use smartredis_mod  ! Hypothetical module for SmartRedis Fortran API
  use Comm_Interface  ! Ensures adherence to the core interfaces
  implicit none
  private
  public :: init_comm, finalize_comm
  public :: send_tensor, receive_tensor

contains

  subroutine init_comm(client, comm)
    !---------------------------------------------------------------------
    ! Initialize the SmartRedis client and MPI environment
    !---------------------------------------------------------------------
    type(SmartRedisClient), intent(out) :: client
    integer, intent(in) :: comm  ! MPI communicator
    integer :: ierr
    logical :: mpi_initialized

    ! Check if MPI is already initialized
    call MPI_Initialized(mpi_initialized, ierr)
    if (.not. mpi_initialized) then
      call MPI_Init(ierr)
    endif

    ! Initialize the SmartRedis client
    call SmartRedisClient_create(client)
  end subroutine init_comm

  subroutine finalize_comm(client)
    !---------------------------------------------------------------------
    ! Finalize the SmartRedis client and MPI environment
    !---------------------------------------------------------------------
    type(SmartRedisClient), intent(inout) :: client
    integer :: ierr
    logical :: mpi_finalized

    ! Finalize the SmartRedis client
    call SmartRedisClient_destroy(client)

    ! Finalize MPI if it was initialized here
    call MPI_Finalized(mpi_finalized, ierr)
    if (.not. mpi_finalized) then
      call MPI_Finalize(ierr)
    endif
  end subroutine finalize_comm

  subroutine send_tensor(client, local_tensor, global_name, comm)
    !---------------------------------------------------------------------
    ! Sends a tensor to the ML component via SmartRedis
    !---------------------------------------------------------------------
    type(SmartRedisClient), intent(inout) :: client
    real, intent(in) :: local_tensor(:)
    character(len=*), intent(in) :: global_name
    integer, intent(in) :: comm
    integer :: ierr, rank, size
    integer :: local_size, total_size
    integer, allocatable :: recv_counts(:), displs(:)
    real, allocatable :: global_tensor(:)

    ! Get MPI rank and size
    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, size, ierr)

    local_size = size(local_tensor)
    allocate(recv_counts(size), displs(size))

    ! Gather the sizes of local tensors from all ranks
    call MPI_Gather(local_size, 1, MPI_INTEGER, &
                    recv_counts, 1, MPI_INTEGER, 0, comm, ierr)

    ! Calculate displacements and total size on root
    if (rank == 0) then
      displs(1) = 0
      do i = 2, size
        displs(i) = displs(i-1) + recv_counts(i-1)
      end do
      total_size = displs(size) + recv_counts(size)
      allocate(global_tensor(total_size))
    end if

    ! Gather all local tensors to the root rank
    call MPI_Gatherv(local_tensor, local_size, MPI_REAL, &
                     global_tensor, recv_counts, displs, MPI_REAL, &
                     0, comm, ierr)

    ! Root rank sends the global tensor to the ML component
    if (rank == 0) then
      call client%put_tensor(global_name, global_tensor, (/total_size/))
      deallocate(global_tensor)
    end if

    deallocate(recv_counts, displs)
  end subroutine send_tensor

  subroutine receive_tensor(client, local_tensor, global_name, comm)
    !---------------------------------------------------------------------
    ! Receives a tensor from the ML component via SmartRedis
    !---------------------------------------------------------------------
    type(SmartRedisClient), intent(inout) :: client
    real, intent(out) :: local_tensor(:)
    character(len=*), intent(in) :: global_name
    integer, intent(in) :: comm
    integer :: ierr, rank, size
    integer :: local_size, total_size
    integer, allocatable :: send_counts(:), displs(:)
    real, allocatable :: global_tensor(:)

    ! Get MPI rank and size
    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, size, ierr)

    local_size = size(local_tensor)
    allocate(send_counts(size), displs(size))

    ! Gather the sizes of local tensors from all ranks
    call MPI_Gather(local_size, 1, MPI_INTEGER, &
                    send_counts, 1, MPI_INTEGER, 0, comm, ierr)

    ! Calculate displacements and total size on root
    if (rank == 0) then
      displs(1) = 0
      do i = 2, size
        displs(i) = displs(i-1) + send_counts(i-1)
      end do
      total_size = displs(size) + send_counts(size)
      allocate(global_tensor(total_size))

      ! Root rank retrieves the global tensor from the ML component
      call client%get_tensor(global_name, global_tensor, (/total_size/))
    end if

    ! Scatter the global tensor to all local tensors
    call MPI_Scatterv(global_tensor, send_counts, displs, MPI_REAL, &
                      local_tensor, local_size, MPI_REAL, &
                      0, comm, ierr)

    if (rank == 0) then
      deallocate(global_tensor)
    end if

    deallocate(send_counts, displs)
  end subroutine receive_tensor

end module Comm_SmartRedis