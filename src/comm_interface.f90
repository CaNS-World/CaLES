!=======================================================================
! comm_interface.f90
! Defines abstract interfaces for core communication functions
!=======================================================================
module Comm_Interface
  implicit none
  private
  public :: init_comm, finalize_comm
  public :: send_tensor, receive_tensor

  ! Abstract interfaces for core communication functions
  abstract interface

    subroutine init_comm(client, comm)
      type(SmartRedisClient), intent(out) :: client
      integer, intent(in) :: comm  ! MPI communicator
    end subroutine init_comm

    subroutine finalize_comm(client)
      type(SmartRedisClient), intent(inout) :: client
    end subroutine finalize_comm

    subroutine send_tensor(client, local_tensor, global_name, comm)
      type(SmartRedisClient), intent(inout) :: client
      real, intent(in) :: local_tensor(:)
      character(len=*), intent(in) :: global_name
      integer, intent(in) :: comm  ! MPI communicator
    end subroutine send_tensor

    subroutine receive_tensor(client, local_tensor, global_name, comm)
      type(SmartRedisClient), intent(inout) :: client
      real, intent(out) :: local_tensor(:)
      character(len=*), intent(in) :: global_name
      integer, intent(in) :: comm  ! MPI communicator
    end subroutine receive_tensor

  end interface

end module Comm_Interface