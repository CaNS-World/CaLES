! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_comm_manager
  use mpi
  use mod_channel_comm
  implicit none! -
  !
  ! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
  ! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
  ! SPDX-License-Identifier: MIT
  !
  ! -
  private
  public :: define_parameters_comm_manager
  public :: init_comm_manager
  public :: finalize_comm_manager
  public :: exchange_data_comm_manager
  public :: analyze_comm_manager

  ! Parameters for communication types
  integer, parameter :: PRM_COMM_NONE     = 0
  integer, parameter :: PRM_COMM_CHANNEL  = 1
  ! Add more parameters as needed (e.g., PRM_COMM_HIT = 2)

  ! Variables
  integer :: Comm_Type
  logical :: doComm
  type(SmartRedisClient) :: Client
  integer :: comm
  integer :: ierr  ! For error handling

contains

  subroutine define_parameters()
    !---------------------------------------------------------------------
    ! Define communication parameters based on configuration or input
    !---------------------------------------------------------------------
    ! Example assignment; replace with actual configuration logic
    Comm_Type = PRM_COMM_CHANNEL  ! Set to PRM_COMM_HIT or other as needed
    doComm = (Comm_Type /= PRM_COMM_NONE)
  end subroutine define_parameters

  subroutine initialize()
    !---------------------------------------------------------------------
    ! Initialize communication client and set up communication
    !---------------------------------------------------------------------
    if (doComm) then
      call init_comm(Client, comm)
    endif
  end subroutine initialize

  subroutine finalize()
    !---------------------------------------------------------------------
    ! Finalize communication client
    !---------------------------------------------------------------------
    if (doComm) then
      call finalize_comm(Client)
    endif
  end subroutine finalize

  subroutine exchange_data(U, firstTimeStep, lastTimeStep)
    !---------------------------------------------------------------------
    ! Coordinates data exchange based on the communication type
    !---------------------------------------------------------------------
    real, intent(inout) :: U(:,:,:,:,:)
    logical, intent(in) :: firstTimeStep, lastTimeStep

    if (.not. doComm) return

    select case (Comm_Type)
    case (PRM_COMM_CHANNEL)
      call exchange_data_channel(Client, U, firstTimeStep, lastTimeStep, comm)
    case (PRM_COMM_HIT)
      call exchange_data_hit(Client, U, firstTimeStep, lastTimeStep, comm)
    case default
      print *, 'Unknown communication type: ', Comm_Type
      call MPI_ABORT(comm, 1, ierr)
    end select
  end subroutine exchange_data

  subroutine analyze()
    !---------------------------------------------------------------------
    ! Implement analysis routines related to communication if necessary
    !---------------------------------------------------------------------
    ! Placeholder for analysis logic
    ! ...
  end subroutine analyze

end module mod_comm_manager