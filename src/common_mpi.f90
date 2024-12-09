! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_common_mpi
  implicit none
  public
  integer :: nprocs,myid,ierr
  integer :: halo(3)
  integer :: ipencil_axis
end module mod_common_mpi
