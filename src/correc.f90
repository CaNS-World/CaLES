! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_correc
  use mod_precision, only: rp,sp,dp,i8,MPI_REAL_RP
  implicit none
  private
  public correc
  contains
  subroutine correc(n,dli,dzci,dt,p,u,v,w)
    !
    ! corrects the velocity so that it is divergence free
    ! when homogeneous Neumann boundary condition is imposed for pressure, the
    ! correction precedure does not change the wall-normal velocity at a wall,
    ! so the wall-normal velocity keeps its zero value from the prediction step,
    ! considering that impose_norm_bc=0 for the correction step.
    ! although the wall-parallel velocities at the ghost cells are also corrected,
    ! the correction is superfluous since the ghost cell values are replaced
    ! in the bounduvw step. Hence, the correction step does not use the ghost cells
    ! in essence, which implies that the ghost cell values of wall-paralel
    ! velocities are used only for computing the viscous terms and the velocity
    ! divergence at the right-hand side of the pressure Poisson equation.
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(3 ) :: dli
    real(rp), intent(in), dimension(0:) :: dzci
    real(rp), intent(in) :: dt
    real(rp), intent(in   ), dimension(0:,0:,0:) :: p
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    real(rp) :: factori,factorj
    !real(rp), dimension(0:n(3)+1) :: factork
    integer :: i,j,k
    !
    !factor = rkcoeffab(rkiter)*dt
    !
    factori = dt*dli(1)
    factorj = dt*dli(2)
    !factork = dt*dzci!dli(3)
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=0,n(3)+1
      do j=0,n(2)+1
        do i=0,n(1)
          u(i,j,k) = u(i,j,k) - factori*(   p(i+1,j,k)-p(i,j,k))
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=0,n(3)+1
      do j=0,n(2)
        do i=0,n(1)+1
          v(i,j,k) = v(i,j,k) - factorj*(   p(i,j+1,k)-p(i,j,k))
        end do
      end do
    end do
    !$acc parallel loop collapse(3) default(present) async(1)
    do k=0,n(3)
      do j=0,n(2)+1
        do i=0,n(1)+1
          w(i,j,k) = w(i,j,k) - dt*dzci(k)*(p(i,j,k+1)-p(i,j,k))
        end do
      end do
    end do
  end subroutine correc
end module mod_correc
