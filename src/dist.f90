! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_dist
  use mod_precision, only: rp,sp,dp,i8,MPI_REAL_RP
  use mod_params, only: big
  implicit none
  private
  public wall_dist
  contains
  subroutine wall_dist(cbc,n,is_bound,l,dl,zc,dzc,dw)
    !
    ! distance to the nearest wall. Identification of walls is based on
    ! non-penetration boundary conditions, which applies to both no-slip
    ! and wall model cases.
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in ), dimension(3)        :: n
    logical , intent(in ), dimension(0:1,3)    :: is_bound
    real(rp), intent(in ), dimension(3)        :: l,dl
    real(rp), intent(in ), dimension(0:)       :: zc,dzc
    real(rp), intent(out), dimension(0:,0:,0:) :: dw
    real(rp) :: this_dw
    integer  :: i,j,k
    !
    dw = big
    !
    if(is_bound(0,1).and.cbc(0,1,1)=='D') then
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            this_dw = dl(1)*(i-0.5)
            dw(i,j,k) = min(dw(i,j,k),this_dw)
          end do
        end do
      end do
    end if
    if(is_bound(1,1).and.cbc(1,1,1)=='D') then
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            this_dw = dl(1)*(n(1)-i+0.5)
            dw(i,j,k) = min(dw(i,j,k),this_dw)
          end do
        end do
      end do
    end if
    !
    if(is_bound(0,2).and.cbc(0,2,2)=='D') then
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            this_dw = dl(2)*(j-0.5)
            dw(i,j,k) = min(dw(i,j,k),this_dw)
          end do
        end do
      end do
    end if
    if(is_bound(1,2).and.cbc(1,2,2)=='D') then
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            this_dw = dl(2)*(n(2)-j+0.5)
            dw(i,j,k) = min(dw(i,j,k),this_dw)
          end do
        end do
      end do
    end if
    !
    if(is_bound(0,3).and.cbc(0,3,3)=='D') then
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            this_dw = zc(k)
            dw(i,j,k) = min(dw(i,j,k),this_dw)
          end do
        end do
      end do
    end if
    if(is_bound(1,3).and.cbc(1,3,3)=='D') then
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            this_dw = l(3)-zc(k)
            dw(i,j,k) = min(dw(i,j,k),this_dw)
          end do
        end do
      end do
    end if
  end subroutine wall_dist
end module mod_dist