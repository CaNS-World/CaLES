! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_output
  use mpi
  use decomp_2d_io
  use mod_common_mpi, only:ierr,myid
  use mod_precision, only: rp,sp,dp,i8,MPI_REAL_RP
  implicit none
  private
  public out0d,gen_alias,out1d,out1d_chan,out2d,out3d,write_log_output,write_visu_2d,write_visu_3d,out2d_duct
  public out1d_single_point_chan
  contains
  subroutine out0d(fname,n,var)
    !
    ! appends the first n entries of an array
    ! var to a file
    ! fname -> name of the file
    ! n     -> number of entries
    ! var   -> input array of real values
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in) :: n
    real(rp), intent(in), dimension(:) :: var
    integer :: iunit
    !
    if (myid == 0) then
      open(newunit=iunit,file=fname,position='append')
      write(iunit,'(*(E16.7e3))') var(1:n)
      close(iunit)
    end if
  end subroutine out0d
  !
  subroutine gen_alias(myid,datadir,fname,fname_alias)
    !
    ! this subroutine generates a symlink with name `fname_alias`, pointing to
    ! file `datadir//fname` using the `execute_command_line` intrinsic;
    ! it is called by task `myid`
    !
    integer, intent(in) :: myid
    character(len=*), intent(in) :: datadir,fname,fname_alias
    if(myid == 0) call execute_command_line('ln -sf '//trim(fname)//' '//trim(datadir)//fname_alias)
  end subroutine gen_alias
  !
  subroutine out1d(fname,ng,lo,hi,idir,l,dl,z_g,dz,p)
    !
    ! writes the profile of a variable averaged
    ! over two domain directions
    !
    ! fname -> name of the file
    ! ng    -> global domain sizes
    ! lo,hi -> upper and lower extents of the input array
    ! idir  -> direction of the profile
    ! dl,l  -> uniform grid spacing and length arrays
    ! z_g   -> global z coordinate array (grid is non-uniform in z)
    ! dz    -> local z grid spacing array (should work also with the global one)
    ! p     -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: z_g
    real(rp), intent(in), dimension(lo(3)-1:) :: dz
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(dp), allocatable, dimension(:) :: p1d
    integer :: i,j,k
    integer :: iunit
    real(dp) :: grid_area_ratio,p1d_s
    !
    allocate(p1d(ng(idir)))
    select case(idir)
    case(3)
      grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
      !$acc data copyout(p1d) async(1)
      !$acc kernels default(present) async(1)
      p1d(:) = 0._rp
      !$acc end kernels
      !$acc parallel loop gang default(present) private(p1d_s) async(1)
      do k=lo(3),hi(3)
        p1d_s = 0._rp
        !$acc loop vector collapse(2) reduction(+:p1d_s)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            p1d_s = p1d_s + p(i,j,k)
          end do
        end do
        p1d(k) = p1d_s*grid_area_ratio
      end do
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          write(iunit,'(2E16.7e3)') z_g(k),p1d(k)
        end do
        close(iunit)
      end if
    case(2)
      grid_area_ratio = dl(1)/(l(1)*l(3))
      !$acc data copyout(p1d) async(1)
      !$acc kernels default(present) async(1)
      p1d(:) = 0._rp
      !$acc end kernels
      !$acc parallel loop gang default(present) private(p1d_s) async(1)
      do j=lo(2),hi(2)
        p1d_s = 0._rp
        !$acc loop vector collapse(2) reduction(+:p1d_s)
        do k=lo(3),hi(3)
          do i=lo(1),hi(1)
            p1d_s = p1d_s + p(i,j,k)*dz(k)
          end do
        end do
        p1d(j) = p1d_s*grid_area_ratio
      end do
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(2),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do j=1,ng(2)
          write(iunit,'(2E16.7e3)') (j-.5)*dl(2),p1d(j)
        end do
        close(iunit)
      end if
    case(1)
      grid_area_ratio = dl(2)/(l(2)*l(3))
      !$acc data copyout(p1d) async(1)
      !$acc kernels default(present) async(1)
      p1d(:) = 0._rp
      !$acc end kernels
      !$acc parallel loop gang default(present) private(p1d_s) async(1)
      do i=lo(1),hi(1)
        p1d_s = 0._rp
        !$acc loop vector collapse(2) reduction(+:p1d_s)
        do k=lo(3),hi(3)
          do j=lo(2),hi(2)
            p1d_s = p1d_s + p(i,j,k)*dz(k)
          end do
        end do
        p1d(i) = p1d_s*grid_area_ratio
      end do
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1),ng(1),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do i=1,ng(1)
          write(iunit,'(2E16.7e3)') (i-.5)*dl(1),p1d(i)
        end do
        close(iunit)
      end if
    end select
  end subroutine out1d
  !
  subroutine out2d(fname,inorm,islice,p)
    use mod_common_mpi, only: ipencil => ipencil_axis
    !
    ! saves a planar slice of a scalar field into a binary file
    !
    ! fname  -> name of the output file
    ! inorm  -> plane is perpendicular to direction
    !           inorm (1,2,3)
    ! islice -> plane is of constant index islice
    !           in direction inorm
    ! p      -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in) :: inorm,islice
    real(rp), intent(in), dimension(:,:,:) :: p
    !
    select case(inorm)
    case(1) !normal to x --> yz plane
       call decomp_2d_write_plane(ipencil,p,inorm,islice,'.',fname,'dummy')
    case(2) !normal to y --> zx plane
       call decomp_2d_write_plane(ipencil,p,inorm,islice,'.',fname,'dummy')
    case(3) !normal to z --> xy plane
       call decomp_2d_write_plane(ipencil,p,inorm,islice,'.',fname,'dummy')
    end select
  end subroutine out2d
  !
  subroutine out3d(fname,nskip,p)
    use mod_common_mpi, only: ipencil => ipencil_axis
    !
    ! saves a 3D scalar field into a binary file
    !
    ! fname  -> name of the output file
    ! nskip  -> array with the step size for which the
    !           field is written; i.e.: [1,1,1]
    !           writes the full field
    ! p      -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: nskip
    real(rp), intent(in), dimension(:,:,:) :: p
    integer :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    !
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fname, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierr)
    disp = 0_MPI_OFFSET_KIND
#if !defined(_OPENACC)
    call decomp_2d_write_every(ipencil,p,nskip(1),nskip(2),nskip(3),fname,.true.)
#else
    !
    ! alternative way of writing a full 3D field that was needed in the past
    !
    block
      use decomp_2d
      use mod_load, only: io_field
      integer, dimension(3) :: ng,lo,hi
      ng(:) = [nx_global,ny_global,nz_global]
      select case(ipencil)
      case(1)
        lo(:) = xstart(:)
        hi(:) = xend(:)
      case(2)
        lo(:) = ystart(:)
        hi(:) = yend(:)
      case(3)
        lo(:) = zstart(:)
        hi(:) = zend(:)
      end select
      if(any(nskip /= 1) .and. myid == 0) &
        print*, 'Warning: `nskip` should be `[1,1,1]` if `io_field()` is used to output 3D field data'
      call io_field('w',fh,ng,[0,0,0],lo,hi,disp,p)
    end block
#endif
    call MPI_FILE_CLOSE(fh,ierr)
  end subroutine out3d
  !
  subroutine write_log_output(fname,fname_fld,varname,nmin,nmax,nskip,time,istep)
    !
    ! appends information about a saved binary file to a log file
    ! this file is used to generate a xdmf file for visualization of field data
    !
    ! fname     -> name of the output log file
    ! fname_fld -> name of the saved binary file (excluding the directory)
    ! varname   -> name of the variable that is saved
    ! nmin      -> first element of the field that is saved in each direction, e.g. [1,1,1]
    ! nmax      -> last  element of the field that is saved in each direction, e.g. [ng(1),ng(2),ng(3)]
    ! nskip     -> step size between nmin and nmax, e.g. [1,1,1] if the whole array is saved
    ! time      -> physical time
    ! istep     -> time step number
    !
    implicit none
    character(len=*), intent(in) :: fname,fname_fld,varname
    integer , intent(in), dimension(3) :: nmin,nmax,nskip
    real(rp), intent(in)               :: time
    integer , intent(in)               :: istep
    character(len=100) :: cfmt
    integer :: iunit
    !
    write(cfmt, '(A)') '(A30,A15,9I5,E16.7E3,I7)'
    if (myid  ==  0) then
      open(newunit=iunit,file=fname,position='append')
      write(iunit,trim(cfmt)) trim(fname_fld),trim(varname),nmin,nmax,nskip,time,istep
      close(iunit)
    end if
  end subroutine write_log_output
  !
  subroutine write_visu_3d(datadir,fname_bin,fname_log,varname,nmin,nmax,nskip,time,istep,p)
    !
    ! wraps the calls of out3d and write_log_output into the same subroutine
    !
    implicit none
    character(len=*), intent(in)          :: datadir,fname_bin,fname_log,varname
    integer , intent(in), dimension(3)    :: nmin,nmax,nskip
    real(rp), intent(in)                  :: time
    integer , intent(in)                  :: istep
    real(rp), intent(in), dimension(:,:,:) :: p
    !
    call out3d(trim(datadir)//trim(fname_bin),nskip,p)
    call write_log_output(trim(datadir)//trim(fname_log),trim(fname_bin),trim(varname),nmin,nmax,nskip,time,istep)
  end subroutine write_visu_3d
  !
  subroutine write_visu_2d(datadir,fname_bin,fname_log,varname,inorm,nslice,ng,time,istep,p)
    !
    ! wraps the calls of out2d and write-log_output into the same subroutine
    !
    implicit none
    character(len=*), intent(in)          :: datadir,fname_bin,fname_log,varname
    integer , intent(in)                  :: inorm,nslice
    integer , intent(in), dimension(3)    :: ng
    real(rp), intent(in)                  :: time
    integer , intent(in)                  :: istep
    real(rp), intent(in), dimension(:,:,:) :: p
    integer , dimension(3) :: nmin_2d,nmax_2d
    !
    call out2d(trim(datadir)//trim(fname_bin),inorm,nslice,p)
    select case(inorm)
    case(1)
      nmin_2d(:) = [nslice,1    ,1    ]
      nmax_2d(:) = [nslice,ng(2),ng(3)]
    case(2)
      nmin_2d(:) = [1    ,nslice,1    ]
      nmax_2d(:) = [ng(1),nslice,ng(3)]
    case(3)
      nmin_2d(:) = [1    ,1    ,nslice]
      nmax_2d(:) = [ng(1),ng(2),nslice]
    end select
    call write_log_output(trim(datadir)//trim(fname_log),trim(fname_bin),trim(varname),nmin_2d,nmax_2d,[1,1,1],time,istep)
  end subroutine write_visu_2d
  !
  subroutine out1d_chan(fname,ng,lo,hi,idir,l,dl,z_g,u,v,w)
    !
    ! for a channel with streamwise dir in x
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: z_g
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    real(dp), allocatable, dimension(:) :: um,vm,wm,u2,v2,w2,uw
    real(dp) :: buf01,buf02,buf03,buf04,buf05,buf06,buf07
    integer :: i,j,k
    integer :: iunit
    integer :: q
    real(dp) :: grid_area_ratio
    !
    q = ng(idir)
    select case(idir)
    case(3)
      grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
      allocate(um(q),vm(q),wm(q),u2(q),v2(q),w2(q),uw(q))
      !$acc data copyout(um,vm,wm,u2,v2,w2,uw) async(1)
      !$acc kernels default(present) async(1)
      um(:) = 0._rp
      vm(:) = 0._rp
      wm(:) = 0._rp
      u2(:) = 0._rp
      v2(:) = 0._rp
      w2(:) = 0._rp
      uw(:) = 0._rp
      !$acc end kernels
      !$acc parallel loop gang default(present) async(1) &
      !$acc private(buf01,buf02,buf03,buf04,buf05,buf06,buf07)
      do k=lo(3),hi(3)
        buf01 = 0._rp
        buf02 = 0._rp
        buf03 = 0._rp
        buf04 = 0._rp
        buf05 = 0._rp
        buf06 = 0._rp
        buf07 = 0._rp
        !$acc loop vector collapse(2) &
        !$acc reduction(+:buf01,buf02,buf03,buf04,buf05,buf06,buf07)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            ! located at zc
            buf01 = buf01 + u(i,j,k)
            buf02 = buf02 + v(i,j,k)
            buf03 = buf03 + 0.50*(w(i,j,k-1) + w(i,j,k))
            buf04 = buf04 + u(i,j,k)**2
            buf05 = buf05 + v(i,j,k)**2
            buf06 = buf06 + 0.50*(w(i,j,k)**2+w(i,j,k-1)**2)
            buf07 = buf07 + 0.25*(u(i-1,j,k) + u(i,j,k))* &
                                 (w(i,j,k-1) + w(i,j,k))
          end do
        end do
        um(k) = buf01*grid_area_ratio
        vm(k) = buf02*grid_area_ratio
        wm(k) = buf03*grid_area_ratio
        u2(k) = buf04*grid_area_ratio
        v2(k) = buf05*grid_area_ratio
        w2(k) = buf06*grid_area_ratio
        uw(k) = buf07*grid_area_ratio
      end do
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,um(1),ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vm(1),ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wm(1),ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,u2(1),ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,v2(1),ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,w2(1),ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uw(1),ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          write(iunit,'(8E16.7e3)') z_g(k),um(k),vm(k),wm(k), &
                                           u2(k),v2(k),w2(k), &
                                           uw(k)
        end do
        close(iunit)
      end if
    case(2)
    case(1)
    end select
  end subroutine out1d_chan
  !
  subroutine out2d_duct(fname,ng,lo,hi,idir,l,dl,z_g,u,v,w)
    !
    ! for a duct, cell-centered values
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: z_g
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    real(dp), allocatable, dimension(:,:) :: um,vm,wm,u2,v2,w2,uv,uw,vw
    real(dp) :: buf01,buf02,buf03,buf04,buf05,buf06,buf07,buf08,buf09
    integer :: i,j,k
    integer :: iunit
    integer :: p,q
    real(dp) :: x_g,y_g,grid_area_ratio
    !
    select case(idir) ! streamwise direction
    case(3)
    case(2)
    case(1)
      grid_area_ratio = dl(1)/l(1)
      p = ng(2)
      q = ng(3)
      allocate(um(p,q),vm(p,q),wm(p,q),u2(p,q),v2(p,q),w2(p,q),uv(p,q),uw(p,q),vw(p,q))
      !$acc data copyout(um,vm,wm,u2,v2,w2,uv,uw,vw) async(1)
      !$acc kernels default(present) async(1)
      um(:,:) = 0._rp
      vm(:,:) = 0._rp
      wm(:,:) = 0._rp
      u2(:,:) = 0._rp
      v2(:,:) = 0._rp
      w2(:,:) = 0._rp
      uv(:,:) = 0._rp
      uw(:,:) = 0._rp
      vw(:,:) = 0._rp
      !$acc end kernels
      !$acc parallel loop gang collapse(2) default(present) async(1) &
      !$acc private(buf01,buf02,buf03,buf04,buf05,buf06,buf07,buf08,buf09)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          buf01 = 0._rp
          buf02 = 0._rp
          buf03 = 0._rp
          buf04 = 0._rp
          buf05 = 0._rp
          buf06 = 0._rp
          buf07 = 0._rp
          buf08 = 0._rp
          buf09 = 0._rp
          !$acc loop vector reduction(+:buf01,buf02,buf03,buf04,buf05,buf06,buf07,buf08,buf09)
          do i=lo(1),hi(1)
            buf01 = buf01 + u(i,j,k)
            buf02 = buf02 + 0.5*(v(i,j-1,k)+v(i,j,k))
            buf03 = buf03 + 0.5*(w(i,j,k-1)+w(i,j,k))
            buf04 = buf04 + u(i,j,k)**2
            buf05 = buf05 + 0.5*(v(i,j-1,k)**2+v(i,j,k)**2)
            buf06 = buf06 + 0.5*(w(i,j,k-1)**2+w(i,j,k)**2)
            buf07 = buf07 + 0.25*(u(i-1,j,k) + u(i,j,k))* &
                                 (v(i,j-1,k) + v(i,j,k))
            buf08 = buf08 + 0.25*(u(i-1,j,k) + u(i,j,k))* &
                                 (w(i,j,k-1) + w(i,j,k))
            buf09 = buf09 + 0.25*(v(i,j-1,k) + v(i,j,k))* &
                                 (w(i,j,k-1) + w(i,j,k))
          end do
          um(j,k) = buf01*grid_area_ratio
          vm(j,k) = buf02*grid_area_ratio
          wm(j,k) = buf03*grid_area_ratio
          u2(j,k) = buf04*grid_area_ratio
          v2(j,k) = buf05*grid_area_ratio
          w2(j,k) = buf06*grid_area_ratio
          uv(j,k) = buf07*grid_area_ratio
          uw(j,k) = buf08*grid_area_ratio
          vw(j,k) = buf09*grid_area_ratio
        end do
      end do
      !$acc end data
      !$acc wait(1)
      call MPI_ALLREDUCE(MPI_IN_PLACE,um(1,1),ng(2)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vm(1,1),ng(2)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,wm(1,1),ng(2)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,u2(1,1),ng(2)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,v2(1,1),ng(2)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,w2(1,1),ng(2)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uv(1,1),ng(2)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,uw(1,1),ng(2)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,vw(1,1),ng(2)*ng(3),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        open(newunit=iunit,file=fname)
        do k=1,ng(3)
          do j=1,ng(2)
            y_g = (j-0.5)*dl(2)
            write(iunit,'(11E16.7e3)') y_g,z_g(k),um(j,k),vm(j,k),wm(j,k), &
                                                  u2(j,k),v2(j,k),w2(j,k), &
                                                  uv(j,k),uw(j,k),vw(j,k)
          end do
        end do
        close(iunit)
      end if
    end select
  end subroutine out2d_duct
  !
  subroutine out1d_single_point_chan(fname,ng,lo,hi,idir,l,dl,dzc_g,dzf_g,zc_g,zf_g,u,v,w,p,visct)
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: dzc_g,dzf_g,zc_g,zf_g
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w,p,visct
    real(dp), allocatable, dimension(:,:) :: buf
    real(dp) :: tmp_x,tmp_y,tmp_z
    integer :: i,j,k,q
    integer :: iunit
    integer :: nn,nvars
    character(len=30) :: cfmt
    real(dp) :: grid_area_ratio
    real(dp) :: buf01,buf02,buf03,buf04,buf05,buf06,buf07,buf08,buf09,buf10, &
                buf11,buf12,buf13,buf14,buf15,buf16,buf17,buf18,buf19,buf20, &
                buf21,buf22,buf23,buf24,buf25,buf26,buf27,buf28,buf29,buf30, &
                buf31,buf32,buf33,buf34,buf35,buf36,buf37,buf38
    real(dp) :: dudx_ip,dudx_im,dvdy_jp,dvdy_jm,dwdz_kp,dwdz_km,dudz,dwdx
    real(dp) :: s_ccc,s_pcc,s_cpc,s_ccp,s_pcp
    real(dp) :: div
    !
    nn = ng(idir)
    nvars = 38
    select case(idir)
    case(3)
      grid_area_ratio = dl(1)*dl(2)/(l(1)*l(2))
      allocate(buf(nvars,nn))
      !$acc enter data create(buf) async(1)
      !$acc kernels default(present) async(1)
      buf(:,:) = 0._rp
      !$acc end kernels
      !$acc parallel loop gang default(present) async(1) &
      !$acc private(buf01,buf02,buf03,buf04,buf05,buf06,buf07,buf08,buf09,buf10) &
      !$acc private(buf11,buf12,buf13,buf14,buf15,buf16,buf17,buf18,buf19,buf20) &
      !$acc private(buf21,buf22,buf23,buf24,buf25,buf26,buf27,buf28,buf29,buf30) &
      !$acc private(buf31,buf32,buf33,buf34,buf35,buf36,buf37,buf38) &
      !$acc private(tmp_x,tmp_y,tmp_z) &
      !$acc private(s_ccc,s_pcc,s_cpc,s_ccp,s_pcp) &
      !$acc private(dudx_ip,dudx_im,dvdy_jp,dvdy_jm,dwdz_kp,dwdz_km,dudz,dwdx)
      do k=lo(3),hi(3)
        buf01 = 0._rp
        buf02 = 0._rp
        buf03 = 0._rp
        buf04 = 0._rp
        buf05 = 0._rp
        buf06 = 0._rp
        buf07 = 0._rp
        buf08 = 0._rp
        buf09 = 0._rp
        buf10 = 0._rp
        buf11 = 0._rp
        buf12 = 0._rp
        buf13 = 0._rp
        buf14 = 0._rp
        buf15 = 0._rp
        buf16 = 0._rp
        buf17 = 0._rp
        buf18 = 0._rp
        buf19 = 0._rp
        buf20 = 0._rp
        buf21 = 0._rp
        buf22 = 0._rp
        buf23 = 0._rp
        buf24 = 0._rp
        buf25 = 0._rp
        buf26 = 0._rp
        buf27 = 0._rp
        !$acc loop vector collapse(2) &
        !$acc reduction(+:buf01,buf02,buf03,buf04,buf05,buf06,buf07,buf08,buf09,buf10) &
        !$acc reduction(+:buf11,buf12,buf13,buf14,buf15,buf16,buf17,buf18,buf19,buf20) &
        !$acc reduction(+:buf21,buf22,buf23,buf24,buf25,buf26,buf27)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            !
            ! velocity
            !
            buf01 = buf01  + u(i,j,k)
            buf02 = buf02  + v(i,j,k)
            buf03 = buf03  + w(i,j,k)
            buf04 = buf04  + u(i,j,k)**2
            buf05 = buf05  + v(i,j,k)**2
            buf06 = buf06  + w(i,j,k)**2
            buf07 = buf07  + 0.25_rp*(u(i,j,k+1) + u(i  ,j,k))* &
                                     (w(i,j,k  ) + w(i+1,j,k)) ! cell edge
            buf08 = buf08 + u(i,j,k)**3
            buf09 = buf09 + v(i,j,k)**3
            buf10 = buf10 + w(i,j,k)**3
            buf11 = buf11 + u(i,j,k)**4
            buf12 = buf12 + v(i,j,k)**4
            buf13 = buf13 + w(i,j,k)**4
            !
            ! pressure
            !
            buf14 = buf14 + p(i,j,k)
            buf15 = buf15 + p(i,j,k)**2
            !
            ! vorticity
            !
            ! x component
            !
            tmp_x = (w(i,j+1,k)-w(i,j,k))/dl(2) - (v(i,j,k+1)-v(i,j,k))/dzc_g(k)
            !
            ! y component
            !
            tmp_y = (u(i,j,k+1)-u(i,j,k))/dzc_g(k) - (w(i+1,j,k)-w(i,j,k))/dl(1)
            !
            ! z component
            !
            tmp_z = (v(i+1,j,k)-v(i,j,k))/dl(1) - (u(i,j+1,k)-u(i,j,k))/dl(2)
            !
            buf16 = buf16 + tmp_x
            buf17 = buf17 + tmp_y
            buf18 = buf18 + tmp_z
            buf19 = buf19 + tmp_x**2
            buf20 = buf20 + tmp_y**2
            buf21 = buf21 + tmp_z**2
            !
            ! modelled (subgrid) stress
            ! stored at the same locations as resolved <uu>, <vv>, <ww>, <uw>
            !
            s_ccc = visct(i  ,j  ,k  )
            s_pcc = visct(i+1,j  ,k  )
            s_cpc = visct(i  ,j+1,k  )
            s_ccp = visct(i  ,j  ,k+1)
            s_pcp = visct(i+1,j  ,k+1)
            !
            dudx_ip = (u(i+1,j  ,k  )-u(i  ,j  ,k  ))/dl(1)
            dudx_im = (u(i  ,j  ,k  )-u(i-1,j  ,k  ))/dl(1)
            dvdy_jp = (v(i  ,j+1,k  )-v(i  ,j  ,k  ))/dl(2)
            dvdy_jm = (v(i  ,j  ,k  )-v(i  ,j-1,k  ))/dl(2)
            dwdz_kp = (w(i  ,j  ,k+1)-w(i  ,j  ,k  ))/dzf_g(k+1)
            dwdz_km = (w(i  ,j  ,k  )-w(i  ,j  ,k-1))/dzf_g(k  )
            dudz    = (u(i  ,j  ,k+1)-u(i  ,j  ,k  ))/dzc_g(k  )
            dwdx    = (w(i+1,j  ,k  )-w(i  ,j  ,k  ))/dl(1)
            !
            buf22 = buf22 - 0.5_rp*(s_pcc*(dudx_ip+dudx_ip) + s_ccc*(dudx_im+dudx_im))
            buf23 = buf23 - 0.5_rp*(s_cpc*(dvdy_jp+dvdy_jp) + s_ccc*(dvdy_jm+dvdy_jm))
            buf24 = buf24 - 0.5_rp*(s_ccp*(dwdz_kp+dwdz_kp) + s_ccc*(dwdz_km+dwdz_km))
            buf25 = buf25 - 0.25_rp*(s_ccc+s_pcc+s_ccp+s_pcp)*(dudz+dwdx) ! cell edge
            buf26 = buf26 + visct(i,j,k)
            !
            ! viscous shear stress, uw
            !
            buf27 = buf27 + dudz ! cell edge
          end do
        end do
        buf( 1,k) = buf01*grid_area_ratio
        buf( 2,k) = buf02*grid_area_ratio
        buf( 3,k) = buf03*grid_area_ratio
        buf( 4,k) = buf04*grid_area_ratio
        buf( 5,k) = buf05*grid_area_ratio
        buf( 6,k) = buf06*grid_area_ratio
        buf( 7,k) = buf07*grid_area_ratio
        buf( 8,k) = buf08*grid_area_ratio
        buf( 9,k) = buf09*grid_area_ratio
        buf(10,k) = buf10*grid_area_ratio
        buf(11,k) = buf11*grid_area_ratio
        buf(12,k) = buf12*grid_area_ratio
        buf(13,k) = buf13*grid_area_ratio
        buf(14,k) = buf14*grid_area_ratio
        buf(15,k) = buf15*grid_area_ratio
        buf(16,k) = buf16*grid_area_ratio
        buf(17,k) = buf17*grid_area_ratio
        buf(18,k) = buf18*grid_area_ratio
        buf(19,k) = buf19*grid_area_ratio
        buf(20,k) = buf20*grid_area_ratio
        buf(21,k) = buf21*grid_area_ratio
        buf(22,k) = buf22*grid_area_ratio
        buf(23,k) = buf23*grid_area_ratio
        buf(24,k) = buf24*grid_area_ratio
        buf(25,k) = buf25*grid_area_ratio
        buf(26,k) = buf26*grid_area_ratio
        buf(27,k) = buf27*grid_area_ratio
      end do
      !$acc update self(buf) async(1)
      !$acc wait(1)
      call mpi_allreduce(MPI_IN_PLACE,buf(1,1),size(buf),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        nvars = 27
        write(cfmt,'(A)') '(*(es24.16e3,1x))'
        open(newunit=iunit,file=fname//'.out')
        do k=1,nn
          write(iunit,trim(cfmt)) zc_g(k),zf_g(k),(buf(i,k),i=1,nvars),dzc_g(k),dzf_g(k)
        end do
        close(iunit)
        open(10,file=trim(fname)//'.bin',access='stream')
        write(10) buf(1:nvars,1:nn)
        close(10)
      end if
      !
      ! MKE and Reynolds shear stresses budgets
      !
      !$acc parallel loop gang default(present) async(1) &
      !$acc private(buf01,buf02,buf03,buf04,buf05,buf06,buf07,buf08,buf09,buf10) &
      !$acc private(buf11,buf12,buf13,buf14,buf15,buf16,buf17,buf18,buf19,buf20) &
      !$acc private(buf21,buf22,buf23,buf24,buf25,buf26,buf27,buf28,buf29,buf30) &
      !$acc private(buf31,buf32,buf33,buf34,buf35,buf36,buf37,buf38)
      do k=lo(3),hi(3)
        buf01 = 0._rp
        buf02 = 0._rp
        buf03 = 0._rp
        buf04 = 0._rp
        buf05 = 0._rp
        buf06 = 0._rp
        buf07 = 0._rp
        buf08 = 0._rp
        buf09 = 0._rp
        buf10 = 0._rp
        buf11 = 0._rp
        buf12 = 0._rp
        buf13 = 0._rp
        buf14 = 0._rp
        buf15 = 0._rp
        buf16 = 0._rp
        buf17 = 0._rp
        buf18 = 0._rp
        buf19 = 0._rp
        buf20 = 0._rp
        buf21 = 0._rp
        buf22 = 0._rp
        buf23 = 0._rp
        buf24 = 0._rp
        buf25 = 0._rp
        buf26 = 0._rp
        buf27 = 0._rp
        buf28 = 0._rp
        buf29 = 0._rp
        buf30 = 0._rp
        buf31 = 0._rp
        buf32 = 0._rp
        buf33 = 0._rp
        buf34 = 0._rp
        buf35 = 0._rp
        buf36 = 0._rp
        buf37 = 0._rp
        buf38 = 0._rp
        !$acc loop vector collapse(2) &
        !$acc reduction(+:buf01,buf02,buf03,buf04,buf05,buf06,buf07,buf08,buf09,buf10) &
        !$acc reduction(+:buf11,buf12,buf13,buf14,buf15,buf16,buf17,buf18,buf19,buf20) &
        !$acc reduction(+:buf21,buf22,buf23,buf24,buf25,buf26,buf27,buf28,buf29,buf30) &
        !$acc reduction(+:buf31,buf32,buf33,buf34,buf35,buf36,buf37,buf38)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            !
            ! terms needed for MKE (excl. mean pressure gradient)
            !
            buf01 = buf01 + u(i,j,k)                           ! cell center -> for MKE work by pressure gradient
            buf02 = buf02 + 0.5_rp*(u(i,j,k)+u(i,j,k+1))       ! cell edge   -> for MKE turbulent transport
            !
            buf03 = buf03 + (u(i,j,k+1)   -u(i,j,k)   )/dzc_g(k) ! cell edge -> for MKE viscous transport
            buf04 = buf04 + (u(i,j,k+1)**2-u(i,j,k)**2)/dzc_g(k) ! cell edge -> for MKE viscous transport
            !
            buf05 = buf05 + 0.25_rp*(u(i,j,k+1) + u(i  ,j,k))* &
                                    (w(i,j,k  ) + w(i+1,j,k))  ! cell edge   -> for MKE turbulent transport
            buf06 = buf06 + 0.25_rp*(u(i-1,j,k) + u(i,j,k  ))* &
                                    (w(i  ,j,k) + w(i,j,k-1))  ! cell center -> for MKE turbulent dissipation
            !
            buf07 = buf07 + 0.25_rp*( &
                                   ((u(i  ,j,k+1)-u(i  ,j,k  ))/dzc_g(k  )) + &
                                   ((u(i  ,j,k  )-u(i  ,j,k-1))/dzc_g(k-1)) + &
                                   ((u(i-1,j,k+1)-u(i-1,j,k  ))/dzc_g(k  )) + &
                                   ((u(i-1,j,k  )-u(i-1,j,k-1))/dzc_g(k-1)) &
                                  ) ! cell center for production and viscous dissipation
            !
            ! terms needed for uu (excl. those computed above)
            !
            ! for transport
            !
            buf08 = buf08 + 0.125_rp*(u(i,j,k+1) + u(i  ,j,k))**2* &
                                   (w(i,j,k  ) + w(i+1,j,k))   ! cell edge
            !
            ! for pressure-strain
            !
            buf09 = buf09 + p(i,j,k); q = q + 1 ! cell center
            buf10 = buf10 + (u(i,j,k)-u(i-1,j,k))/dl(1)*p(i,j,k) ! cell center
            !
            ! for dissipation
            !
            buf11 = buf11 + &
                                      ((u(i  ,j  ,k  )-u(i-1,j  ,k  ))/dl(1))**2 + &
                             0.25_rp*( &
                                      ((u(i  ,j+1,k  )-u(i  ,j  ,k  ))/dl(2))**2 + &
                                      ((u(i  ,j  ,k  )-u(i  ,j-1,k  ))/dl(2))**2 + &
                                      ((u(i-1,j+1,k  )-u(i-1,j  ,k  ))/dl(2))**2 + &
                                      ((u(i-1,j  ,k  )-u(i-1,j-1,k  ))/dl(2))**2 &
                                     ) + &
                             0.25_rp*( &
                                      ((u(i  ,j  ,k+1)-u(i  ,j  ,k  ))/dzc_g(k  ))**2 + &
                                      ((u(i  ,j  ,k  )-u(i  ,j  ,k-1))/dzc_g(k-1))**2 + &
                                      ((u(i-1,j  ,k+1)-u(i-1,j  ,k  ))/dzc_g(k  ))**2 + &
                                      ((u(i-1,j  ,k  )-u(i-1,j  ,k-1))/dzc_g(k-1))**2 &
                                     ) ! cell center
            !
            ! terms needed for vv (excl. those computed above)
            !
            ! for transport
            !
            buf12 = buf12  + (v(i,j,k+1)**2-v(i,j,k)**2)/dzc_g(k)    ! cell edge
            buf13 = buf13  + 0.125_rp*(v(i,j,k+1) + v(i,j  ,k))**2* &
                                      (w(i,j,k  ) + w(i,j+1,k))    ! cell edge
            !
            ! for pressure-strain
            !
            buf14 = buf14  + (v(i,j,k)-v(i,j-1,k))/dl(2)*p(i,j,k) ! cell center
            !
            ! for dissipation
            !
            buf15 = buf15  + &
                              0.25_rp*( &
                                       ((v(i+1,j  ,k  )-v(i  ,j  ,k  ))/dl(1))**2 + &
                                       ((v(i  ,j  ,k  )-v(i-1,j  ,k  ))/dl(1))**2 + &
                                       ((v(i+1,j-1,k  )-v(i  ,j-1,k  ))/dl(1))**2 + &
                                       ((v(i  ,j-1,k  )-v(i-1,j-1,k  ))/dl(1))**2 &
                                      ) + &
                                       ((v(i  ,j  ,k  )-v(i  ,j-1,k  ))/dl(2))**2 + &
                              0.25_rp*( &
                                       ((v(i  ,j  ,k+1)-v(i  ,j  ,k  ))/dzc_g(k  ))**2 + &
                                       ((v(i  ,j  ,k  )-v(i  ,j  ,k-1))/dzc_g(k-1))**2 + &
                                       ((v(i  ,j-1,k+1)-v(i  ,j-1,k  ))/dzc_g(k  ))**2 + &
                                       ((v(i  ,j-1,k  )-v(i  ,j-1,k-1))/dzc_g(k-1))**2 &
                                      ) ! cell center
            !
            ! terms needed for ww (excl. those computed above)
            !
            ! for transport
            !
            buf16 = buf16 + 0.5_rp*( (w(i,j,k+1)**2-w(i,j,k  )**2)/dzf_g(k+1) + &
                                     (w(i,j,k  )**2-w(i,j,k-1)**2)/dzf_g(k  ) ) ! cell edge
            buf17 = buf17 + w(i,j,k)**3                                       ! cell edge
            buf18 = buf18 + w(i,j,k)*0.5_rp*(p(i,j,k+1)+p(i,j,k))             ! cell edge
            !
            ! for pressure-strain
            !
            buf19 = buf19 + (w(i,j,k)-w(i,j,k-1))/dzf_g(k)*p(i,j,k)
            !
            ! for dissipation
            !
            buf20 = buf20 + &
                             0.25_rp*( &
                                      ((w(i+1,j,k  )-w(i  ,j,k  ))/dl(1))**2 + &
                                      ((w(i  ,j,k  )-w(i-1,j,k  ))/dl(1))**2 + &
                                      ((w(i+1,j,k-1)-w(i  ,j,k-1))/dl(1))**2 + &
                                      ((w(i  ,j,k-1)-w(i-1,j,k-1))/dl(1))**2 &
                                     ) + &
                             0.25_rp*( &
                                      ((w(i,j+1,k  )-w(i,j  ,k  ))/dl(2))**2 + &
                                      ((w(i,j  ,k  )-w(i,j-1,k  ))/dl(2))**2 + &
                                      ((w(i,j+1,k-1)-w(i,j  ,k-1))/dl(2))**2 + &
                                      ((w(i,j  ,k-1)-w(i,j-1,k-1))/dl(2))**2 &
                                     ) + &
                                      ((w(i,j,k)-w(i,j,k-1))/dzf_g(k))**2
            !
            ! terms needed for uw (excl. those computed above)
            !
            ! for production
            !
            buf21 = buf21 + 0.5_rp*(w(i,j,k)**2+w(i,j,k-1)**2)
            !
            ! for transport
            !
            buf22 = buf22 + ( 0.25_rp*( w(i,j,k) + w(i,j,k+1) + w(i+1,j,k+1) + w(i+1,j,k) )*u(i,j,k+1) - &
                              0.25_rp*( w(i,j,k) + w(i,j,k-1) + w(i+1,j,k-1) + w(i+1,j,k) )*u(i,j,k  ) )/dzc_g(k) ! cell edge
            buf23 = buf23 + w(i,j,k)**2 ! cell edge
            buf24 = buf24  + 0.125_rp*(u(i,j,k+1) + u(i  ,j,k))* &
                                      (w(i,j,k  ) + w(i+1,j,k))**2 ! cell edge
            buf25 = buf25 + 0.5_rp*(p(i,j,k+1)+p(i,j,k)) ! cell edge
            buf26 = buf26 + 0.25_rp*( u(i,j,k) + u(i,j,k+1) + u(i-1,j,k+1) + u(i-1,j,k) )* &
                            0.50_rp*( p(i,j,k+1) + p(i,j,k) ) ! cell edge
            !
            ! for pressure-strain
            !
            buf27 = buf27  + &
                              0.25_rp*( &
                                       (u(i  ,j,k+1)-u(i  ,j,k  ))/dzc_g(k  ) + &
                                       (u(i  ,j,k  )-u(i  ,j,k-1))/dzc_g(k-1) + &
                                       (u(i-1,j,k+1)-u(i-1,j,k  ))/dzc_g(k  ) + &
                                       (u(i-1,j,k  )-u(i-1,j,k-1))/dzc_g(k-1) &
                                      )*p(i,j,k) + &
                              0.25_rp*( &
                                       (w(i+1,j,k  )-w(i  ,j,k  ))/dl(1) + &
                                       (w(i  ,j,k  )-w(i-1,j,k  ))/dl(1) + &
                                       (w(i+1,j,k-1)-w(i  ,j,k-1))/dl(1) + &
                                       (w(i  ,j,k-1)-w(i-1,j,k-1))/dl(1) &
                                      )*p(i,j,k)
            !
            ! for dissipation
            !
            buf28 = buf28  + &
                                       (u(i  ,j,k  )-u(i-1,j,k  ))/dl(1)* &
                              0.25_rp*( &
                                       (w(i+1,j,k  )-w(i  ,j,k  ))/dl(1) + &
                                       (w(i  ,j,k  )-w(i-1,j,k  ))/dl(1) + &
                                       (w(i+1,j,k-1)-w(i  ,j,k-1))/dl(1) + &
                                       (w(i  ,j,k-1)-w(i-1,j,k-1))/dl(1) &
                                      ) + &
                               0.25_rp*( &
                                        (u(i  ,j+1,k)-u(i  ,j  ,k))/dl(2) + &
                                        (u(i  ,j  ,k)-u(i  ,j-1,k))/dl(2) + &
                                        (u(i-1,j+1,k)-u(i-1,j  ,k))/dl(2) + &
                                        (u(i-1,j  ,k)-u(i-1,j-1,k))/dl(2) &
                                       ) * &
                               0.25_rp*( &
                                        (w(i,j+1,k  )-w(i,j  ,k  ))/dl(2) + &
                                        (w(i,j  ,k  )-w(i,j-1,k  ))/dl(2) + &
                                        (w(i,j+1,k-1)-w(i,j  ,k-1))/dl(2) + &
                                        (w(i,j  ,k-1)-w(i,j-1,k-1))/dl(2) &
                                       ) + &
                                      0.25_rp*( &
                                        (u(i  ,j,k+1)-u(i  ,j,k  ))/dzc_g(k  ) + &
                                        (u(i  ,j,k  )-u(i  ,j,k-1))/dzc_g(k-1) + &
                                        (u(i-1,j,k+1)-u(i-1,j,k  ))/dzc_g(k  ) + &
                                        (u(i-1,j,k  )-u(i-1,j,k-1))/dzc_g(k-1) &
                                       )* &
                                        (w(i  ,j,k  )-w(i  ,j,k-1))/dzf_g(k)
            !
            ! split dissipation contributions
            !
            ! MKE
            !
            buf29 = buf29 + ((u(i  ,j  ,k+1)-u(i  ,j  ,k  ))/dzc_g(k))
            !
            ! uu
            !
            buf30 = buf30 + ((u(i  ,j  ,k  )-u(i-1,j  ,k  ))/dl(1)   )**2
            buf31 = buf31 + ((u(i  ,j+1,k  )-u(i  ,j  ,k  ))/dl(2)   )**2
            buf32 = buf32 + ((u(i  ,j  ,k+1)-u(i  ,j  ,k  ))/dzc_g(k))**2
            !
            ! vv
            !
            buf33 = buf33 + ((v(i+1,j  ,k  )-v(i  ,j  ,k  ))/dl(1)   )**2
            buf34 = buf34 + ((v(i  ,j  ,k  )-v(i  ,j-1,k  ))/dl(2)   )**2
            buf35 = buf35 + ((v(i  ,j  ,k+1)-v(i  ,j  ,k  ))/dzc_g(k))**2
            !
            ! ww
            !
            buf36 = buf36 + ((w(i+1,j  ,k  )-w(i  ,j  ,k  ))/dl(1)   )**2
            buf37 = buf37 + ((w(i  ,j+1,k  )-w(i  ,j  ,k  ))/dl(2)   )**2
            buf38 = buf38 + ((w(i  ,j  ,k  )-w(i  ,j  ,k-1))/dzf_g(k))**2
          end do
        end do
        buf( 1,k) = buf01*grid_area_ratio
        buf( 2,k) = buf02*grid_area_ratio
        buf( 3,k) = buf03*grid_area_ratio
        buf( 4,k) = buf04*grid_area_ratio
        buf( 5,k) = buf05*grid_area_ratio
        buf( 6,k) = buf06*grid_area_ratio
        buf( 7,k) = buf07*grid_area_ratio
        buf( 8,k) = buf08*grid_area_ratio
        buf( 9,k) = buf09*grid_area_ratio
        buf(10,k) = buf10*grid_area_ratio
        buf(11,k) = buf11*grid_area_ratio
        buf(12,k) = buf12*grid_area_ratio
        buf(13,k) = buf13*grid_area_ratio
        buf(14,k) = buf14*grid_area_ratio
        buf(15,k) = buf15*grid_area_ratio
        buf(16,k) = buf16*grid_area_ratio
        buf(17,k) = buf17*grid_area_ratio
        buf(18,k) = buf18*grid_area_ratio
        buf(19,k) = buf19*grid_area_ratio
        buf(20,k) = buf20*grid_area_ratio
        buf(21,k) = buf21*grid_area_ratio
        buf(22,k) = buf22*grid_area_ratio
        buf(23,k) = buf23*grid_area_ratio
        buf(24,k) = buf24*grid_area_ratio
        buf(25,k) = buf25*grid_area_ratio
        buf(26,k) = buf26*grid_area_ratio
        buf(27,k) = buf27*grid_area_ratio
        buf(28,k) = buf28*grid_area_ratio
        buf(29,k) = buf29*grid_area_ratio
        buf(30,k) = buf30*grid_area_ratio
        buf(31,k) = buf31*grid_area_ratio
        buf(32,k) = buf32*grid_area_ratio
        buf(33,k) = buf33*grid_area_ratio
        buf(34,k) = buf34*grid_area_ratio
        buf(35,k) = buf35*grid_area_ratio
        buf(36,k) = buf36*grid_area_ratio
        buf(37,k) = buf37*grid_area_ratio
        buf(38,k) = buf38*grid_area_ratio
      end do
      !$acc update self(buf) async(1)
      !$acc wait(1)
      call mpi_allreduce(MPI_IN_PLACE,buf(1,1),size(buf),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        nvars = 38
        write(cfmt,'(A)') '(*(es24.16e3,1x))'
        open(newunit=iunit,file=fname//'_reystr_budget.out')
        do k=1,nn
          write(iunit,trim(cfmt)) zc_g(k),zf_g(k),(buf(i,k),i=1,nvars),dzc_g(k),dzf_g(k)
        end do
        close(iunit)
        open(10,file=trim(fname)//'_reystr_budget.bin',access='stream')
        write(10) buf(1:nvars,1:nn)
        close(10)
      end if
      block
      use mod_params, only:dx,dy
      !$acc parallel loop gang default(present) async(1) &
      !$acc private(buf01,buf02,buf03,buf04,buf05,buf06) &
      !$acc private(div)
      do k=lo(3),hi(3)
        buf01 = 0._rp
        buf02 = 0._rp
        buf03 = 0._rp
        buf04 = 0._rp
        buf05 = 0._rp
        buf06 = 0._rp
        !$acc loop vector collapse(2) reduction(+:buf02,buf03,buf05,buf06) reduction(max:buf01,buf04)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            !
            ! velocity
            !
            div = (w(i,j,k)-w(i,j,k-1))/dzf_g(k) + &
                  (v(i,j,k)-v(i,j-1,k))/dy       + &
                  (u(i,j,k)-u(i-1,j,k))/dx
            buf01 = max(abs(div),buf01)
            buf02 = buf02 + abs(div)
            buf03 = buf03 + div
            buf04 = max(abs(div)*dzf_g(k),buf04)
            buf05 = buf05 + abs(div)*dzf_g(k)
            buf06 = buf06 + div*dzf_g(k)
          end do
        end do
        buf(1,k) = buf01
        buf(2,k) = buf02*grid_area_ratio
        buf(3,k) = buf03*grid_area_ratio
        buf(4,k) = buf04
        buf(5,k) = buf05*grid_area_ratio
        buf(6,k) = buf06*grid_area_ratio
      end do
      !$acc update self(buf) async(1)
      !$acc wait(1)
      call mpi_allreduce(MPI_IN_PLACE,buf(1,1),size(buf),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myid == 0) then
        nvars = 6
        write(cfmt,'(A)') '(*(es24.16e3,1x))'
        open(newunit=iunit,file=fname//'_leakage.out')
        do k=1,nn
          write(iunit,trim(cfmt)) zc_g(k),zf_g(k),(buf(i,k),i=1,nvars),dzc_g(k),dzf_g(k)
        end do
        close(iunit)
        open(10,file=trim(fname)//'_leakage.bin',access='stream')
        write(10) buf(1:nvars,1:nn)
        close(10)
      end if
      end block
    case(2)
    case(1)
    end select
    !$acc exit data delete(buf) async(1)
  end subroutine out1d_single_point_chan
end module mod_output
