! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_initmpi
  use mpi
  use decomp_2d
  use mod_common_mpi, only: myid,ierr,halo,ipencil => ipencil_axis
  use mod_precision
  !@acc use openacc
  !@acc use cudecomp
  !@cuf use cudafor, only: cudaGetDeviceCount,cudaSetDevice
#if defined(_OPENACC)
  use mod_common_cudecomp, only: cudecomp_real_rp, &
                                 ch => handle,gd => gd_halo,gd_poi, &
                                 ap_x,ap_y,ap_z,ap_x_poi,ap_y_poi,ap_z_poi
  use mod_param, only: cudecomp_t_comm_backend     ,cudecomp_h_comm_backend    , &
                       cudecomp_is_t_comm_autotune ,cudecomp_is_h_comm_autotune, &
                       cudecomp_is_t_enable_nccl   ,cudecomp_is_h_enable_nccl  , &
                       cudecomp_is_t_enable_nvshmem,cudecomp_is_h_enable_nvshmem
#endif
  implicit none
  private
  public initmpi
  !@acc integer, parameter :: CUDECOMP_RANK_NULL = -1
  contains
  subroutine initmpi(ng,dims,sgstype,cbcvel,cbcpre,lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z,nb,is_bound)
    implicit none
    integer, intent(in   ), dimension(3) :: ng
    integer, intent(inout), dimension(2) :: dims
    character(len=*), intent(in) :: sgstype
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3) :: cbcpre
    integer, intent(out), dimension(3    ) :: lo,hi,n,n_x_fft,n_y_fft,lo_z,hi_z,n_z
    integer, intent(out), dimension(0:1,3) :: nb
    logical, intent(out), dimension(0:1,3) :: is_bound
    logical, dimension(3) :: periods
    integer :: l,ipencil_t(2),nproc
#if defined(_OPENACC)
    integer(acc_device_kind) ::dev_type
    integer :: local_comm,mydev,ndev
    integer :: istat
    type(cudecompGridDescConfig)          :: conf,conf_poi
    type(cudecompGridDescAutotuneOptions) :: atune_conf
#else
    integer :: comm_cart
#endif
    !
#if !defined(_DECOMP_Y) && !defined(_DECOMP_Z)
    ipencil=1
#elif defined(_DECOMP_Y)
    ipencil=2
#elif defined(_DECOMP_Z)
    ipencil=3
#endif
    ipencil_t(:) = pack([1,2,3],[1,2,3] /= ipencil)
    is_bound(:,:) = .false.
    !
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
    if(product(dims)==0.and.nproc>=2.and.trim(sgstype)=='smag') then
      call calc_dims(cbcvel,sgstype,ipencil_t,nproc,dims)
      if(myid == 0) then
        print*, 'In auto-tuning mode......'
        print*, 'p_row x p_col',dims(1),dims(2)
      end if
    end if
    !
    periods(:) = .false.
    where(cbcpre(0,:)//cbcpre(1,:) == 'PP') periods(:) = .true.
    !
#if defined(_OPENACC)
    call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,local_comm,ierr)
    call MPI_COMM_RANK(local_comm,mydev,ierr)
    dev_type = acc_get_device_type()
#if 1
    istat = cudaGetDeviceCount(ndev)      ! may be tweaked with environment variable CUDA_VISIBLE_DEVICES
    mydev = mod(mydev,ndev)
    istat = cudaSetDevice(mydev)
    if(istat /= 0) print*,'MPI rank: ',myid,' error assigning GPU: ',mydev
#else
    ndev  = acc_get_num_devices(dev_type) ! may be tweaked with environment variable ACC_DEVICE_NUM
    mydev = mod(mydev,ndev)
#endif
    call acc_set_device_num(mydev,dev_type)
    call acc_init(dev_type)
    !
    istat = cudecompInit(ch,MPI_COMM_WORLD)
    !
    ! setup descriptor for the Poisson solver
    !
    istat = cudecompGridDescConfigSetDefaults(conf)
    conf%transpose_comm_backend = cudecomp_t_comm_backend
    conf%transpose_axis_contiguous(:) = [.true.,.true.,.false.]
    conf%gdims(:)      = [2*(ng(1)/2+1),2*(ng(2)/2+1),ng(3)]
    conf%pdims(:)      = dims(1:2)
    conf%gdims_dist(:) = ng(:)
    istat = cudecompGridDescAutotuneOptionsSetDefaults(atune_conf)
    if(rp == dp) then
      cudecomp_real_rp = CUDECOMP_DOUBLE
    else
      cudecomp_real_rp = CUDECOMP_FLOAT
    end if
    atune_conf%dtype = cudecomp_real_rp
    atune_conf%grid_mode = CUDECOMP_AUTOTUNE_GRID_TRANSPOSE
    atune_conf%autotune_transpose_backend = cudecomp_is_t_comm_autotune
    atune_conf%disable_nccl_backends    = .not.cudecomp_is_t_enable_nccl
    atune_conf%disable_nvshmem_backends = .not.cudecomp_is_t_enable_nvshmem
    istat = cudecompGridDescCreate(ch,gd_poi,conf,atune_conf)
    conf_poi = conf
    dims(:) = conf%pdims
    !
    ! setup descriptor for halo exchanges
    !
    istat = cudecompGridDescConfigSetDefaults(conf)
    conf%gdims(:)      = ng(:)
    conf%pdims(:)      = dims(1:2)
    conf%halo_comm_backend = cudecomp_h_comm_backend
    conf%transpose_axis_contiguous(:) = .false.
    istat = cudecompGridDescAutotuneOptionsSetDefaults(atune_conf)
    atune_conf%halo_extents(:) = 1
    atune_conf%halo_periods(:) = periods(:)
    atune_conf%dtype = cudecomp_real_rp
    atune_conf%autotune_halo_backend = cudecomp_is_h_comm_autotune
    atune_conf%disable_nccl_backends    = .not.cudecomp_is_t_enable_nccl
    atune_conf%disable_nvshmem_backends = .not.cudecomp_is_t_enable_nvshmem
    if(all(conf_poi%transpose_comm_backend /= [CUDECOMP_TRANSPOSE_COMM_NVSHMEM,CUDECOMP_TRANSPOSE_COMM_NVSHMEM_PL])) then
      !
      ! disable NVSHMEM halo backend autotuning when NVSHMEM is NOT used for transposes
      !
      atune_conf%disable_nvshmem_backends = .true.
    end if
    istat = cudecompGridDescCreate(ch,gd,conf,atune_conf)
#endif
    call decomp_2d_init(ng(1),ng(2),ng(3),dims(1),dims(2),periods)
#if defined(_OPENACC)
    !
    ! fetch lo(:), hi(:), n(:) and n_z(:) from cuDecomp (should match the modified 2decomp one)
    !
    istat = cudecompGetPencilInfo(ch,gd    ,ap_x    ,1)
    istat = cudecompGetPencilInfo(ch,gd    ,ap_y    ,2)
    istat = cudecompGetPencilInfo(ch,gd    ,ap_z    ,3)
    istat = cudecompGetPencilInfo(ch,gd_poi,ap_x_poi,1)
    istat = cudecompGetPencilInfo(ch,gd_poi,ap_y_poi,2)
    istat = cudecompGetPencilInfo(ch,gd_poi,ap_z_poi,3)
    select case(ipencil)
    case(1)
      lo(:) = ap_x%lo(:)
      hi(:) = ap_x%hi(:)
    case(2)
      lo(:) = ap_y%lo(:)
      hi(:) = ap_y%hi(:)
    case(3)
      lo(:) = ap_z%lo(:)
      hi(:) = ap_z%hi(:)
    end select
    n(:)       = hi(:)-lo(:)+1
    n_x_fft(:) = ap_x_poi%shape(:)
    n_y_fft(:) = ap_y_poi%shape(:)
    lo_z(:)    = ap_z%lo(:)
    hi_z(:)    = ap_z%hi(:)
    n_z(:)     = ap_z%shape(:)
    nb(:,ipencil) = CUDECOMP_RANK_NULL
    associate(ip_t => ipencil_t)
      istat = cudecompGetShiftedRank(ch,gd,ipencil,ip_t(1),-1,periods(ip_t(1)),nb(0,ip_t(1)))
      istat = cudecompGetShiftedRank(ch,gd,ipencil,ip_t(1), 1,periods(ip_t(1)),nb(1,ip_t(1)))
      istat = cudecompGetShiftedRank(ch,gd,ipencil,ip_t(2),-1,periods(ip_t(2)),nb(0,ip_t(2)))
      istat = cudecompGetShiftedRank(ch,gd,ipencil,ip_t(2), 1,periods(ip_t(2)),nb(1,ip_t(2)))
    end associate
    where(nb(:,:) == CUDECOMP_RANK_NULL) is_bound(:,:) = .true.
#else
    select case(ipencil)
    case(1)
      lo(:) = xstart(:)
      hi(:) = xend(:)
      comm_cart = DECOMP_2D_COMM_CART_X
    case(2)
      lo(:) = ystart(:)
      hi(:) = yend(:)
      comm_cart = DECOMP_2D_COMM_CART_Y
    case(3)
      lo(:) = zstart(:)
      hi(:) = zend(:)
      comm_cart = DECOMP_2D_COMM_CART_Z
    end select
    n(:)       = hi(:)-lo(:)+1
    n_x_fft(:) = xsize(:)
    n_y_fft(:) = ysize(:)
    lo_z(:)    = zstart(:)
    hi_z(:)    = zend(:)
    n_z(:)     = zsize(:)
    do l=1,3
      call makehalo(l,1,n(:),halo(l))
    end do
    nb(:,ipencil) = MPI_PROC_NULL
    call MPI_CART_SHIFT(comm_cart,0,1,nb(0,ipencil_t(1)),nb(1,ipencil_t(1)),ierr)
    call MPI_CART_SHIFT(comm_cart,1,1,nb(0,ipencil_t(2)),nb(1,ipencil_t(2)),ierr)
    where(nb(:,:) == MPI_PROC_NULL) is_bound(:,:) = .true.
#endif
  end subroutine initmpi
  !
  subroutine makehalo(idir,nh,n,halo)
    implicit none
    integer, intent(in ) :: idir,nh
    integer, intent(in ), dimension(3) :: n
    integer, intent(out) :: halo
    integer, dimension(3) :: nn
    nn(:) = n(:) + 2*nh
    select case(idir)
    case(1)
      call MPI_TYPE_VECTOR(nn(2)*nn(3),nh            ,nn(1)            ,MPI_REAL_RP,halo,ierr)
    case(2)
      call MPI_TYPE_VECTOR(      nn(3),nh*nn(1)      ,nn(1)*nn(2)      ,MPI_REAL_RP,halo,ierr)
    case(3)
      call MPI_TYPE_VECTOR(          1,nh*nn(1)*nn(2),nn(1)*nn(2)*nn(3),MPI_REAL_RP,halo,ierr)
    end select
    call MPI_TYPE_COMMIT(halo,ierr)
    !number of blocks, number of elements in each block, displacement between the blocks
    !exchange (n(2)+2)*(n(3)+2) cells in the x direction
    !exchange (n(1)+2)*(n(3)+2) cells in the y direction
    !exchange (n(1)+2)*(n(2)+2) cells in the z direction
  end subroutine makehalo
  !
  subroutine calc_dims(cbcvel,sgstype,ipencil_t,nproc,dims)
    !
    ! calculate dims(1:2) to ensure <= 2 subdomains between two opposite walls.
    ! The splitting could be more general, but unnecessary for the current purpose.
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=*), intent(in) :: sgstype
    integer, intent(in),dimension(2) :: ipencil_t(2)
    integer, intent(in) :: nproc
    integer, intent(inout), dimension(2) :: dims
    character(len=2) :: bc01v
    integer :: i,idir,ivel
    !
    do i = 1,2
      idir = ipencil_t(i)
      ivel = ipencil_t(i)
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      if(bc01v == 'DD') then
        if(mod(nproc,2) == 1) then
          dims(i)   = 1
          dims(3-i) = nproc
        else
          dims(i)   = 2
          dims(3-i) = nproc/2
        end if
        exit
      end if
    end do
  end subroutine calc_dims
end module mod_initmpi
