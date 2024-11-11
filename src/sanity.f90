! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors.
! SPDX-FileCopyrightText: Modifications Copyright (c) 2023-2024 Maochao Xiao and the CaLES contributors.
! SPDX-License-Identifier: MIT
!
! -
module mod_sanity
  use, intrinsic :: iso_c_binding, only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound     , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv    , only: chkdiv
  use mod_common_mpi, only: myid,ierr,ipencil_axis
  use mod_correc    , only: correc
  use mod_fft       , only: fftend
  use mod_fillps    , only: fillps
  use mod_initflow  , only: add_noise
  use mod_initmpi   , only: initmpi
  use mod_param     , only: small
#if !defined(_OPENACC)
  use mod_solver    , only: solver
#else
  use mod_solver_gpu, only: solver => solver_gpu
#endif
  use mod_wallmodel    , only: updt_wallmodelbc
  use mod_typedef   , only: bound
  use mod_precision
  implicit none
  private
  public test_sanity_input
  contains
  subroutine test_sanity_input(ng,dims,sgstype,stop_type,cbcvel,cbcpre,cbcsgs,bcvel,bcpre,bcsgs, &
                               n,is_bound,lwm,l,zc,dl,h,is_forced)
    !
    ! performs some a priori checks of the input files before the calculation starts
    !
    implicit none
    integer , intent(in), dimension(3) :: ng
    integer , intent(in), dimension(2) :: dims
    character(len=*), intent(in) :: sgstype
    logical , intent(in), dimension(3) :: stop_type
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3)   :: cbcpre,cbcsgs
    real(rp), intent(in), dimension(0:1,3,3) :: bcvel
    real(rp), intent(in), dimension(0:1,3) :: bcpre,bcsgs
    integer, intent(in), dimension(3) :: n
    logical, intent(in), dimension(0:1,3) :: is_bound
    integer, intent(in), dimension(0:1,3) :: lwm
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in) :: h
    logical, intent(in), dimension(3) :: is_forced
    logical :: passed
    !
    call chk_dims(ng,dims,cbcvel,sgstype,passed); if(.not.passed) call abortit
    call chk_stop_type(stop_type,passed); if(.not.passed) call abortit
    call chk_bc(cbcvel,cbcpre,cbcsgs,bcvel,bcpre,bcsgs,n,is_bound,lwm,l,zc,dl,h,passed); if(.not.passed) call abortit
    call chk_forcing(cbcpre,is_forced,passed); if(.not.passed) call abortit
#if defined(_IMPDIFF_1D) && !defined(_IMPDIFF)
    if(myid == 0)  print*, 'ERROR: `_IMPDIFF_1D` cpp macro requires building with `_IMPDIFF` too.'; call abortit
#endif
#if defined(_IMPDIFF_1D) && !defined(_DECOMP_Z)
    if(myid == 0)  print*, 'WARNING: a run with implicit Z diffusion (`_IMPDIFF_1D`) is much more efficient &
                                   & when combined with a Z-pencils parallelization (`_DECOMP_Z`).'
#endif
  end subroutine test_sanity_input
  !
  subroutine chk_stop_type(stop_type,passed)
    implicit none
    logical, intent(in), dimension(3) :: stop_type
    logical, intent(out) :: passed
    passed = .true.
    if(.not.any(stop_type(:))) then
      if(myid == 0) print*, 'ERROR: stopping criterion not chosen.'
      passed = .false.
    end if
  end subroutine chk_stop_type
  !
  subroutine chk_dims(ng,dims,cbcvel,sgstype,passed)
    implicit none
    integer, intent(in), dimension(3) :: ng
    integer, intent(in), dimension(2) :: dims
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=*), intent(in) :: sgstype
    logical, intent(out) :: passed
    integer, dimension(2) :: ii
    logical :: passed_loc
    character(len=2) :: bc01v
    integer :: idir,ivel,i
    passed = .true.
    ii = pack([1,2,3],[1,2,3] /= ipencil_axis)
    passed_loc = all(dims(:)<=ng(ii)).and.all(dims(:)>=1)
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: 1 <= dims(:) <= [itot,jtot], or [itot,ktot], or [jtot ktot] depending on the decomposition.'
    passed = passed.and.passed_loc
    !
    if(trim(sgstype)=='smag') then
      passed_loc = .true.
      do i = 1,2
        idir = ii(i)
        ivel = ii(i)
        bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
        if(bc01v == 'DD') then
          passed_loc = passed_loc.and.(dims(i)<=2)
        end if
      end do
      if(myid == 0.and.(.not.passed_loc)) &
        print*, 'ERROR: more than two subdomains between two opposite walls.'
      passed = passed.and.passed_loc
    end if
    !
  end subroutine chk_dims
  !
  subroutine chk_bc(cbcvel,cbcpre,cbcsgs,bcvel,bcpre,bcsgs,n,is_bound,lwm,l,zc,dl,h,passed)
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3) :: cbcpre,cbcsgs
    real(rp)        , intent(in), dimension(0:1,3,3) :: bcvel
    real(rp)        , intent(in), dimension(0:1,3) :: bcpre,bcsgs
    integer, intent(in), dimension(3) :: n
    logical, intent(in), dimension(0:1,3) :: is_bound
    integer, intent(in), dimension(0:1,3) :: lwm
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(0:) :: zc
    real(rp), intent(in) :: h
    logical, intent(out) :: passed
    character(len=2) :: bc01v,bc01p,bc01s
    integer :: i,ivel,idir
    logical :: passed_loc
    passed = .true.
    !
    ! check validity of pressure and velocity BCs
    !
    passed_loc = .true.
    do ivel = 1,3
      do idir=1,3
        bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
        passed_loc = passed_loc.and.( (bc01v == 'PP').or. &
                                      (bc01v == 'ND').or. &
                                      (bc01v == 'DN').or. &
                                      (bc01v == 'NN').or. &
                                      (bc01v == 'DD') )
      end do
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: velocity BCs not valid.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01p == 'PP').or. &
                                    (bc01p == 'ND').or. &
                                    (bc01p == 'DN').or. &
                                    (bc01p == 'NN').or. &
                                    (bc01p == 'DD') )
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: pressure BCs not valid.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      ivel = idir
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01v == 'PP'.and.bc01p == 'PP').or. &
                                    (bc01v == 'ND'.and.bc01p == 'DN').or. &
                                    (bc01v == 'DN'.and.bc01p == 'ND').or. &
                                    (bc01v == 'DD'.and.bc01p == 'NN').or. &
                                    (bc01v == 'NN'.and.bc01p == 'DD') )
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: velocity and pressure BCs not compatible.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      bc01s = cbcsgs(0,idir)//cbcsgs(1,idir)
      passed_loc = passed_loc.and.( (bc01s == 'PP').or. &
                                    (bc01s == 'ND').or. &
                                    (bc01s == 'DN').or. &
                                    (bc01s == 'DD').or. &
                                    (bc01s == 'NN') )
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: sgs BCs not valid.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      ivel = idir
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      bc01s = cbcsgs(0,idir)//cbcsgs(1,idir)
      passed_loc = passed_loc.and.( (bc01v == 'PP'.and.bc01s == 'PP').or. &
                                    (bc01v == 'ND'.and.bc01s == 'DD').or. &
                                    (bc01v == 'DN'.and.bc01s == 'DD').or. &
                                    (bc01v == 'DD'.and.bc01s == 'DD').or. &
                                    (bc01v == 'NN'.and.bc01s == 'DD') )
    end do
    if(myid == 0.and.(.not.passed_loc)) print*, 'ERROR: velocity and sgs BCs not compatible.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,2
      passed_loc = passed_loc.and.((bcpre(0,idir) == 0.).and.(bcpre(1,idir) == 0.))
    end do
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: pressure BCs in directions x and y must be homogeneous (value = 0.).'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir = 1,3
      do i = 0,1
        if(lwm(i,idir)/=0) then
          do ivel = 1,3
            passed_loc = passed_loc.and.cbcvel(i,idir,ivel)=='D'
          end do
        end if
      end do
    end do
    if(myid == 0.and.(.not.passed_loc)) &
    print*, 'ERROR: wall model BCs must be Dirichlet.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    if(is_bound(0,1).and.lwm(0,1)/=0) passed_loc = passed_loc.and.(h>0.5_rp*dl(1) .and.h<(n(1)-0.5_rp)*dl(1))
    if(is_bound(1,1).and.lwm(1,1)/=0) passed_loc = passed_loc.and.(h>0.5_rp*dl(1) .and.h<(n(1)-0.5_rp)*dl(1))
    if(is_bound(0,2).and.lwm(0,2)/=0) passed_loc = passed_loc.and.(h>0.5_rp*dl(2) .and.h<(n(2)-0.5_rp)*dl(2))
    if(is_bound(1,2).and.lwm(1,2)/=0) passed_loc = passed_loc.and.(h>0.5_rp*dl(2) .and.h<(n(2)-0.5_rp)*dl(2))
    if(is_bound(0,3).and.lwm(0,3)/=0) passed_loc = passed_loc.and.(h>zc(1)        .and.h<zc(n(3))  )
    if(is_bound(1,3).and.lwm(1,3)/=0) passed_loc = passed_loc.and.(h>l(3)-zc(n(3)).and.h<l(3)-zc(1))
    if(.not.passed_loc) print*, 'ERROR: invalid wall model height.'
    passed = passed.and.passed_loc
    !
#if defined(_IMPDIFF) && !defined(_IMPDIFF_1D)
    ! implicit: bcs in x and y must be 'D0' or 'P',
    ! so implicit mode cannot be used for WMLES of square ducts
    passed_loc = .true.
    do ivel = 1,3
      do idir=1,2
        bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
        passed_loc = passed_loc.and.(bc01v /= 'NN')
      end do
    end do
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: Neumann-Neumann velocity BCs with implicit diffusion currently not supported in x and y; only in z.'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do ivel = 1,3
      do idir=1,2
        passed_loc = passed_loc.and.((bcvel(0,idir,ivel) == 0.).and.(bcvel(1,idir,ivel) == 0.))
      end do
    end do
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: velocity BCs with implicit diffusion in directions x and y must be homogeneous (value = 0.).'
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    passed_loc = passed_loc.and.( lwm(0,1)==0.and.lwm(1,1)==0.and. &
                                  lwm(0,2)==0.and.lwm(1,2)==0 )
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: wall model BCs cannot be used in x and y directions when 3D implicit diffusion is applied.'
    passed = passed.and.passed_loc
    !
#endif
#if defined(_OPENACC)
    do idir=1,2
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and..not.( (bc01p == 'DN').or. &
                                         (bc01p == 'ND') )
    end do
    if(myid == 0.and.(.not.passed_loc)) &
      print*, 'ERROR: pressure BCs "ND" or "DN" along x or y not implemented on GPUs yet.'
#endif
  end subroutine chk_bc
  !
  subroutine chk_forcing(cbcpre,is_forced,passed)
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbcpre
    logical         , intent(in), dimension(3) :: is_forced
    logical         , intent(out) :: passed
    integer :: idir
    passed = .true.
    !
    ! 1) check for compatibility between pressure BCs and flow forcing
    !
    do idir=1,3
      if(is_forced(idir)) then
        passed = passed.and.(cbcpre(0,idir)//cbcpre(1,idir) == 'PP')
      end if
    end do
    if(myid == 0.and.(.not.passed)) &
    print*, 'ERROR: Flow cannot be forced in a non-periodic direction; check the BCs and is_forced in `input.nml`.'
  end subroutine chk_forcing
  !
  subroutine abortit
    implicit none
    if(myid == 0) print*, ''
    if(myid == 0) print*, '*** Simulation aborted due to errors in the input file ***'
    if(myid == 0) print*, '    check `input.nml`.'
    call decomp_2d_finalize
    call MPI_FINALIZE(ierr)
    error stop
  end subroutine abortit
end module mod_sanity
