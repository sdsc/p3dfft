
! This sample program illustrates the
! use of P3DFFT library for highly scalable parallel 3D FFT.
!
! This program initializes a 3D array with a 3D sine wave, then
! performs forward cosine transform, differentiation by Chebyshev method, 
! backward transform, and checks that
! the results are correct. It can be used both as a correctness
! test and for timing the library functions.
!
! Processor grid in this test is chosen to be square or close to square.
! For better performance, experiment with this setting, varying
! iproc and jproc. In many cases, minimizing iproc gives best results.
! Setting it to 1 corresponds to one-dimensional decomposition, which
! is the best in case it's feasible.
!
! To build this program, use one of the provided makefiles, modifying
! it as needed for your compilers, flags and library locations.
! Be sure to link the program with both the P3DFFT library and the underlying
! 1D FFT library such as ESSL or FFTW.
!
! If you have questions please contact Dmitry Pekurovsky, dmitry@sdsc.edu

      program cheby3d

      use p3dfft
      implicit none
      include 'mpif.h'

      integer i,j,k,n,m,x,y,z, msize(3)
      integer fstatus
      integer*8 mem_rsize1d,mem_csize1d
      real(mytype) pi
      logical iex
      real(8) timer,rtime2

      complex(mytype), dimension(:,:,:), allocatable, target::mem

      integer*8 Ntot
      integer ierr,nu,ndim,dims(2),nproc,proc_id
      integer iproc,jproc
      integer memsize(3)

      integer nx,ny,nz
      common /grid/ nx,ny,nz

      integer istart(3),iend(3),isize(3)
      integer fstart(3),fend(3),fsize(3)
      common /p3dfftsets/ istart,iend,isize,fstart,fend,fsize,proc_id


      real(mytype), dimension(:), allocatable :: coordX, coordY, coordZ
      real(mytype) :: Lx,Ly,Lz
      common /Lxyz/ Lx,Ly,Lz

      call MPI_INIT (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)

      pi =  4.d0 *datan(1.d0)
         n    = 4
         ndim = 2
         nx   = 32; Lx = 2.d0*pi
         ny   = 32; Ly = 4.d0*pi
         nz   = 33; Lz = 2.d0

      if (proc_id.eq.0) then 
	print *,'P3DFFT test, Chebyshev'
         open (unit=3,file='stdin',status='old',  &
              access='sequential',form='formatted', iostat=fstatus)
         if (fstatus .eq. 0) then
            write(*, *) ' Reading from input file stdin'
         endif 
         ndim = 2

        read (3,*) nx, ny, nz, ndim,n
        write (*,*) "procs=",nproc," nx=",nx,  &
               " ny=", ny," nz=", nz,"ndim=",ndim," repeat=", n
        if(mytype .eq. 4) then
           print *,'Single precision version'
        else if(mytype .eq. 8) then
           print *,'Double precision version'
        endif
       endif

      call MPI_Bcast(nx,1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ny,1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(nz,1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(n,1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ndim,1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      	allocate(coordX(nx))
      	do i=1,nx
      		coordX(i) = Lx*(2*i-1)/dble(2*nx)
      	enddo
      	allocate(coordY(ny))
      	do j=1,ny
      		coordY(j) = Ly*(2*j-1)/dble(2*ny)
      	enddo
      	allocate(coordZ(nz))
      	do k=1,nz
#ifndef SINGLE_PREC
      		coordZ(k) = dcos( pi*(k-1)/dble(nz-1) ) *2.d0/Lz
#else
      		coordZ(k) = cos( pi*(k-1)/dble(nz-1) ) *2.d0/Lz
#endif
      	enddo
      	if(mod(nz,2) .eq. 0) then
      		if(proc_id .eq. 0) write(*,*) 'ERROR: Chebyshev expects nz to be odd!'
      		stop
      	endif

!    nproc is devided into a iproc x jproc stencle
!
      if(ndim .eq. 1) then
         dims(1) = 1
         dims(2) = nproc
      else if(ndim .eq. 2) then
	inquire(file='dims',exist=iex)
	if (iex) then
           if (proc_id.eq.0) print *, 'Reading proc. grid from file dims'
           open (999,file='dims')
           read (999,*) dims(1), dims(2)
           close (999)
           if(dims(1) * dims(2) .ne. nproc) then
              dims(2) = nproc / dims(1)
           endif
	else
           if (proc_id.eq.0) print *, 'Creating proc. grid with mpi_dims_create'
           dims(1) = 0
           dims(2) = 0
           call MPI_Dims_create(nproc,2,dims,ierr)
           if(dims(1) .gt. dims(2)) then
              dims(1) = dims(2)
              dims(2) = nproc / dims(1)
           endif
        endif
      endif

      iproc = dims(1)
      jproc = dims(2)

      if(proc_id .eq. 0) then
         print *,'Using processor grid ',iproc,' x ',jproc
      endif

      call p3dfft_setup (dims,nx,ny,nz,MPI_COMM_WORLD,nx,ny,nz,.false.,memsize)
      call p3dfft_get_dims(istart,iend,isize,1)
      call p3dfft_get_dims(fstart,fend,fsize,2)

      mem_rsize1d = int(isize(1),8) *int(isize(2),8) *int(isize(3),8)
      mem_csize1d = int(fsize(1),8) *int(fsize(2),8) *int(fsize(3),8) *2
      if(mem_rsize1d .gt. mem_csize1d) then
      	msize = isize
      else
      	msize = (/2*fsize(1),fsize(2),fsize(3)/)
      endif
      allocate(mem(msize(1),msize(2),msize(3)),stat=ierr)
      if(ierr .ne. 0) then
         print *,'Error ',ierr,' allocating array mem3d'
      endif

      timer = 0.
      do i=1,n
      	call cheby_test(mem,mem,coordX,coordY,coordZ,timer)
      enddo

      call MPI_Reduce(timer,rtime2,1,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      if (proc_id.eq.0) write(6,*)'time per loop', rtime2/dble(n)

      call MPI_FINALIZE (ierr)

      stop

      end

!=========================================================

!     !========================================
!     !
!     !		cheby test
!     !			calc partial(X)/partial(Z) = derivative in chebyshev-direction
!     !
!     !========================================
      subroutine cheby_test(rmem,cmem,coordX,coordY,coordZ,timer)
      	use p3dfft
      	use mpi
      	implicit none

!    !	common block
      	integer proc_id
      	integer istart(3),iend(3),isize(3)
      	integer fstart(3),fend(3),fsize(3)
      	common /p3dfftsets/istart,iend,isize,fstart,fend,fsize,proc_id
        real(8) timer

      	real(mytype) :: Lx,Ly,Lz
      	common /Lxyz/ Lx,Ly,Lz

      	integer nx,ny,nz
      	common /grid/ nx,ny,nz

!     !	function args
      	real(mytype), dimension(istart(1):iend(1), &
                               istart(2):iend(2), &
                               istart(3):iend(3)) ::	rmem
      	complex(mytype), dimension(fstart(1):fend(1), &
                               fstart(2):fend(2),  &
                              fstart(3):fend(3)) ::	cmem
      	real(mytype), intent(in) :: coordX(nx), coordY(ny), coordZ(nz)

      	double precision  :: delta,prec
      	double precision  :: maxdelta1, maxdelta01
      	integer            :: i,j,k,ierr

      	do k=istart(3),iend(3)
#ifndef SINGLE_PREC
           rmem(:,:,k) = dsin(coordZ(k))
#else
           rmem(:,:,k) = sin(coordZ(k))
#endif
        enddo

         call MPI_Barrier(MPI_COMM_WORLD,ierr)
        timer = timer - MPI_Wtime() 

        call p3dfft_cheby(rmem,cmem,Lz)

        timer = timer + MPI_Wtime()

        if(proc_id .eq. 0) then
           print *,'After cheby transform'
        endif
        call print_buf(cmem,fsize(1),fsize(2),fsize(3))

#ifdef STRIDE1

        do i=fstart(3),fend(3)
           do j=fstart(2),fend(2)
              cmem(1,j,i) = cmem(1,j,i) * 2.0d0
              cmem(nz,j,i) = cmem(nz,j,i) * 2.0d0
           enddo
        enddo
#else
        do j=fstart(2),fend(2)
           do i=fstart(1),fend(1)
              cmem(i,j,1) = cmem(i,j,1) * 2.0d0
              cmem(i,j,nz) = cmem(i,j,nz) * 2.0d0
           enddo
        enddo
#endif
        cmem = cmem * 0.5d0

      	call p3dfft_btran_c2r(cmem, rmem,'cff')

      	maxdelta1=0.d0
        do k=istart(3),iend(3)
      	  do j=istart(2),iend(2)
             do i=istart(1),iend(1)
#ifndef SINGLE_PREC
     		  delta = rmem(i,j,k) -dcos(coordZ(k))
#else
     		  delta = rmem(i,j,k) -cos(coordZ(k))
#endif
     		  if(abs(delta) .gt. maxdelta1) maxdelta1 = abs(delta)
     		enddo
     	  enddo
     	enddo

      	call MPI_reduce(maxdelta1,maxdelta01,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      	if(proc_id .eq. 0) then
           write(*,*) '  FFT(chebyshev) partial3D error in Z : ',maxdelta01
         if(mytype .eq. 8) then
            prec = 1e-14
         else
            prec = 1e-5
         endif
         if(maxdelta01 .gt. prec * Nx * Ny * Nz*0.25) then
            print *,'Results are incorrect'
         else
            print *,'Results are correct'
         endif
      endif

      return
      end subroutine cheby_test



