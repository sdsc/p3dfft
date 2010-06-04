! This sample program illustrates the 
! use of P3DFFT library for highly scalable parallel 3D FFT. 
!
! This program initializes a 3D array with RANDOM NUMBERS, then 
! performs a forward FFT transform, and computes POWER SPECTRUM.
!
! The program expects 'stdin' file in the working directory, with 
! a single line of numbers : Nx,Ny,Nz,Ndim,Nrep. Here Nx,Ny,Nz
! are box dimensions, Ndim is the dimentionality of processor grid
! (1 or 2), and Nrep is the number of repititions. Optionally
! a file named 'dims' can also be provided to guide in the choice 
! of processor geometry in case of 2D decomposition. It should contain 
! two numbers in a line, with their product equal to the total number
! of tasks. Otherwise processor grid geometry is chosen automatically.
! For better performance, experiment with this setting, varying 
! iproc and jproc. In many cases, minimizing iproc gives best results. 
! Setting it to 1 corresponds to one-dimensional decomposition.
!
! If you have questions please contact Dmitry Pekurovsky, dmitry@sdsc.edu

      program fft3d_spec

      use p3dfft
      implicit none
      include 'mpif.h'

      integer i,n,nx,ny,nz
      integer m,x,y,z
      integer fstatus
      logical flg_inplace

      real(mytype), dimension(:,:,:),  allocatable :: BEG
      complex(mytype), dimension(:,:,:),  allocatable :: AEND

      integer(8) Ntot
      real(mytype) factor
      real(8) rtime1,rtime2
      real(8) gt(4,3),gtcomm(3),tc
      integer ierr,nu,ndim,dims(2),nproc,proc_id
      integer istart(3),iend(3),isize(3)
      integer fstart(3),fend(3),fsize(3)
      integer iproc,jproc,ng(3),kmax,k
      real(mytype), allocatable :: E(:)
      logical iex

      call MPI_INIT (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)

      timers = 0.0
      gt=0.0
      gtcomm=0.0

      if (proc_id.eq.0) then 
         open (unit=3,file='stdin',status='old', &
               access='sequential',form='formatted', iostat=fstatus)
         if (fstatus .eq. 0) then
            write(*, *) ' Reading from input file stdin'
         endif 
         ndim = 2

        read (3,*) nx, ny, nz, ndim
        write (*,*) "procs=",nproc," nx=",nx, &
                " ny=", ny," nz=", nz,"ndim=",ndim
        if(mytype .eq. 4) then
           print *,'Single precision version'
        else if(mytype .eq. 8) then
           print *,'Double precision version'
        endif
       endif

      call MPI_Bcast(nx,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(ny,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(nz,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(ndim,1, MPI_INTEGER,0,mpi_comm_world,ierr)

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

! Initialize P3DFFT

      call p3dfft_setup (dims,nx,ny,nz,.false.)

! Get dimensions for initial and outcoming arrays

      call get_dims(istart,iend,isize,1)
      call get_dims(fstart,fend,fsize,2)

! Allocate array for initial data

      allocate (BEG(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)), stat=ierr)
      if(ierr .ne. 0) then
         print *,'Error ',ierr,' allocating array BEG'
      endif

! Allocate the (complex) array for transformed data

      allocate (AEND(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)), stat=ierr)
      if(ierr .ne. 0) then
         print *,'Error ',ierr,' allocating array AEND'
      endif

! Initialize with random numbers
!
      do z=istart(3),iend(3)
         do y=istart(2),iend(2)
            do x=istart(1),iend(1)
               call random_number(BEG(x,y,z))
            enddo
         enddo
      enddo

!
! transform from physical space to wavenumber space
! (XgYiZj to XiYjZg)
! then transform back to physical space
! (XiYjZg to XgYiZj)
!

      factor = 1.0d0/(nx*ny)
      factor = factor / nz
      Ntot = fsize(1)*fsize(2)*fsize(3)

! Barrier for correct timing
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      rtime1 = - MPI_wtime()

! Do forward Fourier transform
      call p3dfft_ftran_r2c (BEG,AEND)
         
      rtime1 = rtime1 + MPI_wtime()
         
! Normalize
      call mult_array(AEND, Ntot,factor)

      kmax = sqrt(real(nx)*nx + ny*ny + nz*nz)*0.5 + 0.5
      allocate(E(0:kmax))
      ng(1) = nx
      ng(2) = ny
      ng(3) = nz

! Compute power spectrum
      call compute_spectrum(AEND,fstart,fsize,fend,ng,E,kmax,0)

      if(proc_id .eq. 0) then
         print *,'Power spectrum:'
         do k=0,kmax
             print *,k,E(k)
         enddo
      endif

! Free work space
      call p3dfft_clean

! Process timing statistics
      call MPI_Reduce(rtime1,rtime2,1,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      if (proc_id.eq.0) write(6,*)'proc_id, cpu time', &
         proc_id,rtime2

      call MPI_Reduce(timers,gt(1,1),4,mpi_real8,MPI_SUM,0, &
        MPI_COMM_WORLD,ierr)

      call MPI_Reduce(timers,gt(1,2),4,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      call MPI_Reduce(timers,gt(1,3),4,mpi_real8,MPI_MIN,0, &
        MPI_COMM_WORLD,ierr)

      tc = (timers(1)+timers(2)+timers(3)+timers(4))
      call MPI_Reduce(tc,gtcomm(1),1,mpi_real8,MPI_SUM,0, &
        MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tc,gtcomm(2),1,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tc,gtcomm(3),1,mpi_real8,MPI_MIN,0, &
        MPI_COMM_WORLD,ierr)

      gt(1:4,1) = gt(1:4,1) / dble(nproc)
      gtcomm(1) = gtcomm(1) / dble(nproc)

      if(proc_id .eq. 0) then
         do i=1,4
            print *,'timer',i,' (avg/max/min): ',gt(i,:)
         enddo
         print *,'Total comm (avg/max/min): ',gtcomm
      endif


      call MPI_FINALIZE (ierr)

      contains 
!=========================================================

      subroutine compute_spectrum(B,st,sz,en,ng,E,kmax,root)

      use p3dfft
      
      integer st(3),sz(3),en(3),ng(3),root
      complex(mytype) B(sz(1),sz(2),sz(3))
      real(mytype) E(0:kmax),el(0:kmax),k2
      integer kmax,nl(3),x,y,z,i,ik,kx,ky,kz,cx,cy,cz

      if(4.0 * kmax**2 .lt. (en(1)-1)**2 +(en(2)-1)**2+(en(3)-1)**2) then
         print *,'Error: kmax is less than it should be (',kmax,'),vs. ',en
         call abort
      endif

      do ik=0,kmax
         el(ik) = 0.0
         E(ik) = 0.0
      enddo

      do z=st(3),en(3)
         if(z .gt. ng(3)/2) then
            kz = ng(3) - z +1
         else
            kz = z-1
         endif

         do y=st(2),en(2)
            if(y .gt. ng(2)/2) then
               ky = ng(2) - y +1
            else
               ky = y-1
            endif

            do x=st(1),en(1)
               kx = x-1

               k2 = kx**2 + ky**2 + kz**2
               ik = sqrt(k2) + 0.5
               el(ik) = el(ik) + k2 * real(B(x-st(1)+1,y-st(2)+1,z-st(3)+1) * &
                   conjg(B(x-st(1)+1,y-st(2)+1,z-st(3)+1))) 

            enddo
         enddo
      enddo

      call mpi_reduce(el,E,kmax+1,mpireal,MPI_SUM,root,MPI_COMM_WORLD,ierr)

      return 
      end subroutine


!=========================================================
      subroutine mult_array(X,nar,f)

      use p3dfft

      integer(8) nar,i
      complex(mytype) X(nar)
      real(mytype) f

      do i=1,nar
         X(i) = X(i) * f
      enddo

      return
      end subroutine

      end
