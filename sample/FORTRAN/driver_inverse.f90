! This sample program illustrates the 
! use of P3DFFT library for highly scalable parallel 3D FFT. 
!
! This program initializes a 3D array of complex numbers with a 
! 3D sine/cosine wave, then performs inverse FFT transform, and checks that 
! the results are correct. This sample program also demonstrates 
! how to work with complex arrays in wavenumber space, declared as real.
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

      program fft3d

      use p3dfft
      implicit none
      include 'mpif.h'

      integer i,n,nx,ny,nz
      integer m,x,y,z
      integer fstatus
      logical flg_inplace

      real(mytype), dimension(:,:,:),  allocatable :: Fin
      real(mytype), dimension(:,:,:),  allocatable :: Beg
      real(mytype),dimension(:),allocatable:: sinx,siny,sinz
      real(mytype),dimension(:),allocatable:: cosx,cosy,cosz
      real(mytype) pi,twopi,sinyz,diff,cdiff,ccdiff,ans

      integer(8) Ntot
      real(mytype) factor
      real(8) rtime1,rtime2,Nglob,prec
      real(8) gt(12,3),gtcomm(3),tc
      integer ierr,nu,ndim,dims(2),nproc,proc_id
      integer istart(3),iend(3),isize(3)
      integer fstart(3),fend(3),fsize(3)
      integer iproc,jproc
      logical iex

      call MPI_INIT (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)

      print *,'Mytype=',mytype

      twopi=atan(1.0d0)*8.0d0

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

        read (3,*) nx, ny, nz, ndim,n
        write (*,*) "procs=",nproc," nx=",nx, &
                " ny=", ny," nz=", nz,"ndim=",ndim," repeat=", n
        if(mytype .eq. 4) then
           print *,'Single precision version'
        else if(mytype .eq. 8) then
           print *,'Double precision version'
        endif
       endif

      call MPI_Bcast(nx,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(ny,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(nz,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(n,1, MPI_INTEGER,0,mpi_comm_world,ierr)
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

! Set up work structures for P3DFFT
      call p3dfft_setup (dims,nx,ny,nz,.true.)

! Get dimensions for the original array of complex numbers, (X- or Z-pencils
! depending on how the library was compiled)
      call get_dims(istart,iend,isize,2)

! Get dimensions for the C2R-transformed array of real numbers in X-pencils 
! 
      call get_dims(fstart,fend,fsize,1)


! Allocate initial and output arrays

! BEG is real posing as complex
      allocate (BEG(2*istart(1)-1:2*iend(1),istart(2):iend(2),istart(3):iend(3)), stat=ierr)
      if(ierr .ne. 0) then
         print *,'Error ',ierr,' allocating array BEG'
      endif

      allocate (FIN(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)), stat=ierr)
      if(ierr .ne. 0) then
         print *,'Error ',ierr,' allocating array Fin'
      endif

!
! Initialize the array to be transformed (Z-pencils)
!
      allocate (sinx(nx))
      allocate (siny(ny))
      allocate (sinz(nz))
      allocate (cosx(nx))
      allocate (cosy(ny))
      allocate (cosz(nz))

      if(isize(1) .eq. nz) then
! i.e. the library was compiled with --enable-stride1 option and Z is the 
! fastest running index.

         call init_wave1(BEG)
         
      else

         call init_wave2(BEG)

      endif
!
! transform from wavenumber space to physical space
! Repeat n times

      Ntot = fsize(1)*fsize(2)*fsize(3)
      Nglob = nx * ny 
      Nglob = Nglob * nz

      rtime1 = 0.0               

      do  m=1,n
         if(proc_id .eq. 0) then
            print *,'Iteration ',m
         endif
         
         FIN = 0.0d0

! Barrier for correct timing
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         rtime1 = rtime1 - MPI_wtime()

! Call Inverse transform (call a wrapper routine since first argument is
! a real array posing as complex)
         call btran_c2r (BEG,FIN)
         
         rtime1 = rtime1 + MPI_wtime()
         
         if(proc_id .eq. 0) then
            print *,'Result of inverse transform:'
         endif
         call print_all_real(FIN,Ntot,proc_id,Nglob)
                  
      end do

! Free work space
      call p3dfft_clean

! Check results: we expect four non-zero values

      cdiff = 0.0;
      do z=fstart(3),fend(3)
         do y=fstart(2),fend(2)
            do x=fstart(1),fend(1)
               if(x .eq. Nx) then
                  if(y .eq. 3 .and. z .eq. 4) then
                     cdiff = max(abs(FIN(x,y,z)+Nx*Ny*Nz*0.25d0),cdiff)
                  else if(y .eq. 3 .and. z .eq. Nz-2) then
                     cdiff = max(abs(FIN(x,y,z)-Nx*Ny*Nz*0.25d0),cdiff)
                  else if(y .eq. Ny-1 .and. z .eq. 4) then
                     cdiff = max(abs(FIN(x,y,z)-Nx*Ny*Nz*0.25d0),cdiff)
                  else if(y .eq. Ny-1 .and. z .eq. Nz-2) then
                     cdiff = max(abs(FIN(x,y,z)+Nx*Ny*Nz*0.25d0),cdiff)
                  else
                     cdiff = max(abs(FIN(x,y,z)),cdiff)
                  endif
               else
                  cdiff = max(abs(FIN(x,y,z)),cdiff)
               endif

            enddo
         enddo
      enddo

      call MPI_Reduce(cdiff,ccdiff,1,mpireal,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      if(proc_id .eq. 0) then
         if(mytype .eq. 8) then
            prec = 1e-14
         else
            prec = 1e-5
         endif
         if(ccdiff .gt. prec * Nx*Ny*Nz*0.25) then
            print *,'Results are incorrect'
         else
            print *,'Results are correct'
         endif
         print *,'Max diff =',ccdiff
      endif
      
      deallocate(sinx,siny,sinz,cosx,cosy,cosz)
      deallocate(BEG,FIN)

! Gather timing statistics
      call MPI_Reduce(rtime1,rtime2,1,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      if (proc_id.eq.0) write(6,*)'proc_id, cpu time per loop', &
         proc_id,rtime2/dble(n)

      timers = timers / dble(n)

      call MPI_Reduce(timers,gt(1,1),12,mpi_real8,MPI_SUM,0, &
        MPI_COMM_WORLD,ierr)

      call MPI_Reduce(timers,gt(1,2),12,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      call MPI_Reduce(timers,gt(1,3),12,mpi_real8,MPI_MIN,0, &
        MPI_COMM_WORLD,ierr)

      tc = (timers(1)+timers(2)+timers(3)+timers(4))
      call MPI_Reduce(tc,gtcomm(1),1,mpi_real8,MPI_SUM,0, &
        MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tc,gtcomm(2),1,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tc,gtcomm(3),1,mpi_real8,MPI_MIN,0, &
        MPI_COMM_WORLD,ierr)

      gt(1:12,1) = gt(1:12,1) / dble(nproc)
      gtcomm(1) = gtcomm(1) / dble(nproc)

      if(proc_id .eq. 0) then
         do i=1,12
            print *,'timer',i,' (avg/max/min): ',gt(i,:)
         enddo
         print *,'Total comm (avg/max/min): ',gtcomm
      endif


      call MPI_FINALIZE (ierr)

      contains 
!=========================================================
! Initialize complex array passed as real

        subroutine init_wave1(A)

           implicit none

          real(mytype) A((2*istart(1)-1):(2*iend(1)),istart(2):iend(2),istart(3):iend(3))
          integer x,y,z

         do z=istart(1),iend(1)
            sinz(z)=sin((z-1)*3.0d0 * twopi/nz)
         enddo
         do y=istart(2),iend(2)
            siny(y)= sin((y-1)*2.0d0 * twopi/ny)
         enddo
         do x=istart(3),iend(3)
            sinx(x)=sin((x-1)*twopi/nx)
            cosx(x)=cos((x-1)*twopi/nx)
         enddo

         do x=istart(3),iend(3)
            do y=istart(2),iend(2)
               do z=istart(1),iend(1)
                  A(2*z-1,y,x)=cosx(x)*siny(y)*sinz(z)
                  A(2*z,y,x)=sinx(x)*siny(y)*sinz(z)
               enddo
            enddo
         enddo

         end subroutine

         subroutine init_wave2(A)
           
           implicit none

          real(mytype) A((2*istart(1)-1):(2*iend(1)),istart(2):iend(2),istart(3):iend(3))
          integer x,y,z

         do x=istart(1),iend(1)
            sinx(x)=sin((x-1)*twopi/nx)
            cosx(x)=cos((x-1)*twopi/nx)
         enddo
         do y=istart(2),iend(2)
            siny(y)= sin((y-1)*2.0d0 * twopi/ny)
         enddo
         do z=istart(3),iend(3)
            sinz(z)=sin((z-1)*3.0d0 * twopi/nz)
         enddo

         do z=istart(3),iend(3)
            do y=istart(2),iend(2)
               do x=istart(1),iend(1)
                  A(2*x-1,y,z)=cosx(x)*siny(y)*sinz(z)
                  A(2*x,y,z)=sinx(x)*siny(y)*sinz(z)
               enddo
            enddo
         enddo


         end subroutine

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

!=========================================================
! Translate one-dimensional index into three dimensions,
!    print out significantly non-zero values
!

      subroutine print_all_real(Ar,Nar,proc_id,Nglob)

      use p3dfft

      integer x,y,z,proc_id
      integer(8) i,Nar
      real(8) Nglob
      real(mytype), target :: Ar(1,1,*)
      integer Fstart(3),Fend(3),Fsize(3)

      call get_dims(Fstart,Fend,Fsize,1)
      do i=1,Nar
         if(abs(Ar(1,1,i)) .gt. Nglob *1.25e-6) then
            z = (i-1)/(Fsize(1)*Fsize(2))
            y = (i-1 - z * Fsize(1)*Fsize(2))/(Fsize(1))
            x = i -1 - z*Fsize(1)*Fsize(2) - y*Fsize(1)
            print *,'(',x+Fstart(1),y+Fstart(2),z+Fstart(3),') ',Ar(1,1,i)
         endif
      enddo

      return
      end subroutine


      subroutine print_all(Ar,Nar,proc_id,Nglob)

      use p3dfft

      integer x,y,z,proc_id
      integer(8) i,Nar
      complex(mytype) Ar(1,1,*)
      integer Fstart(3),Fend(3),Fsize(3)
      real(8) Nglob

      call get_dims(Fstart,Fend,Fsize,2)

      do i=1,Nar
         if(abs(Ar(1,1,i)) .gt. Nglob *1.25e-4) then
            z = (i-1)/(Fsize(1)*Fsize(2))
            y = (i-1 - z * Fsize(1)*Fsize(2))/Fsize(1)
            x = i-1-z*Fsize(1)*Fsize(2) - y*Fsize(1)
            print *,'(',x+Fstart(1),y+Fstart(2),z+Fstart(3),') ',Ar(1,1,i)
         endif
      enddo

      return
      end subroutine

      end
