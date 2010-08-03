! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2010 Dmitry Pekurovsky
!    Copyright (C) 2006-2010 University of California
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!
!----------------------------------------------------------------------------

      subroutine p3dfft_ftran_r2c (XgYZ,XYZg)
!========================================================

      use fft_spec
      implicit none

      real(mytype), TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(mytype), TARGET :: XYZg(nz_fft,iistart:iiend,jjstart:jjend)
#else
      complex(mytype), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nz_fft)
#endif
      integer x,y,z,i,err,nx,ny,nz,ithr,OMP_GET_THREAD_NUM,ierr,id
      integer(i8) Nl
      
      if(.not. mpi_set) then
         print *,'P3DFFT error: call setup before other routines'
         return
      endif

      nx = nx_fft
      ny = ny_fft
      nz = nz_fft

! For FFT libraries that require explicit allocation of work space,
! such as ESSL, initialize here

!	print *,taskid,': Enter ftran'

      call init_work(nx,ny,nz)

! FFT transform (R2C) in X for all z and y


      if(jisize * kjsize .gt. 0) then
         call init_f_r2c(XgYZ,nx,buf,nxhp,nx,jisize*kjsize)

         timers(5) = timers(5) - MPI_Wtime()

         call exec_f_r2c(XgYZ,nx,buf,nxhp,nx,jisize*kjsize)

         timers(5) = timers(5) + MPI_Wtime()

      endif


! Exchange data in rows 

      if(iproc .gt. 1) then 

!	print *,taskid,': Calling fcomm1'
         call fcomm1(buf,buf,timers(1),timers(6))
         

#ifdef STRIDE1
      else
         call reorder_f1(buf,buf)
#endif
      endif

! FFT transform (C2C) in Y for all x and z, one Z plane at a time


!	print *,taskid,': Transforming in Y'

      if(iisize * kjsize .gt. 0) then
#ifdef STRIDE1
         call init_f_c(buf,1,ny,buf,1,ny,ny,iisize*kjsize)

         timers(7) = timers(7) - MPI_Wtime()

         call exec_f_c1(buf,1,ny,buf,1,ny,ny,iisize*kjsize)
         timers(7) = timers(7) + MPI_Wtime()


#else
         call init_f_c(buf,iisize,1,buf,iisize,1,ny,iisize)


         timers(7) = timers(7) - MPI_Wtime()

         do z=1,kjsize
            call ftran_y_zplane(buf,z-1,iisize,kjsize,iisize,1, buf,z-1,iisize,kjsize,iisize,1,ny,iisize)
         enddo
         timers(7) = timers(7) + MPI_Wtime()

#endif
      endif

!	print *,taskid,': Calling fcomm2'


! Exchange data in columns
      if(jproc .gt. 1) then

#ifdef STRIDE1
! For stride1 option combine second transpose with transform in Z
         call init_f_c(buf,1,nz, buf,1,nz,nz,jjsize)
         call fcomm2_trans(buf,XYZg,timers(2),timers(8))
#else

! FFT Transform (C2C) in Z for all x and y

! Transpose y-z

         call fcomm2(buf,XYZg,timers(2),timers(8))

! In forward transform we can safely use output array as one of the buffers
! This speeds up FFTW since it is non-stride-1 transform and it is 
! faster than done in-place

!	print *,taskid,': Transforming in Z'

         if(iisize * jjsize .gt. 0) then
            call init_f_c(XYZg,iisize*jjsize, 1, XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
         
            timers(8) = timers(8) - MPI_Wtime()
            call exec_f_c2(XYZg,iisize*jjsize, 1,XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
            timers(8) = timers(8) + MPI_Wtime()
         endif
#endif

      else

#ifdef STRIDE1
         call reorder_trans_f2(buf,XYZg)
#else
         Nl = iisize*jjsize*nz
         call ar_copy(buf,XYZg,Nl)
         call init_f_c(XYZg,iisize*jjsize, 1, XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
         call exec_f_c2(XYZg,iisize*jjsize, 1, XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
#endif
      endif


      call free_work

      return
      end subroutine
