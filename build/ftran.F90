! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2010 Dmitry Pekurovsky
!    Copyright (C) 2006-2010 University of California
!    Copyright (C) 2010-2011 Jens Henrik Goebbert
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

      subroutine p3dfft_ftran_r2c (XgYZ,XYZg,op)
!========================================================

      use fft_spec
      implicit none

      real(mytype), TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(mytype), TARGET :: XYZg(nz_fft,iistart:iiend,jjstart:jjend)
#else
      complex(mytype), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nz_fft)
#endif
      integer x,y,z,i,nx,ny,nz,ierr
      integer(i8) Nl
      character(len=3) op
      
      if(.not. mpi_set) then
         print *,'P3DFFT error: call setup before other routines'
         return
      endif

      nx = nx_fft
      ny = ny_fft
      nz = nz_fft

! For FFT libraries that require explicit allocation of work space,
! such as ESSL, initialize here

#ifdef DEBUG
	print *,taskid,': Enter ftran'
#endif
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

#ifdef DEBUG
	print *,taskid,': Calling fcomm1'
#endif
         call fcomm1(buf,buf,timers(1),timers(6))
         

#ifdef STRIDE1
      else
         call reorder_f1(buf,buf)
#endif
      endif

! FFT transform (C2C) in Y for all x and z, one Z plane at a time


#ifdef DEBUG
	print *,taskid,': Transforming in Y'
#endif

      if(iisize * kjsize .gt. 0) then
#ifdef STRIDE1
         call init_f_c(buf,1,ny,XYZg,1,ny,ny,iisize*kjsize)

         timers(7) = timers(7) - MPI_Wtime()

         call exec_f_c1(buf,1,ny,XYZg,1,ny,ny,iisize*kjsize)
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

#ifdef DEBUG
	print *,taskid,': Calling fcomm2'
#endif


! Exchange data in columns
      if(jproc .gt. 1) then

#ifdef STRIDE1
! For stride1 option combine second transpose with transform in Z
         call init_f_c(buf,1,nz, XYZg,1,nz,nz,jjsize,op)
         call fcomm2_trans(XYZg,XYZg,buf,op,timers(2),timers(8))
#else

! FFT Transform (C2C) in Z for all x and y

! Transpose y-z

         call fcomm2(buf,XYZg,timers(2),timers(8))

! In forward transform we can safely use output array as one of the buffers
! This speeds up FFTW since it is non-stride-1 transform and it is 
! faster than done in-place

#ifdef DEBUG
	print *,taskid,': Transforming in Z'
#endif
         if(iisize * jjsize .gt. 0) then
	    if(op(3:3) == 't' .or. op(3:3) == 'f') then
               call init_f_c(XYZg,iisize*jjsize, 1, XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_f_c2(XYZg,iisize*jjsize, 1,XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 'c') then
               call init_ctrans_r2(XYZg,2*iisize*jjsize, 1, XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_ctrans_r2(XYZg,2*iisize*jjsize, 1,XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 's') then
               call init_strans_r2(XYZg,2*iisize*jjsize, 1, XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_strans_r2(XYZg,2*iisize*jjsize, 1,XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()
	    else if(op(3:3) /= 'n' .and. op(3:3) /= '0') then
		print *,'Unknown transform type: ',op(3:3)
		call MPI_Abort(MPI_COMM_WORLD,ierr)
            endif
	endif	
#endif

      else

         timers(8) = timers(8) - MPI_Wtime()	

#ifdef STRIDE1
         call reorder_trans_f2(XYZg,XYZg,buf,op)
#else
         Nl = iisize*jjsize*nz
         call ar_copy(buf,XYZg,Nl)
	    if(op(3:3) == 't' .or. op(3:3) == 'f') then
               call init_f_c(XYZg,iisize*jjsize, 1, XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_f_c2(XYZg,iisize*jjsize, 1,XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 'c') then
               call init_ctrans_r2(XYZg,2*iisize*jjsize, 1, XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_ctrans_r2(XYZg,2*iisize*jjsize, 1,XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 's') then
               call init_strans_r2(XYZg,2*iisize*jjsize, 1, XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_strans_r2(XYZg,2*iisize*jjsize, 1,XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()
	    else if(op(3:3) /= 'n' .and. op(3:3) /= '0') then
		print *,'Unknown transform type: ',op(3:3)
		call MPI_Abort(MPI_COMM_WORLD,ierr)
            endif
#endif

         timers(8) = timers(8) + MPI_Wtime()

      endif


      call free_work

      return
      end subroutine

!---------------------------------------------------------------
	subroutine p3dfft_cheby(in,out,Lz)

!     !	function args
      	real(mytype), dimension(nx_fft,     &
                                jisize,     &
                                kjsize), target	::	in

#ifdef STRIDE1
      	complex(mytype), dimension(nz_fft, &
                                jjsize,&
                                iisize), target	::	out
      	complex(mytype) Old, New, Tmp
#else
      	complex(mytype), dimension(iisize, &
                                jjsize,    &
                                nz_fft), target	::	out
      	complex(mytype),dimension(:,:),pointer :: ptrOld, ptrNew, ptrTmp
#endif
      	real(mytype) :: Lfactor,Lz
	integer k,nz,i,j

	nz = nz_fft


    	call p3dfft_ftran_r2c(in,out,'ffc')
     	out = out *(1.d0/(dble(nx_fft*ny_fft)*(nz_fft-1)))

! less tmp-memory version (but difficult to read)

     	Lfactor = 4.d0/Lz

#ifdef STRIDE1
! first and last cheby-coeff needs to gets multiplied by factor 0.5
! because of relation between cheby and discrete cosinus transforms

	do i=1,iisize
	   do j=1,jjsize
    	     Old = out(nz-1,j,i)
    	     out(nz-1,j,i) = Lfactor *dble(nz-1) *out(nz,j,i) * 0.5d0
    	     out(nz,j,i) = cmplx(0.d0,0.d0)
    	     do k = nz-2, 1, -1
    		New = out(k,j,i)
    		out(k,j,i) = Lfactor *dble(k) *Old +out(k+2,j,i)
    		Tmp = New
    		New = Old
    		Old = Tmp
      	     enddo
    	     out(1,j,i) = out(1,j,i) *0.5d0
	   enddo
	enddo
#else
   	allocate(ptrOld(iisize,jjsize))
      	allocate(ptrNew(iisize,jjsize))

    	ptrOld(:,:) = out(:,:,nz-1)
    	out(:,:,nz-1) = Lfactor *dble(nz-1) *out(:,:,nz) *0.5d0
    	out(:,:,nz  ) = cmplx(0.d0,0.d0)
    	do k = nz-2, 1, -1
    		ptrNew(:,:) = out(:,:,k)
    		out(:,:,k) = Lfactor *dble(k) *ptrOld(:,:) +out(:,:,k+2)
    		ptrTmp => ptrNew
    		ptrNew => ptrOld
    		ptrOld => ptrTmp
    	enddo
    	out(:,:,1) = out(:,:,1) *0.5d0
      	deallocate(ptrOld)
      	deallocate(ptrNew)
#endif

!     		! easy to read version (but more tmp-memory needed)
!     		Lfactor = 4.d0/Lz
!     		ctest10 = out
!    		out(:,:,n3  ) = cmplx(0.d0,0.d0)					!ok
!    		out(:,:,n3-1) = Lfactor *dble(n3-1) *ctest10(:,:,n3)	!ok
!    		do k = n3-2, 1, -1
!    			out(:,:,k) = Lfactor *dble(k) *ctest10(:,:,k+1) +out(:,:,k+2)
!    		enddo
!    			out(:,:,1) = out(:,:,1) *0.5d0

	return
	end subroutine p3dfft_cheby

! --------------------------------------
!
!  p3dfft_ftran_r2c_1d(..)
!
! --------------------------------------
subroutine p3dfft_ftran_r2c_1d (rXgYZ, cXgYZ)
  use fft_spec
  implicit none

  real (mytype), target :: rXgYZ (NX_fft, jistart:jiend, kjstart:kjend)
  real (mytype), target :: cXgYZ (NX_fft+2, jistart:jiend, kjstart:kjend)

!      complex(mytype), allocatable :: XYgZ(:,:,:)
  integer x, y, z, i, err, nx, ny, nz

  if ( .not. mpi_set) then
    print *, 'P3DFFT error: call setup before other routines'
    return
  end if

  nx = NX_fft
  ny = NY_fft
  nz = NZ_fft

!
! FFT transform (R2C) in X for all z and y
!
  if (jisize*kjsize > 0) then
    call init_f_r2c (rXgYZ, nx, cXgYZ, nxhp, nx, jisize*kjsize)
    call exec_f_r2c (rXgYZ, nx, cXgYZ, nxhp, nx, jisize*kjsize)
  end if

end subroutine p3dfft_ftran_r2c_1d
