! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2014 Dmitry Pekurovsky
!    Copyright (C) 2006-2014 University of California
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

! This is a C wrapper routine
!========================================================
      subroutine p3dfft_ftran_r2c_many_w (XgYZ,dim_in,XYZg,dim_out,nv,op) BIND(C,NAME='p3dfft_ftran_r2c_many')
!========================================================
      use, intrinsic :: iso_c_binding
      real(mytype), TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(mytype), TARGET :: XYZg(nzc,iistart:iiend,jjstart:jjend)
#else
      complex(mytype), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nzc)
#endif
      integer dim_in,dim_out,nv
      character, dimension(*), target :: op
      character(4), pointer :: lcl_op
      call c_f_pointer(c_loc(op), lcl_op)

      call p3dfft_ftran_r2c_many (XgYZ,dim_in,XYZg,dim_out,nv,lcl_op) 

      end subroutine

! Forward R2C transform of multiple variables (nv)
!========================================================
      subroutine p3dfft_ftran_r2c_many (XgYZ,dim_in,XYZg,dim_out,nv,op) 
!========================================================

      use fft_spec
      implicit none

      integer dim_in,dim_out

      real(mytype), TARGET :: XgYZ(dim_in,nv)
      complex(mytype), TARGET :: XYZg(dim_out,nv)

      integer x,y,z,i,nx,ny,nz,ierr,dnz,nv,j,err,n1,n2
      integer(i8) Nl
      character(len=3) op
      if(.not. mpi_set) then
         print *,'P3DFFT error: call setup before other routines'
         return
      endif
      
      if(dim_in .lt. nx_fft*jisize*kjsize) then
         print *,taskid,': ftran error: input array dimensions are too low: ',dim_in,' while expecting ',nx_fft*jisize*kjsize
      endif	 

      if(dim_out .lt. nzc*jjsize*iisize) then
         print *,taskid,': ftran error: output array dimensions are too low: ',dim_out,' while expecting ',nzc*jjsize*iisize
      endif	 

!     preallocate memory for FFT-Transforms

      if(nv .gt. nv_preset) then
        nv_preset = nv
        deallocate(buf1,buf2,buf)

        allocate (buf(nxhp*jisize*(kjsize+padd)*nv), stat=err)
        if (err /= 0) then
          print *, 'Error ', err, ' allocating array buf'
        end if
!     initialize buf to avoid "floating point invalid" errors in debug mode
       buf = 0.d0

#ifdef USE_EVEN
        n1 = nv * IfCntMax * iproc /(mytype*2)
        n2 = nv * KfCntMax * jproc / (mytype*2)
        n1 = max(n1,n2)
        allocate(buf1(n1))
        allocate(buf2(n1))
#else
	if(taskid .eq. 0) then
	  print *,'nm=',nm
	endif
        allocate(buf1(nm*nv))
        allocate(buf2(nm*nv))        	
#endif
      endif

      nx = nx_fft
      ny = ny_fft
      nz = nz_fft

! For FFT libraries that require explicit allocation of work space,
! such as ESSL, initialize here

#ifdef DEBUG
	print *,taskid,': Enter ftran',nv,nv_preset
#endif

! FFT transform (R2C) in X for all z and y

      if(jisize * kjsize .gt. 0) then
         call init_f_r2c(XgYZ,nx,buf,nxhp,nx,jisize*kjsize)

         timers(5) = timers(5) - MPI_Wtime()
	 call f_r2c_many(XgYZ,nx,buf,nxhp,nx,jisize*kjsize,dim_in,nv)
         timers(5) = timers(5) + MPI_Wtime()

      endif

! Exchange data in rows 

      if(iproc .gt. 1) then 

#ifdef DEBUG
	print *,taskid,': Calling fcomm1'
#endif
         call fcomm1_many(buf,buf,nv,timers(1),timers(6))
         

#ifdef STRIDE1
      else
         call reorder_f1_many(buf,buf,buf1,nv)
#endif
      endif

! FFT transform (C2C) in Y for all x and z, one Z plane at a time


#ifdef DEBUG
	print *,taskid,': Transforming in Y'
#endif

      if(iisize * kjsize .gt. 0) then
#ifdef STRIDE1
         call init_f_c(buf,1,ny,buf,1,ny,ny,iisize*kjsize)

         timers(7) = timers(7) - MPI_Wtime()
	 call f_c1_many(buf,1,ny,ny,iisize*kjsize,iisize*kjsize*ny,nv)
	 timers(7) = timers(7) + MPI_Wtime()

#else
         call init_f_c(buf,iisize,1,buf,iisize,1,ny,iisize)


         timers(7) = timers(7) - MPI_Wtime()

         do z=1,kjsize*nv
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
         call fcomm2_trans_many(buf,XYZg,buf,dim_out,nv,op,timers(2),timers(8))
#else

! FFT Transform (C2C) in Z for all x and y

         dnz = nz - nzc
	 if(dnz .gt. 0) then

! Transpose y-z

         call fcomm2_many(buf,buf,iisize*jjsize*nz,nv,timers(2),timers(8))

! In forward transform we can safely use output array as one of the buffers
! This speeds up FFTW since it is non-stride-1 transform and it is 
! faster than done in-place

#ifdef DEBUG
	print *,taskid,': Transforming in Z'
#endif
         if(iisize * jjsize .gt. 0) then
	    call ztran_f_same_many(buf,iisize*jjsize,1,nz,iisize*jjsize,iisize*jjsize*nz,nv,op)
            call seg_copy_z_f_many(buf,XYZg,1,iisize,1,jjsize,1,nzhc,0,iisize,jjsize,nz,dim_out,nv)
            call seg_copy_z_f_many(buf,XYZg,1,iisize,1,jjsize,nzhc+1,nzc,dnz,iisize,jjsize,nz,dim_out,nv)
	endif	
      else
        call fcomm2_many(buf,XYZg,dim_out,nv,timers(2),timers(8))

! In forward transform we can safely use output array as one of the buffers
! This speeds up FFTW since it is non-stride-1 transform and it is 
! faster than done in-place

#ifdef DEBUG
        print *,taskid,': Transforming in Z'
#endif
         if(iisize * jjsize .gt. 0) then
	    call ztran_f_same_many(XYZg,iisize*jjsize,1,nz,iisize*jjsize,dim_out,nv,op)
        endif
     endif	

#endif

      else

         timers(8) = timers(8) - MPI_Wtime()	

#ifdef STRIDE1
         call reorder_trans_f2_many(buf,XYZg,buf1,dim_out,nv,op)
#else
         Nl = iisize*jjsize*nz
         dnz = nz - nzc
	 if(dnz .gt. 0) then

	    call ztran_f_same_many(buf,iisize*jjsize,1,nz,iisize*jjsize,iisize*jjsize*nz,nv,op)
            call seg_copy_z_f_many(buf,XYZg,1,iisize,1,jjsize,1,nzhc,0,iisize,jjsize,nz,dim_out,nv)
            call seg_copy_z_f_many(buf,XYZg,1,iisize,1,jjsize,nzhc+1,nzc,dnz,iisize,jjsize,nz,dim_out,nv)
	else

         call ar_copy_many(buf,iisize*jjsize*nz,XYZg,dim_out,Nl,nv)
	 call ztran_f_same_many(XYZg,iisize*jjsize,1,nz,iisize*jjsize,dim_out,nv,op)
        endif
#endif

        timers(8) = timers(8) + MPI_Wtime()

      endif

!      deallocate(buf)

      return
      end subroutine

! This is a C wrapper routine
!========================================================
      subroutine p3dfft_ftran_cheby_many_w (XgYZ,dim_in,XYZg,dim_out,nv,Lz) BIND(C,NAME='p3dfft_cheby_many')
!========================================================

      real(mytype), TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(mytype), TARGET :: XYZg(nzc,iistart:iiend,jjstart:jjend)
#else
      complex(mytype), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nzc)
#endif
      integer dim_in,dim_out,nv
      real(mytype) Lz

      call p3dfft_cheby_many (XgYZ,dim_in,XYZg,dim_out,nv,Lz) 

      end subroutine

! Chebyshev transform (2D R2C forward FFT + Chebyshev) for multiple variables 
!========================================================
	subroutine p3dfft_cheby_many(in,dim_in,out,dim_out,nv,Lz) 
!========================================================

	integer dim_in,dim_out,nv,j
	real(mytype) Lz
      	real(mytype), dimension(dim_in,nv), target	::	in
        complex(mytype), dimension (dim_out,nv), target :: out				

    	call p3dfft_ftran_r2c_many(in,dim_in,out,dim_out,nv,'ffc')

	do j=1,nv
	   call p3dfft_cheby(in(1,j),out(1,j),Lz)
        enddo

	return
	end subroutine


! This is a C wrapper routine
!========================================================
	subroutine p3dfft_cheby_w(in,out,Lz) BIND(C,NAME='p3dfft_cheby')
!========================================================

      	real(mytype), dimension(nx_fft,     &
                                jisize,     &
                                kjsize), target	::	in

#ifdef STRIDE1
      	complex(mytype), dimension(nzc, &
                                jjsize,&
                                iisize), target	::	out
#else
      	complex(mytype), dimension(iisize, &
                                jjsize,    &
                                nzc), target	::	out
#endif
	real(mytype) Lz

	call p3dfft_cheby(in,out,Lz) 

	end subroutine

! Chebyshev transform (2D R2C forward FFT + Chebyshev) for a single variable 
!==============================================================
	subroutine p3dfft_cheby(in,out,Lz) 
!========================================================

!     !	function args
      	real(mytype), dimension(nx_fft,     &
                                jisize,     &
                                kjsize), target	::	in

#ifdef STRIDE1
      	complex(mytype), dimension(nzc, &
                                jjsize,&
                                iisize), target	::	out
      	complex(mytype) Old, New, Tmp
#else
      	complex(mytype), dimension(iisize, &
                                jjsize,    &
                                nzc), target	::	out
      	complex(mytype),dimension(:,:),pointer :: ptrOld, ptrNew, ptrTmp
#endif
      	real(mytype) :: Lfactor,Lz
	integer k,nz,i,j

	nz = nzc

    	call p3dfft_ftran_r2c(in,out,'ffc')
     	out = out *(1.d0/(dble(nx_fft*ny_fft)*(nzc-1)))

! less tmp-memory version (but difficult to read)

     	Lfactor = 4.d0/Lz

#ifdef STRIDE1
! first and last cheby-coeff needs to gets multiplied by factor 0.5
! because of relation between cheby and discrete cosinus transforms

  
	do i=1,iisize
	   do j=1,jjsize
    	     Old = out(nzc-1,j,i)
    	     out(nzc-1,j,i) = Lfactor *dble(nzc-1) *out(nzc,j,i) * 0.5d0
    	     out(nzc,j,i) = cmplx(0.d0,0.d0)
    	     do k = nzc-2, 1, -1
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

    	ptrOld(:,:) = out(:,:,nzc-1)
    	out(:,:,nzc-1) = Lfactor *dble(nzc-1) *out(:,:,nzc) *0.5d0
    	out(:,:,nzc  ) = cmplx(0.d0,0.d0)
    	do k = nzc-2, 1, -1
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


! This is a C wrapper routine
!========================================================
      subroutine p3dfft_ftran_r2c_w (XgYZ,XYZg,op) BIND(C,NAME='p3dfft_ftran_r2c')
!========================================================
      use, intrinsic :: iso_c_binding
      real(mytype), TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(mytype), TARGET :: XYZg(nzc,iistart:iiend,jjstart:jjend)
#else
      complex(mytype), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nzc)
#endif
      integer dim_in,dim_out,nv
      character, dimension(*), target :: op
      character(4), pointer :: lcl_op
      call c_f_pointer(c_loc(op), lcl_op)

      call p3dfft_ftran_r2c (XgYZ,XYZg,lcl_op) 

      end subroutine

! Forward R2C transform of 1 variable
!========================================================
      subroutine p3dfft_ftran_r2c (XgYZ,XYZg,op) 
!========================================================

      use fft_spec
      implicit none

      real(mytype), TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(mytype), TARGET :: XYZg(nzc,iistart:iiend,jjstart:jjend)
#else
      complex(mytype), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nzc)
#endif

      integer x,y,z,i,nx,ny,nz,ierr,dnz
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
         call reorder_f1(buf,buf,buf1)
#endif
      endif

! FFT transform (C2C) in Y for all x and z, one Z plane at a time


#ifdef DEBUG
	print *,taskid,': Transforming in Y'
#endif

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

#ifdef DEBUG
	print *,taskid,': Calling fcomm2'
#endif


! Exchange data in columns
      if(jproc .gt. 1) then

#ifdef STRIDE1
! For stride1 option combine second transpose with transform in Z
         call init_f_c(buf,1,nz, XYZg,1,nz,nz,jjsize,op)
         call fcomm2_trans(buf,XYZg,buf,op,timers(2),timers(8))
#else

! FFT Transform (C2C) in Z for all x and y

         dnz = nz - nzc
	 if(dnz .gt. 0) then

! Transpose y-z

         call fcomm2(buf,buf,timers(2),timers(8))

! In forward transform we can safely use output array as one of the buffers
! This speeds up FFTW since it is non-stride-1 transform and it is 
! faster than done in-place

#ifdef DEBUG
	print *,taskid,': Transforming in Z'
#endif
         if(iisize * jjsize .gt. 0) then
	    if(op(3:3) == 't' .or. op(3:3) == 'f') then
               call init_f_c(buf,iisize*jjsize, 1, buf,iisize*jjsize, 1,nz,iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_f_c2_same(buf,iisize*jjsize, 1,buf,iisize*jjsize, 1,nz,iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 'c') then
               call init_ctrans_r2(buf,2*iisize*jjsize, 1, buf,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_ctrans_r2_same(buf,2*iisize*jjsize, 1,buf,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 's') then
               call init_strans_r2(buf,2*iisize*jjsize, 1, buf,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_strans_r2_same(buf,2*iisize*jjsize, 1,buf,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()
	    else if(op(3:3) .ne. 'n' .and. op(3:3) .ne. '0') then
		print *,'Unknown transform type: ',op(3:3)
		call MPI_Abort(MPI_COMM_WORLD,ierr)
            endif

	    call seg_copy_z(buf,XYZg,1,iisize,1,jjsize,1,nzhc,0,iisize,jjsize,nz)
	    call seg_copy_z(buf,XYZg,1,iisize,1,jjsize,nzhc+1,nzc,dnz,iisize,jjsize,nz)

	endif	
      else
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
              call exec_f_c2_same(XYZg,iisize*jjsize, 1,XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

            else if(op(3:3) == 'c') then
               call init_ctrans_r2(XYZg,2*iisize*jjsize, 1, XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_ctrans_r2_same(XYZg,2*iisize*jjsize, 1,XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

            else if(op(3:3) == 's') then
               call init_strans_r2(XYZg,2*iisize*jjsize, 1, XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_strans_r2_same(XYZg,2*iisize*jjsize, 1,XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()
            else if(op(3:3) .ne. 'n' .and. op(3:3) .ne. '0') then
                print *,'Unknown transform type: ',op(3:3)
                call MPI_Abort(MPI_COMM_WORLD,ierr)
            endif

        endif
     endif	


#endif

      else

         timers(8) = timers(8) - MPI_Wtime()	

#ifdef STRIDE1
         call reorder_trans_f2(buf,XYZg,buf1,op)
#else
         Nl = iisize*jjsize*nz
         dnz = nz - nzc
	 if(dnz .gt. 0) then

	    if(op(3:3) == 't' .or. op(3:3) == 'f') then
               call init_f_c(buf,iisize*jjsize, 1,buf,iisize*jjsize, 1,nz,iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_f_c2_same(buf,iisize*jjsize, 1,buf,iisize*jjsize, 1,nz,iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 'c') then
               call init_ctrans_r2(buf,2*iisize*jjsize, 1,buf,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_ctrans_r2_same(buf,2*iisize*jjsize, 1,buf,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 's') then
               call init_strans_r2(buf,2*iisize*jjsize, 1,buf,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_strans_r2_same(buf,2*iisize*jjsize, 1,buf,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()
	    else if(op(3:3) /= 'n' .and. op(3:3) /= '0') then
		print *,'Unknown transform type: ',op(3:3)
		call MPI_Abort(MPI_COMM_WORLD,ierr)
            endif

	   call seg_copy_z(buf,XYZg,1,iisize,1,jjsize,1,nzhc,0,iisize,jjsize,nz)
	   call seg_copy_z(buf,XYZg,1,iisize,1,jjsize,nzhc+1,nzc,dnz,iisize,jjsize,nz)

	else

         call ar_copy(buf,XYZg,Nl)
            if(op(3:3) == 't' .or. op(3:3) == 'f') then
               call init_f_c(XYZg,iisize*jjsize, 1, XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_f_c2_same(XYZg,iisize*jjsize, 1,XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

            else if(op(3:3) == 'c') then
               call init_ctrans_r2(XYZg,2*iisize*jjsize, 1, XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_ctrans_r2_same(XYZg,2*iisize*jjsize, 1,XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

            else if(op(3:3) == 's') then
               call init_strans_r2(XYZg,2*iisize*jjsize, 1, XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
         
              timers(8) = timers(8) - MPI_Wtime()
              call exec_strans_r2_same(XYZg,2*iisize*jjsize, 1,XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()
            else if(op(3:3) /= 'n' .and. op(3:3) /= '0') then
                print *,'Unknown transform type: ',op(3:3)
                call MPI_Abort(MPI_COMM_WORLD,ierr)
            endif


        endif
#endif

        timers(8) = timers(8) + MPI_Wtime()

      endif

      return
      end subroutine

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

!	 call f_r2c(XgYZ,nx,buf,nxhp,nx,jisize*kjsize,nv)

subroutine f_r2c_many(source,str1,dest,str2,n,m,dim,nv)

  integer str1,str2,n,m,nv,j,dim
  real(mytype) source(dim,nv)
  complex(mytype) dest(n/2+1,m,nv)

  do j=1,nv
    call exec_f_r2c(source(1,j),str1,dest(1,1,j),str2,n,m)
  enddo

  return
  end subroutine

         subroutine f_c1_many(A,str1,str2,n,m,dim,nv)

	   integer n,m,nv,j,str1,str2,dim
	   complex(mytype) A(dim,nv)

	 do j=1,nv
           call exec_f_c1(A(1,j),str1,str2,A(1,j),str1,str2,n,m)
         enddo

	 return
	 end subroutine

         subroutine ztran_f_same_many(A,str1,str2,n,m,dim,nv,op)
	
	   integer str1,str2,n,m,nv,j,ierr,dim
	   complex(mytype) A(dim,nv)
	   character(len=3) op

	    if(op(3:3) == 't' .or. op(3:3) == 'f') then
               call init_f_c(A,str1,str2,A,str1,str2,n,m)
         
              timers(8) = timers(8) - MPI_Wtime()
	      do j=1,nv
                 call exec_f_c2_same(A(1,j),str1,str2,A(1,j),str1,str2,n,m)
              enddo
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 'c') then
               call init_ctrans_r2(A,str1,str2,A,str1,str2,n,m)
         
              timers(8) = timers(8) - MPI_Wtime()
	      do j=1,nv
                 call exec_ctrans_r2_same(A(1,j),2*str1,str2,A(1,j),2*str1,str2,n,2*m)
	      enddo
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 's') then
               call init_strans_r2(A,str1,str2,A,str1,str2,n,m)
         
              timers(8) = timers(8) - MPI_Wtime()
	      do j=1,nv
                 call exec_strans_r2_same(A(1,j),2*str1,str2,A(1,j),2*str1,str2,n,2*m)
              enddo
              timers(8) = timers(8) + MPI_Wtime()
	    else if(op(3:3) .ne. 'n' .and. op(3:3) .ne. '0') then
		print *,'Unknown transform type: ',op(3:3)
		call MPI_Abort(MPI_COMM_WORLD,ierr)
            endif

	    return
	    end subroutine


