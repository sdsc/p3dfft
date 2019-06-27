! This file is part of P3DFFT library

! Title: P3DFFT library

! Authors: Dmitry Pekurovsky

! Copyright (c) 2006-2019 

! The Regents of the University of California.

! All Rights Reserved.                        

 

!    Permission to use, copy, modify and  distribute  any part

!    of this software for  educational,  research  and  non-profit

!    purposes, by individuals or non-profit organizations,

!    without fee,  and  without a written  agreement is

!    hereby granted,  provided  that the  above  copyright notice,

!    this paragraph  and the following  three  paragraphs appear in

!    all copies.       

 

!    For-profit organizations desiring to use this software and others

!    wishing to incorporate this  software into commercial

!    products or use it for  commercial  purposes should contact the:    

!          Office of Innovation & Commercialization 

!          University of California San Diego

!          9500 Gilman Drive,  La Jolla,  California, 92093-0910        

!          Phone: (858) 534-5815

!          E-mail: innovation@ucsd.edu

 

!    IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE

!    TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR    

!    CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT

!    OF THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY OF

!    CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH

!    DAMAGE.

 

!    THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND

!    THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE        

!    MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. 

!    THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND    

!    EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR

!    IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES

!    OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR

!    THAT THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT,        

!    TRADEMARK OR OTHER RIGHTS.


! This is a C wrapper routine
!========================================================
      subroutine p3dfft_ftran_r2c_many_w (XgYZ,dim_in,XYZg,dim_out,nv,op) BIND(C,NAME='p3dfft_ftran_r2c_many')
!========================================================
      use, intrinsic :: iso_c_binding
      real(p3dfft_type), TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(p3dfft_type), TARGET :: XYZg(nzc,jjstart:jjend,iistart:iiend)
#else
      complex(p3dfft_type), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nzc)
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

      real(p3dfft_type), TARGET :: XgYZ(dim_in,nv)
      complex(p3dfft_type), TARGET :: XYZg(dim_out,nv)

      integer x,y,z,i,nx,ny,nz,ierr,dnz,nv,j,err,n1,n2,dny
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

        allocate (buf(nxhp*jisize*(kjsize+padi)*nv), stat=err)
        if (err /= 0) then
          print *, 'Error ', err, ' allocating array buf'
        end if
!     initialize buf to avoid "floating point invalid" errors in debug mode
       buf = 0.d0

#ifdef USE_EVEN
        n1 = nv * IfCntMax * iproc /(p3dfft_type*2)
        n2 = nv * KfCntMax * jproc / (p3dfft_type*2)
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
         call init_f_r2c(XgYZ,nx,buf2,nxhp,nx,jisize*kjsize)

         timers(5) = timers(5) - MPI_Wtime()
	 call f_r2c_many(XgYZ,nx,buf2,nxhp,nx,jisize*kjsize,dim_in,nv)
         timers(5) = timers(5) + MPI_Wtime()

      endif

! Exchange data in rows

      if(iproc .gt. 1) then

#ifdef DEBUG
	print *,taskid,': Calling fcomm1'
#endif
         call fcomm1_many(buf2,buf,nv,timers(1),timers(6))

      else

#ifdef STRIDE1
         call reorder_f1_many(buf2,buf,buf1,nv)
#else
 	 call seg_copy_x_f_many(buf2,buf,1,nxhpc,0,nxhp,nxhpc,jisize,kjsize,nxhpc*jisize*kjsize,nv)
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
         call init_f_c(buf,1,nz, XYZg,1,nz,nz,jjsize)
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

	    dny = ny - nyc
	    call seg_copy_y_f_many(buf,buf1,1,nyhc,0,iisize,ny,nyc,nz,iisize*nyc*nz,nv)
	    call seg_copy_y_f_many(buf,buf1,nycph+1,nyc,dny,iisize,ny,nyc,nz,iisize*nyc*nz,nv)

	    call ztran_f_same_many(buf1,iisize*jjsize,1,nz,iisize*jjsize,iisize*jjsize*nz,nv,op)
            call seg_copy_z_f_many(buf1,XYZg,1,iisize,1,jjsize,1,nzhc,0,iisize,jjsize,nz,dim_out,nv)
            call seg_copy_z_f_many(buf1,XYZg,1,iisize,1,jjsize,nzhc+1,nzc,dnz,iisize,jjsize,nz,dim_out,nv)
	else

	    dny = ny - nyc
	    call seg_copy_y_f_many(buf,XYZg,1,nycph,0,iisize,ny,nyc,nz,iisize*nyc*nz,nv)
	    call seg_copy_y_f_many(buf,XYZg,nycph+1,nyc,dny,iisize,ny,nyc,nz,iisize*nyc*nz,nv)

	    call ztran_f_same_many(XYZg,iisize*jjsize,1,nz,iisize*jjsize,dim_out,nv,op)
        endif
#endif

        timers(8) = timers(8) + MPI_Wtime()

      endif

!      deallocate(buf)

!      call mpi_barrier(mpi_comm_world,ierr)

     return
      end subroutine

! This is a C wrapper routine
!========================================================
      subroutine p3dfft_ftran_cheby_many_w (XgYZ,dim_in,XYZg,dim_out,nv,Lz) BIND(C,NAME='p3dfft_cheby_many')
!========================================================

      real(p3dfft_type), TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(p3dfft_type), TARGET :: XYZg(nzc,jjstart:jjend,iistart:iiend)
#else
      complex(p3dfft_type), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nzc)
#endif
      integer dim_in,dim_out,nv
      real(p3dfft_type) Lz

      call p3dfft_cheby_many (XgYZ,dim_in,XYZg,dim_out,nv,Lz)

      end subroutine

! Chebyshev transform (2D R2C forward FFT + Chebyshev) for multiple variables
!========================================================
	subroutine p3dfft_cheby_many(in,dim_in,out,dim_out,nv,Lz)
!========================================================

	integer dim_in,dim_out,nv,j
	real(p3dfft_type) Lz
      	real(p3dfft_type), dimension(dim_in,nv), target	::	in
        complex(p3dfft_type), dimension (dim_out,nv), target :: out

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

      	real(p3dfft_type), dimension(nx_fft,     &
                                jisize,     &
                                kjsize), target	::	in

#ifdef STRIDE1
      	complex(p3dfft_type), dimension(nzc, &
                                jjsize,&
                                iisize), target	::	out
#else
      	complex(p3dfft_type), dimension(iisize, &
                                jjsize,    &
                                nzc), target	::	out
#endif
	real(p3dfft_type) Lz

	call p3dfft_cheby(in,out,Lz)

	end subroutine

! Chebyshev transform (2D R2C forward FFT + Chebyshev) for a single variable
!==============================================================
	subroutine p3dfft_cheby(in,out,Lz)
!========================================================

!     !	function args
      	real(p3dfft_type), dimension(nx_fft,     &
                                jisize,     &
                                kjsize), target	::	in

#ifdef STRIDE1
      	complex(p3dfft_type), dimension(nzc, &
                                jjsize,&
                                iisize), target	::	out
      	complex(p3dfft_type) Old, New, Tmp
#else
      	complex(p3dfft_type), dimension(iisize, &
                                jjsize,    &
                                nzc), target	::	out
      	complex(p3dfft_type),dimension(:,:),pointer :: ptrOld, ptrNew, ptrTmp
#endif
      	real(p3dfft_type) :: Lfactor,Lz
	integer k,nz,i,j

	nz = nzc

    	call p3dfft_ftran_r2c(in,out,'ffc')
     	out = out *(1.d0/(dble(nx_fft*ny_fft)*dble(nzc-1)))

! less tmp-memory version (but difficult to read)

     	Lfactor = 4.d0/dble(Lz)

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
      real(p3dfft_type), TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(p3dfft_type), TARGET :: XYZg(nzc,jjstart:jjend,iistart:iiend)
#else
      complex(p3dfft_type), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nzc)
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

      real(p3dfft_type), TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(p3dfft_type), TARGET :: XYZg(nzc,jjstart:jjend,iistart:iiend)
#else
      complex(p3dfft_type), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nzc)
#endif

      integer x,y,z,i,nx,ny,nz,ierr,dnz,dny
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
         call init_f_r2c(XgYZ,nx,buf2,nxhp,nx,jisize*kjsize)

         timers(5) = timers(5) - MPI_Wtime()

         call exec_f_r2c(XgYZ,nx,buf2,nxhp,nx,jisize*kjsize)

         timers(5) = timers(5) + MPI_Wtime()

      endif


! Exchange data in rows

      if(iproc .gt. 1) then

#ifdef DEBUG
	print *,taskid,': Calling fcomm1'
#endif
         call fcomm1(buf2,buf,timers(1),timers(6))

      else

#ifdef STRIDE1
#ifdef DEBUG
	print *,taskid,': Calling reorder_f1'
#endif
         call reorder_f1(buf2,buf,buf1)
#else
	call seg_copy_x(buf2,buf,1,nxhpc,0,nxhp,nxhpc,jisize,kjsize)
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
         call init_f_c(buf,1,nz, XYZg,1,nz,nz,jjsize)
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
              call exec_ctrans_r2_complex_same(buf,2*iisize*jjsize, 1,buf,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 's') then
               call init_strans_r2(buf,2*iisize*jjsize, 1, buf,2*iisize*jjsize, 1,nz,2*iisize*jjsize)

              timers(8) = timers(8) - MPI_Wtime()
              call exec_strans_r2_complex_same(buf,2*iisize*jjsize, 1,buf,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()
	    else if(op(3:3) .ne. 'n' .and. op(3:3) .ne. '0') then
		print *,'Unknown transform type: ',op(3:3)
		call MPI_Abort(MPI_COMM_WORLD,ierr)
            endif

	    call seg_copy_z(buf,XYZg,1,iisize,1,jjsize,1,nzcph,0,iisize,jjsize,nz)
	    call seg_copy_z(buf,XYZg,1,iisize,1,jjsize,nzcph+1,nzc,dnz,iisize,jjsize,nz)

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
              call exec_ctrans_r2_complex_same(XYZg,2*iisize*jjsize, 1,XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

            else if(op(3:3) == 's') then
               call init_strans_r2(XYZg,2*iisize*jjsize, 1, XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)

              timers(8) = timers(8) - MPI_Wtime()
              call exec_strans_r2_complex_same(XYZg,2*iisize*jjsize, 1,XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
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

	    dny = ny - nyc
	    call seg_copy_y(buf,buf1,1,nycph,0,iisize,ny,nyc,nz)
	    call seg_copy_y(buf,buf1,nycph+1,nyc,dny,iisize,ny,nyc,nz)

	    if(op(3:3) == 't' .or. op(3:3) == 'f') then
               call init_f_c(buf1,iisize*jjsize, 1,buf1,iisize*jjsize, 1,nz,iisize*jjsize)

              timers(8) = timers(8) - MPI_Wtime()
              call exec_f_c2_same(buf1,iisize*jjsize, 1,buf1,iisize*jjsize, 1,nz,iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 'c') then
               call init_ctrans_r2(buf1,2*iisize*jjsize, 1,buf1,2*iisize*jjsize, 1,nz,2*iisize*jjsize)

              timers(8) = timers(8) - MPI_Wtime()
              call exec_ctrans_r2_complex_same(buf1,2*iisize*jjsize, 1,buf1,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 's') then
               call init_strans_r2(buf1,2*iisize*jjsize, 1,buf1,2*iisize*jjsize, 1,nz,2*iisize*jjsize)

              timers(8) = timers(8) - MPI_Wtime()
              call exec_strans_r2_complex_same(buf1,2*iisize*jjsize, 1,buf1,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()
	    else if(op(3:3) /= 'n' .and. op(3:3) /= '0') then
		print *,'Unknown transform type: ',op(3:3)
		call MPI_Abort(MPI_COMM_WORLD,ierr)
            endif

	   call seg_copy_z(buf1,XYZg,1,iisize,1,jjsize,1,nzcph,0,iisize,jjsize,nz)
	   call seg_copy_z(buf1,XYZg,1,iisize,1,jjsize,nzcph+1,nzc,dnz,iisize,jjsize,nz)

	else

	    dny = ny - nyc
	    call seg_copy_y(buf,XYZg,1,nycph,0,iisize,ny,nyc,nz)
	    call seg_copy_y(buf,XYZg,nycph+1,nyc,dny,iisize,ny,nyc,nz)

            if(op(3:3) == 't' .or. op(3:3) == 'f') then
               call init_f_c(XYZg,iisize*jjsize, 1, XYZg,iisize*jjsize, 1,nz,iisize*jjsize)

              timers(8) = timers(8) - MPI_Wtime()
              call exec_f_c2_same(XYZg,iisize*jjsize, 1,XYZg,iisize*jjsize, 1,nz,iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

            else if(op(3:3) == 'c') then
               call init_ctrans_r2(XYZg,2*iisize*jjsize, 1, XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)

              timers(8) = timers(8) - MPI_Wtime()
              call exec_ctrans_r2_complex_same(XYZg,2*iisize*jjsize, 1,XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()

            else if(op(3:3) == 's') then
               call init_strans_r2(XYZg,2*iisize*jjsize, 1, XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)

              timers(8) = timers(8) - MPI_Wtime()
              call exec_strans_r2_complex_same(XYZg,2*iisize*jjsize, 1,XYZg,2*iisize*jjsize, 1,nz,2*iisize*jjsize)
              timers(8) = timers(8) + MPI_Wtime()
            else if(op(3:3) /= 'n' .and. op(3:3) /= '0') then
                print *,'Unknown transform type: ',op(3:3)
                call MPI_Abort(MPI_COMM_WORLD,ierr)
            endif


        endif
#endif

        timers(8) = timers(8) + MPI_Wtime()

      endif

!#ifdef DEBUG
!      print *,taskid,': Waiting at barrier'
!#endif
!
!      call mpi_barrier(mpi_comm_world,ierr)

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

  real (p3dfft_type), target :: rXgYZ (NX_fft, jistart:jiend, kjstart:kjend)
  real (p3dfft_type), target :: cXgYZ (NX_fft+2, jistart:jiend, kjstart:kjend)

!      complex(p3dfft_type), allocatable :: XYgZ(:,:,:)
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
  real(p3dfft_type) source(dim,nv)
  complex(p3dfft_type) dest(n/2+1,m,nv)

  do j=1,nv
    call exec_f_r2c(source(1,j),str1,dest(1,1,j),str2,n,m)
  enddo

  return
  end subroutine

         subroutine f_c1_many(A,str1,str2,n,m,dim,nv)

	   integer n,m,nv,j,str1,str2,dim
	   complex(p3dfft_type) A(dim,nv)

	 do j=1,nv
           call exec_f_c1(A(1,j),str1,str2,A(1,j),str1,str2,n,m)
         enddo

	 return
	 end subroutine

         subroutine ztran_f_same_many(A,str1,str2,n,m,dim,nv,op)

	   integer str1,str2,n,m,nv,j,ierr,dim
	   complex(p3dfft_type) A(dim,nv)
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
                 call exec_ctrans_r2_complex_same(A(1,j),2*str1,str2,A(1,j),2*str1,str2,n,2*m)
	      enddo
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(3:3) == 's') then
               call init_strans_r2(A,str1,str2,A,str1,str2,n,m)

              timers(8) = timers(8) - MPI_Wtime()
	      do j=1,nv
                 call exec_strans_r2_complex_same(A(1,j),2*str1,str2,A(1,j),2*str1,str2,n,2*m)
              enddo
              timers(8) = timers(8) + MPI_Wtime()
	    else if(op(3:3) .ne. 'n' .and. op(3:3) .ne. '0') then
		print *,'Unknown transform type: ',op(3:3)
		call MPI_Abort(MPI_COMM_WORLD,ierr)
            endif

	    return
	    end subroutine


