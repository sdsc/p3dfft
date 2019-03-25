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


! this is a C wrapper routine
!========================================================
      subroutine p3dfft_btran_c2r_many_w (XYZg,dim_in,XgYZ,dim_out,nv,op) BIND(C,NAME='p3dfft_btran_c2r_many')
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

      call p3dfft_btran_c2r_many (XYZg,dim_in,XgYZ,dim_out,nv,lcl_op)

      end subroutine

! Inverse C2R 3D FFT transform of multiple variables (nv)
!----------------------------------------------------------------------------
      subroutine p3dfft_btran_c2r_many (XYZg,dim_in,XgYZ,dim_out,nv,op)
!========================================================

      use fft_spec
      implicit none

      integer x,y,z,i,k,nx,ny,nz,ierr,dnz,nv,j,n1,n2,dim_in,dim_out,dny
      real(p3dfft_type),TARGET :: XgYZ(dim_out,nv)
#ifdef STRIDE1
      complex(p3dfft_type), TARGET :: XYZg(dim_in,nv)
#else
      complex(p3dfft_type), TARGET :: XYZg(dim_in,nv)
#endif

!      real(p3dfft_type),TARGET :: XgYZ(nx_fft,jisize,kjsize,nv)
!#ifdef STRIDE1
!      complex(p3dfft_type), TARGET :: XYZg(nzc,iisize,jjsize,nv)
!#else
!      complex(p3dfft_type), TARGET :: XYZg(iisize,jjsize,nzc,nv)
!#endif

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

! Allocate work array

      if(nv .gt. nv_preset) then
        nv_preset = nv
	if(allocated(buf)) deallocate(buf)
        allocate(buf(nxhp*jisize*(kjsize+padi)*nv))
        deallocate(buf1,buf2)
#ifdef USE_EVEN
        n1 = nv * IfCntMax * iproc /(p3dfft_type*2)
        n2 = nv * KfCntMax * jproc / (p3dfft_type*2)
        n1 = max(n1,n2)
        allocate(buf1(n1))
        allocate(buf2(n1))
#else
        allocate(buf1(nm*nv))
        allocate(buf2(nm*nv))
#endif
      endif

! FFT Tranform (C2C) in Z for all x and y

      if(jproc .gt. 1) then

#ifdef STRIDE1
         call init_b_c(XYZg, 1,nz, buf, 1, nz,nz,jjsize)
         call bcomm1_trans_many(XYZg,buf,dim_in,nv,op,timers(3),timers(9))
#else

         if(OW .and. nz .eq. nzc) then

            if(iisize*jjsize .gt. 0) then
	       call ztran_b_same_many(XYZg,iisize*jjsize,1,nz,iisize*jjsize,dim_in,nv,op)
            endif
            call bcomm1_many(XYZg,buf,dim_in,nv,timers(3),timers(9))

         else

           if(iisize*jjsize .gt. 0) then
              if(op(1:1) == 'n' .or. op(1:1) == '0') then
                  call bcomm1_many(XYZg,buf,dim_in,nv,timers(3),timers(9))
	      else

	        dnz = nz - nzc
		call seg_copy_z_b_many(XYZg,buf,1,iisize,1,jjsize,1,nzhc,0,iisize,jjsize,nz,dim_in,nv)
		call seg_copy_z_b_many(XYZg,buf,1,iisize,1,jjsize,nz-nzhc+1,nz,-dnz,iisize,jjsize,nz,dim_in,nv)
		call seg_zero_z_many(buf,iisize,jjsize,nzhc+1,nz-nzhc,nz,iisize*jjsize*nz,nv)

		call ztran_b_same_many(buf,iisize*jjsize,1,nz,iisize*jjsize,iisize*jjsize*nz,nv,op)
                call bcomm1_many(buf,buf,iisize*jjsize*nz,nv,timers(3),timers(9))
              endif

            endif
         endif

#endif

      else
            timers(9) = timers(9) - MPI_Wtime()

#ifdef STRIDE1
         call reorder_trans_b1_many(XYZg,buf,buf2,dim_in,nv,op)
#else
         Nl = iisize*jjsize*nzc
         if(OW .and. nz .eq. nzc) then

            call ztran_b_same_many(XYZg,iisize*jjsize,1,nz,iisize*jjsize,dim_in,nv,op)
            call ar_copy_many(XYZg,dim_in,buf,iisize*jjsize*nz,Nl,nv)

         else
           if(iisize*jjsize .gt. 0) then

	      if(op(1:1) == 'n' .or. op(1:1) == '0') then

  	         call seg_copy_z_b_many(XYZg,buf,1,iisize,1,jjsize,1,nz,0,iisize,jjsize,nz,dim_in,nv)
              else

	        dnz = nz - nzc
	        call seg_copy_z_b_many(XYZg,buf,1,iisize,1,jjsize,1,nzhc,0,iisize,jjsize,nz,dim_in,nv)
		call seg_copy_z_b_many(XYZg,buf,1,iisize,1,jjsize,nz-nzhc+1,nz,-dnz,iisize,jjsize,nz,dim_in,nv)
		call seg_zero_z_many(buf,iisize,jjsize,nzhc+1,nz-nzhc,nz,iisize*jjsize*nz,nv)

		call ztran_b_same_many(buf,iisize*jjsize,1,nz,iisize*jjsize,iisize*jjsize*nz,nv,op)

	        dny = ny - nyc
	        call seg_copy_y_b_many(buf1,buf,1,nycph,0,iisize,nyc,ny,nz,iisize*nyc*nz,nv)
		call seg_copy_y_b_many(buf1,buf,ny-nycph+1,ny,-dny,iisize,nyc,ny,nz,iisize*nyc*nz,nv)
		call seg_zero_y_many(buf,nycph+1,ny-nycph,iisize,ny,nz,nv)

              endif
	    endif
         endif
#endif

            timers(9) = timers(9) + MPI_Wtime()

      endif

! Exhange in columns if needed

!
! FFT Transform (C2C) in y dimension for all x, one z-plane at a time
!


      if(iisize * kjsize .gt. 0) then

#ifdef STRIDE1
         call init_b_c(buf,1,ny,buf,1,ny,ny,iisize*kjsize)

         timers(10) = timers(10) - MPI_Wtime()

         call b_c1_many(buf,1,ny,ny,iisize*kjsize,iisize*kjsize*ny,nv)

         timers(10) = timers(10) + MPI_Wtime()


#else
         call init_b_c(buf,iisize,1,buf,iisize,1,ny,iisize)

         timers(10) = timers(10) - MPI_Wtime()
         do z=1,kjsize * nv

            call btran_y_zplane(buf,z-1,iisize,kjsize,iisize,1, &
                                buf,z-1,iisize,kjsize,iisize,1,ny,iisize)

         enddo
         timers(10) = timers(10) + MPI_Wtime()
#endif
      endif

      if(iproc .gt. 1) then
         call bcomm2_many(buf,buf1,nv,timers(4),timers(11))
      else
#ifdef STRIDE1
         call reorder_b2_many(buf,buf1,nv)
#else
	Nl = jisize*kjsize*nxhp
	call seg_copy_x_b_many(buf,buf1,1,nxhpc,0,nxhpc,nxhp,jisize,kjsize,nxhpc*jisize*kjsize,nv)
	call seg_zero_x_many(buf1,nxhpc+1,nxhp,nxhp,jisize,kjsize,nv)
#endif
      endif


! Perform Complex-to-real FFT in x dimension for all y and z
       if(jisize * kjsize .gt. 0) then

          call init_b_c2r(buf1,nxhp,XgYZ,nx,nx,jisize*kjsize)

          timers(12) = timers(12) - MPI_Wtime()
          call b_c2r_many(buf1,nxhp,XgYZ,nx,nx,jisize*kjsize,dim_out,nv)
          timers(12) = timers(12) + MPI_Wtime()
       endif

!	deallocate(buf)

!      do j=1,nv
!        print *,'Exiting btran: j=',j
!        call print_buf_real(XgYZ(1,j),nx,jisize,kjsize)
!      enddo

!      call mpi_barrier(mpi_comm_world,ierr)

      return
      end subroutine

!========================================================
subroutine b_c2r_many(A,str1,B,str2,n,m,dim,nv)
!========================================================

	   integer str1,str2,n,m,nv,j,dim
	   complex(p3dfft_type) A(n/2+1,m,nv)
	   real(p3dfft_type) B(dim,nv)

	   do j=1,nv
	      call exec_b_c2r(A(1,1,j),str1,B(1,j),str2,n,m)
           enddo

	   return
	   end subroutine

!========================================================
subroutine ztran_b_same_many(A,str1,str2,n,m,dim,nv,op)
!========================================================

	   integer str1,str2,n,m,nv,j,ierr,dim
	   complex(p3dfft_type) A(dim,nv)
	   character(len=3) op

	    if(op(1:1) == 't' .or. op(1:1) == 'f') then
               call init_b_c(A,str1,str2,A,str1,str2,n,m)

              timers(8) = timers(8) - MPI_Wtime()
	      do j=1,nv
                 call exec_b_c2_same(A(1,j),str1,str2,A(1,j),str1,str2,n,m)
              enddo
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(1:1) == 'c') then
               call init_ctrans_r2(A,str1,str2,A,str1,str2,n,m)

              timers(8) = timers(8) - MPI_Wtime()
	      do j=1,nv
                 call exec_ctrans_r2_complex_same(A(1,j),2*str1,str2,A(1,j),2*str1,str2,n,2*m)
	      enddo
              timers(8) = timers(8) + MPI_Wtime()

	    else if(op(1:1) == 's') then
               call init_strans_r2(A,str1,str2,A,str1,str2,n,m)

              timers(8) = timers(8) - MPI_Wtime()
	      do j=1,nv
                 call exec_strans_r2_complex_same(A(1,j),2*str1,str2,A(1,j),2*str1,str2,n,2*m)
              enddo
              timers(8) = timers(8) + MPI_Wtime()
	    else if(op(1:1) .ne. 'n' .and. op(1:1) .ne. '0') then
		print *,'Unknown transform type: ',op(1:1)
		call MPI_Abort(MPI_COMM_WORLD,ierr)
            endif

	    return
	    end subroutine

	    subroutine b_c1_many(A,str1,str2,n,m,dim,nv)

	   integer n,m,nv,j,str1,str2,dim
	   complex(p3dfft_type) A(dim,nv)

  	   do j=1,nv
             call exec_b_c1(A(1,j),str1,str2,A(1,j),str1,str2,n,m)
           enddo

	   return
	   end subroutine

! This is a C wrapper routine
!========================================================
      subroutine p3dfft_btran_c2r_w (XYZg,XgYZ,op) BIND(C,NAME='p3dfft_btran_c2r')
!========================================================
      use, intrinsic :: iso_c_binding
      real(p3dfft_type), TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(p3dfft_type), TARGET :: XYZg(nzc,iistart:iiend,jjstart:jjend)
#else
      complex(p3dfft_type), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nzc)
#endif
      character, dimension(*), target :: op
      character(4), pointer :: lcl_op
      call c_f_pointer(c_loc(op), lcl_op)

      call p3dfft_btran_c2r (XYZg,XgYZ,lcl_op)

      end subroutine

! Inverse C2R 3D FFT transform of a single variable
!----------------------------------------------------------------------------
      subroutine p3dfft_btran_c2r (XYZg,XgYZ,op)
!========================================================

      use fft_spec
      implicit none

      real(p3dfft_type),TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(p3dfft_type), TARGET :: XYZg(nzc,jjstart:jjend,iistart:iiend)
#else
      complex(p3dfft_type), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nzc)
#endif

      integer x,y,z,i,k,nx,ny,nz,ierr,dnz,dny
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

! Allocate work array

!      allocate(buf(nxhp,jistart:jiend,kjstart:kjend+padd))

! FFT Tranform (C2C) in Z for all x and y

      if(jproc .gt. 1) then

#ifdef STRIDE1
         call init_b_c(XYZg, 1,nz, buf, 1, nz,nz,jjsize)
         call bcomm1_trans(XYZg,buf,op,timers(3),timers(9))
#else

         if(OW .and. nz .eq. nzc) then

            if(iisize*jjsize .gt. 0) then
		if(op(1:1) == 't' .or. op(1:1) == 'f') then
                   call init_b_c(XYZg, iisize*jjsize, 1, &
                                 XYZg, iisize*jjsize, 1,nz,iisize*jjsize)

                   timers(9) = timers(9) - MPI_Wtime()
                   call exec_b_c2_same(XYZg, iisize*jjsize,1, XYZg, &
				iisize*jjsize, 1,nz,iisize*jjsize)
                   timers(9) = timers(9) + MPI_Wtime()
 		else if(op(1:1) == 'c') then
	           call init_ctrans_r2 (XYZg, 2*iisize*jjsize, 1, &
					XYZg, 2*iisize*jjsize, 1, &
					nz, 2*iisize*jjsize)
                   call exec_ctrans_r2_complex_same (XYZg, 2*iisize*jjsize, 1, &
 					XYZg, 2*iisize*jjsize, 1, &
					nz, 2*iisize*jjsize)
 		else if(op(1:1) == 's') then
	           call init_strans_r2 (XYZg, 2*iisize*jjsize, 1, &
					XYZg, 2*iisize*jjsize, 1, &
					nz, 2*iisize*jjsize)
                   call exec_strans_r2_complex_same (XYZg, 2*iisize*jjsize, 1, &
 					XYZg, 2*iisize*jjsize, 1, &
					nz, 2*iisize*jjsize)
	        else if(op(1:1) /= 'n' .and. op(1:1) /= '0') then
		   print *,taskid,'Unknown transform type: ',op(1:1)
		   call MPI_abort(MPI_COMM_WORLD,ierr)
		endif
            endif
            call bcomm1(XYZg,buf,timers(3),timers(9))

         else

           if(iisize*jjsize .gt. 0) then
              if(op(1:1) == 'n' .or. op(1:1) == '0') then
                  call bcomm1(XYZg,buf,timers(3),timers(9))
	      else

	        dnz = nz - nzc
		call seg_copy_z(XYZg,buf,1,iisize,1,jjsize,1,nzcph,0,iisize,jjsize,nz)
		call seg_copy_z(XYZg,buf,1,iisize,1,jjsize,nz-nzcph+1,nz,-dnz,iisize,jjsize,nz)
		call seg_zero_z(buf,iisize,jjsize,nzcph+1,nz-nzcph,nz)

		 if(op(1:1) == 't' .or. op(1:1) == 'f') then
                   call init_b_c(buf, iisize*jjsize, 1, &
				 buf, iisize*jjsize, 1,nz,iisize*jjsize)

                   timers(9) = timers(9) - MPI_Wtime()
    	           call exec_b_c2_same(buf, iisize*jjsize,1, &
                                 buf, iisize*jjsize, 1,nz,iisize*jjsize)
                   timers(9) = timers(9) + MPI_Wtime()
		 else if(op(1:1) == 'c') then
	           call init_ctrans_r2 (buf, 2*iisize*jjsize, 1, &
					buf, 2*iisize*jjsize, 1, &
					nz, 2*iisize*jjsize)
                   call exec_ctrans_r2_complex_same (buf, 2*iisize*jjsize, 1, &
 					buf, 2*iisize*jjsize, 1, &
					nz, 2*iisize*jjsize)
		 else if(op(1:1) == 's') then
	           call init_strans_r2 (buf, 2*iisize*jjsize, 1, &
					buf, 2*iisize*jjsize, 1, &
					nz, 2*iisize*jjsize)
                   call exec_strans_r2_complex_same (buf, 2*iisize*jjsize, 1, &
 					buf, 2*iisize*jjsize, 1, &
					nz, 2*iisize*jjsize)
	         else if(op(1:1) /= 'n' .and. op(1:1) /= '0') then
		   print *,taskid,'Unknown transform type: ',op(1:1)
		   call MPI_abort(MPI_COMM_WORLD,ierr)
		 endif

                 call bcomm1(buf,buf,timers(3),timers(9))
              endif

            endif
         endif

#endif


      else

            timers(9) = timers(9) - MPI_Wtime()

#ifdef STRIDE1
         call reorder_trans_b1(XYZg,buf,buf2,op)
#else
         Nl = iisize*jjsize*nzc
         if(OW .and. nz .eq. nzc) then

  	    if(op(1:1) == 't' .or. op(1:1) == 'f') then
               call init_b_c(XYZg, iisize*jjsize, 1, &
			     XYZg, iisize*jjsize, 1,nz,iisize*jjsize)
               call exec_b_c2_same(XYZg, iisize*jjsize, 1, &
			      XYZg, iisize*jjsize, 1,nz,iisize*jjsize)
   	    else if(op(1:1) == 'c') then
	       call init_ctrans_r2 (XYZg, 2*iisize*jjsize, 1, &
			            XYZg, 2*iisize*jjsize, 1, &
				    nz, 2*iisize*jjsize)
               call exec_ctrans_r2_complex_same (XYZg, 2*iisize*jjsize, 1, &
 				    XYZg, 2*iisize*jjsize, 1, &
				    nz, 2*iisize*jjsize)
     	    else if(op(1:1) == 's') then
	       call init_strans_r2 (XYZg, 2*iisize*jjsize, 1, &
		    	            XYZg, 2*iisize*jjsize, 1, &
			            nz, 2*iisize*jjsize)
               call exec_strans_r2_complex_same (XYZg, 2*iisize*jjsize, 1, &
				    XYZg, 2*iisize*jjsize, 1, &
				    nz, 2*iisize*jjsize)
	    else if(op(1:1) /= 'n' .and. op(1:1) /= '0') then
		print *,taskid,'Unknown transform type: ',op(1:1)
   	        call MPI_abort(MPI_COMM_WORLD,ierr)
	    endif
            call ar_copy(XYZg,buf,Nl)

         else
           if(iisize*jjsize .gt. 0) then

	        dnz = nz - nzc
		call seg_copy_z(XYZg,buf1,1,iisize,1,jjsize,1,nzcph,0,iisize,jjsize,nz)
		call seg_copy_z(XYZg,buf1,1,iisize,1,jjsize,nz-nzcph+1,nz,-dnz,iisize,jjsize,nz)
		call seg_zero_z(buf1,iisize,jjsize,nzcph+1,nz-nzcph,nz)

    	         if(op(1:1) == 't' .or. op(1:1) == 'f') then
                    call init_b_c(buf1, iisize*jjsize, 1,  &
			     buf1, iisize*jjsize, 1,nz,iisize*jjsize)
                    call exec_b_c2_same(buf1, iisize*jjsize, 1, &
			      buf1, iisize*jjsize, 1,nz,iisize*jjsize)
   	         else if(op(1:1) == 'c') then
	            call init_ctrans_r2 (buf1, 2*iisize*jjsize, 1, &
			            buf1, 2*iisize*jjsize, 1, &
				    nz, 2*iisize*jjsize)
                    call exec_ctrans_r2_complex_same (buf1, 2*iisize*jjsize, 1, &
 				    buf1, 2*iisize*jjsize, 1, &
				    nz, 2*iisize*jjsize)
     	         else if(op(1:1) == 's') then
	            call init_strans_r2 (buf1, 2*iisize*jjsize, 1, &
		    	            buf1, 2*iisize*jjsize, 1, &
			            nz, 2*iisize*jjsize)
                    call exec_strans_r2_complex_same (buf1, 2*iisize*jjsize, 1, &
				    buf1, 2*iisize*jjsize, 1, &
				    nz, 2*iisize*jjsize)
                 else if(op(1:1) /= 'n' .and. op(1:1) /= '0') then
		    print *,taskid,'Unknown transform type: ',op(1:1)
	            call MPI_abort(MPI_COMM_WORLD,ierr)
 	         endif

 		 dny = ny - nyc
		 call seg_copy_y(buf1,buf,1,nycph,0,iisize,nyc,ny,nz)
		 call seg_copy_y(buf1,buf,ny-nycph+1,ny,-dny,iisize,nyc,ny,nz)
		 call seg_zero_y(buf,nycph+1,ny-nycph,iisize,ny,nz)
	    endif
         endif
#endif

            timers(9) = timers(9) + MPI_Wtime()

      endif

! Exhange in columns if needed

!
! FFT Transform (C2C) in y dimension for all x, one z-plane at a time
!


      if(iisize * kjsize .gt. 0) then

#ifdef STRIDE1
         call init_b_c(buf,1,ny,buf,1,ny,ny,iisize*kjsize)

         timers(10) = timers(10) - MPI_Wtime()

         call exec_b_c1(buf,1,ny,buf,1,ny,ny,iisize*kjsize)

         timers(10) = timers(10) + MPI_Wtime()

#else
         call init_b_c(buf,iisize,1,buf,iisize,1,ny,iisize)

         timers(10) = timers(10) - MPI_Wtime()
         do z=kjstart,kjend

            call btran_y_zplane(buf,z-kjstart,iisize,kjsize,iisize,1, &
                                buf,z-kjstart,iisize,kjsize,iisize,1,ny,iisize)

         enddo
         timers(10) = timers(10) + MPI_Wtime()
#endif
      endif

#ifdef STRIDE1
      if(iproc .gt. 1) then
         call bcomm2(buf,buf1,timers(4),timers(11))
      else
         call reorder_b2(buf,buf1)
      endif
! Perform Complex-to-real FFT in x dimension for all y and z
      if(jisize * kjsize .gt. 0) then

         call init_b_c2r(buf1,nxhp,XgYZ,nx,nx,jisize*kjsize)

         timers(12) = timers(12) - MPI_Wtime()

         call exec_b_c2r(buf1,nxhp,XgYZ,nx,nx,jisize*kjsize)
         timers(12) = timers(12) + MPI_Wtime()

      endif
#else
      if(iproc .gt. 1) then
         call bcomm2(buf,buf,timers(4),timers(11))
! Perform Complex-to-real FFT in x dimension for all y and z
         if(jisize * kjsize .gt. 0) then

            call init_b_c2r(buf,nxhp,XgYZ,nx,nx,jisize*kjsize)

            timers(12) = timers(12) - MPI_Wtime()

            call exec_b_c2r(buf,nxhp,XgYZ,nx,nx,jisize*kjsize)
            timers(12) = timers(12) + MPI_Wtime()
         endif
      else
! Perform Complex-to-real FFT in x dimension for all y and z
         if(jisize * kjsize .gt. 0) then

	    call seg_copy_x(buf,buf1,1,nxhpc,0,nxhpc,nxhp,jisize,kjsize)
	    call seg_zero_x(buf1,nxhpc+1,nxhp,nxhp,jisize,kjsize)

            call init_b_c2r(buf1,nxhp,XgYZ,nx,nx,jisize*kjsize)

            timers(12) = timers(12) - MPI_Wtime()

            call exec_b_c2r(buf1,nxhp,XgYZ,nx,nx,jisize*kjsize)
            timers(12) = timers(12) + MPI_Wtime()
         endif
       endif

#endif

!      call mpi_barrier(mpi_comm_world,ierr)

      return
      end subroutine


