! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2013 Dmitry Pekurovsky
!    Copyright (C) 2006-2013 University of California
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
      subroutine p3dfft_btran_c2r (XYZg,XgYZ,op)
!========================================================

      use fft_spec
      implicit none

      real(mytype),TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(mytype), TARGET :: XYZg(nzc,iistart:iiend,jjstart:jjend)
#else
      complex(mytype), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nzc)
#endif

      integer x,y,z,i,k,nx,ny,nz,ierr,dnz
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
         call bcomm1_trans(XYZg,buf,buf2,op,timers(3),timers(9))
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
                   call exec_ctrans_r2_same (XYZg, 2*iisize*jjsize, 1, &
 					XYZg, 2*iisize*jjsize, 1, &
					nz, 2*iisize*jjsize)
 		else if(op(1:1) == 's') then	
	           call init_strans_r2 (XYZg, 2*iisize*jjsize, 1, &
					XYZg, 2*iisize*jjsize, 1, &
					nz, 2*iisize*jjsize)
                   call exec_strans_r2_same (XYZg, 2*iisize*jjsize, 1, &
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
		call seg_copy_z(XYZg,buf,1,iisize,1,jjsize,1,nzhc,0,iisize,jjsize,nz)
		call seg_copy_z(XYZg,buf,1,iisize,1,jjsize,nz-nzhc+1,nz,-dnz,iisize,jjsize,nz)
		call seg_zero_z(buf,iisize,jjsize,nzhc+1,nz-nzhc,nz)

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
                   call exec_ctrans_r2_same (buf, 2*iisize*jjsize, 1, & 
 					buf, 2*iisize*jjsize, 1, &
					nz, 2*iisize*jjsize)
		 else if(op(1:1) == 's') then
	           call init_strans_r2 (buf, 2*iisize*jjsize, 1, &
					buf, 2*iisize*jjsize, 1, &
					nz, 2*iisize*jjsize)
                   call exec_strans_r2_same (buf, 2*iisize*jjsize, 1, & 
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
               call exec_ctrans_r2_same (XYZg, 2*iisize*jjsize, 1, &
 				    XYZg, 2*iisize*jjsize, 1, &
				    nz, 2*iisize*jjsize)
     	    else if(op(1:1) == 's') then
	       call init_strans_r2 (XYZg, 2*iisize*jjsize, 1, &
		    	            XYZg, 2*iisize*jjsize, 1, & 
			            nz, 2*iisize*jjsize)
               call exec_strans_r2_same (XYZg, 2*iisize*jjsize, 1, &
				    XYZg, 2*iisize*jjsize, 1, &
				    nz, 2*iisize*jjsize)
	    else if(op(1:1) /= 'n' .and. op(1:1) /= '0') then
		print *,taskid,'Unknown transform type: ',op(1:1)
   	        call MPI_abort(MPI_COMM_WORLD,ierr)
	    endif
            call ar_copy(XYZg,buf,Nl)

         else
           if(iisize*jjsize .gt. 0) then

	      if(op(1:1) == 'n' .or. op(1:1) == '0') then

		call seg_copy_z(XYZg,buf,1,iisize,1,jjsize,1,nz,0,iisize,jjsize,nz)
              else

	        dnz = nz - nzc
		call seg_copy_z(XYZg,buf,1,iisize,1,jjsize,1,nzhc,0,iisize,jjsize,nz)
		call seg_copy_z(XYZg,buf,1,iisize,1,jjsize,nz-nzhc+1,nz,-dnz,iisize,jjsize,nz)
		call seg_zero_z(buf,iisize,jjsize,nzhc+1,nz-nzhc,nz)

 
    	         if(op(1:1) == 't' .or. op(1:1) == 'f') then
                    call init_b_c(buf, iisize*jjsize, 1,  &
			     buf, iisize*jjsize, 1,nz,iisize*jjsize)
                    call exec_b_c2_same(buf, iisize*jjsize, 1, &
			      buf, iisize*jjsize, 1,nz,iisize*jjsize)
   	         else if(op(1:1) == 'c') then
	            call init_ctrans_r2 (buf, 2*iisize*jjsize, 1, &
			            buf, 2*iisize*jjsize, 1, &
				    nz, 2*iisize*jjsize)
                    call exec_ctrans_r2_same (buf, 2*iisize*jjsize, 1, & 
 				    buf, 2*iisize*jjsize, 1, &
				    nz, 2*iisize*jjsize)
     	         else if(op(1:1) == 's') then
	            call init_strans_r2 (buf, 2*iisize*jjsize, 1, &
		    	            buf, 2*iisize*jjsize, 1, &
			            nz, 2*iisize*jjsize)
                    call exec_strans_r2_same (buf, 2*iisize*jjsize, 1, & 
				    buf, 2*iisize*jjsize, 1, &
				    nz, 2*iisize*jjsize) 
                 else if(op(1:1) /= 'n' .and. op(1:1) /= '0') then
		    print *,taskid,'Unknown transform type: ',op(1:1)
	            call MPI_abort(MPI_COMM_WORLD,ierr)
	         endif
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
         call bcomm2(buf,buf,timers(4),timers(11))
      else
         call reorder_b2(buf,buf)
      endif
! Perform Complex-to-real FFT in x dimension for all y and z
      if(jisize * kjsize .gt. 0) then

         call init_b_c2r(buf,nxhp,XgYZ,nx,nx,jisize*kjsize)

         timers(12) = timers(12) - MPI_Wtime()

         call exec_b_c2r(buf,nxhp,XgYZ,nx,nx,jisize*kjsize)
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

            call init_b_c2r(buf,nxhp,XgYZ,nx,nx,jisize*kjsize)

            timers(12) = timers(12) - MPI_Wtime()

            call exec_b_c2r(buf,nxhp,XgYZ,nx,nx,jisize*kjsize)
            timers(12) = timers(12) + MPI_Wtime()
         endif
       endif	

#endif

      return
      end subroutine


