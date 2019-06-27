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
!----------------------------------------------------------------------

! This file contains routines for transposing data locally in memory
! They are used then either iproc=1 or jproc=1 (or both)


! This routine is called only when jproc=1, and only when stride1 is used
! transform backward in Z and transpose array in memory

!=============================================================
      subroutine reorder_trans_b1_many(A,B,C,dim,nv,op)
!=============================================================

      use fft_spec

      character(len=3) op
      integer x,y,z,iy,iz,y2,z2,ierr,dnz,dny,nv,j,dim
      complex(p3dfft_type) A(dim,nv)
      complex(p3dfft_type) B(ny_fft,iisize,nz_fft,nv)
      complex(p3dfft_type) C(nz_fft,nyc)

      do j=1,nv
         call reorder_trans_b1(A(1,j),B(1,1,1,j),C,op)
      enddo

      return
      end subroutine

!=============================================================
      subroutine reorder_trans_b1(A,B,C,op)
!=============================================================

      use fft_spec

      complex(p3dfft_type) A(nzc,nyc,iisize)
      complex(p3dfft_type) B(ny_fft,iisize,nz_fft)
      complex(p3dfft_type) C(nz_fft,nyc)
      integer x,y,z,iy,iz,y2,z2,ierr,dnz,dny
      character(len=3) op

      dny = ny_fft - nyc
      dnz = nz_fft - nzc
      if(op(1:1) == '0' .or. op(1:1) == 'n') then

!$OMP parallel do private(x,y,z,y2,z2,iy,iz,C)
         do x=1,iisize
            do y=1,nycph,NBy2
               y2 = min(y+NBy2-1,nycph)

   	       do z=1,nzcph,NBz
	          z2 = min(z+NBz-1,nzcph)
                    do iy=y,y2
                        do iz=z,z2
			   B(iy,x,iz) = A(iz,iy,x)
                        enddo
                     enddo
                enddo

   	       do z=nzcph+dnz+1,nz_fft,NBz
	          z2 = min(z+NBz-1,nz_fft)
                    do iy=y,y2
                        do iz=z,z2
			   B(iy,x,iz) = A(iz-dnz,iy,x)
                        enddo
                     enddo
                enddo

              enddo

            do y=nycph+1,nyc,NBy2
               y2 = min(y+NBy2-1,nyc)

   	       do z=1,nzcph,NBz
	          z2 = min(z+NBz-1,nzcph)
                    do iy=y,y2
                        do iz=z,z2
			   B(iy+dny,x,iz) = A(iz,iy,x)
                        enddo
                     enddo
                enddo

   	       do z=nzcph+dnz+1,nz_fft,NBz
	          z2 = min(z+NBz-1,nz_fft)
                    do iy=y,y2
                        do iz=z,z2
			   B(iy+dny,x,iz) = A(iz-dnz,iy,x)
                        enddo
                     enddo
                enddo

              enddo

          enddo

	  do z=nzcph+1,nzcph+dnz
	     do x=1,iisize
                do y=1,ny_fft
	 	   B(y,x,z) = 0
                enddo
             enddo
          enddo
	else

!$OMP parallel do private(x,y,z,y2,z2,iy,iz,C)
           do x=1,iisize
	      do y=1,nyc
	         do z=1,nzcph
		    C(z,y) = A(z,y,x)
		 enddo
	         do z=nzcph+1,nzcph+dnz
		    C(z,y) = 0.
		 enddo
	         do z=nzcph+dnz+1,nz_fft
		    C(z,y) = A(z-dnz,y,x)
		 enddo
	      enddo

	      if(op(1:1) == 't' .or. op(1:1) == 'f') then
                 call exec_b_c2_same_serial(C, 1,nz_fft, &
				  C, 1,nz_fft,nz_fft,nyc)
 	      else if(op(1:1) == 'c') then
                 call exec_ctrans_r2_complex_same(C, 2,2*nz_fft, &
				  C, 2,2*nz_fft,nz_fft,nyc)
 	      else if(op(1:1) == 's') then
                 call exec_strans_r2_complex_same(C, 2,2*nz_fft, &
				  C, 2,2*nz_fft,nz_fft,nyc)
              else
	         print *,taskid,'Unknown transform type: ',op(1:1)
	         call MPI_abort(MPI_COMM_WORLD,ierr)
	      endif

              do y=1,nycph,NBy2
                 y2 = min(y+NBy2-1,nycph)
     	         do z=1,nz_fft,NBz
	            z2 = min(z+NBz-1,nz_fft)
                       do iy=y,y2
                          do iz=z,z2
			     B(iy,x,iz) = C(iz,iy)
                          enddo
                       enddo
                  enddo
              enddo
              do y=nycph+1,nyc,NBy2
                 y2 = min(y+NBy2-1,nyc)
     	         do z=1,nz_fft,NBz
	            z2 = min(z+NBz-1,nz_fft)
                       do iy=y,y2
                          do iz=z,z2
			     B(iy+dny,x,iz) = C(iz,iy)
                          enddo
                       enddo
                  enddo
              enddo
          enddo
     endif

     do z=1,nz_fft
        do x=1,iisize
           do y=nycph+1,nycph+dny
   	      B(y,x,z) = 0.
	   enddo
        enddo
      enddo

      return
      end subroutine

! This routine is called only when iproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Y

!=============================================================
      subroutine reorder_b2_many(A,B,nv)
!=============================================================

      implicit none

      complex(p3dfft_type) B(nxhp,ny_fft,kjsize,nv)
      complex(p3dfft_type) A(ny_fft,nxhpc,kjsize,nv)
      complex(p3dfft_type) tmp(nxhpc,ny_fft)
      integer x,y,z,iy,x2,ix,y2,nv,j

!$OMP parallel do private(x,y,z,y2,x2,iy,ix,tmp) collapse(2)

      do j=1,nv
      do z=1,kjsize
         do x=1,nxhpc,nbx
            x2 = min(x+nbx-1,nxhpc)
            do y=1,ny_fft,nby1
               y2 = min(y+nby1-1,ny_fft)
               do ix=x,x2
                  do iy = y,y2
                     tmp(ix,iy) = A(iy,ix,z,j)
                  enddo
               enddo
            enddo
         enddo

         do y=1,ny_fft
            do x=1,nxhpc
               B(x,y,z,j) = tmp(x,y)
            enddo
	    do x=nxhpc+1,nxhp
	       B(x,y,z,j) = 0.
	    enddo
         enddo
      enddo
      enddo

      return
      end subroutine

! This routine is called only when iproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Y

!=============================================================
      subroutine reorder_f1_many(A,B,tmp,nv)
!=============================================================

      implicit none

      complex(p3dfft_type) A(nxhp,ny_fft,kjsize,nv)
      complex(p3dfft_type) B(ny_fft,nxhpc,kjsize,nv)
      integer x,y,z,iy,x2,ix,y2,dnx,nv,j
      complex(p3dfft_type) tmp(ny_fft,nxhpc,kjsize)
!      complex(p3dfft_type), allocatable :: tmp(:,:)

!      allocate(tmp(ny_fft,nxhpc))

!!$OMP parallel do private(x,y,z,y2,x2,iy,ix,tmp) collapse(2)
       do j=1,nv
       do z=1,kjsize
         do y=1,ny_fft,nby1
            y2 = min(y+nby1-1,ny_fft)
            do x=1,nxhpc,nbx
               x2 = min(x+nbx-1,nxhpc)
               do iy = y,y2
                  do ix=x,x2
                     tmp(iy,ix,z) = A(ix,iy,z,j)
                  enddo
               enddo
            enddo
         enddo

         do x=1,nxhpc
            do y=1,ny_fft
               B(y,x,z,j) = tmp(y,x,z)
            enddo
         enddo
      enddo

      enddo

      return
      end subroutine

! This routine is called only when jproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Z

!=============================================================
      subroutine reorder_trans_f2_many(A,B,C,dim,nv,op)
!=============================================================

      use fft_spec
      implicit none

      integer x,y,z,iy,iz,y2,z2,ierr,dnz,dny,nv,j,dim
      complex(p3dfft_type) B(dim,nv)
      complex(p3dfft_type) A(ny_fft,iisize,nz_fft,nv)
      complex(p3dfft_type) C(nz_fft,nyc)
      character(len=3) op

      do j=1,nv
         call reorder_trans_f2(A(1,1,1,j),B(1,j),C,op)
      enddo

      return
      end subroutine

!=============================================================
      subroutine reorder_trans_f2(A,B,C,op)
!=============================================================

      use fft_spec
      implicit none

      integer x,y,z,iy,iz,y2,z2,ierr,dnz,dny
      complex(p3dfft_type) A(ny_fft,iisize,nz_fft)
      complex(p3dfft_type) B(nzc,nyc,iisize)
      complex(p3dfft_type) C(nz_fft,nyc)
      character(len=3) op

      dnz = nz_fft - nzc
      dny = ny_fft - nyc
      if(op(3:3) == '0' .or. op(3:3) == 'n') then

!$OMP parallel do private(x,y,z,y2,z2,iy,iz,C)
         do x=1,iisize
            do z=1,nzcph,NBz
	       z2 = min(z+NBz-1,nzcph)

               do y=1,nycph,NBy2
                  y2 = min(y+NBy2-1,nycph)
                  do iz=z,z2
                     do iy=y,y2
  		        B(iz,iy,x) = A(iy,x,iz)
                     enddo
                  enddo
               enddo

               do y=nycph+1,nyc,NBy2
                  y2 = min(y+NBy2-1,ny_fft)
                  do iz=z,z2
                     do iy=y,y2
  		        B(iz,iy,x) = A(iy+dny,x,iz)
                     enddo
                  enddo
               enddo
            enddo
            do z=nzcph+1,nzc,NBz
	       z2 = min(z+NBz-1,nz_fft)

               do y=1,nycph,NBy2
                  y2 = min(y+NBy2-1,nycph)
                  do iz=z,z2
                     do iy=y,y2
		        B(iz,iy,x) = A(iy,x,iz+dnz)
                     enddo
                  enddo
               enddo
               do y=nycph+1,nyc,NBy2
                  y2 = min(y+NBy2-1,ny_fft)
                  do iz=z,z2
                     do iy=y,y2
		        B(iz,iy,x) = A(iy+dny,x,iz+dnz)
                     enddo
                  enddo
               enddo
            enddo
	  enddo

	else

!$OMP parallel do private(x,y,z,y2,z2,iy,iz,C)
           do x=1,iisize
              do z=1,nz_fft,NBz
	         z2 = min(z+NBz-1,nz_fft)

                 do y=1,nycph,NBy2
                    y2 = min(y+NBy2-1,nycph)
                    do iz=z,z2
                        do iy=y,y2
			   C(iz,iy) = A(iy,x,iz)
                        enddo
                     enddo
                  enddo
                 do y=nycph+1,nyc,NBy2
                    y2 = min(y+NBy2-1,nyc)
                    do iz=z,z2
                        do iy=y,y2
			   C(iz,iy) = A(iy+dny,x,iz)
                        enddo
                     enddo
                  enddo
	       enddo

	      if(dnz .gt. 0) then

   	         if(op(3:3) == 't' .or. op(3:3) == 'f') then
                    call exec_f_c2_same(C, 1,nz_fft, &
			  C, 1,nz_fft,nz_fft,nyc)
	         else if(op(3:3) == 'c') then
                   call exec_ctrans_r2_complex_same(C, 2,2*nz_fft, &
				  C, 2,2*nz_fft,nz_fft,nyc)
 	         else if(op(3:3) == 's') then
                   call exec_strans_r2_complex_same(C, 2,2*nz_fft, &
				  C, 2,2*nz_fft,nz_fft,nyc)
                 else
	           print *,taskid,'Unknown transform type: ',op(3:3)
	           call MPI_abort(MPI_COMM_WORLD,ierr)
	         endif
	 	   do y=1,nyc
		      do z=1,nzcph
			B(z,y,x) = C(z,y)
		      enddo
		      do z=nzcph+1,nzc
			B(z,y,x) = C(z+dnz,y)
		      enddo
		   enddo
	      else
   	         if(op(3:3) == 't' .or. op(3:3) == 'f') then
                    call exec_f_c2_dif(C, 1,nz_fft, &
			  B(1,1,x), 1,nz_fft,nz_fft,nyc)
	         else if(op(3:3) == 'c') then
                   call exec_ctrans_r2_complex_dif(C, 2,2*nz_fft, &
				  B(1,1,x), 2,2*nz_fft,nz_fft,nyc)
 	         else if(op(3:3) == 's') then
                   call exec_strans_r2_complex_dif(C, 2,2*nz_fft, &
				  B(1,1,x), 2,2*nz_fft,nz_fft,nyc)
                 else
	           print *,taskid,'Unknown transform type: ',op(3:3)
	           call MPI_abort(MPI_COMM_WORLD,ierr)
	         endif
	      endif
         enddo

      endif

      return
      end subroutine


! This routine is called only when iproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Y

!=============================================================
      subroutine reorder_b2(A,B)
!=============================================================

      implicit none

      complex(p3dfft_type) B(nxhp,ny_fft,kjsize)
      complex(p3dfft_type) A(ny_fft,nxhpc,kjsize)
      integer x,y,z,iy,x2,ix,y2
      complex(p3dfft_type) tmp(nxhpc,ny_fft)

!$OMP parallel do private(x,y,z,y2,x2,iy,ix,tmp)
      do z=1,kjsize
         do x=1,nxhpc,nbx
            x2 = min(x+nbx-1,nxhpc)
            do y=1,ny_fft,nby1
               y2 = min(y+nby1-1,ny_fft)
               do ix=x,x2
                  do iy = y,y2
                     tmp(ix,iy) = A(iy,ix,z)
                  enddo
               enddo
            enddo
         enddo

         do y=1,ny_fft
            do x=1,nxhpc
               B(x,y,z) = tmp(x,y)
            enddo
	    do x=nxhpc+1,nxhp
	       B(x,y,z) = 0.
	    enddo
         enddo
      enddo

      return
      end subroutine

! This routine is called only when iproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Y

!=============================================================
      subroutine reorder_f1(A,B,tmp)
!=============================================================

      implicit none

      complex(p3dfft_type) A(nxhp,ny_fft,kjsize)
      complex(p3dfft_type) B(ny_fft,nxhpc,kjsize)
      integer x,y,z,iy,x2,ix,y2,dnx
      complex(p3dfft_type) tmp(ny_fft,nxhpc)
!      complex(p3dfft_type), allocatable :: tmp(:,:)

!      allocate(tmp(ny_fft,nxhpc))


!$OMP parallel do private(x,y,z,y2,x2,iy,ix,tmp)
      do z=1,kjsize
         do y=1,ny_fft,nby1
            y2 = min(y+nby1-1,ny_fft)
            do x=1,nxhpc,nbx
               x2 = min(x+nbx-1,nxhpc)
               do iy = y,y2
                  do ix=x,x2
                     tmp(iy,ix) = A(ix,iy,z)
                  enddo
               enddo
            enddo
         enddo

         do x=1,nxhpc
            do y=1,ny_fft
               B(y,x,z) = tmp(y,x)
            enddo
         enddo
      enddo

      return
      end subroutine

