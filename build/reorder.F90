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
      complex(mytype) A(dim,nv)
      complex(mytype) B(ny_fft,iisize,nz_fft,nv)
      complex(mytype) C(nz_fft,nyc)

      do j=1,nv
         call reorder_trans_b1(A(1,j),B(1,1,1,j),C,op)
      enddo

      return
      end subroutine

!=============================================================
      subroutine reorder_trans_b1(A,B,C,op)
!=============================================================

      use fft_spec

      complex(mytype) A(nzc,nyc,iisize)
      complex(mytype) B(ny_fft,iisize,nz_fft)
      complex(mytype) C(nz_fft,nyc)
      integer x,y,z,iy,iz,y2,z2,ierr,dnz,dny
      character(len=3) op

      dny = ny_fft - nyc
      dnz = nz_fft - nzc
      if(op(1:1) == '0' .or. op(1:1) == 'n') then

         do x=1,iisize
            do y=1,nyhc,NBy2
               y2 = min(y+NBy2-1,nyhc)

   	       do z=1,nzhc,NBz
	          z2 = min(z+NBz-1,nzhc)
                    do iy=y,y2
                        do iz=z,z2
			   B(iy,x,iz) = A(iz,iy,x)
                        enddo
                     enddo
                enddo

   	       do z=nzhc+dnz+1,nz_fft,NBz
	          z2 = min(z+NBz-1,nz_fft)
                    do iy=y,y2
                        do iz=z,z2
			   B(iy,x,iz) = A(iz-dnz,iy,x)
                        enddo
                     enddo
                enddo

              enddo 

            do y=nyhc+1,nyc,NBy2
               y2 = min(y+NBy2-1,nyc)

   	       do z=1,nzhc,NBz
	          z2 = min(z+NBz-1,nzhc)
                    do iy=y,y2
                        do iz=z,z2
			   B(iy+dny,x,iz) = A(iz,iy,x)
                        enddo
                     enddo
                enddo

   	       do z=nzhc+dnz+1,nz_fft,NBz
	          z2 = min(z+NBz-1,nz_fft)
                    do iy=y,y2
                        do iz=z,z2
			   B(iy+dny,x,iz) = A(iz-dnz,iy,x)
                        enddo
                     enddo
                enddo

              enddo 

          enddo

	  do z=nzhc+1,nzhc+dnz
	     do x=1,iisize
                do y=1,ny_fft
	 	   B(y,x,z) = 0
                enddo
             enddo
          enddo

	else

           do x=1,iisize
	      do y=1,nyc
	         do z=1,nzhc
		    C(z,y) = A(z,y,x)
		 enddo
	         do z=nzhc+1,nzhc+dnz
		    C(z,y) = 0.
		 enddo
	         do z=nzhc+dnz+1,nz_fft
		    C(z,y) = A(z-dnz,y,x)
		 enddo
	      enddo 

	      if(op(1:1) == 't' .or. op(1:1) == 'f') then
                 call exec_b_c2_same(C, 1,nz_fft, &
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

              do y=1,nyhc,NBy2
                 y2 = min(y+NBy2-1,nyhc)
     	         do z=1,nz_fft,NBz
	            z2 = min(z+NBz-1,nz_fft)
                       do iy=y,y2
                          do iz=z,z2
			     B(iy,x,iz) = C(iz,iy)
                          enddo
                       enddo
                  enddo
              enddo 
              do y=nyhc+1,nyc,NBy2
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
           do y=nyhc+1,nyhc+dny
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

      complex(mytype) B(nxhp,ny_fft,kjsize,nv)
      complex(mytype) A(ny_fft,nxhpc,kjsize,nv)
      integer x,y,z,iy,x2,ix,y2,nv,j
      complex(mytype), allocatable :: tmp(:,:)


      allocate(tmp(nxhpc,ny_fft))

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

      deallocate(tmp)

      return
      end subroutine

! This routine is called only when iproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Y

!=============================================================
      subroutine reorder_f1_many(A,B,tmp,nv)
!=============================================================

      implicit none

      complex(mytype) A(nxhp,ny_fft,kjsize,nv)
      complex(mytype) B(ny_fft,nxhpc,kjsize,nv)
      integer x,y,z,iy,x2,ix,y2,dnx,nv,j
      complex(mytype) tmp(ny_fft,nxhpc,kjsize)
!      complex(mytype), allocatable :: tmp(:,:)

!      allocate(tmp(ny_fft,nxhpc))

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
      complex(mytype) B(dim,nv)
      complex(mytype) A(ny_fft,iisize,nz_fft,nv)
      complex(mytype) C(nz_fft,ny_fft)
      character(len=3) op

      do j=1,nv
         call reorder_trans_f2(A(1,1,1,j),B(1,j),C,op)
      enddo

      return
      end subroutine

      subroutine reorder_trans_f2(A,B,C,op)

      use fft_spec
      implicit none

      integer x,y,z,iy,iz,y2,z2,ierr,dnz,dny
      complex(mytype) A(ny_fft,iisize,nz_fft)
      complex(mytype) B(nzc,nyc,iisize)
      complex(mytype) C(nz_fft,ny_fft)
      character(len=3) op

      dnz = nz_fft - nzc
      dny = ny_fft - nyc
      if(op(3:3) == '0' .or. op(3:3) == 'n') then
	
         do x=1,iisize
            do z=1,nzhc,NBz
	       z2 = min(z+NBz-1,nzhc)
            
               do y=1,nyhc,NBy2
                  y2 = min(y+NBy2-1,nyhc)
                  do iz=z,z2
                     do iy=y,y2
  		        B(iz,iy,x) = A(iy,x,iz)
                     enddo
                  enddo
               enddo

               do y=nyhc+dny+1,ny_fft,NBy2
                  y2 = min(y+NBy2-1,ny_fft)
                  do iz=z,z2
                     do iy=y,y2
  		        B(iz,iy,x) = A(iy-dny,x,iz)
                     enddo
                  enddo
               enddo
            enddo
            do z=nzhc+dnz+1,nz_fft,NBz
	       z2 = min(z+NBz-1,nz_fft)
            
               do y=1,nyhc,NBy2
                  y2 = min(y+NBy2-1,nyhc)
                  do iz=z,z2
                     do iy=y,y2
		        B(iz,iy,x) = A(iy,x,iz-dnz)
                     enddo
                  enddo
               enddo
               do y=nyhc+dny+1,ny_fft,NBy2
                  y2 = min(y+NBy2-1,ny_fft)
                  do iz=z,z2
                     do iy=y,y2
		        B(iz,iy,x) = A(iy-dny,x,iz-dnz)
                     enddo
                  enddo
               enddo
            enddo 
	  enddo

	else

           do x=1,iisize
              do z=1,nz_fft,NBz
	         z2 = min(z+NBz-1,nz_fft)
            
                 do y=1,nyhc,NBy2
                    y2 = min(y+NBy2-1,nyhc)
                    do iz=z,z2
                        do iy=y,y2
			   C(iz,iy) = A(iy,x,iz)
                        enddo
                     enddo
                  enddo
                 do y=nyhc+dny+1,ny_fft,NBy2
                    y2 = min(y+NBy2-1,ny_fft)
                    do iz=z,z2
                        do iy=y,y2
			   C(iz,iy) = A(iy-dny,x,iz)
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
		      do z=1,nzhc
			B(z,y,x) = C(z,y)
		      enddo	
		      do z=nzhc+1,nzc
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

      complex(mytype) B(nxhp,ny_fft,kjsize)
      complex(mytype) A(ny_fft,nxhpc,kjsize)
      integer x,y,z,iy,x2,ix,y2
      complex(mytype), allocatable :: tmp(:,:)


      allocate(tmp(nxhpc,ny_fft))

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


      deallocate(tmp)

      return
      end subroutine

! This routine is called only when iproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Y

!=============================================================
      subroutine reorder_f1(A,B,tmp)
!=============================================================

      implicit none

      complex(mytype) A(nxhp,ny_fft,kjsize)
      complex(mytype) B(ny_fft,nxhpc,kjsize)
      integer x,y,z,iy,x2,ix,y2,dnx
      complex(mytype) tmp(ny_fft,nxhpc,kjsize)
!      complex(mytype), allocatable :: tmp(:,:)

!      allocate(tmp(ny_fft,nxhpc))


      do z=1,kjsize
         do y=1,ny_fft,nby1
            y2 = min(y+nby1-1,ny_fft)
            do x=1,nxhpc,nbx
               x2 = min(x+nbx-1,nxhpc)
               do iy = y,y2
                  do ix=x,x2
                     tmp(iy,ix,z) = A(ix,iy,z)
                  enddo
               enddo
            enddo
         enddo

         do x=1,nxhpc
            do y=1,ny_fft
               B(y,x,z) = tmp(y,x,z)
            enddo
         enddo
      enddo

      return
      end subroutine

