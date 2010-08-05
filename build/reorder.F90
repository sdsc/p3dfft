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

! This routine is called only when jproc=1, and only when stride1 is used
! transform backward in Z and transpose array in memory 

      subroutine reorder_trans_b1(A,B,C)

      use fft_spec

      complex(mytype) B(ny_fft,iisize,nz_fft)
      complex(mytype) A(nz_fft,ny_fft,iisize)
      complex(mytype) C(nz_fft,ny_fft)
      integer x,y,z,iy,iz,y2,z2

      if(OW) then

         do x=1,iisize

            call exec_b_c2(A(1,1,x),1,nz_fft,A(1,1,x),1,nz_fft, nz_fft,ny_fft)            
            do y=1,ny_fft,NBy2
               y2 = min(y+NBy2-1,ny_fft)
   	       do z=1,nz_fft,NBz
	          z2 = min(z+NBz-1,nz_fft)
                    do iy=y,y2
                        do iz=z,z2
			   B(iy,x,iz) = A(iz,iy,x)
                        enddo
                     enddo
                  enddo
             enddo 

          enddo

	else

         do x=1,iisize

            call exec_f_c2(A(1,1,x),1,nz_fft,C,1,nz_fft, nz_fft,ny_fft)            
            do y=1,ny_fft,NBy2
               y2 = min(y+NBy2-1,ny_fft)
   	       do z=1,nz_fft,NBz
	          z2 = min(z+NBz-1,nz_fft)
                    do iy=y,y2
                        do iz=z,z2
			   B(iy,x,iz) = C(iz,iy)
                        enddo
                     enddo
                enddo
              enddo 
          enddo
      endif

      return
      end subroutine

! This routine is called only when iproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Y

      subroutine reorder_b2(A,B)

      implicit none

      complex(mytype) B(nxhp,ny_fft,kjsize)
      complex(mytype) A(ny_fft,nxhp,kjsize)
      integer x,y,z,iy,x2,ix,y2
      complex(mytype), allocatable :: tmp(:,:)


      allocate(tmp(nxhp,0:ny_fft-1))

      do z=1,kjsize
         do x=1,nxhp,nbx
            x2 = min(x+nbx-1,nxhp)
            do y=0,ny_fft-1,nby1
               y2 = min(y+nby1-1,ny_fft-1)
               do ix=x,x2
                  do iy = y,y2
                     tmp(ix,iy) = A(iy+1,ix,z)
                  enddo
               enddo
            enddo
         enddo

         do y=0,ny_fft-1
            do x=1,nxhp
               B(x,y+1,z) = tmp(x,y)
            enddo
         enddo
      enddo

      deallocate(tmp)

      return
      end subroutine

! This routine is called only when iproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Y

      subroutine reorder_f1(A,B)

      implicit none

      complex(mytype) A(nxhp,ny_fft,kjsize)
      complex(mytype) B(ny_fft,nxhp,kjsize)
      integer x,y,z,iy,x2,ix,y2
      complex(mytype), allocatable :: tmp(:,:)

      allocate(tmp(0:ny_fft-1,nxhp))


      do z=1,kjsize
         do y=0,ny_fft-1,nby1
            y2 = min(y+nby1-1,ny_fft-1)
            do x=1,nxhp,nbx
               x2 = min(x+nbx-1,nxhp)
               do iy = y,y2
                  do ix=x,x2
                     tmp(iy,ix) = A(ix,iy+1,z)
                  enddo
               enddo
            enddo
         enddo

         do x=1,nxhp
            do y=0,ny_fft-1
               B(y+1,x,z) = tmp(y,x)
            enddo
         enddo
      enddo

      deallocate(tmp)

      return
      end subroutine

! This routine is called only when jproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Z

      subroutine reorder_trans_f2(A,B)

      use fft_spec
      implicit none

      complex(mytype) A(ny_fft,iisize,nz_fft)
      complex(mytype) B(nz_fft,ny_fft,iisize)
      integer x,y,z,iy,iz,y2,z2


      do x=1,iisize
         do z=1,nz_fft,NBz
	    z2 = min(z+NBz-1,nz_fft)
            
            do y=1,ny_fft,NBy2
               y2 = min(y+NBy2-1,ny_fft)

                  do iz=z,z2
                     do iy=y,y2
			B(iz,iy,x) = A(iy,x,iz)
                     enddo
                  enddo
               enddo
           enddo 

         call exec_f_c2(B(1,1,x),1,nz_fft,B(1,1,x),1,nz_fft, nz_fft,ny_fft)
      enddo

      return
      end subroutine
