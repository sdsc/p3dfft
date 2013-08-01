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

!========================================================
! Transpose  Z-pencils  into Y-pencils

      subroutine bcomm1 (source,dest,t,tc)
!========================================================
      implicit none

      complex(mytype) source(iisize,jjsize,nz_fft)
      complex(mytype) dest(iisize,ny_fft,kjsize)
      real(r8) t,tc
      integer x,y,z,i,ierr,xs,ys,iy,iz,y2,z2,dny
      integer(i8) position,pos1
      
!     Pack the data for sending

#ifdef USE_EVEN

      if(KfCntUneven) then
         tc = tc - MPI_Wtime()
         position = 1
         do i=0,jproc-1
            do z=kjst(i),kjen(i)
               do y=1,jjsize
                  do x=1,iisize
                     buf1(position) = source(x,y,z)
                     position = position+1
                  enddo
               enddo
            enddo
            position = (i+1)*KfCntMax/(mytype*2)+1
         enddo
         tc = tc + MPI_Wtime()
  
         t = t - MPI_Wtime() 
         call mpi_alltoall(buf1,KfCntMax,mpi_byte, &
              buf2,KfCntMax,mpi_byte,mpi_comm_col,ierr)
                  
      else
         t = t - MPI_Wtime() 
         call mpi_alltoall(source,KfCntMax,mpi_byte, &
              buf2,KfCntMax,mpi_byte,mpi_comm_col,ierr)
      endif
#else      
 
!     Exchange data in columns
      t = t - MPI_Wtime() 
      call mpi_alltoallv(source,JrSndCnts, JrSndStrt,mpi_byte, &
           buf2,JrRcvCnts, JrRcvStrt,mpi_byte,mpi_comm_col,ierr)
#endif      


      t = t + MPI_Wtime() 

! Unpack receive buffers into dest

      position=1
      dny = ny_fft - nyc
      do i=0,jproc-1
! If clearly in the first half of ny
         if(jjen(i) .le. nyhc) then
            do z=1,kjsize
               do y=jjst(i),jjen(i)
                  do x=1,iisize
                     dest(x,y,z) = buf2(position)
                     position = position+1
                  enddo
               enddo
	    enddo
! If clearly in the second half of ny
         else if (jjst(i) .ge. nyhc+1) then
            do z=1,kjsize
               do y=jjst(i)+dny,jjen(i)+dny
                  do x=1,iisize
                     dest(x,y,z) = buf2(position)
                     position = position +1
                  enddo
               enddo
            enddo         

! If spanning the first and second half of nz (i.e. jproc is odd)  
         else
	    do z=1,kjsize
               do y=jjst(i),nyhc
                  do x=1,iisize
                     dest(x,y,z) = buf2(position)
                     position = position +1
                  enddo
		enddo
                do y=ny_fft-nyhc+1,jjen(i)+dny
                  do x=1,iisize
                     dest(x,y,z) = buf2(position)
                     position = position +1
                  enddo
               enddo
            enddo
         endif

! Fill center with zeros
         do z=1,kjsize
            do y=nyhc+1,ny_fft-nyhc
               do x=1,iisize
	          dest(x,y,z) = 0.0
               enddo
            enddo
         enddo

#ifdef USE_EVEN
         position = (i+1)*KfCntMax/(mytype*2)+1
#endif
      enddo
      
      
      return
      end subroutine

