! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2013 Dmitry Pekurovsky
!    Copyright (C) 2006-2013 University of California
!    Copyright (C) 2010-2011 Jens Henrik Goebbert
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
! Transpose an array in Y=pencils into Z-pencils
! Uses MPI_Alltoall(v)
!
      subroutine fcomm2(source,dest,t,tc)
!========================================================

      implicit none

      complex(mytype) source(iisize,ny_fft,kjsize)
      complex(mytype) dest(iisize,jjsize,nz_fft)

      real(r8) t,tc
      integer x,z,y,i,ierr,xs,ys,y2,z2,iy,iz,dny
      integer(i8) position,pos1,pos0

!	if(taskid .eq. 0) then
!	  print *,'Entring fcomm2'
!        endif	
!	call print_buf(source,iisize,ny_fft,kjsize)

      dny = ny_fft-nyc

! Pack send buffers for exchanging y and z for all x at once 

      position = 1
      do i=0,jproc-1
#ifdef USE_EVEN
         pos0 = i*KfCntMax/(mytype*2)  + 1 
#else
         pos0 = KfSndStrt(i)/(mytype*2)+ 1 
#endif

! If clearly in the first half of ny

         if(jjen(i) .le. nyhc) then
     	    do z=1,kjsize
               position = pos0 +(z-1)*jjsz(i)*iisize
               do y=jjst(i),jjen(i)
   		  do x=1,iisize
                     buf1(position) = source(x,y,z)
                     position = position+1
                  enddo
               enddo	
            enddo

! If clearly in the second half of ny
         else if (jjst(i) .ge. nyhc+1) then
     	    do z=1,kjsize
               position = pos0 +(z-1)*jjsz(i)*iisize
               do y=jjst(i)+dny,jjen(i)+dny
                  do x=1,iisize
                     buf1(position) = source(x,y,z)
                     position = position+1
                  enddo
               enddo	
            enddo



! If spanning the first and second half of ny (e.g. iproc is odd)
         else
     	    do z=1,kjsize
               position = pos0 +(z-1)*jjsz(i)*iisize
               do y=jjst(i),nyhc
                  do x=1,iisize
                     buf1(position) = source(x,y,z)
                     position = position+1
                  enddo
	       enddo	
               do y=ny_fft-nyhc+1,jjen(i)+dny
                  do x=1,iisize
                     buf1(position) = source(x,y,z)
                     position = position+1
                  enddo
               enddo	
            enddo
         endif

      enddo
      
! Exchange y-z buffers in columns of processors

      t = t - MPI_Wtime()

#ifdef USE_EVEN
! Use MPI_Alltoall

      if(KfCntUneven) then

         call mpi_alltoall(buf1,KfCntMax, mpi_byte, &
           buf2,KfCntMax, mpi_byte,mpi_comm_col,ierr)

         t = MPI_Wtime() + t

         tc = tc - MPI_Wtime()

         position = 1
         do i=0,jproc-1
            do z=kjst(i),kjen(i)
               do y=1,jjsize
                  do x=1,iisize
                     dest(x,y,z) = buf2(position)
                     position = position +1
                  enddo
               enddo
            enddo
            position = (i+1)*KfCntMax/(mytype*2)+1
         enddo

         tc = tc + MPI_Wtime()

      else

         call mpi_alltoall(buf1,KfCntMax, mpi_byte, &
           dest,KfCntMax, mpi_byte,mpi_comm_col,ierr)
         t = MPI_Wtime() + t

      endif

#else
! Use MPI_Alltoallv

      call mpi_alltoallv(buf1,KfSndCnts, KfSndStrt,mpi_byte, &
           dest,KfRcvCnts, KfRcvStrt,mpi_byte,mpi_comm_col,ierr)
      t = MPI_Wtime() + t
         
#endif
      return
      end subroutine
