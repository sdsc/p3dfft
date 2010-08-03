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

!========================================================
! Transpose back Z to Y pencils
! Assumes stride1 data structure

      subroutine bcomm1_trans (source,buf3,dest,t,tc)
!========================================================

      use fft_spec
      implicit none

! Assume STRIDE1
      complex(mytype) source(nz_fft,jjsize,iisize)
      complex(mytype) buf3(nz_fft,jjsize)
      complex(mytype) dest(ny_fft,iisize,kjsize)

      real(r8) t,tc
      integer x,y,z,i,ierr,xs,ys,iy,y2,z2,ix,x2,n,iz
      integer(i8) position,pos1,pos0


!     Pack the data for sending

      tc = tc - MPI_Wtime()

#ifdef USE_EVEN
! Use MPI_Alltoall

      if(OW) then

         do x=1,iisize
             call exec_b_c2(source(1,1,x),1,nz_fft,source(1,1,x),1,nz_fft,nz_fft,jjsize)
            
            do i=0,jproc-1
               
               pos0 = i * KfCntMax / (mytype*2) + (x-1)*jjsize 

               do z=kjst(i),kjen(i),NBz
                  z2 = min(z+NBz-1,kjen(i))
                  
                  do y=1,jjsize,NBy2
                     y2 = min(y+NBy2-1,jjsize)
                     
                     pos1 = pos0 +y
                     
                     do iz=z,z2
                        position = pos1 
                        do iy=y,y2
! Here we are sure that dest and buf are different
                           buf1(position) = source(iz,iy,x)
                           position = position + 1
                        enddo
                        pos1 = pos1 + iisize * jjsize
                     enddo
                  enddo
                  pos0 = pos0 + iisize*jjsize*Nbz
               enddo
               
            enddo
         enddo
         
      else

         do x=1,iisize
            call exec_b_c2(source(1,1,x),1,nz_fft,buf3,1,nz_fft,nz_fft,jjsize)
         
            do i=0,jproc-1
            
               pos0 = i * KfCntMax / (mytype*2) + (x-1)*jjsize 

               do z=kjst(i),kjen(i),NBz
                  z2 = min(z+NBz-1,kjen(i))

                  do y=1,jjsize,NBy2
                     y2 = min(y+NBy2-1,jjsize)
                     
                     pos1 = pos0 +y
                     
                     do iz=z,z2
                        position = pos1 
                        do iy=y,y2
! Here we are sure that dest and buf are different
                           buf1(position) = buf3(iz,iy)
                           position = position + 1
                        enddo
                        pos1 = pos1 + iisize * jjsize
                     enddo
                  enddo
               pos0 = pos0 + iisize*jjsize*NBz
               enddo

            enddo
         enddo
      endif

      tc = tc + MPI_Wtime()
      t = t - MPI_Wtime() 
      call mpi_alltoall(buf1,KfCntMax, mpi_byte, buf2,KfCntMax,mpi_byte,mpi_comm_col,ierr)
      t = t + MPI_Wtime() 
      tc = tc - MPI_Wtime()
         
#else
! Use MPI_Alltoallv

      if(OW) then
         do x=1,iisize
            call exec_b_c2(source(1,1,x),1,nz_fft,source(1,1,x),1,nz_fft,nz_fft,jjsize)
            
            pos0 = (x-1)*jjsize
            do z=1,nz_fft,NBz
               z2 = min(z+NBz-1,nz_fft)
               
               do y=1,jjsize,NBy2
                  y2 = min(y+NBy2-1,jjsize)
                  
                  pos1 = pos0 + y
                  
                  do iz=z,z2
                     position = pos1
                     do iy=y,y2
                        buf1(position) = source(iz,iy,x)
                        position = position +1
                     enddo
                     pos1 = pos1 + iisize * jjsize
                  enddo
               enddo
               pos0 = pos0 + iisize*jjsize*NBz
            enddo
            
         enddo

         else

            do x=1,iisize
               call exec_b_c2(source(1,1,x),1,nz_fft,buf3,1,nz_fft,nz_fft,jjsize)
            
               pos0 = (x-1)*jjsize 

               do z=1,nz_fft,NBz
                  z2 = min(z+NBz,nz_fft)
               
                  do y=1,jjsize,NBy2
                     y2 = min(y+NBy2-1,jjsize)
                  
                     pos1 = pos0 + y
                  
                     do iz=z,z2
                        position = pos1
                        do iy=y,y2
                           buf1(position) = buf3(iz,iy)
                           position = position +1
                        enddo
                        pos1 = pos1 + iisize * jjsize
                     enddo
                  enddo
                  pos0 = pos0 + iisize*jjsize*NBz
               enddo
            
            enddo

      endif
         
      tc = tc + MPI_Wtime()
      t = t - MPI_Wtime() 
      call mpi_alltoallv(buf1,JrSndCnts, JrSndStrt,mpi_byte, buf2,JrRcvCnts, JrRcvStrt,mpi_byte,mpi_comm_col,ierr)
      t = t + MPI_Wtime() 
      tc = tc - MPI_Wtime()

#endif

! Unpack receive buffers into dest

         do i=0,jproc-1
#ifdef USE_EVEN
            pos0 = i*KfCntMax/(mytype*2) +  1 
#else
            pos0 = KfSndStrt(i)/(mytype*2)+ 1 
#endif
            do z=1,kjsize
               pos1 = pos0
               do x=1,iisize
                  position = pos1
                  do y=jjst(i),jjen(i)
                     dest(y,x,z) = buf2(position) 
                     position = position+1
                  enddo
                  pos1 = pos1 + jjsz(i)
               enddo
               pos0 = pos0 + jjsz(i)*iisize
            enddo
         enddo
      
      tc = tc + MPI_Wtime()
      
      return
      end subroutine
