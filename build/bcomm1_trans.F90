! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2010 Dmitry Pekurovsky
!    Copyright (C) 2006-2010 University of California
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
! Transpose back Z to Y pencils
! Assumes stride1 data structure

      subroutine bcomm1_trans_many (source,dest,dim,nv,op,t,tc)
!========================================================

      use fft_spec
      implicit none

      integer j,nv,nz,dim
! Assume STRIDE1
!      complex(mytype) source(nzc,jjsize,iisize,nv)
      complex(mytype) source(dim,nv)
      complex(mytype), allocatable :: buf3(:,:)
      complex(mytype) dest(ny_fft,iisize,kjsize,nv)

      real(r8) t,tc
      integer x,y,z,i,ierr,xs,ys,iy,y2,z2,ix,x2,n,iz,dny,dnz
      integer(i8) position,pos1,pos0
      character(len=3) op
      integer sndcnts(0:jproc-1)
      integer rcvcnts(0:jproc-1)
      integer sndstrt(0:jproc-1)
      integer rcvstrt(0:jproc-1)

      nz = nz_fft

!     Pack the data for sending
      
      tc = tc - MPI_Wtime()

      allocate(buf3(nz_fft,jjsize))

     if(jjsize .gt. 0) then
        do j=1,nv
          call pack_bcomm1_trans(buf1,source(1,j),buf3,j,nv,op)     
	enddo
     endif

      tc = tc + MPI_Wtime()
      t = t - MPI_Wtime() 

#ifdef USE_EVEN
      call mpi_alltoall(buf1,KfCntMax*nv, mpi_byte, buf2,KfCntMax*nv,mpi_byte,mpi_comm_col,ierr)
#else
      sndcnts = JrSndCnts * nv
      sndstrt = JrSndStrt * nv
      rcvcnts = JrRcvCnts * nv
      rcvstrt = JrRcvStrt * nv

      call mpi_alltoallv(buf1,SndCnts, SndStrt,mpi_byte, buf2,RcvCnts, RcvStrt,mpi_byte,mpi_comm_col,ierr)
#endif

      t = t + MPI_Wtime() 
      tc = tc - MPI_Wtime()

! Unpack receive buffers into dest

      call unpack_bcomm1_trans_many(dest,buf2,nv)


      deallocate(buf3)

      tc = tc + MPI_Wtime()
      
      return
      end subroutine

      subroutine unpack_bcomm1_trans_many(dest,buf2,nv)

      complex(mytype) dest(ny_fft,iisize,kjsize,nv)
      complex(mytype) buf2(ny_fft*iisize*kjsize*nv)
      integer i,j,nv,position,pos0,pos1,x,y,z,dny      

      dny = ny_fft - nyc
      do i=0,jproc-1
#ifdef USE_EVEN
         pos0 = i*KfCntMax *nv/(mytype*2) +  1 
#else
         pos0 = JrRcvStrt(i) *nv/(mytype*2)+ 1 
#endif
	do j=1,nv
	 do z=1,kjsize
            pos1 = pos0
            do x=1,iisize
               position = pos1
! If clearly in the first half of ny
               if(jjen(i) .le. nyhc) then
                  do y=jjst(i),jjen(i)
                     dest(y,x,z,j) = buf2(position) 
                     position = position+1
                  enddo
! If clearly in the second half of ny
               else if (jjst(i) .ge. nyhc+1) then
                  do y=jjst(i)+dny,jjen(i)+dny
                     dest(y,x,z,j) = buf2(position)
                     position = position +1
                  enddo

! If spanning the first and second half of nz (i.e. jproc is odd)  
              else
                  do y=jjst(i),nyhc
                     dest(y,x,z,j) = buf2(position)
                     position = position +1
   		  enddo
                  do y=ny_fft-nyhc+1,jjen(i)+dny
                     dest(y,x,z,j) = buf2(position)
                     position = position +1
                  enddo
               endif
               pos1 = pos1 + jjsz(i)
            enddo
            pos0 = pos0 + jjsz(i)*iisize
         enddo
	 enddo
      enddo


! Fill center in Y with zeros
      if(dny .ne. 0) then	
        do j=1,nv
         do z=1,kjsize
            do x=1,iisize
               do y=nyhc+1,ny_fft-nyhc
	          dest(y,x,z,j) = 0.0
               enddo
            enddo
         enddo
        enddo
      endif

      return
      end subroutine


      subroutine pack_bcomm1_trans(sendbuf,source,buf3,j,nv,op)

      implicit none

      complex(mytype) source(nzc,jjsize,iisize)
#ifdef USE_EVEN
      complex(mytype) sendbuf(KfCntMax*nv*jproc/(mytrype*2))
#else
      complex(mytype) sendbuf(nzc*jjsize*iisize*nv)
#endif
      complex(mytype) buf3(nz_fft,jjsize)
      integer nz,dnz,i,j,x,y,z,iz,iy,z2,y2,ierr,nv
      integer*8 position,pos0,pos1,pos2
      character(len=3) op

      nz = nz_fft

      dnz = nz - nzc
      if(op(1:1) == '0' .or. op(1:1) == 'n') then

         do x=1,iisize
            do i=0,jproc-1
            
#ifdef USE_EVEN
               pos0 = (i*nv + j-1) * KfCntMax/ (mytype*2) + (x-1)*jjsize 
#else
               pos0 = JrSndStrt(i) *nv/ (mytype*2) + (x-1)*jjsize 
               pos0 = pos0 +(j-1)*kjsz(i)*iisize*jjsize
#endif
	       
	       pos1 = pos0

	       if(kjen(i) .lt. nzhc .or. kjst(i) .gt. nzhc+1) then
! Just copy the data
	           do z=kjst(i),kjen(i),NBz
        	      z2 = min(z+NBz-1,kjen(i))

              		do y=1,jjsize,NBy2
                	   y2 = min(y+NBy2-1,jjsize)
                    
                 	   pos2 = pos1 +y
                     
	                   do iz=z,z2
           	              position = pos2 
                 	      do iy=y,y2
! Here we are sure that dest and buf are different
                      	         sendbuf(position) = source(iz,iy,x)
                       	         position = position + 1
                    	      enddo
                    	   pos2 = pos2 + iisize * jjsize
                        enddo
                     enddo
                     pos1 = pos1 + iisize*jjsize*NBz
	          enddo
	    else 
! Copy some data, then insert zeros to restore full dimension, then again 
! copy some data if needed	
	           do z=kjst(i),nzhc,NBz
        	      z2 = min(z+NBz-1,nzhc)

              		do y=1,jjsize,NBy2
                	   y2 = min(y+NBy2-1,jjsize)
                    
                 	   pos2 = pos1 +y
                     
	                   do iz=z,z2
           	              position = pos2 
                 	      do iy=y,y2
! Here we are sure that dest and buf are different
                      	         sendbuf(position) = source(iz,iy,x)
                       	         position = position + 1
                    	      enddo
                    	   pos2 = pos2 + iisize * jjsize
                        enddo
                     enddo
                     pos1 = pos1 + iisize*jjsize*NBz
	          enddo

                  pos1 = pos0 + (nzhc-kjst(i)+1)*iisize*jjsize
		  do z=1,dnz
		    position = pos1
		    do y=1,jjsize
			sendbuf(position) = 0.
		        position = position +1
		    enddo
		    pos1 = pos1 + iisize*jjsize
		  enddo

	          do z=nzhc+1,kjen(i),NBz
        	      z2 = min(z+NBz-1,kjen(i))

              		do y=1,jjsize,NBy2
                	   y2 = min(y+NBy2-1,jjsize)
                    
                 	   pos2 = pos1 +y
                     
	                   do iz=z,z2
           	              position = pos2 
                 	      do iy=y,y2
! Here we are sure that dest and buf are different
                      	         sendbuf(position) = source(iz,iy,x)
                       	         position = position + 1
                    	      enddo
                    	   pos2 = pos2 + iisize * jjsize
                           enddo
                        enddo
                        pos1 = pos1 + iisize*jjsize*NBz
	           enddo
	        endif
	        pos0 = pos0 + iisize*jjsize*kjsz(i)
	     enddo
         enddo

      else
         do x=1,iisize

	    if(nz .ne. nzc) then

	       do y=1,jjsize
                  do z=1,nzhc
		     buf3(z,y) = source(z,y,x)
	          enddo
	          do z=nzhc+1,nzhc+dnz
		     buf3(z,y) = 0.
	          enddo
	          do z=nzhc+dnz+1,nz_fft
		     buf3(z,y) = source(z-dnz,y,x)
 	          enddo
	       enddo

   	       if(op(1:1) == 't' .or. op(1:1) == 'f') then
                call exec_b_c2_same(buf3, 1,nz_fft, &
				  buf3, 1,nz_fft,nz_fft,jjsize)
 	       else if(op(1:1) == 'c') then	
                   call exec_ctrans_r2_complex_same(buf3, 2,2*nz_fft, & 
				  buf3, 2,2*nz_fft,nz_fft,jjsize)
 	       else if(op(1:1) == 's') then	
                   call exec_strans_r2_complex_same(buf3, 2,2*nz_fft, & 
				  buf3, 2,2*nz_fft,nz_fft,jjsize)
	       else
		   print *,taskid,'Unknown transform type: ',op(1:1)
		   call MPI_abort(MPI_COMM_WORLD,ierr)
	       endif

	    else

     	       if(op(1:1) == 't' .or. op(1:1) == 'f') then
                  call exec_b_c2_dif(source(1,1,x), 1,nz_fft, &
				  buf3, 1,nz_fft,nz_fft,jjsize)
 	       else if(op(1:1) == 'c') then	
                   call exec_ctrans_r2_complex_dif(source(1,1,x), 2,2*nz_fft, & 
				  buf3, 2,2*nz_fft,nz_fft,jjsize)
 	       else if(op(1:1) == 's') then	
                   call exec_strans_r2_complex_dif(source(1,1,x), 2,2*nz_fft, & 
				  buf3, 2,2*nz_fft,nz_fft,jjsize)
	       else
		   print *,taskid,'Unknown transform type: ',op(1:1)
		   call MPI_abort(MPI_COMM_WORLD,ierr)
	       endif

	    endif
         
            do i=0,jproc-1
            

#ifdef USE_EVEN
               pos0 = (i*nv + j-1) * KfCntMax/ (mytype*2) + (x-1)*jjsize 
#else
               pos0 = JrSndStrt(i) *nv/ (mytype*2) + (x-1)*jjsize 
               pos0 = pos0 +(j-1)*kjsz(i)*iisize*jjsize
#endif
	       pos1 = pos0 
               do z=kjst(i),kjen(i),NBz
                  z2 = min(z+NBz-1,kjen(i))

                  do y=1,jjsize,NBy2
                     y2 = min(y+NBy2-1,jjsize)
                     
                     pos2 = pos1 +y
                     
                     do iz=z,z2
                        position = pos2 
                        do iy=y,y2
! Here we are sure that dest and buf are different
                           sendbuf(position) = buf3(iz,iy)
                           position = position + 1
                        enddo
                        pos2 = pos2 + iisize * jjsize
                     enddo
                  enddo
                  pos1 = pos1 + iisize*jjsize*NBz
               enddo
!	       pos0 = pos0 + iisize*jjsize*kjsz(i)
            enddo

         enddo

	endif

	return
	end subroutine

      subroutine bcomm1_trans (source,dest,op,t,tc)
!========================================================

      use fft_spec
      implicit none

! Assume STRIDE1
      complex(mytype) source(nzc,jjsize,iisize)
      complex(mytype), allocatable :: buf3(:,:)
      complex(mytype) dest(ny_fft,iisize,kjsize)

      real(r8) t,tc
      integer x,y,z,i,ierr,xs,ys,iy,y2,z2,ix,x2,n,iz,dny,dnz
      integer(i8) position,pos1,pos0
      character(len=3) op


!     Pack the data for sending

      tc = tc - MPI_Wtime()

      allocate(buf3(nz_fft,jjsize))

      call pack_bcomm1_trans(buf1,source,buf3,1,1,op)     


      tc = tc + MPI_Wtime()
      t = t - MPI_Wtime() 
#ifdef USE_EVEN
      call mpi_alltoall(buf1,KfCntMax, mpi_byte, buf2,KfCntMax,mpi_byte,mpi_comm_col,ierr)
#else
! Use MPI_Alltoallv

      t = t + MPI_Wtime() 
      tc = tc - MPI_Wtime()
         
      call mpi_alltoallv(buf1,JrSndCnts, JrSndStrt,mpi_byte, buf2,JrRcvCnts, JrRcvStrt,mpi_byte,mpi_comm_col,ierr)
#endif

      t = t + MPI_Wtime() 
      tc = tc - MPI_Wtime()
      
      call unpack_bcomm1_trans(dest,buf2)


! Unpack receive buffers into dest

      deallocate(buf3)

      tc = tc + MPI_Wtime()
      
      return
      end subroutine

      subroutine unpack_bcomm1_trans(dest,buf2)

      complex(mytype) dest(ny_fft,iisize,kjsize)
      complex(mytype) buf2(ny_fft*iisize*kjsize)
      integer i,dny,position,pos0,pos1,x,y,z

      dny = ny_fft - nyc
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
! If clearly in the first half of ny
               if(jjen(i) .le. nyhc) then
                  do y=jjst(i),jjen(i)
                     dest(y,x,z) = buf2(position) 
                     position = position+1
                  enddo
! If clearly in the second half of ny
               else if (jjst(i) .ge. nyhc+1) then
                  do y=jjst(i)+dny,jjen(i)+dny
                     dest(y,x,z) = buf2(position)
                     position = position +1
                  enddo

! If spanning the first and second half of nz (i.e. jproc is odd)  
              else
                  do y=jjst(i),nyhc
                     dest(y,x,z) = buf2(position)
                     position = position +1
   		  enddo
                  do y=ny_fft-nyhc+1,jjen(i)+dny
                     dest(y,x,z) = buf2(position)
                     position = position +1
                  enddo
               endif
               pos1 = pos1 + jjsz(i)
            enddo
            pos0 = pos0 + jjsz(i)*iisize
         enddo
      enddo

! Fill center in Y with zeros
      if(dny .ne. 0) then	
         do z=1,kjsize
            do x=1,iisize
               do y=nyhc+1,ny_fft-nyhc
	          dest(y,x,z) = 0.0
               enddo
            enddo
         enddo
      endif

      return 
      end subroutine
