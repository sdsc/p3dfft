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

      subroutine bcomm1_trans (source,dest,buf3,op,t,tc)
!========================================================

      use fft_spec
      implicit none

! Assume STRIDE1
      complex(mytype) source(nzc,jjsize,iisize)
      complex(mytype) buf3(nz_fft,jjsize)
      complex(mytype) dest(ny_fft,iisize,kjsize)

      real(r8) t,tc
      integer x,y,z,i,ierr,xs,ys,iy,y2,z2,ix,x2,n,iz,dny,dnz
      integer(i8) position,pos1,pos0
      character(len=3) op


!     Pack the data for sending


      tc = tc - MPI_Wtime()

#ifdef USE_EVEN
! Use MPI_Alltoall

     if(jjsize .gt. 0) then

      dnz = nz - nzc
      if(op(1:1) == '0' .or. op(1:1) == 'n') then
         do x=1,iisize
            do i=0,jproc-1
            
               pos0 = i * KfCntMax / (mytype*2) + (x-1)*jjsize 
	       if(kjen(i) .lt. nzhc .or. kjst(i) .gt. nzhc+1) then
! Just copy the data
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
                     pos0 = pos0 + iisize*jjsize*NBz
	          enddo
	    else then
! Copy some data, then insert zeros to restore full dimension, then again 
! copy some data if needed	
	           do z=kjst(i),nzhc,NBz
        	      z2 = min(z+NBz-1,nzhc)

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
                     pos0 = pos0 + iisize*jjsize*NBz
	          enddo

                  pos0 = i * KfCntMax / (mytype*2) + (x-1)*jjsize  &
		     +(nzhc-kjst(i)+1)*iisize*jjsize
		  do z=1,dnz
		    position = pos0
		    do y=1,jjsize
			buf1(position) = 0.
		        position = position +1
		    enddo
		    pos0 = pos0 + iisize*jjsize
		  enddo

	           do z=nzhc+1,kjen(i),NBz
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
                     pos0 = pos0 + iisize*jjsize*NBz
	          enddo
	     endif
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
     endif

      tc = tc + MPI_Wtime()
      t = t - MPI_Wtime() 
      call mpi_alltoall(buf1,KfCntMax, mpi_byte, buf2,KfCntMax,mpi_byte,mpi_comm_col,ierr)
      t = t + MPI_Wtime() 
      tc = tc - MPI_Wtime()
         
#else
! Use MPI_Alltoallv

     if(jjsize .gt. 0) then

      dnz = nz_fft - nzc
      if(op(1:1) == '0' .or. op(1:1) == 'n') then

	 do x=1,iisize
           pos0 = (x-1)*jjsize 

           do z=1,nzhc,NBz
              z2 = min(z+NBz-1,nzhc)

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
              pos0 = pos0 + iisize*jjsize*NBz
	   enddo

           pos1 = (x-1)*jjsize + iisize*jjsize*nzhc +1
           do z=nzhc+1,nz_fft-nzhc
	      position = pos1
   	      do y=1,jjsize
		 buf1(position) = 0.
                 position = position + 1
              enddo
              pos1 = pos1 + iisize * jjsize
	   enddo

	   pos0 = pos1 -1
           do z=nz_fft-nzhc+1,nz_fft,NBz
              z2 = min(z+NBz-1,nz_fft)

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
              pos0 = pos0 + iisize*jjsize*NBz
	   enddo

        enddo

      else
         do x=1,iisize
           pos0 = (x-1)*jjsize 

	    do y=1,jjsize
               do z=1,nzhc
		  buf3(z,y) = source(z,y,x)
	       enddo
	       do z=nzhc+1,nz_fft-nzhc
		  buf3(z,y) = 0.
	       enddo
	       do z=nz_fft-nzhc+1,nz_fft
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
            
            do z=1,nz_fft,NBz
              z2 = min(z+NBz-1,nz_fft)

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
	endif
     endif
         
      tc = tc + MPI_Wtime()
      t = t - MPI_Wtime() 
      call mpi_alltoallv(buf1,JrSndCnts, JrSndStrt,mpi_byte, buf2,JrRcvCnts, JrRcvStrt,mpi_byte,mpi_comm_col,ierr)
      t = t + MPI_Wtime() 
      tc = tc - MPI_Wtime()

#endif

! Unpack receive buffers into dest

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

      tc = tc + MPI_Wtime()
      
      return
      end subroutine
