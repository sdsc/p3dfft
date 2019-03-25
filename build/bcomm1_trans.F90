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
!      complex(p3dfft_type) source(nzc,jjsize,iisize,nv)
      complex(p3dfft_type) source(dim,nv)
      complex(p3dfft_type), allocatable :: buf3(:,:)
      complex(p3dfft_type) dest(ny_fft,iisize,kjsize,nv)

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

      complex(p3dfft_type) dest(ny_fft,iisize,kjsize,nv)
      complex(p3dfft_type) buf2(ny_fft*iisize*kjsize*nv)
      integer i,j,nv,position,pos0,pos1,x,y,z,dny

      dny = ny_fft - nyc
!$OMP PARALLEL DO private(i,j,pos0,pos1,position,x,y,z)
      do i=0,jproc-1
#ifdef USE_EVEN
         pos0 = i*KfCntMax *nv/(p3dfft_type*2) +  1
#else
         pos0 = JrRcvStrt(i) *nv/(p3dfft_type*2)+ 1
#endif
	do j=1,nv
	 do z=1,kjsize
            pos1 = pos0
            do x=1,iisize
               position = pos1

! First set of significant Fourier modes (zero and positive, total a little more than half)
               do y=jjst(i),min(jjen(i),nycph)
                  dest(y,x,z,j) = buf2(position)
                  position = position+1
               enddo
! Second set of Fourier modes (negative)
               do y=max(jjst(i)+dny,ny_fft-nyhc+1),jjen(i)+dny
                  dest(y,x,z,j) = buf2(position)	
                  position = position +1
               enddo

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
               do y=nycph+1,ny_fft-nyhc
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

      complex(p3dfft_type) source(nzc,jjsize,iisize)
#ifdef USE_EVEN
      complex(p3dfft_type) sendbuf(KfCntMax*nv*jproc/(p3dfft_type*2))
#else
      complex(p3dfft_type) sendbuf(nz_fft*jjsize*iisize*nv)
#endif
      complex(p3dfft_type) buf3(nz_fft,jjsize)
      integer nz,dnz,i,j,x,y,z,iz,iy,z2,y2,ierr,nv
      integer*8 position,pos0,pos1,pos2
      character(len=3) op

      nz = nz_fft

      dnz = nz - nzc
      if(op(1:1) == '0' .or. op(1:1) == 'n') then

!$OMP PARALLEL DO private(i,pos0,pos1,pos2,position,x,y,z,iy,y2,iz,z2) collapse(2)
         do x=1,iisize
            do i=0,jproc-1

#ifdef USE_EVEN
               pos0 = i*nv * KfCntMax/ (p3dfft_type*2) + (x-1)*jjsize
#else
               pos0 = JrSndStrt(i) *nv/ (p3dfft_type*2) + (x-1)*jjsize
#endif
               pos0 = pos0 +(j-1)*kjsz(i)*iisize*jjsize

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
!$OMP PARALLEL DO private(i,pos0,pos1,pos2,position,x,y,z,iy,y2,iz,z2,buf3)
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
               pos0 = i*nv * KfCntMax/ (p3dfft_type*2) + (x-1)*jjsize
#else
               pos0 = JrSndStrt(i) *nv/ (p3dfft_type*2) + (x-1)*jjsize
#endif
               pos0 = pos0 +(j-1)*kjsz(i)*iisize*jjsize
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
      complex(p3dfft_type) source(nzc,jjsize,iisize)
      complex(p3dfft_type), allocatable :: buf3(:,:)
      complex(p3dfft_type) dest(ny_fft,iisize,kjsize)

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

      complex(p3dfft_type) dest(ny_fft,iisize,kjsize)
      complex(p3dfft_type) buf2(ny_fft*iisize*kjsize)
      integer i,dny,position,pos0,pos1,x,y,z

      dny = ny_fft - nyc
!$OMP PARALLEL DO private(i,pos0,pos1,position,x,y,z)
      do i=0,jproc-1
#ifdef USE_EVEN
         pos0 = i*KfCntMax/(p3dfft_type*2) +  1
#else
         pos0 = KfSndStrt(i)/(p3dfft_type*2)+ 1
#endif
	 do z=1,kjsize
            pos1 = pos0
            do x=1,iisize
               position = pos1

! First set of significant Fourier modes (zero and positive, total a little more than half)
               do y=jjst(i),min(jjen(i),nycph)
                  dest(y,x,z) = buf2(position)
                  position = position+1
               enddo
! Second set of Fourier modes (negative)
               do y=max(jjst(i)+dny,ny_fft-nyhc+1),jjen(i)+dny
                  dest(y,x,z) = buf2(position)	
                  position = position +1
               enddo

               pos1 = pos1 + jjsz(i)
            enddo
            pos0 = pos0 + jjsz(i)*iisize
         enddo
      enddo

! Fill center in Y with zeros
      if(dny .ne. 0) then
         do z=1,kjsize
            do x=1,iisize
               do y=nycph+1,ny_fft-nyhc
	          dest(y,x,z) = 0.0
               enddo
            enddo
         enddo
      endif

      return
      end subroutine
