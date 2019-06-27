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
! Transpose  Z-pencils  into Y-pencils

      subroutine bcomm1_many (source,dest,dim,nv,t,tc)
!========================================================
      implicit none

      integer x,y,z,i,ierr,xs,ys,iy,iz,y2,z2,dny,j,nv,dim
!      complex(p3dfft_type) source(iisize,jjsize,nz_fft,nv)
      complex(p3dfft_type) source(dim,nv)
      complex(p3dfft_type) dest(iisize,ny_fft,kjsize,nv)
      real(r8) t,tc
      integer(i8) position,pos1
      integer sndcnts(0:jproc-1)
      integer rcvcnts(0:jproc-1)
      integer sndstrt(0:jproc-1)
      integer rcvstrt(0:jproc-1)


!     Pack the data for sending

      tc = tc - MPI_Wtime()
      do j=1,nv
	 call pack_bcomm1(source(1,j),j,nv)
      enddo
      tc = tc + MPI_Wtime()

!     Exchange data in columns
      t = t - MPI_Wtime()

#ifdef USE_EVEN
      call mpi_alltoall(buf1,KfCntMax *nv,mpi_byte, &
              buf2,KfCntMax *nv,mpi_byte,mpi_comm_col,ierr)

#else
      sndcnts = JrSndCnts * nv
      sndstrt = JrSndStrt * nv
      rcvcnts = JrRcvCnts * nv
      rcvstrt = JrRcvStrt * nv

      call mpi_alltoallv(buf1,SndCnts, SndStrt,mpi_byte, &
           buf2,RcvCnts, RcvStrt,mpi_byte,mpi_comm_col,ierr)
#endif

      t = t + MPI_Wtime()

! Unpack receive buffers into dest

      call unpack_bcomm1_many(dest,buf2,nv)

      return
      end subroutine

!========================================================
      subroutine unpack_bcomm1_many(dest,buf2,nv)
!========================================================

      complex(p3dfft_type) dest(iisize,ny_fft,kjsize,nv)
      complex(p3dfft_type) buf2(iisize*ny_fft*kjsize*nv)
      integer i,j,position,pos0,pos1,x,y,z,dny,nv

      position=1
      dny = ny_fft - nyc
!$OMP PARALLEL DO private(i,j,pos0,position,x,y,z) collapse(2)
      do j=1,nv
      do i=0,jproc-1
#ifdef USE_EVEN
         pos0 = i * nv * KfCntMax/(p3dfft_type*2)
#else
         pos0 = nv * KfSndStrt(i)/(p3dfft_type*2)
#endif
 	pos0 = pos0 + (j-1)*KfSndCnts(i)/(p3dfft_type*2)+ 1

! If clearly in the first half of ny
         if(jjen(i) .le. nyhc) then
            do z=1,kjsize
               position = pos0 +(z-1)*jjsz(i)*iisize
               do y=jjst(i),jjen(i)
                  do x=1,iisize
                     dest(x,y,z,j) = buf2(position)
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
                     dest(x,y,z,j) = buf2(position)
                     position = position +1
                  enddo
               enddo
            enddo

! If spanning the first and second half of nz (i.e. jproc is odd)
         else
	    do z=1,kjsize
               position = pos0 +(z-1)*jjsz(i)*iisize
               do y=jjst(i),nyhc
                  do x=1,iisize
                     dest(x,y,z,j) = buf2(position)
                     position = position +1
                  enddo
		enddo
                do y=ny_fft-nyhc+1,jjen(i)+dny
                  do x=1,iisize
                     dest(x,y,z,j) = buf2(position)
                     position = position +1
                  enddo
               enddo
            enddo
         endif

! Fill center with zeros
         do z=1,kjsize
            do y=nyhc+1,ny_fft-nyhc
               do x=1,iisize
	          dest(x,y,z,j) = 0.0
               enddo
            enddo
         enddo
         enddo
      enddo

      end subroutine


!========================================================
      subroutine pack_bcomm1(A,j,nv)
!========================================================

      integer j,nv
      complex(p3dfft_type) A(iisize,jjsize,nz_fft)
      integer x,y,z,i,ierr,xs,ys,iy,iz,y2,z2,dny
      integer(i8) position,pos1
      integer sndcnts(0:jproc-1)
      integer rcvcnts(0:jproc-1)
      integer sndstrt(0:jproc-1)
      integer rcvstrt(0:jproc-1)

!$OMP PARALLEL DO private(i,position,x,y,z)
      do i=0,jproc-1
#ifdef USE_EVEN
         position = i*KfCntMax*nv/(p3dfft_type*2)+1
#else
         position = JrSndStrt(i)*nv/(p3dfft_type*2)+1
#endif
         position = position + (j-1)*iisize*jjsize*kjsz(i)

         do z=kjst(i),kjen(i)
            do y=1,jjsize
               do x=1,iisize
                  buf1(position) = A(x,y,z)
                  position = position+1
               enddo
            enddo
         enddo
      enddo

      return
      end subroutine

!========================================================
      subroutine bcomm1 (source,dest,t,tc)
!========================================================
      implicit none

      complex(p3dfft_type) source(iisize,jjsize,nz_fft)
      complex(p3dfft_type) dest(iisize,ny_fft,kjsize)
      real(r8) t,tc
      integer x,y,z,i,j,ierr,xs,ys,iy,iz,y2,z2,dny
      integer(i8) position,pos1,pos0

!     Pack the data for sending

#ifdef USE_EVEN

      if(KfCntUneven) then
         tc = tc - MPI_Wtime()
         position = 1
!$OMP PARALLEL DO private(i,j,pos0,position,x,y,z)
         do i=0,jproc-1
            position = i*KfCntMax/(p3dfft_type*2)  + 1
            do z=kjst(i),kjen(i)
               do y=1,jjsize
                  do x=1,iisize
                     buf1(position) = source(x,y,z)
                     position = position+1
                  enddo
               enddo
            enddo
            position = (i+1)*KfCntMax/(p3dfft_type*2)+1
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

      call unpack_bcomm1(dest,buf2)

      return
      end subroutine

      subroutine unpack_bcomm1(dest,buf2)

      complex(p3dfft_type) dest(iisize,ny_fft,kjsize)
      complex(p3dfft_type) buf2(iisize*ny_fft*kjsize)
      integer i,position,pos0,x,y,z,dny

      position=1
      dny = ny_fft - nyc
!$OMP PARALLEL DO private(i,pos0,position,x,y,z)
      do i=0,jproc-1
#ifdef USE_EVEN
         pos0 = i*KfCntMax/(p3dfft_type*2)  + 1
#else
         pos0 = KfSndStrt(i)/(p3dfft_type*2)+ 1
#endif

! If clearly in the first half of ny
         if(jjen(i) .le. nyhc) then
            do z=1,kjsize
               position = pos0 +(z-1)*jjsz(i)*iisize
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
               position = pos0 +(z-1)*jjsz(i)*iisize
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
               position = pos0 +(z-1)*jjsz(i)*iisize
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

      enddo


      return
      end subroutine

