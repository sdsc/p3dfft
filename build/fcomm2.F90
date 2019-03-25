!! This file is part of P3DFFT library

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
!----------------------------------------------------------------------------

!========================================================
! Transpose an array in Y=pencils into Z-pencils
! Uses MPI_Alltoall(v)
!
      subroutine fcomm2_many(source,dest,dim_out,nv,t,tc)
!========================================================

      implicit none

      real(r8) t,tc
      integer x,z,y,i,ierr,xs,ys,y2,z2,iy,iz,dny,nv,j,dim_out
      integer(i8) position,pos1,pos0
      integer sndcnts(0:jproc-1)
      integer rcvcnts(0:jproc-1)
      integer sndstrt(0:jproc-1)
      integer rcvstrt(0:jproc-1)
      complex(p3dfft_type) source(iisize,ny_fft,kjsize,nv)
      complex(p3dfft_type) dest(dim_out,nv)


! Pack send buffers for exchanging y and z for all x at once

      call pack_fcomm2_many(buf1,source,nv)

! Exchange y-z buffers in columns of processors

      t = t - MPI_Wtime()

#ifdef USE_EVEN
! Use MPI_Alltoall

         call mpi_alltoall(buf1,KfCntMax * nv, mpi_byte, &
           buf2,KfCntMax * nv, mpi_byte,mpi_comm_col,ierr)

#else
! Use MPI_Alltoallv

      sndcnts = KfSndCnts * nv
      sndstrt = KfSndStrt * nv
      rcvcnts = KfRcvCnts * nv
      rcvstrt = KfRcvStrt * nv

      call mpi_alltoallv(buf1,SndCnts, SndStrt,mpi_byte, &
           buf2,RcvCnts, RcvStrt,mpi_byte,mpi_comm_col,ierr)
#endif

      t = MPI_Wtime() + t

         tc = tc - MPI_Wtime()
	 do j=1,nv
	    call unpack_fcomm2(dest(1,j),j,nv)
         enddo

         tc = tc + MPI_Wtime()



      return
      end subroutine

      subroutine pack_fcomm2_many(sndbuf,source,nv)

      use fft_spec
      implicit none

      complex(p3dfft_type) source(iisize,ny_fft,kjsize,nv)
      complex(p3dfft_type) sndbuf(iisize*ny_fft*kjsize*nv)
      integer nv,j,i,position,pos0,pos1,x,y,z,dny

      dny = ny_fft-nyc
      position = 1

!$OMP PARALLEL DO private(i,j,pos0,position,x,y,z) collapse(2)
      do j=1,nv

      do i=0,jproc-1
#ifdef USE_EVEN
        pos0 = i * nv * KfCntMax/(p3dfft_type*2)
#else
        pos0 = nv * KfSndStrt(i) /(p3dfft_type*2)
#endif
 	pos0 = pos0 + (j-1)*KfSndCnts(i)/(p3dfft_type*2)+ 1


! If clearly in the first half of ny

         if(jjen(i) .le. nyhc) then
     	    do z=1,kjsize
               position = pos0 +(z-1)*jjsz(i)*iisize
               do y=jjst(i),jjen(i)
   		  do x=1,iisize
                     sndbuf(position) = source(x,y,z,j)
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
                     sndbuf(position) = source(x,y,z,j)
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
                     sndbuf(position) = source(x,y,z,j)
                     position = position+1
                  enddo
	       enddo
               do y=ny_fft-nyhc+1,jjen(i)+dny
                  do x=1,iisize
                     sndbuf(position) = source(x,y,z,j)
                     position = position+1
                  enddo
               enddo
            enddo
         endif

      enddo
      enddo

      end subroutine


      subroutine unpack_fcomm2(dest,j,nv)

      implicit none
      integer j,i,x,y,z,nv
      integer*8 position
      complex(p3dfft_type) dest(iisize,jjsize,nz_fft)


!$OMP PARALLEL DO private(i,position,x,y,z)
         do i=0,jproc-1
#ifdef USE_EVEN
            position = i*KfCntMax*nv/(p3dfft_type*2)+1
#else
            position = KfRcvStrt(i)*nv/(p3dfft_type*2)+1
#endif
 	    position = position + (j-1)*iisize*jjsize*kjsz(i)

            do z=kjst(i),kjen(i)
               do y=1,jjsize
                  do x=1,iisize
                     dest(x,y,z) = buf2(position)
                     position = position +1
                  enddo
               enddo
            enddo
         enddo

	 return
	 end subroutine

      subroutine fcomm2(source,dest,t,tc)
!========================================================

      implicit none

      complex(p3dfft_type) source(iisize,ny_fft,kjsize)
      complex(p3dfft_type) dest(iisize,jjsize,nz_fft)

      real(r8) t,tc
      integer x,z,y,i,j,ierr,xs,ys,y2,z2,iy,iz,dny
      integer(i8) position,pos1,pos0



! Pack send buffers for exchanging y and z for all x at once

     call pack_fcomm2(buf1,source)

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
!$OMP PARALLEL DO private(i,j,pos0,position,x,y,z)
         do i=0,jproc-1
	    position = i*KfCntMax/(p3dfft_type*2)  + 1
            do z=kjst(i),kjen(i)
               do y=1,jjsize
                  do x=1,iisize
                     dest(x,y,z) = buf2(position)
                     position = position +1
                  enddo
               enddo
            enddo
!            position = (i+1)*KfCntMax/(p3dfft_type*2)+1
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

      subroutine pack_fcomm2(buf1,source)

      complex(p3dfft_type) source(iisize,ny_fft,kjsize)
      complex(p3dfft_type) buf1(iisize*ny_fft*kjsize)
      integer i,dny,position,pos0,x,y,z

      dny = ny_fft-nyc
      position = 1
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

      return
      end subroutine
