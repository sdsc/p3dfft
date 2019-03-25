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
! Transpose back Y to X pencils

      subroutine bcomm2_many(source,dest,nv,t,tc)
!========================================================

      implicit none

      complex(p3dfft_type) dest(nxhp,jisize,kjsize,nv)
#ifdef STRIDE1
      complex(p3dfft_type) source(ny_fft,iisize,kjsize,nv)
#else
      complex(p3dfft_type) source(iisize,ny_fft,kjsize,nv)
#endif
      real(r8) t,tc
      integer x,y,z,i,ierr,ix,iy,x2,y2,l,j,nv,dim
      integer(i8) position,pos1,pos0,pos2
      integer sndcnts(0:iproc-1)
      integer rcvcnts(0:iproc-1)
      integer sndstrt(0:iproc-1)
      integer rcvstrt(0:iproc-1)

      tc = tc - MPI_Wtime()

! Pack and exchange x-z buffers in rows

!$OMP PARALLEL DO private(i,j,pos0,pos1,pos2,ix,iy,x2,y2,position,x,y,z) collapse(2)
      do i=0,iproc-1

	do j=1,nv

#ifdef USE_EVEN
         pos0 = (IfCntMax*nv*i +(j-1) * KrSndCnts(i))/(p3dfft_type*2) +1
#else
         pos0 = (KrSndStrt(i) * nv + (j-1) * KrSndCnts(i))/(p3dfft_type*2) +1
#endif

         do z=1,kjsize

#ifdef STRIDE1
            pos1 = pos0
            do y=jist(i),jien(i),nby1
               y2 = min(y+nby1-1,jien(i))
               do x=1,iisize,nbx
                  x2 = min(x+nbx-1,iisize)
                  pos2 = pos1 + x-1
                  do iy = y,y2
                     position = pos2
                     do ix=x,x2
                        buf1(position) = source(iy,ix,z,j)
                        position = position + 1
                     enddo
                     pos2 = pos2 + iisize
                  enddo
               enddo
               pos1 = pos1 + iisize*nby1
            enddo

#else
            position = pos0
            do y=jist(i),jien(i)
               do x=1,iisize
                  buf1(position) = source(x,y,z,j)
                  position = position + 1
               enddo
            enddo
#endif

            pos0 = pos0 + iisize*jisz(i)
         enddo
      enddo
      enddo

      tc = tc + MPI_Wtime()
      t = t - MPI_Wtime()

#ifdef USE_EVEN
      call mpi_alltoall (buf1,IfCntMax*nv,mpi_byte, buf2,IfCntMax*nv,mpi_byte,mpi_comm_row,ierr)
#else
      sndcnts = KrSndCnts * nv
      sndstrt = KrSndStrt * nv
      rcvcnts = KrRcvCnts * nv
      rcvstrt = KrRcvStrt * nv

      call mpi_alltoallv (buf1,SndCnts, SndStrt, mpi_byte, buf2,RcvCnts,RcvStrt,mpi_byte,mpi_comm_row,ierr)
#endif

      t = t + MPI_Wtime()
      tc = tc - MPI_Wtime()

! Unpack receive buffers into dest

!$OMP PARALLEL DO private(i,j,pos0,position,x,y,z)
      do i=0,iproc-1

#ifdef USE_EVEN
         pos0 = i*IfCntMax*nv/(p3dfft_type*2) + 1
#else
         pos0 = IfSndStrt(i)*nv/(p3dfft_type*2) + 1
#endif

	do j=1,nv
         do z=1,kjsize
            position = pos0
            do y=1,jisize
               do x=iist(i),iien(i)
                  dest(x,y,z,j) = buf2(position)
                  position = position +1
               enddo
            enddo
            pos0 = pos0 + iisz(i)*jisize
         enddo
	 enddo
      enddo
      do j=1,nv
      do z=1,kjsize
         do y=1,jisize
	    do x=nxhpc+1,nxhp
	       dest(x,y,z,j) = 0.
	    enddo
	 enddo
      enddo
      enddo

      tc = tc + MPI_Wtime()

      return
      end subroutine

      subroutine bcomm2(source,dest,t,tc)
!========================================================

      implicit none

      complex(p3dfft_type) dest(nxhp,jisize,kjsize)
#ifdef STRIDE1
      complex(p3dfft_type) source(ny_fft,iisize,kjsize)
#else
      complex(p3dfft_type) source(iisize,ny_fft,kjsize)
#endif
      real(r8) t,tc
      integer x,y,z,i,ierr,ix,iy,x2,y2,l
      integer(i8) position,pos1,pos0,pos2

      tc = tc - MPI_Wtime()

! Pack and exchange x-z buffers in rows

!$OMP PARALLEL DO private(i,pos0,pos1,pos2,ix,iy,x2,y2,position,x,y,z) collapse(2)
      do i=0,iproc-1
         do z=1,kjsize

#ifdef STRIDE1
#ifdef USE_EVEN
         pos1 = i*IfCntMax/(p3dfft_type*2)+1 +(z-1)*iisize*jisz(i)
#else
         pos1 = IfRcvStrt(i)/(p3dfft_type*2)+1 +(z-1)*iisize*jisz(i)
#endif
            do y=jist(i),jien(i),nby1
               y2 = min(y+nby1-1,jien(i))
               do x=1,iisize,nbx
                  x2 = min(x+nbx-1,iisize)
                  pos2 = pos1 + x-1
                  do iy = y,y2
                     position = pos2
                     do ix=x,x2
                        buf1(position) = source(iy,ix,z)
                        position = position + 1
                     enddo
                     pos2 = pos2 + iisize
                  enddo
               enddo
               pos1 = pos1 + iisize*nby1
            enddo

#else
#ifdef USE_EVEN
         position = i*IfCntMax/(p3dfft_type*2)+1 + (z-1)*iisize*jisz(i)
#else
         position = IfRcvStrt(i)/(p3dfft_type*2)+1+(z-1)*iisize*jisz(i)
#endif
            do y=jist(i),jien(i)
               do x=1,iisize
                  buf1(position) = source(x,y,z)
                  position = position + 1
               enddo
            enddo
#endif
         enddo
      enddo

      tc = tc + MPI_Wtime()
      t = t - MPI_Wtime()

#ifdef USE_EVEN
      call mpi_alltoall (buf1,IfCntMax,mpi_byte, buf2,IfCntMax,mpi_byte,mpi_comm_row,ierr)
#else
      call mpi_alltoallv (buf1,KrSndCnts, KrSndStrt, mpi_byte, buf2,KrRcvCnts,KrRcvStrt,mpi_byte,mpi_comm_row,ierr)
#endif

      t = t + MPI_Wtime()
      tc = tc - MPI_Wtime()

! Unpack receive buffers into dest

!$OMP PARALLEL DO private(i,pos0,position,x,y,z)
      do i=0,iproc-1
#ifdef USE_EVEN
         pos0 = i*IfCntMax/(p3dfft_type*2) + 1
#else
         pos0 = IfSndStrt(i)/(p3dfft_type*2) + 1
#endif
         do z=1,kjsize
            position = pos0
            do y=1,jisize
               do x=iist(i),iien(i)
                  dest(x,y,z) = buf2(position)
                  position = position +1
               enddo
            enddo
            pos0 = pos0 + iisz(i)*jisize
         enddo
      enddo
      do z=1,kjsize
         do y=1,jisize
	    do x=nxhpc+1,nxhp
	       dest(x,y,z) = 0.
	    enddo
	 enddo
      enddo

      tc = tc + MPI_Wtime()

      return
      end subroutine
