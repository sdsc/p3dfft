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
! Transpose X and Y pencils

      subroutine fcomm1_many(source,dest,nv,t,tc)
!========================================================

      implicit none

      complex(p3dfft_type) source(nxhp,jisize,kjsize,nv)
#ifdef STRIDE1
      complex(p3dfft_type) dest(ny_fft,iisize,kjsize,nv)
#else
      complex(p3dfft_type) dest(iisize,ny_fft,kjsize,nv)
#endif

      real(r8) t,tc
      integer x,y,i,ierr,z,xs,j,n,ix,iy,y2,x2,l,nv
      integer(i8) position,pos1,pos0,pos2
      integer sndcnts(0:iproc-1)
      integer rcvcnts(0:iproc-1)
      integer sndstrt(0:iproc-1)
      integer rcvstrt(0:iproc-1)

!	if(taskid .eq. 0) then
!	  print *,'Entering fcomm1'
!        endif
!	call print_buf(source,nxhp,jisize,kjsize)

! Pack the send buffer for exchanging y and x (within a given z plane ) into sendbuf

      tc = tc - MPI_Wtime()

!$OMP PARALLEL DO private(i,position,x,y,z)
      do i=0,iproc-1
#ifdef USE_EVEN
         position = i*IfCntMax * nv/(p3dfft_type*2) + 1
#else
         position = IfSndStrt(i) *nv/(p3dfft_type*2) + 1
#endif
	do j=1,nv
           do z=1,kjsize
              do y=1,jisize
                 do x=iist(i),iien(i)
                    buf1(position) = source(x,y,z,j)
                    position = position +1
                 enddo
              enddo
            enddo
	 enddo
      enddo
      tc = tc + MPI_Wtime()
      t = t - MPI_Wtime()

#ifdef USE_EVEN

! Use MPI_Alltoall
! Exchange the y-x buffers (in rows of processors)

      call mpi_alltoall(buf1,IfCntMax*nv, mpi_byte, buf2,IfCntMax*nv, mpi_byte,mpi_comm_row,ierr)

#else
! Use MPI_Alltoallv
! Exchange the y-x buffers (in rows of processors)
      sndcnts = IfSndCnts * nv
      sndstrt = IfSndStrt * nv
      rcvcnts = IfRcvCnts * nv
      rcvstrt = IfRcvStrt * nv
      call mpi_alltoallv(buf1,SndCnts, SndStrt,mpi_byte, buf2,RcvCnts, RcvStrt,mpi_byte,mpi_comm_row,ierr)
#endif

      t = MPI_Wtime() + t
      tc = - MPI_Wtime() + tc

! Unpack the data

!$OMP PARALLEL DO private(i,j,pos0,pos1,pos2,ix,iy,x2,y2,position,x,y,z) collapse(2)
      do i=0,iproc-1

	do j=1,nv
#ifdef USE_EVEN
         pos0 = (IfCntMax*nv*i +(j-1) * IfRcvCnts(i))/(p3dfft_type*2) +1
#else
         pos0 = (IfRcvStrt(i) * nv + (j-1) * IfRcvCnts(i))/(p3dfft_type*2) +1
#endif
         do z=1,kjsize
            pos1 = pos0 + (z-1)*iisize*jisz(i)

#ifdef STRIDE1
            do y=jist(i),jien(i),nby1
               y2 = min(y+nby1-1,jien(i))
               do x=1,iisize,nbx
                  x2 = min(x+nbx-1,iisize)
                  pos2 = pos1 + x-1
                  do iy = y,y2
                     position = pos2
                     do ix=x,x2
                        dest(iy,ix,z,j) = buf2(position)
                        position = position + 1
                     enddo
                     pos2 = pos2 + iisize
                  enddo
               enddo
               pos1 = pos1 + iisize*nby1
            enddo
#else
            position = pos1
            do y=jist(i),jien(i)
               do x=1,iisize
                  dest(x,y,z,j) = buf2(position)
                  position = position + 1
               enddo
            enddo
#endif

         enddo
      enddo
      enddo

      tc = tc + MPI_Wtime()


      return
      end subroutine


      subroutine fcomm1(source,dest,t,tc)
!========================================================

      implicit none

      complex(p3dfft_type) source(nxhp,jisize,kjsize)
#ifdef STRIDE1
      complex(p3dfft_type) dest(ny_fft,iisize,kjsize)
#else
      complex(p3dfft_type) dest(iisize,ny_fft,kjsize)
#endif

      real(r8) t,tc
      integer x,y,i,ierr,z,xs,j,n,ix,iy,y2,x2,l
      integer(i8) position,pos1,pos0,pos2

!	if(taskid .eq. 0) then
!	  print *,'Entering fcomm1'
!        endif
!	call print_buf(source,nxhp,jisize,kjsize)

! Pack the send buffer for exchanging y and x (within a given z plane ) into sendbuf

#ifdef DEBUG
       print *,taskid,': fcomm1: packing data'
#endif

      tc = tc - MPI_Wtime()

!$OMP PARALLEL DO private(i,position,x,y,z)
      do i=0,iproc-1
#ifdef USE_EVEN
         position = i*IfCntMax/(p3dfft_type*2) + 1
#else
         position = IfSndStrt(i)/(p3dfft_type*2) + 1
#endif
         do z=1,kjsize
            do y=1,jisize
               do x=iist(i),iien(i)
                  buf1(position) = source(x,y,z)
                  position = position +1
               enddo
            enddo
         enddo
      enddo
      tc = tc + MPI_Wtime()
      t = t - MPI_Wtime()

#ifdef DEBUG
       print *,taskid,': fcomm1: initiating exchange'
#endif

#ifdef USE_EVEN

! Use MPI_Alltoall
! Exchange the y-x buffers (in rows of processors)

      call mpi_alltoall(buf1,IfCntMax, mpi_byte, buf2,IfCntMax, mpi_byte,mpi_comm_row,ierr)

#else
! Use MPI_Alltoallv
! Exchange the y-x buffers (in rows of processors)
      call mpi_alltoallv(buf1,IfSndCnts, IfSndStrt,mpi_byte, buf2,IfRcvCnts, IfRcvStrt,mpi_byte,mpi_comm_row,ierr)
#endif

      t = MPI_Wtime() + t
      tc = - MPI_Wtime() + tc

! Unpack the data
#ifdef DEBUG
       print *,taskid,': fcomm1: unpacking data'
#endif


!$OMP PARALLEL DO private(i,position,x,y,z,pos0,pos1,pos2,iy,y2,ix,x2) collapse(2)
      do i=0,iproc-1
         do z=1,kjsize
#ifdef USE_EVEN
            pos1 = i*IfCntMax/(p3dfft_type*2)+1 + (z-1)*iisize*jisz(i)
#else
	    pos1 = IfRcvStrt(i)/(p3dfft_type*2)+1 + (z-1)*iisize*jisz(i)
#endif

#ifdef STRIDE1
            do y=jist(i),jien(i),nby1
               y2 = min(y+nby1-1,jien(i))
               do x=1,iisize,nbx
                  x2 = min(x+nbx-1,iisize)
                  pos2 = pos1 + x-1
                  do iy = y,y2
                     position = pos2
                     do ix=x,x2
                        dest(iy,ix,z) = buf2(position)
                        position = position + 1
                     enddo
                     pos2 = pos2 + iisize
                  enddo
               enddo
               pos1 = pos1 + iisize*nby1
            enddo
#else
            position = pos1
            do y=jist(i),jien(i)
               do x=1,iisize
                  dest(x,y,z) = buf2(position)
                  position = position + 1
               enddo
            enddo
#endif

         enddo
      enddo

      tc = tc + MPI_Wtime()


      return
      end subroutine

