! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2014 Dmitry Pekurovsky
!    Copyright (C) 2006-2014 University of California
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

