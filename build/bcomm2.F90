!========================================================
! Transpose back Y to X pencils

      subroutine bcomm2(source,dest,t,tc)
!========================================================

      implicit none

      complex(mytype) dest(nxhp,jisize,kjsize)
#ifdef STRIDE1
      complex(mytype) source(ny_fft,iisize,kjsize)
#else
      complex(mytype) source(iisize,ny_fft,kjsize)
#endif
      real(r8) t,tc
      integer x,y,z,i,ierr,ix,iy,x2,y2,l
      integer(i8) position,pos1,pos0,pos2

      tc = tc - MPI_Wtime()

! Pack and exchange x-z buffers in rows
      
      do i=0,iproc-1


#ifdef USE_EVEN
         pos0 = i*IfCntMax/(mytype*2)+1
#else
         pos0 = IfRcvStrt(i)/(mytype*2)+1
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
                        buf1(position) = source(iy,ix,z)
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
                  buf1(position) = source(x,y,z)
                  position = position + 1
               enddo
            enddo
#endif

            pos0 = pos0 + iisize*jisz(i)
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

      do i=0,iproc-1
#ifdef USE_EVEN
         pos0 = i*IfCntMax/(mytype*2) + 1 
#else
         pos0 = IfSndStrt(i)/(mytype*2) + 1
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

      tc = tc + MPI_Wtime() 
      
      return
      end subroutine
