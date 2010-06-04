!========================================================
! Transpose  Z-pencils  into Y-pencils

      subroutine bcomm1 (source,dest,t,tc)
!========================================================

      implicit none

      complex(mytype) source(iisize,jjsize,nz_fft)
      complex(mytype) dest(iisize,ny_fft,kjsize)
      real(8) t,tc
      integer x,y,z,i,ierr,xs,ys,iy,iz,y2,z2
      integer(8) position,pos1
      
!     Pack the data for sending


#ifdef USE_EVEN

      if(KfCntUneven) then
         tc = tc - MPI_Wtime()
         position = 1
         do i=0,jproc-1
            do z=kjst(i),kjen(i)
               do y=1,jjsize
                  do x=1,iisize
                     buf1(position) = source(x,y,z)
                     position = position+1
                  enddo
               enddo
            enddo
            position = (i+1)*KfCntMax/(mytype*2)+1
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

      position=1
      do i=0,jproc-1
         do z=1,kjsize
#ifdef STRIDE12
            do x=1,iisize
               do y=jjst(i),jjen(i)
                  dest(y,x,z) = buf2(position)
                  position = position+1
               enddo
            enddo
#else
            do y=jjst(i),jjen(i)
               do x=1,iisize
                  dest(x,y,z) = buf2(position)
                  position = position+1
               enddo
            enddo
#endif
         enddo
#ifdef USE_EVEN
         position = (i+1)*KfCntMax/(mytype*2)+1
#endif
      enddo
      
      
      return
      end subroutine
