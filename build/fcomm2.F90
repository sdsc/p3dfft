!========================================================
! Transpose an array in Y=pencils into Z-pencils
! Uses MPI_Alltoall(v)
!
      subroutine fcomm2(source,dest,t,tc)
!========================================================

      implicit none

      complex(mytype) source(iisize,ny_fft,kjsize)
      complex(mytype) dest(iisize,jjsize,nz_fft)

      real(r8) t,tc
      integer x,z,y,i,ierr,xs,ys,y2,z2,iy,iz
      integer(i8) position,pos1

!	if(taskid .eq. 0) then
!	  print *,'Entring fcomm2'
!        endif	
!	call print_buf(source,iisize,ny_fft,kjsize)


! Pack send buffers for exchanging y and z for all x at once 

      position = 1
 
      do i=0,jproc-1
         do z=1,kjsize
            do y=jjst(i),jjen(i)
               do x=1,iisize
                  buf1(position) = source(x,y,z)
                  position = position+1
               enddo
            enddo
         enddo

#ifdef USE_EVEN
         position = position + (KfCntMax/(mytype*2) - jjsz(i)*iisize*kjsize)
#endif
      enddo
      
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
         do i=0,jproc-1
            do z=kjst(i),kjen(i)
               do y=1,jjsize
                  do x=1,iisize
                     dest(x,y,z) = buf2(position)
                     position = position +1
                  enddo
               enddo
            enddo
            position = (i+1)*KfCntMax/(mytype*2)+1
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
