!========================================================
! Transpose Y and Z pencils
! Assume Stride1 data structure

      subroutine fcomm2_trans(source,dest,t,tc)
!========================================================

      use fft_spec
      implicit none

! Assume stride12
      complex(mytype) source(ny_fft,iisize,kjsize)
!      complex(mytype) buf(nz_fft*iisize*jjsize)
      complex(mytype) dest(nz_fft,jjsize,iisize)

      real(8) t,tc
      integer x,z,y,i,ierr,xs,ys,y2,z2,iy,iz,ix,x2,n,sz,l
      integer*8 position,pos1,pos0

!	if(taskid .eq. 0) then
!	  print *,'Entring fcomm2_trans'
!        endif	
!	call print_buf(source,ny_fft,iisize,kjsize)

! Pack send buffers for exchanging y and z for all x at once 

      tc = tc - MPI_Wtime()

      do i=0,jproc-1
#ifdef USE_EVEN
         pos0 = i*KfCntMax/(mytype*2)  + 1 
#else
         pos0 = KfSndStrt(i)/(mytype*2)+ 1 
#endif
         do z=1,kjsize
            position = pos0 +(z-1)*jjsz(i)*iisize
            do x=1,iisize
               do y=jjst(i),jjen(i)
                  buf1(position) = source(y,x,z)
                  position = position+1
               enddo
            enddo
         enddo
      enddo
      
      tc = tc + MPI_Wtime()
      t = t - MPI_Wtime()

! Exchange y-z buffers in columns of processors

#ifdef USE_EVEN
      call mpi_alltoall(buf1,KfCntMax, mpi_byte, buf2,KfCntMax, mpi_byte,mpi_comm_col,ierr)

      t = MPI_Wtime() + t
      tc = - MPI_Wtime() + tc

      do x=1,iisize

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
                        dest(iz,iy,x) = buf2(position)
                        position = position + 1
                     enddo
                     pos1 = pos1 + iisize * jjsize
                  enddo
               enddo
               pos0 = pos0 + jjsize*iisize*NBz
           enddo 
         enddo


         call exec_f_c2(dest(1,1,x),1,nz_fft,dest(1,1,x),1,nz_fft, nz_fft,jjsize)
      enddo

      tc = tc + MPI_Wtime()

#else
      call mpi_alltoallv(buf1,KfSndCnts, KfSndStrt,mpi_byte,buf2,KfRcvCnts, KfRcvStrt,mpi_byte,mpi_comm_col,ierr)

      t = MPI_Wtime() + t
      tc = -MPI_Wtime() + tc

      do x=1,iisize

         pos0 = (x-1)*jjsize 

         do z=1,nz_fft,NBz
            z2 = min(z+NBz-1,nz_fft)
            
            do y=1,jjsize,NBy2
               y2 = min(y+NBy2-1,jjsize)
               
               pos1 = pos0 + y 
               
               do iz=z,z2
                  position = pos1
                  do iy=y,y2
! Here we are sure that dest and buf are different
                     dest(iz,iy,x) = buf2(position)
                     position = position +1
                  enddo
                  pos1 = pos1 + iisize * jjsize
               enddo
            enddo
            pos0 = pos0 + iisize*jjsize*NBz
         enddo

         call exec_f_c2(dest(1,1,x),1,nz_fft,dest(1,1,x),1,nz_fft,nz_fft,jjsize)

      enddo
      tc = tc + MPI_Wtime()
         
#endif

      return
      end subroutine

