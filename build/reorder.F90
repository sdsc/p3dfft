! This routine is called only when jproc=1, and only when stride1 is used
! transform backward in Z and transpose array in memory 

      subroutine reorder_trans_b1(A,B,C)

      use fft_spec

      complex(mytype) B(iisize,ny_fft,nz_fft)
      complex(mytype) A(nz_fft,iisize,ny_fft)
      complex(mytype) C(nz_fft,iisize,NBy2*iisize*nz_fft)
      integer x,y,z,iy,iz,y2,z2

      if(OW) then
         do y=1,ny_fft,nby2
            y2 = min(y+nby2-1,ny_fft)
            call exec_b_c2(A(1,1,y),1,nz_fft,A(1,1,y),1,nz_fft,nz_fft,nby2*iisize)
            do z=1,nz_fft,nbz
               z2 = min(z+nbz-1,nz_fft)
               do iz=z,z2
                  do iy=y,y2
                     do x=1,iisize
                        B(x,iy,iz) = A(iz,x,iy)
                     enddo
                  enddo
               enddo
               
            enddo                  
         enddo

      else
         do y=1,ny_fft,nby2
            y2 = min(y+nby2-1,ny_fft)
            call exec_b_c2(A(1,1,y),1,nz_fft,C(1,1,y),1,nz_fft,nz_fft,nby2*iisize)
            do z=1,nz_fft,nbz
               z2 = min(z+nbz-1,nz_fft)
               do iz=z,z2
                  do iy=y,y2
                     do x=1,iisize
                        B(x,iy,iz) = C(iz,x,iy)
                     enddo
                  enddo
               enddo
               
            enddo                  
         enddo
         
      endif

      return
      end subroutine

! This routine is called only when iproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Y

      subroutine reorder_b2(A,B)

      implicit none

      complex(mytype) B(nxhp,ny_fft,kjsize)
      complex(mytype) A(ny_fft,nxhp,kjsize)
      integer x,y,z,iy,x2,ix,y2
      complex(mytype), allocatable :: tmp(:,:)


      allocate(tmp(nxhp,0:ny_fft-1))

      do z=1,kjsize
         do x=1,nxhp,nbx
            x2 = min(x+nbx-1,nxhp)
            do y=0,ny_fft-1,nby1
               y2 = min(y+nby1-1,ny_fft-1)
               do ix=x,x2
                  do iy = y,y2
                     tmp(ix,iy) = A(iy+1,ix,z)
                  enddo
               enddo
            enddo
         enddo

         do y=0,ny_fft-1
            do x=1,nxhp
               B(x,y+1,z) = tmp(x,y)
            enddo
         enddo
      enddo

      deallocate(tmp)

      return
      end subroutine

! This routine is called only when iproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Y

      subroutine reorder_f1(A,B)

      implicit none

      complex(mytype) A(nxhp,ny_fft,kjsize)
      complex(mytype) B(ny_fft,nxhp,kjsize)
      integer x,y,z,iy,x2,ix,y2
      complex(mytype), allocatable :: tmp(:,:)

      allocate(tmp(ny_fft-1,nxhp))


      do z=1,kjsize
         do y=0,ny_fft-1,nby1
            y2 = min(y+nby1-1,ny_fft-1)
            do x=1,nxhp,nbx
               x2 = min(x+nbx-1,nxhp)
               do iy = y,y2
                  do ix=x,x2
                     tmp(iy,ix) = A(ix,iy+1,z)
                  enddo
               enddo
            enddo
         enddo

         do x=1,nxhp
            do y=0,ny_fft-1
               B(y+1,x,z) = tmp(y,x)
            enddo
         enddo
      enddo

      deallocate(tmp)

      return
      end subroutine

! This routine is called only when jproc=1, and only when stride1 is used
! Transpose array in memory and transform forward in Z

      subroutine reorder_trans_f2(A,B)

      use fft_spec
      implicit none

      complex(mytype) A(iisize,ny_fft,nz_fft)
      complex(mytype) B(nz_fft,iisize,ny_fft)
      integer x,y,z,iy,iz,y2,z2

      do y=1,ny_fft,nby2
         y2 = min(y+nby2-1,ny_fft)
         do z=1,nz_fft,nbz
            z2 = min(z+nbz-1,nz_fft)
            do iz=z,z2
               do iy=y,y2
                  do x=1,iisize
                     B(iz,x,iy) = A(x,iy,iz)
                  enddo
               enddo
            enddo
         enddo                  
         call exec_f_c2(B(1,1,y),1,nz_fft,B(1,1,y),1,nz_fft,nz_fft,nby2*iisize)
      enddo

      return
      end subroutine
