      subroutine p3dfft_btran_c2r (XYZg,XgYZ)
!========================================================

      use fft_spec
      implicit none

      real(mytype),TARGET :: XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
#ifdef STRIDE1
      complex(mytype), TARGET :: XYZg(nz_fft,iistart:iiend,jjstart:jjend)
#else
      complex(mytype), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nz_fft)
#endif

      integer x,y,z,i,k,nx,ny,nz
      integer(8) Nl

      if(.not. mpi_set) then
         print *,'P3DFFT error: call setup before other routines'
         return
      endif
      nx = nx_fft
      ny = ny_fft
      nz = nz_fft

      call init_work(nx,ny,nz)

! For FFT libraries that require explicit allocation of work space,
! such as ESSL, initialize here

! Allocate work array

!      allocate(buf(nxhp,jistart:jiend,kjstart:kjend+padd))

! FFT Tranform (C2C) in Z for all x and y 

      if(jproc .gt. 1) then

#ifdef STRIDE1
         call init_b_c(buf, 1,nz, buf, 1, nz,nz,jjsize)
         call bcomm1_trans(XYZg,buf2,buf,timers(3),timers(9))
#else

         if(OW) then

            if(iisize*jjsize .gt. 0) then
               call init_b_c(XYZg, iisize*jjsize, 1, XYZg, iisize*jjsize, 1,nz,iisize*jjsize)

               timers(9) = timers(9) - MPI_Wtime()
               call exec_b_c2(XYZg, iisize*jjsize,1, XYZg, iisize*jjsize, 1,nz,iisize*jjsize)
               timers(9) = timers(9) + MPI_Wtime()
            endif
            call bcomm1(XYZg,buf,timers(3),timers(9))
   
         else
            if(iisize*jjsize .gt. 0) then
               call init_b_c(buf, iisize*jjsize, 1, buf, iisize*jjsize, 1,nz,iisize*jjsize)
            
               timers(9) = timers(9) - MPI_Wtime()
               call exec_b_c2(XYZg, iisize*jjsize,1, buf, iisize*jjsize, 1,nz,iisize*jjsize)
               timers(9) = timers(9) + MPI_Wtime()
               call bcomm1(buf,buf,timers(3),timers(9))
            endif
         endif

#endif


      else
#ifdef STRIDE1
         call reorder_trans_b1(XYZg,buf,buf2)
#else
         if(OW) then   
            Nl = iisize*jjsize*nz
            call init_b_c(XYZg, iisize*jjsize, 1,XYZg, iisize*jjsize, 1,nz,iisize*jjsize)
            call exec_b_c2(XYZg, iisize*jjsize, 1,XYZg, iisize*jjsize, 1,nz,iisize*jjsize)
            call ar_copy(XYZg,buf,Nl)
         else
            call init_b_c(XYZg, iisize*jjsize, 1, buf, iisize*jjsize, 1,nz,iisize*jjsize)
            call exec_b_c2(XYZg, iisize*jjsize, 1, buf, iisize*jjsize, 1,nz,iisize*jjsize)
         endif
#endif
      endif

! Exhange in columns if needed

!
! FFT Transform (C2C) in y dimension for all x, one z-plane at a time
!

      
      if(iisize * kjsize .gt. 0) then

#ifdef STRIDE1
         call init_b_c(buf,1,ny,buf,1,ny,ny,iisize*kjsize)

         timers(10) = timers(10) - MPI_Wtime()

         call exec_b_c1(buf,1,ny,buf,1,ny,ny,iisize*kjsize)

         timers(10) = timers(10) + MPI_Wtime()

#else
         call init_b_c(buf,iisize,1,buf,iisize,1,ny,iisize)
         
         timers(10) = timers(10) - MPI_Wtime()
         do z=kjstart,kjend
               
            call btran_y_zplane(buf,z-kjstart,iisize,kjsize,iisize,1, buf,z-kjstart,iisize,kjsize,iisize,1,ny,iisize)
            
         enddo
         timers(10) = timers(10) + MPI_Wtime()
#endif
      endif

      if(iproc .gt. 1) then 
         call bcomm2(buf,buf,timers(4),timers(11))
#ifdef STRIDE1
      else
         call reorder_b2(buf,buf)
#endif
      endif

! Perform Complex-to-real FFT in x dimension for all y and z
      if(jisize * kjsize .gt. 0) then

         call init_b_c2r(buf,nxhp,XgYZ,nx,nx,jisize*kjsize)

         timers(12) = timers(12) - MPI_Wtime()

         call exec_b_c2r(buf,nxhp,XgYZ,nx,nx,jisize*kjsize)
         timers(12) = timers(12) + MPI_Wtime()

      endif


      call free_work


      return
      end subroutine

      subroutine wrap_exec_b_c2(A,strideA,B,strideB, N,m,L,k)

      complex(mytype) A(L,N),B(L,N)
      integer strideA,strideB,N,m,L,k
 

      call exec_b_c2(A(k,1),strideA,1,B(k,1),strideB,1,N,m)

      return
      end subroutine
