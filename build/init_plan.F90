!========================================================

      subroutine init_plan(A,B,n1)

      use fft_spec
      implicit none

      complex(mytype) A(n1)
      real(mytype) B(n1*2)
      integer(8) n1

      
      call init_work(nx_fft,ny_fft,nz_fft)
      call plan_f_r2c(B,nx_fft,A,nxhp,nx_fft,jisize*kjsize,.false.) 
      call plan_b_c2r(A,nxhp,B,nx_fft,nx_fft,jisize*kjsize,.false.) 
#ifdef STRIDE1
      call plan_f_c1(A,1,ny_fft,A,1,ny_fft,ny_fft,iisize*kjsize,.false.)
      call plan_b_c1(A,1,ny_fft,A,1,ny_fft,ny_fft,iisize*kjsize,.false.)
         call plan_f_c2(A,1,nz_fft, &
           A,1,nz_fft,nz_fft,jjsize,.false.)
      if(OW) then
         call plan_b_c2(A,1,nz_fft,A,1,nz_fft,nz_fft,jjsize,.false.)
      else
         call plan_b_c2(A,1,nz_fft,B,1,nz_fft,nz_fft,jjsize,.false.)
      endif
#else
      call plan_f_c1(A,iisize,1,A,iisize,1,ny_fft,iisize,.false.)
      call plan_b_c1(A,iisize,1,A,iisize,1,ny_fft,iisize,.false.)
      call plan_f_c2(A,iisize*jjsize, 1, A,iisize*jjsize, 1,nz_fft,iisize*jjsize,.false.)
      if(OW) then
         call plan_b_c2(A,iisize*jjsize, 1, A,iisize*jjsize, 1,nz_fft,iisize*jjsize,.false.)
      else
         call plan_b_c2(A,iisize*jjsize, 1, B,iisize*jjsize, 1,nz_fft,iisize*jjsize,.false.)
      endif
#endif

      return
      end subroutine
