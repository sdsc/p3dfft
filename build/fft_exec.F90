! This file is part of P3DFFT library
!
! Version 2.3
!
! Copyright (C) 2006-2008 Dmitry Pekurovsky
!
!    P3DFFT is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

!----------------------------------------------------------------------------

! This file contains routines for executing pre-initialized 1D FFT operations

      subroutine ftran_y_zplane(In,z_in,xsize_in,zsize_in,stride1_in, &
         stride2_in, Out,z_out,xsize_out,zsize_out,stride1_out,stride2_out, &
         ny,m)

      use p3dfft
      implicit none

      integer xsize_in,zsize_in,xsize_out,zsize_out,stride1_in,stride2_in
      integer stride1_out,stride2_out,ny,m,z_in,z_out
      complex(mytype) In(xsize_in,ny,0:zsize_in-1),Out(xsize_out,ny,0:zsize_out-1)

      call exec_f_c1(In(1,1,z_in),stride1_in,stride2_in,Out(1,1,z_out),stride1_out,stride2_out,ny,m)

      return
      end


      subroutine btran_y_zplane(In,z_in,xsize_in,zsize_in,stride1_in, &
         stride2_in, Out,z_out,xsize_out,zsize_out,stride1_out,stride2_out, &
         ny,m)

      use p3dfft
      implicit none

      integer xsize_in,zsize_in,xsize_out,zsize_out,stride1_in,stride2_in
      integer stride1_out,stride2_out,ny,m,z_in,z_out
      complex(mytype) In(xsize_in,ny,0:zsize_in-1),Out(xsize_out,ny,0:zsize_out-1)

      call exec_b_c1(In(1,1,z_in),stride1_in,stride2_in,Out(1,1,z_out),stride1_out,stride2_out,ny,m)

      return
      end


! Execute backward complex-to-complex 1D FFT

      subroutine exec_b_c1(X,stride_x1,stride_x2,Y,stride_y1, &
          stride_y2,N,m)

      use fft_spec
      use p3dfft
      implicit none

      integer stride_x1,stride_x2,stride_y1,stride_y2,N,m
      complex(mytype) X(N*stride_x1+m*stride_x2),Y(N*stride_y1+m*stride_y2)

#ifdef FFTW
#ifndef SINGLE_PREC
         call dfftw_execute_dft(plan1_bc,X,Y)
#else
         call sfftw_execute_dft(plan1_bc,X,Y)
#endif
#elif defined ESSL

#ifndef SINGLE_PREC
      call dcft (0,X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m,-1,1.0d0, &
              caux1,cnaux,caux2,cnaux)       
#else
      call scft (0,X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m,-1,1.0, &
              caux1,cnaux,caux2,cnaux)       
#endif

#else
      Error: undefined FFT library
#endif
      return
      end

      subroutine exec_b_c2(X,stride_x1,stride_x2,Y,stride_y1, &
          stride_y2,N,m)

      use fft_spec
      use p3dfft
      implicit none

      integer stride_x1,stride_x2,stride_y1,stride_y2,N,m
      complex(mytype) X(N*stride_x1+m*stride_x2),Y(N*stride_y1+m*stride_y2)

#ifdef FFTW
#ifndef SINGLE_PREC
         call dfftw_execute_dft(plan2_bc,X,Y)
#else
         call sfftw_execute_dft(plan2_bc,X,Y)
#endif
#elif defined ESSL

#ifndef SINGLE_PREC
      call dcft (0,X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m,-1,1.0d0, &
              caux1,cnaux,caux2,cnaux)       
#else
      call scft (0,X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m,-1,1.0, &
              caux1,cnaux,caux2,cnaux)       
#endif

#else
      Error: undefined FFT library
#endif
      return
      end

! Execute backward complex-to-real 1D FFT

      subroutine exec_b_c2r(X,dimx,Y,dimy,N,m)

      use fft_spec
      use p3dfft
      implicit none

      integer dimx,dimy,N,m
      complex(mytype) X((N/2+1)*m)
      real(mytype) Y(N*m)

#ifdef FFTW
#ifndef SINGLE_PREC
      call dfftw_execute_dft_c2r(plan1_bcr,X,Y)
#else
      call sfftw_execute_dft_c2r(plan1_bcr,X,Y)
#endif

#elif defined ESSL

#ifndef SINGLE_PREC
      call dcrft(0,X,dimx,Y,dimy,N,m,-1,1.0d0, &
            raux1, rnaux1,raux2,rnaux2,raux3,1)  
#else
      call scrft(0,X,dimx,Y,dimy,N,m,-1,1.0, &
            raux1, rnaux1,raux2,rnaux2,raux3,1)  
#endif
#else
      Error: unknown FFT library
#endif

      return
      end

! Execute forward complex-to-complex 1D FFT

      subroutine exec_f_c1(X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m) 

      use fft_spec
      use p3dfft
      implicit none

      integer N,m,stride_x1,stride_x2,stride_y1,stride_y2
      complex(mytype) X(N*stride_x1+m*stride_x2),Y(N*stride_y1+m*stride_y2)


#ifdef FFTW
#ifndef SINGLE_PREC
      call dfftw_execute_dft(plan1_fc,X,Y)
#else
      call sfftw_execute_dft(plan1_fc,X,Y)
#endif

#elif defined ESSL

#ifndef SINGLE_PREC
      call dcft(0,X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m,1,1.0d0,&
              caux1,cnaux,caux2,cnaux)
#else
      call scft(0,X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m,1,1.0,  &
              caux1,cnaux,caux2,cnaux)
#endif
#else
      Error: undefined FFT library 
#endif

      return
      end

      subroutine exec_f_c2(X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m) 

      use fft_spec
      use p3dfft
      implicit none

      integer N,m,stride_x1,stride_x2,stride_y1,stride_y2
      complex(mytype) X(N*stride_x1+m*stride_x2),Y(N*stride_y1+m*stride_y2)

#ifdef FFTW
#ifndef SINGLE_PREC
      call dfftw_execute_dft(plan2_fc,X,Y)
#else
      call sfftw_execute_dft(plan2_fc,X,Y)
#endif

#elif defined ESSL

#ifndef SINGLE_PREC
      call dcft(0,X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m,1,1.0d0, &
              caux1,cnaux,caux2,cnaux)
#else
      call scft(0,X,stride_x1,stride_x2,Y,stride_y1,stride_y2,N,m,1,1.0,  &
              caux1,cnaux,caux2,cnaux)
#endif
#else
      Error: undefined FFT library 
#endif

      return
      end

! Execute forward real-to-complex 1D FFT

      subroutine exec_f_r2c(X,dimx,Y,dimy,N,m)

      use fft_spec
      use p3dfft

      integer dimx,dimy,N,m
      real(mytype) X(N*m)
      complex(mytype) Y((N/2+1)*m)

#ifdef FFTW

#ifndef SINGLE_PREC
      call dfftw_execute_dft_r2c(plan1_frc,X,Y)
#else
      call sfftw_execute_dft_r2c(plan1_frc,X,Y)
#endif

#elif defined ESSL

#ifndef SINGLE_PREC
      call drcft(0,X,dimx,Y,dimy,N,m,1,1.0d0, &
            raux1, rnaux1,raux2,rnaux2,raux3,1)     
#else
      call srcft(0,X,dimx,Y,dimy,N,m,1,1.0, &
            raux1, rnaux1,raux2,rnaux2,raux3,1)     
#endif
#else
      Error: undefined FFT library
#endif

      return
      end

