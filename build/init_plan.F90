! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2010 Dmitry Pekurovsky
!    Copyright (C) 2006-2010 University of California
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

      subroutine init_plan(A,B,n2)

      use fft_spec
      implicit none

      integer(i8) n2
      complex(mytype) A(n2)
      real(mytype) B(n2*2)

      
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
