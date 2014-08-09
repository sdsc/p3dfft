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

      subroutine init_plan(A,B,C,n2)

      use fft_spec
      implicit none

      integer(i8) n2
      complex(mytype) A(n2),C(n2)
      real(mytype) B(n2*2)

      
      call init_work(nx_fft,ny_fft,nz_fft)

      if(jisize*kjsize .gt. 0) then

#ifdef DEBUG
	print *,taskid,': doing plan_f_r2c'
#endif
        call plan_f_r2c(B,nx_fft,A,nxhp,nx_fft,jisize*kjsize) 

#ifdef DEBUG
	print *,taskid,': doing plan_b_c2r'
#endif
      call plan_b_c2r(A,nxhp,B,nx_fft,nx_fft,jisize*kjsize) 

     endif

#ifdef STRIDE1

     if(iisize*kjsize .gt. 0) then 
     
#ifdef DEBUG
	print *,taskid, ': doing plan_f_c1'
#endif
      call plan_f_c1(A,1,ny_fft,A,1,ny_fft,ny_fft,iisize*kjsize)
#ifdef DEBUG
	print *,taskid,': doing plan_b_c1'
#endif
      call plan_b_c1(A,1,ny_fft,A,1,ny_fft,ny_fft,iisize*kjsize)

     endif	

     if(jjsize .gt. 0) then
#ifdef DEBUG
	print *,taskid,': doing plan_f_c2 & plan_b_c2'
#endif

        call plan_b_c2_same(A,1,nz_fft,A,1,nz_fft,nz_fft,jjsize)
        call plan_b_c2_dif(A,1,nz_fft,C,1,nz_fft,nz_fft,jjsize)
        call plan_f_c2_same(A,1,nz_fft, A,1,nz_fft,nz_fft,jjsize)
        call plan_f_c2_dif(A,1,nz_fft, C,1,nz_fft,nz_fft,jjsize)
        call plan_ctrans_r2_same (A, 2,2*nz_fft, &
                         A, 2,2*nz_fft, NZ_fft, jjsize)
        call plan_strans_r2_same (A, 2,2*nz_fft, &
                         A, 2,2*nz_fft, NZ_fft, jjsize)
        call plan_ctrans_r2_dif (A, 2,2*nz_fft, &
                         C, 2,2*nz_fft, NZ_fft, jjsize)
        call plan_strans_r2_dif (A, 2,2*nz_fft, &
                         C, 2,2*nz_fft, NZ_fft, jjsize)

     endif
#else
     if(iisize .gt. 0) then
#ifdef DEBUG
      print *,taskid,': doing plan_f_c1'
#endif
      call plan_f_c1(A,iisize,1,A,iisize,1,ny_fft,iisize)
#ifdef DEBUG
      print *,taskid,': doing plan_b_c1'
#endif
      call plan_b_c1(A,iisize,1,A,iisize,1,ny_fft,iisize)

       if(jjsize .gt. 0) then

#ifdef DEBUG
          print *,taskid,': doing plan_f_c2'
#endif
          call plan_f_c2_same(A,iisize*jjsize, 1, A,iisize*jjsize, 1,nz_fft,iisize*jjsize)
#ifdef DEBUG
          print *,taskid,': doing plan_b_c2'
#endif
        call plan_b_c2_same(A,iisize*jjsize, 1, A,iisize*jjsize, 1,nz_fft,iisize*jjsize)

       endif

!cccccccccccccccccccccc added chebyshev cccccccccccccccccccccccccccccccc
    call plan_ctrans_r2_same (A, 2*iisize*jjsize, 1, &
                         A, 2 * iisize * jjsize, 1, NZ_fft, 2 * iisize * jjsize)
    call plan_strans_r2_same (A, 2*iisize*jjsize, 1, &
                         A, 2 * iisize * jjsize, 1, NZ_fft, 2 * iisize * jjsize)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     endif
#endif


#ifdef DEBUG
    print *,taskid,': Finished init_plan'
#endif

      return
      end subroutine
