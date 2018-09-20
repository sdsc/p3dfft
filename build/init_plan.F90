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

      subroutine init_plan

      use fft_spec
      implicit none

      integer(i8) n2
!      complex(p3dfft_type) A(n2),C(n2)
!      real(p3dfft_type) B(n2*2)
       real(p3dfft_type), allocatable :: B(:)
       complex(p3dfft_type), allocatable :: A(:),C(:)
       integer omp_get_num_threads,omp_get_thread_num,l,m,tid,ierr

#ifdef OPENMP

!$OMP PARALLEL
!$OMP MASTER
      num_thr = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL

#else
      num_thr = 1
#endif

      if(taskid .eq. 0) then
         print *,'Running on ',num_thr,'threads'
      endif

      call init_work(nx_fft,ny_fft,nz_fft)

#ifdef FFTW
        allocate(plan1_frc(0:num_thr-1),plan1_bcr(0:num_thr-1),plan1_fc(0:num_thr-1),plan1_bc(0:num_thr-1))
	allocate(plan_ctrans_same(0:num_thr-1), plan_strans_same(0:num_thr-1))
	allocate(plan_ctrans_dif(0:num_thr-1), plan_strans_dif(0:num_thr-1))
	allocate(plan2_bc_same(0:num_thr-1),plan2_fc_same(0:num_thr-1),plan2_bc_dif(0:num_thr-1),plan2_fc_dif(0:num_thr-1))
	allocate(startx_frc(0:num_thr-1),startx_bcr(0:num_thr-1),startx_f_c1(0:num_thr-1),startx_b_c1(0:num_thr-1))
	allocate(startx_ctrans_same(0:num_thr-1), startx_strans_same(0:num_thr-1))
	allocate(startx_ctrans_dif(0:num_thr-1), startx_strans_dif(0:num_thr-1))
	allocate(startx_b_c2_same(0:num_thr-1),startx_f_c2_same(0:num_thr-1),startx_b_c2_dif(0:num_thr-1),startx_f_c2_dif(0:num_thr-1))
	allocate(starty_frc(0:num_thr-1),starty_bcr(0:num_thr-1),starty_f_c1(0:num_thr-1),starty_b_c1(0:num_thr-1))
	allocate(starty_ctrans_same(0:num_thr-1), starty_strans_same(0:num_thr-1))
	allocate(starty_ctrans_dif(0:num_thr-1), starty_strans_dif(0:num_thr-1))
	allocate(starty_b_c2_same(0:num_thr-1),starty_f_c2_same(0:num_thr-1),starty_b_c2_dif(0:num_thr-1),starty_f_c2_dif(0:num_thr-1))

      if(jisize*kjsize .gt. 0) then

      l = mod(jisize*kjsize,num_thr)
      m = jisize*kjsize/num_thr
      startx_frc(0) = 1
      starty_frc(0) = 1
      startx_bcr(0) = 1
      starty_bcr(0) = 1
      do tid=0,l-1
         startx_frc(tid+1) = startx_frc(tid) + (m+1)*nx_fft
         starty_frc(tid+1) = starty_frc(tid) + (m+1)*nxhp
         starty_bcr(tid+1) = starty_bcr(tid) + (m+1)*nx_fft
         startx_bcr(tid+1) = startx_bcr(tid) + (m+1)*nxhp
      enddo
      do tid=l,num_thr-2
         startx_frc(tid+1) = startx_frc(tid) + m*nx_fft
         starty_frc(tid+1) = starty_frc(tid) + m*nxhp
         starty_bcr(tid+1) = starty_bcr(tid) + m*nx_fft
         startx_bcr(tid+1) = startx_bcr(tid) + m*nxhp
      enddo

#ifdef DEBUG
	print *,taskid,': doing plan_f_r2c and plan_b_c2r'
#endif

!!$OMP PARALLEL private(tid,A,B) shared(nx_fft,m,l,plan1_frc,plan1_bcr)
!      tid = omp_get_thread_num()
      allocate(A(nxhp*(m+1)),B(nx_fft*(m+1)))

      do tid=0,num_thr-1
#ifndef SINGLE_PREC
      if(tid .lt. l) then
         call dfftw_plan_many_dft_r2c(plan1_frc(tid),1,nx_fft,m+1, &
           B,NULL,1,nx_fft, A,NULL,1,nxhp,fftw_flag)
         call dfftw_plan_many_dft_c2r(plan1_bcr(tid),1,nx_fft,m+1, &
              A,NULL,1,nxhp, B,NULL,1,nx_fft,fftw_flag)
      else
         call dfftw_plan_many_dft_r2c(plan1_frc(tid),1,nx_fft,m, &
           B,NULL,1,nx_fft, A,NULL,1,nxhp,fftw_flag)
         call dfftw_plan_many_dft_c2r(plan1_bcr(tid),1,nx_fft,m, &
              A,NULL,1,nxhp, B,NULL,1,nx_fft,fftw_flag)
      endif
#else
      if(tid .lt. l) then
         call sfftw_plan_many_dft_r2c(plan1_frc(tid),1,nx_fft,m+1, &
           B,NULL,1,nx_fft, A,NULL,1,nxhp,fftw_flag)
         call sfftw_plan_many_dft_c2r(plan1_bcr(tid),1,nx_fft,m+1, &
              A,NULL,1,nxhp, B,NULL,1,nx_fft,fftw_flag)
      else
         call sfftw_plan_many_dft_r2c(plan1_frc(tid),1,nx_fft,m, &
           B,NULL,1,nx_fft, A,NULL,1,nxhp,fftw_flag)
         call sfftw_plan_many_dft_c2r(plan1_bcr(tid),1,nx_fft,m, &
              A,NULL,1,nxhp, B,NULL,1,nx_fft,fftw_flag)
      endif
#endif
      enddo

      deallocate(A,B)

! !$OMP END PARALLEL
     endif

#ifdef STRIDE1

#ifdef DEBUG
	print *,taskid, ': doing plan_f_c1'
#endif

     if(iisize*kjsize .gt. 0) then

      l = mod(iisize*kjsize,num_thr)
      m = iisize*kjsize/num_thr

      startx_f_c1(0) = 1
      starty_f_c1(0) = 1
      startx_b_c1(0) = 1
      starty_b_c1(0) = 1
      do tid=0,l-1
         startx_f_c1(tid+1) = startx_f_c1(tid) + (m+1)*ny_fft
         starty_f_c1(tid+1) = starty_f_c1(tid) + (m+1)*ny_fft
         starty_b_c1(tid+1) = starty_b_c1(tid) + (m+1)*ny_fft
         startx_b_c1(tid+1) = startx_b_c1(tid) + (m+1)*ny_fft
      enddo
      do tid=l,num_thr-2
         startx_f_c1(tid+1) = startx_f_c1(tid) + m*ny_fft
         starty_f_c1(tid+1) = starty_f_c1(tid) + m*ny_fft
         starty_b_c1(tid+1) = starty_b_c1(tid) + m*ny_fft
         startx_b_c1(tid+1) = startx_b_c1(tid) + m*ny_fft
      enddo

!!$OMP PARALLEL private(tid,A)
!      tid = omp_get_thread_num()
      allocate(A(ny_fft*(m+1)))

      do tid=0,num_thr-1
#ifndef SINGLE_PREC
      if(tid .lt. l) then
         call dfftw_plan_many_dft(plan1_fc(tid),1,ny_fft,m+1, A,NULL,1,ny_fft, &
           A,NULL,1,ny_fft,FFTW_FORWARD,fftw_flag)
      else
         call dfftw_plan_many_dft(plan1_fc(tid),1,ny_fft,m, A,NULL,1,ny_fft, &
           A,NULL,1,ny_fft,FFTW_FORWARD,fftw_flag)
      endif
#else
      if(tid .lt. l) then
         call sfftw_plan_many_dft(plan1_fc(tid),1,ny_fft,m+1, A,NULL,1,ny_fft, &
           A,NULL,1,ny_fft,FFTW_FORWARD,fftw_flag)
      else
         call sfftw_plan_many_dft(plan1_fc(tid),1,ny_fft,m, A,NULL,1,ny_fft, &
           A,NULL,1,ny_fft,FFTW_FORWARD,fftw_flag)
      endif
#endif

!#ifdef DEBUG
!!$OMP MASTER
!	print *,taskid,': doing plan_b_c1'
!!$OMP END MASTER
!#endif


#ifndef SINGLE_PREC
      if(tid .lt. l) then
         call dfftw_plan_many_dft(plan1_bc(tid),1,ny_fft,m+1, A,NULL,1,ny_fft, &
           A,NULL,1,ny_fft,FFTW_BACKWARD,fftw_flag)
      else
         call dfftw_plan_many_dft(plan1_bc(tid),1,ny_fft,m, A,NULL,1,ny_fft, &
           A,NULL,1,ny_fft,FFTW_BACKWARD,fftw_flag)
      endif
#else
      if(tid .lt. l) then
         call sfftw_plan_many_dft(plan1_bc(tid),1,ny_fft,m+1, A,NULL,1,ny_fft, &
           A,NULL,1,ny_fft,FFTW_BACKWARD,fftw_flag)
      else
         call sfftw_plan_many_dft(plan1_bc(tid),1,ny_fft,m, A,NULL,1,ny_fft, &
           A,NULL,1,ny_fft,FFTW_BACKWARD,fftw_flag)
      endif
#endif

     enddo

     deallocate(A)
! !$OMP END PARALLEL
     endif

     if(jjsize .gt. 0) then

      startx_b_c2_same = 1
      startx_b_c2_dif = 1
      startx_f_c2_same = 1
      startx_f_c2_dif = 1
      starty_b_c2_same = 1
      starty_b_c2_dif = 1
      starty_f_c2_same = 1
      starty_f_c2_dif = 1
      startx_ctrans_same = 1
      startx_strans_same = 1
      startx_ctrans_dif = 1
      startx_strans_dif = 1
      starty_ctrans_same = 1
      starty_strans_same = 1
      starty_ctrans_dif = 1
      starty_strans_dif = 1

!      l = mod(jjsize,num_thr)
      m = jjsize

#ifdef DEBUG
	print *,taskid,': doing plan_f_c2 & plan_b_c2'
#endif

!!$OMP PARALLEL private(tid,A,C)
!      tid = omp_get_thread_num()
      allocate(A(nz_fft*m),C(nz_fft*m))

      do tid=0,num_thr-1
#ifndef SINGLE_PREC
#ifdef DEBUG
	print *,taskid,': doing plan_b_c2'
#endif
         call dfftw_plan_many_dft(plan2_bc_same(tid),1,nz_fft,m, A,NULL,1, &
	    nz_fft, A,NULL,1,nz_fft,FFTW_BACKWARD,fftw_flag)
         call dfftw_plan_many_dft(plan2_bc_dif(tid),1,nz_fft,m, A,NULL,1, &
	    nz_fft, C,NULL,1,nz_fft,FFTW_BACKWARD,fftw_flag)
#ifdef DEBUG
	print *,taskid,': doing plan_f_c2'
#endif
         call dfftw_plan_many_dft(plan2_fc_same(tid),1,nz_fft,m, A,NULL,1, &
	    nz_fft, A,NULL,1,nz_fft,FFTW_FORWARD,fftw_flag)
         call dfftw_plan_many_dft(plan2_fc_dif(tid),1,nz_fft,m, A,NULL,1, &
	    nz_fft, C,NULL,1,nz_fft,FFTW_FORWARD,fftw_flag)
#ifdef DEBUG
	print *,taskid,': doing plan_ctrans'
#endif
         call dfftw_plan_many_r2r (plan_ctrans_same(tid), 1, nz_fft, m, &
                              A, NULL, 2,2*nz_fft, &
                              A, NULL, 2,2*nz_fft, FFTW_REDFT00, fftw_flag)
         call dfftw_plan_many_r2r (plan_ctrans_dif(tid), 1, nz_fft, m, &
                              A, NULL, 2,2*nz_fft, &
                              C, NULL, 2,2*nz_fft, FFTW_REDFT00, fftw_flag)
#ifdef DEBUG
	print *,taskid,': doing plan_strans'
#endif
         call dfftw_plan_many_r2r (plan_strans_same(tid), 1, nz_fft, m, &
                              A, NULL, 2,2*nz_fft, &
                              A, NULL, 2,2*nz_fft, FFTW_RODFT00, fftw_flag)
         call dfftw_plan_many_r2r (plan_strans_dif(tid), 1, nz_fft, m, &
                              A, NULL, 2,2*nz_fft, &
                              C, NULL, 2,2*nz_fft, FFTW_RODFT00, fftw_flag)
#else
         call sfftw_plan_many_dft(plan2_bc_same(tid),1,nz_fft,m, A,NULL,1, &
	    nz_fft, A,NULL,1,nz_fft,FFTW_BACKWARD,fftw_flag)
         call sfftw_plan_many_dft(plan2_bc_dif(tid),1,nz_fft,m, A,NULL,1, &
	    nz_fft, C,NULL,1,nz_fft,FFTW_BACKWARD,fftw_flag)
         call sfftw_plan_many_dft(plan2_fc_same(tid),1,nz_fft,m, A,NULL,1, &
	    nz_fft, A,NULL,1,nz_fft,FFTW_FORWARD,fftw_flag)
         call sfftw_plan_many_dft(plan2_fc_dif(tid),1,nz_fft,m, A,NULL,1, &
	    nz_fft, C,NULL,1,nz_fft,FFTW_FORWARD,fftw_flag)
         call sfftw_plan_many_r2r (plan_ctrans_same(tid), 1, nz_fft, m, &
                              A, NULL, 2,2*nz_fft, &
                              A, NULL, 2,2*nz_fft, FFTW_REDFT00, fftw_flag)
         call sfftw_plan_many_r2r (plan_strans_same(tid), 1, nz_fft, m, &
                              A, NULL, 2,2*nz_fft, &
                              A, NULL, 2,2*nz_fft, FFTW_RODFT00, fftw_flag)
         call sfftw_plan_many_r2r (plan_ctrans_dif(tid), 1, nz_fft, m, &
                              A, NULL, 2,2*nz_fft, &
                              C, NULL, 2,2*nz_fft, FFTW_REDFT00, fftw_flag)
         call sfftw_plan_many_r2r (plan_strans_dif(tid), 1, nz_fft, m, &
                              A, NULL, 2,2*nz_fft, &
                              C, NULL, 2,2*nz_fft, FFTW_RODFT00, fftw_flag)
#endif
     enddo

     deallocate(A,C)
! !$OMP END PARALLEL
     endif

#else   ! STRIDE != 1

     if(iisize .gt. 0) then

      l = mod(iisize,num_thr)
      m = iisize/num_thr
#ifdef DEBUG
      print *,taskid,': doing plan_f_c1 and plan_b_c1',m,l
#endif

      startx_f_c1(0) = 1
      starty_f_c1(0) = 1
      startx_b_c1(0) = 1
      starty_b_c1(0) = 1
      do tid=0,l-1
         startx_f_c1(tid+1) = startx_f_c1(tid) + (m+1)
         starty_f_c1(tid+1) = starty_f_c1(tid) + (m+1)
         starty_b_c1(tid+1) = starty_b_c1(tid) + (m+1)
         startx_b_c1(tid+1) = startx_b_c1(tid) + (m+1)
      enddo
      do tid=l,num_thr-2
         startx_f_c1(tid+1) = startx_f_c1(tid) + m
         starty_f_c1(tid+1) = starty_f_c1(tid) + m
         starty_b_c1(tid+1) = starty_b_c1(tid) + m
         startx_b_c1(tid+1) = startx_b_c1(tid) + m
      enddo

!!$OMP PARALLEL private(tid,A)
!      tid = omp_get_thread_num()
      allocate(A(ny_fft*(iisize+1)))
      A = 0.

      do tid=0,num_thr-1
#ifndef SINGLE_PREC
      if(tid .lt. l) then
         call dfftw_plan_many_dft(plan1_fc(tid),1,ny_fft,m+1, A,NULL,iisize,1, &
           A,NULL,iisize,1,FFTW_FORWARD,fftw_flag)
         call dfftw_plan_many_dft(plan1_bc(tid),1,ny_fft,m+1, A,NULL,iisize,1, &
           A,NULL,iisize,1,FFTW_BACKWARD,fftw_flag)
      else
         call dfftw_plan_many_dft(plan1_fc(tid),1,ny_fft,m, A,NULL,iisize,1, &
           A,NULL,iisize,1,FFTW_FORWARD,fftw_flag)
         call dfftw_plan_many_dft(plan1_bc(tid),1,ny_fft,m, A,NULL,iisize,1, &
           A,NULL,iisize,1,FFTW_BACKWARD,fftw_flag)
      endif
#else
      if(tid .lt. l) then
         call sfftw_plan_many_dft(plan1_fc(tid),1,ny_fft,m+1, A,NULL,iisize,1, &
           A,NULL,iisize,1,FFTW_FORWARD,fftw_flag)
         call sfftw_plan_many_dft(plan1_bc(tid),1,ny_fft,m+1, A,NULL,iisize,1, &
           A,NULL,iisize,1,FFTW_BACKWARD,fftw_flag)
      else
         call sfftw_plan_many_dft(plan1_fc(tid),1,ny_fft,m, A,NULL,iisize,1, &
           A,NULL,iisize,1,FFTW_FORWARD,fftw_flag)
         call sfftw_plan_many_dft(plan1_bc(tid),1,ny_fft,m, A,NULL,iisize,1, &
           A,NULL,iisize,1,FFTW_BACKWARD,fftw_flag)
      endif
#endif
	enddo

#ifdef DEBUG
       print *,'plan1_fc=',plan1_fc
       print *,'plan1_bc=',plan1_bc
       print *,'plan1_frc=',plan1_frc
#endif
	deallocate(A)

! !$OMP END PARALLEL

       if(jjsize .gt. 0) then

#ifdef DEBUG
          print *,taskid,': doing plan_f_c2 and plan_b_c2'
#endif

      l = mod(iisize*jjsize,num_thr)
      m = iisize*jjsize/num_thr
      startx_f_c2_same(0) = 1
      starty_f_c2_same(0) = 1
      startx_b_c2_same(0) = 1
      starty_b_c2_same(0) = 1
      startx_ctrans_same(0) = 1
      startx_strans_same(0) = 1
      starty_ctrans_same(0) = 1
      starty_strans_same(0) = 1
      do tid=0,l-1
         startx_f_c2_same(tid+1) = startx_f_c2_same(tid) + (m+1)
         starty_f_c2_same(tid+1) = starty_f_c2_same(tid) + (m+1)
         starty_b_c2_same(tid+1) = starty_b_c2_same(tid) + (m+1)
         startx_b_c2_same(tid+1) = startx_b_c2_same(tid) + (m+1)
	 startx_ctrans_same(tid+1) = startx_ctrans_same(tid) + (m+1)*2
	 startx_strans_same(tid+1) = startx_strans_same(tid) + (m+1)*2
	 starty_ctrans_same(tid+1) = starty_ctrans_same(tid) + (m+1)*2
	 starty_strans_same(tid+1) = starty_strans_same(tid) + (m+1)*2
      enddo
      do tid=l,num_thr-2
         startx_f_c2_same(tid+1) = startx_f_c2_same(tid) + m
         starty_f_c2_same(tid+1) = starty_f_c2_same(tid) + m
         starty_b_c2_same(tid+1) = starty_b_c2_same(tid) + m
         startx_b_c2_same(tid+1) = startx_b_c2_same(tid) + m
	 startx_ctrans_same(tid+1) = startx_ctrans_same(tid) + m*2
	 startx_strans_same(tid+1) = startx_strans_same(tid) + m*2
	 starty_ctrans_same(tid+1) = starty_ctrans_same(tid) + m*2
	 starty_strans_same(tid+1) = starty_strans_same(tid) + m*2
      enddo

!!$OMP PARALLEL private(tid,A)
!      tid = omp_get_thread_num()
      allocate(A((iisize*jjsize+1)*nz_fft))

      do tid=0,num_thr-1
#ifndef SINGLE_PREC
      if(tid .lt. l) then
         call dfftw_plan_many_dft(plan2_bc_same(tid),1,nz_fft,m+1, A,NULL, &
	    iisize*jjsize,1, A,NULL,iisize*jjsize,1,FFTW_BACKWARD,fftw_flag)
         call dfftw_plan_many_dft(plan2_fc_same(tid),1,nz_fft,m+1, A,NULL, &
	    iisize*jjsize,1, A,NULL,iisize*jjsize,1,FFTW_FORWARD,fftw_flag)
         call dfftw_plan_many_r2r (plan_ctrans_same(tid), 1, nz_fft, m+1, &
                              A, NULL, 2*iisize*jjsize,2, &
                              A, NULL, 2*iisize*jjsize,2, FFTW_REDFT00, fftw_flag)
         call dfftw_plan_many_r2r (plan_strans_same(tid), 1, nz_fft, m+1, &
                              A, NULL, 2*iisize*jjsize,2, &
                              A, NULL, 2*iisize*jjsize,2, FFTW_RODFT00, fftw_flag)
      else
         call dfftw_plan_many_dft(plan2_bc_same(tid),1,nz_fft,m, A,NULL, &
	    iisize*jjsize,1, A,NULL,iisize*jjsize,1,FFTW_BACKWARD,fftw_flag)
         call dfftw_plan_many_dft(plan2_fc_same(tid),1,nz_fft,m, A,NULL, &
	    iisize*jjsize,1, A,NULL,iisize*jjsize,1,FFTW_FORWARD,fftw_flag)
         call dfftw_plan_many_r2r (plan_ctrans_same(tid), 1, nz_fft, m, &
                              A, NULL, 2*iisize*jjsize,2, &
                              A, NULL, 2*iisize*jjsize,2, FFTW_REDFT00, fftw_flag)
         call dfftw_plan_many_r2r (plan_strans_same(tid), 1, nz_fft, m, &
                              A, NULL, 2*iisize*jjsize,2, &
                              A, NULL, 2*iisize*jjsize,2, FFTW_RODFT00, fftw_flag)

      endif
#else
      if(tid .lt. l) then
         call sfftw_plan_many_dft(plan2_bc_same(tid),1,nz_fft,m+1, A,NULL, &
	    iisize*jjsize,1, A,NULL,iisize*jjsize,1,FFTW_BACKWARD,fftw_flag)
         call sfftw_plan_many_dft(plan2_fc_same(tid),1,nz_fft,m+1, A,NULL, &
	    iisize*jjsize,1, A,NULL,iisize*jjsize,1,FFTW_FORWARD,fftw_flag)
         call sfftw_plan_many_r2r (plan_ctrans_same(tid), 1, nz_fft, m+1, &
                              A, NULL, 2*iisize*jjsize,2, &
                              A, NULL, 2*iisize*jjsize,2, FFTW_REDFT00, fftw_flag)
         call sfftw_plan_many_r2r (plan_strans_same(tid), 1, nz_fft, m+1, &
                              A, NULL, 2*iisize*jjsize,2, &
                              A, NULL, 2*iisize*jjsize,2, FFTW_RODFT00, fftw_flag)
      else
         call sfftw_plan_many_dft(plan2_bc_same(tid),1,nz_fft,m, A,NULL, &
	    iisize*jjsize,1, A,NULL,iisize*jjsize,1,FFTW_BACKWARD,fftw_flag)
         call sfftw_plan_many_dft(plan2_fc_same(tid),1,nz_fft,m, A,NULL, &
	    iisize*jjsize,1, A,NULL,iisize*jjsize,1,FFTW_FORWARD,fftw_flag)
         call sfftw_plan_many_r2r (plan_ctrans_same(tid), 1, nz_fft, m, &
                              A, NULL, 2*iisize*jjsize,2, &
                              A, NULL, 2*iisize*jjsize,2, FFTW_REDFT00, fftw_flag)
         call sfftw_plan_many_r2r (plan_strans_same(tid), 1, nz_fft, m, &
                              A, NULL, 2*iisize*jjsize,2, &
                              A, NULL, 2*iisize*jjsize,2, FFTW_RODFT00, fftw_flag)

      endif
#endif
	enddo
	deallocate(A)
! !$OMP END PARALLEL

       endif

     endif
#endif


#ifdef DEBUG
    print *,taskid,': Finished init_plan'
#endif

#endif


      return
      end subroutine
