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

      module p3dfft

      implicit none

      include 'mpif.h'

      private

! Set precision

#ifndef SINGLE_PREC
       integer, parameter,public :: mytype=8
       integer, parameter,public:: mpireal = MPI_DOUBLE_PRECISION
       integer,parameter,public:: mpicomplex = MPI_DOUBLE_COMPLEX
#else
       integer, parameter,public :: mytype=4
       integer, parameter,public:: mpireal = MPI_REAL
       integer,parameter,public:: mpicomplex = MPI_COMPLEX
#endif

! global variables

      integer, save,public :: padd,num_thr
      real(8), save,public :: timers(12)

      integer,save :: NX_fft,NY_fft,NZ_fft,numtasks,iproc,jproc
      integer,save :: ipid,jpid,taskid
      integer,save :: iistart,iiend,iisize,jjstart,jjsize,jjend
      integer,save ::jistart,kjstart,jisize,kjsize,jiend,kjend

      integer,save ::  nxh,nxhp

! mpi process info
!
      logical :: mpi_set=.false.
      integer, save :: mpi_comm_cart      
      integer, save :: mpi_comm_row, mpi_comm_col
      integer,save, dimension(:), allocatable :: iist,iien,iisz
      integer,save, dimension(:), allocatable :: jist,jien,jisz
      integer,save, dimension(:), allocatable :: jjst,jjen,jjsz
      integer,save, dimension(:), allocatable :: kjst,kjen,kjsz

! mpi derived data types for implementing alltoallv using send-recvs
      integer,save,dimension(:),allocatable:: IfSndCnts,IfSndStrt
      integer,save,dimension(:),allocatable:: IfRcvCnts,IfRcvStrt
      integer,save,dimension(:),allocatable:: KfSndCnts,KfSndStrt
      integer,save,dimension(:),allocatable:: KfRcvCnts,KfRcvStrt
      integer,save,dimension(:),allocatable:: JrSndCnts,JrSndStrt
      integer,save,dimension(:),allocatable:: JrRcvCnts,JrRcvStrt
      integer,save,dimension(:),allocatable:: KrSndCnts,KrSndStrt
      integer,save,dimension(:),allocatable:: KrRcvCnts,KrRcvStrt
      integer,save,dimension(:,:),allocatable:: status
      complex(mytype), save, allocatable :: buf(:),buf1(:),buf2(:)
      logical :: OW = .false.
#ifdef USE_EVEN
      integer(8),save :: IfCntMax,KfCntMax
      logical KfCntUneven
#endif

#ifdef STRIDE1
      integer CB,NBx,NBy1,NBy2,NBz
#endif


#ifdef STRIDE1
#ifdef NBL_X
      integer, parameter :: NB1 = NBL_X
#else
      integer, parameter :: NB1=4
#endif
#ifdef NBL_Y
      integer, parameter :: NB=NBL_Y
#else
      integer, parameter :: NB=1
#endif
#endif

! ghost-cell support using p3dfft_init_ghosts(), p3dfft_update_ghosts()
      logical :: ghosts_set=.false.
      integer,save, dimension(:),   allocatable :: proc_id2coords
      integer,save, dimension(:,:), allocatable :: proc_coords2id
      integer,save, dimension(:), allocatable :: gneighb_r,gneighb_c
      real(kind=mytype),save, dimension(:), allocatable :: gbuf_snd,gbuf_recv
      integer,save :: goverlap
      integer,save :: gmess_rsize(3), gmess_csize(3)
      integer,save :: gproc_rsize(3), gproc_csize(3)
      integer,save :: gproc_rstart(3), gproc_cstart(3)
      integer,save :: gproc_rend(3), gproc_cend(3)
      integer,save :: gslab_rsize(3,3), gslab_rstart(3,6), gslab_rend(3,6)
      integer,save :: gslab_csize(3,3), gslab_cstart(3,6), gslab_cend(3,6)
      integer(kind=8),save :: gmem_rstart(6), gmem_cstart(6)

      public :: get_dims,p3dfft_setup,p3dfft_ftran_r2c,p3dfft_btran_c2r, &
                 p3dfft_clean,print_buf,&
                 p3dfft_init_ghosts, update_rghosts,  &
                 update_cghosts, gr_ijk2i, &
                 proc_id2coords, proc_coords2id, &
                 gmem_rstart, gmem_cstart

!-------------------
      contains
!-------------------

#include "setup.F90"
#include "init_plan.F90"
#include "ftran.F90"
#include "btran.F90"
#ifdef STRIDE1
#include "reorder.F90"
#endif
#include "fcomm1.F90"
#ifdef STRIDE1
#include "fcomm2_trans.F90"
#include "bcomm1_trans.F90"
#else
#include "fcomm2.F90"
#include "bcomm1.F90"
#endif
#include "bcomm2.F90"
#include "ghost_cell.F90"

!=====================================================
! Return array dimensions for either real-space (conf=1) or wavenumber-space(conf=2)
! 
      subroutine get_dims(istart,iend,isize,conf)
!=====================================================

      integer istart(3),iend(3),isize(3),conf
      
      if(.not. mpi_set) then
         print *,'P3DFFT error: call setup before other routines'
      else

      if(conf .eq. 1) then
         istart(1) = 1
         iend(1) = NX_fft
         isize(1) = NX_fft
         istart(2) = jistart
         iend(2) = jiend
         isize(2) = jisize
         istart(3) = kjstart
         iend(3) = kjend
         isize(3) = kjsize
      else if(conf .eq. 2) then
#ifdef STRIDE1
         istart(3) = iistart
         iend(3) = iiend
         isize(3) = iisize
         istart(2) = jjstart
         iend(2) = jjend
         isize(2) = jjsize
         istart(1) = 1
         iend(1) = NZ_fft
         isize(1) = NZ_fft
#else
         istart(1) = iistart
         iend(1) = iiend
         isize(1) = iisize
         istart(2) = jjstart
         iend(2) = jjend
         isize(2) = jjsize
         istart(3) = 1
         iend(3) = NZ_fft
         isize(3) = NZ_fft
#endif
      endif

      endif
      end subroutine get_dims


!========================================================
      subroutine p3dfft_clean

!!--------------------------------------------------------------
! Clean-up routines for FFTW

      use fft_spec

#ifdef FFTW
#ifndef SINGLE_PREC
      call dfftw_destroy_plan(plan1_frc)      
      call dfftw_destroy_plan(plan1_bcr)      
      call dfftw_destroy_plan(plan1_fc)      
      call dfftw_destroy_plan(plan2_fc)      
      call dfftw_destroy_plan(plan1_bc)      
      call dfftw_destroy_plan(plan2_bc)      
#else
      call sfftw_destroy_plan(plan1_frc)      
      call sfftw_destroy_plan(plan1_bcr)      
      call sfftw_destroy_plan(plan1_fc)      
      call sfftw_destroy_plan(plan2_fc)      
      call sfftw_destroy_plan(plan1_bc)      
      call sfftw_destroy_plan(plan2_bc)      
#endif

#elif defined ESSL
      deallocate(caux1)      
      deallocate(caux2)
      deallocate(raux1)
      deallocate(raux2)
#endif

      deallocate(buf1)
      deallocate(buf2)
      deallocate(buf)

      return
      end subroutine

!========================================================
      subroutine ar_copy(A,B,nar)

      integer(8) nar,i
      complex(mytype) A(nar,1,1),B(nar,1,1)

      do i=1,nar
         B(i,1,1)=A(i,1,1)
      enddo

      return
      end subroutine



!========================================================
      subroutine print_buf(A,lx,ly,lz)

      complex(mytype) A(lx,ly,lz)
      integer lx,ly,lz,i,j,k
      
      do k=1,lz
         do j=1,ly
            do i=1,lx
               if(abs(A(i,j,k)) .gt. 0.0000005) then
                  print *,taskid,': (',i,j,k,') =',A(i,j,k)
               endif
            enddo
         enddo
      enddo

      end subroutine


      end module


