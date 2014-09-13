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

      module p3dfft

      implicit none

      include 'mpif.h'

      private

! Set precision

#ifndef SINGLE_PREC
       integer, parameter,public :: mytype=KIND(1.0d0)
       integer, parameter,public:: mpireal = MPI_DOUBLE_PRECISION
       integer,parameter,public:: mpicomplex = MPI_DOUBLE_COMPLEX
#else
       integer, parameter,public :: mytype=KIND(1.0)
       integer, parameter,public:: mpireal = MPI_REAL
       integer,parameter,public:: mpicomplex = MPI_COMPLEX
#endif

! global variables

      integer, parameter, public :: r8 = KIND(1.0d0)
      integer, parameter, public :: i8 = SELECTED_INT_KIND(18)
      integer, save, public :: num_thr,padi
      real(r8), save,public :: timers(12)
      real(r8), save :: timer(12)
       integer, public :: real_size,complex_size

      integer,save :: NX_fft,NY_fft,NZ_fft,nxh,nxhp,nv_preset
      integer,save :: nxc,nyc,nzc,nxhc,nxhpc,nyh,nzh,nyhc,nzhc,nyhcp,nzhcp	
      integer,save :: ipid,jpid,taskid,numtasks,iproc,jproc
      integer,save :: iistart,iiend,iisize,jjstart,jjsize,jjend
      integer,save ::jistart,kjstart,jisize,kjsize,jiend,kjend
    integer, save :: ijstart, ijsize, ijend, iiistart, iiisize, iiiend

    integer, save :: maxisize, maxjsize, maxksize


! mpi process info
!
      logical :: mpi_set=.false.
      integer, save :: mpi_comm_cart,mpicomm      
      integer, save :: mpi_comm_row, mpi_comm_col
      integer,save, dimension(:), allocatable :: iist,iien,iisz
      integer,save, dimension(:), allocatable :: jist,jien,jisz
      integer,save, dimension(:), allocatable :: jjst,jjen,jjsz
      integer,save, dimension(:), allocatable :: kjst,kjen,kjsz
    integer, save, dimension (:), allocatable :: iiist, iiien, iiisz
    integer, save, dimension (:), allocatable :: ijst, ijen, ijsz

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
      integer, save, dimension (:), allocatable :: IiCnts, IiStrt
      integer, save, dimension (:), allocatable :: IjCnts, IjStrt
      integer, save, dimension (:), allocatable :: JiCnts, JiStrt
      integer, save, dimension (:), allocatable :: KjCnts, KjStrt
 
#ifdef USE_EVEN
    integer * 8, save :: IfCntMax, KfCntMax
    integer * 8, save :: IJCntMax, JICntMax
    integer * 8, save :: IKCntMax, KICntMax
    integer * 8, save :: IiCntMax, KjCntMax
    logical KfCntUneven
#endif
    integer*8 :: nm

#ifdef STRIDE1
      integer CB,NBx,NBy1,NBy2,NBz
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

! trans2proc support
    integer, save, dimension (:), allocatable :: proc_id2coords
    integer, save, dimension (:, :), allocatable :: proc_coords2id
    integer, save, dimension (:, :, :), allocatable :: proc_dims
    integer, save, dimension (:, :), allocatable :: proc_parts
 
    public :: p3dfft_get_dims, p3dfft_get_mpi_info, p3dfft_setup, &
		p3dfft_ftran_r2c, p3dfft_btran_c2r, p3dfft_cheby, &
		p3dfft_ftran_r2c_many, p3dfft_btran_c2r_many, p3dfft_cheby_many, &
		get_timers,set_timers,&
              p3dfft_clean, print_buf, print_buf_real, &
              proc_id2coords, proc_coords2id, &
              proc_dims, proc_parts, get_proc_parts, &
              rtran_x2y, rtran_y2x, rtran_x2z, rtran_z2x, &
              p3dfft_ftran_r2c_1d

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
!#include "wrap.F90"
!#include "ghost_cell.F90"


!=====================================================
! this is a C wrapper routine
      subroutine p3dfft_get_dims_w(istart,iend,isize,conf) BIND(C,NAME='p3dfft_get_dims')
!=====================================================
      integer istart(3),iend(3),isize(3),conf

      call p3dfft_get_dims(istart,iend,isize,conf) 

      end subroutine

!=====================================================
! Return array dimensions for either real-space (conf=1) or wavenumber-space(conf=2)
! 
      subroutine p3dfft_get_dims(istart,iend,isize,conf) 
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
         iend(1) = NZc
         isize(1) = NZc
#else
         istart(1) = iistart
         iend(1) = iiend
         isize(1) = iisize
         istart(2) = jjstart
         iend(2) = jjend
         isize(2) = jjsize
         istart(3) = 1
         iend(3) = NZc
         isize(3) = NZc
#endif
        else if (conf == 3) then
          istart = (/ 0, 0, 0 /)
          iend = (/ maxisize, maxjsize, maxksize /)
          isize = (/ maxisize, maxjsize, maxksize /)
      endif

      endif
      end subroutine p3dfft_get_dims

! --------------------------------------
!
!  get_mpi_info
!
! --------------------------------------
    subroutine p3dfft_get_mpi_info (mpi_taskid, mpi_tasks, mpi_comm)
      implicit none
 
!         function args
      integer, intent (out) :: mpi_taskid
      integer, intent (out) :: mpi_tasks
      integer, intent (out) :: mpi_comm
 
      if ( .not. mpi_set) then
        print *, 'P3DFFT error: call setup before other routines'
        return
      end if
 
      mpi_taskid = taskid
      mpi_tasks = numtasks
      mpi_comm = mpi_comm_cart
 
    end subroutine 	

! this is a C wrapper routine
!========================================================
      subroutine p3dfft_clean_w() BIND(C,NAME='p3dfft_clean')
!========================================================

      call p3dfft_clean

      end subroutine 

!========================================================
      subroutine p3dfft_clean() 
!========================================================

!!--------------------------------------------------------------
! Clean-up routines for FFTW

      use fft_spec

#ifdef FFTW
#ifndef SINGLE_PREC
      call dfftw_destroy_plan(plan1_frc)      
      call dfftw_destroy_plan(plan1_bcr)      
      call dfftw_destroy_plan(plan1_fc)      
      call dfftw_destroy_plan(plan2_fc_same)      
      call dfftw_destroy_plan(plan2_fc_dif)      
      call dfftw_destroy_plan(plan1_bc)      
      call dfftw_destroy_plan(plan2_bc_same)      
      call dfftw_destroy_plan(plan2_bc_dif)      
      call dfftw_destroy_plan(plan_ctrans_same)      
      call dfftw_destroy_plan(plan_ctrans_dif)      
      call dfftw_destroy_plan(plan_strans_same)      
      call dfftw_destroy_plan(plan_strans_dif)      
#else
      call sfftw_destroy_plan(plan1_frc)      
      call sfftw_destroy_plan(plan1_bcr)      
      call sfftw_destroy_plan(plan1_fc)      
      call sfftw_destroy_plan(plan1_bc)      
      call sfftw_destroy_plan(plan2_fc_same)      
      call sfftw_destroy_plan(plan2_fc_dif)      
      call sfftw_destroy_plan(plan2_bc_same)      
      call sfftw_destroy_plan(plan2_bc_dif)      
      call sfftw_destroy_plan(plan_ctrans_same)      
      call sfftw_destroy_plan(plan_ctrans_dif)      
      call sfftw_destroy_plan(plan_strans_same)      
      call sfftw_destroy_plan(plan_strans_dif)      
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

    deallocate (iist)
    deallocate (iisz)
    deallocate (iien)
    deallocate (jjst)
    deallocate (jjsz)
    deallocate (jjen)
    deallocate (jist)
    deallocate (jisz)
    deallocate (jien)
    deallocate (kjst)
    deallocate (kjsz)
    deallocate (kjen)
  
    deallocate (IfSndStrt)
    deallocate (IfSndCnts)
    deallocate (IfRcvStrt)
    deallocate (IfRcvCnts)
    deallocate (KfSndStrt)
    deallocate (KfSndCnts)
    deallocate (KfRcvStrt)
    deallocate (KfRcvCnts)
    deallocate (JrSndStrt)
    deallocate (JrSndCnts)
    deallocate (JrRcvStrt)
    deallocate (JrRcvCnts)
    deallocate (KrSndStrt)
    deallocate (KrSndCnts)
    deallocate (KrRcvStrt)
    deallocate (KrRcvCnts)
 
    deallocate (IiCnts)
    deallocate (IiStrt)
    deallocate (IjCnts)
    deallocate (IjStrt)
    deallocate (JiCnts)
    deallocate (JiStrt)
    deallocate (KjCnts)
    deallocate (KjStrt)
 
    deallocate (proc_id2coords)
    deallocate (proc_coords2id)
    deallocate (proc_dims)
    deallocate (proc_parts)
 
    mpi_set = .false.

      return
      end subroutine

!========================================================
      subroutine ar_copy_many(A,dim_a,B,dim_b,nar,nv)
!========================================================

      integer(i8) nar,i
      integer dim_a,dim_b,nv,j
      complex(mytype) A(dim_a,nv),B(dim_b,nv)

      do j=1,nv
         call ar_copy(A(1,j),B(1,j),nar)
      enddo
      
      return
      end subroutine


!========================================================
      subroutine ar_copy(A,B,nar)
!========================================================

      integer(i8) nar,i
      complex(mytype) A(nar,1,1),B(nar,1,1)

      do i=1,nar
         B(i,1,1)=A(i,1,1)
      enddo

      return
      end subroutine

!========================================================
      subroutine seg_zero_z_many(A,xdim,ydim,z1,z2,zdim,dim,nv)
!========================================================

      implicit none
      integer x,y,z,xdim,ydim,zdim,z1,z2,nv,j,dim
      complex(mytype) A(dim,nv)

      do j=1,nv
         call seg_zero_z(A(1,j),xdim,ydim,z1,z2,zdim)
      enddo

      return
      end subroutine

!========================================================
      subroutine seg_zero_z(A,xdim,ydim,z1,z2,zdim)
!========================================================

      implicit none
      integer x,y,z,xdim,ydim,zdim,z1,z2
      complex(mytype) A(xdim,ydim,zdim)

      do z=z1,z2
         do y=1,ydim
	    do x=1,xdim
	       A(x,y,z) = 0.
	    enddo
	 enddo
      enddo

      return
      end subroutine seg_zero_z

!========================================================
    subroutine seg_copy_z_f_many(in,out,x1,x2,y1,y2,z1,z2,shift_z,xdim,ydim,zdim,dim,nv)
!========================================================

    implicit none
    integer x1,x2,y1,y2,z1,z2,xdim,ydim,zdim,shift_z,x,y,z,nv,j,dim
    complex(mytype) in(xdim,ydim,zdim,nv), out(dim,nv)

    do j=1,nv
      call seg_copy_z(in(1,1,1,j),out(1,j),x1,x2,y1,y2,z1,z2,shift_z,xdim,ydim,zdim)    
    enddo

    return
    end subroutine

!========================================================
    subroutine seg_copy_z_b_many(in,out,x1,x2,y1,y2,z1,z2,shift_z,xdim,ydim,zdim,dim,nv)
!========================================================

    implicit none
    integer x1,x2,y1,y2,z1,z2,xdim,ydim,zdim,shift_z,x,y,z,nv,j,dim
    complex(mytype) out(xdim,ydim,zdim,nv), in(dim,nv)

    do j=1,nv
      call seg_copy_z(in(1,j),out(1,1,1,j),x1,x2,y1,y2,z1,z2,shift_z,xdim,ydim,zdim)    
    enddo

    return
    end subroutine

!========================================================
    subroutine seg_copy_z(in,out,x1,x2,y1,y2,z1,z2,shift_z,xdim,ydim,zdim)
!========================================================

    implicit none
    integer x1,x2,y1,y2,z1,z2,xdim,ydim,zdim,shift_z,x,y,z
    complex(mytype) in(xdim,ydim,zdim), out(xdim,ydim,zdim)


    do z=z1,z2
       do y=y1,y2
          do x=x1,x2
	     out(x,y,z) = in(x,y,z+shift_z)
	  enddo	
        enddo	
    enddo

    return
    end subroutine seg_copy_z


!========================================================
      subroutine get_timers_w(timer) BIND(C,name='get_timers')
!========================================================

      real(r8) timer(12)

      call get_timers(timer)

      return
      end subroutine

      subroutine set_timers_w() BIND(C,name='set_timers')

      call set_timers

      return
      end subroutine

      subroutine get_timers(timer)
         real(r8) timer(12)
         timer(1) = timers(1)
         timer(2) = timers(2)
         timer(3) = timers(3)
         timer(4) = timers(4)
         timer(5) = timers(5)
         timer(6) = timers(6)
         timer(7) = timers(7)
         timer(8) = timers(8)
         timer(9) = timers(9)
         timer(10) = timers(10)
         timer(11) = timers(11)
         timer(12) = timers(12)
      end subroutine

      subroutine set_timers()
         timers = 0
      end subroutine

!========================================================
      subroutine print_buf(A,lx,ly,lz)
!========================================================

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
!========================================================
      subroutine print_buf_real(A,lx,ly,lz)
!========================================================

      real(mytype) A(lx,ly,lz)
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

! --------------------------------------
!
!  proc_neighb(..)
!
! --------------------------------------
  function proc_neighb (base_proc_id, orient, direction)
    implicit none
 
!         function args
    integer, intent (in) :: base_proc_id, orient, direction
    integer :: proc_neighb
 
    integer :: coord_i, coord_j
 
    proc_neighb = - 1
 
!         check input
    if (base_proc_id < 0 .or. base_proc_id >= iproc*jproc) then
      return
    end if
    if (orient /= 1 .and. orient /=-1) then
      return
    end if
    if (direction /= 1 .and. direction /= 2) then
      return
    end if
 
    coord_i = proc_id2coords (base_proc_id*2)
    coord_j = proc_id2coords (base_proc_id*2+1)
 
!         direction=1:i, direction=2:j, orient=-1 or 1
    if (direction == 1 .and. &
       (coord_i+orient) >= 0 .and. &
       (coord_i+orient) <= iproc) then
    proc_neighb = proc_coords2id (coord_i+orient, coord_j)
  else if (direction == 2 .and. &
       (coord_j+orient) >= 0 .and. &
       (coord_j+orient) <= jproc) then
    proc_neighb = proc_coords2id (coord_i, coord_j+orient)
  end if
 
  return
end function proc_neighb
 
! --------------------------------------
!
!  search_proc(..)
!
! --------------------------------------
function search_proc (point_i, point_j, dimI_offset, dimJ_offset, conf)
  implicit none
 
!         function args
  integer, intent (in) :: point_i, point_j, dimI_offset, dimJ_offset
  integer, intent (in) :: conf
  integer :: search_proc
 
  integer :: proc_id
 
  search_proc = - 1
 
!         check input
  if (dimI_offset < 1 .or. dimI_offset > 3 .or. &
      dimJ_offset < 1 .or. dimJ_offset > 3 .or. &
      conf < 1 .or. conf > 2) then
    return
  end if
 
  proc_id = proc_coords2id (0, 0)
 
!         search in direction i
  do
    if (point_i < proc_dims(conf, dimI_offset, proc_id)+proc_dims(conf, dimI_offset+6, proc_id)) then
      exit
    else
!             search for the next proc in i-direction
      proc_id = proc_neighb (proc_id, 1, 1)
      if (proc_id < 0) then
        return
      end if
    end if
  end do
 
!         search in direction j
  do
    if (point_j < proc_dims(conf, dimJ_offset, proc_id)+proc_dims(conf, dimJ_offset+6, proc_id)) then
      exit
    else
!             search for the next proc in j-direction
      proc_id = proc_neighb (proc_id, 1, 2)
      if (proc_id < 0) then
        return
      end if
    end if
  end do
 
  search_proc = proc_id
  return
end function search_proc
 
! --------------------------------------
!
!  get_proc_parts(..)
!
! --------------------------------------
subroutine get_proc_parts (base_x, base_y, base_z, size_x, size_y, size_z, conf, ierr)
  implicit none
 
!         function args
  integer, intent (in) :: base_x, base_y, base_z
  integer, intent (in) :: size_x, size_y, size_z
  integer, intent (in) :: conf
  integer, intent (out) :: ierr
 
!         other vars
  integer :: base_i, base_j, base_k
  integer :: size_i, size_j, size_k
  integer :: dimI_offset, dimJ_offset
  integer :: no_parts, p
 
  integer :: found_proc_id
  integer :: end_of_dim_i, end_of_dim_j
  integer :: size_i_exist
  integer :: start_base_j, start_size_j, start_id_j
 
  ierr = 0
  proc_parts = - 1
 
!    int s,t;
!         pencil in X-dimension (physical space)
  if (conf == 1) then
    base_i = base_y
    base_j = base_z
    base_k = base_x
    size_i = size_y
    size_j = size_z
    size_k = size_x
    dimI_offset = 2
    dimJ_offset = 3
 
!         pencil in Z-dimension (fourier space)
  else if (conf == 2) then
    base_i = base_x
    base_j = base_y
    base_k = base_z
    size_i = size_x
    size_j = size_y
    size_k = size_z
    dimI_offset = 1
    dimJ_offset = 2
 
  else
    ierr = 1
    return
  end if
 
!         find the proc, which contain the base point, this proc will be as start search point
  found_proc_id = search_proc (base_i, base_j, dimI_offset, dimJ_offset, conf)
  if (found_proc_id < 0) then
    ierr = - 1
    return
  end if
 
  no_parts = 0
  end_of_dim_i = 0
  end_of_dim_j = 0
 
!         start to check which proc is contained, write the result in proc_partsbackup in i,j,k system
  do while (end_of_dim_i ==  0)
    no_parts = no_parts + 1
    proc_parts (no_parts, 1) = found_proc_id
 
    start_id_j = found_proc_id
    start_base_j = base_j
    start_size_j = size_j
 
!           dim-i
!           check if it is contained in one proc on i dim
    if (base_i+size_i <= proc_dims(conf, dimI_offset, found_proc_id) &
                        +proc_dims(conf, dimI_offset+6, found_proc_id)) then
      proc_parts (no_parts, 2) = base_i
      proc_parts (no_parts, 5) = size_i
      end_of_dim_i = 1
    else
      proc_parts (no_parts, 2) = base_i
      proc_parts (no_parts, 5) = proc_dims (conf, dimI_offset, found_proc_id) &
                               + proc_dims (conf, dimI_offset+6, found_proc_id) -base_i
!                 compute the new base point and length in the proc
      size_i = base_i + size_i - proc_dims (conf, dimI_offset, found_proc_id) - proc_dims (conf, dimI_offset+6, found_proc_id)
      base_i = proc_dims (conf, dimI_offset, found_proc_id) + proc_dims (conf, dimI_offset+6, found_proc_id)
    end if
 
!           dim-k
    proc_parts (no_parts, 4) = base_k
    proc_parts (no_parts, 7) = size_k
 
!           dim-j
    end_of_dim_j = 0
    size_i_exist = 1
    base_j = start_base_j
    size_j = start_size_j
    do while (end_of_dim_j ==  0)
      if (size_i_exist == 1) then
!               use the old values of i,k stuff
        size_i_exist = 0
      else
        no_parts = no_parts + 1
 
!                the values are copied from last proc data in array
        proc_parts (no_parts, 1) = found_proc_id
        proc_parts (no_parts, 4) = proc_parts (no_parts-1, 4)
        proc_parts (no_parts, 5) = proc_parts (no_parts-1, 5)
        proc_parts (no_parts, 7) = proc_parts (no_parts-1, 7)
      end if
 
      if (base_j+size_j <= proc_dims(conf, dimJ_offset, found_proc_id) &
                         + proc_dims(conf, dimJ_offset+6, found_proc_id)) then
        proc_parts (no_parts, 3) = base_j
        proc_parts (no_parts, 6) = size_j
        end_of_dim_j = 1
      else
        proc_parts (no_parts, 3) = base_j
        proc_parts (no_parts, 6) = proc_dims (conf, dimJ_offset, found_proc_id) &
                                 + proc_dims (conf, dimJ_offset+6, found_proc_id) - base_j
!                     compute the new base point and length in the proc
        size_j = base_j + size_j - proc_dims (conf, dimJ_offset, found_proc_id) - proc_dims (conf, dimJ_offset+6, found_proc_id)
        base_j = proc_dims (conf, dimJ_offset, found_proc_id) + proc_dims (conf, dimJ_offset+6, found_proc_id)
      end if
 
      if (end_of_dim_j == 0) then
        found_proc_id = proc_neighb (found_proc_id, 1, 2)
        if (found_proc_id < 0) exit
      end if
    end do
 
    if (end_of_dim_i == 0) then
      found_proc_id = proc_neighb (start_id_j, 1, 1)
      if (found_proc_id < 0) exit
    end if
 
  end do
 
!         transform the result from i,j,k system to x,y,z system according to the style
  do p = 0, (iproc*jproc)
    base_i = proc_parts (p, 2)
    base_j = proc_parts (p, 3)
    base_k = proc_parts (p, 4)
    size_i = proc_parts (p, 5)
    size_j = proc_parts (p, 6)
    size_k = proc_parts (p, 7)
 
!           physical space
    if (conf == 1) then
      proc_parts (p, 2) = base_k
      proc_parts (p, 3) = base_i
      proc_parts (p, 4) = base_j
      proc_parts (p, 5) = size_k
      proc_parts (p, 6) = size_i
      proc_parts (p, 7) = size_j
 
!           fourier space (no translation nessessary)
    else if (conf .eq. 2) then
!              proc_parts(p,2) = base_i
!              proc_parts(p,3) = base_j
!              proc_parts(p,4) = base_k
!              proc_parts(p,5) = size_i
!              proc_parts(p,6) = size_j
!              proc_parts(p,7) = size_k
    end if
  end do
 
end subroutine get_proc_parts

! --------------------------------------
!
!     Transpose X and Y pencils (real data)
!
! --------------------------------------
subroutine rtran_x2y (source, dest, rbuf1, rbuf2, dstart, dend, dsize, t)
  implicit none
 
  real (mytype), intent (in) :: source (NX_fft, jisize, kjsize)
  real (mytype), intent (out) :: dest (iiisize, NY_fft, kjsize)
  real (mytype), intent (inout) :: rbuf1 (*), rbuf2 (*)
  integer, intent (out) :: dstart (3), dend (3), dsize (3)
 
  real (8) :: t
  integer :: x, y, i, ierr, z, xs, j, N, ix, iy
  integer * 8 :: position
 
! Pack the send buffer for exchanging y and x (within a given z plane ) into sendbuf
  position = 1
 
  do i = 0, iproc - 1
 
    do z = 1, kjsize
      do y = 1, jisize
        do x = iiist (i), iiien (i)
          rbuf1 (position) = source (x, y, z)
          position = position + 1
        end do
      end do
    end do
#ifdef USE_EVEN
    position = (i+1) * IiCntMax / (mytype) + 1
#endif
  end do
 
  t = t - MPI_Wtime ()
#ifdef USE_EVEN
! Use MPI_Alltoall
! Exchange the y-x buffers (in rows of processors)
  call mpi_alltoall (rbuf1, IiCntMax, mpi_byte, rbuf2, IiCntMax, mpi_byte, mpi_comm_row, ierr)
#else
! Use MPI_Alltoallv
! Exchange the y-x buffers (in rows of processors)
  call mpi_alltoallv (rbuf1, IiCnts, IiStrt, mpi_byte, rbuf2, JiCnts, JiStrt, mpi_byte, mpi_comm_row, ierr)
#endif
  t = t + MPI_Wtime ()
 
! Unpack the data
  position = 1
  do i = 0, iproc - 1
    do z = 1, kjsize
      do y = jist (i), jien (i)
        do x = 1, iiisize
          dest (x, y, z) = rbuf2 (position)
          position = position + 1
        end do
      end do
    end do
#ifdef USE_EVEN
    position = (i+1) * IiCntMax / (mytype) + 1
#endif
  end do
 
  dstart (1) = iiistart
  dend (1) = iiiend
  dsize (1) = iiisize
  dstart (2) = 1
  dend (2) = NY_fft
  dsize (2) = NY_fft
  dstart (3) = kjstart
  dend (3) = kjend
  dsize (3) = kjsize
 
  return
end subroutine
 
! --------------------------------------
!
!     Transpose Y and X pencils (real data)
!
! --------------------------------------
subroutine rtran_y2x (source, dest, rbuf1, rbuf2, dstart, dend, dsize, t)
  implicit none
 
  real (mytype), intent (in) :: source (iiisize, NY_fft, kjsize)
  real (mytype), intent (out) :: dest (NX_fft, jisize, kjsize)
  real (mytype), intent (inout) :: rbuf1 (*), rbuf2 (*)
  integer, intent (out) :: dstart (3), dend (3), dsize (3)
 
  real (8) :: t
  integer :: x, y, i, ierr, z, xs, j, N, ix, iy
  integer * 8 :: position, pos1
 
! Pack the send buffer for exchanging y and x (within a given z plane ) into sendbuf
  position = 1
 
  do i = 0, iproc - 1
    do z = 1, kjsize
      do y = jist (i), jien (i)
        do x = 1, iiisize
          rbuf1 (position) = source (x, y, z)
          position = position + 1
        end do
      end do
    end do
#ifdef USE_EVEN
    position = (i+1) * JICntMax / (mytype) + 1
#endif
  end do
 
  t = t - MPI_Wtime ()
#ifdef USE_EVEN
! Use MPI_Alltoall
! Exchange the y-x buffers (in rows of processors)
  call mpi_alltoall (rbuf1, JICntMax, mpi_byte, &
                     rbuf2, JICntMax, mpi_byte, mpi_comm_row, ierr)
#else
! Use MPI_Alltoallv
! Exchange the y-x buffers (in rows of processors)
  call mpi_alltoallv (rbuf1, JiCnts, JiStrt, mpi_byte, &
                      rbuf2, IiCnts, IiStrt, mpi_byte, mpi_comm_row, ierr)
#endif
  t = t + MPI_Wtime ()
 
! Unpack the data
  position = 1
  do i = 0, iproc - 1
    do z = 1, kjsize
      do y = 1, jisize
        do x = iiist (i), iiien (i)
          dest (x, y, z) = rbuf2 (position)
          position = position + 1
        end do
      end do
    end do
#ifdef USE_EVEN
    position = (i+1) * JICntMax / (mytype) + 1
#endif
  end do
 
  dstart (1) = 1
  dend (1) = NX_fft
  dsize (1) = NX_fft
  dstart (2) = jistart
  dend (2) = jiend
  dsize (2) = jisize
  dstart (3) = kjstart
  dend (3) = kjend
  dsize (3) = kjsize
 
  return
end subroutine
 
! --------------------------------------
!
!     Transpose X and Z pencils (real data)
!
! --------------------------------------
subroutine rtran_x2z (source, dest, rbuf1, rbuf2, dstart, dend, dsize, t)
  implicit none
 
  real (mytype), intent (in) :: source (NX_fft, jisize, kjsize)
  real (mytype), intent (out) :: dest (ijsize, jisize, NZ_fft)
  real (mytype), intent (inout) :: rbuf1 (*), rbuf2 (*)
  integer, intent (out) :: dstart (3), dend (3), dsize (3)
 
  real (8) :: t
  integer :: x, y, i, ierr, z, xs, j, N, ix, iy
  integer * 8 :: position
 
! Pack the send buffer for exchanging x and z (within a given y plane ) into sendbuf
  position = 1
 
  do i = 0, jproc - 1
 
    do z = 1, kjsize
      do y = 1, jisize
        do x = ijst (i), ijen (i)
          rbuf1 (position) = source (x, y, z)
          position = position + 1
        end do
      end do
    end do
#ifdef USE_EVEN
    position = (i+1) * IJCntMax / (mytype) + 1
#endif
  end do
 
  t = t - MPI_Wtime ()
#ifdef USE_EVEN
! Use MPI_Alltoall
! Exchange the x-z buffers (in rows of processors)
  call mpi_alltoall (rbuf1, IJCntMax, mpi_byte, &
                     rbuf2, IJCntMax, mpi_byte, mpi_comm_col, ierr)
#else
! Use MPI_Alltoallv
! Exchange the x-z buffers (in rows of processors)
  call mpi_alltoallv (rbuf1, IjCnts, IjStrt, mpi_byte, &
                      rbuf2, KjCnts, KjStrt, mpi_byte, mpi_comm_col, ierr)
#endif
  t = t + MPI_Wtime ()
 
! Unpack the data
  position = 1
  do i = 0, jproc - 1
    do z = kjst (i), kjen (i)
      do y = 1, jisize
        do x = 1, ijsize
          dest (x, y, z) = rbuf2 (position)
          position = position + 1
        end do
      end do
    end do
#ifdef USE_EVEN
    position = (i+1) * IJCntMax / (mytype) + 1
#endif
  end do
 
  dstart (1) = ijstart
  dend (1) = ijend
  dsize (1) = ijsize
  dstart (2) = jistart
  dend (2) = jiend
  dsize (2) = jisize
  dstart (3) = 1
  dend (3) = NZ_fft
  dsize (3) = NZ_fft
 
  return
end subroutine
 
! --------------------------------------
!
!     Transpose Z and X pencils (real data)
!
! --------------------------------------
subroutine rtran_z2x (source, dest, rbuf1, rbuf2, dstart, dend, dsize, t)
  implicit none
 
  real (mytype), intent (in) :: source (ijsize, jisize, NZ_fft)
  real (mytype), intent (out) :: dest (NX_fft, jisize, kjsize)
  real (mytype), intent (inout) :: rbuf1 (*), rbuf2 (*)
  integer, intent (out) :: dstart (3), dend (3), dsize (3)
 
  real (8) :: t
  integer :: x, y, i, ierr, z, xs, j, N, ix, iy
  integer * 8 :: position
 
! Pack the send buffer for exchanging y and x (within a given z plane ) into sendbuf
  position = 1
 
  do i = 0, jproc - 1
 
    do z = kjst (i), kjen (i)
      do y = 1, jisize
        do x = 1, ijsize
          rbuf1 (position) = source (x, y, z)
          position = position + 1
        end do
      end do
    end do
#ifdef USE_EVEN
    position = (i+1) * KjCntMax / (mytype) + 1
#endif
  end do
 
  t = t - MPI_Wtime ()
#ifdef USE_EVEN
! Use MPI_Alltoall
! Exchange the z-x buffers (in rows of processors)
  call mpi_alltoall (rbuf1, KjCntMax, mpi_byte, rbuf2, KjCntMax, mpi_byte, mpi_comm_col, ierr)
#else
! Use MPI_Alltoallv
! Exchange the z-x buffers (in rows of processors)
  call mpi_alltoallv (rbuf1, KjCnts, KjStrt, mpi_byte, rbuf2, IjCnts, IjStrt, mpi_byte, mpi_comm_col, ierr)
#endif
  t = t + MPI_Wtime ()
 
! Unpack the data
  position = 1
  do i = 0, jproc - 1
    do z = 1, kjsize
      do y = 1, jisize
        do x = ijst (i), ijen (i)
          dest (x, y, z) = rbuf2 (position)
          position = position + 1
        end do
      end do
    end do
#ifdef USE_EVEN
    position = (i+1) * KjCntMax / (mytype) + 1
#endif
  end do
 
  dstart (1) = 1
  dend (1) = NX_fft
  dsize (1) = NX_fft
  dstart (2) = jistart
  dend (2) = jiend
  dsize (2) = jisize
  dstart (3) = kjstart
  dend (3) = kjend
  dsize (3) = kjsize
 
  return
end subroutine
 

 

 

 

 


 
!========================================================
 
!      subroutine p3dfft_rot_x180 (XYZg)
!      use fft_spec
!      implicit none
!
!      complex(mytype), TARGET :: XYZg(iistart:iiend,jjstart:jjend,nz_fft)
!
!
!      integer x,y,z,i,k,nx,ny,nz
!
!      if(.not. mpi_set) then
!         print *,'P3DFFT error: call setup before other routines'
!         return
!      endif
!      nx = nx_fft
!      ny = ny_fft
!      nz = nz_fft
!
!     ! flip data in z-direction
!
!     ! z-pencil to y-pencil
!      if(jproc > 1) then
!          call bcomm1(XYZg,XYgZ,timers(3))
!      endif
!
!     ! flip data in y-direction
!      do k=1,nz_fft
!        do j=jjstart,jjend
!            do i=iistart,iiend
!            enddo
!          enddo
!      enddo
!
!     ! y-pencil to z-pencil
!      if(jproc > 1) then
!          call fcomm2(XYgZ,XYZg,timers(2))
!      endif
!
!      end subroutine p3dfft_rot_x180

      end module


