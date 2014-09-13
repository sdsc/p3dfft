! This file is part of P3DFFT library
!
! Version 2.5
!
! Copyright (C) 2011 Jens Henrik Goebbert
!
!    P3DFFT is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!----------------------------------------------------------------------------

!=======================================================================
!                 ghost cell support
!               ----------------------
!> \page ghosts Ghost Cell Support
!!
!! \section use How to use
!!
!!    - 1) init ghost-support:
!!       ==>> call ghosts_init(goverlap)
!!       goverlap = overlap of ghost-cells in all directions (incl. edges and corners)
!!                  (limitations: only one goverlap for all directions supported)
!!       (we are assuming periodic boundary conditions)
!!
!!    - 2) allocation as much memory as value of "gsize_1d" for the data array
!!       ==>> call ghosts_mem(gsize_1d)
!!       - the ghostcells as added at the _end_ of the memory used by p3dfft for fourier transformations
!!       - the min. size: (psize(1)+2*goverlap) *(psize(2)+2*goverlap) *(psize(3)+2*goverlap)
!!           with psize(:) the max. dimension in x,y,z of real or complex
!!                   ==  psize(1) *psize(2) *psize(3)   |xx|yy|zz| = slab 0           <BR>
!!                    faces:                                                          <BR>
!!                      + psize(2)*psize(3)*goverlap    |-x|yy|zz| = slab 1,1 = -X    <BR>
!!                      + psize(2)*psize(3)*goverlap    |+x|yy|zz| = slab 1,2 = +X    <BR>
!!                      + psize(1)*psize(3)*goverlap    |xx|-y|zz| = slab 2,1 = -Y    <BR>
!!                      + psize(1)*psize(3)*goverlap    |xx|+y|zz| = slab 2,2 = +Y    <BR>
!!                      + psize(1)*psize(3)*goverlap    |xx|yy|-z| = slab 3,1 = -Z    <BR>
!!                      + psize(1)*psize(3)*goverlap    |xx|yy|+z| = slab 3,2 = +Z    <BR>
!!                    x-edges on y-planes:                                            <BR>
!!                      + psize(1)*goverlap*goverlap    |xx|-y|-z| = slab 4,1         <BR>
!!                      + psize(1)*goverlap*goverlap    |xx|+y|-z| = slab 4,2         <BR>
!!                      + psize(1)*goverlap*goverlap    |xx|-y|+z| = slab 5,1         <BR>
!!                      + psize(1)*goverlap*goverlap    |xx|+y|+z| = slab 5,2         <BR>
!!                    y-edges on x-planes:                                            <BR>
!!                      + psize(2)*goverlap*goverlap    |-x|yy|-z| = slab 6,1         <BR>
!!                      + psize(2)*goverlap*goverlap    |+x|yy|-z| = slab 6,2         <BR>
!!                      + psize(2)*goverlap*goverlap    |-x|yy|+z| = slab 7,1         <BR>
!!                      + psize(2)*goverlap*goverlap    |+x|yy|+z| = slab 7,2         <BR>
!!                    z-edges on x-planes:                                            <BR>
!!                      + psize(3)*goverlap*goverlap    |-x|-y|zz| = slab 8,1         <BR>
!!                      + psize(3)*goverlap*goverlap    |+x|-y|zz| = slab 8,2         <BR>
!!                      + psize(3)*goverlap*goverlap    |-x|+y|zz| = slab 9,1         <BR>
!!                      + psize(3)*goverlap*goverlap    |+x|+y|zz| = slab 9,2         <BR>
!!                    corners on x-planes                                             <BR>
!!                      + goverlap*goverlap*goverlap    |-x|-y|-z| = slab 10,1        <BR>
!!                      + goverlap*goverlap*goverlap    |+x|-y|-z| = slab 10,2        <BR>
!!                      + goverlap*goverlap*goverlap    |-x|-y|+z| = slab 11,1        <BR>
!!                      + goverlap*goverlap*goverlap    |+x|-y|+z| = slab 11,2        <BR>
!!                      + goverlap*goverlap*goverlap    |-x|+y|-z| = slab 12,1        <BR>
!!                      + goverlap*goverlap*goverlap    |+x|+y|-z| = slab 12,2        <BR>
!!                      + goverlap*goverlap*goverlap    |-x|+y|+z| = slab 13,1        <BR>
!!                      + goverlap*goverlap*goverlap    |+x|+y|+z| = slab 13,2        <BR>
!!
!!     - 3) exchange ghost-cells among processors
!!       to update the ghost-cells among the proceesors
!!       ==>> call ghosts_update(array,glevel,timer).
!!       This will assume data represents physical space!
!!       For best performance update only the ghostcells which are really needed.
!!       Therefore setting glevel has the following effect:
!!         =1  exchange faces
!!         =2  exchange faces and edges
!!         =3  exchange faces,edges and corners
!!
!!     - 4) use/read ghost-cells
!!       if you want to access/read ghost-cells you have functions for different scenarios
!!       which returns the position in a 1D-array as an integer*8
!!        performance = max. 2xif
!!        - search x-faces only(j,k-coordinate does not leave pencil) ==> ghosts_i002i(i,j,k)
!!        - search y-faces only(i,k-coordinate does not leave pencil) ==> ghosts_0j02i(i,j,k)
!!        - search z-faces only(i,j-coordinate does not leave pencil) ==> ghosts_00k2i(i,j,k)
!!        performance = face:5or6xif, edges:5or4xif, corners:4or3xif
!!        - search in this order corner,edges,faces ==> ghosts_ijk2i(i,j,k)
!!       if you do not want to handle 1d-arrays the following functions
!!       return the value of the array at the coordinate (just pass the 3d-array)
!!        performance = max. 2xif
!!        - search x-faces only(j,k-coordinate does not leave pencil) ==> ghosts_i002(array,i,j,k)
!!        - search y-faces only(i,k-coordinate does not leave pencil) ==> ghosts_0j02(array,i,j,k)
!!        - search z-faces only(i,j-coordinate does not leave pencil) ==> ghosts_00k2(array,i,j,k)
!!        performance = face:5or6xif, edges:5or4xif, corners:4or3xif
!!        - search in this order corner,edges,faces ==> ghosts_ijk(array,i,j,k)
!!
!!
!!  \section mem array memory:
!! ==============
!!     - 1) size if complex-type: (psize(1)+2*goverlap)
!!                             *(psize(2)+2*goverlap)
!!                             *(psize(3)+2*goverlap)
!!     - 2a) slabs in memory:
!!      |000000000000000000000000000|11111|22222|33333|44444|55555|66666|77777| ....
!!                 ^^^^              ^     ^     ^     ^     ^     ^
!!          input/output of FFT      A     B     C     D     E     F
!!
!!      ghost-cell-slab at each pencil-side:
!!      A:  -X,Y*Z	A(istart(1)-i, j,k)		=> gmem_start(1,1)
!!      B:  +X,Y*Z	B(i,j,k)				=> gmem_start(1,2)
!!      C:  -Y,X*Z	C(i,j,k)				=> gmem_start(2,1)
!!      D:  +Y,X*Z	D(i,j,k)				=> gmem_start(2,2)
!!      E:  -Z,X*Y	E(i,j,k)				=> gmem_start(3,1)
!!      F:  +Z,X*Y	F(i,j,k)				=> gmem_start(3,2)
!!      G:  ...etc.
!!
!!  \section vars important variables:
!! =====================
!!     - 1) goverlap
!!       the data overlap between two processors - equal to the thickness
!!       of a ghost-slab
!!       This version supports only one fixed overlap value, because it simplifies
!!       the coding. If different overlaps for different arrays are needed the
!!       following varables must be dependent on the overlap-value and therefor
!!       should have one additional dimension.
!!       A good idea would be to make p3dfft_init_ghosts(..) return an id,which
!!       has to get set when calling update_ghosts(data,id).
!!
!!     - 2) glevel
!!       =1 : send face-ghostcells only
!!       =2 : send face- and edge-ghostcells
!!       =3 : send face-,edge- and corner-ghostcells
!!
!!     - 3) gmem_rstart(1..13,1..2) (13 ghost-slab pairs: faces, edges, corners)
!!       start of data in memory for ghost-cells of certain ghost-slab
!!
!!     - 4) gslab_rstart(1..13, 1..2)
!!       relative start coord. of ghost-slab in pencil ijk(1..3) for faces and
!!       (0,0,0) for edges and corners
!!
!!     - 5) gslab_size(1..13, 1..2)
!!       size of ghost-slab for faces,edges,corners
!!
!!     - 6) gneighb_r(mpi_tasks*(3*2))
!!       extended processor grid with all processor neighbours in -/+XYZ(1..6)
!!       including periodic boundary neighbours
!!        example for 2 processors with gneighb_r(:) = 0,0,0,0,1,1,1,1,1,1,0,0
!!         pencil-side: -X  +X  -Y  +Y  -Z  +Z
!!         taskid 0:     0   0   0   0   1   1
!!         taskid 1:     1   1   1   1   0   0
!!
!!  \section ghostdeps dependencies to p3dfft:
!! ========================
!!   The following informations are borrowed from p3dfft...
!!     - p3dfft_get_mpi_info(mpi_taskid, mpi_tasks, mpi_comm)
!!     - p3dfft_get_dims(gproc_rstart, gproc_rend, gproc_rsize, 1)
!!     - p3dfft_get_dims(mstart, mend, msize, 3)
!!     - proc_id2coords
!!     - proc_coords2id
!!
!=======================================================================

      module p3dfft_ghosts
      use p3dfft
      implicit none

      include 'mpif.h'

      private

! ghost-cell support using p3dfft_init_ghosts(), p3dfft_update_ghosts()
      integer,save :: mpi_taskid
      integer,save :: mpi_tasks
      integer,save :: mpi_comm
      integer,save :: mpi_mytype

      integer,save :: goverlap
      logical :: ghosts_set = .false.

      integer,save :: gproc_rsize(3), gproc_rstart(3), gproc_rend(3)

      integer,save, dimension(:), allocatable :: gneighb_r
      integer,save, dimension(:), allocatable :: gneighb_c

      integer(kind=1),save, dimension(:), allocatable :: gslabbuf_snd
      integer(kind=1),save, dimension(:), allocatable :: gslabbuf_recv
      integer,save :: gslabbuf_bytesize

      integer,save :: gslab_sndnghb (13,2)
      integer,save :: gslab_recvnghb(13,2)

      integer,save :: gslab_rsize (3,13)

      integer,save :: mpi_mytype_gslab_rrecv(13,2)
      integer,save :: mpi_mytype_gslab_rsnd (13,2)

      integer(kind=8),save :: gmem_rstart(13,2)
      integer(kind=8),save :: gmem_rsndstart(13,2)

      public :: ghosts_init
      public :: ghosts_mem
      public :: ghosts_update
      public :: ghosts_check
      public :: ghosts_ijk2i, ghosts_i002i, ghosts_0j02i, ghosts_00k2i
      public :: ghosts_ijk,   ghosts_i00,   ghosts_0j0,   ghosts_00k
      public :: ghosts_ijk_debug

!-------------------
      contains
!-------------------

!    !========================================
!    !
!    !		init ghosts
!    !
!    !========================================
      subroutine ghosts_init(goverlap_in)
      	implicit none

!     	function args
      	integer, intent(in) :: goverlap_in

!     	other vars
      	integer, dimension(:,:), allocatable :: coords2id
      	integer :: iproc,jproc
      	integer :: pneighb_r(6)
      	integer :: pCoo(2)
      	integer :: msize(3), mstart(3), mend(3)
      	integer :: i,p,d,mp,mp1,ierr

        integer :: gslabbuf_rsndsize (13), gslabbuf_rsndsize_tmp, gslabbuf_rsndsize_max
        integer :: gslabbuf_rrecvsize(13), gslabbuf_rrecvsize_tmp, gslabbuf_rrecvsize_max

      	integer(kind=8) :: gmess_rsize(7)
      	integer :: gslab_rparentsize(3,13)
      	integer :: gslab_rsndstart(3,13,2)

      	integer :: start000(3)

      	goverlap = goverlap_in

      	call get_mpi_info(mpi_taskid, mpi_tasks, mpi_comm)

      	call get_dims(gproc_rstart, gproc_rend, gproc_rsize, 1)
      	call get_dims(mstart, mend, msize, 3)

      	if(     goverlap .gt. gproc_rsize(1)			 &
      	   .or. goverlap .gt. gproc_rsize(2)			 &
      	   .or. goverlap .gt. gproc_rsize(3) ) then
      			print *,'P3DFFT error: in one direction goverlap is greater then pencil size'
      			ierr = 1
      			call abort
      	endif

!     !	build help-array with coordinates to proc-ids
      	iproc = size(proc_coords2id,1)
      	jproc = size(proc_coords2id,2)

      	allocate(coords2id(-1:iproc,-1:jproc))
      	coords2id(:,:) = -1
      	do i=0,mpi_tasks-1
      		coords2id( proc_id2coords(2*i), proc_id2coords(2*i+1)) = i
      	enddo

!     !	extend array with one line more on each side (periodicy)
      	coords2id(-1,   0:jproc-1) = coords2id(iproc-1,0:jproc-1)
      	coords2id(iproc,0:jproc-1) = coords2id(0,      0:jproc-1)
      	coords2id(0:iproc-1,   -1) = coords2id(0:iproc-1,jproc-1)
      	coords2id(0:iproc-1,jproc) = coords2id(0:iproc-1,      0)

!      	! build pneighb_r
      	pCoo(1) = proc_id2coords(mpi_taskid*2)
      	pCoo(2) = proc_id2coords(mpi_taskid*2+1)
      	do d=0, 2	! loop over X,Y,Z sides
      	  do mp=1,2		! loop over direction (-/+)
      		mp1 = mp*2-3	! 1=>-1 , 2=>+1

!      		! build pneighb_r
      		if(	d .eq. 0 .or.							 & ! no decomp. in X
      			d .eq. 1 .and. iproc .eq. 1 .or.		 & ! no decomp. in Y
      			d .eq. 2 .and. jproc .eq. 1) then		  ! no decomp. in Z
      				pneighb_r(d*2+mp) = mpi_taskid
      		else if(d .eq. 1) then					! decomp. in Y
      			pneighb_r(d*2+mp) = coords2id(pCoo(1)+mp1, pCoo(2))
      		else if(d .eq. 2) then					! decomp. in Z
      			pneighb_r(d*2+mp) = coords2id(pCoo(1), pCoo(2)+mp1)
      		endif

      	  enddo
      	enddo

      	deallocate(coords2id)

!     ! distribute pneigh_r
      	allocate(gneighb_r(mpi_tasks*6))
      	call MPI_Allgather(	pneighb_r,	6, MPI_INTEGER,	 &
      						gneighb_r,	6, MPI_INTEGER,	 &
      						mpi_comm, ierr)

!     	set mpi message sizes for exchange of ghost-cells between processors
      	gmess_rsize(1) =  goverlap       *gproc_rsize(2) *gproc_rsize(3)		! |-x|yy|zz|, |+x|yy|zz|
      	gmess_rsize(2) =  gproc_rsize(1) *goverlap       *gproc_rsize(3)		! |xx|-y|zz|, |xx|+y|zz|
      	gmess_rsize(3) =  gproc_rsize(1) *gproc_rsize(2) *goverlap				! |xx|yy|-z|, |xx|yy|+z|
      	gmess_rsize(4) =  gproc_rsize(1) *goverlap *goverlap					! |xx|-y|-z|, |xx|+y|-z|, |xx|-y|+z|,  |xx|+y|+z|
      	gmess_rsize(5) =  goverlap *gproc_rsize(2) *goverlap					! |-x|yy|-z|, |+x|yy|-z|, |-x|yy|+z|,  |+x|yy|+z|
      	gmess_rsize(6) =  goverlap *goverlap *gproc_rsize(3)					! |-x|-y|zz|, |+x|-y|zz|, |-x|+y|zz|,  |+x|+y|zz|
      	gmess_rsize(7) =  goverlap *goverlap *goverlap							! |-x|-y|-z|, |+x|-y|-z|, |-x|-y|+z|, |+x|-y|+z|, |-x|+y|-z|, |+x|+y|-z|, |-x|+y|+z|, |+x|+y|+z|

!     	precalc ghost-slab start-position in memory of in/out array
!     	faces
      	gmem_rstart(1,1)  = int(msize(1),8) *int(msize(2),8) *int(msize(3),8) +1	! |-x|yy|zz|
      	gmem_rstart(1,2)  = gmem_rstart(1,1)  +gmess_rsize(1)						! |+x|yy|zz|
     	gmem_rstart(2,1)  = gmem_rstart(1,2)  +gmess_rsize(1)						! |xx|-y|zz|
      	gmem_rstart(2,2)  = gmem_rstart(2,1)  +gmess_rsize(2)						! |xx|+y|zz|
      	gmem_rstart(3,1)  = gmem_rstart(2,2)  +gmess_rsize(2)						! |xx|yy|-z|
      	gmem_rstart(3,2)  = gmem_rstart(3,1)  +gmess_rsize(3)						! |xx|yy|+z|
!     	x-edges swapped in y
      	gmem_rstart(4,1)  = gmem_rstart(3,2)  +gmess_rsize(3)						! |xx|-y|-z|
      	gmem_rstart(4,2)  = gmem_rstart(4,1)  +gmess_rsize(4)						! |xx|+y|-z|
      	gmem_rstart(5,1)  = gmem_rstart(4,2)  +gmess_rsize(4)						! |xx|-y|+z|
      	gmem_rstart(5,2)  = gmem_rstart(5,1)  +gmess_rsize(4)						! |xx|+y|+z|
!     	y-edges swapped in x
      	gmem_rstart(6,1)  = gmem_rstart(5,2)  +gmess_rsize(4)						! |-x|yy|-z|
      	gmem_rstart(6,2)  = gmem_rstart(6,1)  +gmess_rsize(5)						! |+x|yy|-z|
      	gmem_rstart(7,1)  = gmem_rstart(6,2)  +gmess_rsize(5)						! |-x|yy|+z|
      	gmem_rstart(7,2)  = gmem_rstart(7,1)  +gmess_rsize(5)						! |+x|yy|+z|
!     	z-edges swapped in x
      	gmem_rstart(8,1)  = gmem_rstart(7,2)  +gmess_rsize(5)						! |-x|-y|zz|
      	gmem_rstart(8,2)  = gmem_rstart(8,1)  +gmess_rsize(6)						! |+x|-y|zz|
      	gmem_rstart(9,1)  = gmem_rstart(8,2)  +gmess_rsize(6)						! |-x|+y|zz|
      	gmem_rstart(9,2)  = gmem_rstart(9,1)  +gmess_rsize(6)						! |+x|+y|zz|
!     	corners swapped in x
      	gmem_rstart(10,1) = gmem_rstart(9,2)  +gmess_rsize(6)						! |-x|-y|-z|
      	gmem_rstart(10,2) = gmem_rstart(10,1) +gmess_rsize(7)						! |+x|-y|-z|
      	gmem_rstart(11,1) = gmem_rstart(10,2) +gmess_rsize(7)						! |-x|-y|+z|
      	gmem_rstart(11,2) = gmem_rstart(11,1) +gmess_rsize(7)						! |+x|-y|+z|
      	gmem_rstart(12,1) = gmem_rstart(11,2) +gmess_rsize(7)						! |-x|+y|-z|
      	gmem_rstart(12,2) = gmem_rstart(12,1) +gmess_rsize(7)						! |+x|+y|-z|
      	gmem_rstart(13,1) = gmem_rstart(12,2) +gmess_rsize(7)						! |-x|+y|+z|
      	gmem_rstart(13,2) = gmem_rstart(13,1) +gmess_rsize(7)						! |+x|+y|+z|

!     	precalc ghost-slab start-position in memory of in/out array
!     	faces
      	gmem_rsndstart(1,1)  = 1										! |-x|yy|zz|
      	gmem_rsndstart(1,2)  = 1										! |+x|yy|zz|
     	gmem_rsndstart(2,1)  = 1										! |xx|-y|zz|
      	gmem_rsndstart(2,2)  = 1										! |xx|+y|zz|
      	gmem_rsndstart(3,1)  = 1										! |xx|yy|-z|
      	gmem_rsndstart(3,2)  = 1										! |xx|yy|+z|
!     	x-edges swapped in y
      	gmem_rsndstart(4,1)  = gmem_rstart(3,1)							! |xx|-y|-z|
      	gmem_rsndstart(4,2)  = gmem_rstart(3,1)							! |xx|+y|-z|
      	gmem_rsndstart(5,1)  = gmem_rstart(3,2)							! |xx|-y|+z|
      	gmem_rsndstart(5,2)  = gmem_rstart(3,2)							! |xx|+y|+z|
!     	y-edges swapped in x
      	gmem_rsndstart(6,1)  = gmem_rstart(3,1)							! |-x|yy|-z|
      	gmem_rsndstart(6,2)  = gmem_rstart(3,1)							! |+x|yy|-z|
      	gmem_rsndstart(7,1)  = gmem_rstart(3,2)							! |-x|yy|+z|
      	gmem_rsndstart(7,2)  = gmem_rstart(3,2)							! |+x|yy|+z|
!     	z-edges swapped in x
      	gmem_rsndstart(8,1)  = gmem_rstart(2,1)							! |-x|-y|zz|
      	gmem_rsndstart(8,2)  = gmem_rstart(2,1)							! |+x|-y|zz|
      	gmem_rsndstart(9,1)  = gmem_rstart(2,2)							! |-x|+y|zz|
      	gmem_rsndstart(9,2)  = gmem_rstart(2,2)							! |+x|+y|zz|
!     	corners swapped in x
      	gmem_rsndstart(10,1) = gmem_rstart(4,1)							! |-x|-y|-z|
      	gmem_rsndstart(10,2) = gmem_rstart(4,1)							! |+x|-y|-z|
      	gmem_rsndstart(11,1) = gmem_rstart(5,1)							! |-x|-y|+z|
      	gmem_rsndstart(11,2) = gmem_rstart(5,1)							! |+x|-y|+z|
      	gmem_rsndstart(12,1) = gmem_rstart(4,2)							! |-x|+y|-z|
      	gmem_rsndstart(12,2) = gmem_rstart(4,2)							! |+x|+y|-z|
      	gmem_rsndstart(13,1) = gmem_rstart(5,2)							! |-x|+y|+z|
      	gmem_rsndstart(13,2) = gmem_rstart(5,2)							! |+x|+y|+z|

!     	precalc ghost-slab size
!     	faces
      	gslab_rsize(:,1)  = (/goverlap,       gproc_rsize(2), gproc_rsize(3)/)	! |-x|yy|zz|, |+x|yy|zz|
      	gslab_rsize(:,2)  = (/gproc_rsize(1), goverlap,       gproc_rsize(3)/)	! |xx|-y|zz|, |xx|+y|zz|
      	gslab_rsize(:,3)  = (/gproc_rsize(1), gproc_rsize(2), goverlap/)		! |xx|yy|-z|, |xx|yy|+z|
!     	x-edges swapped in y
      	gslab_rsize(:,4)  = (/gproc_rsize(1), goverlap,       goverlap      /)	! |xx|-y|-z|, |xx|+y|-z|
      	gslab_rsize(:,5)  = (/gproc_rsize(1), goverlap,       goverlap      /)	! |xx|-y|+z|, |xx|+y|+z|
!     	y-edges swapped in x
      	gslab_rsize(:,6)  = (/goverlap,       gproc_rsize(2), goverlap      /)	! |-x|yy|-z|, |+x|yy|-z|
      	gslab_rsize(:,7)  = (/goverlap,       gproc_rsize(2), goverlap      /)	! |-x|yy|+z|, |+x|yy|+z|
!     	z-edges swapped in x
      	gslab_rsize(:,8)  = (/goverlap,       goverlap,       gproc_rsize(3)/)	! |-x|-y|zz|, |+x|-y|zz|
      	gslab_rsize(:,9)  = (/goverlap,       goverlap,       gproc_rsize(3)/)	! |-x|+y|zz|, |+x|+y|zz|
!     	corners swapped in x
      	gslab_rsize(:,10) = (/goverlap,       goverlap,       goverlap      /)	! |-x|-y|-z|, |+x|-y|-z|
      	gslab_rsize(:,11) = (/goverlap,       goverlap,       goverlap      /)	! |-x|-y|+z|, |+x|-y|+z|
      	gslab_rsize(:,12) = (/goverlap,       goverlap,       goverlap      /)	! |-x|+y|-z|, |+x|+y|-z|
      	gslab_rsize(:,13) = (/goverlap,       goverlap,       goverlap      /)	! |-x|+y|+z|, |+x|+y|+z|

!     	precalc size of array the ghost-slab gets copied from (its parent)
!     	faces
      	gslab_rparentsize(:,1)  = gproc_rsize							! |-x|yy|zz|, |+x|yy|zz|
      	gslab_rparentsize(:,2)  = gproc_rsize							! |xx|-y|zz|, |xx|+y|zz|
      	gslab_rparentsize(:,3)  = gproc_rsize							! |xx|yy|-z|, |xx|yy|+z|
!     	x-edges swapped in y
      	gslab_rparentsize(:,4)  = gslab_rsize(:,3)						! |xx|-y|-z|, |xx|+y|-z|
      	gslab_rparentsize(:,5)  = gslab_rsize(:,3)						! |xx|-y|+z|, |xx|+y|+z|
!     	y-edges swapped in x
      	gslab_rparentsize(:,6)  = gslab_rsize(:,3)						! |-x|yy|-z|, |+x|yy|-z|
      	gslab_rparentsize(:,7)  = gslab_rsize(:,3)						! |-x|yy|+z|, |+x|yy|+z|
!     	z-edges swapped in x
      	gslab_rparentsize(:,8)  = gslab_rsize(:,2)						! |-x|-y|zz|, |+x|-y|zz|
      	gslab_rparentsize(:,9)  = gslab_rsize(:,2)						! |-x|+y|zz|, |+x|+y|zz|
!     	corners swapped in x
      	gslab_rparentsize(:,10) = gslab_rsize(:,4)						! |-x|-y|-z|, |+x|-y|-z|
      	gslab_rparentsize(:,11) = gslab_rsize(:,5)						! |-x|-y|+z|, |+x|-y|+z|
      	gslab_rparentsize(:,12) = gslab_rsize(:,4)						! |-x|+y|-z|, |+x|+y|-z|
      	gslab_rparentsize(:,13) = gslab_rsize(:,5)						! |-x|+y|+z|, |+x|+y|+z|

!     	precalc ghost-slab start in memory
!     	faces
      	gslab_rsndstart(:,1,1)  = (/0, 0, 0/)							! snd to |-x|yy|zz|
      	gslab_rsndstart(:,1,2)  = (/gproc_rsize(1)-gslab_rsize(1,1),	 &
      							    gproc_rsize(2)-gslab_rsize(2,1),	 &
      							    gproc_rsize(3)-gslab_rsize(3,1)/)	! snd to |+x|yy|zz|
      	gslab_rsndstart(:,2,1)  = (/0, 0, 0/)							! snd to |xx|-y|zz|
      	gslab_rsndstart(:,2,2)  = (/gproc_rsize(1)-gslab_rsize(1,2),	 &
      							    gproc_rsize(2)-gslab_rsize(2,2),	 &
      							    gproc_rsize(3)-gslab_rsize(3,2)/)	! snd to |xx|+y|zz|
      	gslab_rsndstart(:,3,1)  = (/0, 0, 0/)							! snd to |xx|yy|-z|
      	gslab_rsndstart(:,3,2)  = (/gproc_rsize(1)-gslab_rsize(1,3),	 &
      							    gproc_rsize(2)-gslab_rsize(2,3),	 &
      							    gproc_rsize(3)-gslab_rsize(3,3)/)	! snd to |xx|yy|+z|
!     	x-edges swapped in y
      	gslab_rsndstart(:,4,1)  = (/0, 0, 0/)							! snd to |xx|-y|-z| from 3,1
      	gslab_rsndstart(:,4,2)  = (/0, gproc_rsize(2)-goverlap, 0/)		! snd to |xx|+y|-z| from 3,1
      	gslab_rsndstart(:,5,1)  = (/0, 0, 0/)							! snd to |xx|-y|+z| from 3,2
      	gslab_rsndstart(:,5,2)  = (/0, gproc_rsize(2)-goverlap, 0/)		! snd to |xx|+y|+z| from 3,2
!     	y-edges swapped in x
      	gslab_rsndstart(:,6,1)  = (/0, 0, 0/)							! snd to |-x|yy|-z| from 3,1
      	gslab_rsndstart(:,6,2)  = (/gproc_rsize(1)-goverlap, 0, 0/)		! snd to |+x|yy|-z| from 3,1
      	gslab_rsndstart(:,7,1)  = (/0, 0, 0/)							! snd to |-x|yy|+z| from 3,2
      	gslab_rsndstart(:,7,2)  = (/gproc_rsize(1)-goverlap, 0, 0/)		! snd to |+x|yy|+z| from 3,2
!     	z-edges swapped in x
      	gslab_rsndstart(:,8,1)  = (/0, 0, 0/)							! snd to |-x|-y|zz| from 2,1
      	gslab_rsndstart(:,8,2)  = (/gproc_rsize(1)-goverlap, 0, 0/)		! snd to |+x|-y|zz| from 2,1
      	gslab_rsndstart(:,9,1)  = (/0, 0, 0/)							! snd to |-x|+y|zz| from 2,2
      	gslab_rsndstart(:,9,2)  = (/gproc_rsize(1)-goverlap, 0, 0/)		! snd to |+x|+y|zz| from 2,2
!     	corners swapped in x
      	gslab_rsndstart(:,10,1) = (/0, 0, 0/)							! |-x|-y|-z|
      	gslab_rsndstart(:,10,2) = (/gproc_rsize(1)-goverlap, 0, 0/)		! |+x|-y|-z|
      	gslab_rsndstart(:,11,1) = (/0, 0, 0/)							! |-x|-y|+z|
      	gslab_rsndstart(:,11,2) = (/gproc_rsize(1)-goverlap, 0, 0/)		! |+x|-y|+z|
      	gslab_rsndstart(:,12,1) = (/0, 0, 0/)							! |-x|+y|-z|
      	gslab_rsndstart(:,12,2) = (/gproc_rsize(1)-goverlap, 0, 0/)		! |+x|+y|-z|
      	gslab_rsndstart(:,13,1) = (/0, 0, 0/)							! |-x|+y|+z|
      	gslab_rsndstart(:,13,2) = (/gproc_rsize(1)-goverlap, 0, 0/)		! |+x|+y|+z|

      	call MPI_Type_contiguous(mytype, MPI_BYTE, mpi_mytype, ierr)
      	call MPI_Type_Commit(mpi_mytype, ierr)

        mpi_mytype_gslab_rsnd(:,:) = 0
        mpi_mytype_gslab_rrecv(:,:) = 0
      	do i=1, 13
      	  do p=1,2
      		call MPI_Type_Create_Subarray( 3, gslab_rparentsize(:,i), gslab_rsize(:,i), gslab_rsndstart(:,i,p),					 &
      									MPI_ORDER_FORTRAN,																		 &
      									mpi_mytype,																				 &
      									mpi_mytype_gslab_rsnd(i,p),																 &
      									ierr )
      		call MPI_Type_Commit(mpi_mytype_gslab_rsnd(i,p), ierr)
      		start000= (/0,0,0/)
      		call MPI_Type_Create_Subarray( 3, gslab_rsize(:,i), gslab_rsize(:,i), start000,										 &
      									MPI_ORDER_FORTRAN,																		 &
      									mpi_mytype,																				 &
      									mpi_mytype_gslab_rrecv(i,p),															 &
      									ierr )
      		call MPI_Type_Commit(mpi_mytype_gslab_rrecv(i,p), ierr)
      	  enddo
      	enddo

!     !	calc send/recieve buffer in bytes and allocate
        gslabbuf_rsndsize_max = 0
        gslabbuf_rrecvsize_max = 0
      	do i=1, 13
      	  do p=1,2
      		call MPI_Pack_Size(1, mpi_mytype_gslab_rsnd(i,p), mpi_comm, gslabbuf_rsndsize(i), ierr)
      		call MPI_Pack_Size(1, mpi_mytype_gslab_rrecv(i,p), mpi_comm, gslabbuf_rrecvsize(i), ierr)
      	  enddo
      	enddo

        ! buffer size for (glevel > 1)
        gslabbuf_rsndsize_tmp = 0
        gslabbuf_rrecvsize_tmp = 0
      	do i=4,5
      	    gslabbuf_rsndsize_tmp  = gslabbuf_rsndsize_tmp  +gslabbuf_rsndsize(i)
      	    gslabbuf_rrecvsize_tmp = gslabbuf_rrecvsize_tmp +gslabbuf_rrecvsize(i)
      	end do
        if(gslabbuf_rsndsize_tmp  > gslabbuf_rsndsize_max) gslabbuf_rsndsize_max  = gslabbuf_rsndsize_tmp
        if(gslabbuf_rrecvsize_tmp > gslabbuf_rrecvsize_max)gslabbuf_rrecvsize_max = gslabbuf_rrecvsize_tmp

        ! buffer size for (glevel > 2)
        gslabbuf_rsndsize_tmp = 0
        gslabbuf_rrecvsize_tmp = 0
        do i=6,9
            gslabbuf_rsndsize_tmp  = gslabbuf_rsndsize_tmp  +gslabbuf_rsndsize(i)
            gslabbuf_rrecvsize_tmp = gslabbuf_rrecvsize_tmp +gslabbuf_rrecvsize(i)
        end do
        do i=10,13
            gslabbuf_rsndsize_tmp  = gslabbuf_rsndsize_tmp  +gslabbuf_rsndsize(i)
            gslabbuf_rrecvsize_tmp = gslabbuf_rrecvsize_tmp +gslabbuf_rrecvsize(i)
        end do
        if(gslabbuf_rsndsize_tmp  > gslabbuf_rsndsize_max) gslabbuf_rsndsize_max  = gslabbuf_rsndsize_tmp
        if(gslabbuf_rrecvsize_tmp > gslabbuf_rrecvsize_max)gslabbuf_rrecvsize_max = gslabbuf_rrecvsize_tmp

        gslabbuf_bytesize = max(gslabbuf_rsndsize_max,gslabbuf_rrecvsize_max)

      	if(allocated(gslabbuf_snd)) deallocate(gslabbuf_snd)
      	if(allocated(gslabbuf_recv)) deallocate(gslabbuf_recv)
      	allocate(gslabbuf_snd(gslabbuf_bytesize),  stat=ierr)
      	allocate(gslabbuf_recv(gslabbuf_bytesize), stat=ierr)

!     	precalc snd neighb (send_neighb = (taskid*6) +(d-1)*2 +1|2))
!     	faces
      	gslab_sndnghb(1,1)  = 1	! |-x|yy|zz|
      	gslab_sndnghb(1,2)  = 1	! |+x|yy|zz|
      	gslab_sndnghb(2,1)  = 2	! |xx|-y|zz|
      	gslab_sndnghb(2,2)  = 2	! |xx|+y|zz|
      	gslab_sndnghb(3,1)  = 3	! |xx|yy|-z|
      	gslab_sndnghb(3,2)  = 3	! |xx|yy|+z|
!     	x-edges swapped in y
      	gslab_sndnghb(4,1)  = 2	! |xx|-y|-z|
      	gslab_sndnghb(4,2)  = 2	! |xx|+y|-z|
      	gslab_sndnghb(5,1)  = 2	! |xx|-y|+z|
      	gslab_sndnghb(5,2)  = 2	! |xx|+y|+z|
!     	y-edges swapped in x
      	gslab_sndnghb(6,1)  = 1	! |-x|yy|-z|
      	gslab_sndnghb(6,2)  = 1	! |+x|yy|-z|
      	gslab_sndnghb(7,1)  = 1	! |-x|yy|+z|
      	gslab_sndnghb(7,2)  = 1	! |+x|yy|+z|
!     	z-edges swapped in x
      	gslab_sndnghb(8,1)  = 1	! |-x|-y|zz|
      	gslab_sndnghb(8,2)  = 1	! |+x|-y|zz|
      	gslab_sndnghb(9,1)  = 1	! |-x|+y|zz|
      	gslab_sndnghb(9,2)  = 1	! |+x|+y|zz|
!     	corners swapped in x
      	gslab_sndnghb(10,1) = 1	! |-x|-y|-z|
      	gslab_sndnghb(10,2) = 1	! |+x|-y|-z|
      	gslab_sndnghb(11,1) = 1	! |-x|-y|+z|
      	gslab_sndnghb(11,2) = 1	! |+x|-y|+z|
      	gslab_sndnghb(12,1) = 1	! |-x|+y|-z|
      	gslab_sndnghb(12,2) = 1	! |+x|+y|-z|
      	gslab_sndnghb(13,1) = 1	! |-x|+y|+z|
      	gslab_sndnghb(13,2) = 1	! |+x|+y|+z|
      	do i=1, 13
      	  gslab_sndnghb(i,1) = (mpi_taskid*6) +(gslab_sndnghb(i,1)-1)*2 +1
      	  gslab_sndnghb(i,2) = (mpi_taskid*6) +(gslab_sndnghb(i,2)-1)*2 +2
      	enddo

!     	precalc recv neighb (recv_neighb = (taskid*6) +(d-1)*2 +1|2)
!     	faces
      	gslab_recvnghb(1,1)  = 1	! |-x|yy|zz|
      	gslab_recvnghb(1,2)  = 1	! |+x|yy|zz|
      	gslab_recvnghb(2,1)  = 2	! |xx|-y|zz|
      	gslab_recvnghb(2,2)  = 2	! |xx|+y|zz|
      	gslab_recvnghb(3,1)  = 3	! |xx|yy|-z|
      	gslab_recvnghb(3,2)  = 3	! |xx|yy|+z|
!     	x-edges swapped in y
      	gslab_recvnghb(4,1)  = 2	! |xx|-y|-z|
      	gslab_recvnghb(4,2)  = 2	! |xx|+y|-z|
      	gslab_recvnghb(5,1)  = 2	! |xx|-y|+z|
      	gslab_recvnghb(5,2)  = 2	! |xx|+y|+z|
!     	y-edges swapped in x
      	gslab_recvnghb(6,1)  = 1	! |-x|yy|-z|
      	gslab_recvnghb(6,2)  = 1	! |+x|yy|-z|
      	gslab_recvnghb(7,1)  = 1	! |-x|yy|+z|
      	gslab_recvnghb(7,2)  = 1	! |+x|yy|+z|
!     	z-edges swapped in x
      	gslab_recvnghb(8,1)  = 1	! |-x|-y|zz|
      	gslab_recvnghb(8,2)  = 1	! |+x|-y|zz|
      	gslab_recvnghb(9,1)  = 1	! |-x|+y|zz|
      	gslab_recvnghb(9,2)  = 1	! |+x|+y|zz|
!     	corners swapped in x
      	gslab_recvnghb(10,1) = 1	! |-x|-y|-z|
      	gslab_recvnghb(10,2) = 1	! |+x|-y|-z|
      	gslab_recvnghb(11,1) = 1	! |-x|-y|+z|
      	gslab_recvnghb(11,2) = 1	! |+x|-y|+z|
      	gslab_recvnghb(12,1) = 1	! |-x|+y|-z|
      	gslab_recvnghb(12,2) = 1	! |+x|+y|-z|
      	gslab_recvnghb(13,1) = 1	! |-x|+y|+z|
      	gslab_recvnghb(13,2) = 1	! |+x|+y|+z|
      	do i=1, 13
      	  gslab_recvnghb(i,1) = (mpi_taskid*6) +(gslab_recvnghb(i,1)-1)*2 +1
      	  gslab_recvnghb(i,2) = (mpi_taskid*6) +(gslab_recvnghb(i,2)-1)*2 +2
      	enddo

      	ghosts_set = .true.

      end subroutine ghosts_init

!    !========================================
!    !
!    !		memory of ghost-pencils
!    !
!    !========================================
      subroutine ghosts_mem(gsize_mem)
      	implicit none

!     	function args
      	integer(kind=8), intent(out) :: gsize_mem

!     	other vars
      	integer :: msize(3), mstart(3), mend(3)

      	if(.not. ghosts_set) then
      		print *,'P3DFFT error: call ghosts_init() before'
      		stop
      	endif

      	call get_dims(mstart, mend, msize, 3)
      	gsize_mem = int(msize(1)+2*goverlap,8)	 &
      	           *int(msize(2)+2*goverlap,8)	 &
      	           *int(msize(3)+2*goverlap,8)

      end subroutine ghosts_mem

!    !========================================
!    !  ghosts_update(XYZ, goverlap)
!    !     # XYZ - array of real made out of two parts
!    !              1. part = input data = FFT-pencil
!    !                        [(fsize(1)*2) *(fsize(2)*2) *(fsize(3)*2)]
!    !              2. part = output data = ghost-cells
!    !                         [(fsize(1)*2) *(fsize(2)*2) *goverlap]*2	2* XY-plane
!    !                        +[(fsize(1)*2) *(fsize(3)*2) *goverlap]*2	2* XZ-plane
!    !                        +[(fsize(2)*2) *(fsize(3)*2) *goverlap]*2	2* YZ-plane
!    !     # gneighb    - neighbours of all pencils in all directions
!    !     # glevel     - 1 == update only faces
!    !                  - 2 == update faces and edges
!    !                  - 3 == update faces,edges and corners
!    !
!    !     TODO FOR SPEEDUP: use mem-copy for x-direction
!    !
!    !========================================
      subroutine ghosts_update(XYZ, glevel, timer)
      	implicit none

!     !	function args
      	real(mytype) ,target :: XYZ(*)
      	integer,intent(in) :: glevel
      	double precision, intent(out) :: timer

!     !	usual stuff
     	integer :: i, p
     	integer :: sndid, recvid
     	integer :: request, status(mpi_status_size), ierr
     	integer :: gslabbuf_pos

      	if(.not. ghosts_set) then
      		print *,'P3DFFT error: call ghosts_init() before'
      		stop
      	endif

      	timer = -MPI_Wtime()

!     !	exchange ghost-cells faces
      	if(glevel .gt. 0) then
      	  do i=1,3
      		do p=1,2
      		  if(p .eq. 1) then; sndid=1; recvid=2;					! send from task with higher abs. coords to task with lower abs. coords.
      		  else;              sndid=2; recvid=1; endif		    ! send from task with lower abs. coords to task with higher abs. coords.

      		  call MPI_IRecv( XYZ(gmem_rstart(i,recvid)),																		 &
      		                 1, mpi_mytype_gslab_rrecv(i,recvid),																 &
      		                 gneighb_r(gslab_recvnghb(i,recvid)), i,															 &
      		                 mpi_comm, request, ierr)

      		  call MPI_Send( XYZ(gmem_rsndstart(i,sndid)),																		 &
      		                1, mpi_mytype_gslab_rsnd(i,sndid),																	 &
      		                gneighb_r(gslab_sndnghb(i,sndid)), i,																 &
      		                mpi_comm, ierr)

      		  call MPI_Wait( request, status, ierr )
      		enddo
      	  enddo
      	endif

!     !	exchange ghost-cell-edges swapped in y
      	if(glevel .gt. 1) then
      	  do p=1,2
      		if(p .eq. 1) then; sndid=1; recvid=2;
      		else;               sndid=2; recvid=1; endif

!     ! 	pack edges
      		gslabbuf_pos = 0
      		do i=4,5
      		  call MPI_Pack( XYZ(gmem_rsndstart(i,sndid)),																		 &
      		                 1, mpi_mytype_gslab_rsnd(i,sndid),																	 &
      		                 gslabbuf_snd, gslabbuf_bytesize, gslabbuf_pos,														 &
      		                 mpi_comm, ierr)
      		enddo

!     !		snd/recv data
      		call MPI_IRecv( gslabbuf_recv, gslabbuf_bytesize,																	 &
      		                MPI_PACKED,																							 &
      		                gneighb_r(gslab_recvnghb(4,recvid)), p,																 &
      		                mpi_comm, request, ierr)

      		call MPI_Send( gslabbuf_snd, gslabbuf_pos,																			 &
      		               MPI_PACKED,																							 &
      		               gneighb_r(gslab_sndnghb(4,sndid)), p,																 &
      		               mpi_comm, ierr)

      		call MPI_Wait( request, status, ierr )

!     !		unpack edges
      		gslabbuf_pos = 0
      		do i=4,5
      		  call MPI_Unpack( gslabbuf_recv, gslabbuf_bytesize, gslabbuf_pos,													 &
      		                   XYZ(gmem_rstart(i,recvid)),																		 &
      		                   1, mpi_mytype_gslab_rrecv(i,recvid),																 &
      		                   mpi_comm, ierr)
      		enddo
      	  enddo
      	endif

!     	exchange ghost-cell-edges and ghost-cell-corners swapped in x
      	if(glevel .gt. 1) then
      	  do p=1,2
      		if(p .eq. 1) then; sndid=1; recvid=2;
      		else;               sndid=2; recvid=1; endif

!     ! 	pack edges
      		gslabbuf_pos = 0
      		do i=6,9
      		  call MPI_Pack( XYZ(gmem_rsndstart(i,sndid)),																		 &
      		                 1, mpi_mytype_gslab_rsnd(i,sndid),																	 &
      		                 gslabbuf_snd, gslabbuf_bytesize, gslabbuf_pos,														 &
      		                 mpi_comm, ierr)
      		enddo

!     		pack corners
      		if(glevel .gt. 2) then
      		  do i=10,13
      			call MPI_Pack( XYZ(gmem_rsndstart(i,sndid)),																	 &
      		                   1, mpi_mytype_gslab_rsnd(i,sndid),																 &
      		                   gslabbuf_snd, gslabbuf_bytesize, gslabbuf_pos,													 &
      		                   mpi_comm, ierr)
      		  enddo
      		endif

!     !		snd/recv data
      		call MPI_IRecv( gslabbuf_recv, gslabbuf_bytesize,																	 &
      		                MPI_PACKED,																							 &
      		                gneighb_r(gslab_recvnghb(6,recvid)), p,																 &
      		                mpi_comm, request, ierr)

      		call MPI_Send( gslabbuf_snd, gslabbuf_pos,																			 &
      		               MPI_PACKED,																							 &
      		               gneighb_r(gslab_sndnghb(6,sndid)), p,																 &
      		               mpi_comm, ierr)

      		call MPI_Wait( request, status, ierr )

!     !		unpack data
      		gslabbuf_pos = 0
      		do i=6,9
      		  call MPI_Unpack( gslabbuf_recv, gslabbuf_bytesize, gslabbuf_pos,											         &
      		                   XYZ(gmem_rstart(i,recvid)),																		 &
      		                   1, mpi_mytype_gslab_rrecv(i,recvid),																 &
      		                   mpi_comm, ierr)
      		enddo

!     		unpack corners
      		if(glevel .gt. 2) then
      		  do i=10,13
      			call MPI_Unpack( gslabbuf_recv, gslabbuf_bytesize, gslabbuf_pos,											     &
      		                     XYZ(gmem_rstart(i,recvid)),																	 &
      		                     1, mpi_mytype_gslab_rrecv(i,recvid),															 &
      		                     mpi_comm, ierr)
      		  enddo
      		endif

      	  enddo
      	endif

      	timer = timer +MPI_Wtime()

      end subroutine ghosts_update

!    !========================================
!    !
!    !	ghosts_ijk2i() - calc 1d-array-position
!    !
!    !
!    !     TODO FOR SPEEDUP: precalc values in ghosts_init()
!    !
!    !========================================
      function ghosts_ijk2i(i,j,k)
      	implicit none

!     	function args
      	integer, intent(in) :: i,j,k
      	integer(kind=8) :: ghosts_ijk2i

!       -X plane (incl. X,Y edges)
      	if(i .lt. gproc_rstart(1) ) then								! |-x|
      	  if(j .lt. gproc_rstart(2) ) then								! |-x|-y|
      		if(k .lt. gproc_rstart(3) ) then							! |-x|-y|-z|
      		  ghosts_ijk2i = gmem_rstart(10,1)																					 &
      			+int( gslab_rsize(1,10)*gslab_rsize(2,10),8) *int( k-gproc_rstart(3) +gslab_rsize(3,10),8)						 &
      			+int( gslab_rsize(1,10)                  ,8) *int( j-gproc_rstart(2) +gslab_rsize(2,10),8)						 &
      			+                                             int( i-gproc_rstart(1) +gslab_rsize(1,10),8)
      		  return
      		else if(k .gt. gproc_rend(3) ) then						! |-x|-y|+z|
      		  ghosts_ijk2i = gmem_rstart(11,1)																					 &
      			+int( gslab_rsize(1,11)*gslab_rsize(2,11),8) *int( k  -gproc_rend(3)-1,8)										 &
      			+int( gslab_rsize(1,11)                  ,8) *int( j-gproc_rstart(2) +gslab_rsize(2,11),8)						 &
      			+                                             int( i-gproc_rstart(1) +gslab_rsize(1,11),8)
      		  return
      		else														! |-x|-y|zz|
      		  ghosts_ijk2i = gmem_rstart(8,1)																					 &
      			+int( gslab_rsize(1,8)*gslab_rsize(2,8),8) *int( k-gproc_rstart(3),8)											 &
      			+int( gslab_rsize(1,8)                 ,8) *int( j-gproc_rstart(2) +gslab_rsize(2,8),8)							 &
      			+                                           int( i-gproc_rstart(1) +gslab_rsize(1,8),8)
      		  return
      		endif
      	  else if(j .gt. gproc_rend(2) ) then							! |-x|+y|
      		if(k .lt. gproc_rstart(3) ) then							! |-x|+y|-z|
      		  ghosts_ijk2i = gmem_rstart(12,1)																					 &
      			+int( gslab_rsize(1,12)*gslab_rsize(2,12),8) *int( k-gproc_rstart(3) +gslab_rsize(3,12),8)						 &
      			+int( gslab_rsize(1,12)                  ,8) *int( j  -gproc_rend(2)-1,8)										 &
      			+                                             int( i-gproc_rstart(1) +gslab_rsize(1,12),8)
      		  return
      		else if(k .gt. gproc_rend(3) ) then						! |-x|+y|+z|
      		  ghosts_ijk2i = gmem_rstart(13,1)																					 &
      			+int( gslab_rsize(1,13)*gslab_rsize(2,13),8) *int( k  -gproc_rend(3)-1,8)										 &
      			+int( gslab_rsize(1,13)                  ,8) *int( j  -gproc_rend(2)-1,8)										 &
      			+                                             int( i-gproc_rstart(1) +gslab_rsize(1,13),8)
      		  return
      		else														! |-x|+y|zz|
      		  ghosts_ijk2i = gmem_rstart(9,1)																					 &
      			+int( gslab_rsize(1,9)*gslab_rsize(2,9),8) *int( k-gproc_rstart(3),8)											 &
      			+int( gslab_rsize(1,9)                 ,8) *int( j  -gproc_rend(2)-1,8)											 &
      			+                                           int( i-gproc_rstart(1) +gslab_rsize(1,9),8)
      		  return
      		endif
      	  else															! |-x|yy|
      		if(k .lt. gproc_rstart(3) ) then							! |-x|yy|-z|
      		  ghosts_ijk2i = gmem_rstart(6,1)																					 &
      		  +int( gslab_rsize(1,6)*gslab_rsize(2,6),8) *int( k-gproc_rstart(3)+gslab_rsize(3,6),8)							 &
      		  +int( gslab_rsize(1,6)                 ,8) *int( j-gproc_rstart(2),8)												 &
      		  +                                           int( i-gproc_rstart(1)+gslab_rsize(1,6),8)
      		  return
      		else if(k .gt. gproc_rend(3) ) then						! |-x|yy|+z|
      		  ghosts_ijk2i = gmem_rstart(7,1)																					 &
      		  +int( gslab_rsize(1,7)*gslab_rsize(2,7),8) *int( k  -gproc_rend(3)-1,8)											 &
      		  +int( gslab_rsize(1,7)                 ,8) *int( j-gproc_rstart(2),8)												 &
      		  +                                           int( i-gproc_rstart(1)+gslab_rsize(1,7),8)
      		  return
      		else														! |-x|yy|zz|
      		  ghosts_ijk2i = gmem_rstart(1,1)																					 &
      			+int( gslab_rsize(1,1)*gslab_rsize(2,1),8) *int( k-gproc_rstart(3),8)											 &
      			+int( gslab_rsize(1,1)                 ,8) *int( j-gproc_rstart(2),8)											 &
      			+                                           int( i-gproc_rstart(1) +gslab_rsize(1,1),8)
      		  return
      		endif
      	  endif
!       +X plane (incl. X,Y edges)
      	else if(i .gt. gproc_rend(1) ) then							! |+x|
      	  if(j .lt. gproc_rstart(2) ) then								! |+x|-y|
      		if(k .lt. gproc_rstart(3) ) then							! |+x|-y|-z|
      		  ghosts_ijk2i = gmem_rstart(10,2)																					 &
      			+int( gslab_rsize(1,10)*gslab_rsize(2,10),8) *int( k-gproc_rstart(3) +gslab_rsize(3,10),8)						 &
      			+int( gslab_rsize(1,10)                  ,8) *int( j-gproc_rstart(2) +gslab_rsize(2,10),8)						 &
      			+                                             int( i  -gproc_rend(1)-1,8)
      		  return
      		else if(k .gt. gproc_rend(3) ) then						! |+x|-y|+z|
      		  ghosts_ijk2i = gmem_rstart(11,2)																					 &
      			+int( gslab_rsize(1,11)*gslab_rsize(2,11),8) *int( k  -gproc_rend(3)-1,8)										 &
      			+int( gslab_rsize(1,11)                  ,8) *int( j-gproc_rstart(2) +gslab_rsize(2,11),8)						 &
      			+                                             int( i  -gproc_rend(1)-1,8)
      		  return
      		else														! |+x|-y|zz|
      		  ghosts_ijk2i = gmem_rstart(8,2)																					 &
      			+int( gslab_rsize(1,8)*gslab_rsize(2,8),8) *int( k-gproc_rstart(3),8)											 &
      			+int( gslab_rsize(1,8)                 ,8) *int( j-gproc_rstart(2) +gslab_rsize(2,8),8)							 &
      			+                                           int( i  -gproc_rend(1)-1,8)
      		  return
      		endif
      	  else if(j .gt. gproc_rend(2) ) then							! |+x|+y|
      		if(k .lt. gproc_rstart(3) ) then							! |+x|+y|-z|
      		  ghosts_ijk2i = gmem_rstart(12,2)																					 &
      			+int( gslab_rsize(1,12)*gslab_rsize(2,12),8) *int( k-gproc_rstart(3) +gslab_rsize(3,12),8)						 &
      			+int( gslab_rsize(1,12)                  ,8) *int( j  -gproc_rend(2)-1,8)										 &
      			+                                             int( i  -gproc_rend(1)-1,8)
      		  return
      		else if(k .gt. gproc_rend(3) ) then						! |+x|+y|+z|
      		  ghosts_ijk2i = gmem_rstart(13,2)																					 &
      			+int( gslab_rsize(1,13)*gslab_rsize(2,13),8) *int( k  -gproc_rend(3)-1,8)										 &
      			+int( gslab_rsize(1,13)                  ,8) *int( j  -gproc_rend(2)-1,8)										 &
      			+                                             int( i  -gproc_rend(1)-1,8)
      		  return
      		else														! |+x|+y|zz|
      		  ghosts_ijk2i = gmem_rstart(9,2)																					 &
      			+int( gslab_rsize(1,9)*gslab_rsize(2,9),8) *int( k-gproc_rstart(3),8)											 &
      			+int( gslab_rsize(1,9)                 ,8) *int( j  -gproc_rend(2)-1,8)											 &
      			+                                           int( i  -gproc_rend(1)-1,8)
      		  return
      		endif
      	  else															! |+x|yy|
      		if(k .lt. gproc_rstart(3) ) then							! |+x|yy|-z|
      		  ghosts_ijk2i = gmem_rstart(6,2)																					 &
      		  +int( gslab_rsize(1,6)*gslab_rsize(2,6),8) *int( k-gproc_rstart(3)+gslab_rsize(3,6),8)							 &
      		  +int( gslab_rsize(1,6)                 ,8) *int( j-gproc_rstart(2),8)												 &
      		  +                                           int( i  -gproc_rend(1)-1,8)
      		  return
      		else if(k .gt. gproc_rend(3) ) then						! |+x|yy|+z|
      		  ghosts_ijk2i = gmem_rstart(7,2)																					 &
      		  +int( gslab_rsize(1,7)*gslab_rsize(2,7),8) *int( k  -gproc_rend(3)-1,8)											 &
      		  +int( gslab_rsize(1,7)                 ,8) *int( j-gproc_rstart(2),8)												 &
      		  +                                           int( i  -gproc_rend(1)-1,8)
      		  return
      		else														! |+x|yy|zz|
      		  ghosts_ijk2i = gmem_rstart(1,2)																					 &
      			+int( gslab_rsize(1,1)*gslab_rsize(2,1),8) *int( k-gproc_rstart(3),8)											 &
      			+int( gslab_rsize(1,1)                 ,8) *int( j-gproc_rstart(2),8)											 &
      			+                                           int( i  -gproc_rend(1)-1,8)
      		  return
      		endif
      	  endif

!       -Y plane (incl. X edges)
!       (-k,+k was checked before: gproc_rstart(3) <= k <= gproc_rend(3))
      	else if(j .lt. gproc_rstart(2) ) then							! |xx|-y|
      	  if(k .lt. gproc_rstart(3) ) then								! |xx|-y|-z|
      	    ghosts_ijk2i = gmem_rstart(4,1)																						 &
      		  +int( gslab_rsize(1,4)*gslab_rsize(2,4),8) *int( k-gproc_rstart(3)+gslab_rsize(3,4),8)							 &
      		  +int( gslab_rsize(1,4)                 ,8) *int( j-gproc_rstart(2)+gslab_rsize(2,4),8)							 &
      		  +                                           int( i-gproc_rstart(1),8)
      	    return
      	  else if(k .gt. gproc_rend(3) ) then							! |xx|-y|+z|
      	    ghosts_ijk2i = gmem_rstart(5,1)																						 &
      		  +int( gslab_rsize(1,5)*gslab_rsize(2,5),8) *int( k  -gproc_rend(3)-1,8)											 &
      		  +int( gslab_rsize(1,5)                 ,8) *int( j-gproc_rstart(2)+gslab_rsize(2,5),8)							 &
      		  +                                           int( i-gproc_rstart(1),8)
      	    return
      	  else															! |xx|-y|zz|
      		ghosts_ijk2i = gmem_rstart(2,1)																						 &
      		  +int( gslab_rsize(1,2)*gslab_rsize(2,2),8) *int( k-gproc_rstart(3),8)												 &
      		  +int( gslab_rsize(1,2)                 ,8) *int( j-gproc_rstart(2) +gslab_rsize(2,2),8)							 &
      		  +                                           int( i-gproc_rstart(1),8)
      		  return
      	  endif

!       +Y plane (incl. X edges)
      	else if(j .gt. gproc_rend(2) ) then							! |xx|+y|
      	  if(k .lt. gproc_rstart(3) ) then								! |xx|+y|-z|
      	    ghosts_ijk2i = gmem_rstart(4,2)																						 &
      		  +int( gslab_rsize(1,4)*gslab_rsize(2,4),8) *int( k-gproc_rstart(3)+gslab_rsize(3,4),8)							 &
      		  +int( gslab_rsize(1,4)                 ,8) *int( j  -gproc_rend(2)-1,8)											 &
      		  +                                           int( i-gproc_rstart(1),8)
      	    return
      	  else if(k .gt. gproc_rend(3) ) then							! |xx|+y|+z|
      	    ghosts_ijk2i = gmem_rstart(5,2)							 &
      		  +int( gslab_rsize(1,5)*gslab_rsize(2,5),8) *int( k  -gproc_rend(3)-1,8)											 &
      		  +int( gslab_rsize(1,5)                 ,8) *int( j  -gproc_rend(2)-1,8)											 &
      		  +                                           int( i-gproc_rstart(1),8)
      	    return
      	  else															! |xx|+y|zz|
      		ghosts_ijk2i = gmem_rstart(2,2)																						 &
      		  +int( gslab_rsize(1,2)*gslab_rsize(2,2),8) *int( k-gproc_rstart(3),8)												 &
      		  +int( gslab_rsize(1,2)                 ,8) *int( j  -gproc_rend(2)-1,8)											 &
      		  +                                           int( i-gproc_rstart(1),8)
      		return
      	  endif

!       -Z plane
      	else if(k .lt. gproc_rstart(3) ) then							! |xx|yy|-z|
      	  ghosts_ijk2i = gmem_rstart(3,1)																						 &
      		+int( gslab_rsize(1,3)*gslab_rsize(2,3),8) *int( k-gproc_rstart(3) +gslab_rsize(3,3),8)								 &
      		+int( gslab_rsize(1,3)                 ,8) *int( j-gproc_rstart(2),8)												 &
      		+                                           int( i-gproc_rstart(1),8)
      	  return

!       +Z plane
      	else if(k .gt. gproc_rend(3) ) then							! |xx|yy|+z|
      	  ghosts_ijk2i = gmem_rstart(3,2)																						 &
      		+int( gslab_rsize(1,3)*gslab_rsize(2,3),8) *int( k  -gproc_rend(3) -1,8)											 &
      		+int( gslab_rsize(1,3)                 ,8) *int( j-gproc_rstart(2),8)												 &
      		+                                           int( i-gproc_rstart(1),8)
      	  return

!       no ghostcell
      	else															! |xx|yy|zz|
      	  ghosts_ijk2i = 1																										 &
      		+int(gproc_rsize(1)*gproc_rsize(2),8) *int(k-gproc_rstart(3),8)														 &
      		+int(gproc_rsize(1)               ,8) *int(j-gproc_rstart(2),8)														 &
      		+                                      int(i-gproc_rstart(1),8)
      	  return
      	endif

      end function ghosts_ijk2i

!    !========================================
!    !
!    !	ghosts_i002i() - calc 1d-array-position if only i might leave pencil
!    !
!    !     TODO FOR SPEEDUP: precalc values in ghosts_init()
!    !
!    !========================================
      function ghosts_i002i(i,j,k)
      	implicit none

!     	function args
      	integer, intent(in) :: i,j,k
      	integer(kind=8) :: ghosts_i002i

      	if(i .lt. gproc_rstart(1) ) then								! |-x|yy|zz|
      		  ghosts_i002i = gmem_rstart(1,1)																					 &
      		+int( gslab_rsize(1,1)*gslab_rsize(2,1),8) *int( k-gproc_rstart(3),8)												 &
      		+int( gslab_rsize(1,1)                 ,8) *int( j-gproc_rstart(2),8)												 &
      		+                                           int( i-gproc_rstart(1) +gslab_rsize(1,1),8)
      		  return

      	else if(i .gt. gproc_rend(1) ) then							! |+x|yy|zz|
      		  ghosts_i002i = gmem_rstart(1,2)																					 &
      		+int( gslab_rsize(1,1)*gslab_rsize(2,1),8) *int( k-gproc_rstart(3),8)												 &
      		+int( gslab_rsize(1,1)                 ,8) *int( j-gproc_rstart(2),8)												 &
      		+                                           int( i-gproc_rend(1) -1,8)
      		  return

!       no ghostcell
      	else															! |xx|yy|zz|
      	  ghosts_i002i = 1																										 &
      		+int(gproc_rsize(1)*gproc_rsize(2),8) *int(k-gproc_rstart(3),8)														 &
      		+int(gproc_rsize(1)               ,8) *int(j-gproc_rstart(2),8)														 &
      		+                                      int(i-gproc_rstart(1),8)
      	  return
      	endif

      end function ghosts_i002i

!    !========================================
!    !
!    !	ghosts_0j02i() - calc 1d-array-position if only j might leave pencil
!    !
!    !     TODO FOR SPEEDUP: precalc values in ghosts_init()
!    !
!    !========================================
      function ghosts_0j02i(i,j,k)
      	implicit none

!     	function args
      	integer, intent(in) :: i,j,k
      	integer(kind=8) :: ghosts_0j02i

      	if(j .lt. gproc_rstart(2) ) then								! |xx|-y|zz|
      		ghosts_0j02i = gmem_rstart(2,1)																						 &
      		  +int( gslab_rsize(1,2)*gslab_rsize(2,2),8) *int( k-gproc_rstart(3),8)												 &
      		  +int( gslab_rsize(1,2)                 ,8) *int( j-gproc_rstart(2) +gslab_rsize(2,2),8)							 &
      		  +                                           int( i-gproc_rstart(1),8)
      		  return

      	else if(j .gt. gproc_rend(2) ) then							! |xx|+y|zz|
      		ghosts_0j02i = gmem_rstart(2,2)																						 &
      		  +int( gslab_rsize(1,2)*gslab_rsize(2,2),8) *int( k-gproc_rstart(3),8)												 &
      		  +int( gslab_rsize(1,2)                 ,8) *int( j-gproc_rend(2)-1,8)												 &
      		  +                                           int( i-gproc_rstart(1),8)
      		return

!       no ghostcell
      	else															! |xx|yy|zz|
      	  ghosts_0j02i = 1																										 &
      		+int(gproc_rsize(1)*gproc_rsize(2),8) *int(k-gproc_rstart(3),8)														 &
      		+int(gproc_rsize(1)               ,8) *int(j-gproc_rstart(2),8)														 &
      		+                                      int(i-gproc_rstart(1),8)
      	  return
      	endif

      end function ghosts_0j02i

!    !========================================
!    !
!    !	ghosts_00k2i() - calc 1d-array-position if only k might leave pencil
!    !
!    !     TODO FOR SPEEDUP: precalc values in ghosts_init()
!    !
!    !========================================
      function ghosts_00k2i(i,j,k)
      	implicit none

!     	function args
      	integer, intent(in) :: i,j,k
      	integer(kind=8) :: ghosts_00k2i

      	if(k .lt. gproc_rstart(3) ) then								! |xx|yy|-z|
      	  ghosts_00k2i = gmem_rstart(3,1)																						 &
      			+int( gslab_rsize(1,3)*gslab_rsize(2,3),8) *int( k-gproc_rstart(3) +gslab_rsize(3,3),8)							 &
      			+int( gslab_rsize(1,3)                 ,8) *int( j-gproc_rstart(2),8)											 &
      			+                                           int( i-gproc_rstart(1),8)
      	  return

      	else if(k .gt. gproc_rend(3) ) then							! |xx|yy|+z|
      	  ghosts_00k2i = gmem_rstart(3,2)																						 &
      			+int( gslab_rsize(1,3)*gslab_rsize(2,3),8) *int( k-gproc_rend(3)-1,8)											 &
      			+int( gslab_rsize(1,3)                 ,8) *int( j-gproc_rstart(2),8)											 &
      			+                                           int( i-gproc_rstart(1),8)
      	  return

!       no ghostcell
      	else															! |xx|yy|zz|
      	  ghosts_00k2i = 1																										 &
      		+int(gproc_rsize(1)*gproc_rsize(2),8) *int(k-gproc_rstart(3),8)														 &
      		+int(gproc_rsize(1)               ,8) *int(j-gproc_rstart(2),8)														 &
      		+                                      int(i-gproc_rstart(1),8)
      	  return
      	endif

      end function ghosts_00k2i

!    !========================================
!    !
!    !	ghosts_ijk()
!    !
!    !========================================
      function ghosts_ijk(XgYZ,i,j,k)
      	implicit none

!    !	function vars
      	real(mytype) ,intent(in), target :: XgYZ(*)
      	integer, intent(in) :: i,j,k
      	real(mytype) :: ghosts_ijk

      	ghosts_ijk = XgYZ(ghosts_ijk2i(i,j,k))

      end function ghosts_ijk

!    !========================================
!    !
!    !  ghosts_ijk_debug()
!    !
!    !========================================
      function ghosts_ijk_debug(XgYZ,i,j,k,XgYZsize)
        implicit none

!    !  function vars
        integer(kind=8),intent(in) :: XgYZsize
        real(mytype) ,intent(in), target :: XgYZ(XgYZsize)
        integer, intent(in) :: i,j,k
        real(mytype) :: ghosts_ijk_debug

     !  other vars
        integer(kind=8) :: mempos

        mempos = ghosts_ijk2i(i,j,k)

        if( mempos > XgYZsize .or. mempos < 1) then
          write(*,*) 'ERROR: ghosts_ijk2i(i,j,k) returns mem-position larger than memory.'
        end if

        ghosts_ijk_debug = XgYZ(mempos)

      end function ghosts_ijk_debug

!    !========================================
!    !
!    !	ghosts_i00()
!    !
!    !========================================
      function ghosts_i00(XgYZ,i,j,k)
      	implicit none

!    !	function vars
      	real(mytype) ,intent(in), target :: XgYZ(*)
      	integer, intent(in) :: i,j,k
      	real(mytype) :: ghosts_i00

      	ghosts_i00 = XgYZ(ghosts_i002i(i,j,k))

      end function ghosts_i00

!    !========================================
!    !
!    !	ghosts_0j0()
!    !
!    !========================================
      function ghosts_0j0(XgYZ,i,j,k)
      	implicit none

!    !	function vars
      	real(mytype) ,intent(in), target :: XgYZ(*)
      	integer, intent(in) :: i,j,k
      	real(mytype) :: ghosts_0j0

      	ghosts_0j0 = XgYZ(ghosts_0j02i(i,j,k))

      end function ghosts_0j0

!    !========================================
!    !
!    !	ghosts_00k()
!    !
!    !========================================
      function ghosts_00k(XgYZ,i,j,k)
      	implicit none

!    !	function vars
      	real(mytype) ,intent(in), target :: XgYZ(*)
      	integer, intent(in) :: i,j,k
      	real(mytype) :: ghosts_00k

      	ghosts_00k = XgYZ(ghosts_00k2i(i,j,k))

      end function ghosts_00k

!    !========================================
!    !
!    !		check ghosts
!    !
!    !========================================
      subroutine ghosts_check(mem3d,istart,iend,n1,n2,n3,maxdelta,ierr)
      	implicit none

!!    !	function args
      	integer :: istart(3), iend(3)
      	real(mytype),	dimension(	istart(1):iend(1),					 &
      								istart(2):iend(2),					 &
      								istart(3):iend(3)),					 &
      						target, intent(inout) :: mem3d
      	integer :: n1,n2,n3
      	double precision, intent(out) :: maxdelta
      	integer, intent(out) :: ierr

!!    !	other vars
      	integer :: i,j,k
      	real(mytype) :: sin3_1, sin3_2, sin3_3, sin3_4
      	real(mytype) :: sinyz
      	real(mytype) :: delta1, delta2, delta3, delta4
      	double precision :: maxdelta1,maxdelta2,maxdelta3,maxdelta4
      	double precision :: pi2
      	double precision :: t

      	double precision, dimension(:), allocatable :: coordX
      	double precision, dimension(:), allocatable :: coordY
      	double precision, dimension(:), allocatable :: coordZ
      	double precision :: Lx,Ly,Lz

      	ierr = 0
      	maxdelta  = 0.d0
      	maxdelta1 = 0.d0
      	maxdelta2 = 0.d0
      	maxdelta3 = 0.d0
      	maxdelta4 = 0.d0

      	Lx = 2.d0
      	Ly = 2.d0
      	Lz = 2.d0
      	pi2 = 8.d0 *datan(1.d0)

      	if(.not. ghosts_set) then
      		print *,'P3DFFT error: call p3dfft_init_ghosts() before'
      		return
      	endif

      	allocate(coordX(-goverlap:n1+goverlap), stat=ierr)
      	allocate(coordY(-goverlap:n2+goverlap), stat=ierr)
      	allocate(coordZ(-goverlap:n3+goverlap), stat=ierr)
      	do i=-goverlap,n1+goverlap
      		coordX(i) = Lx *((2*i-1)/dble(2*n1))
      	enddo
      	do j=-goverlap,n2+goverlap
      		coordY(j) = Ly *((2*j-1)/dble(2*n2))
      	enddo
      	do k=-goverlap,n3+goverlap
      		coordZ(k) = Lz *((2*k-1)/dble(2*n3))
      	enddo

!    !	fill array (pencil) with sinus-wave
      	do k=istart(3),iend(3)
      	   do j=istart(2),iend(2)
      	      sinyz = sin(pi2/Ly *coordY(j)) *sin(pi2/Lz *coordZ(k))
      	      do i=istart(1),iend(1)
      	         mem3d(i,j,k) = sin(pi2/Lx *coordX(i)) *sinyz +5.d0
      	      enddo
      	   enddo
      	enddo

!     !	update ghostcells
     	call ghosts_update(mem3d,3,t)

!    !	check array (pencil) faces
      	do k=istart(3),iend(3)
      	   do j=istart(2),iend(2)
      	      do i=istart(1),iend(1)
      	        sin3_1 = sin(pi2/Lx *coordX(i-goverlap)) *sin(pi2/Ly *coordY(j)      ) *sin(pi2/Lz *coordZ(k)      ) +5.d0
      	        sin3_2 = sin(pi2/Lx *coordX(i)      ) *sin(pi2/Ly *coordY(j-goverlap)) *sin(pi2/Lz *coordZ(k)      ) +5.d0
      	        sin3_3 = sin(pi2/Lx *coordX(i)      ) *sin(pi2/Ly *coordY(j)      ) *sin(pi2/Lz *coordZ(k-goverlap)) +5.d0
      	        delta1 = abs( ghosts_i00(mem3d,i-goverlap,j,k) -sin3_1 )
      	        delta2 = abs( ghosts_0j0(mem3d,i,j-goverlap,k) -sin3_2 )
      	        delta3 = abs( ghosts_00k(mem3d,i,j,k-goverlap) -sin3_3 )
      	        if( delta1 .gt. maxdelta1) then
      	            maxdelta1 = delta1
      	            if(delta1 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |-x|yy|zz| error ', delta1, ' at ',i,j,k
      	        endif
      	        if( delta2 .gt. maxdelta2) then
      	            maxdelta2 = delta2
                   if(delta2 .gt. 1e-4)	 		 &
                     write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |xx|-y|zz| error ', delta2, ' at ',i,j,k
      	        endif
      	        if( delta3 .gt. maxdelta3) then
      	            maxdelta3 = delta3
                   if(delta3 .gt. 1e-4)			 &
                     write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |xx|yy|-z| error ', delta3, ' at ',i,j,k
      	       endif
      	       sin3_1 = sin(pi2/Lx *coordX(i+goverlap)) *sin(pi2/Ly *coordY(j)      ) *sin(pi2/Lz *coordZ(k)      ) +5.d0
      	       sin3_2 = sin(pi2/Lx *coordX(i)      ) *sin(pi2/Ly *coordY(j+goverlap)) *sin(pi2/Lz *coordZ(k)      ) +5.d0
      	       sin3_3 = sin(pi2/Lx *coordX(i)      ) *sin(pi2/Ly *coordY(j)      ) *sin(pi2/Lz *coordZ(k+goverlap)) +5.d0
      	       delta1 = abs( ghosts_i00(mem3d,i+goverlap,j,k) -sin3_1 )
      	       delta2 = abs( ghosts_0j0(mem3d,i,j+goverlap,k) -sin3_2 )
      	       delta3 = abs( ghosts_00k(mem3d,i,j,k+goverlap) -sin3_3 )
      	       if( delta1 .gt. maxdelta1) then
      	            maxdelta1 = delta1
      	            if(delta1 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |+x|yy|zz| error ', delta1, ' at ',i,j,k
      	       endif
      	       if( delta2 .gt. maxdelta2) then
      	            maxdelta2 = delta2
      	            if(delta2 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |xx|+y|zz| error ', delta2, ' at ',i,j,k
      	       endif
      	       if( delta3 .gt. maxdelta3) then
      	            maxdelta3 = delta3
      	            if(delta3 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |xx|yy|+z| error ', delta3, ' at ',i,j,k
      	       endif
      	     enddo
      	   enddo
      	enddo

!    !	check array (pencil) edges in x-direction
      	do k=istart(3),iend(3)
      	   do j=istart(2),iend(2)
      	      do i=istart(1),iend(1)
      	        sin3_1 = sin(pi2/Lx *coordX(i)) *sin(pi2/Ly *coordY(j-goverlap)) *sin(pi2/Lz *coordZ(k-goverlap)) +5.d0
      	        sin3_2 = sin(pi2/Lx *coordX(i)) *sin(pi2/Ly *coordY(j-goverlap)) *sin(pi2/Lz *coordZ(k+goverlap)) +5.d0
      	        sin3_3 = sin(pi2/Lx *coordX(i)) *sin(pi2/Ly *coordY(j+goverlap)) *sin(pi2/Lz *coordZ(k-goverlap)) +5.d0
      	        sin3_4 = sin(pi2/Lx *coordX(i)) *sin(pi2/Ly *coordY(j+goverlap)) *sin(pi2/Lz *coordZ(k+goverlap)) +5.d0
      	        delta1 = abs( ghosts_ijk(mem3d,i,j-goverlap,k-goverlap) -sin3_1 )
      	        delta2 = abs( ghosts_ijk(mem3d,i,j-goverlap,k+goverlap) -sin3_2 )
      	        delta3 = abs( ghosts_ijk(mem3d,i,j+goverlap,k-goverlap) -sin3_3 )
      	        delta4 = abs( ghosts_ijk(mem3d,i,j+goverlap,k+goverlap) -sin3_4 )
      	        if( delta1 .gt. maxdelta1) then
      	            maxdelta1 = delta1
      	            if(delta1 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |xx|-y|-z| error ', delta1, ' at ',i,j-goverlap,k-goverlap
      	        endif
      	        if( delta2 .gt. maxdelta2) then
      	            maxdelta2 = delta2
      	            if(delta2 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |xx|-y|+z| error ', delta2, ' at ',i,j-goverlap,k+goverlap
      	        endif
      	        if( delta3 .gt. maxdelta3) then
      	            maxdelta3 = delta3
      	            if(delta3 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |xx|+y|-z| error ', delta3, ' at ',i,j+goverlap,k-goverlap
      	        endif
      	        if( delta4 .gt. maxdelta4) then
      	            maxdelta4 = delta4
      	            if(delta4 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |xx|+y|+z| error ', delta4, ' at ',i,j+goverlap,k+goverlap
      	        endif
      	      enddo
      	   enddo
      	enddo

!    !	check array (pencil) edges in y-direction
      	do k=istart(3),iend(3)
      	   do j=istart(2),iend(2)
      	      do i=istart(1),iend(1)
      	        sin3_1 = sin(pi2/Lx *coordX(i-goverlap)) *sin(pi2/Ly *coordY(j)) *sin(pi2/Lz *coordZ(k-goverlap)) +5.d0
      	        sin3_2 = sin(pi2/Lx *coordX(i-goverlap)) *sin(pi2/Ly *coordY(j)) *sin(pi2/Lz *coordZ(k+goverlap)) +5.d0
      	        sin3_3 = sin(pi2/Lx *coordX(i+goverlap)) *sin(pi2/Ly *coordY(j)) *sin(pi2/Lz *coordZ(k-goverlap)) +5.d0
      	        sin3_4 = sin(pi2/Lx *coordX(i+goverlap)) *sin(pi2/Ly *coordY(j)) *sin(pi2/Lz *coordZ(k+goverlap)) +5.d0
      	        delta1 = abs( ghosts_ijk(mem3d,i-goverlap,j,k-goverlap) -sin3_1 )
      	        delta2 = abs( ghosts_ijk(mem3d,i-goverlap,j,k+goverlap) -sin3_2 )
      	        delta3 = abs( ghosts_ijk(mem3d,i+goverlap,j,k-goverlap) -sin3_3 )
      	        delta4 = abs( ghosts_ijk(mem3d,i+goverlap,j,k+goverlap) -sin3_4 )
      	        if( delta1 .gt. maxdelta1) then
      	            maxdelta1 = delta1
      	            if(delta1 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |-x|yy|-z| error ', delta1, ' at ',i-goverlap,j,k-goverlap
      	        endif
      	        if( delta2 .gt. maxdelta2) then
      	            maxdelta2 = delta2
      	            if(delta2 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |-x|yy|+z| error ', delta2, ' at ',i-goverlap,j,k+goverlap
      	        endif
      	        if( delta3 .gt. maxdelta3) then
      	            maxdelta3 = delta3
      	            if(delta3 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |+x|yy|-z| error ', delta3, ' at ',i+goverlap,j,k-goverlap
      	        endif
      	        if( delta4 .gt. maxdelta4) then
      	            maxdelta4 = delta4
      	            if(delta4 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |+x|yy|+z| error ', delta4, ' at ',i+goverlap,j,k+goverlap
      	        endif
      	      enddo
      	   enddo
      	enddo

!    !	check array (pencil) edges in z-direction
      	do k=istart(3),iend(3)
      	   do j=istart(2),iend(2)
      	      do i=istart(1),iend(1)
      	        sin3_1 = sin(pi2/Lx *coordX(i-goverlap)) *sin(pi2/Ly *coordY(j-goverlap)) *sin(pi2/Lz *coordZ(k)) +5.d0
      	        sin3_2 = sin(pi2/Lx *coordX(i-goverlap)) *sin(pi2/Ly *coordY(j+goverlap)) *sin(pi2/Lz *coordZ(k)) +5.d0
      	        sin3_3 = sin(pi2/Lx *coordX(i+goverlap)) *sin(pi2/Ly *coordY(j-goverlap)) *sin(pi2/Lz *coordZ(k)) +5.d0
      	        sin3_4 = sin(pi2/Lx *coordX(i+goverlap)) *sin(pi2/Ly *coordY(j+goverlap)) *sin(pi2/Lz *coordZ(k)) +5.d0
      	        delta1 = abs( ghosts_ijk(mem3d,i-goverlap,j-goverlap,k) -sin3_1 )
      	        delta2 = abs( ghosts_ijk(mem3d,i-goverlap,j+goverlap,k) -sin3_2 )
      	        delta3 = abs( ghosts_ijk(mem3d,i+goverlap,j-goverlap,k) -sin3_3 )
      	        delta4 = abs( ghosts_ijk(mem3d,i+goverlap,j+goverlap,k) -sin3_4 )
      	        if( delta1 .gt. maxdelta1) then
      	            maxdelta1 = delta1
      	            if(delta1 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |-x|-y|zz| error ', delta1, ' at ',i-goverlap,j-goverlap,k
      	        endif
      	        if( delta2 .gt. maxdelta2) then
      	            maxdelta2 = delta2
      	            if(delta2 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |-x|+y|zz| error ', delta2, ' at ',i-goverlap,j+goverlap,k
      	        endif
      	        if( delta3 .gt. maxdelta3) then
      	            maxdelta3 = delta3
      	            if(delta3 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |+x|-y|zz| error ', delta3, ' at ',i+goverlap,j-goverlap,k
      	        endif
      	        if( delta4 .gt. maxdelta4) then
      	            maxdelta4 = delta4
      	            if(delta4 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |+x|+y|zz| error ', delta4, ' at ',i+goverlap,j+goverlap,k
      	        endif
      	      enddo
      	   enddo
      	enddo

!    !	check array (pencil) corner in minus-z-direction
      	do k=istart(3),iend(3)
      	   do j=istart(2),iend(2)
      	      do i=istart(1),iend(1)
      	        sin3_1 = sin(pi2/Lx *coordX(i-goverlap)) *sin(pi2/Ly *coordY(j-goverlap)) *sin(pi2/Lz *coordZ(k-goverlap)) +5.d0
      	        sin3_2 = sin(pi2/Lx *coordX(i-goverlap)) *sin(pi2/Ly *coordY(j+goverlap)) *sin(pi2/Lz *coordZ(k-goverlap)) +5.d0
      	        sin3_3 = sin(pi2/Lx *coordX(i+goverlap)) *sin(pi2/Ly *coordY(j-goverlap)) *sin(pi2/Lz *coordZ(k-goverlap)) +5.d0
      	        sin3_4 = sin(pi2/Lx *coordX(i+goverlap)) *sin(pi2/Ly *coordY(j+goverlap)) *sin(pi2/Lz *coordZ(k-goverlap)) +5.d0
      	        delta1 = abs( ghosts_ijk(mem3d,i-goverlap,j-goverlap,k-goverlap) -sin3_1 )
      	        delta2 = abs( ghosts_ijk(mem3d,i-goverlap,j+goverlap,k-goverlap) -sin3_2 )
      	        delta3 = abs( ghosts_ijk(mem3d,i+goverlap,j-goverlap,k-goverlap) -sin3_3 )
      	        delta4 = abs( ghosts_ijk(mem3d,i+goverlap,j+goverlap,k-goverlap) -sin3_4 )
      	        if( delta1 .gt. maxdelta1) then
      	            maxdelta1 = delta1
      	            if(delta1 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |-x|-y|-z| error ', delta1, ' at ',i-goverlap,j-goverlap,k-goverlap
      	        endif
      	        if( delta2 .gt. maxdelta2) then
      	            maxdelta2 = delta2
      	            if(delta2 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |-x|+y|-z| error ', delta2, ' at ',i-goverlap,j+goverlap,k-goverlap
      	        endif
      	        if( delta3 .gt. maxdelta3) then
      	            maxdelta3 = delta3
      	            if(delta3 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |+x|-y|-z| error ', delta3, ' at ',i+goverlap,j-goverlap,k-goverlap
      	        endif
      	        if( delta4 .gt. maxdelta4) then
      	            maxdelta4 = delta4
      	            if(delta4 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |+x|+y|-z| error ', delta4, ' at ',i+goverlap,j+goverlap,k-goverlap
      	        endif
      	      enddo
      	   enddo
      	enddo

!    !	check array (pencil) corner in plus-z-direction
      	do k=istart(3),iend(3)
      	   do j=istart(2),iend(2)
      	      do i=istart(1),iend(1)
      	        sin3_1 = sin(pi2/Lx *coordX(i-goverlap)) *sin(pi2/Ly *coordY(j-goverlap)) *sin(pi2/Lz *coordZ(k+goverlap)) +5.d0
      	        sin3_2 = sin(pi2/Lx *coordX(i-goverlap)) *sin(pi2/Ly *coordY(j+goverlap)) *sin(pi2/Lz *coordZ(k+goverlap)) +5.d0
      	        sin3_3 = sin(pi2/Lx *coordX(i+goverlap)) *sin(pi2/Ly *coordY(j-goverlap)) *sin(pi2/Lz *coordZ(k+goverlap)) +5.d0
      	        sin3_4 = sin(pi2/Lx *coordX(i+goverlap)) *sin(pi2/Ly *coordY(j+goverlap)) *sin(pi2/Lz *coordZ(k+goverlap)) +5.d0
      	        delta1 = abs( ghosts_ijk(mem3d,i-goverlap,j-goverlap,k+goverlap) -sin3_1 )
      	        delta2 = abs( ghosts_ijk(mem3d,i-goverlap,j+goverlap,k+goverlap) -sin3_2 )
      	        delta3 = abs( ghosts_ijk(mem3d,i+goverlap,j-goverlap,k+goverlap) -sin3_3 )
      	        delta4 = abs( ghosts_ijk(mem3d,i+goverlap,j+goverlap,k+goverlap) -sin3_4 )
      	        if( delta1 .gt. maxdelta1) then
      	            maxdelta1 = delta1
      	            if(delta1 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |-x|-y|+z| error ', delta1, ' at ',i-goverlap,j-goverlap,k+goverlap
      	        endif
      	        if( delta2 .gt. maxdelta2) then
      	            maxdelta2 = delta2
      	            if(delta2 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |-x|+y|+z| error ', delta2, ' at ',i-goverlap,j+goverlap,k+goverlap
      	        endif
      	        if( delta3 .gt. maxdelta3) then
      	            maxdelta3 = delta3
      	            if(delta3 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |+x|-y|+z| error ', delta3, ' at ',i+goverlap,j-goverlap,k+goverlap
      	        endif
      	        if( delta4 .gt. maxdelta4) then
      	            maxdelta4 = delta4
      	            if(delta4 .gt. 1e-4)		 &
      	              write(*,fmt='(a,f24.16,a,i4,i4,i4)') '  |+x|+y|+z| error ', delta4, ' at ',i+goverlap,j+goverlap,k+goverlap
      	        endif
      	      enddo
      	   enddo
      	enddo

      	deallocate(coordX)
      	deallocate(coordY)
      	deallocate(coordZ)

      	maxdelta = max(maxdelta1, maxdelta2, maxdelta3, maxdelta4)

      end subroutine ghosts_check

      end module
