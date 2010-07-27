!=======================================================================
!                 ghost cell support
!               ----------------------
!
!  how to use:
! ==============
!     1) init ghost-support:
!       ==>> call p3dfft_init_ghosts(goverlap)
!       goverlap = overlap of ghost-cell-slab
!                  (limitations: only one goverlap supported)
!       (of cause we are still assuming periodic boundary conditions)
!
!     2) allocation of memory
!       each array which needs ghost-cell-support has to have
!       the min. size:	fsize(1)+2*goverlap *
!                       fsize(2)+2*goverlap *
!                       fsize(3)+2*goverlap
!
!     3) exchange ghost-cells among processors
!       to update the ghost-cells among the proceesors
!       ==>> call update_rghosts(array).
!       This will assume data represents physical space!
!
!     4) use/read ghost-cells
!       if you want to access/read ghost-cells
!       ==>> call ijk2i(i,j,k)
!       which returns the position in a 1D-array as an integer*8
!
!  array-memory:
! ==============
!     1) size if complex-type: fsize(1)+2*goverlap *
!                              fsize(2)+2*goverlap *
!                              fsize(3)+2*goverlap
!     2a) structure in memory:
!      |000000000000000000000000000|11111|22222|33333|44444|55555|66666|
!                 ^^^^              ^     ^     ^     ^     ^     ^
!          input/output of FFT      A     B     C     D     E     F
!
!      ghost-cell-slab at each pencil-side:
!      A:  -X,Y*Z	A(istart(1)-i, j,k)		=> gmem_start(1)
!      B:  +X,Y*Z	B(i,j,k)				=> gmem_start(2)
!      C:  -Y,X*Z	C(i,j,k)				=> gmem_start(3)
!      D:  +Y,X*Z	D(i,j,k)				=> gmem_start(4)
!      E:  -Z,X*Y	E(i,j,k)				=> gmem_start(5)
!      F:  +Z,X*Y	F(i,j,k)				=> gmem_start(6)
!
!     2b) a more closer look at "A":
!      ....|aaaaaabbbbbbccccccdddddd|....
!           ^     ^     ^     ^
!           a1    b1    c1    d1
!
!      each ghost-cell-slab is build up by no. goverlap planes:
!      a1:  Y*Z-plane with lowest X-coord
!       |
!      d1:  Y*Z-plane with highest X-coord
!
!  important variables:
! =====================
!     1) goverlap
!       the data overlap between two processors - equal to the thickness
!       of a ghost-slab
!       This version supports only one fixed overlap value, because it simplifies
!       the coding. If different overlaps for different arrays are needed the
!       following varables must be dependent on the overlap-value and therefor
!       should have one additional dimension.
!       A good idea would be to make p3dfft_init_ghosts(..) return an id,which
!       has to get set when calling update_ghosts(data,id).
!
!     2) gmem_start(1..6) (-X,-Y,-Z,+X,+Y,+Z)
!       start of data in memory for ghost-cells of certain side
!
!     3) gmess_size(1..3) (X,Y,Z)
!       size of mpi message needed to send/recieve ghost-slab each dimension(1..3)
!
!     4a) gslab_start(1..3, 1..6)
!       relative start coord. of ghost-slab in pencil ijk(1..3), [neg./pos side(1..2) for each dimension(1..3)](1..6)
!           -X   -Y   -Z   +X   +Y   +Z
!         i
!         j
!         k
!
!     4b) gslab_end(1..3, 1..6)
!       relative end coord. of ghost-slab in pencil ijk(1..3), [neg./pos side(1..2) for each dimension(1..3)](1..6)
!           -X   -Y   -Z   +X   +Y   +Z
!         i
!         j
!         k
!
!     4c) gslab_size(1..3, 1..3)
!       size of ghost-slab for each dimension(1..3), ijk(1..3)
!
!     6) gneighb_(numtasks*(3*2))
!       extended processor grid with all processor neighbours in -/+XYZ(1..6)
!       including periodic boundary neighbours
!        example for 2 processors with gneighb_r(:) = 0,0,0,0,1,1,1,1,1,1,0,0
!         pencil-side: -X  +X  -Y  +Y  -Z  +Z
!         taskid 0:     0   0   0   0   1   1
!         taskid 1:     1   1   1   1   0   0
!=======================================================================

      subroutine p3dfft_init_ghosts(goverlap_in)
      	implicit none

!     	function args
      	integer, intent(in) :: goverlap_in

!     	other vars
      	integer, dimension(:,:), allocatable :: coords2id
      	integer :: pneighb_r(6), pneighb_c(6)
      	integer :: pCoo(2)
      	integer :: rstart(3), rend(3), cstart(3), cend(3)
      	integer :: i,d,mp,mp1,msize,ierr,dmp,m,p

      	if(.not. mpi_set) then
      		print *,'P3DFFT error: call setup before other routines'
      		ierr = 1
      		call abort
      	endif
      	
      	goverlap = goverlap_in

      	call get_dims(gproc_rstart, gproc_rend, gproc_rsize, 1)
      	call get_dims(gproc_cstart, gproc_cend, gproc_csize, 2)

      	if(     goverlap .gt. gproc_csize(1) &
      	   .or. goverlap .gt. gproc_csize(2) &
      	   .or. goverlap .gt. gproc_csize(3) ) then
      			print *,'P3DFFT error: in one direction goverlap is greater then complex pencil size'
      			ierr = 1
      			call abort
      	endif

!     	set mpi message sizes for exchange of ghost-cells between processors
      	gmess_rsize(1) = goverlap       *gproc_rsize(2) *gproc_rsize(3)
      	gmess_rsize(2) = gproc_rsize(1) *goverlap       *gproc_rsize(3)
      	gmess_rsize(3) = gproc_rsize(1) *gproc_rsize(2) *goverlap

      	gmess_csize(1) = goverlap       *gproc_csize(2) *gproc_csize(3) *2
      	gmess_csize(2) = gproc_csize(1) *goverlap       *gproc_csize(3) *2
      	gmess_csize(3) = gproc_csize(1) *gproc_csize(2) *goverlap       *2

!     	precalc ghost-slab start-position in memory of in/out array (-X,-Y
      	gmem_rstart(1) = int(gproc_csize(1),8)  *int(gproc_csize(2),8) *int(gproc_csize(3),8) *2 +1
     	gmem_rstart(2) = gmem_rstart(1) +int(goverlap      *gproc_rsize(2) *gproc_rsize(3) ,8)
      	gmem_rstart(3) = gmem_rstart(2) +int(gproc_rsize(1) *goverlap      *gproc_rsize(3) ,8)
      	gmem_rstart(4) = gmem_rstart(3) +int(gproc_rsize(1) *gproc_rsize(2) *goverlap      ,8)
      	gmem_rstart(5) = gmem_rstart(4) +int(goverlap       *gproc_rsize(2) *gproc_rsize(3),8)
      	gmem_rstart(6) = gmem_rstart(5) +int(gproc_rsize(1) *goverlap       *gproc_rsize(3),8)

      	gmem_cstart(1) = int(gproc_csize(1),8)  *int(gproc_csize(2),8) *int(gproc_csize(3),8) *2 +1
     	gmem_cstart(2) = gmem_cstart(1) +int(goverlap      *gproc_csize(2) *gproc_csize(3) ,8)
      	gmem_cstart(3) = gmem_cstart(2) +int(gproc_csize(1) *goverlap      *gproc_csize(3) ,8)
      	gmem_cstart(4) = gmem_cstart(3) +int(gproc_csize(1) *gproc_csize(2) *goverlap      ,8)
      	gmem_cstart(5) = gmem_cstart(4) +int(goverlap       *gproc_csize(2) *gproc_csize(3),8)
      	gmem_cstart(6) = gmem_cstart(5) +int(gproc_csize(1) *goverlap       *gproc_csize(3),8)

!     	write(*,*) 'gmem_start: ', gmem_rstart(:)

!     	allocate send/recieve buffer
      	msize = gmess_rsize(1)
      	if(gmess_rsize(2) .gt. msize) msize = gmess_rsize(2)
      	if(gmess_rsize(3) .gt. msize) msize = gmess_rsize(3)
      	if(gmess_csize(1) .gt. msize) msize = gmess_csize(1)
      	if(gmess_csize(2) .gt. msize) msize = gmess_csize(2)
      	if(gmess_csize(3) .gt. msize) msize = gmess_csize(3)
      	allocate(gbuf_snd(msize),  stat=ierr)
      	allocate(gbuf_recv(msize), stat=ierr)

!     	precalc ghost-slab size,start,end (checked)
      	gslab_rsize(:,1) = gproc_rsize(:)
      	gslab_rsize(:,2) = gproc_rsize(:)
      	gslab_rsize(:,3) = gproc_rsize(:)
      	gslab_csize(:,1) = gproc_csize(:)
      	gslab_csize(:,2) = gproc_csize(:)
      	gslab_csize(:,3) = gproc_csize(:)

      	gslab_rsize(1,1) = goverlap
      	gslab_rsize(2,2) = goverlap
      	gslab_rsize(3,3) = goverlap
      	gslab_csize(1,1) = goverlap
      	gslab_csize(2,2) = goverlap
      	gslab_csize(3,3) = goverlap
      	do d=1, 3	! loop over X,Y,Z sides
      	  do mp=1,2		! loop over direction (-/+)
      		dmp = (mp-1)*3 +d
      		if(mp .eq. 1) then; m=1; p=0; else; m=0; p=1; endif
      		
      		gslab_rstart(1,dmp) = m*1 +p*(gproc_rsize(1)-gslab_rsize(1,d)+1)
      		gslab_rstart(2,dmp) = m*1 +p*(gproc_rsize(2)-gslab_rsize(2,d)+1)
      		gslab_rstart(3,dmp) = m*1 +p*(gproc_rsize(3)-gslab_rsize(3,d)+1)

      		gslab_rend(1,dmp) = gslab_rstart(1,dmp) +gslab_rsize(1,d) -1
      		gslab_rend(2,dmp) = gslab_rstart(2,dmp) +gslab_rsize(2,d) -1
      		gslab_rend(3,dmp) = gslab_rstart(3,dmp) +gslab_rsize(3,d) -1

      		gslab_cstart(1,dmp) = m*1 +p*(gproc_csize(1)-gslab_csize(1,d)+1)
      		gslab_cstart(2,dmp) = m*1 +p*(gproc_csize(2)-gslab_csize(2,d)+1)
      		gslab_cstart(3,dmp) = m*1 +p*(gproc_csize(3)-gslab_csize(3,d)+1)

      		gslab_cend(1,dmp) = gslab_rstart(1,dmp) +gslab_rsize(1,d) -1
      		gslab_cend(2,dmp) = gslab_rstart(2,dmp) +gslab_rsize(2,d) -1
      		gslab_cend(3,dmp) = gslab_rstart(3,dmp) +gslab_rsize(3,d) -1
      	  enddo
      	enddo

!      	print *, 'gproc_rstart ',taskid, gproc_rstart(:)
!      	call MPI_Barrier(mpi_comm_cart,ierr)
!      	print *, 'gproc_rend ',taskid, gproc_rend(:)
!      	call MPI_Barrier(mpi_comm_world,ierr)
!      	do d=1, 3
!      		write(*,*) 'gslab_rsize ',taskid,d, gslab_rsize(d,:)
!      	enddo
!      	call MPI_Barrier(mpi_comm_cart,ierr)
!      	do d=1, 3
!      		write(*,*) 'gslab_rstart',taskid,d, gslab_rstart(d,:)
!      	enddo
!      	call MPI_Barrier(mpi_comm_cart,ierr)
!      	do d=1, 3
!      		write(*,*) 'gslab_rend  ',taskid,d, gslab_rend(d,:)
!      	enddo
!      	call MPI_Barrier(mpi_comm_cart,ierr)
      

!     !	build help-array with coordinates to proc-ids
      	allocate(coords2id(-1:iproc,-1:jproc))
      	coords2id(:,:) = -1
      	do i=0,numtasks-1
      		coords2id( proc_id2coords(2*i), proc_id2coords(2*i+1)) = i
      	enddo
!     !	extend array with one line more on each side (periodicy)
      	coords2id(-1,   0:jproc-1) = coords2id(iproc-1,0:jproc-1)
      	coords2id(iproc,0:jproc-1) = coords2id(0,      0:jproc-1)
      	coords2id(0:iproc-1,   -1) = coords2id(0:iproc-1,jproc-1)
      	coords2id(0:iproc-1,jproc) = coords2id(0:iproc-1,      0)

      	pCoo(1) = proc_id2coords(taskid*2)
      	pCoo(2) = proc_id2coords(taskid*2+1)
      	do d=0, 2	! loop over X,Y,Z sides
      	  do mp=1,2		! loop over direction (-/+)
      		mp1 = mp*2-3	! 1=>-1 , 2=>+1

!      		! build pneighb_r
      		if(	d .eq. 0 .or.			 &	! no decomp. in X  
      			d .eq. 1 .and. iproc .eq. 1 .or. &	! no decomp. in Y  
      			d .eq. 2 .and. jproc .eq. 1) then	! no decomp. in Z
      				pneighb_r(d*2+mp) = taskid
      		else if(d .eq. 1) then					! decomp. in Y
      			pneighb_r(d*2+mp) = coords2id(pCoo(1)+mp1, pCoo(2))
      		else if(d .eq. 2) then					! decomp. in Z
      			pneighb_r(d*2+mp) = coords2id(pCoo(1), pCoo(2)+mp1)
      		endif

!      		! build pneighb_c
      		if(	d .eq. 0 .and. iproc .eq. 1 .or. &	! no decomp. in X  
      			d .eq. 1 .and. jproc .eq. 1 .or. &	! no decomp. in Y  
      			d .eq. 2) then						! no decomp. in Z
      				pneighb_c(d*2+mp) = taskid
      		else if(d .eq. 0) then					! decomp. in X
      			pneighb_c(d*2+mp) = coords2id(pCoo(1)+mp1, pCoo(2))
      		else if(d .eq. 1) then					! decomp. in Y
      			pneighb_c(d*2+mp) = coords2id(pCoo(1), pCoo(2)+mp1)
      		endif

      	  enddo
      	enddo

!     ! distribute pneigh_r
      	allocate(gneighb_r(numtasks*6))
      	call MPI_Allgather(	pneighb_r,	6, MPI_INTEGER, &
      						gneighb_r,	6, MPI_INTEGER, &
      						mpi_comm_cart,ierr)

!     ! distribute pneigh_c
      	allocate(gneighb_c(numtasks*6))
      	call MPI_Allgather(	pneighb_c,	6, MPI_INTEGER, &
      						gneighb_c,	6, MPI_INTEGER, &
      						mpi_comm_cart,ierr)

!      	if(taskid .eq. 0) write(*,*) 'gneighb: ',gneighb_r(:)

      	ghosts_set = .true.

      end subroutine p3dfft_init_ghosts

!     ------------------------------------------------------------------
!     !  p3dfft_update_ghosts(XYZ, goverlap)
!     !     # XYZ - array of real made out of two parts
!     !              1. part = input data = FFT-pencil
!     !                        [(fsize(1)*2) *(fsize(2)*2) *(fsize(3)*2)]
!     !              2. part = output data = ghost-cells
!     !                         [(fsize(1)*2) *(fsize(2)*2) *goverlap]*2	2* XY-plane
!     !                        +[(fsize(1)*2) *(fsize(3)*2) *goverlap]*2	2* XZ-plane
!     !                        +[(fsize(2)*2) *(fsize(3)*2) *goverlap]*2	2* YZ-plane
!     !     # gneighb    - neighbours of all pencils in all directions

      subroutine p3dfft_update_ghosts(XYZ, gneighb, gproc_size, &
                                     gmem_start, gslab_start, gslab_end)
      	implicit none

!     !	function args
      	real(mytype) ,TARGET :: XYZ(*)
      	integer, intent(in) :: gneighb(*)
      	integer, intent(in) :: gproc_size(*)
      	integer(i8), intent(in) :: gmem_start(*)
      	integer, intent(in) :: gslab_start(3,6)
      	integer, intent(in) :: gslab_end(3,6)

!     !	usual stuff
     	integer,save :: d,mp, dmp, pmd
     	integer,save :: send_neighb, recv_neighb
     	integer(i8),save :: i,j,k,jk
     	integer,save :: g, buf_size
     	integer(i8),save :: jpos,kpos
     	integer, save :: request, status(mpi_status_size), ierr

!     !	exchange ghost-cells sides of pencil/slide
      	do d=1, 3	! X,Y,Z sides
      		do mp=1,2
      			if(mp .eq. 1) then	! from - to +
      				send_neighb = (taskid*6) +(d-1)*2 +2
      				recv_neighb = (taskid*6) +(d-1)*2 +1
      				dmp = d
      				pmd = 3+d
      			else				! from + to -
      				send_neighb = (taskid*6) +(d-1)*2 +1
      				recv_neighb = (taskid*6) +(d-1)*2 +2
      				dmp = 3+d
      				pmd = d
      			endif

!      			dimension d is contained entirely within processors memory
      			if(gneighb(send_neighb) .eq. taskid) then
!      			   .and. gneighb(recv_neighb) .eq. taskid) should always be the case if gneighb(send_neighb)=taskid

!      				just copy data in local memory
      				g=0
      				do k=gslab_start(3,dmp), gslab_end(3,dmp)
      				  kpos = (k-1)*gproc_size(1)*gproc_size(2)
      				  do j=gslab_start(2,dmp), gslab_end(2,dmp)
      					jpos = kpos +(j-1)*gproc_size(1)
      					do i=gslab_start(1,dmp), gslab_end(1,dmp)
      					  XYZ(gmem_start(pmd) +g) = XYZ(jpos +i)
      					  g = g+1
      					enddo
      				  enddo
      				enddo

!     			dimension d is block-distributed among processors
      			else

!      				pack send-buffer
      				g=0
      				do k=gslab_start(3,dmp), gslab_end(3,dmp)
      				  kpos = (k-1)*gproc_size(1)*gproc_size(2)
      				  do j=gslab_start(2,dmp), gslab_end(2,dmp)
      					jpos = kpos +(j-1)*gproc_size(1)
      					do i=gslab_start(1,dmp), gslab_end(1,dmp)
      					  g = g+1
      					  gbuf_snd(g) = XYZ(jpos +i)
      					enddo
      				  enddo
      				enddo

!     				start a receive, send ghost-cells then wait
      				buf_size = g*mytype
      				call MPI_IRecv(	gbuf_recv, buf_size, MPI_BYTE, &
      						gneighb(send_neighb), d, mpi_comm_cart, &
      						request, ierr) 

      				call MPI_Send(	gbuf_snd, buf_size, MPI_BYTE, &
      						gneighb(recv_neighb), d, mpi_comm_cart, ierr)

      				call MPI_Wait( request, status, ierr )

!      				unpack recieve-buffer
      				do i=1, g
      					XYZ(gmem_start(pmd)+i-1) = gbuf_recv(i)
      				enddo

      			endif

      		enddo
      	enddo

      end subroutine p3dfft_update_ghosts

!     ------------------------------------------------------------------
!
!     	ijk2i() - calc 1d-array-position
!
      function gr_ijk2i(i,j,k)
      	implicit none

!     	function args
      	integer, intent(in) :: i,j,k
      	integer(i8) :: gr_ijk2i

      	if(.not. ghosts_set) then
      		print *,'P3DFFT error: call p3dfft_init_ghosts before'
      		return
      	endif

!     	ghost-cell in i-dim
      	if(i .lt. gproc_rstart(1)) then
      	  gr_ijk2i = gmem_rstart(1) &
      		+int( gslab_rsize(1,1)*gslab_rsize(2,1),8) *int( k-gproc_rstart(3),8) &
      		+int( gslab_rsize(1,1)                 ,8) *int( j-gproc_rstart(2),8) &
      		+                                           int( i-gproc_rstart(1) +gslab_rsize(1,1),8)
      		return
      	else if(i .gt. gproc_rend(1)  ) then
      	  gr_ijk2i = gmem_rstart(4) &
      		+int( gslab_rsize(1,1)*gslab_rsize(1,2),8) *int( k-gproc_rstart(3),8) &
      		+int( gslab_rsize(1,1)                 ,8) *int( j-gproc_rstart(2),8) &
      		+                                           int( i-gproc_rend(1)-1,8)
      		return

!     	ghost-cell in j-dim
      	else if(j .lt. gproc_rstart(2)) then
      	  gr_ijk2i = gmem_rstart(2) &
      		+int( gslab_rsize(1,2)*gslab_rsize(2,2),8) *int( k-gproc_rstart(3),8) &
      		+int( gslab_rsize(1,2)                 ,8) *int( j-gproc_rstart(2) +gslab_rsize(2,2),8) &
      		+                                           int( i-gproc_rstart(1),8)
      		return
      	else if(j .gt. gproc_rend(2)  ) then
      	  gr_ijk2i = gmem_rstart(5) &
      		+int( gslab_rsize(1,2)*gslab_rsize(2,2),8) *int( k-gproc_rstart(3),8) &
      		+int( gslab_rsize(1,2)                 ,8) *int( j-gproc_rend(2)-1,8) &
      		+                                           int( i-gproc_rstart(1),8)
      		return

!     	ghost-cell in k-dim
      	else if(k .lt. gproc_rstart(3)) then
      	  gr_ijk2i = gmem_rstart(3) &
      		+int( gslab_rsize(1,3)*gslab_rsize(2,3),8) *int( k-gproc_rstart(3) +gslab_rsize(3,3),8) &
      		+int( gslab_rsize(1,3)                 ,8) *int( j-gproc_rstart(2),8) &
      		+                                           int( i-gproc_rstart(1),8)
      	else if(k .gt. gproc_rend(3)  ) then
      	  gr_ijk2i = gmem_rstart(6) &
      		+int( gslab_rsize(1,3)*gslab_rsize(2,3),8) *int( k-gproc_rend(3)-1,8) &
      		+int( gslab_rsize(1,3)                 ,8) *int( j-gproc_rstart(2),8) &
      		+                                           int( i-gproc_rstart(1),8)
      		return

!     	no ghost-cell
      	else
      	  gr_ijk2i = 1 &
      		+int(gproc_rsize(1)*gproc_rsize(2),8) *int(k-gproc_rstart(3),8) &
      		+int(gproc_rsize(1)               ,8) *int(j-gproc_rstart(2),8) &
      		+                                      int(i-gproc_rstart(1),8)
      	endif

      end function gr_ijk2i

!     ------------------------------------------------------------------
      subroutine update_rghosts(XgYZ)
      	real(mytype) ,TARGET :: XgYZ(1,1,*)

      	if(.not. ghosts_set) then
      		print *,'P3DFFT error: call p3dfft_init_ghosts before'
      		return
      	endif
      	call p3dfft_update_ghosts(XgYZ, gneighb_r, gproc_rsize, &
                                   gmem_rstart, gslab_rstart, gslab_rend)

      end subroutine update_rghosts

!     ------------------------------------------------------------------
      subroutine update_cghosts(XYZg)
      	real(mytype) ,TARGET :: XYZg(1,1,*)

      	print *,'P3DFFT error: support for complex ghost-cells not implemented'
      	print *,'              in p3dfft_update_ghosts() and ijk2i() yet,'
      	print *,'              perhaps you like to do the job (?)'
      	return

      	if(.not. ghosts_set) then
      		print *,'P3DFFT error: call p3dfft_init_ghosts before'
      		return
      	endif
      	call p3dfft_update_ghosts(XYZg, gneighb_c, gproc_csize, &
                                   gmem_cstart, gslab_cstart, gslab_cend)

      end subroutine update_cghosts
