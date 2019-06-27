! This file is part of P3DFFT library

! Title: P3DFFT library

! Authors: Dmitry Pekurovsky

! Copyright (c) 2006-2019 

! The Regents of the University of California.

! All Rights Reserved.                        

 

!    Permission to use, copy, modify and  distribute  any part

!    of this software for  educational,  research  and  non-profit

!    purposes, by individuals or non-profit organizations,

!    without fee,  and  without a written  agreement is

!    hereby granted,  provided  that the  above  copyright notice,

!    this paragraph  and the following  three  paragraphs appear in

!    all copies.       

 

!    For-profit organizations desiring to use this software and others

!    wishing to incorporate this  software into commercial

!    products or use it for  commercial  purposes should contact the:    

!          Office of Innovation & Commercialization 

!          University of California San Diego

!          9500 Gilman Drive,  La Jolla,  California, 92093-0910        

!          Phone: (858) 534-5815

!          E-mail: innovation@ucsd.edu

 

!    IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE

!    TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR    

!    CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT

!    OF THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY OF

!    CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH

!    DAMAGE.

 

!    THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND

!    THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE        

!    MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. 

!    THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND    

!    EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR

!    IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES

!    OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR

!    THAT THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT,        

!    TRADEMARK OR OTHER RIGHTS.
!----------------------------------------------------------------------

! =========================================================
      subroutine p3dfft_setup_c(dims,nx,ny,nz,mpi_comm_in,nxcut,nycut,nzcut,OW,memsize) BIND(C,NAME='p3dfft_setup')
!========================================================

      use iso_c_binding
      implicit none

      integer nx,ny,nz,mpi_comm_in,dims(2)
      integer,intent (out) :: memsize (3)
      integer,intent (in) :: nxcut,nycut,nzcut
      integer,intent(in) :: OW
      logical overwrite

      if(OW .ne. 0) then
         overwrite = .true.
      else
         overwrite = .false.
      endif

      call p3dfft_setup(dims,nx,ny,nz,mpi_comm_in,nxcut,nycut,nzcut,overwrite,memsize)

      return
      end subroutine

! =========================================================
      subroutine p3dfft_setup(dims,nx,ny,nz,mpi_comm_in,nxcut,nycut,nzcut,overwrite,memsize)
!========================================================

      implicit none

      integer i,j,k,nx,ny,nz,err,mpi_comm_in
      integer ierr, dims(2),  cartid(2),mydims(2)
      logical periodic(2),remain_dims(2)
      integer impid, ippid, jmpid, jppid
      integer(i8) n1,n2,pad1,padd
      real(p3dfft_type), allocatable :: R(:)
      integer, optional, intent (out) :: memsize (3)
      integer, optional, intent (in) :: nxcut,nycut,nzcut
      logical, optional, intent(in) :: overwrite

      integer my_start (3), my_end (3), my_size (3)
      integer my_proc_dims (2, 9)

      if(nx .le. 0 .or. ny .le. 0 .or. nz .le. 0) then
         print *,'Invalid dimensions :',nx,ny,nz
         call MPI_ABORT(mpicomm, 0)
      endif

      if(mpi_set) then
         print *,'P3DFFT Setup error: the problem is already initialized. '
         print *,'Currently multiple setups not supported.'
         print *,'Quit the library using p3dfft_clean before initializing another setup'
         call MPI_ABORT(mpicomm, 0)
      endif

      if(present(overwrite)) then
         OW = overwrite
      else
         OW = .true.
      endif


      timers = 0.0

      mpi_set = .true.
      mpicomm = mpi_comm_in
      nx_fft = nx
      ny_fft = ny
      nz_fft = nz
      if(present(nxcut)) then
	nxc = nxcut
      else
        nxc = nx
      endif
      if(present(nycut)) then
	nyc = nycut
      else
        nyc = ny
      endif
      if(present(nzcut)) then
	nzc = nzcut
      else
        nzc = nz
      endif

      nxh=nx/2
      nxhp=nxh+1
	nxhc = nxc/2
	nxhpc = nxhc + 1
	nyh = ny/2
	nzh = nz/2
	nyhc = nyc / 2
	nzhc = nzc / 2
	nycph = (nyc + 1)/2
	nzcph = (nzc + 1)/2

      call MPI_COMM_SIZE (mpicomm,numtasks,ierr)
      call MPI_COMM_RANK (mpicomm,taskid,ierr)

      if(dims(1) .le. 0 .or. dims(2) .le. 0 .or.  dims(1)*dims(2) .ne. numtasks) then
         print *,'Invalid processor geometry: ',dims,' for ',numtasks, 'tasks'
         call MPI_ABORT(mpicomm, 0)
      endif

#ifdef STRIDE1
      if(taskid .eq. 0) then
         print *,'Using stride-1 layout'
      endif
#endif

      iproc = dims(1)
      jproc = dims(2)

#ifndef DIMS_C
!       i = dims(1)
!       dims(1) = dims(2)
!       dims(2) = i
	mydims(1) = dims(2)
	mydims(2) = dims(1)
#else
	mydims = dims
#endif

      periodic(1) = .false.
      periodic(2) = .false.
! creating cartesian processor grid
      call MPI_Cart_create(mpicomm,2,mydims,periodic,.false.,mpi_comm_cart,ierr)
! Obtaining process ids with in the cartesian grid
      call MPI_Cart_coords(mpi_comm_cart,taskid,2,cartid,ierr)
! process with a linear id of 5 may have cartid of (3,1)

#ifdef DIMS_C
      ipid = cartid(1)
      jpid = cartid(2)
#else
      ipid = cartid(2)
      jpid = cartid(1)
#endif

! store processor-grid-informations
    cartid (1) = ipid
    cartid (2) = jpid
    allocate (proc_id2coords(0:(iproc*jproc)*2-1))
    call MPI_Allgather (cartid, 2, MPI_INTEGER, proc_id2coords, 2, MPI_INTEGER, mpi_comm_cart, ierr)
    allocate (proc_coords2id(0:iproc-1, 0:jproc-1))
    do i = 0, (iproc*jproc) - 1
      proc_coords2id (proc_id2coords(2*i), &
                      proc_id2coords(2*i+1)) = i
    end do

! here i is east-west j is north-south
! impid is west neighbour ippid is east neighbour and so on
      impid = ipid - 1
      ippid = ipid + 1
      jmpid = jpid - 1
      jppid = jpid + 1
!boundary processes
      if (ipid.eq.0) impid = MPI_PROC_NULL
      if (jpid.eq.0) jmpid = MPI_PROC_NULL
      if (ipid.eq.iproc-1) ippid = MPI_PROC_NULL
      if (jpid.eq.jproc-1) jppid = MPI_PROC_NULL
! using cart comworld create east-west(row) sub comworld

#ifdef DIMS_C
      remain_dims(1) = .true.
      remain_dims(2) = .false.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_row,ierr)
! using cart comworld create north-south(column) sub comworld
      remain_dims(1) = .false.
      remain_dims(2) = .true.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_col,ierr)
#else
      remain_dims(1) = .true.
      remain_dims(2) = .false.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_col,ierr)
! using cart comworld create north-south(column) sub comworld
      remain_dims(1) = .false.
      remain_dims(2) = .true.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_row,ierr)
#endif

      allocate (iist(0:iproc-1))
      allocate (iisz(0:iproc-1))
      allocate (iien(0:iproc-1))
      allocate (jjst(0:jproc-1))
      allocate (jjsz(0:jproc-1))
      allocate (jjen(0:jproc-1))
      allocate (jist(0:iproc-1))
      allocate (jisz(0:iproc-1))
      allocate (jien(0:iproc-1))
      allocate (kjst(0:jproc-1))
      allocate (kjsz(0:jproc-1))
      allocate (kjen(0:jproc-1))
!
!Mapping 3-D data arrays onto 2-D process grid
! (nx+2,ny,nz) => (iproc,jproc)
!
      call MapDataToProc(nxhpc,iproc,iist,iien,iisz)
      call MapDataToProc(ny,iproc,jist,jien,jisz)
      call MapDataToProc(nyc,jproc,jjst,jjen,jjsz)
      call MapDataToProc(nz,jproc,kjst,kjen,kjsz)

! These are local array indices for each processor

      iistart = iist(ipid)
      jjstart = jjst(jpid)
      jistart = jist(ipid)
      kjstart = kjst(jpid)
      iisize= iisz(ipid)
      jjsize= jjsz(jpid)
      jisize= jisz(ipid)
      kjsize= kjsz(jpid)
      iiend = iien(ipid)
      jjend = jjen(jpid)
      jiend = jien(ipid)
      kjend = kjen(jpid)

    allocate (iiist(0:iproc-1))
    allocate (iiisz(0:iproc-1))
    allocate (iiien(0:iproc-1))
    allocate (ijst(0:jproc-1))
    allocate (ijsz(0:jproc-1))
    allocate (ijen(0:jproc-1))
    call MapDataToProc (nx, iproc, iiist, iiien, iiisz)
    call MapDataToProc (nx, jproc, ijst, ijen, ijsz)
    iiistart = iiist (ipid)
    iiisize = iiisz (ipid)
    iiiend = iiien (ipid)
    ijstart = ijst (jpid)
    ijsize = ijsz (jpid)
    ijend = ijen (jpid)

#ifdef USE_EVEN

      IfCntMax = iisz(iproc-1)*jisz(iproc-1)*kjsize*p3dfft_type*2
      KfCntMax = iisize * jjsz(jproc-1) * kjsz(jproc-1)*p3dfft_type*2
      if(mod(ny,jproc) .ne. 0 .or. mod(nz,jproc) .ne. 0) then
         KfCntUneven = .true.
      else
         KfCntUneven = .false.
      endif
    IiCntMax = iiisz (iproc-1) * jisize * kjsize * p3dfft_type
    IJCntMax = ijsz (jproc-1) * jisize * kjsize * p3dfft_type
    JICntMax = jisz (iproc-1) * iiisize * kjsize * p3dfft_type
    KjCntMax = kjsz (jproc-1) * jisize * ijsize * p3dfft_type
#endif


#ifdef STRIDE1
#ifdef CACHE_BL
      CB = CACHE_BL
#else
      CB = 32768
#endif

#ifdef NBL_X
      NBx = NBL_X
#else
      NBx=CB/(4*p3dfft_type*Ny)
#endif

#ifdef NBL_Y1
      NBy1=NBL_Y1
#else
      NBy1 = CB/(4*p3dfft_type*iisize)
#endif

#ifdef NBL_Y2
      NBy2=NBL_Y2
#else
      NBy2 = CB/(4*p3dfft_type*Nz)
#endif

#ifdef NBL_Z
      NBz=NBL_Z
#else
      NBz = (CB/Nx)*(numtasks/(2*p3dfft_type*Ny))
#endif

      if(NBx .eq. 0) then
         NBx = 1
      endif
      if(NBy1 .eq. 0) then
         NBy1 = 1
      endif
      if(NBy2 .eq. 0) then
         NBy2 = 1
      endif
      if(NBz .eq. 0) then
         NBz = 1
      endif

      if(taskid .eq. 0) then
         print *,'Using loop block sizes ',NBx,NBy1,NBy2,NBz
      endif

#endif


! We may need to pad arrays due to uneven size
      padd = max(iisize*jjsize*nz_fft,iisize*ny_fft*kjsize) - nxhp*jisize*kjsize
      padi = padd
      if(padi .le. 0) then
         padi=0
      else
         if(mod(padi,nxhp*jisize) .eq. 0) then
            padi = padi / (nxhp*jisize)
         else
            padi = padi / (nxhp*jisize)+1
         endif

      endif

!      print *,'padi=',padi

! Initialize FFTW and allocate buffers for communication
      nm = nxhp * jisize * (kjsize+padi)
      nv_preset = 1
      if(nm .gt. 0) then
        allocate(buf1(nm),stat=err)
        if(err .ne. 0) then
           print *,'p3dfft_setup: Error allocating buf1 (',nm
        endif
        allocate(buf2(nm),stat=err)
        if(err .ne. 0) then
           print *,'p3dfft_setup: Error allocating buf2 (',nm
        endif
1        allocate(R(nm*2),stat=err)
!        if(err .ne. 0) then
!           print *,'p3dfft_setup: Error allocating R (',nm*2
!        endif
        buf1 = 0.0
        R = 0.0

! For FFT libraries that allocate work space implicitly such as through
! plans (e.g. FFTW) initialize here

        call init_plan
!(buf1,R,buf2,nm)

!        deallocate(R)
        allocate(buf(nm),stat=err)
        if(err .ne. 0) then
           print *,'p3dfft_setup: Error allocating buf (',nm
        endif


     endif

#ifdef USE_EVEN
      n1 = IfCntMax * iproc /(p3dfft_type*2)
      n2 = KfCntMax * jproc / (p3dfft_type*2)
      n1 = max(n1,n2)
      if(n1 .gt. nm) then
         deallocate(buf1)
         allocate(buf1(n1))
         deallocate(buf2)
         allocate(buf2(n1))
      endif
#endif

!#ifdef USE_EVEN
!
!     n1 = IfCntMax * iproc / (p3dfft_type*2)
!     n2 = KfCntMax * jproc / (p3dfft_type*2)
!     allocate(buf_x(n1))
!     allocate(buf_z(n2))
!     n1 = max(n1,n2)
!     allocate(buf_y(n1))
!#else
!     allocate(buf_x(nxhp*jisize*kjsize))
!     allocate(buf_y(ny*iisize*kjsize))
!     allocate(buf_z(nz*iisize*jjsize))
!#endif

! Displacements and buffer counts for mpi_alltoallv

      allocate (IfSndStrt(0:iproc-1))
      allocate (IfSndCnts(0:iproc-1))
      allocate (IfRcvStrt(0:iproc-1))
      allocate (IfRcvCnts(0:iproc-1))

      allocate (KfSndStrt(0:jproc-1))
      allocate (KfSndCnts(0:jproc-1))
      allocate (KfRcvStrt(0:jproc-1))
      allocate (KfRcvCnts(0:jproc-1))

      allocate (JrSndStrt(0:jproc-1))
      allocate (JrSndCnts(0:jproc-1))
      allocate (JrRcvStrt(0:jproc-1))
      allocate (JrRcvCnts(0:jproc-1))

      allocate (KrSndStrt(0:iproc-1))
      allocate (KrSndCnts(0:iproc-1))
      allocate (KrRcvStrt(0:iproc-1))
      allocate (KrRcvCnts(0:iproc-1))


!   start pointers and types of send  for the 1st forward transpose
      do i=0,iproc-1
         IfSndStrt(i) = (iist(i) -1)* jisize*kjsize*p3dfft_type*2
         IfSndCnts(i) = iisz(i) * jisize*kjsize*p3dfft_type*2

!   start pointers and types of recv for the 1st forward transpose
         IfRcvStrt(i) = (jist(i) -1) * iisize*kjsize*p3dfft_type*2
         IfRcvCnts(i) = jisz(i) * iisize*kjsize*p3dfft_type*2
      end do

!   start pointers and types of send  for the 2nd forward transpose
      do i=0,jproc-1
         KfSndStrt(i) = (jjst(i) -1)*iisize*kjsize*p3dfft_type*2
         KfSndCnts(i) = iisize*kjsize*jjsz(i)*p3dfft_type*2

!   start pointers and types of recv for the 2nd forward transpose
         KfRcvStrt(i) = (kjst(i) -1) * iisize * jjsize*p3dfft_type*2
         KfRcvCnts(i) = iisize*jjsize*kjsz(i)*p3dfft_type*2
      end do

!   start pointers and types of send  for the 1st inverse transpose
      do i=0,jproc-1
         JrSndStrt(i) = (kjst(i) -1) * iisize * jjsize*p3dfft_type*2
         JrSndCnts(i) = iisize*jjsize*kjsz(i)*p3dfft_type*2

!   start pointers and types of recv for the 1st inverse transpose
         JrRcvStrt(i) = (jjst(i) -1)*iisize*kjsize*p3dfft_type*2
         JrRcvCnts(i) = jjsz(i) * iisize * kjsize*p3dfft_type*2
      end do

!   start pointers and types of send  for the 2nd inverse transpose
      do i=0,iproc-1
         KrSndStrt(i) = (jist(i) -1) * iisize*kjsize*p3dfft_type*2
         KrSndCnts(i) = jisz(i) * iisize*kjsize*p3dfft_type*2

!   start pointers and types of recv for the 2nd inverse transpose
         KrRcvStrt(i) = (iist(i) -1) * jisize*kjsize*p3dfft_type*2
         KrRcvCnts(i) = jisize*iisz(i)*kjsize*p3dfft_type*2
      enddo

! Displacements and buffer counts for mpi_alltoallv in transpose-functions(..)
    allocate (IiStrt(0:iproc-1))
    allocate (IiCnts(0:iproc-1))
    allocate (JiStrt(0:iproc-1))
    allocate (JiCnts(0:iproc-1))

    allocate (IjStrt(0:jproc-1))
    allocate (IjCnts(0:jproc-1))
    allocate (KjStrt(0:jproc-1))
    allocate (KjCnts(0:jproc-1))

!   start pointers and size for the x<->y transpose
    do i = 0, iproc - 1
!        x->y
      IiStrt (i) = (iiist(i)-1) * jisize * kjsize * p3dfft_type
      IiCnts (i) = iiisz (i) * jisize * kjsize * p3dfft_type

      JiStrt (i) = (jist(i)-1) * iiisize * kjsize * p3dfft_type
      JiCnts (i) = jisz (i) * iiisize * kjsize * p3dfft_type
    end do

!   start pointers and size for the x<->z transpose
    do i = 0, jproc - 1
!        x->z
      IjStrt (i) = (ijst(i)-1) * jisize * kjsize * p3dfft_type
      IjCnts (i) = ijsz (i) * jisize * kjsize * p3dfft_type

      KjStrt (i) = (kjst(i)-1) * jisize * ijsize * p3dfft_type
      KjCnts (i) = kjsz (i) * jisize * ijsize * p3dfft_type
    end do

!   create proc_dims = send information of each pencil-dimensions to all procs
    call p3dfft_get_dims (my_start, my_end, my_size, 1)
    my_proc_dims (1, 1) = my_start (1)
    my_proc_dims (1, 2) = my_start (2)
    my_proc_dims (1, 3) = my_start (3)
    my_proc_dims (1, 4) = my_end (1)
    my_proc_dims (1, 5) = my_end (2)
    my_proc_dims (1, 6) = my_end (3)
    my_proc_dims (1, 7) = my_size (1)
    my_proc_dims (1, 8) = my_size (2)
    my_proc_dims (1, 9) = my_size (3)
    call p3dfft_get_dims (my_start, my_end, my_size, 2)
    my_proc_dims (2, 1) = my_start (1)
    my_proc_dims (2, 2) = my_start (2)
    my_proc_dims (2, 3) = my_start (3)
    my_proc_dims (2, 4) = my_end (1)
    my_proc_dims (2, 5) = my_end (2)
    my_proc_dims (2, 6) = my_end (3)
    my_proc_dims (2, 7) = my_size (1)
    my_proc_dims (2, 8) = my_size (2)
    my_proc_dims (2, 9) = my_size (3)

    allocate (proc_dims(2, 9, 0:(iproc*jproc)-1))
    call MPI_Allgather (my_proc_dims, 2*9, MPI_INTEGER,proc_dims, 2 * 9, MPI_INTEGER,mpi_comm_cart, ierr)

    allocate (proc_parts((iproc*jproc), 7))
    proc_parts = - 1

!     calc max. needed memory (attention: cast to integer8 included)
      pad1 = 2* max(nz*jjsize*iisize,ny*kjsize*iisize) - nx*jisize*kjsize
!      print *,taskid,': pad1=',pad1

      if(pad1 .le. 0) then
        pad1 = 0
      endif

      padi=pad1
      if(mod(padi,nx*jisize) .ne. 0) then
         padi = padi / (nx*jisize) + 1
      else
         padi = padi / (nx*jisize)
      endif


	maxisize = nx
	maxjsize = jisize
	maxksize = kjsize + padi

	if(present(memsize)) then
	  memsize(1) = maxisize
	  memsize(2) = maxjsize
	  memsize(3) = maxksize
	endif

      end subroutine p3dfft_setup

!==================================================================
      subroutine MapDataToProc (data,proc,st,en,sz)
!========================================================
!
       implicit none
       integer data,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
       integer i,size,nl,nu

       size=data/proc
       nu = data - size * proc
       nl = proc - nu
       st(0) = 1
       sz(0) = size
       en(0) = size
       do i=1,nl-1
         st(i) = st(i-1) + size
         sz(i) = size
         en(i) = en(i-1) + size
      enddo
      size = size + 1
      do i=nl,proc-1
         st(i) = en(i-1) + 1
         sz(i) = size
         en(i) = en(i-1) + size
      enddo
      en(proc-1)= data
      sz(proc-1)= data-st(proc-1)+1

      end subroutine

