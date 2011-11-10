! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2010 Dmitry Pekurovsky
!    Copyright (C) 2006-2010 University of California
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

! =========================================================
      subroutine p3dfft_setup(dims,nx,ny,nz,overwrite,memsize)
!========================================================

      implicit none

      integer i,j,k,nx,ny,nz,err
      integer ierr, dims(2),  cartid(2)
      logical periodic(2),remain_dims(2),overwrite
      integer impid, ippid, jmpid, jppid
      integer(i8) nm,n1,n2
      real(mytype), allocatable :: R(:)
      integer, optional, intent (out) :: memsize (3)

      integer my_start (3), my_end (3), my_size (3)
      integer my_proc_dims (2, 9)

      if(nx .le. 0 .or. ny .le. 0 .or. nz .le. 0) then
         print *,'Invalid dimensions :',nx,ny,nz
         call MPI_ABORT(MPI_COMM_WORLD, 0)
      endif

      if(mpi_set) then
         print *,'P3DFFT Setup error: the problem is already initialized. '
         print *,'Currently multiple setups not supported.'       
         print *,'Quit the library using p3dfft_clean before initializing another setup'
         call MPI_ABORT(MPI_COMM_WORLD, 0)
      endif

      OW = overwrite

      timers = 0.0

      mpi_set = .true.
      nx_fft = nx
      ny_fft = ny
      nz_fft = nz
      nxh=nx/2
      nxhp=nxh+1

      call MPI_COMM_SIZE (MPI_COMM_WORLD,numtasks,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,taskid,ierr)

      if(dims(1) .le. 0 .or. dims(2) .le. 0 .or.  dims(1)*dims(2) .ne. numtasks) then
         print *,'Invalid processor geometry: ',dims,' for ',numtasks, 'tasks'
         call MPI_ABORT(MPI_COMM_WORLD, 0)
      endif

#ifdef STRIDE1      
      if(taskid .eq. 0) then 
         print *,'Using stride-1 layout'
      endif
#endif

      iproc = dims(1)
      jproc = dims(2)

#ifndef DIMS_C
       i = dims(1)  
       dims(1) = dims(2)
       dims(2) = i
#endif

      periodic(1) = .false.
      periodic(2) = .false.
! creating cartesian processor grid
      call MPI_Cart_create(MPI_COMM_WORLD,2,dims,periodic,.false.,mpi_comm_cart,ierr)
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
      call MapDataToProc(nxhp,iproc,iist,iien,iisz)
      call MapDataToProc(ny,iproc,jist,jien,jisz)
      call MapDataToProc(ny,jproc,jjst,jjen,jjsz)
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

      IfCntMax = iisz(iproc-1)*jisz(iproc-1)*kjsize*mytype*2
      KfCntMax = iisize * jjsz(jproc-1) * kjsz(jproc-1)*mytype*2
      if(mod(ny,jproc) .ne. 0 .or. mod(nz,jproc) .ne. 0) then 
         KfCntUneven = .true.
      else
         KfCntUneven = .false.
      endif
    IiCntMax = iiisz (iproc-1) * jisize * kjsize * mytype
    IJCntMax = ijsz (iproc-1) * jisize * kjsize * mytype
    JICntMax = jisz (iproc-1) * iiisize * kjsize * mytype
    KjCntMax = kjsz (iproc-1) * jisize * ijsize * mytype
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
      NBx=CB/(4*mytype*Ny)
#endif

#ifdef NBL_Y1
      NBy1=NBL_Y1
#else
      NBy1 = CB/(4*mytype*iisize)
#endif

#ifdef NBL_Y2
      NBy2=NBL_Y2
#else
      NBy2 = CB/(4*mytype*Nz)
#endif

#ifdef NBL_Z
      NBz=NBL_Z
#else
      NBz = (CB/Nx)*(numtasks/(2*mytype*Ny))
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
      if(padd .le. 0) then 
         padd=0
      else
         if(mod(padd,nxhp*jisize) .eq. 0) then
            padd = padd / (nxhp*jisize)
         else
            padd = padd / (nxhp*jisize)+1
         endif

      endif

!      print *,taskid,': padd=',padd
! Initialize FFTW and allocate buffers for communication
      nm = nxhp * jisize * (kjsize+padd) 
      if(nm .gt. 0) then       
        allocate(buf1(nm),stat=err)
        if(err .ne. 0) then
           print *,'p3dfft_setup: Error allocating buf1 (',nm
        endif
        allocate(buf2(nm),stat=err)
        if(err .ne. 0) then
           print *,'p3dfft_setup: Error allocating buf2 (',nm
        endif
        allocate(R(nm*2),stat=err)
        if(err .ne. 0) then
           print *,'p3dfft_setup: Error allocating R (',nm*2
        endif
        buf1 = 0.0
        R = 0.0

! For FFT libraries that allocate work space implicitly such as through 
! plans (e.g. FFTW) initialize here

        call init_plan(buf1,R,buf2,nm)

        deallocate(R)
     endif

#ifdef USE_EVEN
      n1 = IfCntMax * iproc /(mytype*2)
      n2 = KfCntMax * jproc / (mytype*2)
      n1 = max(n1,n2)
      if(n1 .gt. nm) then
         deallocate(buf1)
         allocate(buf1(n1))
         deallocate(buf2)
         allocate(buf2(n1))
      endif
#endif

!     preallocate memory for FFT-Transforms
    allocate (buf(nxhp*jisize*(kjsize+padd)), stat=err)
!     initialize buf to avoid "floating point invalid" errors in debug mode
    buf = 0.d0
    if (err /= 0) then
      print *, 'Error ', err, ' allocating array XYgZ'
    end if


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
         IfSndStrt(i) = (iist(i) -1)* jisize*kjsize*mytype*2
         IfSndCnts(i) = iisz(i) * jisize*kjsize*mytype*2

!   start pointers and types of recv for the 1st forward transpose
         IfRcvStrt(i) = (jist(i) -1) * iisize*kjsize*mytype*2
         IfRcvCnts(i) = jisz(i) * iisize*kjsize*mytype*2
      end do

!   start pointers and types of send  for the 2nd forward transpose
      do i=0,jproc-1
         KfSndStrt(i) = (jjst(i) -1)*iisize*kjsize*mytype*2
         KfSndCnts(i) = iisize*kjsize*jjsz(i)*mytype*2

!   start pointers and types of recv for the 2nd forward transpose
         KfRcvStrt(i) = (kjst(i) -1) * iisize * jjsize*mytype*2
         KfRcvCnts(i) = iisize*jjsize*kjsz(i)*mytype*2
      end do

!   start pointers and types of send  for the 1st inverse transpose
      do i=0,jproc-1
         JrSndStrt(i) = (kjst(i) -1) * iisize * jjsize*mytype*2
         JrSndCnts(i) = iisize*jjsize*kjsz(i)*mytype*2

!   start pointers and types of recv for the 1st inverse transpose
         JrRcvStrt(i) = (jjst(i) -1)*iisize*kjsize*mytype*2
         JrRcvCnts(i) = jjsz(i) * iisize * kjsize*mytype*2
      end do

!   start pointers and types of send  for the 2nd inverse transpose
      do i=0,iproc-1
         KrSndStrt(i) = (jist(i) -1) * iisize*kjsize*mytype*2
         KrSndCnts(i) = jisz(i) * iisize*kjsize*mytype*2

!   start pointers and types of recv for the 2nd inverse transpose
         KrRcvStrt(i) = (iist(i) -1) * jisize*kjsize*mytype*2
         KrRcvCnts(i) = jisize*iisz(i)*kjsize*mytype*2
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
      IiStrt (i) = (iiist(i)-1) * jisize * kjsize * mytype
      IiCnts (i) = iiisz (i) * jisize * kjsize * mytype

      JiStrt (i) = (jist(i)-1) * iiisize * kjsize * mytype
      JiCnts (i) = jisz (i) * iiisize * kjsize * mytype
    end do

!   start pointers and size for the x<->z transpose
    do i = 0, jproc - 1
!        x->z
      IjStrt (i) = (ijst(i)-1) * jisize * kjsize * mytype
      IjCnts (i) = ijsz (i) * jisize * kjsize * mytype

      KjStrt (i) = (kjst(i)-1) * jisize * ijsize * mytype
      KjCnts (i) = kjsz (i) * jisize * ijsize * mytype
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
	maxisize = 2 * nxhp
	maxjsize = jisize
	maxksize = kjsize + padd

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

