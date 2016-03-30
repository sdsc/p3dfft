! This file is part of P3DFFT library
!
!    P3DFFT
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2006-2014 Dmitry Pekurovsky
!    Copyright (C) 2006-2014 University of California
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

      subroutine ftran_r2c_many(IN,dim_in,OUT,dim_out,nv,op)

      use p3dfft

      integer dim_in,dim_out
      real(p3dfft_type), TARGET :: IN(dim_in,nv)
      complex(p3dfft_type), TARGET :: OUT(dim_out,nv)
!      real(p3dfft_type) IN(nx_fft,jistart:jiend,kjstart:kjend)
!      complex(p3dfft_type) OUT(iistart:iiend,jjstart:jjend,nz_fft)
!      logical flg_inplace
      character(len=3) op
      integer nv

      call p3dfft_ftran_r2c_many(IN,dim_in,OUT,dim_out,nv,op)

      return
      end

      subroutine btran_c2r_many(IN,dim_in,OUT,dim_out,nv,op)

      use p3dfft

      integer dim_in,dim_out
      real(p3dfft_type), TARGET :: OUT(dim_out,nv)
      complex(p3dfft_type), TARGET :: IN(dim_in,nv)
!      real(p3dfft_type) OUT(nx_fft,jistart:jiend,kjstart:kjend)
!      complex(p3dfft_type) IN(iistart:iiend,jjstart:jjend,nz_fft)
!      logical flg_inplace
      character(len=3) op
      integer nv

      call p3dfft_btran_c2r_many(IN,dim_in,OUT,dim_out,nv,op)

      return
      end

      subroutine ftran_r2c(IN,OUT,op)

      use p3dfft

      real(p3dfft_type), TARGET :: IN(1,1,*)
      complex(p3dfft_type), TARGET :: OUT(1,1,*)
!      real(p3dfft_type) IN(nx_fft,jistart:jiend,kjstart:kjend)
!      complex(p3dfft_type) OUT(iistart:iiend,jjstart:jjend,nz_fft)
!      logical flg_inplace
      character(len=3) op
      integer nv

      call p3dfft_ftran_r2c(IN,OUT,op)

      return
      end

      subroutine btran_c2r(IN,OUT,op)

      use p3dfft

      real(p3dfft_type), TARGET :: OUT(1,1,*)
      complex(p3dfft_type), TARGET :: IN(1,1,*)
!      real(p3dfft_type) OUT(nx_fft,jistart:jiend,kjstart:kjend)
!      complex(p3dfft_type) IN(iistart:iiend,jjstart:jjend,nz_fft)
!      logical flg_inplace
      character(len=3) op
      integer nv

      call p3dfft_btran_c2r(IN,OUT,op)

      return
      end
