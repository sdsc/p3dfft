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
