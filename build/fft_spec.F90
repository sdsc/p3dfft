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

      module fft_spec

      integer, parameter, public :: r8 = KIND(1.0d0)
      integer, parameter, public :: i8 = SELECTED_INT_KIND(18)

#ifdef FFTW
      include "fftw3.f"

      integer(i8), allocatable, dimension(:) :: plan1_frc,plan1_bcr,plan1_fc,plan1_bc
      integer(i8), allocatable, dimension(:) :: plan_ctrans_same, plan_strans_same,  plan_ctrans_dif, plan_strans_dif
      integer(i8), allocatable, dimension(:) :: plan2_bc_same,plan2_fc_same,plan2_bc_dif,plan2_fc_dif
      integer(i8), allocatable, dimension(:) :: startx_frc,startx_bcr,startx_f_c1,startx_b_c1
      integer(i8), allocatable, dimension(:) :: startx_ctrans_same, startx_strans_same,  startx_ctrans_dif, startx_strans_dif
      integer(i8), allocatable, dimension(:) :: startx_b_c2_same,startx_f_c2_same,startx_b_c2_dif,startx_f_c2_dif
      integer(i8), allocatable, dimension(:) :: starty_frc,starty_bcr,starty_f_c1,starty_b_c1
      integer(i8), allocatable, dimension(:) :: starty_ctrans_same, starty_strans_same, &
         starty_ctrans_dif, starty_strans_dif
      integer(i8), allocatable, dimension(:) :: starty_b_c2_same,starty_f_c2_same,starty_b_c2_dif,starty_f_c2_dif

      integer fftw_flag,NULL
#ifdef ESTIMATE
      parameter(fftw_flag = FFTW_ESTIMATE,NULL=0)
#elif defined PATIENT
      parameter(fftw_flag = FFTW_PATIENT,NULL=0)
#else
      parameter(fftw_flag = FFTW_MEASURE,NULL=0)
#endif

#endif

#ifdef ESSL
      integer :: cnaux,rnaux1,rnaux2
      real(r8),save,allocatable :: caux1(:),caux2(:),raux1(:),raux2(:)
      real(r8),save :: raux3(1)
#endif

      end module
