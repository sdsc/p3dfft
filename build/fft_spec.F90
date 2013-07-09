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

      module fft_spec

      integer, parameter, public :: r8 = KIND(1.0d0)
      integer, parameter, public :: i8 = SELECTED_INT_KIND(18)

#ifdef FFTW
      include "fftw3.f"
      integer(i8) plan1_frc,plan1_bcr,plan1_fc,plan1_bc
      integer(i8) plan_ctrans_same, plan_strans_same, &
         plan_ctrans_dif, plan_strans_dif
      integer(i8) plan2_bc_same,plan2_fc_same,plan2_bc_dif,plan2_fc_dif
!      integer(i8) plan1,plan2,plan3      
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
