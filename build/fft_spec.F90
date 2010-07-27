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

      module fft_spec

#ifdef FFTW
      include "fftw3.f"
      integer(SELECTED_INT_KIND(16)) plan1_frc,plan1_bcr,plan1_fc,plan2_fc,plan1_bc,plan2_bc
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
