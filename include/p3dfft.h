/*
! This file is part of P3DFFT library
!
! Version 2.7
!
! Copyright (C) 2006-2014 Dmitry Pekurovsky
! Copyright (C) 2006-2014 UC San Diego
!
! P3DFFT is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! P3DFFT is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with P3DFFT.  If not, see <http://www.gnu.org/licenses/>.
!
*/

#include <stdlib.h>

#define FORT_MOD_NAME(NAME) NAME

#ifdef IBM

#define FORTNAME(NAME) NAME

#elif defined INTEL

#define FORTNAME(NAME) NAME##_

#elif defined PGI 

/* #define FORT_MOD_NAME(NAME) p3dfft_##NAME##_ */

#define FORTNAME(NAME) NAME##_

#elif defined CRAY

#define FORTNAME(NAME) NAME

#elif defined GNU
#include "gnu.h"
//#define FORTNAME(NAME) NAME

#else
#define FORTNAME(NAME) NAME

#endif


extern void FORT_MOD_NAME(p3dfft_setup)(int *dims,int *nx,int *ny,int *nz, int * comm, int *nxc, int *nyc, int *nzc, int *ow, int *memsize);
extern void FORT_MOD_NAME(p3dfft_get_dims)(int *,int *,int *,int *);
extern void FORT_MOD_NAME(get_timers)(double *timers);
extern void FORT_MOD_NAME(set_timers)();

#ifndef SINGLE_PREC
extern void FORT_MOD_NAME(p3dfft_ftran_r2c)(double *A,double *B, unsigned char *op);
extern void FORT_MOD_NAME(p3dfft_btran_c2r)(double *A,double *B, unsigned char *op);
#else
extern void FORT_MOD_NAME(p3dfft_ftran_r2c)(float *A,float *B, unsigned char *op);
extern void FORT_MOD_NAME(p3dfft_btran_c2r)(float *A,float *B, unsigned char *op);
#endif

extern void FORT_MOD_NAME(p3dfft_clean)();


extern void Cp3dfft_setup(int *dims,int nx,int ny,int nz,int comm, int nxc,int nyc, int nzc, int ovewrite, int *memsize);
extern void Cp3dfft_clean();

extern void Cp3dfft_get_dims(int *,int *,int *,int );

#ifndef SINGLE_PREC
extern void Cp3dfft_ftran_r2c(double *A,double *B, unsigned char *op);
extern void Cp3dfft_btran_c2r(double *A,double *B, unsigned char *op);
#else
extern void Cp3dfft_ftran_r2c(float *A,float *B, unsigned char *op);
extern void Cp3dfft_btran_c2r(float *A,float *B, unsigned char *op);
#endif

extern void Cget_timers(double *timers);
extern void Cset_timers();


inline void Cp3dfft_setup(int *dims,int nx,int ny,int nz, int comm, int nxc, int nyc, int nzc, int overwrite, int * memsize)
{
  FORT_MOD_NAME(p3dfft_setup)(dims,&nx,&ny,&nz,&comm, &nxc, &nyc, &nzc, &overwrite,memsize);
}

inline void Cp3dfft_clean()
{
  FORT_MOD_NAME(p3dfft_clean)();
}

inline void Cp3dfft_get_dims(int *start,int *end,int *size,int conf)
{
  FORT_MOD_NAME(p3dfft_get_dims)(start,end,size,&conf);
}

inline void Cget_timers(double *timers)
{
  FORT_MOD_NAME(get_timers)(timers);
}

inline void Cset_timers()
{
  FORT_MOD_NAME(set_timers)();
}


#ifndef SINGLE_PREC
inline void Cp3dfft_ftran_r2c(double *A,double *B, unsigned char *op)
{
  FORT_MOD_NAME(p3dfft_ftran_r2c)(A,B,op);
}
#else
inline void Cp3dfft_ftran_r2c(float *A,float *B, unsigned char *op)
{
  FORT_MOD_NAME(p3dfft_ftran_r2c)(A,B,op);
}

#endif

#ifndef SINGLE_PREC
inline void Cp3dfft_btran_c2r(double *A,double *B, unsigned char *op)
{
  FORT_MOD_NAME(p3dfft_btran_c2r)(A,B,op);
}
#else
inline void Cp3dfft_btran_c2r(float *A,float *B, unsigned char *op)
{
  FORT_MOD_NAME(p3dfft_btran_c2r)(A,B,op);
}
#endif

