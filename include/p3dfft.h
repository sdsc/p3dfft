/*
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

#ifdef __cplusplus
extern "C"
{
#endif

extern void FORT_MOD_NAME(p3dfft_setup)(int *dims,int *nx,int *ny,int *nz, int * comm, int *nxc, int *nyc, int *nzc, int *ow, int *memsize);
extern void FORT_MOD_NAME(p3dfft_get_dims)(int *,int *,int *,int *);
extern void FORT_MOD_NAME(get_timers)(double *timers);
extern void FORT_MOD_NAME(set_timers)();

#ifndef SINGLE_PREC
extern void FORT_MOD_NAME(p3dfft_ftran_r2c)(double *A,double *B, unsigned char *op);
extern void FORT_MOD_NAME(p3dfft_ftran_r2c_many)(double *A, int *dim_in, double *B, int *dim_out, int *nv, unsigned char *op);
extern void FORT_MOD_NAME(p3dfft_btran_c2r)(double *A,double *B, unsigned char *op);
extern void FORT_MOD_NAME(p3dfft_btran_c2r_many)(double *A, int *dim_in, double *B, int *dim_out, int *nv, unsigned char *op);
extern void FORT_MOD_NAME(p3dfft_cheby)(double *A,double *B, double *Lz);
extern void FORT_MOD_NAME(p3dfft_cheby_many)(double *A, int *dim_in, double *B, int *dim_out, int *nv, double *Lz);
#else
extern void FORT_MOD_NAME(p3dfft_ftran_r2c)(float *A,float *B, unsigned char *op);
extern void FORT_MOD_NAME(p3dfft_ftran_r2c_many)(float *A, int *dim_in, float *B, int *dim_out, int *nv, unsigned char *op);
extern void FORT_MOD_NAME(p3dfft_btran_c2r)(float *A,float *B, unsigned char *op);
extern void FORT_MOD_NAME(p3dfft_btran_c2r_many)(float *A, int *dim_in, float *B, int *dim_out, int *nv, unsigned char *op);
extern void FORT_MOD_NAME(p3dfft_cheby)(float *A,float *B, float *Lz);
extern void FORT_MOD_NAME(p3dfft_cheby_many)(float *A, int *dim_in, float *B, int *dim_out, int *nv, float *Lz);
#endif

extern void FORT_MOD_NAME(p3dfft_clean)();


extern void Cp3dfft_setup(int *dims,int nx,int ny,int nz,int comm, int nxc,int nyc, int nzc, int ovewrite, int *memsize);
extern void Cp3dfft_clean();

extern void Cp3dfft_get_dims(int *,int *,int *,int );

#ifndef SINGLE_PREC
extern void Cp3dfft_ftran_r2c(double *A,double *B, unsigned char *op);
extern void Cp3dfft_btran_c2r(double *A,double *B, unsigned char *op);
extern void Cp3dfft_ftran_r2c_many(double *A,int dim_in,double *B, int dim_out,int nv, unsigned char *op);
extern void Cp3dfft_btran_c2r_many(double *A,int dim_in,double *B, int dim_out,int nv, unsigned char *op);
#else
extern void Cp3dfft_ftran_r2c(float *A,float *B, unsigned char *op);
extern void Cp3dfft_btran_c2r(float *A,float *B, unsigned char *op);
  extern void Cp3dfft_ftran_r2c_many(float *A,int dim_in,float *B, int dim_out,int nv,unsigned char *op);
  extern void Cp3dfft_btran_c2r_many(float *A,int dim_in,float *B, int dim_out,int nv, unsigned char *op);
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


#ifndef SINGLE_PREC
inline void Cp3dfft_cheby(double *A,double *B, double Lz)
{
  FORT_MOD_NAME(p3dfft_cheby)(A,B,&Lz);
}
#else
inline void Cp3dfft_cheby(float *A,float *B, float Lz)
{
  FORT_MOD_NAME(p3dfft_cheby)(A,B,&Lz);
}
#endif

#ifndef SINGLE_PREC
  inline void Cp3dfft_cheby_many(double *A,int dim_in,double *B, int dim_out, int nv, double Lz)
{
  FORT_MOD_NAME(p3dfft_cheby_many)(A,&dim_in,B,&dim_out,&nv,&Lz);
}
#else
  inline void Cp3dfft_cheby_many(float *A, int dim_in, float *B, int dim_out, int nv, float Lz)
{
  FORT_MOD_NAME(p3dfft_cheby_many)(A,&dim_in,B,&dim_out,&nv,&Lz);
}
#endif

#ifndef SINGLE_PREC
inline void Cp3dfft_ftran_r2c_many(double *A, int dim_in, double *B, int dim_out, int nv, unsigned char *op)
{
  FORT_MOD_NAME(p3dfft_ftran_r2c_many)(A,&dim_in,B,&dim_out,&nv,op);
}
#else
inline void Cp3dfft_ftran_r2c_many(float *A, int dim_in, float *B, int dim_out, int nv, unsigned char *op)
{
  FORT_MOD_NAME(p3dfft_ftran_r2c_many)(A,&dim_in,B,&dim_out,&nv,op);
}
#endif

#ifndef SINGLE_PREC
inline void Cp3dfft_btran_c2r_many(double *A, int dim_in, double *B, int dim_out, int nv, unsigned char *op)
{
  FORT_MOD_NAME(p3dfft_btran_c2r_many)(A,&dim_in,B,&dim_out,&nv,op);
}
#else
inline void Cp3dfft_btran_c2r_many(float *A, int dim_in, float *B, int dim_out, int nv, unsigned char *op)
{
  FORT_MOD_NAME(p3dfft_btran_c2r_many)(A,&dim_in,B,&dim_out,&nv,op);
}
#endif

#ifdef __cplusplus
}
#endif
