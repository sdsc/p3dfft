/*
! This sample program illustrates the 
! use of P3DFFT library for highly scalable parallel 3D FFT. 
!
! This program initializes a 3D array of complex numbers with a 
! 3D sine/cosine wave, then performs inverse FFT transform, and checks that 
! the results are correct. This sample program also demonstrates 
! how to work with complex arrays in wavenumber space, declared as real.
!
! The program expects 'stdin' file in the working directory, with 
! a single line of numbers : Nx,ny,nz,Ndim,Nrep. Here Nx,ny,nz
! are box dimensions, Ndim is the dimentionality of processor grid
! (1 or 2), and Nrep is the number of repititions. Optionally
! a file named 'dims' can also be provided to guide in the choice 
! of processor geometry in case of 2D decomposition. It should contain 
! two numbers in a line, with their product equal to the total number
! of tasks. Otherwise processor grid geometry is chosen automatically.
! For better performance, experiment with this setting, varying 
! iproc and jproc. In many cases, minimizing iproc gives best results. 
! Setting it to 1 corresponds to one-dimensional decomposition.
!
! If you have questions please contact Dmitry Pekurovsky, dmitry@sdsc.edu
*/

#include "p3dfft.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define max(a,b) ((a>b)?a:b)

double FORTNAME(t1),FORTNAME(t2),FORTNAME(t3),FORTNAME(t4),FORTNAME(tp1);
/* double t1,t2,t3,t4,tp1; */

int main(int argc,char **argv)
{
#ifndef SINGLE_PREC
   double *A,*B,*p,*C;
#else
   float *A,*B,*p,*C;
#endif
   int i,j,k,x,y,z,nx,ny,nz,proc_id,nproc,dims[2],ndim,nu;
   int istart[3],isize[3],iend[3];
   int fstart[3],fsize[3],fend[3];
   int iproc,jproc,ng[3],iex,conf,m,n;
   long int Nglob,Ntot;
   double rtime1,rtime2,gt[12],gt1[12],gt2[12],timers[12];
#ifndef SINGLE_PREC
   double cdiff,ccdiff,prec;
#else
   float cdiff,ccdiff,prec;
#endif
   FILE *fp;

#ifndef SINGLE_PREC
   void print_all(double *,int,long int),mult_array(double *,long int,double);
   void print_all_init(double *,int,long int);
   void init_wave1(double *,int *,int *,int,int,int);
   void init_wave2(double *,int *,int *,int,int,int);
#else
   void print_all(float *,int,long int),mult_array(float *,long int,double);
   void print_all_init(float *,int,long int);
   void init_wave1(float *,int *,int *,int,int,int);
   void init_wave2(float *,int *,int *,int,int,int);
#endif

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&nproc);
   MPI_Comm_rank(MPI_COMM_WORLD,&proc_id);

   for(i=0; i< 12; i++) {
     gt[i] = gt1[i] = 0.0;
     gt2[i] = 1E10;
   }
        
   set_timers();
 
   if(proc_id == 0) {
     if((fp=fopen("stdin", "r"))==NULL){
        printf("Cannot open file. Setting to default nx=ny=nz=128, ndim=2, n=1.\n");
        nx=ny=nz=128; n=1;
     } else {
        fscanf(fp,"%d %d %d %d %d\n",&nx,&ny,&nz,&ndim,&n);
        fclose(fp);
     }
#ifndef SINGLE_PREC
     printf("Double precision\n (%d %d %d) grid\n %d proc. dimensions\n%d repetitions\n",nx,ny,nz,ndim,n);
#else
     printf("Single precision\n (%d %d %d) grid\n %d proc. dimensions\n%d repetitions\n",nx,ny,nz,ndim,n);
#endif
   }
   MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ny,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nz,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ndim,1,MPI_INT,0,MPI_COMM_WORLD);
   
   if(ndim == 1) {
     dims[0] = 1; dims[1] = nproc;
   }
   else if(ndim == 2) {
     fp = fopen("dims","r");
     if(fp != NULL) {
       if(proc_id == 0)
         printf("Reading proc. grid from file dims\n");
       fscanf(fp,"%d %d\n",dims,dims+1);
       fclose(fp);
       if(dims[0]*dims[1] != nproc) 
          dims[1] = nproc / dims[0];
     }
     else {
       if(proc_id == 0) 
          printf("Creating proc. grid with mpi_dims_create\n");
       dims[0]=dims[1]=0;
       MPI_Dims_create(nproc,2,dims);
       if(dims[0] > dims[1]) {
          dims[0] = dims[1];
          dims[1] = nproc/dims[0];
       }
     }
   }

   if(proc_id == 0) 
      printf("Using processor grid %d x %d\n",dims[0],dims[1]);

   /* Initialize P3DFFT */
   p3dfft_setup(dims,nx,ny,nz,1);
   /* Get dimensions for input array - complex numbers, Z-pencil shape.
      Stride-1 dimension could be X or Z, depending on how the library 
      was compiled (stride1 option).
      Note : returns Fortran-style dimensions, i.e. starting with 1. */
   conf = 2;
   p3dfft_get_dims(istart,iend,isize,conf);
   /* Get dimensions for input array - real numbers, X-pencil shape.
      Note that we are following the Fortran ordering, i.e. 
      the dimension  with stride-1 is X. */
   conf = 1;
   p3dfft_get_dims(fstart,fend,fsize,conf);

   /* Allocate input and output arrays; note the extra factor of 2 for complex numbers*/
#ifndef SINGLE_PREC
   A = (double *) malloc(sizeof(double) * isize[0]*isize[1]*isize[2]*2);
   B = (double *) malloc(sizeof(double) * fsize[0]*fsize[1]*fsize[2]);
#else
   A = (float *) malloc(sizeof(float) * isize[0]*isize[1]*isize[2]*2);
   B = (float *) malloc(sizeof(float) * fsize[0]*fsize[1]*fsize[2]);
#endif

   if(A == NULL) 
     printf("%d: Error allocating array A (%d)\n",proc_id,isize[0]*isize[1]*isize[2]*2);

   if(B == NULL) 
     printf("%d: Error allocating array B (%d)\n",proc_id,fsize[0]*fsize[1]*fsize[2]);

#ifdef STRIDE1
   init_wave1(A,isize,istart,nx,ny,nz);
#else
   init_wave2(A,isize,istart,nx,ny,nz);
#endif

   Nglob = nx * ny;
   Nglob *= nz;

  /*
   if(proc_id == 0) printf("Initial array:\n");
   print_all_init(A,proc_id,Nglob);
   */

   rtime1 = 0.0;
   for(m=0;m < n;m++) {

     if(proc_id == 0) 
        printf("Iteration %d\n",m);

     /* Compute backward transform on A, store results in B */
     MPI_Barrier(MPI_COMM_WORLD);
     rtime1 = rtime1 - MPI_Wtime();
     p3dfft_btran_c2r(A,B);
     rtime1 = rtime1 + MPI_Wtime();

     if(proc_id == 0) 
       printf("Results of inverse transform: \n");
     print_all(B,proc_id,Nglob);

   } 
   
   /* free work space */
  p3dfft_clean();
  
  /* Check results */
  cdiff = 0.0; p = B;
  for(z=0;z < fsize[2];z++)
    for(y=0;y < fsize[1];y++)  
      for(x=0;x < fsize[0];x++) {
	 if(x +fstart[0] == nx) 
	   if(y +fstart[1] == 3 && z+fstart[2] == 4) 
	     cdiff = max(fabs((*p) + nx*ny*nz*0.25),cdiff);
	   else if(y +fstart[1]== 3 && z+fstart[2] == nz-2)
	     cdiff = max(fabs((*p) - nx*ny*nz*0.25),cdiff);
	   else if(y +fstart[1]== ny-1 && z+fstart[2] == 4)
	     cdiff = max(fabs((*p) - nx*ny*nz*0.25),cdiff);
	   else if(y +fstart[1]== ny-1 && z+fstart[2] == nz-2)
	     cdiff = max(fabs((*p) + nx*ny*nz*0.25),cdiff);
	   else
	     cdiff = max(fabs(*p),cdiff);
	 else
	     cdiff = max(fabs(*p),cdiff);
	 p++;
      }


   MPI_Reduce(&cdiff,&ccdiff,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  if(proc_id == 0) {
#ifndef SINGLE_PREC
    prec = 1.0e-14;
#else
    prec = 1.0e-5;
#endif
    if(ccdiff > prec * nx*ny*nz*0.25)
      printf("Results are incorrect\n");
    else
      printf("Results are correct\n");

    printf("max diff =%g\n",ccdiff);
  }

  get_timers(timers);

  /* Gather timing statistics */
  MPI_Reduce(&rtime1,&rtime2,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

  for (i=0;i < 12;i++) {
    timers[i] = timers[i] / ((double) n);
  }

  MPI_Reduce(&timers,&gt,12,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&timers,&gt1,12,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&timers,&gt2,12,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

  for (i=0;i < 12;i++) {
    gt[i] = gt[i]/ ((double) nproc);
  }

  if(proc_id == 0) { 
     printf("Time per loop=%lg\n",rtime2/((double) n));
     for(i=0;i < 12;i++) {
       printf("timer[%d] (avg/max/min): %lE %lE %lE\n",i+1,gt[i],gt1[i],gt2[i]);
     }
  }

  MPI_Finalize();

}

#ifndef SINGLE_PREC
void init_wave1(double *A,int isize[3],int istart[3],int nx,int ny,int nz) {
#else
  void init_wave1(float *A,int isize[3],int istart[3],int nx,int ny,int nz) {
#endif
  
  int x,y,z;
  double sinyz,*sinx,*cosx,*siny,*sinz,twopi;

   cosx = malloc(sizeof(double)*nx);
   sinx = malloc(sizeof(double)*nx);
   siny = malloc(sizeof(double)*ny);
   sinz = malloc(sizeof(double)*nz);

   twopi = atan(1.0)*8.0;

   for(z=0;z < isize[0];z++)
     sinz[z] = sin(3.0*(z+istart[0]-1)*twopi/nz);
   for(y=0;y < isize[1];y++)
     siny[y] = sin(2.0*(y+istart[1]-1)*twopi/ny);
   for(x=0;x < isize[2];x++) {
     sinx[x] = sin((x+istart[2]-1)*twopi/nx);
     cosx[x] = cos((x+istart[2]-1)*twopi/nx);
   }

   for(x=0;x < isize[2];x++)
     for(y=0;y < isize[1];y++) 
       for(z=0;z < isize[0];z++) {
	 *A++ = cosx[x]*siny[y]*sinz[z];
	 *A++ = sinx[x]*siny[y]*sinz[z];
       }

}

#ifndef SINGLE_PREC
  void init_wave2(double *A,int isize[3],int istart[3],int nx,int ny,int nz) {
#else
    void init_wave2(float *A,int isize[3],int istart[3],int nx,int ny,int nz) {
#endif
  
  int x,y,z;
  double sinyz,*sinx,*cosx,*siny,*sinz,twopi;

   cosx = malloc(sizeof(double)*nx);
   sinx = malloc(sizeof(double)*nx);
   siny = malloc(sizeof(double)*ny);
   sinz = malloc(sizeof(double)*nz);

   twopi = atan(1.0)*8.0;

   for(z=0;z < isize[2];z++)
     sinz[z] = sin(3.0*(z+istart[2]-1)*twopi/nz);
   for(y=0;y < isize[1];y++)
     siny[y] = sin(2.0*(y+istart[1]-1)*twopi/ny);
   for(x=0;x < isize[0];x++) {
     sinx[x] = sin((x+istart[0]-1)*twopi/nx);
     cosx[x] = cos((x+istart[0]-1)*twopi/nx);
   }


   for(z=0;z < isize[2];z++)
     for(y=0;y < isize[1];y++) {
       sinyz = siny[y]*sinz[z];
       for(x=0;x < isize[0];x++) {
          *A++ = cosx[x]*sinyz;
          *A++ = sinx[x]*sinyz;
       }
     }

}

#ifndef SINGLE_PREC
void mult_array(double *A,long int nar,double f)
#else
void mult_array(float *A,long int nar,double f)
#endif
{
  long int i;

  for(i=0;i < nar;i++)
    A[i] *= f;
}

#ifndef SINGLE_PREC
void print_all(double *A,int proc_id,long int Nglob)
#else
void print_all(float *A,int proc_id,long int Nglob)
#endif
{
  int x,y,z,conf,Fstart[3],Fsize[3],Fend[3];
  long int i,nar;

  conf = 1;
  p3dfft_get_dims(Fstart,Fend,Fsize,conf);
  nar = Fsize[0]*Fsize[1]*Fsize[2];

  for(i=0;i < nar;i++)
    if(fabs(A[i]) > Nglob *1.25e-4) {
      z = i/(Fsize[0]*Fsize[1]);
      y = i/(Fsize[0]) - z*Fsize[1];
      x = i-z*Fsize[0]*Fsize[1] - y*Fsize[0];
      printf("(%d,%d,%d) %lf\n",x+Fstart[0],y+Fstart[1],z+Fstart[2],A[i]);
    }
}

#ifndef SINGLE_PREC
void print_all_init(double *A,int proc_id,long int Nglob)
#else
void print_all_init(float *A,int proc_id,long int Nglob)
#endif
{
  int x,y,z,conf,Fstart[3],Fsize[3],Fend[3];
  long int i,nar;

  conf = 2;
  p3dfft_get_dims(Fstart,Fend,Fsize,conf);
  nar = Fsize[0]*Fsize[1]*Fsize[2]*2;

  for(i=0;i < nar;i+=2)
    if(fabs(A[i]) + fabs(A[i+1])> 1.25e-4) {
      z = i/(2*Fsize[0]*Fsize[1]);
      y = i/(2*Fsize[0]) - z*Fsize[1];
      x = i/2-z*Fsize[0]*Fsize[1] - y*Fsize[0];
      printf("(%d,%d,%d) %lf %lf\n",x+Fstart[0],y+Fstart[1],z+Fstart[2],A[i],A[i+1]);
    }
}

  
