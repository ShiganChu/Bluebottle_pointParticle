/*
extern "C"
{
#include "bluebottle.h"
#include "stdio.h"
#include "cuda.h"
//#include "cuda_scalar.h"
}
*/
#include "cuda_scalar.h"
#include "cuda_point.h"

__global__ void array_init(real *A,dom_struct *dom, int n, real a)
{
int indx =  threadIdx.x + blockIdx.x*blockDim.x;
int indy =  threadIdx.y + blockIdx.y*blockDim.y;
int pp=indx+indy*gridDim.x*blockDim.x;
if(pp<n) A[pp]=a;
}



/*
A is fluid property located at cell center or face center depending on dirc 
dirc=0 coresponds to cell-center; dirc=1 x-face-center;dirc=2 y-face-center;dirc=3 z-face-center
dir
*/
__global__ void lpt_mollify_sc(int npoints, real *A, point_struct *points, dom_struct *dom,real *Ksi, int dirc, int dirc2)
{

int indx =  threadIdx.x + blockIdx.x*blockDim.x;
int indy =  threadIdx.y + blockIdx.y*blockDim.y;
int pp=indx+indy*gridDim.x*blockDim.x;

//if(pp<npoints*STENCIL3)  Ksi[pp]=0.0;
//__sycthreads();

if(pp<npoints)
{

////printf("\npp indx %d %d %d \n",pp, indx, indy);
////printf("\nblockDim.x,gridDim.x %d %d %d %d\n",blockDim.x,gridDim.x,blockDim.y,gridDim.y );

//Define the STENCIL of the Gausian filter, the spread length of Gaussian kernel
real ksi[STENCIL][STENCIL][STENCIL];
real xs=dom->xs;real ys=dom->ys;real zs=dom->zs;
real dx=dom->dx;real dy=dom->dy;real dz=dom->dz;
  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;
real cellVol=dx*dy*dz;

real  xp =  points[pp].x;
real  yp =  points[pp].y;
real  zp =  points[pp].z;

// i,j,k should be the start number of cell center, so that xm(i)<=xp<xm(i+1)
//xm(i)= (i-0.5)*dom->dx +dom->xs
//i=floor((xm - dom->xs) * ddx+ 0.5);
//i=floor((xm - dom->xs) * ddx- 0.5) + DOM_BUF;
  int ip,jp,kp;
  int is,js,ks;

//The cell-center or face-center index to locate the particle 
switch(dirc)
{
case 0:
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for xm[i]
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;//for ym[j]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;//for zm[k]
  break;
case 1:
  ip = floor((xp - dom->xs) * ddx + 0.5) + DOM_BUF;//for xm[i+1], ip should be i+1
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;	   //for y[i]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;	   //for z[i]
  break;
case 2:
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for x[i] 
  jp = floor((yp - dom->ys) * ddy + 0.5) + DOM_BUF;	   //for ym[i+1]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;	   //for z[i]
  break;
case 3:
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for x[i] 
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;	   //for y[i]
  kp = floor((zp - dom->zs) * ddz+ 0.5) + DOM_BUF;  //for zm[i+1]
  break;
default:
  printf("Wrong direction in lpt_mollify");
}

points[pp].i=ip;
points[pp].j=jp;
points[pp].k=kp;

real buf=0;
int began=- round((STENCIL-1)/2.0f);//equal to -1 if STENCIL=3
int end=1+ ceil((STENCIL-1)/2.0f);//equal to 2 if STENCIL=3


//Calculate the filter strength from Gaussian kernel
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
{
is=di-began;js=dj-began;ks=dk-began;
ksi[is][js][ks]=lpt_integrate_mol(ip+di,jp+dj,kp+dk,xp,yp,zp,dx,dy,dz,xs,ys,zs,dirc);
buf+=ksi[is][js][ks];
}

//TODO add mask infomation as in lpt_mollify_sc in lpt_interpolator.f90

real Ap;
switch(dirc)
{
case 0:
    switch(dirc2)
   { case 0:   
//Exchange rate of soluble mass from particles. All the following source should divide cellVol to ensure conservation law
     Ap=-points[pp].msdot/cellVol;break;
//     Ap=1/cellVol;break;

     case 1:
//Particle volume filter
     Ap=PI*4/3*points[pp].r*points[pp].r*points[pp].r/cellVol;break;
     default: Ap=0;break;
	}
break;
case 1:
//Particle x-momentum reaction to fluid
Ap=-points[pp].Fx/cellVol;
break;
case 2:
Ap=-points[pp].Fy/cellVol;
break;
case 3:
//Particle z-momentum reaction to fluid
Ap=-points[pp].Fz/cellVol;
break;
default:
break;
	}		


// Normalize  ksi = ksi/buf
if(buf>0)
	{
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
							{
is=di-began;js=dj-began;ks=dk-began;
ksi[is][js][ks]=ksi[is][js][ks]/buf;
Ksi[is+js*STENCIL+ks*STENCIL2+pp*STENCIL3]=ksi[is][js][ks]*Ap; }
		}


	}

}


//Add the point source to fluid source weighted by the Gaussian filter coefficient Ksi
__global__ void lpt_sum_ksi(int npoints, real *A, point_struct *points, dom_struct *dom,real *Ksi, int dirc,int dirc2)
{
int indx =  threadIdx.x + blockIdx.x*blockDim.x;
int indy =  threadIdx.y + blockIdx.y*blockDim.y;
int np=indx+indy*gridDim.x*blockDim.x;
if(np<1)
{

int began=-round((STENCIL-1)/2.0);//equal to -1 if STENCIL=3
int end=1+ ceil((STENCIL-1)/2.0);//equal to 2 if STENCIL=3

//printf("\nST %d %d\n",end-began+1,STENCIL);
  int ip,jp,kp;
  int is,js,ks;
  int ic,jc,kc;
real ksi;

for(int pp=0;pp<npoints;pp++)
{
ip =  points[pp].i;
jp =  points[pp].j;
kp =  points[pp].k;

// Perform the actual extrapolation on A
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
{
  ic=di+ip;jc=dj+jp;kc=dk+kp;
  is=di-began;js=dj-began;ks=dk-began;
  ksi=Ksi[is+js*STENCIL+ks*STENCIL2+pp*STENCIL3];
switch(dirc)
{
case 0:
//Make the index inside the domain before we change the property at this location
periodic_grid_index(ic,jc,kc,dom,0);
A[ic+jc*dom->Gcc._s1b+kc*dom->Gcc._s2b]+=ksi;
break;
case 1:
periodic_grid_index(ic,jc,kc,dom,1);
A[ic+jc*dom->Gfx._s1b+kc*dom->Gfx._s2b]+=ksi;
break;
case 2:
periodic_grid_index(ic,jc,kc,dom,2);
A[ic+jc*dom->Gfy._s1b+kc*dom->Gfy._s2b]+=ksi;
break;
case 3:
periodic_grid_index(ic,jc,kc,dom,3);
A[ic+jc*dom->Gfz._s1b+kc*dom->Gfz._s2b]+=ksi;
break;
default:
break;
			}		

		}
	}
  }
}



/*
__global__ void lpt_mollify_sc(int npoints, real *A, point_struct *points, dom_struct *dom,real *Ksi, int dirc, int dirc2)
{

int indx =  threadIdx.x + blockIdx.x*blockDim.x;
int indy =  threadIdx.y + blockIdx.y*blockDim.y;
int pp=indx+indy*gridDim.x*blockDim.x;

//if(pp<npoints*STENCIL3)  Ksi[pp]=0.0;
//__sycthreads();

if(pp<npoints)
{

////printf("\npp indx %d %d %d \n",pp, indx, indy);
////printf("\nblockDim.x,gridDim.x %d %d %d %d\n",blockDim.x,gridDim.x,blockDim.y,gridDim.y );

//Define the STENCIL of the Gausian filter, the spread length of Gaussian kernel
real ksi[STENCIL][STENCIL][STENCIL];
real xs=dom->xs;real ys=dom->ys;real zs=dom->zs;
real dx=dom->dx;real dy=dom->dy;real dz=dom->dz;
  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;
//real cellVol=dx*dy*dz;

real  xp =  points[pp].x;
real  yp =  points[pp].y;
real  zp =  points[pp].z;

// i,j,k should be the start number of cell center, so that xm(i)<=xp<xm(i+1)
//xm(i)= (i-0.5)*dom->dx +dom->xs
//i=floor((xm - dom->xs) * ddx+ 0.5);
//i=floor((xm - dom->xs) * ddx- 0.5) + DOM_BUF;
  int ip,jp,kp;
  int is,js,ks;

//The cell-center or face-center index to locate the particle 
switch(dirc)
{
case 0:
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for xm[i]
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;//for ym[j]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;//for zm[k]
  break;
case 1:
  ip = floor((xp - dom->xs) * ddx + 0.5) + DOM_BUF;//for xm[i+1], ip should be i+1
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;	   //for y[i]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;	   //for z[i]
  break;
case 2:
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for x[i] 
  jp = floor((yp - dom->ys) * ddy + 0.5) + DOM_BUF;	   //for ym[i+1]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;	   //for z[i]
  break;
case 3:
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for x[i] 
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;	   //for y[i]
  kp = floor((zp - dom->zs) * ddz+ 0.5) + DOM_BUF;  //for zm[i+1]
  break;
default:
  printf("Wrong direction in lpt_mollify");
}

points[pp].i=ip;
points[pp].j=jp;
points[pp].k=kp;

real buf=0;
int began=-(STENCIL-1)/2;//equal to -1 if STENCIL=3
int end=1+(STENCIL-1)/2;//equal to 2 if STENCIL=3


//Calculate the filter strength from Gaussian kernel
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
{
is=di-began;js=dj-began;ks=dk-began;
ksi[is][js][ks]=lpt_integrate_mol(ip+di,jp+dj,kp+dk,xp,yp,zp,dx,dy,dz,xs,ys,zs,dirc);
buf+=ksi[is][js][ks];
}

//TODO add mask infomation as in lpt_mollify_sc in lpt_interpolator.f90

// Normalize  ksi = ksi/buf
if(buf>0)
{
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
	{
is=di-began;js=dj-began;ks=dk-began;
ksi[is][js][ks]=ksi[is][js][ks]/buf;
Ksi[is+js*STENCIL+ks*STENCIL2+pp*STENCIL3]=ksi[is][js][ks];
			}
		}


	}

}


//Add the point source to fluid source weighted by the Gaussian filter coefficient Ksi
__global__ void lpt_sum_ksi(int npoints, real *A, point_struct *points, dom_struct *dom,real *Ksi, int dirc,int dirc2)
{
int indx =  threadIdx.x + blockIdx.x*blockDim.x;
int indy =  threadIdx.y + blockIdx.y*blockDim.y;
int np=indx+indy;
if(np<1)
{

int began=-(STENCIL-1)/2;//equal to -1 if STENCIL=3
int end=1+(STENCIL-1)/2;//equal to 2 if STENCIL=3
  int ip,jp,kp;
  int is,js,ks;
  int ic,jc,kc;
real Ap,ksi;
real cellVol=dom->dx*dom->dy*dom->dz;

for(int pp=0;pp<npoints;pp++)
{

ip =  points[pp].i;
jp =  points[pp].j;
kp =  points[pp].k;

// Perform the actual extrapolation on A
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
{
  ic=di+ip;jc=dj+jp;kc=dk+kp;
  is=di-began;js=dj-began;ks=dk-began;
  ksi=Ksi[is+js*STENCIL+ks*STENCIL2+pp*STENCIL3];
switch(dirc)
{
case 0:
//Make the index inside the domain before we change the property at this location
periodic_grid_index(ic,jc,kc,dom,0);
    switch(dirc2)
   { case 0:   
//Exchange rate of soluble mass from particles. All the following source should divide cellVol to ensure conservation law
     Ap=-points[pp].msdot/cellVol;break;
//    Ap=1/cellVol;break;
     case 1:
//Particle volume filter
     Ap=PI*4/3*points[pp].r*points[pp].r*points[pp].r/cellVol;break;
     default:
printf("Wrong direction in cell center mollify");
	}
 A[ic+jc*dom->Gcc._s1b+kc*dom->Gcc._s2b]+=ksi*Ap;
break;
case 1:
periodic_grid_index(ic,jc,kc,dom,1);

//Particle x-momentum reaction to fluid
Ap=-points[pp].Fx/cellVol;
//Ap=1/cellVol;
//Ap=0;
A[ic+jc*dom->Gfx._s1b+kc*dom->Gfx._s2b]+=ksi*Ap;
break;
case 2:
periodic_grid_index(ic,jc,kc,dom,2);
//Particle y-momentum reaction to fluid
Ap=-points[pp].Fy/cellVol;
//Ap=1/cellVol;
//Ap=0;
A[ic+jc*dom->Gfy._s1b+kc*dom->Gfy._s2b]+=ksi*Ap;
break;
case 3:
periodic_grid_index(ic,jc,kc,dom,3);
//Particle z-momentum reaction to fluid
Ap=-points[pp].Fz/cellVol;
//Ap=1/cellVol;
//Ap=0;
A[ic+jc*dom->Gfz._s1b+kc*dom->Gfz._s2b]+=ksi*Ap;

break;
default:
break;
			}		

		}
	}
  }
}

*/


__device__ void periodic_grid_index(int ic,int jc,int kc,dom_struct *dom, int dirc)
{
//int tag=0;
switch(dirc)
{
case 0:
if(ic<dom->Gcc.is||ic > dom->Gcc.ie-1)  ic=(int) fmodf(ic,dom->Gcc.in);
if(jc<dom->Gcc.js||jc > dom->Gcc.je-1)  jc=(int) fmodf(jc,dom->Gcc.jn);
if(kc<dom->Gcc.ks||kc > dom->Gcc.ke-1)  kc=(int) fmodf(kc,dom->Gcc.kn);
//after fmod, the value could still be negative
  if(ic < dom->Gcc.is) ic = ic +dom->Gcc.in;
  if(jc < dom->Gcc.js) jc = jc +dom->Gcc.jn;
  if(kc < dom->Gcc.ks) kc = kc +dom->Gcc.kn;
 break;
case 1:
if(ic<dom->Gfx.is||ic > dom->Gfx.ie-2)  ic=(int) fmodf(ic,dom->Gfx.in-1);
if(jc<dom->Gfx.js||jc > dom->Gfx.je-1)  jc=(int) fmodf(jc,dom->Gfx.jn);
if(kc<dom->Gfx.ks||kc > dom->Gfx.ke-1)  kc=(int) fmodf(kc,dom->Gfx.kn);
  if(ic < dom->Gfx.is) ic = ic+(dom->Gfx.in-1);
  if(jc < dom->Gfx.js) jc = jc+dom->Gfx.jn;
  if(kc < dom->Gfx.ks) kc = kc+dom->Gfx.kn;
break;
case 2:
if(ic<dom->Gfy.is||ic > dom->Gfy.ie-1)  ic=(int) fmodf(ic,dom->Gfy.in);
if(jc<dom->Gfy.js||jc > dom->Gfy.je-2)  jc=(int) fmodf(jc,dom->Gfy.jn-1);
if(kc<dom->Gfy.ks||kc > dom->Gfy.ke-1)  kc=(int) fmodf(kc,dom->Gfy.kn);
  if(ic < dom->Gfy.is) ic = ic+dom->Gfy.in;
  if(jc < dom->Gfy.js) jc = jc+(dom->Gfy.jn-1);
  if(kc < dom->Gfy.ks) kc = kc+dom->Gfy.kn;
break;
case 3:
if(ic<dom->Gfz.is||ic > dom->Gfz.ie-1)  ic=(int) fmodf(ic,dom->Gfz.in);
if(jc<dom->Gfz.js||jc > dom->Gfz.je-1)  jc=(int) fmodf(jc,dom->Gfz.jn);
if(kc<dom->Gfz.ks||kc > dom->Gfz.ke-2)  kc=(int) fmodf(kc,dom->Gfz.kn-1);
  if(ic < dom->Gfz.is) ic = ic+dom->Gfz.in;
  if(jc < dom->Gfz.js) jc = jc+dom->Gfz.jn;
  if(kc < dom->Gfz.ks) kc = kc+(dom->Gfz.kn-1);
break;
default:
break;
}

}

__device__ void periodic_grid_position(real x,real y,real z,dom_struct *dom)
{
if(x<dom->xs||x>dom->xe)  x=fmodf(x,dom->xl);
if(y<dom->ys||y>dom->ye)  y=fmodf(y,dom->yl);
if(z<dom->zs||z>dom->ze)  z=fmodf(z,dom->zl);

  if(x < dom->xs ) x = x + dom->xl;
  if(y < dom->ys ) y = y + dom->yl;
  if(z < dom->zs ) z = z + dom->zl;
}

//Gausian kernel to calculate weight coefficient of the filter
__device__ real lpt_integrate_mol(int ic,int jc,int kc,real xp,real yp,real zp, real dx,real dy,real dz,real xs,real ys,real zs, int dirc)
{

real xm,ym,zm;
real x,y,z,r2;
switch(dirc)
{
case 0:
xm = (ic-DOM_BUF+0.5) * dx + xs;
ym = (jc-DOM_BUF+0.5) * dy + ys;
zm = (kc-DOM_BUF+0.5) * dz + zs;
r2 = (xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm);
break;
case 1:
x =  (ic-DOM_BUF) * dx + xs;
ym = (jc-DOM_BUF+0.5) * dy + ys;
zm = (kc-DOM_BUF+0.5) * dz + zs;
r2 = (xp-x)*(xp-x)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm);
break;
case 2:
xm = (ic-DOM_BUF+0.5) * dx + xs;
y = (jc-DOM_BUF) * dy + ys;
zm = (kc-DOM_BUF+0.5) * dz + zs;
r2 = (xp-xm)*(xp-xm)+(yp-y)*(yp-y)+(zp-zm)*(zp-zm);
break;
case 3:
xm = (ic-DOM_BUF+0.5) * dx + xs;
ym = (jc-DOM_BUF+0.5) * dy + ys;
z = (kc-DOM_BUF) * dz + zs;
r2 = (xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-z)*(zp-z);
break;
default:
printf("Wrong dirc in lpt_integrate_mol");
}

//TODO make this as defined value avaible from host and device!!!
real min_meshsize=min(min(dx,dy),dz);
//2.0f*sqrt(2.0f*log(2.0f)) = 2.355
real sig= KERNEL_WIDTH *min_meshsize/(2.0f*sqrt(2.0f*log(2.0f)));

//r should be 3 times larger than sig to make val smaller than 0.01!!
real val = exp(-r2/(2.0f*sig*sig));
return val;
}




// xp~zp is the particle location, x(1~N) is the grid face position
// if x(ip)<=xp<x(ip+1) points[pp].i=ip, this is the grid face number rather than grid center number
__global__ void lpt_localize(int npoints, point_struct *points, dom_struct *dom, BC bc)
{
 int pp =  threadIdx.x + blockIdx.x*blockDim.x;

if(pp<npoints)
{
  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;

// Cartesian location of node
real  xp =  points[pp].x;
real  yp =  points[pp].y;
real  zp =  points[pp].z;

//TODO whether periodic BC for point particle need to be determined in future
periodic_grid_position(xp,yp,zp,dom);
 
  int ip,jp,kp;
//x(i)=(i-DOM_BUF)*dom->dx+dom->xs
//i=floor((x(i) - dom->xs) * ddx) + DOM_BUF;
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for x[i]
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;//for y[j]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;//for z[k]

points[pp].i= ip;
points[pp].j= jp;
points[pp].k= kp;
}

}


__global__ void stress_u(real rho_f, real nu, real *u0,real *p,real *p0,real *stress, dom_struct *dom,int *flag_u, real dt, real dt0)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_u0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u back
  __shared__ real s_u1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u center
  __shared__ real s_u2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u forward

  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  
  __shared__ real grad_P[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real grad_P0[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
 
// working constants

/*
  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??
  real ab0= 0.5 * a ;
  real ab=1. + ab0 ;
  */

//Weight for Adam-Bashforth interpolation
  real a = dt0/dt;
  a = (a + 2.)/(a + 1.);
  real ab0= a-1.0 ;
  real ab=a;

  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed


  // loop over u-planes
  for(int i = dom->Gfx._is; i < dom->Gfx._ie; i++) {

 // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int j = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int k = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int tj = threadIdx.x;
    int tk = threadIdx.y;

 // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations
    // TODO: THIS CAN BE FIXED BY PADDING ALL OF THESE ARRAYS WHEN COPYING FROM
    // HOST TO DEVICE
    if((k >= dom->Gfx._ksb && k < dom->Gfx._keb)
      && (j >= dom->Gfx._jsb && j < dom->Gfx._jeb)) {
      s_u0[tj + tk*blockDim.x] = u0[(i-1) + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u1[tj + tk*blockDim.x] = u0[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u2[tj + tk*blockDim.x] = u0[(i+1) + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
   
    __syncthreads();


    // if off the shared memory block boundary
    if((tj > 0 && tj < blockDim.x-1) && (tk > 0 && tk < blockDim.y-1)) {
      real u011 = s_u0[tj + tk*blockDim.x];
      real u111 = s_u1[tj + tk*blockDim.x];
      real u211 = s_u2[tj + tk*blockDim.x];

      real u101 = s_u1[(tj-1) + tk*blockDim.x];
      real u121 = s_u1[(tj+1) + tk*blockDim.x];

      real u110 = s_u1[tj + (tk-1)*blockDim.x];
      real u112 = s_u1[tj + (tk+1)*blockDim.x];

 // compute diffusion term (Adams-Bashforth stepping)
      real dud1 = (u211 - u111) * ddx;
      real dud0 = (u111 - u011) * ddx;
      real ddudxx = (dud1 - dud0) * ddx;

      dud1 = (u121 - u111) * ddy;
      dud0 = (u111 - u101) * ddy;
      real ddudyy = (dud1 - dud0) * ddy;

      dud1 = (u112 - u111) * ddz;
      dud0 = (u111 - u110) * ddz;
      real ddudzz = (dud1 - dud0) * ddz;

      s_d[tj + tk*blockDim.x] = nu * (ddudxx + ddudyy + ddudzz);
/*
      grad_P[tj + tk*blockDim.x]=abs(flag_u[i + tj*dom->Gfx._s1b
        + tk*dom->Gfx._s2b])
        * ddx * (p[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p[(i-1) + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]);
  */    
      grad_P0[tj + tk*blockDim.x]=abs(flag_u[i + tj*dom->Gfx._s1b
        + tk*dom->Gfx._s2b])
        * ddx * (p0[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p0[(i-1) + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]);

    }
   // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((k >= dom->Gfx._ks && k < dom->Gfx._ke)
      && (j >= dom->Gfx._js && j < dom->Gfx._je)
      && (tj > 0 && tj < (blockDim.x-1))
      && (tk > 0 && tk < (blockDim.y-1))) { 
//     stress[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] =rho_f* s_d[tj + tk*blockDim.x]-ab*grad_P[tj + tk*blockDim.x]+ab0*grad_P0[tj + tk*blockDim.x]-gradP_x;
     stress[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] =rho_f* s_d[tj + tk*blockDim.x]-ab*grad_P[tj + tk*blockDim.x]+ab0*grad_P0[tj + tk*blockDim.x];
// stress[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] =rho_f* s_d[tj + tk*blockDim.x]-grad_P[tj + tk*blockDim.x];
    }
  }
} 


__global__ void stress_v(real rho_f, real nu, real *v0,real *p,real *p0,real *stress, dom_struct *dom,int *flag_v, real dt, real dt0)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_v0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v back
  __shared__ real s_v1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v center
  __shared__ real s_v2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v forward

  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff

  __shared__ real grad_P[MAX_THREADS_DIM * MAX_THREADS_DIM];       // y-force
 __shared__ real grad_P0[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
 
// working constants
/*  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??
  real ab0= 0.5 * a ;
  real ab=1. + ab0 ;
  */
 real a = dt0/dt;
  a = (a + 2.)/(a + 1.);
  real ab0= a-1.0 ;
  real ab=a;

  // working constants
 
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

  // loop over u-planes
  for(int j = dom->Gfy._js; j < dom->Gfy._je; j++) {
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int k = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int i = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int tk = threadIdx.x;
    int ti = threadIdx.y;

    // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations
    if((i >= dom->Gfy._isb && i < dom->Gfy._ieb)
      && (k >= dom->Gfy._ksb && k < dom->Gfy._keb)) {
     
      s_v0[tk + ti*blockDim.x] = v0[i + (j-1)*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v1[tk + ti*blockDim.x] = v0[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v2[tk + ti*blockDim.x] = v0[i + (j+1)*dom->Gfy._s1b + k*dom->Gfy._s2b];
     }
   
    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute right-hand side
    // if off the shared memory block boundary
    if((tk > 0 && tk < blockDim.x-1) && (ti > 0 && ti < blockDim.y-1)) {
      // grab the required data points for calculations
      real v101 = s_v0[tk + ti*blockDim.x];
      real v111 = s_v1[tk + ti*blockDim.x];
      real v121 = s_v2[tk + ti*blockDim.x];

      real v110 = s_v1[(tk-1) + ti*blockDim.x];
      real v112 = s_v1[(tk+1) + ti*blockDim.x];
   
      real v011 = s_v1[tk + (ti-1)*blockDim.x];
      real v211 = s_v1[tk + (ti+1)*blockDim.x];

      // compute diffusive term
      real dvd1 = (v211 - v111) * ddx;
      real dvd0 = (v111 - v011) * ddx;
      real ddvdxx = (dvd1 - dvd0) * ddx;

      dvd1 = (v121 - v111) * ddy;
      dvd0 = (v111 - v101) * ddy;
      real ddvdyy = (dvd1 - dvd0) * ddy;

      dvd1 = (v112 - v111) * ddz;
      dvd0 = (v111 - v110) * ddz;
      real ddvdzz = (dvd1 - dvd0) * ddz;

      s_d[tk + ti*blockDim.x] = nu * (ddvdxx + ddvdyy + ddvdzz);
/* 
     grad_P[tk + ti*blockDim.x]=abs(flag_v[ti + j*dom->Gfy._s1b
        + tk*dom->Gfy._s2b])
        * ddy * (p[ti + j*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p[ti + (j-1)*dom->Gcc._s1b + tk*dom->Gcc._s2b]);
*/
      grad_P0[tk + ti*blockDim.x]=abs(flag_v[ti + j*dom->Gfy._s1b
        + tk*dom->Gfy._s2b])
        * ddy * (p0[ti + j*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p0[ti + (j-1)*dom->Gcc._s1b + tk*dom->Gcc._s2b]);

     
	}

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((i >= dom->Gfy._is && i < dom->Gfy._ie)
      && (k >= dom->Gfy._ks && k < dom->Gfy._ke)
      && (tk > 0 && tk < (blockDim.x-1))
      && (ti > 0 && ti < (blockDim.y-1))) {
  //    stress[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b] = rho_f*s_d[tk + ti*blockDim.x]-ab*grad_P[tk + ti*blockDim.x]+ab0*grad_P0[tk + ti*blockDim.x]-gradP_y;
//     stress[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b] = rho_f*s_d[tk + ti*blockDim.x]-grad_P[tk + ti*blockDim.x];
    stress[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b] = rho_f*s_d[tk + ti*blockDim.x]-ab*grad_P[tk + ti*blockDim.x]+ab0*grad_P0[tk + ti*blockDim.x];
      }
  }
}

__global__ void stress_w(real rho_f, real nu, real *w0,real *p,real *p0,real *stress, dom_struct *dom,int *flag_w, real dt, real dt0)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_w2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward
  
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff0
  
  __shared__ real grad_P[MAX_THREADS_DIM * MAX_THREADS_DIM];  // pressure gradient
  __shared__ real grad_P0[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
 
// working constants
/*  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??
  real ab0= 0.5 * a ;
  real ab=1. + ab0 ;
  */
  real a = dt0/dt;
  a = (a + 2.)/(a + 1.);
  real ab0= a-1.0 ;
  real ab=a;

  // working constants
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

  // loop over w-planes
  for(int k = dom->Gfz._ks; k < dom->Gfz._ke; k++) {
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int i = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int ti = threadIdx.x;
    int tj = threadIdx.y;

    // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations
    if((j >= dom->Gfz._jsb && j < dom->Gfz._jeb)
      && (i >= dom->Gfz._isb && i < dom->Gfz._ieb)) {
      
      s_w0[ti + tj*blockDim.x] = w0[i + j*dom->Gfz._s1b + (k-1)*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w0[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w2[ti + tj*blockDim.x] = w0[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
      
    }
     // make sure all threads complete shared memory copy
    __syncthreads();

    // compute right-hand side
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) {
      // grab the required data points for calculations
      real w110 = s_w0[ti + tj*blockDim.x];
      real w111 = s_w1[ti + tj*blockDim.x];
      real w112 = s_w2[ti + tj*blockDim.x];

      real w011 = s_w1[(ti-1) + tj*blockDim.x];
      real w211 = s_w1[(ti+1) + tj*blockDim.x];
   
      real w101 = s_w1[ti + (tj-1)*blockDim.x];
      real w121 = s_w1[ti + (tj+1)*blockDim.x];

      // compute diffusive term
      real dwd1 = (w211 - w111) * ddx;
      real dwd0 = (w111 - w011) * ddx;
      real ddwdxx = (dwd1 - dwd0) * ddx;

      dwd1 = (w121 - w111) * ddy;
      dwd0 = (w111 - w101) * ddy;
      real ddwdyy = (dwd1 - dwd0) * ddy;

      dwd1 = (w112 - w111) * ddz;
      dwd0 = (w111 - w110) * ddz;
      real ddwdzz = (dwd1 - dwd0) * ddz;

      s_d[ti + tj*blockDim.x] = nu * (ddwdxx + ddwdyy + ddwdzz);
/*
      grad_P[ti + tj*blockDim.x] = abs(flag_w[ti + tj*dom->Gfz._s1b
        + k*dom->Gfz._s2b])
        * ddz * (p[ti + tj*dom->Gcc._s1b + k*dom->Gcc._s2b]
        - p[ti + tj*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b]);
*/	
      grad_P0[ti + tj*blockDim.x] = abs(flag_w[ti + tj*dom->Gfz._s1b
        + k*dom->Gfz._s2b])
        * ddz * (p0[ti + tj*dom->Gcc._s1b + k*dom->Gcc._s2b]
        - p0[ti + tj*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b]);

 
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((j >= dom->Gfz._js && j < dom->Gfz._je)
      && (i >= dom->Gfz._is && i < dom->Gfz._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
//      stress[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b]= rho_f* s_d[ti + tj*blockDim.x]-ab*grad_P[ti + tj*blockDim.x]+ab0*grad_P0[ti + tj*blockDim.x]-gradP_z;
//stress[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b]= rho_f* s_d[ti + tj*blockDim.x]-grad_P[ti + tj*blockDim.x];
      
stress[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b]= rho_f* s_d[ti + tj*blockDim.x]-ab*grad_P[ti + tj*blockDim.x]+ab0*grad_P0[ti + tj*blockDim.x];
    }
  }
}



//Interpolate cell-center flow field data A to the particle position, and get Ag at each particle position
__global__ void interpolate_point_scalar_Lag2(int npoints,real *A,real *Ag, point_struct *points, dom_struct *dom)
{
  int pp =  threadIdx.x + blockIdx.x*blockDim.x;
  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;
 
if(pp<npoints)
{

real  x =  points[pp].x;
real  y =  points[pp].y;
real  z =  points[pp].z;
//real Ap=points[pp].mdot;
 
//TODO threat BC in the future
//periodic_grid_position(x,y,z,dom);

 __syncthreads();
int i,j,k;   
real x1,x2,y1,y2,z1,z2;
int ic,jc,kc;
real wx[2],wy[2],wz[2];
real weight[2][2][2];

  i = floor((x - dom->xs) * ddx- 0.5) + DOM_BUF;//for xm[i]
  j = floor((y - dom->ys) * ddy- 0.5) + DOM_BUF;//for ym[j]
  k = floor((z - dom->zs) * ddz- 0.5) + DOM_BUF;//for zm[k]

//periodic_grid_index(i,j,k,dom,0);

  x1 = (i-0.5) * dom->dx + dom->xs;
  x2 = (i+0.5) * dom->dx + dom->xs;
  y1=  (j-0.5) * dom->dy + dom->ys;
  y2=  (j+0.5) * dom->dy + dom->ys;
  z1=  (k-0.5) * dom->dz + dom->zs;
  z2=  (k+0.5) * dom->dz + dom->zs;

//2nd order lagragian interpolation
 wx[0]=(x2-x)/(x2-x1);
 wx[1]=(x-x1)/(x2-x1);
 wy[0]=(y2-y)/(y2-y1);
 wy[1]=(y-y1)/(y2-y1);
 wz[0]=(z2-z)/(z2-z1);
 wz[1]=(z-z1)/(z2-z1);

real buf=0.f;
//Normalize the weight
 for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
      { 
	weight[ii][jj][kk]=wx[ii]*wy[jj]*wz[kk];	
	buf+=weight[ii][jj][kk];
	}

//Normalize
if(fabs(buf-1)>EPSILON)
{
 for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
	weight[ii][jj][kk]=weight[ii][jj][kk]/buf;	
}

Ag[pp]=0;
//Add them to the scalar field source
 for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
      { 
  ic=i+ii;jc=j+jj;kc=k+kk;
periodic_grid_index(ic,jc,kc,dom,0);
  Ag[pp]+=weight[ii][jj][kk]*A[ic +jc*dom->Gcc.s1b + kc*dom->Gcc.s2b];
 }

//TODO add mask infomation as in lpt_mollify_sc in lpt_interpolator.f90
  }
}

//Bilinear interpolation!  This skeme is accurate enough, which has been tested by oscllating flows.
//TODO treat boudary specification,such as mask information
__global__ void interpolate_point_vel_Lag2(real *u, real *v, real *w,int npoints,real rho_f, real nu,real *ug,real *vg,real *wg, point_struct *points, dom_struct *dom, BC bc)
{
  // int node = threadIdx.x;
  int pp =  threadIdx.x + blockIdx.x*blockDim.x;
  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;


if(pp<npoints)
{
// Cartesian location of node
real  x =  points[pp].x;
real  y =  points[pp].y;
real  z =  points[pp].z;

// periodic_grid_position(x,y,z,dom);

  __syncthreads();

int i,j,k;   
real x1,x2,y1,y2,z1,z2;
int ic,jc,kc;
real wx[2],wy[2],wz[2];
real weight[2][2][2];

  // interpolate velocities
  // interpolate u-velocity
  i = floor((x - dom->xs) * ddx) + DOM_BUF; 	//for x[i]
  j = floor((y - dom->ys) * ddy- 0.5) + DOM_BUF;//for ym[j]
  k = floor((z - dom->zs) * ddz- 0.5) + DOM_BUF;//for zm[k]
//periodic_grid_index(ic,jc,kc,dom,1);

  x1 = (i-DOM_BUF) * dom->dx + dom->xs;
  x2 = (i+1-DOM_BUF) * dom->dx + dom->xs;
  y1=  (j-0.5) * dom->dy + dom->ys;
  y2=  (j+0.5) * dom->dy + dom->ys;
  z1=  (k-0.5) * dom->dz + dom->zs;
  z2=  (k+0.5) * dom->dz + dom->zs;

 
 wx[0]=(x2-x)/(x2-x1);
 wx[1]=(x-x1)/(x2-x1);
 wy[0]=(y2-y)/(y2-y1);
 wy[1]=(y-y1)/(y2-y1);
 wz[0]=(z2-z)/(z2-z1);
 wz[1]=(z-z1)/(z2-z1);

ug[pp]=0;
  for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
      { 
	
weight[ii][jj][kk]=wx[ii]*wy[jj]*wz[kk];	
ic=i+ii;
jc=j+jj;
kc=k+kk;
periodic_grid_index(ic,jc,kc,dom,1);
	ug[pp]+=weight[ii][jj][kk]*u[ic +jc*dom->Gfx.s1b + kc*dom->Gfx.s2b];
	}

  // interpolate V-velocity
  i = floor((x - dom->xs) * ddx- 0.5) + DOM_BUF;
  j = floor((y - dom->ys) * ddy) + DOM_BUF;
  k = floor((z - dom->zs) * ddz- 0.5) + DOM_BUF;

//periodic_grid_index(ic,jc,kc,dom,2);

  x1 = (i-0.5) * dom->dx + dom->xs;
  x2 = (i+0.5) * dom->dx + dom->xs;
  y1=  (j-DOM_BUF) * dom->dy + dom->ys;
  y2=  (j+1-DOM_BUF) * dom->dy + dom->ys;
  z1=  (k-0.5) * dom->dz + dom->zs;
  z2=  (k+0.5) * dom->dz + dom->zs;


 wx[0]=(x2-x)/(x2-x1);
 wx[1]=(x-x1)/(x2-x1);
 wy[0]=(y2-y)/(y2-y1);
 wy[1]=(y-y1)/(y2-y1);
 wz[0]=(z2-z)/(z2-z1);
 wz[1]=(z-z1)/(z2-z1);

vg[pp]=0;
  for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
      { 
	weight[ii][jj][kk]=wx[ii]*wy[jj]*wz[kk];
ic=i+ii;
jc=j+jj;
kc=k+kk;
periodic_grid_index(ic,jc,kc,dom,2);
	vg[pp]+=weight[ii][jj][kk]*v[ic +jc*dom->Gfy.s1b + kc*dom->Gfy.s2b];
	}


  // interpolate W-velocity
  i = floor((x - dom->xs) * ddx- 0.5) + DOM_BUF;
  j = floor((y - dom->ys) * ddy- 0.5) + DOM_BUF;
  k = floor((z - dom->zs) * ddz) + DOM_BUF;
//periodic_grid_index(ic,jc,kc,dom,3);

  x1 = (i-0.5) * dom->dx + dom->xs;
  x2 = (i+0.5) * dom->dx + dom->xs;
  y1=  (j-0.5) * dom->dy + dom->ys;
  y2=  (j+0.5) * dom->dy + dom->ys;
  z1=  (k-DOM_BUF) * dom->dz + dom->zs;
  z2=  (k+1-DOM_BUF) * dom->dz + dom->zs;


 wx[0]=(x2-x)/(x2-x1);
 wx[1]=(x-x1)/(x2-x1);
 wy[0]=(y2-y)/(y2-y1);
 wy[1]=(y-y1)/(y2-y1);
 wz[0]=(z2-z)/(z2-z1);
 wz[1]=(z-z1)/(z2-z1);

wg[pp]=0;

  for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
      { 
	weight[ii][jj][kk]=wx[ii]*wy[jj]*wz[kk];
ic=i+ii;
jc=j+jj;
kc=k+kk;
periodic_grid_index(ic,jc,kc,dom,3);
	wg[pp]+=weight[ii][jj][kk]*w[ic + jc*dom->Gfz.s1b + kc*dom->Gfz.s2b];
	}
 }

}


__global__ void point_interp_init(int npoints,point_struct *points,real *ug,real *vg,real *wg,real *lpt_stress_u,real *lpt_stress_v,real *lpt_stress_w,real *scg)
{
 int pp = threadIdx.x + blockIdx.x*blockDim.x;
  if(pp < npoints) {
//real rad=  points[pp].r;
ug[pp]=0;
vg[pp]=0;
wg[pp]=0;

lpt_stress_u[pp]=0;
lpt_stress_v[pp]=0;
lpt_stress_w[pp]=0;

scg[pp]=0;
//volPoint[pp]=PI*4/3*rad*rad*rad;
}
}

__global__ void drag_points(point_struct *points, int npoints,
real *ug,real *vg,real *wg,
real *lpt_stress_u,real *lpt_stress_v,real *lpt_stress_w,real *scg,
real rho_f,real mu, g_struct g,gradP_struct gradP,
real C_add,real C_stress,real C_drag,
real sc_eq,real DIFF)
//gradP serve as bodyforce for the time being
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;

  if(pp < npoints) {
real up=points[pp].u;
real vp=points[pp].v;
real wp=points[pp].w;
real dia=2*points[pp].r;
real rhod=points[pp].rho;//rhod is the particle density

//particle interaction force
real iFx=points[pp].iFx;
real iFy=points[pp].iFy;
real iFz=points[pp].iFz;

//fluid velocity at the particle position
real uf=ug[pp];
real vf=vg[pp];
real wf=wg[pp];

//Fluid stress at the particle position
real stress_x=lpt_stress_u[pp];
real stress_y=lpt_stress_v[pp];
real stress_z=lpt_stress_w[pp];


//realative velocity between point_particle and fluid
real ur=sqrt((up-uf)*(up-uf)+(vp-vf)*(vp-vf)+(wp-wf)*(wp-wf));
real nu=mu/rho_f;

//Particle Reynolds number
real Rep=ur*dia/nu+EPSILON;


//Based on hp*dp/D=2+0.6Re_p^0.5 Sc^{1/3}, ref eq (22) in Oresta&&Prosperetti(2014)
real Nu =2+0.6*sqrt(Rep)*powf(nu/DIFF,1.0/3.0);

//Nu=2;
real hp =Nu*DIFF/dia;
/*
real hh =(0.6*sqrt(Rep)*powf(nu/DIFF,1.0/3.0));
int n3=5;
int n2=1;
real r=n3*DIFF/PI/2;
real strength=r/n2;
if(ur>0.5*strength) printf("\nhh %f %f %f %f\n",hh,Rep,ur,DIFF);
*/

//Modification to stokes theory when Rep>1
real F = 1.0f+0.15f*powf(Rep,0.687f); 
//real F = 1.0f; 
real taud = rhod*dia*dia/(18.0f*mu);
real itau=F/taud;

real volume=1./6. * PI * dia*dia*dia;

//Including soluble part:rhod*volume  and insoluble part: ms
real mp =  rhod *volume + points[pp].ms;
//fluid mass in particle volume
real mf =  rho_f *volume;
real msdot=points[pp].msdot;
real gammar=mp/mf;

//drag force on particle
real drag_x=(uf-up)*itau*mp;
real drag_y=(vf-vp)*itau*mp;
real drag_z=(wf-wp)*itau*mp;


//Total add mass,  -gradP.x~z are the body force on fluid apart from gravity
real add_x=(stress_x/rho_f-gradP.x)*mf;
real add_y=(stress_y/rho_f-gradP.y)*mf;
real add_z=(stress_z/rho_f-gradP.z)*mf;

//Store the fluid force on particle including add mass, drag, fluid force and gravity
//default C_add=0.5; C_stress=1; C_drag=1;
real Fx=(C_add*add_x+C_stress*stress_x*volume+C_drag*drag_x)+(mp-mf)*g.x;
real Fy=(C_add*add_y+C_stress*stress_y*volume+C_drag*drag_y)+(mp-mf)*g.y;
real Fz=(C_add*add_z+C_stress*stress_z*volume+C_drag*drag_z)+(mp-mf)*g.z;

//Store the temp particle acceleration, also including the particle soluble mass change  in the last term!
real udot =(Fx+ iFx -up*msdot)/mp;
real vdot =(Fy+ iFy -vp*msdot)/mp;
real wdot =(Fz+ iFz -wp*msdot)/mp;

//acount for added mass effect since it appears also on the left handside of particle governing equation
if(fabs(C_add-0.5)<EPSILON)
{
      udot = udot/(1+C_add/gammar);
      vdot = vdot/(1+C_add/gammar);
      wdot = wdot/(1+C_add/gammar);
}

//particle acceleration
      points[pp].udot = udot;
      points[pp].vdot = vdot;
      points[pp].wdot = wdot;


/*
Store the fluid force on particle including add mass, drag, fluid force and gravity
added mass effect on particle-fluid force interaction
gravity is not implemented on fluid, but we need to add -mf*g to reaction force
F=d(mp*u)/dt -mp*g, according to eq 6 from dropDiffusionForceImpleBlue.pdf
*/
points[pp].Fx=Fx-C_add*udot*mf-mp*g.x;
points[pp].Fy=Fy-C_add*vdot*mf-mp*g.y;
points[pp].Fz=Fz-C_add*wdot*mf-mp*g.z;


/*
Exchange rate of soluble mass into scalar field
dms/dt = pi *dp^2*hp*(rho_s-rho_{sat})   ref eq (20) in Oresta&&Prosperetti(2014)
*/
if(points[pp].ms>0)  points[pp].msdot= PI*dia*dia*hp*(scg[pp]- sc_eq);
else   points[pp].msdot= 0;

//printf("\nmsdot %f %f %f %f\n",points[pp].msdot,hp,scg[pp],sc_eq);

  }
}



//update point velocity -> 1st step of Eulerian prediction
__global__ void move_points_a(point_struct *points, int npoints,
  real dt)
{

  int pp = threadIdx.x + blockIdx.x*blockDim.x; // point_point_particle number
  
  //real m = 4./3. * PI * points[pp].rho * points[pp].r*points[pp].r*points[pp].r;
  //real dT = dt / dt0;

  if(pp < npoints) {
    // update position
      points[pp].x = points[pp].x0 + 0.5*points[pp].u * dt;
      points[pp].y = points[pp].y0 + 0.5*points[pp].v * dt;
      points[pp].z = points[pp].z0 + 0.5*points[pp].w * dt;

     // update linear velocities
      points[pp].u = points[pp].u0 + 0.5*points[pp].udot * dt;
      points[pp].v = points[pp].v0 + 0.5*points[pp].vdot * dt;
      points[pp].w = points[pp].w0 + 0.5*points[pp].wdot * dt;

     // update soluble mass
      points[pp].ms = points[pp].ms0 + 0.5*points[pp].msdot * dt;
  }
}

//update point velocity -> 2nd step of Eulerian prediction
__global__ void move_points_b(dom_struct *dom,point_struct *points, int npoints,
  real dt)
{

  int pp = threadIdx.x + blockIdx.x*blockDim.x; // point_point_particle number
  
  //real m = 4./3. * PI * points[pp].rho * points[pp].r*points[pp].r*points[pp].r;
  //real dT = dt / dt0;

  if(pp < npoints) {
    // update position
      points[pp].x = points[pp].x0 + points[pp].u * dt;
      points[pp].y = points[pp].y0 + points[pp].v * dt;
      points[pp].z = points[pp].z0 + points[pp].w * dt;

     // update linear velocities
      points[pp].u = points[pp].u0 + points[pp].udot * dt;
      points[pp].v = points[pp].v0 + points[pp].vdot * dt;
      points[pp].w = points[pp].w0 + points[pp].wdot * dt;

     // update soluble mass
      points[pp].ms = points[pp].ms0 + points[pp].msdot * dt;
  

//TODO periodic BC for particles, may need to change in future
periodic_grid_position(points[pp].x,points[pp].y,points[pp].z,dom);

//update old values
      points[pp].x0 = points[pp].x;
      points[pp].y0 = points[pp].y;
      points[pp].z0 = points[pp].z;

      points[pp].u0 = points[pp].u;
      points[pp].v0 = points[pp].v;
      points[pp].w0 = points[pp].w;

      points[pp].ms0 = points[pp].ms;
  }
}







__global__ void reset_flag_u(int *flag_u, dom_struct *dom, BC bc)
{
  int i;    // iterator
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  // flag everything as fluid
  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    for(i = dom->Gfx._isb; i < dom->Gfx._ieb; i++) {
      flag_u[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = 1;
    }
  }

  __syncthreads();

  // flag external boundaries
  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    if(bc.uW != PERIODIC)
      flag_u[dom->Gfx._is + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = 0;
    if(bc.uE != PERIODIC)
      flag_u[dom->Gfx._ie-1 + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = 0;
  }
}

__global__ void reset_flag_v(int *flag_v, dom_struct *dom, BC bc)
{
  int j;    // iterator
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  // flag everything as fluid
  if((tk < dom->Gfy._knb) && (ti < dom->Gfy._inb)) {
    for(j = dom->Gfy._jsb; j < dom->Gfy._jeb; j++) {
      flag_v[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b] = 1;
    }
  }

  __syncthreads();

  // flag external boundaries
  if((tk < dom->Gfy._knb) && (ti < dom->Gfy._inb)) {
    if(bc.vS != PERIODIC)
      flag_v[ti + dom->Gfy._js*dom->Gfy._s1b + tk*dom->Gfy._s2b] = 0;
    if(bc.vN != PERIODIC)
      flag_v[ti + (dom->Gfy._je-1)*dom->Gfy._s1b + tk*dom->Gfy._s2b] = 0;
  }
}

__global__ void reset_flag_w(int *flag_w, dom_struct *dom, BC bc)
{
  int k;    // iterator
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  // flag everything as fluid
  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    for(k = dom->Gfz._ksb; k < dom->Gfz._keb; k++) {
      flag_w[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b] = 1;
    }
  }

  __syncthreads();

  // flag external boundaries
  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    if(bc.wB != PERIODIC)
      flag_w[ti + tj*dom->Gfz._s1b + dom->Gfz._ks*dom->Gfz._s2b] = 0;
    if(bc.wT != PERIODIC)
      flag_w[ti + tj*dom->Gfz._s1b + (dom->Gfz._ke-1)*dom->Gfz._s2b] = 0;
  }
}


__global__ void cage_flag_u_periodic_W(int *flag_u, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfx.jnb && tk < dom->Gfx.knb) {
    flag_u[dom->Gfx.isb + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b]
      = flag_u[(dom->Gfx.ie-2) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_u_periodic_E(int *flag_u, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfx.jnb && tk < dom->Gfx.knb) {
    flag_u[(dom->Gfx.ieb-1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b]
      = flag_u[(dom->Gfx.is+1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_u_periodic_S(int *flag_u, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfx.knb && ti < dom->Gfx.inb) {
    flag_u[ti + dom->Gfx.jsb*dom->Gfx.s1b + tk*dom->Gfx.s2b]
      = flag_u[ti + (dom->Gfx.je-1)*dom->Gfx.s1b + tk*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_u_periodic_N(int *flag_u, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfx.knb && ti < dom->Gfx.inb) {
    flag_u[ti + (dom->Gfx.jeb-1)*dom->Gfx.s1b + tk*dom->Gfx.s2b]
      = flag_u[ti + dom->Gfx.js*dom->Gfx.s1b + tk*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_u_periodic_B(int *flag_u, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfx.inb && tj < dom->Gfx.jnb) {
    flag_u[ti + tj*dom->Gfx.s1b + dom->Gfx.ksb*dom->Gfx.s2b]
      = flag_u[ti + tj*dom->Gfx.s1b + (dom->Gfx.ke-1)*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_u_periodic_T(int *flag_u, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfx.inb && tj < dom->Gfx.jnb) {
    flag_u[ti + tj*dom->Gfx.s1b + (dom->Gfx.keb-1)*dom->Gfx.s2b]
      = flag_u[ti + tj*dom->Gfx.s1b + dom->Gfx.ks*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_v_periodic_W(int *flag_v, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfy.jnb && tk < dom->Gfy.knb) {
    flag_v[dom->Gfy.isb + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b]
      = flag_v[(dom->Gfy.ie-1) + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_v_periodic_E(int *flag_v, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfy.jnb && tk < dom->Gfy.knb) {
    flag_v[(dom->Gfy.ieb-1) + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b]
      = flag_v[dom->Gfy.is + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_v_periodic_S(int *flag_v, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfy.knb && ti < dom->Gfy.inb) {
    flag_v[ti + dom->Gfy.jsb*dom->Gfy.s1b + tk*dom->Gfy.s2b]
      = flag_v[ti + (dom->Gfy.je-2)*dom->Gfy.s1b + tk*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_v_periodic_N(int *flag_v, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfy.knb && ti < dom->Gfy.inb) {
    flag_v[ti + (dom->Gfy.jeb-1)*dom->Gfy.s1b + tk*dom->Gfy.s2b]
      = flag_v[ti + (dom->Gfy.js+1)*dom->Gfy.s1b + tk*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_v_periodic_B(int *flag_v, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfy.inb && tj < dom->Gfy.jnb) {
    flag_v[ti + tj*dom->Gfy.s1b + dom->Gfy.ksb*dom->Gfy.s2b]
      = flag_v[ti + tj*dom->Gfy.s1b + (dom->Gfy.ke-1)*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_v_periodic_T(int *flag_v, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfy.inb && tj < dom->Gfy.jnb) {
    flag_v[ti + tj*dom->Gfy.s1b + (dom->Gfy.keb-1)*dom->Gfy.s2b]
      = flag_v[ti + tj*dom->Gfy.s1b + dom->Gfy.ks*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_w_periodic_W(int *flag_w, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfz.jnb && tk < dom->Gfz.knb) {
    flag_w[dom->Gfz.isb + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b]
      = flag_w[(dom->Gfz.ie-1)+ tj*dom->Gfz.s1b + tk*dom->Gfz.s2b];
  }
}

__global__ void cage_flag_w_periodic_E(int *flag_w, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfz.jnb && tk < dom->Gfz.knb) {
    flag_w[(dom->Gfz.ieb-1) + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b]
      = flag_w[dom->Gfz.is + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b];
  }
}

__global__ void cage_flag_w_periodic_S(int *flag_w, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfz.knb && ti < dom->Gfz.inb) {
    flag_w[ti + dom->Gfz.jsb*dom->Gfz.s1b + tk*dom->Gfz.s2b]
      = flag_w[ti + (dom->Gfz.je-1)*dom->Gfz.s1b + tk*dom->Gfz.s2b];
  }
}

__global__ void cage_flag_w_periodic_N(int *flag_w, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfz.knb && ti < dom->Gfz.inb) {
    flag_w[ti + (dom->Gfz.jeb-1)*dom->Gfz.s1b + tk*dom->Gfz.s2b]
      = flag_w[ti + dom->Gfz.js*dom->Gfz.s1b + tk*dom->Gfz.s2b];
  }
}

__global__ void cage_flag_w_periodic_B(int *flag_w, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfz.inb && tj < dom->Gfz.jnb) {
    flag_w[ti + tj*dom->Gfz.s1b + dom->Gfz.ksb*dom->Gfz.s2b]
      = flag_w[ti + tj*dom->Gfz.s1b + (dom->Gfz.ke-2)*dom->Gfz.s2b];
  }
}

__global__ void cage_flag_w_periodic_T(int *flag_w, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfz.inb && tj < dom->Gfz.jnb) {
    flag_w[ti + tj*dom->Gfz.s1b + (dom->Gfz.keb-1)*dom->Gfz.s2b]
      = flag_w[ti + tj*dom->Gfz.s1b + (dom->Gfz.ks+1)*dom->Gfz.s2b];
  }
}





