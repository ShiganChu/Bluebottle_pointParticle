#include "cuda_bicgstab.h"
#include "bluebottle.h"
#include "shigan.h"
#include "cuda.h"
#include <cusp/dia_matrix.h>

//#include "cuPrintf.cu"
#include <stdio.h>

/*
//TODO verify the process of omega
__global__ void Omega_x(real rho_f, real nu, real *u,real *v,real *w,dom_struct *dom,real *omega_x, real dt, real dt0)
{
// create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  
//  __shared__ real grad_P[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  
// working constants
  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??
  real ab0= 0.5 * a ;
  real ab=1. + ab0 ;
  
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
omega_x[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b]=ddy * (w[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - w[i+(tj-1)*dom->Gcc._s1b + tk*dom->Gcc._s2b])
					      +ddz * (u[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - u[i+tj*dom->Gcc._s1b + (tk-1)*dom->Gcc._s2b]);
  }
}


__global__ void Omega_y(real rho_f, real nu, real *u,real *v,real *w,dom_struct *dom,real *omega_x, real dt, real dt0)
{
// create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  
//  __shared__ real grad_P[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  
// working constants
  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??
  real ab0= 0.5 * a ;
  real ab=1. + ab0 ;
  

  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed


  // loop over u-planes
  for(int j = dom->Gfy._js; j < dom->Gfy._je; j++) {
    int i = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int k = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int ti = threadIdx.x;
    int tk = threadIdx.y;
omega_y[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b]=ddz * (u[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - u[i+tj*dom->Gcc._s1b + (tk-1)*dom->Gcc._s2b])
					      -ddx * (w[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - w[(i-1)+tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]);
}
}


*/


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
  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??
  real ab0= 0.5 * a ;
  real ab=1. + ab0 ;
  

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
      grad_P[tj + tk*blockDim.x]=abs(flag_u[i + tj*dom->Gfx._s1b
        + tk*dom->Gfx._s2b])
        * ddx * (p[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p[(i-1) + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]);
      
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
    }
//if(i==16&&j==8&&k==8) printf("diffusion_u %f, grad_Px %f\n",s_d[tj + tk*blockDim.x],grad_P[tj + tk*blockDim.x]);
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
  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??
  real ab0= 0.5 * a ;
  real ab=1. + ab0 ;
  
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
      grad_P[tk + ti*blockDim.x]=abs(flag_v[ti + j*dom->Gfy._s1b
        + tk*dom->Gfy._s2b])
        * ddy * (p[ti + j*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p[ti + (j-1)*dom->Gcc._s1b + tk*dom->Gcc._s2b]);

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
  //    stress[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] = rho_f*s_d[tk + ti*blockDim.x]-ab*grad_P[tk + ti*blockDim.x]+ab0*grad_P0[tk + ti*blockDim.x]-gradP_y;
      stress[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] = rho_f*s_d[tk + ti*blockDim.x]-ab*grad_P[tk + ti*blockDim.x]+ab0*grad_P0[tk + ti*blockDim.x];
  //  stress[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] = rho_f*s_d[tk + ti*blockDim.x];
//   stress[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] =-ab*grad_P[tk + ti*blockDim.x]+ab0*grad_P0[tk + ti*blockDim.x]-gradP_y;
//if(i==16&&j==8&&k==8) printf("diffusion_v %f, grad_Py %f\n",s_d[tk + ti*blockDim.x],grad_P[tk + ti*blockDim.x]);
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
  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??
  real ab0= 0.5 * a ;
  real ab=1. + ab0 ;
  
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
      grad_P[ti + tj*blockDim.x] = abs(flag_w[ti + tj*dom->Gfz._s1b
        + k*dom->Gfz._s2b])
        * ddz * (p[ti + tj*dom->Gcc._s1b + k*dom->Gcc._s2b]
        - p[ti + tj*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b]);
	
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
//      stress[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b]= rho_f* s_d[ti + tj*blockDim.x]-ab*grad_P[ti + tj*blockDim.x]+ab0*grad_P0[ti + tj*blockDim.x]-gradP_z;
      stress[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b]= rho_f* s_d[ti + tj*blockDim.x]-ab*grad_P[ti + tj*blockDim.x]+ab0*grad_P0[ti + tj*blockDim.x];
//if(i==16&&j==8&&k==8) printf("diffusion_w %f, grad_Pz %f\n",s_d[ti + tj*blockDim.x],grad_P[ti + tj*blockDim.x]);
    }
  }
}




//TODO need to change "C" to "points" and add "npoints" in  Pressure interpolation
//interpolate pressure from grid to point
__global__ void interpolate_point_pre(real *p0, real *p,real *pg,point_struct *points, dom_struct *dom, real dt0, real dt, BC bc)
{
  //int node = threadIdx.x;
  int point = blockIdx.x;

  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz; 
  // convert node (r, theta, phi) to (x, y, z)
 // real xp, yp, zp;  // Cartesian radial vector
 
// Cartesian location of node
real  x =  points[point].x;
real  y =  points[point].y;
real  z =  points[point].z;

  if(x < dom->xs && bc.uW == PERIODIC) x = x + dom->xl;
  else if(x > dom->xe && bc.uE == PERIODIC) x = x - dom->xl;
  if(y < dom->ys && bc.vS == PERIODIC) y = y + dom->yl;
  else if(y > dom->ye && bc.vN == PERIODIC) y = y - dom->yl;
  if(z < dom->zs && bc.wB == PERIODIC) z = z + dom->zl;
  else if(z > dom->ze && bc.wT == PERIODIC) z = z - dom->zl;

  __syncthreads();

  // find index of cell containing node
  int i = floor((x - dom->xs) * ddx) + DOM_BUF;
  int j = floor((y - dom->ys) * ddy) + DOM_BUF;
  int k = floor((z - dom->zs) * ddz) + DOM_BUF;
  if(i < dom->Gcc.is) i = dom->Gcc.ie-1;
  if(j < dom->Gcc.js) j = dom->Gcc.je-1;
  if(k < dom->Gcc.ks) k = dom->Gcc.ke-1;
  if(i > dom->Gcc.ie-1) i = dom->Gcc.is;
  if(j > dom->Gcc.je-1) j = dom->Gcc.js;
  if(k > dom->Gcc.ke-1) k = dom->Gcc.ks;
  int C = i + j*dom->Gcc.s1b + k*dom->Gcc.s2b;

  // Cartesian location of center of cell
  real xx = (i-0.5) * dom->dx + dom->xs;
  real yy = (j-0.5) * dom->dy + dom->ys;
  real zz = (k-0.5) * dom->dz + dom->zs;

  // interpolate pressure
  // TODO: IMPROVE THIS INTERPOLATION METHOD!!!
  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??
  real pc = (1. + 0.5 * a) * p[C] - 0.5 * a * p0[C];
  real pw = (1. + 0.5 * a) * p[C-1] - 0.5 * a * p0[C-1];
  real pe = (1. + 0.5 * a) * p[C+1] - 0.5 * a * p0[C+1];
  real ps = (1. + 0.5 * a) * p[C-dom->Gcc.s1b] - 0.5 * a * p0[C-dom->Gcc.s1b];
  real pn = (1. + 0.5 * a) * p[C+dom->Gcc.s1b] - 0.5 * a * p0[C+dom->Gcc.s1b];
  real pb = (1. + 0.5 * a) * p[C-dom->Gcc.s2b] - 0.5 * a * p0[C-dom->Gcc.s2b];
  real pt = (1. + 0.5 * a) * p[C+dom->Gcc.s2b] - 0.5 * a * p0[C+dom->Gcc.s2b];
  real dpdx = 0.5 * (pe - pw) * ddx;
  real dpdy = 0.5 * (pn - ps) * ddy;
  real dpdz = 0.5 * (pt - pb) * ddz;
//  pp[node+nnodes*point] = pc + dpdx*(x-xx) + dpdy*(y-yy) + dpdz*(z-zz);
  pg[C] = pc + dpdx*(x-xx) + dpdy*(y-yy) + dpdz*(z-zz);
 
}

//interpolate pressure from grid to point
__global__ void interpolate_point_pre_shigan(real *p0, real *p,real *pg,int npoints,point_struct *points, dom_struct *dom, real dt0, real dt, BC bc)
{
  //int node = threadIdx.x;
  int point = blockIdx.x;
if(point<npoints)
{
  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz; 
  // convert node (r, theta, phi) to (x, y, z)
 // real xp, yp, zp;  // Cartesian radial vector
 
// Cartesian location of node
real  x =  points[point].x;
real  y =  points[point].y;
real  z =  points[point].z;

  if(x < dom->xs && bc.uW == PERIODIC) x = x + dom->xl;
  else if(x > dom->xe && bc.uE == PERIODIC) x = x - dom->xl;
  if(y < dom->ys && bc.vS == PERIODIC) y = y + dom->yl;
  else if(y > dom->ye && bc.vN == PERIODIC) y = y - dom->yl;
  if(z < dom->zs && bc.wB == PERIODIC) z = z + dom->zl;
  else if(z > dom->ze && bc.wT == PERIODIC) z = z - dom->zl;


   real x1,x2,y1,y2,z1,z2;
   int ic,jc,kc;
int i,j,k;

// interpolate u-velocity
  i = floor((x - dom->xs) * ddx- 0.5) + DOM_BUF;//for xm[i]
  j = floor((y - dom->ys) * ddy- 0.5) + DOM_BUF;//for ym[j]
  k = floor((z - dom->zs) * ddz- 0.5) + DOM_BUF;//for zm[k]
  if(i < dom->Gcc.is) i = dom->Gcc.ie-1;
  if(j < dom->Gcc.js) j = dom->Gcc.je-1;
  if(k < dom->Gcc.ks) k = dom->Gcc.ke-1;
  if(i > dom->Gcc.ie-1) i = dom->Gcc.is;
  if(j > dom->Gcc.je-1) j = dom->Gcc.js;
  if(k > dom->Gcc.ke-1) k = dom->Gcc.ks;

  x1 = (i-0.5) * dom->dx + dom->xs;
  x2 = (i+0.5) * dom->dx + dom->xs;
  y1=  (j-0.5) * dom->dy + dom->ys;
  y2=  (j+0.5) * dom->dy + dom->ys;
  z1=  (k-0.5) * dom->dz + dom->zs;
  z2=  (k+0.5) * dom->dz + dom->zs;

real wx[2],wy[2],wz[2];
real weight[2][2][2];
int  CC[2][2][2];

 wx[0]=(x2-x)/(x2-x1);
 wx[1]=(x-x1)/(x2-x1);
 wy[0]=(y2-y)/(y2-y1);
 wy[1]=(y-y1)/(y2-y1);
 wz[0]=(z2-z)/(z2-z1);
 wz[1]=(z-z1)/(z2-z1);

  // interpolate pressure
  // TODO: IMPROVE THIS INTERPOLATION METHOD!!!
  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??


pg[point]=0;
  for(int kk = 0; kk < 2; kk++)
   for(int jj = 0; jj < 2; jj++)
    for(int ii = 0; ii < 2; ii++)
      {

weight[ii][jj][kk]=wx[ii]*wy[jj]*wz[kk];
ic=i+ii;
jc=j+jj;
kc=k+kk;
  if(ic < dom->Gfx.is) ic = dom->Gfx.ie-1;
  if(jc < dom->Gfx.js) jc = dom->Gfx.je-1;
  if(kc < dom->Gfx.ks) kc = dom->Gfx.ke-1;
  if(ic > dom->Gfx.ie-1) ic = dom->Gfx.is;
  if(jc > dom->Gfx.je-1) jc = dom->Gfx.js;
  if(kc > dom->Gfx.ke-1) kc = dom->Gfx.ks;
      CC[ii][jj][kk]=ic +jc*dom->Gfx.s1b + kc*dom->Gfx.s2b;
        pg[point]+=weight[ii][jj][kk]*((1. + 0.5 * a)*p[CC[ii][jj][kk]]-0.5 * a*p0[CC[ii][jj][kk]]);
        }
   }  
}


//Bilinear interpolation!
//TODO treat boudary specification!!! such as mask!!!
__global__ void interpolate_point_vel_shigan(real *u, real *v, real *w,int npoints,real rho_f, real nu,real *ug,real *vg,real *wg, point_struct *points, dom_struct *dom, real dt0, real dt, BC bc)
{
  // int node = threadIdx.x;
  int point =  threadIdx.x + blockIdx.x*blockDim.x;

  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;

if(point<npoints)
{
// Cartesian location of node
real  x =  points[point].x;
real  y =  points[point].y;
real  z =  points[point].z;

  if(x < dom->xs && bc.uW == PERIODIC) x = x + dom->xl;
  else if(x > dom->xe && bc.uE == PERIODIC) x = x - dom->xl;
  if(y < dom->ys && bc.vS == PERIODIC) y = y + dom->yl;
  else if(y > dom->ye && bc.vN == PERIODIC) y = y - dom->yl;
  if(z < dom->zs && bc.wB == PERIODIC) z = z + dom->zl;
  else if(z > dom->ze && bc.wT == PERIODIC) z = z - dom->zl;

  __syncthreads();

int i,j,k;   
  real x1,x2,y1,y2,z1,z2;
  int ic,jc,kc;


  // interpolate velocities
  // interpolate u-velocity
  i = floor((x - dom->xs) * ddx) + DOM_BUF; 	//for x[i]
  j = floor((y - dom->ys) * ddy- 0.5) + DOM_BUF;//for ym[j]
  k = floor((z - dom->zs) * ddz- 0.5) + DOM_BUF;//for zm[k]
  if(i < dom->Gfx.is) i = dom->Gfx.ie-1;
  if(j < dom->Gfx.js) j = dom->Gfx.je-1;
  if(k < dom->Gfx.ks) k = dom->Gfx.ke-1;
  if(i > dom->Gfx.ie-1) i = dom->Gfx.is;
  if(j > dom->Gfx.je-1) j = dom->Gfx.js;
  if(k > dom->Gfx.ke-1) k = dom->Gfx.ks;

  x1 = (i-DOM_BUF) * dom->dx + dom->xs;
  x2 = (i+1-DOM_BUF) * dom->dx + dom->xs;
  y1=  (j-0.5) * dom->dy + dom->ys;
  y2=  (j+0.5) * dom->dy + dom->ys;
  z1=  (k-0.5) * dom->dz + dom->zs;
  z2=  (k+0.5) * dom->dz + dom->zs;

real wx[2],wy[2],wz[2];
real weight[2][2][2];
int  CC[2][2][2];
 
 wx[0]=(x2-x)/(x2-x1);
 wx[1]=(x-x1)/(x2-x1);
 wy[0]=(y2-y)/(y2-y1);
 wy[1]=(y-y1)/(y2-y1);
 wz[0]=(z2-z)/(z2-z1);
 wz[1]=(z-z1)/(z2-z1);

ug[point]=0;
  for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
      { 
	
weight[ii][jj][kk]=wx[ii]*wy[jj]*wz[kk];	
ic=i+ii;
jc=j+jj;
kc=k+kk;
  if(ic < dom->Gfx.is) ic = dom->Gfx.ie-1;
  if(jc < dom->Gfx.js) jc = dom->Gfx.je-1;
  if(kc < dom->Gfx.ks) kc = dom->Gfx.ke-1;
  if(ic > dom->Gfx.ie-1) ic = dom->Gfx.is;
  if(jc > dom->Gfx.je-1) jc = dom->Gfx.js;
  if(kc > dom->Gfx.ke-1) kc = dom->Gfx.ks;
      CC[ii][jj][kk]=ic +jc*dom->Gfx.s1b + kc*dom->Gfx.s2b;
	ug[point]+=weight[ii][jj][kk]*u[CC[ii][jj][kk]];
	}
/*
real x0=x;
real y0=y;
real halfH=2;
real F=1;
real omega=0.5*(PI/halfH)*(PI/halfH);
real nu=1;
real k_r=sqrt(omega/2/nu);
real H=halfH;

real zm=y0;
real PP=cos(k_r*(H+zm))*cosh(k_r*(H-zm))+cos(k_r*(H-zm))*cosh(k_r*(H+zm));
real QQ=sin(k_r*(H+zm))*sinh(k_r*(H-zm))+sin(k_r*(H-zm))*sinh(k_r*(H+zm));
real SS=cos(2*H*k_r)+cosh(2*H*k_r);
real u_r=F/omega/SS*(-sin(omega*dt)*(PP-SS)+cos(omega*dt)*QQ);
*/

//printf("shigan_ug %f %d %d %d %f %f %f %f %f %f \n",ug[point],i,j,k,x1,y1,z1,x,y,z);
//printf("shigan_ug %f %f %f %f %f %f\n",ug[point]-u_r,u_r,z1,z,x,dt);
//printf("shigan_xyz %f %f %f %f %f %f %f %f %f\n",x,x1,x2,y,y1,y2,z,z1,z2);

  // interpolate V-velocity
  i = floor((x - dom->xs) * ddx- 0.5) + DOM_BUF;
  j = floor((y - dom->ys) * ddy) + DOM_BUF;
  k = floor((z - dom->zs) * ddz- 0.5) + DOM_BUF;
  if(i < dom->Gfy.is) i = dom->Gfy.ie-1;
  if(j < dom->Gfy.js) j = dom->Gfy.je-1;
  if(k < dom->Gfy.ks) k = dom->Gfy.ke-1;
  if(i > dom->Gfy.ie-1) i = dom->Gfy.is;
  if(j > dom->Gfy.je-1) j = dom->Gfy.js;
  if(k > dom->Gfy.ke-1) k = dom->Gfy.ks;

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

vg[point]=0;

  for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
      { 
	weight[ii][jj][kk]=wx[ii]*wy[jj]*wz[kk];
ic=i+ii;
jc=j+jj;
kc=k+kk;
  if(ic < dom->Gfy.is) ic = dom->Gfy.ie-1;
  if(jc < dom->Gfy.js) jc = dom->Gfy.je-1;
  if(kc < dom->Gfy.ks) kc = dom->Gfy.ke-1;
  if(ic > dom->Gfy.ie-1) ic = dom->Gfy.is;
  if(jc > dom->Gfy.je-1) jc = dom->Gfy.js;
  if(kc > dom->Gfy.ke-1) kc = dom->Gfy.ks;
        CC[ii][jj][kk]=(i+ii) + (j+jj)*dom->Gfy.s1b + (k+kk)*dom->Gfy.s2b;
	vg[point]+=weight[ii][jj][kk]*v[CC[ii][jj][kk]];
	}










  // interpolate W-velocity
  i = floor((x - dom->xs) * ddx- 0.5) + DOM_BUF;
  j = floor((y - dom->ys) * ddy- 0.5) + DOM_BUF;
  k = floor((z - dom->zs) * ddz) + DOM_BUF;
  if(i < dom->Gfz.is) i = dom->Gfz.ie-1;
  if(j < dom->Gfz.js) j = dom->Gfz.je-1;
  if(k < dom->Gfz.ks) k = dom->Gfz.ke-1;
  if(i > dom->Gfz.ie-1) i = dom->Gfz.is;
  if(j > dom->Gfz.je-1) j = dom->Gfz.js;
  if(k > dom->Gfz.ke-1) k = dom->Gfz.ks;

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

wg[point]=0;

  for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
      { 
	weight[ii][jj][kk]=wx[ii]*wy[jj]*wz[kk];
      ic=i+ii;
jc=j+jj;
kc=k+kk;
  if(ic < dom->Gfz.is) ic = dom->Gfz.ie-1;
  if(jc < dom->Gfz.js) jc = dom->Gfz.je-1;
  if(kc < dom->Gfz.ks) kc = dom->Gfz.ke-1;
  if(ic > dom->Gfz.ie-1) ic = dom->Gfz.is;
  if(jc > dom->Gfz.je-1) jc = dom->Gfz.js;
  if(kc > dom->Gfz.ke-1) kc = dom->Gfz.ks;
  CC[ii][jj][kk]=(i+ii) + (j+jj)*dom->Gfz.s1b + (k+kk)*dom->Gfz.s2b;
	wg[point]+=weight[ii][jj][kk]*w[CC[ii][jj][kk]];
	}


//printf("shigan: ug[point],vg[point],wg[point] %f %f %f \n",ug[point],vg[point],wg[point]);
 }

}
















__global__ void interpolate_point_vel(real *u, real *v, real *w,int npoints,real rho_f, real nu,real *ug,real *vg,real *wg, point_struct *points, dom_struct *dom, real dt0, real dt, BC bc)
{
  // int node = threadIdx.x;
  int point =  threadIdx.x + blockIdx.x*blockDim.x;

  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;


  //real uu, vv, ww;  // temporary nodes for Cartesian result of interpolation

  // convert node (r, theta, phi) to (x, y, z)
  //  real xp, yp, zp;  // Cartesian radial vector
if(point<npoints)
{
// Cartesian location of node
real  x =  points[point].x;
real  y =  points[point].y;
real  z =  points[point].z;

  if(x < dom->xs && bc.uW == PERIODIC) x = x + dom->xl;
  else if(x > dom->xe && bc.uE == PERIODIC) x = x - dom->xl;
  if(y < dom->ys && bc.vS == PERIODIC) y = y + dom->yl;
  else if(y > dom->ye && bc.vN == PERIODIC) y = y - dom->yl;
  if(z < dom->zs && bc.wB == PERIODIC) z = z + dom->zl;
  else if(z > dom->ze && bc.wT == PERIODIC) z = z - dom->zl;

  __syncthreads();

  // find index of cell containing node
  int i = floor((x - dom->xs) * ddx) + DOM_BUF;
  int j = floor((y - dom->ys) * ddy) + DOM_BUF;
  int k = floor((z - dom->zs) * ddz) + DOM_BUF;
  if(i < dom->Gcc.is) i = dom->Gcc.ie-1;
  if(j < dom->Gcc.js) j = dom->Gcc.je-1;
  if(k < dom->Gcc.ks) k = dom->Gcc.ke-1;
  if(i > dom->Gcc.ie-1) i = dom->Gcc.is;
  if(j > dom->Gcc.je-1) j = dom->Gcc.js;
  if(k > dom->Gcc.ke-1) k = dom->Gcc.ks;
  int C = i + j*dom->Gcc.s1b + k*dom->Gcc.s2b;

  // Cartesian location of center of cell
  real xx = (i-0.5) * dom->dx + dom->xs;
  real yy = (j-0.5) * dom->dy + dom->ys;
  real zz = (k-0.5) * dom->dz + dom->zs;


  // interpolate velocities
  // don't work with cell-center anymore;
  // find closest cell face in x-direction

  // interpolate u-velocity
  i = round((x - dom->xs) * ddx - 0.5) + DOM_BUF;
  j = floor((y - dom->ys) * ddy) + DOM_BUF;
  k = floor((z - dom->zs) * ddz) + DOM_BUF;
  if(i < dom->Gfx.is) i = dom->Gfx.ie-1;
  if(j < dom->Gfx.js) j = dom->Gfx.je-1;
  if(k < dom->Gfx.ks) k = dom->Gfx.ke-1;
  if(i > dom->Gfx.ie-1) i = dom->Gfx.is;
  if(j > dom->Gfx.je-1) j = dom->Gfx.js;
  if(k > dom->Gfx.ke-1) k = dom->Gfx.ks;
  xx = (i-DOM_BUF) * dom->dx + dom->xs;
  yy = (j-0.5) * dom->dy + dom->ys;
  zz = (k-0.5) * dom->dz + dom->zs;
  C = i + j*dom->Gfx.s1b + k*dom->Gfx.s2b;
  real dudx = 0.5 * (u[C+1] - u[C-1]) * ddx;
  real dudy = 0.5 * (u[C+dom->Gfx.s1b] - u[C-dom->Gfx.s1b]) * ddy;
  real dudz = 0.5 * (u[C+dom->Gfx.s2b] - u[C-dom->Gfx.s2b]) * ddz;
  ug[point] = u[C] + dudx * (x - xx) + dudy * (y - yy) + dudz * (z - zz);

/*
real x0=x;
real y0=y;
real halfH=2;
real F=1;
real omega=0.5*(PI/halfH)*(PI/halfH);
real nu=1;
real k_r=sqrt(omega/2/nu);
real H=halfH;
 
real zm=y0;
real PP=cos(k_r*(H+zm))*cosh(k_r*(H-zm))+cos(k_r*(H-zm))*cosh(k_r*(H+zm));
real QQ=sin(k_r*(H+zm))*sinh(k_r*(H-zm))+sin(k_r*(H-zm))*sinh(k_r*(H+zm));
real SS=cos(2*H*k_r)+cosh(2*H*k_r);
real u_r=F/omega/SS*(-sin(omega*dt)*(PP-SS)+cos(omega*dt)*QQ);
*/
  
//printf("ug      %f %d %d %d %f %f %f %f %f %f \n",ug[point],i,j,k,xx,yy,zz,x,y,z);
//printf("ug      %f %f %f %f %f %f\n",(ug[point]-u_r)/u_r,u_r,zz,z,x,dt);

 // set equal to the relative velocity between this point_particle and the
  // offending point_particle if this node interferes with another point_particle
  // note that node carries this point_particle if intersecting with a wall, so
  // you get zero for the point_particle, which is what we expect

/*
  uu = (points[point].nodes[node]==-1)*uu + (points[point].nodes[node]>-1)
    *(points[points[point].nodes[node]].u - points[point].u);
*/
  // interpolate v-velocity
  i = floor((x - dom->xs) * ddx) + DOM_BUF;
  j = round((y - dom->ys) * ddy - 0.5) + DOM_BUF;
  k = floor((z - dom->zs) * ddz) + DOM_BUF;
  if(i < dom->Gfy.is) i = dom->Gfy.ie-1;
  if(j < dom->Gfy.js) j = dom->Gfy.je-1;
  if(k < dom->Gfy.ks) k = dom->Gfy.ke-1;
  if(i > dom->Gfy.ie-1) i = dom->Gfy.is;
  if(j > dom->Gfy.je-1) j = dom->Gfy.js;
  if(k > dom->Gfy.ke-1) k = dom->Gfy.ks;
  xx = (i-0.5) * dom->dx + dom->xs;
  yy = (j-DOM_BUF) * dom->dy + dom->ys;
  zz = (k-0.5) * dom->dz + dom->zs;
  C = i + j*dom->Gfy.s1b + k*dom->Gfy.s2b;
  real dvdx = 0.5 * (v[C+1] - v[C-1]) * ddx;
  real dvdy = 0.5 * (v[C+dom->Gfy.s1b] - v[C-dom->Gfy.s1b]) * ddy;
  real dvdz = 0.5 * (v[C+dom->Gfy.s2b] - v[C-dom->Gfy.s2b]) * ddz;
  vg[point] = v[C] + dvdx * (x - xx) + dvdy * (y - yy) + dvdz * (z - zz);

//printf("ug[C],vg[C],wg[C],v[C],C %f %f %f %f %d\n",ug[C],vg[C],wg[C],v[C],C);
//printf("Vg[C],V[C],dvdx,ddx,x,xx %f %f %f %f %f %f\n",vg[C],v[C],dvdx,ddx,x,xx);
  // interpolate w-velocity
  i = floor((x - dom->xs) * ddx) + DOM_BUF;
  j = floor((y - dom->ys) * ddy) + DOM_BUF;
  k = round((z - dom->zs) * ddz - 0.5) + DOM_BUF;
  if(i < dom->Gfz.is) i = dom->Gfz.ie-1;
  if(j < dom->Gfz.js) j = dom->Gfz.je-1;
  if(k < dom->Gfz.ks) k = dom->Gfz.ke-1;
  if(i > dom->Gfz.ie-1) i = dom->Gfz.is;
  if(j > dom->Gfz.je-1) j = dom->Gfz.js;
  if(k > dom->Gfz.ke-1) k = dom->Gfz.ks;
  xx = (i-0.5) * dom->dx + dom->xs;
  yy = (j-0.5) * dom->dy + dom->ys;
  zz = (k-DOM_BUF) * dom->dz + dom->zs;
  C = i + j*dom->Gfz.s1b + k*dom->Gfz.s2b;
  real dwdx = 0.5 * (w[C+1] - w[C-1]) * ddx;
  real dwdy = 0.5 * (w[C+dom->Gfz.s1b] - w[C-dom->Gfz.s1b]) * ddy;
  real dwdz = 0.5 * (w[C+dom->Gfz.s2b] - w[C-dom->Gfz.s2b]) * ddz;
  wg[point] = w[C] + dwdx * (x - xx) + dwdy * (y - yy) + dwdz * (z - zz);

//printf("ug[point],vg[point],wg[point] %f %f %f \n",ug[point],vg[point],wg[point]);
 }

}
__global__ void drag_points(point_struct *points, int npoints,
real *ug,real *vg,real *wg,
real *lpt_stress_u,real *lpt_stress_v,real *lpt_stress_w,
real rho_f,real mu, g_struct g,gradP_struct gradP,
real C_add,real C_stress,real C_drag)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;

  if(i < npoints) {
real up=points[i].u;
real vp=points[i].v;
real wp=points[i].w;
real dia=2*points[i].r;
real rhod=points[i].rho;//rhod is the particle density

real uf=ug[i];
real vf=vg[i];
real wf=wg[i];

//Fluid stress at the particle position
real stress_x=lpt_stress_u[i];
real stress_y=lpt_stress_v[i];
real stress_z=lpt_stress_w[i];


//realative velocity between point_particle and fluid
real ur=sqrt((up-uf)*(up-uf)+(vp-vf)*(vp-vf)+(wp-wf)*(wp-wf));

//printf("up uf %f %f",up,uf);

real Rep=rho_f*ur*dia/mu+EPSILON;

//drag force on the point_particle
real F = 1.0f+0.15f*powf(Rep,0.687f); 
//real F = 1.0f; 


real taud = rhod*dia*dia/(18.0f*mu);
real itau=F/taud;
//printf("rho_f is %f ur is %f, dia %f %f %f",rho_f,ur,dia,mu,EPSILON);

real volume=1./6. * PI * dia*dia*dia;
//real m =  rhod *volume;
real drag_x=(uf-up)*itau*rhod;
real drag_y=(vf-vp)*itau*rhod;
real drag_z=(wf-wp)*itau*rhod;


//real C_drag=1.0;
//real C_drag=0;
//real C_stress=1.0;
//real C_stress=0;

// Add mass effect, by shigan 10_2_2014
// Add mass effect of relative acceleration
//real C_add=0.5;
//real C_add=0.0;


/*
real C1=1;
real C2=10;
real add_x=(C2*stress_x-C1*rho_f*gradP.x);
real add_y=(C2*stress_y-C1*rho_f*gradP.y);
real add_z=(C2*stress_z-C1*rho_f*gradP.z);
*/

//////Total add mass
real add_x=(stress_x-rho_f*gradP.x);
real add_y=(stress_y-rho_f*gradP.y);
real add_z=(stress_z-rho_f*gradP.z);


points[i].Fx=(C_add*add_x+C_stress*stress_x+C_drag*drag_x)*volume;
points[i].Fy=(C_add*add_y+C_stress*stress_y+C_drag*drag_y)*volume;
points[i].Fz=(C_add*add_z+C_stress*stress_z+C_drag*drag_z)*volume;


//printf("stress&position %f %f %f, %f %f %f\n",stress_x,stress_y,stress_z,points[i].x,points[i].y,points[i].z);
//printf("uf,up,wp,itau,taud %f %f %f %f\n",uf,up,wp,itau,taud);
  }
}



/*
__global__ void addMass_points(point_struct *points, int npoints,real *ug,real *vg,real *wg,real rho_f,real mu)
{
int i = threadIdx.x + blockIdx.x*blockDim.x;

  if(i < npoints) {
real up=points[i].u;
real vp=points[i].v;
real wp=points[i].w;
real dia=2*points[i].r;
real rhod=points[i].rho;//rhod is the particle density

real uf=ug[i];
real vf=vg[i];
real wf=wg[i];





real m = 1./6. * PI * rhod * dia*dia*dia;

points[i].Fx=(uf-up)*itau*m;
points[i].Fy=(vf-vp)*itau*m;
points[i].Fz=(wf-wp)*itau*m;

printf("uf,up,itau,taud %f %f %f %f\n",uf,up,itau,taud);
  }
}
*/

//// fluid stress on the particles   du/dt
//__global__ void  stress_points(real *ug, real *vg, real *wg,real *ug0, real *vg0, real *wg0,real *conv_ug, real *conv_vg, real *conv_wg,point_struct *points, int npoints,real rho_f, real dt)
//
//{
//  int i = threadIdx.x + blockIdx.x*blockDim.x;
//
//  if(i < npoints) {
//real rhod=points[i].rho;
//real dia=2*points[i].r;
//real m = 1./6. * PI * rhod * dia*dia*dia;
//real ratio=rho_f/rhod*m;
//
//
//real dudt=(ug[i]-ug0[i])/dt+conv_ug[i];
//real dvdt=(vg[i]-vg0[i])/dt+conv_vg[i];
//real dwdt=(wg[i]-wg0[i])/dt+conv_wg[i];
//
//
////Whether will accumulate??  Fluid stress
//points[i].Fx+=dudt*ratio;
//points[i].Fy+=dvdt*ratio;
//points[i].Fz+=dwdt*ratio;
//
////Added Mass
//real addMassRatio=ratio*0.5;
//points[i].Fx+=(dudt-points[i].udot)*addMassRatio;
//points[i].Fy+=(dvdt-points[i].vdot)*addMassRatio;
//points[i].Fz+=(dwdt-points[i].wdot)*addMassRatio;
//
//  }
//}


//update point velocity
__global__ void move_points_b(dom_struct *dom, point_struct *points, int npoints,
  real dt, real dt0, g_struct g, real rho_f, real ttime,
real C_add,real C_stress,real C_drag)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x; // point_point_particle number
  real rho_p=points[pp].rho;
  real r=points[pp].r;
  real m = 4./3. * PI * rho_p*r*r*r;

  //real m = 4./3. * PI * points[pp].rho * points[pp].r*points[pp].r*points[pp].r;
  //real dT = dt / dt0;

  if(pp < npoints) {

      // update linear accelerations
      points[pp].udot = (points[pp].Fx +  points[pp].iFx ) / m
        + (points[pp].rho - rho_f) / points[pp].rho * g.x;
      points[pp].vdot = (points[pp].Fy + points[pp].iFy  ) / m
        + (points[pp].rho - rho_f) / points[pp].rho * g.y;
      points[pp].wdot = (points[pp].Fz + points[pp].iFz  ) / m
        + (points[pp].rho - rho_f) / points[pp].rho * g.z;


if((C_add-0.5)<EPSILON)
{
real gammar=rho_p/rho_f;
      points[pp].udot = points[pp].udot/(1+C_add/gammar);
      points[pp].vdot = points[pp].vdot/(1+C_add/gammar);
      points[pp].wdot = points[pp].wdot/(1+C_add/gammar);
}
/*
      points[pp].udot = points[pp].Fx + (points[pp].rho - rho_f) / points[pp].rho * g.x;
      points[pp].vdot = points[pp].Fy + (points[pp].rho - rho_f) / points[pp].rho * g.y;
      points[pp].wdot = points[pp].Fz + (points[pp].rho - rho_f) / points[pp].rho * g.z;
*/

//printf("points[pp].Fz,g.z %f %f\n",points[pp].Fz,g.z);
 


     // update linear velocities
      points[pp].u = points[pp].u0 + points[pp].udot * dt;
      points[pp].v = points[pp].v0 + points[pp].vdot * dt;
      points[pp].w = points[pp].w0 + points[pp].wdot * dt;

      points[pp].u0 = points[pp].u;
      points[pp].v0 = points[pp].v;
      points[pp].w0 = points[pp].w;

      // update position
      points[pp].x = points[pp].x + points[pp].u * dt;
      if(points[pp].x < dom->xs) points[pp].x = points[pp].x + dom->xl;
      else if(points[pp].x > dom->xe) points[pp].x = points[pp].x - dom->xl;
      points[pp].y = points[pp].y + points[pp].v * dt;
      if(points[pp].y < dom->ys) points[pp].y = points[pp].y + dom->yl;
      else if(points[pp].y > dom->ye) points[pp].y = points[pp].y - dom->yl;
      points[pp].z = points[pp].z + points[pp].w * dt;
      if(points[pp].z < dom->zs) points[pp].z = points[pp].z + dom->zl;
      else if(points[pp].z > dom->ze) points[pp].z = points[pp].z - dom->zl;
    
/*
    if(points[pp].rotating) {
      // update angular accelerations
      real I = 0.4 * m * points[pp].r*points[pp].r;
      points[pp].oxdot = (points[pp].Lx + points[pp].aLx) / I;
      points[pp].oydot = (points[pp].Ly + points[pp].aLy) / I;
      points[pp].ozdot = (points[pp].Lz + points[pp].aLz) / I;

      // update angular velocities
      points[pp].ox = points[pp].ox + points[pp].oxdot * dt;
      points[pp].oy = points[pp].oy + points[pp].oydot * dt;
      points[pp].oz = points[pp].oz + points[pp].ozdot * dt;

      // update basis vectors
      // calculate rotation magnitude
      real mag = sqrt(points[pp].ox*points[pp].ox + points[pp].oy*points[pp].oy
        + points[pp].oz*points[pp].oz);
      // calculate normalized rotation axis
      real X = 0;
      real Y = 0;
      real Z = 0;
      if(mag > 0) {
        X = points[pp].ox / mag;
        Y = points[pp].oy / mag;
        Z = points[pp].oz / mag;
      }
      // calculate rotation quaternion
      real theta = mag * dt;
      real qr = cos(0.5*theta);
      real qi = X * sin(0.5*theta);
      real qj = Y * sin(0.5*theta);
      real qk = Z * sin(0.5*theta);
      // compute quaternion conjugation to apply rotation to basis vectors
      rotate(qr, qi, qj, qk, &points[pp].axx, &points[pp].axy, &points[pp].axz);
      rotate(qr, qi, qj, qk, &points[pp].ayx, &points[pp].ayy, &points[pp].ayz);
      rotate(qr, qi, qj, qk, &points[pp].azx, &points[pp].azy, &points[pp].azz);
    }
*/
  }
}







