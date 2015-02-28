#include "cuda_scalar.h"
#include "cuda_point.h"


//Using central difference scheme in space, and adam-bashforth scheme in time.
__global__ void diffScalar_explicitD(real *sc,real *sc0,
dom_struct *dom, real DIFF_dt)
{
//store scalar value of difference direction at cell center
__shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM]; // sc at bottom center
__shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM]; // sc at top center
__shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM]; // sc at center
real ddx = 1. / dom->dx; // to limit the number of divisions needed
real ddy = 1. / dom->dy; // to limit the number of divisions needed
real ddz = 1. / dom->dz; // to limit the number of divisions needed
int C;
// loop over z-planes
// for(int k = dom->Gcc._ksb; k < dom->Gcc._keb; k++) {
for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
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
if( j < dom->Gcc._jeb&& i < dom->Gcc._ieb) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];
}
// make sure all threads complete shared memory copy
__syncthreads();
// compute convective term
// if off the shared memory block boundary
if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) {
//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;
real sc_diff= (ddsdxx+ddsdyy+ddsdzz)*DIFF_dt;
sc_c[ti + tj*blockDim.x] +=sc_diff;
}
// copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
if((j >= dom->Gcc._js && j < dom->Gcc._je)
&& (i >= dom->Gcc._is && i < dom->Gcc._ie)
&& (ti > 0 && ti < (blockDim.x-1))
&& (tj > 0 && tj < (blockDim.y-1))) {
C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
sc[C] = sc_c[ti + tj*blockDim.x];
}
}
}


//2nd order for source as well as for convective&diffusion term
__global__ void diffScalar_rhs_CN(real DIFF_dt, 
 real *sc0, real *sc_rhs, dom_struct *dom)
{

//store scalar value of difference direction at cell center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top    center
  __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at        center

  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed
   int C;
  // loop over w-planes
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) { //k from 1~nz
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int i = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int ti = threadIdx.x;
    int tj = threadIdx.y;

   //j range in 0~ny+1; i range in 0~nx+1
    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];//k-1 from 0~nz-1
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];//k+1 from 2~nz+1

 }

  // make sure all threads complete shared memory copy
    __syncthreads();

 
  // compute convective term
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) {
 
//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;

      real diff_sc= (ddsdxx+ddsdyy+ddsdzz)*DIFF_dt;
     // sc_c[ti + tj*blockDim.x] +=(1-theta)*diff_sc;
      sc_c[ti + tj*blockDim.x] +=diff_sc;
//s_rhs[ti + tj*blockDim.x] =sc_c[ti + tj*blockDim.x];
}

// copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
     && (i >= dom->Gcc._is && i < dom->Gcc._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      sc_rhs[C] = sc_c[ti + tj*blockDim.x];
    }
  }
}




//Diff_dt is the length scale square; and it is equall to (delta_f^2-dx^2)/16/ln2
__global__ void diffScalar_coeffs_CN(real thetaD, dom_struct *dom, int pitch,
  real *values,  int *flag_u, int *flag_v, int *flag_w,int coordiSys)
{
  int i;  // iterator
  int C;  // cell locations
 // int W, E, S, N, B, T;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);
  real ddy = 1. / (dom->dy * dom->dy);
  real ddz = 1. / (dom->dz * dom->dz);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;


//get domain start and end index
int incGhost=0;
int is,js,ks,ie,je,ke;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);

int index; 
    index=coordiSys*48        +18;
int s1=tex1Dfetch(texRefDomInfo,index);
int s2=tex1Dfetch(texRefDomInfo,index+1);

  // loop over slices to set values
  if(tj < je && tk < ke) {
    for(i = is; i < ie; i++) {
      C = (i-DOM_BUF) + (tj-DOM_BUF)*s1 + (tk-DOM_BUF)*s2;
//real dx=dom->dx;
//real theta=0.5-1/12.f/DIFF_dt*dx*dx;
//real theta=0.6;
//real theta=0.f;
//real thetaD=theta*DIFF_dt;
real buf=0.f;
      values[C + pitch * 1]  -= thetaD*ddz;
      values[C + pitch * 3]  -= thetaD*ddy;
      values[C + pitch * 5]  -= thetaD*ddx;
      buf  += 1.;
      buf  += 2*thetaD*ddx;
      buf  += 2*thetaD*ddy;
      buf  += 2*thetaD*ddz;
      values[C + pitch * 6]  += buf;
      values[C + pitch * 7]  -= thetaD*ddx;
      values[C + pitch * 9]  -= thetaD*ddy;
      values[C + pitch * 11] -= thetaD*ddz;
/*
      values[C + pitch * 1]  -= 0.5*DIFF_dt*ddz;
      values[C + pitch * 3]  -= 0.5*DIFF_dt*ddy;
      values[C + pitch * 5]  -= 0.5*DIFF_dt*ddx;
      values[C + pitch * 6]  += 1.;
      values[C + pitch * 6]  += DIFF_dt*ddx;
      values[C + pitch * 6]  += DIFF_dt*ddy;
      values[C + pitch * 6]  += DIFF_dt*ddz;
      values[C + pitch * 7]  -= 0.5*DIFF_dt*ddx;
      values[C + pitch * 9]  -= 0.5*DIFF_dt*ddy;
      values[C + pitch * 11] -= 0.5*DIFF_dt*ddz;
*/

    }
  }
 
}

 

__global__ void copy_diffSc_noghost(real *sc_noghost, real *sc_ghost, dom_struct *dom,int coordiSys)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

//get domain start and end index
int incGhost=0;
int is,js,ks,ie,je,ke;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);

int index; 
    index=coordiSys*48        +18;
int s1=tex1Dfetch(texRefDomInfo,index);
int s2=tex1Dfetch(texRefDomInfo,index+1);
    index=coordiSys*48        +21;
int s1b=tex1Dfetch(texRefDomInfo,index);
int s2b=tex1Dfetch(texRefDomInfo,index+1);


  if(tj < je-DOM_BUF && tk < ke-DOM_BUF) {
    for(int i = is-DOM_BUF; i < ie-DOM_BUF; i++) {
      sc_noghost[i + tj*s1 + tk*s2]
        = sc_ghost[(i+DOM_BUF) + (tj+DOM_BUF)*s1b
        + (tk+DOM_BUF)*s2b];
    }
  }
}


 
__global__ void diffScalar_coeffs_init(dom_struct *dom, int pitch, real *values,int coordiSys)
{
  int i;  // iterator
  int C;  // cell

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

//get domain start and end index
int incGhost=0;
int is,js,ks,ie,je,ke;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);

//int in=ie-is;
int jn=je-js;
int kn=ke-ks;

//get s1&s2
int   index=coordiSys*48        +18;
int s1=tex1Dfetch(texRefDomInfo,index);
int s2=tex1Dfetch(texRefDomInfo,index+1);


  // loop over all slices to initialize to zeint is,js,ks,ie,je,ke;
  if(tj < jn && tk < kn) {
    for(i = is-DOM_BUF; i < ie-DOM_BUF; i++) {
      C = i + tj*s1 + tk*s2;
      values[C + 0*pitch]  = 0.;
      values[C + 1*pitch]  = 0.;
      values[C + 2*pitch]  = 0.;
      values[C + 3*pitch]  = 0.;
      values[C + 4*pitch]  = 0.;
      values[C + 5*pitch]  = 0.;
      values[C + 6*pitch]  = 0.;
      values[C + 7*pitch]  = 0.;
      values[C + 8*pitch]  = 0.;
      values[C + 9*pitch]  = 0.;
      values[C + 10*pitch] = 0.;
      values[C + 11*pitch] = 0.;
      values[C + 12*pitch] = 0.;
    }
  }
}




__global__ void diffScalar_coeffs_periodic_W(real DIFF_dt, dom_struct *dom,
  int pitch, real *values,int coordiSys)
{
//get domain start and end index
int incGhost=0;
int is,js,ks,ie,je,ke;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);

 
//get s1&s2
int   index=coordiSys*48        +18;
int s1=tex1Dfetch(texRefDomInfo,index);
int s2=tex1Dfetch(texRefDomInfo,index+1);

  int i = is;  // iterator
  int C;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tj < je) && (tk < ke)) {
    C = (i-DOM_BUF) + (tj-DOM_BUF)*s1 + (tk-DOM_BUF)*s2;
    values[C + pitch * 5] += DIFF_dt*ddx;
    values[C + pitch * 8] -= DIFF_dt*ddx;
  }
}

__global__ void diffScalar_coeffs_periodic_E(real DIFF_dt, dom_struct *dom,
  int pitch, real *values,int coordiSys)
{
//get domain start and end index
int incGhost=0;
int is,js,ks,ie,je,ke;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);

 
//get s1&s2
int   index=coordiSys*48        +18;
int s1=tex1Dfetch(texRefDomInfo,index);
int s2=tex1Dfetch(texRefDomInfo,index+1);

  int i = ie-1;  // iterator
  int C;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tj < je) && (tk < ke)) {
    C = (i-DOM_BUF) + (tj-DOM_BUF)*s1 + (tk-DOM_BUF)*s2;
    values[C + pitch * 4] -= DIFF_dt*ddx;
    values[C + pitch * 7] += DIFF_dt*ddx;
  }
}

__global__ void diffScalar_coeffs_periodic_S(real DIFF_dt, dom_struct *dom,
  int pitch, real *values,int coordiSys)
{
//get domain start and end index
int incGhost=0;
int is,js,ks,ie,je,ke;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);

 
//get s1&s2
int   index=coordiSys*48        +18;
int s1=tex1Dfetch(texRefDomInfo,index);
int s2=tex1Dfetch(texRefDomInfo,index+1);

  int j = js;  // iterator
  int C;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tk < ke) && (ti < ie)) {
    C = (ti-DOM_BUF) + (j-DOM_BUF)*s1 + (tk-DOM_BUF)*s2;
    values[C + pitch * 3]  += DIFF_dt*ddy;
    values[C + pitch * 10] -= DIFF_dt*ddy;
  }
}

__global__ void diffScalar_coeffs_periodic_N(real DIFF_dt, dom_struct *dom,
  int pitch, real *values,int coordiSys)
{
//get domain start and end index
int incGhost=0;
int is,js,ks,ie,je,ke;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);

 
//get s1&s2
int   index=coordiSys*48        +18;
int s1=tex1Dfetch(texRefDomInfo,index);
int s2=tex1Dfetch(texRefDomInfo,index+1);

  int j = je-1;  // iterator
  int C;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tk < ke) && (ti < ie)) {
    C = (ti-DOM_BUF) + (j-DOM_BUF)*s1 + (tk-DOM_BUF)*s2;
    values[C + pitch * 2] -= DIFF_dt*ddy;
    values[C + pitch * 9] += DIFF_dt*ddy;
  }
}

__global__ void diffScalar_coeffs_periodic_B(real DIFF_dt, dom_struct *dom,
  int pitch, real *values,int coordiSys)
{
//get domain start and end index
int incGhost=0;
int is,js,ks,ie,je,ke;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);

 
//get s1&s2
int   index=coordiSys*48        +18;
int s1=tex1Dfetch(texRefDomInfo,index);
int s2=tex1Dfetch(texRefDomInfo,index+1);

  int k = ks;  // iterator
  int C;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((ti < ie) && (tj < je)) {
    C = (ti-DOM_BUF) + (tj-DOM_BUF)*s1 + (k-DOM_BUF)*s2;
    values[C + pitch * 1]  += DIFF_dt*ddz;
    values[C + pitch * 12] -= DIFF_dt*ddz;
  }
}

__global__ void diffScalar_coeffs_periodic_T(real DIFF_dt, dom_struct *dom,
  int pitch, real *values,int coordiSys)
{
//get domain start and end index
int incGhost=0;
int is,js,ks,ie,je,ke;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);

 
//get s1&s2
int   index=coordiSys*48        +18;
int s1=tex1Dfetch(texRefDomInfo,index);
int s2=tex1Dfetch(texRefDomInfo,index+1);

  int k = ke-1;  // iterator
  int C;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((ti < ie) && (tj < je)) {
    C = (ti-DOM_BUF) + (tj-DOM_BUF)*s1 + (k-DOM_BUF)*s2;
    values[C + pitch * 0]  -= DIFF_dt*ddz;
    values[C + pitch * 11] += DIFF_dt*ddz;
  }
}



__global__ void copy_sc_noghost(real *sc_noghost, real *sc_ghost, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gcc.je-DOM_BUF && tk < dom->Gcc.ke-DOM_BUF) {
    for(int i = dom->Gcc.is-DOM_BUF; i < dom->Gcc.ie-DOM_BUF; i++) {
      sc_noghost[i + tj*dom->Gcc._s1 + tk*dom->Gcc._s2]
        = sc_ghost[(i+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc._s1b
        + (tk+DOM_BUF)*dom->Gcc._s2b];
    }
  }
}

__global__ void scalar_coeffs_init(dom_struct *dom, int pitch, real *values)
{
  int i;  // iterator
  int C;  // cell

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  // loop over all slices to initialize to zero
  if(tj < dom->Gcc.jn && tk < dom->Gcc.kn) {
    for(i = dom->Gcc.is-DOM_BUF; i < dom->Gcc.ie-DOM_BUF; i++) {
      C = i + tj*dom->Gcc.s1 + tk*dom->Gcc.s2;
      values[C + 0*pitch]  = 0.;
      values[C + 1*pitch]  = 0.;
      values[C + 2*pitch]  = 0.;
      values[C + 3*pitch]  = 0.;
      values[C + 4*pitch]  = 0.;
      values[C + 5*pitch]  = 0.;
      values[C + 6*pitch]  = 0.;
      values[C + 7*pitch]  = 0.;
      values[C + 8*pitch]  = 0.;
      values[C + 9*pitch]  = 0.;
      values[C + 10*pitch] = 0.;
      values[C + 11*pitch] = 0.;
      values[C + 12*pitch] = 0.;
    }
  }
}

//TODO need to change the BC to scalar BC. Besides, the scalar will also need new periodic_BC func
//This coeff matrix depend on the flow boundary condition rather than the scalar BC!!!
__global__ void scalar_coeffs(real DIFF, real dt, dom_struct *dom, int pitch,
  real *values,  int *flag_u, int *flag_v, int *flag_w,real *epsp)
{
  int i;  // iterator
  int C;  // cell locations
  int W, E, S, N, B, T;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);
  real ddy = 1. / (dom->dy * dom->dy);
  real ddz = 1. / (dom->dz * dom->dz);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  int CC;
  real epsf,iepsf;
//  real epsfe,epsfw,epsfs,epsfn,epsfb,epsft;
//  real iepsfe,iepsfw,iepsfs,iepsfn,iepsfb,iepsft;
  // loop over slices to set values
  if(tj < dom->Gcc.jn + DOM_BUF && tk < dom->Gcc.kn + DOM_BUF) {
    for(i = dom->Gcc.is; i < dom->Gcc.ie; i++) {
      W = i + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;
      E = (i+1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;
      S = i + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b;
      N = i + (tj+1)*dom->Gfy.s1b + tk*dom->Gfy.s2b;
      B = i + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b;
      T = i + tj*dom->Gfz.s1b + (tk+1)*dom->Gfz.s2b;

      C = (i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
     // CC = i + tj*dom->Gcc.s1 + tk*dom->Gcc.s2;
      CC = i + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b;

//TODO the iepsf should be changed to face values correspondingly	
	epsf=1-epsp[CC];
        iepsf=__fdiv_rd(1.f,epsf);
        real buf=0.f;

      values[C + pitch * 1]  -= (real)abs(flag_w[B]) *0.5f*DIFF*dt*ddz*iepsf;
      values[C + pitch * 3]  -= (real)abs(flag_v[S]) *0.5f*DIFF*dt*ddy*iepsf;
      values[C + pitch * 5]  -= (real)abs(flag_u[W]) *0.5f*DIFF*dt*ddx*iepsf;
      buf  += 1.;
      buf  += (real)(abs(flag_u[W]) + abs(flag_u[E]))*0.5f *DIFF*dt*ddx*iepsf;
      buf  += (real)(abs(flag_v[S]) + abs(flag_v[N]))*0.5f *DIFF*dt*ddy*iepsf;
      buf  += (real)(abs(flag_w[B]) + abs(flag_w[T]))*0.5f *DIFF*dt*ddz*iepsf;
      values[C + pitch * 6]  += buf;

      values[C + pitch * 7]  -= (real)abs(flag_u[E]) *0.5f*DIFF*dt*ddx*iepsf;
      values[C + pitch * 9]  -= (real)abs(flag_v[N]) *0.5f*DIFF*dt*ddy*iepsf;
      values[C + pitch * 11] -= (real)abs(flag_w[T]) *0.5f*DIFF*dt*ddz*iepsf;

//if(abs(abs(flag_w[B]*flag_w[T]*flag_u[W]*flag_u[E]*flag_v[N]*flag_v[S])-1)>EPSILON) printf("\ni~k %d %d %d\n",i,tj,tk);
    }
  }
}

/*
__global__ void scalar_coeffs(real DIFF, real dt, dom_struct *dom, int pitch,
  real *values)
{
  int i;  // iterator
  int C;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);
  real ddy = 1. / (dom->dy * dom->dy);
  real ddz = 1. / (dom->dz * dom->dz);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  // loop over slices to set values
  if(tj < dom->Gcc.jn + DOM_BUF && tk < dom->Gcc.kn + DOM_BUF) {
    for(i = dom->Gcc.is; i < dom->Gcc.ie; i++) {
      C = (i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
      values[C + pitch * 1]  -= 0.5*DIFF*dt*ddz;
      values[C + pitch * 3]  -= 0.5*DIFF*dt*ddy;
      values[C + pitch * 5]  -= 0.5*DIFF*dt*ddx;
      values[C + pitch * 6]  += 1.;
      values[C + pitch * 6]  += DIFF*dt*ddx;
      values[C + pitch * 6]  += DIFF*dt*ddy;
      values[C + pitch * 6]  += DIFF*dt*ddz;
      values[C + pitch * 7]  -= 0.5*DIFF*dt*ddx;
      values[C + pitch * 9]  -= 0.5*DIFF*dt*ddy;
      values[C + pitch * 11] -= 0.5*DIFF*dt*ddz;
    }
  }
}
*/


__global__ void copy_sc_ghost(real *sc_ghost, real *sc_noghost, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gcc.je-DOM_BUF && tk < dom->Gcc.ke-DOM_BUF) {
    for(int i = dom->Gcc.is-DOM_BUF; i < dom->Gcc.ie-DOM_BUF; i++) {
      sc_ghost[(i+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc._s1b
        + (tk+DOM_BUF)*dom->Gcc._s2b] = sc_noghost[i + tj*dom->Gcc._s1
        + tk*dom->Gcc._s2];
    }
  }
}

//2nd order for source as well as for convective&diffusion term
__global__ void scalar_rhs_upwind_1st(real rho_f, real DIFF, 
real *u, real *v, real *w, 
real *epsp,real *epsp0,
real *f, real *f0,
real *conv0,real *conv,real *diff,
real *sc0, real *sc_rhs, dom_struct *dom, 
real dt, real dt0)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
 
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_c0[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv0
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // source term
  __shared__ real s_f0[MAX_THREADS_DIM * MAX_THREADS_DIM];       // source term
  __shared__ real s_rhs[MAX_THREADS_DIM * MAX_THREADS_DIM];  // solution

//store scalar value of difference direction at cell center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top    center
  __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at        center

//store point volume percentage in each grid cell
  __shared__ real p_epsp[MAX_THREADS_DIM * MAX_THREADS_DIM];
  __shared__ real p_epsp0[MAX_THREADS_DIM * MAX_THREADS_DIM];

  // working constants
  real ab0 = 0.5 * dt / dt0;   // for Adams-Bashforth stepping
  real ab = 1. + ab0;          // for Adams-Bashforth stepping
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed
   int C;
  // loop over w-planes
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) { //k from 1~nz
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int i = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int ti = threadIdx.x;
    int tj = threadIdx.y;

   //j range in 0~ny+1; i range in 0~nx+1
    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];//k-1 from 0~nz-1
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];//k+1 from 2~nz+1

s_c0[ti + tj*blockDim.x]=conv0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
s_f[ti + tj*blockDim.x]=    f[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
s_f0[ti + tj*blockDim.x]=   f0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];

p_epsp[ti + tj*blockDim.x]= epsp[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
p_epsp0[ti + tj*blockDim.x]= epsp0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
}




  // make sure all threads complete shared memory copy
    __syncthreads();

//No boundary value of u,v,w are used
    if((j >= dom->Gfz._js && j < dom->Gfz._je)
      && (i >= dom->Gfz._is && i < dom->Gfz._ie)) {
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];//k from 1~nz
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];//k from 2~nz+1
    }
    if((j >= dom->Gfx._js && j < dom->Gfx._je)
      && (i >= dom->Gfx._is && i < dom->Gfx._ie)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._js && j < dom->Gfy._je)
      && (i >= dom->Gfy._is && i < dom->Gfy._ie)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }

   s_rhs[ti + tj*blockDim.x] = 0.0;  
    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute convective term
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) {
//scalar on the west and east u face
real u_c=(s_u[ti + tj*blockDim.x]+s_u[(ti+1) + tj*blockDim.x])/2;// at xm(i)
real v_c=(s_v[ti + tj*blockDim.x]+s_v[ti + (tj+1)*blockDim.x])/2;
real w_c=(s_w1[ti + tj*blockDim.x]+s_w0[ti + tj*blockDim.x])/2;//average at cell center, from 1~nz
 
//scalar at face center
real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0f;
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0f;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0f;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0f;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;//average at face center, from 1~nz
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;//average at face center, from 2~nz+1

//1st order upwind scheme
real dsudx,dsvdy,dswdz;
if(u_c>0)	dsudx= (sc_c[ti + tj*blockDim.x]*u_c-sc_uw*s_u[ti + tj*blockDim.x]) *2*ddx;//backward
else		dsudx=-(sc_c[ti + tj*blockDim.x]*u_c-sc_ue*s_u[ti+1 + tj*blockDim.x]) *2* ddx;//forward

if(v_c>0)	dsvdy= (sc_c[ti + tj*blockDim.x]*v_c-sc_vs*s_v[ti + tj*blockDim.x])*2*ddy;//backward
else		dsvdy=-(sc_c[ti + tj*blockDim.x]*v_c-sc_vn*s_v[ti + (tj+1)*blockDim.x])*2*ddy;//forward

if(w_c>0)	dswdz= (sc_c[ti + tj*blockDim.x]*w_c-sc_wb*s_w0[ti + tj*blockDim.x])*2*ddz;//backward
else		dswdz=-(sc_c[ti + tj*blockDim.x]*w_c-sc_wt*s_w1[ti + tj*blockDim.x])*2*ddz;//forward

//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;

     //move convective term to right hand side 
      s_c[ti + tj*blockDim.x]=-(dsudx+dsvdy+dswdz);
     //Diffusion term
      s_d[ti + tj*blockDim.x]= (ddsdxx+ddsdyy+ddsdzz)*DIFF;

//if(isnan(sc_c[ti + tj*blockDim.x])) printf("\nconv_diff %f %f %f %f %f %d %d %d %d %d\n",dsudx,ddsdxx,ddx,ab,sc_c[ti + tj*blockDim.x],ti,tj,i,j,k);
//Adam-bashforth at t+1/2 step
real sc_conv,sc_f,sc_epsp;
if(dt0>0)
{
sc_conv= ab * s_c[ti + tj*blockDim.x]-ab0*s_c0[ti + tj*blockDim.x];
sc_f=    ab * s_f[ti + tj*blockDim.x]-ab0*s_f0[ti + tj*blockDim.x];
sc_epsp= ab * p_epsp[ti + tj*blockDim.x]-ab0*p_epsp0[ti + tj*blockDim.x];
}
else
{
 sc_conv=s_c[ti + tj*blockDim.x];
 sc_f=s_f[ti + tj*blockDim.x];
 sc_epsp=p_epsp[ti + tj*blockDim.x];
}

//Crank-Nicolson method for convective scheme
real sc_diff= 0.5 * s_d[ti + tj*blockDim.x];

//advance scalar. Take the particle volume fraction into consideration
real iepsf=__fdiv_rd(1.f,(1-sc_epsp));
real rhs=(sc_diff+sc_f)*iepsf;
//if(rhs<0) printf("\nrhs error: %d %d %d %f %f %f\n",i,j,k,sc_diff,s_f[ti + tj*blockDim.x],iepsf);

s_d[ti + tj*blockDim.x]=s_d[ti + tj*blockDim.x]*iepsf;

//real rhs=sc_diff+s_f[ti + tj*blockDim.x];
//real rhs=s_f[ti + tj*blockDim.x];
//sc_c[ti + tj*blockDim.x] +=(sc_conv+rhs)*dt;

s_rhs[ti + tj*blockDim.x] =(sc_conv+rhs)*dt;
//if(sc_conv+rhs<0) printf("\nrhs error: %d %d %d %f %f\n",i,j,k,sc_conv,rhs);
s_rhs[ti + tj*blockDim.x] +=sc_c[ti + tj*blockDim.x];
}

// copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
     && (i >= dom->Gcc._is && i < dom->Gcc._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      sc_rhs[C] = s_rhs[ti + tj*blockDim.x];
      conv[C] = s_c[ti + tj*blockDim.x];
      diff[C] =s_d[ti + tj*blockDim.x];
    }
  }
}



__global__ void scalar_rhs_FTCS(real rho_f, real DIFF, real *u, real *v, real *w,  real *epsp, real *f, real *conv0,real *conv,real *diff,real *sc0, real *sc_rhs, dom_struct *dom, real dt, real dt0)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
/*
  __shared__ real s_w2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward
  __shared__ real s_u01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u back
  __shared__ real s_u12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u forward
  __shared__ real s_v01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v back
  __shared__ real s_v12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v forward
*/

  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_c0[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv0
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // source term
  __shared__ real s_rhs[MAX_THREADS_DIM * MAX_THREADS_DIM];  // solution

//store scalar value of difference direction at cell center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top    center
  __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at        center

//store point volume percentage in each grid cell
  __shared__ real p_epsp[MAX_THREADS_DIM * MAX_THREADS_DIM];

  // working constants
  real ab0 = 0.5 * dt / dt0;   // for Adams-Bashforth stepping
  real ab = 1. + ab0;          // for Adams-Bashforth stepping
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed
   int C;
  // loop over w-planes
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int i = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int ti = threadIdx.x;
    int tj = threadIdx.y;


    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];

s_c0[ti + tj*blockDim.x]=conv0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
s_f[ti + tj*blockDim.x]=    f[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];

p_epsp[ti + tj*blockDim.x]= epsp[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
}




  // make sure all threads complete shared memory copy
    __syncthreads();

//No boundary value of u,v,w are used
    if((j >= dom->Gfz._js && j < dom->Gfz._je)
      && (i >= dom->Gfz._is && i < dom->Gfz._ie)) {
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
    }
    if((j >= dom->Gfx._js && j < dom->Gfx._je)
      && (i >= dom->Gfx._is && i < dom->Gfx._ie)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._js && j < dom->Gfy._je)
      && (i >= dom->Gfy._is && i < dom->Gfy._ie)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }

   s_rhs[ti + tj*blockDim.x] = 0.0;  
    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute convective term
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) {
//scalar at face center
real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0f;
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0f;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0f;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0f;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;

//convective term of scalar at cell center
real dsudx=(s_u[(ti+1) + tj*blockDim.x]*sc_ue - s_u[ti + tj*blockDim.x]*sc_uw) * ddx;
real dsvdy=(s_v[ti + (tj+1)*blockDim.x]*sc_vn - s_v[ti + tj*blockDim.x]*sc_vs) * ddy;
real dswdz=(s_w1[ti + tj*blockDim.x]   *sc_wt - s_w0[ti + tj*blockDim.x]*sc_wb) * ddz;

//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;

     //move convective term to right hand side 
      s_c[ti + tj*blockDim.x]=-(dsudx+dsvdy+dswdz);
     //Diffusion term
      s_d[ti + tj*blockDim.x]= (ddsdxx+ddsdyy+ddsdzz)*DIFF;

//if(isnan(sc_c[ti + tj*blockDim.x])) printf("\nconv_diff %f %f %f %f %f %d %d %d %d %d\n",dsudx,ddsdxx,ddx,ab,sc_c[ti + tj*blockDim.x],ti,tj,i,j,k);
//Adam-bashforth at t+1/2 step
real sc_conv;
if(dt0>0)
sc_conv= ab * s_c[ti + tj*blockDim.x]-ab0*s_c0[ti + tj*blockDim.x];
else sc_conv=s_c[ti + tj*blockDim.x];


//Crank-Nicolson method for convective scheme
real sc_diff= 0.5 * s_d[ti + tj*blockDim.x];

//advance scalar. Take the particle volume fraction into consideration
//real rhs=(sc_diff+s_f[ti + tj*blockDim.x])/(1-p_epsp[ti + tj*blockDim.x]);
real rhs=sc_diff+s_f[ti + tj*blockDim.x];
//real rhs=s_f[ti + tj*blockDim.x];
//sc_c[ti + tj*blockDim.x] +=(sc_conv+rhs)*dt;

s_rhs[ti + tj*blockDim.x] =(sc_conv+rhs)*dt;
s_rhs[ti + tj*blockDim.x] +=sc_c[ti + tj*blockDim.x];
}

// copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
     && (i >= dom->Gcc._is && i < dom->Gcc._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      sc_rhs[C] = s_rhs[ti + tj*blockDim.x];
      conv[C] = s_c[ti + tj*blockDim.x];
      diff[C] =2*s_d[ti + tj*blockDim.x];
//if(isnan( s_rhs[ti + tj*blockDim.x])||isinf( s_rhs[ti + tj*blockDim.x])||isnan(s_c[ti + tj*blockDim.x])||isinf(s_c[ti + tj*blockDim.x])||isnan(s_d[ti + tj*blockDim.x])||isinf(s_d[ti + tj*blockDim.x])) printf("\nwrong out %d %d %d %f %f %f\n",i,j,k,s_rhs[ti + tj*blockDim.x],s_c[ti + tj*blockDim.x],s_d[ti + tj*blockDim.x]);

    }
  }
}

   
__global__ void scalar_coeffs_periodic_W(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values,real *epsp,int *flag_u)
{
  int i = dom->Gcc.is;  // iterator
    int C,CC,W;  // cell locations
  real epsf,iepsf;

  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tj < (dom->Gcc.jn + DOM_BUF)) && (tk < (dom->Gcc.kn + DOM_BUF))) {
    C = (i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
   
    //CC = i + tj*dom->Gcc.s1 + tk*dom->Gcc.s2;
    CC = i + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b;
    W = i + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;

    epsf=1-epsp[CC];
    iepsf=__fdiv_rd(1.f,epsf);

    values[C + pitch * 5] +=  (real)abs(flag_u[W])*0.5f*DIFF*iepsf*dt*ddx;
    values[C + pitch * 8] -=  (real)abs(flag_u[W])*0.5f*DIFF*iepsf*dt*ddx;
  }
}

__global__ void scalar_coeffs_periodic_E(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values,real *epsp,int *flag_u)
{
  int i = dom->Gcc.ie-1;  // iterator
    int C,E;  // cell locations
  int CC;
  real epsf,iepsf;

  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tj < (dom->Gcc.jn + DOM_BUF)) && (tk < (dom->Gcc.kn + DOM_BUF))) {
    C = (i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
   // CC = i + tj*dom->Gcc.s1 + tk*dom->Gcc.s2;
    CC = i + tj*dom->Gcc.s1b + tk*dom->Gcc.s2b;
    E = (i+1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;

    epsf=1-epsp[CC];
    iepsf=__fdiv_rd(1.f,epsf);
    values[C + pitch * 4] -=(real)abs(flag_u[E])* 0.5f*DIFF*iepsf*dt*ddx;
    values[C + pitch * 7] +=(real)abs(flag_u[E])* 0.5f*DIFF*iepsf*dt*ddx;
  }
}

__global__ void scalar_coeffs_periodic_S(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values,real *epsp,int *flag_v)
{
  int j = dom->Gcc.js;  // iterator
    int C,S;  // cell locations
  int CC;
  real epsf,iepsf;

  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tk < (dom->Gcc.kn + DOM_BUF)) && (ti < (dom->Gcc.in + DOM_BUF))) {
    C = (ti-DOM_BUF) + (j-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    //CC = ti + j*dom->Gcc.s1 + tk*dom->Gcc.s2;
    CC = ti + j*dom->Gcc.s1b + tk*dom->Gcc.s2b;
   S = ti + j*dom->Gfy.s1b + tk*dom->Gfy.s2b;

    epsf=1-epsp[CC];
    iepsf=__fdiv_rd(1.f,epsf);
    values[C + pitch * 3]  +=(real)abs(flag_v[S])* 0.5f*DIFF*iepsf*dt*ddy;
    values[C + pitch * 10] -=(real)abs(flag_v[S])* 0.5f*DIFF*iepsf*dt*ddy;
  }
}

__global__ void scalar_coeffs_periodic_N(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values,real *epsp,int *flag_v)
{
  int j = dom->Gcc.je-1;  // iterator
    int C,N;  // cell locations
  int CC;
  real epsf,iepsf;

  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tk < (dom->Gcc.kn + DOM_BUF)) && (ti < (dom->Gcc.in + DOM_BUF))) {
    C = (ti-DOM_BUF) + (j-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    //CC = ti + j*dom->Gcc.s1 + tk*dom->Gcc.s2;
    CC = ti + j*dom->Gcc.s1b + tk*dom->Gcc.s2b;
    N = ti + (j+1)*dom->Gfy.s1b + tk*dom->Gfy.s2b;
    epsf=1-epsp[CC];
    iepsf=__fdiv_rd(1.f,epsf);
    values[C + pitch * 2] -=(real)abs(flag_v[N])* 0.5f*DIFF*iepsf*dt*ddy;
    values[C + pitch * 9] +=(real)abs(flag_v[N])* 0.5f*DIFF*iepsf*dt*ddy;
  }
}

__global__ void scalar_coeffs_periodic_B(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values,real *epsp,int *flag_w)
{
  int k = dom->Gcc.ks;  // iterator
    int C,B;  // cell locations
  int CC;
  real epsf,iepsf;

  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((ti < (dom->Gcc.in + DOM_BUF)) && (tj < (dom->Gcc.jn + DOM_BUF))) {
    C = (ti-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (k-DOM_BUF)*dom->Gcc.s2;
    //CC = ti + tj*dom->Gcc.s1 + k*dom->Gcc.s2;
    CC = ti + tj*dom->Gcc.s1b + k*dom->Gcc.s2b;
    B = ti + tj*dom->Gfz.s1b + k*dom->Gfz.s2b;
    epsf=1-epsp[CC];
    iepsf=__fdiv_rd(1.f,epsf);
    values[C + pitch * 1]  +=(real)abs(flag_w[B])* 0.5f*DIFF*iepsf*dt*ddz;
    values[C + pitch * 12] -=(real)abs(flag_w[B])* 0.5f*DIFF*iepsf*dt*ddz;
  }
}

__global__ void scalar_coeffs_periodic_T(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values,real *epsp,int *flag_w)
{
  int k = dom->Gcc.ke-1;  // iterator
    int C,T;  // cell locations
  int CC;
  real epsf,iepsf;

  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((ti < (dom->Gcc.in + DOM_BUF)) && (tj < (dom->Gcc.jn + DOM_BUF))) {
    C = (ti-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (k-DOM_BUF)*dom->Gcc.s2;
    //CC = ti + tj*dom->Gcc.s1 + k*dom->Gcc.s2;
    CC = ti + tj*dom->Gcc.s1b + k*dom->Gcc.s2b;
    T = ti + tj*dom->Gfz.s1b + (k+1)*dom->Gfz.s2b;

    epsf=1-epsp[CC];
    iepsf=__fdiv_rd(1.f,epsf);
    values[C + pitch * 0]  -=(real)abs(flag_w[T])* 0.5f*DIFF*iepsf*dt*ddz;
    values[C + pitch * 11] +=(real)abs(flag_w[T])* 0.5f*DIFF*iepsf*dt*ddz;
  }
}

/*
__global__ void scalar_coeffs_periodic_W(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gcc.is;  // iterator
  int C;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tj < (dom->Gcc.jn + DOM_BUF)) && (tk < (dom->Gcc.kn + DOM_BUF))) {
    C = (i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    values[C + pitch * 5] += 0.5*DIFF*dt*ddx;
    values[C + pitch * 8] -= 0.5*DIFF*dt*ddx;
  }
}

__global__ void scalar_coeffs_periodic_E(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int i = dom->Gcc.ie-1;  // iterator
  int C;  // cell locations
  real ddx = 1. / (dom->dx * dom->dx);

  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tj < (dom->Gcc.jn + DOM_BUF)) && (tk < (dom->Gcc.kn + DOM_BUF))) {
    C = (i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    values[C + pitch * 4] -= 0.5*DIFF*dt*ddx;
    values[C + pitch * 7] += 0.5*DIFF*dt*ddx;
  }
}

__global__ void scalar_coeffs_periodic_S(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gcc.js;  // iterator
  int C;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tk < (dom->Gcc.kn + DOM_BUF)) && (ti < (dom->Gcc.in + DOM_BUF))) {
    C = (ti-DOM_BUF) + (j-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    values[C + pitch * 3]  += 0.5*DIFF*dt*ddy;
    values[C + pitch * 10] -= 0.5*DIFF*dt*ddy;
  }
}

__global__ void scalar_coeffs_periodic_N(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int j = dom->Gcc.je-1;  // iterator
  int C;  // cell locations
  real ddy = 1. / (dom->dy * dom->dy);

  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((tk < (dom->Gcc.kn + DOM_BUF)) && (ti < (dom->Gcc.in + DOM_BUF))) {
    C = (ti-DOM_BUF) + (j-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2;
    values[C + pitch * 2] -= 0.5*DIFF*dt*ddy;
    values[C + pitch * 9] += 0.5*DIFF*dt*ddy;
  }
}

__global__ void scalar_coeffs_periodic_B(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gcc.ks;  // iterator
  int C;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((ti < (dom->Gcc.in + DOM_BUF)) && (tj < (dom->Gcc.jn + DOM_BUF))) {
    C = (ti-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (k-DOM_BUF)*dom->Gcc.s2;
    values[C + pitch * 1]  += 0.5*DIFF*dt*ddz;
    values[C + pitch * 12] -= 0.5*DIFF*dt*ddz;
  }
}

__global__ void scalar_coeffs_periodic_T(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values)
{
  int k = dom->Gcc.ke-1;  // iterator
  int C;  // cell locations
  real ddz = 1. / (dom->dz * dom->dz);

  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if((ti < (dom->Gcc.in + DOM_BUF)) && (tj < (dom->Gcc.jn + DOM_BUF))) {
    C = (ti-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (k-DOM_BUF)*dom->Gcc.s2;
    values[C + pitch * 0]  -= 0.5*DIFF*dt*ddz;
    values[C + pitch * 11] += 0.5*DIFF*dt*ddz;
  }
}

*/



// scalar; west; periodic
__global__ void BC_sc_W_P(real *sc, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    sc[dom->Gcc._isb + tj*s1b + tk*s2b] = sc[(dom->Gcc._ie-1) + tj*s1b + tk*s2b];
}

// scalar; west; Neumann
__global__ void BC_sc_W_N(real *sc, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    sc[dom->Gcc._isb + tj*s1b + tk*s2b] = sc[dom->Gcc._is + tj*s1b + tk*s2b];
}

// scalar; east; periodic
__global__ void BC_sc_E_P(real *sc, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    sc[(dom->Gcc._ieb-1) + tj*s1b + tk*s2b] = sc[dom->Gcc._is + tj*s1b + tk*s2b];
}

// scalar; east; Neumann
__global__ void BC_sc_E_N(real *sc, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    sc[(dom->Gcc._ieb-1) + tj*s1b + tk*s2b] = sc[(dom->Gcc._ie-1)
      + tj*s1b + tk*s2b];
}

// scalar; south; periodic
__global__ void BC_sc_S_P(real *sc, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
    sc[ti + dom->Gcc._jsb*s1b + tk*s2b] = sc[ti + (dom->Gcc._je-1)*s1b + tk*s2b];
}

// scalar; south; Neumann
__global__ void BC_sc_S_N(real *sc, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
    sc[ti + dom->Gcc._jsb*s1b + tk*s2b] = sc[ti + dom->Gcc._js*s1b + tk*s2b];
}

// scalar; north; periodic
__global__ void BC_sc_N_P(real *sc, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
    sc[ti + (dom->Gcc._jeb-1)*s1b + tk*s2b] = sc[ti + dom->Gcc._js*s1b + tk*s2b];
}

// scalar; north; Neumann
__global__ void BC_sc_N_N(real *sc, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
    sc[ti + (dom->Gcc._jeb-1)*s1b + tk*s2b] = sc[ti
      + (dom->Gcc._je-1)*s1b + tk*s2b];
}

// scalar; bottom; periodic
__global__ void BC_sc_B_P(real *sc, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    sc[ti + tj*s1b + dom->Gcc._ksb*s2b] = sc[ti + tj*s1b + (dom->Gcc._ke-1)*s2b];
}

// scalar; bottom; Neumann
__global__ void BC_sc_B_N(real *sc, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b; 
  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    sc[ti + tj*s1b + dom->Gcc._ksb*s2b] = sc[ti + tj*s1b + dom->Gcc._ks*s2b];
}

// scalar; top; periodic
__global__ void BC_sc_T_P(real *sc, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    sc[ti + tj*s1b + (dom->Gcc._keb-1)*s2b] = sc[ti + tj*s1b + dom->Gcc._ks*s2b];
}

// scalar; top; Neumann
__global__ void BC_sc_T_N(real *sc, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    sc[ti + tj*s1b + (dom->Gcc._keb-1)*s2b] = sc[ti
      + tj*s1b + (dom->Gcc._ke-1)*s2b];
}


// scalar; west; Dirichlet
__global__ void BC_sc_W_D(real *sc, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    sc[dom->Gcc._isb + tj*s1b + tk*s2b] = 2*bc-sc[dom->Gcc._is + tj*s1b + tk*s2b];
}


// scalar; east; Dirichlet
__global__ void BC_sc_E_D(real *sc, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    sc[(dom->Gcc._ieb-1) + tj*s1b + tk*s2b] = 2*bc-sc[(dom->Gcc._ie-1) + tj*s1b + tk*s2b];
}

// scalar; south; Dirichlet
__global__ void BC_sc_S_D(real *sc, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
     sc[ti + dom->Gcc._jsb*s1b + tk*s2b] = 2 * bc - sc[ti + dom->Gcc._js*s1b + tk*s2b];
}

// scalar; north; Dirichlet
__global__ void BC_sc_N_D(real *sc, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
     sc[ti + (dom->Gcc._jeb-1)*s1b + tk*s2b] = 2 * bc - sc[ti + (dom->Gcc._je-1)*s1b + tk*s2b];
}

// scalar; bottom; Dirichlet
__global__ void BC_sc_B_D(real *sc, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    sc[ti + tj*s1b + dom->Gcc._ksb*s2b] = 2 * bc - sc[ti + tj*s1b + dom->Gcc._ks*s2b];

}

// scalar; top; Dirichlet
__global__ void BC_sc_T_D(real *sc, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    sc[ti + tj*s1b + (dom->Gcc._keb-1)*s2b] = 2 * bc - sc[ti + tj*s1b + (dom->Gcc._ke-1)*s2b];
}

//QUICK scheme, this one only deal with periodic BC, Doesn't use Adam-Bashforth for the time being
__global__ void advance_sc_QUICK(real DIFF, 
				real *u, real *v, real *w,
				real *f,real *epsp,
				real *diff, real *conv, 
				real *sc,real *sc0,
				dom_struct *dom, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
//store scalar related values such as convective and diffusive term.
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // source term


//store scalar value of difference direction at cell center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top	 center
 __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at  	 center

  __shared__ real sc_ee[MAX_THREADS_DIM * MAX_THREADS_DIM];       
  __shared__ real sc_ww[MAX_THREADS_DIM * MAX_THREADS_DIM];       
  __shared__ real sc_ss[MAX_THREADS_DIM * MAX_THREADS_DIM];       
  __shared__ real sc_nn[MAX_THREADS_DIM * MAX_THREADS_DIM];       
  __shared__ real sc_bb[MAX_THREADS_DIM * MAX_THREADS_DIM];       
  __shared__ real sc_tt[MAX_THREADS_DIM * MAX_THREADS_DIM];       
 
//store point volume percentage in each grid cell
  __shared__ real p_epsp[MAX_THREADS_DIM * MAX_THREADS_DIM];       
 


  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

 
  int C;

  // loop over z-planes
//  for(int k = dom->Gcc._ksb; k < dom->Gcc._keb; k++) {
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
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

    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];

//Peridic BC
int ww=i-2;
int ee=i+2;
if(ww<0) 		ww=dom->Gcc._ie-2;
if(ee>=dom->Gcc._ie)    ee=dom->Gcc._is+1;

int ss=j-2;
int nn=j+2;
if(ss<0) 		ss=dom->Gcc._je-2;
if(nn>=dom->Gcc._je)    nn=dom->Gcc._js+1;

int bb=k-2;
int tt=k+2;
if(bb<0) 		bb=dom->Gcc._ke-2;
if(tt>=dom->Gcc._ke)    tt=dom->Gcc._ks+1;

sc_ee[ti + tj*blockDim.x]=sc0[ee+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_ww[ti + tj*blockDim.x]=sc0[ww+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_ss[ti + tj*blockDim.x]=sc0[i+ss*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_nn[ti + tj*blockDim.x]=sc0[i+nn*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_bb[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + bb*dom->Gcc._s2b];
sc_tt[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + tt*dom->Gcc._s2b];

s_f[ti + tj*blockDim.x]=    f[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
p_epsp[ti + tj*blockDim.x]= epsp[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
}

   // make sure all threads complete shared memory copy
    __syncthreads();

    if((j >= dom->Gfz._js && j < dom->Gfz._je)
      && (i >= dom->Gfz._is && i < dom->Gfz._ie)) {
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
    }
    if((j >= dom->Gfx._js && j < dom->Gfx._je)
      && (i >= dom->Gfx._is && i < dom->Gfx._ie)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._js && j < dom->Gfy._je)
      && (i >= dom->Gfy._is && i < dom->Gfy._ie)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }

  
    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute convective term
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)){
/*
//scalar at face center
real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0f;
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0f;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0f;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0f;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;
*/

//diffusive term of scalar at cell center--CenterScheme
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;

//Convective term using QUICK scheme,  How to treat the ti+2,ti-2 etc in shared memory????? 
real ue_plus=fmax(s_u[ti+1 + tj*blockDim.x],0);
real ue_minus=fmin(s_u[ti+1 + tj*blockDim.x],0);
real uw_plus=fmax(s_u[ti + tj*blockDim.x],0);
real uw_minus=fmin(s_u[ti + tj*blockDim.x],0);

real vn_plus=fmax(s_v[ti + (tj+1)*blockDim.x],0);
real vn_minus=fmin(s_v[ti + (tj+1)*blockDim.x],0);
real vs_plus=fmax(s_v[ti + tj*blockDim.x],0);
real vs_minus=fmin(s_v[ti + tj*blockDim.x],0);

real wt_plus=fmax(s_w1[ti + tj*blockDim.x],0);
real wt_minus=fmin(s_w1[ti + tj*blockDim.x],0);
real wb_plus=fmax(s_w0[ti + tj*blockDim.x],0);
real wb_minus=fmin(s_w0[ti + tj*blockDim.x],0);

real phi_ep= 0.75f*sc_c[ti + tj*blockDim.x]+0.375f*sc_c[ti+1 + tj*blockDim.x]-0.125f*sc_c[ti-1 + tj*blockDim.x];
real phi_wp= 0.75f*sc_c[ti-1 + tj*blockDim.x]+0.375f*sc_c[ti + tj*blockDim.x]-0.125f*sc_ww[ti + tj*blockDim.x];
real phi_wn= 0.75f*sc_c[ti + tj*blockDim.x]+0.375f*sc_c[ti-1 + tj*blockDim.x]-0.125f*sc_c[ti+1 + tj*blockDim.x];
real phi_en= 0.75f*sc_c[ti+1 + tj*blockDim.x]+0.375f*sc_c[ti + tj*blockDim.x]-0.125f*sc_ee[ti + tj*blockDim.x];

real phi_np= 0.75f*sc_c[ti + tj*blockDim.x]+0.375f*sc_c[ti + (tj+1)*blockDim.x]-0.125f*sc_c[ti + (tj-1)*blockDim.x];
real phi_sp= 0.75f*sc_c[ti + (tj-1)*blockDim.x]+0.375f*sc_c[ti + tj*blockDim.x]-0.125f*sc_ss[ti + tj*blockDim.x];
real phi_sn= 0.75f*sc_c[ti + tj*blockDim.x]+0.375f*sc_c[ti + (tj-1)*blockDim.x]-0.125f*sc_c[ti + (tj+1)*blockDim.x];
real phi_nn= 0.75f*sc_c[ti + (tj+1)*blockDim.x]+0.375f*sc_c[ti + tj*blockDim.x]-0.125f*sc_nn[ti + tj*blockDim.x];

real phi_tp= 0.75f*sc_c[ti + tj*blockDim.x]+0.375f*sc_t[ti + tj*blockDim.x]-0.125f*sc_b[ti + tj*blockDim.x];
real phi_bp= 0.75f*sc_b[ti + tj*blockDim.x]+0.375f*sc_c[ti + tj*blockDim.x]-0.125f*sc_bb[ti + tj*blockDim.x];
real phi_bn= 0.75f*sc_c[ti + tj*blockDim.x]+0.375f*sc_b[ti + tj*blockDim.x]-0.125f*sc_t[ti + tj*blockDim.x];
real phi_tn= 0.75f*sc_t[ti + tj*blockDim.x]+0.375f*sc_c[ti + tj*blockDim.x]-0.125f*sc_tt[ti + tj*blockDim.x];

real dsudx=((ue_plus*phi_ep+ue_minus*phi_en)-(uw_plus*phi_wp+uw_minus*phi_wn))*ddx;
real dsvdy=((vn_plus*phi_np+vn_minus*phi_nn)-(vs_plus*phi_sp+vs_minus*phi_sn))*ddy;
real dswdz=((wt_plus*phi_tp+wt_minus*phi_tn)-(wb_plus*phi_bp+wb_minus*phi_bn))*ddz;

     //move convective term to right hand side 
      s_c[ti + tj*blockDim.x]=-(dsudx+dsvdy+dswdz);
     //Diffusion term
      s_d[ti + tj*blockDim.x]= (ddsdxx+ddsdyy+ddsdzz)*DIFF;


real sc_conv= s_c[ti + tj*blockDim.x];
real sc_diff= s_d[ti + tj*blockDim.x];

//TODO should we have s_f at t+1/2 step as well? 
//advance scalar. Take the particle volume fraction into consideration
//real rhs=(sc_diff+s_f[ti + tj*blockDim.x])/(1-p_epsp[ti + tj*blockDim.x]);
real rhs=sc_diff+s_f[ti + tj*blockDim.x];
//real rhs=s_f[ti + tj*blockDim.x];
sc_c[ti + tj*blockDim.x] +=(sc_conv+rhs)*dt;

//if(i==30&&j==30&&k==30) printf("\nsc_diff %f %f %f %f\n",sc_diff,sc_conv,s_f[ti + tj*blockDim.x],dt);
}

// copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
     && (i >= dom->Gcc._is && i < dom->Gcc._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      sc[C] = sc_c[ti + tj*blockDim.x];
      conv[C] = s_c[ti + tj*blockDim.x];
      diff[C] = s_d[ti + tj*blockDim.x];
    }
  }
}

//MacCormac scheme 
__global__ void advance_sc_macCormack(	real DIFF, 
				real *u, real *v, real *w,
				real *f,real *epsp,
				real *diff, real *conv, 
				real *sc,real *sc0,
				dom_struct *dom, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
//store scalar related values such as convective and diffusive term.
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // source term


//store scalar value of difference direction at cell center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top	 center
  __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at  	 center

//store point volume percentage in each grid cell
  __shared__ real p_epsp[MAX_THREADS_DIM * MAX_THREADS_DIM];       
 


  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

 
  int C;

  // loop over z-planes
//  for(int k = dom->Gcc._ksb; k < dom->Gcc._keb; k++) {
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
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

    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];

s_f[ti + tj*blockDim.x]=    f[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
p_epsp[ti + tj*blockDim.x]= epsp[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
}

   // make sure all threads complete shared memory copy
    __syncthreads();

    if((j >= dom->Gfz._js && j < dom->Gfz._je)
      && (i >= dom->Gfz._is && i < dom->Gfz._ie)) {
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
    }
    if((j >= dom->Gfx._js && j < dom->Gfx._je)
      && (i >= dom->Gfx._is && i < dom->Gfx._ie)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._js && j < dom->Gfy._je)
      && (i >= dom->Gfy._is && i < dom->Gfy._ie)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }

  
    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute convective term
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)){
//scalar at face center
real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0f;
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0f;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0f;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0f;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;

//convective term of scalar at cell center
real dsudx=(s_u[(ti+1) + tj*blockDim.x]*sc_ue - s_u[ti + tj*blockDim.x]*sc_uw) * ddx;
real dsvdy=(s_v[ti + (tj+1)*blockDim.x]*sc_vn - s_v[ti + tj*blockDim.x]*sc_vs) * ddy;
real dswdz=(s_w1[ti + tj*blockDim.x]   *sc_wt - s_w0[ti + tj*blockDim.x]*sc_wb) * ddz;


//velocity  at cell center
real u_c=(s_u[(ti+1) + tj*blockDim.x] + s_u[ti + tj*blockDim.x])/2.0f;
real v_c=(s_v[ti + (tj+1)*blockDim.x]+s_v[ti + tj*blockDim.x])/2.0f;
real w_c=(s_w1[ti + tj*blockDim.x]+s_w0[ti + tj*blockDim.x])/2.0f;

//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;

//Artifical diffusivity for MacCormack scheme,  Using the value on x(i),xm(i),x(i+1)
real ddsdxx_num=(sc_uw*s_u[ti + tj*blockDim.x]*s_u[ti + tj*blockDim.x]+sc_ue*s_u[ti+1 + tj*blockDim.x]*s_u[ti+1 + tj*blockDim.x]-2*u_c*u_c*sc_c[ti + tj*blockDim.x]) *ddx* ddx*2;
real ddsdyy_num=(sc_vs*s_v[ti + tj*blockDim.x]*s_v[ti + tj*blockDim.x]+sc_vn*s_v[ti + (tj+1)*blockDim.x]*s_v[ti + (tj+1)*blockDim.x]-2*v_c*v_c*sc_c[ti + tj*blockDim.x]) *ddy* ddy*2;
real ddsdzz_num=(sc_wb*s_w0[ti + tj*blockDim.x]*s_w0[ti + tj*blockDim.x]+sc_wt*s_w1[ti + tj*blockDim.x]*s_w1[ti + tj*blockDim.x]-2*w_c*w_c*sc_c[ti + tj*blockDim.x]) *ddz* ddz*2;

     //move convective term to right hand side 
      s_c[ti + tj*blockDim.x]=-(dsudx+dsvdy+dswdz);
     //Diffusion term
real diff_center=(ddsdxx+ddsdyy+ddsdzz)*DIFF;
real diff_num=  dt*(ddsdxx_num+ddsdyy_num+ddsdzz_num);
 
//s_d[ti + tj*blockDim.x]= diff_center+diff_num;
 s_d[ti + tj*blockDim.x]= ddsdxx*(DIFF+u_c*u_c*dt/2.f)+ddsdyy*(DIFF+v_c*v_c*dt/2.f)+ddsdzz*(DIFF+w_c*w_c*dt/2.f);


real sc_conv= s_c[ti + tj*blockDim.x];
real sc_diff= s_d[ti + tj*blockDim.x];

//TODO should we have s_f at t+1/2 step as well? 
//advance scalar. Take the particle volume fraction into consideration
//real rhs=(sc_diff+s_f[ti + tj*blockDim.x])/(1-p_epsp[ti + tj*blockDim.x]);
real rhs=sc_diff+s_f[ti + tj*blockDim.x];
//real rhs=s_f[ti + tj*blockDim.x];
sc_c[ti + tj*blockDim.x] +=(sc_conv+rhs)*dt;

//if(i==30&&j==30&&k==30) printf("\nsc_diff %f %f %f %f\n",sc_diff,sc_conv,s_f[ti + tj*blockDim.x],dt);
}

// copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
     && (i >= dom->Gcc._is && i < dom->Gcc._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      sc[C] = sc_c[ti + tj*blockDim.x];
      conv[C] = s_c[ti + tj*blockDim.x];
      diff[C] = s_d[ti + tj*blockDim.x];
    }
  }
}

//first order upwind
__global__ void advance_sc_upwind_1st(real DIFF, real *u, real *v, real *w,real *f,real *epsp,
  real *diff0, real *conv0, real *diff, real *conv, real *sc,real *sc0,
  dom_struct *dom, real dt0, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
//store scalar related values such as convective and diffusive term.
  __shared__ real s_d0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // diff0
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // conv0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // source term


//store scalar value of difference direction at cell center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top	 center
  __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at  	 center

//store point volume percentage in each grid cell
  __shared__ real p_epsp[MAX_THREADS_DIM * MAX_THREADS_DIM];       
 

  // working constants
  real ab0 = 0.5 * dt / dt0;  // for Adams-Bashforth stepping
  real ab = 1 + ab0;          // for Adams-Bashforth stepping

  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

  

  int C;

  // loop over z-planes
//  for(int k = dom->Gcc._ksb; k < dom->Gcc._keb; k++) {
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
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

    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];

s_c0[ti + tj*blockDim.x]=conv0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
s_d0[ti + tj*blockDim.x]=diff0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
 s_f[ti + tj*blockDim.x]=    f[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];

p_epsp[ti + tj*blockDim.x]= epsp[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
}

   // make sure all threads complete shared memory copy
    __syncthreads();
 
    if((j >= dom->Gfz._jsb && j < dom->Gfz._jeb)
      && (i >= dom->Gfz._isb && i < dom->Gfz._ieb)) {
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
   }
    if((j >= dom->Gfx._jsb && j < dom->Gfx._jeb)
      && (i >= dom->Gfx._isb && i < dom->Gfx._ieb)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._jsb && j < dom->Gfy._jeb)
      && (i >= dom->Gfy._isb && i < dom->Gfy._ieb)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
 
 
    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute convective term
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) 
 
//if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) 
{
//scalar on the west and east u face
real u_c=(s_u[ti + tj*blockDim.x]+s_u[(ti+1) + tj*blockDim.x])/2;// at xm(i)
real v_c=(s_v[ti + tj*blockDim.x]+s_v[ti + (tj+1)*blockDim.x])/2;
real w_c=(s_w1[ti + tj*blockDim.x]+s_w0[ti + tj*blockDim.x])/2;
 
real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0;//at x(i)
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;

//1st order upwind scheme
real dsudx,dsvdy,dswdz;
if(u_c>0)	dsudx= (sc_c[ti + tj*blockDim.x]*u_c-sc_uw*s_u[ti + tj*blockDim.x]) *2*ddx;//backward
else		dsudx=-(sc_c[ti + tj*blockDim.x]*u_c-sc_ue*s_u[ti+1 + tj*blockDim.x]) *2* ddx;//forward

if(v_c>0)	dsvdy= (sc_c[ti + tj*blockDim.x]*v_c-sc_vs*s_v[ti + tj*blockDim.x])*2*ddy;//backward
else		dsvdy=-(sc_c[ti + tj*blockDim.x]*v_c-sc_vn*s_v[ti + (tj+1)*blockDim.x])*2*ddy;//forward

if(w_c>0)	dswdz= (sc_c[ti + tj*blockDim.x]*w_c-sc_wb*s_w0[ti + tj*blockDim.x])*2*ddz;//backward
else		dswdz=-(sc_c[ti + tj*blockDim.x]*w_c-sc_wt*s_w1[ti + tj*blockDim.x])*2*ddz;//forward


//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;

     //move convective term to right hand side 
//real sc_conv=-(dsudx+dsvdy+dswdz);
      s_c[ti + tj*blockDim.x]=-(dsudx+dsvdy+dswdz);
     //Diffusion term
//real sc_diff=(ddsdxx+ddsdyy+ddsdzz)*DIFF;
      s_d[ti + tj*blockDim.x]= (ddsdxx+ddsdyy+ddsdzz)*DIFF;

//adam-bashforth at t+1/2 step
real sc_conv= ab * s_c[ti + tj*blockDim.x]-ab0*s_c0[ti + tj*blockDim.x];
real sc_diff= ab * s_d[ti + tj*blockDim.x]-ab0*s_d0[ti + tj*blockDim.x];

//TODO should we have s_f at t+1/2 step as well???
//right hand side of scalar equation
//s_rhs[ti + tj*blockDim.x]=sc_conv+sc_diff+s_f[ti + tj*blockDim.x];
//advance scalar
//Take the particle volume fraction into consideration
real rhs=(sc_diff+s_f[ti + tj*blockDim.x])/(1-p_epsp[ti + tj*blockDim.x]);
//real rhs=sc_diff+s_f[ti + tj*blockDim.x];
sc_c[ti + tj*blockDim.x] +=(sc_conv+rhs)*dt;
}

    // copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
     && (i >= dom->Gcc._is && i < dom->Gcc._ie)
  //  if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)
  //    && (i >= dom->Gcc._isb && i < dom->Gcc._ieb)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      sc[C] = sc_c[ti + tj*blockDim.x];
      conv[C] = s_c[ti + tj*blockDim.x];
      diff[C] = s_d[ti + tj*blockDim.x];
    }
  }
}


//first order upwind
__global__ void advance_sc_upwind_1st_init(real DIFF, real *u, real *v, real *w,real *f,real *epsp,
  real *diff0, real *conv0, real *diff, real *conv, real *sc,real *sc0,
  dom_struct *dom, real dt0, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
//store scalar related values such as convective and diffusive term.
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // source term


//store scalar value of difference direction at cell center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top	 center
  __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at  	 center

//store point volume percentage in each grid cell
  __shared__ real p_epsp[MAX_THREADS_DIM * MAX_THREADS_DIM];       
 

 
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

 
 

  int C;

  // loop over z-planes
//  for(int k = dom->Gcc._ksb; k < dom->Gcc._keb; k++) {
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
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

    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];

/*
s_c0[ti + tj*blockDim.x]=conv0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
s_d0[ti + tj*blockDim.x]=diff0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
*/
s_f[ti + tj*blockDim.x]=    f[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];

p_epsp[ti + tj*blockDim.x]= epsp[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
}

   // make sure all threads complete shared memory copy
    __syncthreads();
 
    if((j >= dom->Gfz._jsb && j < dom->Gfz._jeb)
      && (i >= dom->Gfz._isb && i < dom->Gfz._ieb)) {
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
   }
    if((j >= dom->Gfx._jsb && j < dom->Gfx._jeb)
      && (i >= dom->Gfx._isb && i < dom->Gfx._ieb)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._jsb && j < dom->Gfy._jeb)
      && (i >= dom->Gfy._isb && i < dom->Gfy._ieb)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
 
 
    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute convective term
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) 
 
//if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) 
{
//scalar on the west and east u face
real u_c=(s_u[ti + tj*blockDim.x]+s_u[(ti+1) + tj*blockDim.x])/2;// at xm(i)
real v_c=(s_v[ti + tj*blockDim.x]+s_v[ti + (tj+1)*blockDim.x])/2;
real w_c=(s_w1[ti + tj*blockDim.x]+s_w0[ti + tj*blockDim.x])/2;
 
real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0;//at x(i)
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;

//1st order upwind scheme
real dsudx,dsvdy,dswdz;
if(u_c>0)	dsudx= (sc_c[ti + tj*blockDim.x]*u_c-sc_uw*s_u[ti + tj*blockDim.x]) *2*ddx;//backward
else		dsudx=-(sc_c[ti + tj*blockDim.x]*u_c-sc_ue*s_u[ti+1 + tj*blockDim.x]) *2* ddx;//forward

if(v_c>0)	dsvdy= (sc_c[ti + tj*blockDim.x]*v_c-sc_vs*s_v[ti + tj*blockDim.x])*2*ddy;//backward
else		dsvdy=-(sc_c[ti + tj*blockDim.x]*v_c-sc_vn*s_v[ti + (tj+1)*blockDim.x])*2*ddy;//forward

if(w_c>0)	dswdz= (sc_c[ti + tj*blockDim.x]*w_c-sc_wb*s_w0[ti + tj*blockDim.x])*2*ddz;//backward
else		dswdz=-(sc_c[ti + tj*blockDim.x]*w_c-sc_wt*s_w1[ti + tj*blockDim.x])*2*ddz;//forward


//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;

     //move convective term to right hand side 
//real sc_conv=-(dsudx+dsvdy+dswdz);
      s_c[ti + tj*blockDim.x]=-(dsudx+dsvdy+dswdz);
     //Diffusion term
//real sc_diff=(ddsdxx+ddsdyy+ddsdzz)*DIFF;
      s_d[ti + tj*blockDim.x]= (ddsdxx+ddsdyy+ddsdzz)*DIFF;

//adam-bashforth at t+1/2 step
real sc_conv=  s_c[ti + tj*blockDim.x] ;
real sc_diff=  s_d[ti + tj*blockDim.x] ;

//TODO should we have s_f at t+1/2 step as well???
//right hand side of scalar equation
//s_rhs[ti + tj*blockDim.x]=sc_conv+sc_diff+s_f[ti + tj*blockDim.x];
//advance scalar
//Take the particle volume fraction into consideration
real rhs=(sc_diff+s_f[ti + tj*blockDim.x])/(1-p_epsp[ti + tj*blockDim.x]);
//real rhs=sc_diff+s_f[ti + tj*blockDim.x];
sc_c[ti + tj*blockDim.x] +=(sc_conv+rhs)*dt;
}

    // copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
     && (i >= dom->Gcc._is && i < dom->Gcc._ie)
  //  if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)
  //    && (i >= dom->Gcc._isb && i < dom->Gcc._ieb)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      sc[C] = sc_c[ti + tj*blockDim.x];
      conv[C] = s_c[ti + tj*blockDim.x];
      diff[C] = s_d[ti + tj*blockDim.x];
    }
  }
}


//Using central difference scheme in space, and adam-bashforth scheme in time.
__global__ void advance_sc(	real DIFF, 
				real *u, real *v, real *w,
				real *f,real *epsp,
				real *diff0, real *conv0, 
				real *diff, real *conv, 
				real *sc,real *sc0,
				dom_struct *dom, real dt0, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
//store scalar related values such as convective and diffusive term.
  __shared__ real s_d0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // diff0
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // conv0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // source term


//store scalar value of difference direction at cell center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top	 center
  __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at  	 center

//store point volume percentage in each grid cell
  __shared__ real p_epsp[MAX_THREADS_DIM * MAX_THREADS_DIM];       
 

  // working constants
  real ab0 = 0.5 * dt / dt0;  // for Adams-Bashforth stepping
  real ab = 1 + ab0;          // for Adams-Bashforth stepping

  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

 
  int C;

  // loop over z-planes
//  for(int k = dom->Gcc._ksb; k < dom->Gcc._keb; k++) {
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
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

    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];

s_c0[ti + tj*blockDim.x]=conv0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
s_d0[ti + tj*blockDim.x]=diff0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
 s_f[ti + tj*blockDim.x]=    f[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];

p_epsp[ti + tj*blockDim.x]= epsp[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];

//if(isnan(sc_c[ti + tj*blockDim.x])) printf("\nglobal_var %f %f %f %f %d %d %d %d %d\n",sc_b[ti + tj*blockDim.x],sc_t[ti + tj*blockDim.x],s_f[ti + tj*blockDim.x],sc_c[ti + tj*blockDim.x],ti,tj,i,j,k);

}

   // make sure all threads complete shared memory copy
    __syncthreads();

    if((j >= dom->Gfz._js && j < dom->Gfz._je)
      && (i >= dom->Gfz._is && i < dom->Gfz._ie)) {
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
    }
    if((j >= dom->Gfx._js && j < dom->Gfx._je)
      && (i >= dom->Gfx._is && i < dom->Gfx._ie)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._js && j < dom->Gfy._je)
      && (i >= dom->Gfy._is && i < dom->Gfy._ie)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }

  
    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute convective term
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) {
//scalar at face center
real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0f;
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0f;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0f;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0f;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;

//convective term of scalar at cell center
real dsudx=(s_u[(ti+1) + tj*blockDim.x]*sc_ue - s_u[ti + tj*blockDim.x]*sc_uw) * ddx;
real dsvdy=(s_v[ti + (tj+1)*blockDim.x]*sc_vn - s_v[ti + tj*blockDim.x]*sc_vs) * ddy;
real dswdz=(s_w1[ti + tj*blockDim.x]   *sc_wt - s_w0[ti + tj*blockDim.x]*sc_wb) * ddz;

//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;

     //move convective term to right hand side 
      s_c[ti + tj*blockDim.x]=-(dsudx+dsvdy+dswdz);
     //Diffusion term
      s_d[ti + tj*blockDim.x]= (ddsdxx+ddsdyy+ddsdzz)*DIFF;

//if(isnan(sc_c[ti + tj*blockDim.x])) printf("\nconv_diff %f %f %f %f %f %d %d %d %d %d\n",dsudx,ddsdxx,ddx,ab,sc_c[ti + tj*blockDim.x],ti,tj,i,j,k);
//Adam-bashforth at t+1/2 step
real sc_conv= ab * s_c[ti + tj*blockDim.x]-ab0*s_c0[ti + tj*blockDim.x];
real sc_diff= ab * s_d[ti + tj*blockDim.x]-ab0*s_d0[ti + tj*blockDim.x];

//if(isnan(s_f[ti + tj*blockDim.x])) printf("\nsource %f %f %f %f %d %d %d %d %d\n",s_d[ti + tj*blockDim.x],s_c[ti + tj*blockDim.x],s_f[ti + tj*blockDim.x],sc_c[ti + tj*blockDim.x],ti,tj,i,j,k);

//TODO should we have s_f at t+1/2 step as well? 
//advance scalar. Take the particle volume fraction into consideration
//real rhs=(sc_diff+s_f[ti + tj*blockDim.x])/(1-p_epsp[ti + tj*blockDim.x]);
real rhs=sc_diff+s_f[ti + tj*blockDim.x];
//real rhs=s_f[ti + tj*blockDim.x];
sc_c[ti + tj*blockDim.x] +=(sc_conv+rhs)*dt;

//if(i==36&&j==38&&k==35) printf("\nsc_diff %f %f %f %f\n",sc_diff,sc_conv,s_f[ti + tj*blockDim.x],sc_c[ti + tj*blockDim.x]);
}

// copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
     && (i >= dom->Gcc._is && i < dom->Gcc._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      sc[C] = sc_c[ti + tj*blockDim.x];
      conv[C] = s_c[ti + tj*blockDim.x];
      diff[C] = s_d[ti + tj*blockDim.x];
    }
  }
}




__global__ void advance_sc_init(real DIFF, real *u, real *v, real *w,real *f,real *epsp,
  real *diff0, real *conv0, real *diff, real *conv, real *sc,real *sc0,
  dom_struct *dom, real dt0, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
  
//store scalar related values such as convective and diffusive term.
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
//  __shared__ real s_d0[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
//  __shared__ real s_c0[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv0
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // source term


//store scalar value of difference direction at cell center
  __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at  	 center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top	 center

//store point volume percentage in each grid cell
  __shared__ real p_epsp[MAX_THREADS_DIM * MAX_THREADS_DIM];       

  // working constants
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

  int C;

  // loop over z-planes
//  for(int k = dom->Gcc._ksb; k < dom->Gcc._keb; k++) {
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
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

    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];

p_epsp[ti + tj*blockDim.x]= epsp[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];

s_f[ti + tj*blockDim.x]=    f[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
}

  // make sure all threads complete shared memory copy
    __syncthreads();

    if((j >= dom->Gfz._js && j < dom->Gfz._je)
      && (i >= dom->Gfz._is && i < dom->Gfz._ie)) {
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
    }
    if((j >= dom->Gfx._js && j < dom->Gfx._je)
      && (i >= dom->Gfx._is && i < dom->Gfx._ie)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._js && j < dom->Gfy._je)
      && (i >= dom->Gfy._is && i < dom->Gfy._ie)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }

/*
   if((j >= dom->Gfz._jsb && j < dom->Gfz._jeb)
      && (i >= dom->Gfz._isb && i < dom->Gfz._ieb)) {
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
    }
    if((j >= dom->Gfx._jsb && j < dom->Gfx._jeb)
      && (i >= dom->Gfx._isb && i < dom->Gfx._ieb)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._jsb && j < dom->Gfy._jeb)
      && (i >= dom->Gfy._isb && i < dom->Gfy._ieb)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
*/


/*
    //convective term
    s_c[ti + tj*blockDim.x] = 0.0;
    //diffusion term
    s_d[ti + tj*blockDim.x] = 0.0;
*/
   // make sure all threads complete shared memory copy
    __syncthreads();
 
    // compute convective term
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) {
//scalar on the west and east u face
real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0f;
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0f;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0f;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0f;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;

//convective term of scalar at cell center
real dsudx=(s_u[(ti+1) + tj*blockDim.x]*sc_ue - s_u[ti + tj*blockDim.x]*sc_uw) * ddx;
real dsvdy=(s_v[ti + (tj+1)*blockDim.x]*sc_vn - s_v[ti + tj*blockDim.x]*sc_vs) * ddy;
real dswdz=(s_w1[ti + tj*blockDim.x]   *sc_wt - s_w0[ti + tj*blockDim.x]*sc_wb) * ddz;


//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;
 
     //move convective term to right hand side 
      s_c[ti + tj*blockDim.x]=-(dsudx+dsvdy+dswdz);
     //Diffusion term
      s_d[ti + tj*blockDim.x]= (ddsdxx+ddsdyy+ddsdzz)*DIFF;

//advance scalar
//sc_c[ti + tj*blockDim.x] +=(s_c[ti + tj*blockDim.x]+s_d[ti + tj*blockDim.x]+s_f[ti + tj*blockDim.x])*dt;
//real rhs=(s_d[ti + tj*blockDim.x]+s_f[ti + tj*blockDim.x])/(1-p_epsp[ti + tj*blockDim.x]);
real rhs=s_d[ti + tj*blockDim.x]+s_f[ti + tj*blockDim.x];
//real rhs=s_f[ti + tj*blockDim.x];
sc_c[ti + tj*blockDim.x] +=(s_c[ti + tj*blockDim.x]+rhs)*dt;

//if(i==36&&j==38&&k==35) printf("\nsc_diff0 %f %f %f %d\n",s_d[ti + tj*blockDim.x],s_c[ti + tj*blockDim.x],s_f[ti + tj*blockDim.x],k);

}

    __syncthreads();

    // copy shared memory back to global
//sc stores n+1 timestep, conv and diff are n timestep
  if((j >= dom->Gcc._js && j < dom->Gcc._je)
      && (i >= dom->Gcc._is && i < dom->Gcc._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      sc[C] = sc_c[ti + tj*blockDim.x];
      conv[C] = s_c[ti + tj*blockDim.x];
      diff[C] = s_d[ti + tj*blockDim.x];
    }

  }
}


 __global__ void lpt_scalar_source_init(real *A, real *B,dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;
  int ti = blockDim.z*blockIdx.z + threadIdx.z;
/*
  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;
*/
int s1b=tex1Dfetch(texRefDomInfo,21);
int s2b=tex1Dfetch(texRefDomInfo,22);

int ieb=tex1Dfetch(texRefDomInfo,4);
int jeb=tex1Dfetch(texRefDomInfo,10);
int keb=tex1Dfetch(texRefDomInfo,16);


/*
if(tj < dom->Gcc._jeb && tk < dom->Gcc._keb) {
    for(int i = dom->Gcc._isb; i < dom->Gcc._ieb; i++) 
*/

if(ti<ieb &&tj < jeb && tk < keb) {
	{
	A[ti+tj*s1b+tk*s2b]=0.f;
	B[ti+tj*s1b+tk*s2b]=0.f;
	}
 }

}





__global__ void lpt_epsp_clip(real *epsp,dom_struct *dom)
{

  int tj = blockDim.x*blockIdx.x + threadIdx.x- 2*blockIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y- 2*blockIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;
if(tj < dom->Gcc._jeb && tk < dom->Gcc._keb) {
    for(int i = dom->Gcc._isb; i < dom->Gcc._ieb; i++) 
        {
        if(epsp[i+tj*s1b+tk*s2b]>EPSP_CLIP)  epsp[i+tj*s1b+tk*s2b]=EPSP_CLIP;
 //       if(epsp[i+tj*s1b+tk*s2b]>0) printf("\ntest %d %d %d %f\n",i,tj,tk,epsp[i+tj*s1b+tk*s2b]); 
        }
    }
}











//This test case only works with 0 convection!!
//Test scalar source is cos(wx*x)*cos(wy*y)*(cos(wt*t)+(wx*wx+wy*wy)*DIFF*sin(wt*t)/wt); 
//The corresponding solution is cos(wx*x)*cos(wy*y)*DIFF*sin(wt*t)/wt
__global__ void lpt_scalar_source_init_test(real *src,dom_struct *dom,real t, real DIFF)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x- 2*blockIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y- 2*blockIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;
if(tj < dom->Gcc._jeb && tk < dom->Gcc._keb) {
    for(int i = dom->Gcc._isb; i < dom->Gcc._ieb; i++) 
	{
 real x = ((i-0.5) * dom->dx) + dom->xs;
 real y = ((tj-0.5) * dom->dy) + dom->ys;
 real wt=10;
 real wx=1;
 real wy=2;
	src[i+tj*s1b+tk*s2b]=DIFF*cos(wx*x)*cos(wy*y)*(cos(wt*t)+(wx*wx+wy*wy)*DIFF*sin(wt*t)/wt);	
	}
 }

}

//Need to run and modify  forcing_test_taylor<<<numBlocks, dimBlocks>>>(_f_x[dev],_f_y[dev], _dom[dev],ttime,DIFF,_points[dev],n2,n3) in cuda_bluebottle.cu
//suppose u=U0*cos(qx)cos(qy)sin(rt); v=U0*sin(qx)sin(qy)sin(rt)
//Test scalar source is cos(wx*x)*cos(wy*y)*(cos(wt*t)+(wx*wx+wy*wy)*DIFF*sin(wt*t)/wt)-1.f/wt*U0*sin(r*t)*sin(wt*t)*(wx*cos(q*x)*cos(q*y)*cos(wy*y)*sin(wx*x)+wy*cos(wx*x)*sin(q*x)*sin(q*y)*sin(wy*y)); 
//Test solution for scalar is cos(wx*x)*cos(wy*y)*sin(wt*t)
__global__ void lpt_scalar_source_convDiff_test(real *src,dom_struct *dom,real t, real DIFF)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x- 2*blockIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y- 2*blockIdx.y;
  int i;
  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;
/*
 real wt=10*2*PI;
 real wx=6;
 real wy=8;

 real q=6;
 real r=wt;
 real U0=20;
*/

 real wt=10*2*PI;
 real wx=2;
 real wy=2;

 real q=2;
 real r=wt;
 real U0=10;
// real U0=0;



if(tj < dom->Gcc._jeb && tk < dom->Gcc._keb) {
    for( i = dom->Gcc._isb; i < dom->Gcc._ieb; i++) 
	{
//if(tj==dom->Gcc._je&&tk==dom->Gcc._ke&&i==dom->Gcc._ie) printf("\ncorrect block\n");

 real x = ((i-0.5) * dom->dx) + dom->xs;
 real y = ((tj-0.5) * dom->dy) + dom->ys;
	src[i+tj*s1b+tk*s2b]=cos(wx*x)*cos(wy*y)*(cos(wt*t)+(wx*wx+wy*wy)*DIFF*sin(wt*t)/wt)-1.f/wt*U0*sin(r*t)*sin(wt*t)*(wx*cos(q*x)*cos(q*y)*cos(wy*y)*sin(wx*x)+wy*cos(wx*x)*sin(q*x)*sin(q*y)*sin(wy*y));	
	}
 }

}


