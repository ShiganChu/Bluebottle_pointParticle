#include "cuda_scalar.h"
//#include "cuda_point.h"
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
real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0;
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;

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

//Adam-bashforth at t+1/2 step
real sc_conv= ab * s_c[ti + tj*blockDim.x]-ab0*s_c0[ti + tj*blockDim.x];
real sc_diff= ab * s_d[ti + tj*blockDim.x]-ab0*s_d0[ti + tj*blockDim.x];

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

/*
s_c0[ti + tj*blockDim.x]=conv0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
s_d0[ti + tj*blockDim.x]=diff0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
*/
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
real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0;
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;

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

//if(i==30&&j==30&&k>30&&k<31) printf("\nsc_diff0 %f %f %f %d\n",s_d[ti + tj*blockDim.x],s_c[ti + tj*blockDim.x],s_f[ti + tj*blockDim.x],k);

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
  int tj = blockDim.x*blockIdx.x + threadIdx.x- 2*blockIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y- 2*blockIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

if(tj < dom->Gcc._jeb && tk < dom->Gcc._keb) {
    for(int i = dom->Gcc._isb; i < dom->Gcc._ieb; i++) 
	{
	A[i+tj*s1b+tk*s2b]=0.f;
	B[i+tj*s1b+tk*s2b]=0.f;
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
        }
    }
}












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
	src[i+tj*s1b+tk*s2b]=cos(wx*x)*cos(wy*y)*(cos(wt*t)+(wx*wx+wy*wy)*DIFF*sin(wt*t)/wt);	
	}
 }

}


