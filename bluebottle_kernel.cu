#include "cuda_bluebottle.h"

// pressure; west; periodic
__global__ void BC_p_W_P(real *p, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    p[dom->Gcc._isb + tj*s1b + tk*s2b] = p[(dom->Gcc._ie-1) + tj*s1b + tk*s2b];
}

// pressure; west; Neumann
__global__ void BC_p_W_N(real *p, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    p[dom->Gcc._isb + tj*s1b + tk*s2b] = p[dom->Gcc._is + tj*s1b + tk*s2b];
}

// pressure; east; periodic
__global__ void BC_p_E_P(real *p, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    p[(dom->Gcc._ieb-1) + tj*s1b + tk*s2b] = p[dom->Gcc._is + tj*s1b + tk*s2b];
}

// pressure; east; Neumann
__global__ void BC_p_E_N(real *p, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((tj < dom->Gcc._jnb) && (tk < dom->Gcc._knb))
    p[(dom->Gcc._ieb-1) + tj*s1b + tk*s2b] = p[(dom->Gcc._ie-1)
      + tj*s1b + tk*s2b];
}

// pressure; south; periodic
__global__ void BC_p_S_P(real *p, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
    p[ti + dom->Gcc._jsb*s1b + tk*s2b] = p[ti + (dom->Gcc._je-1)*s1b + tk*s2b];
}

// pressure; south; Neumann
__global__ void BC_p_S_N(real *p, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
    p[ti + dom->Gcc._jsb*s1b + tk*s2b] = p[ti + dom->Gcc._js*s1b + tk*s2b];
}

// pressure; north; periodic
__global__ void BC_p_N_P(real *p, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
    p[ti + (dom->Gcc._jeb-1)*s1b + tk*s2b] = p[ti + dom->Gcc._js*s1b + tk*s2b];
}

// pressure; north; Neumann
__global__ void BC_p_N_N(real *p, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tk < dom->Gcc._knb))
    p[ti + (dom->Gcc._jeb-1)*s1b + tk*s2b] = p[ti
      + (dom->Gcc._je-1)*s1b + tk*s2b];
}

// pressure; bottom; periodic
__global__ void BC_p_B_P(real *p, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    p[ti + tj*s1b + dom->Gcc._ksb*s2b] = p[ti + tj*s1b + (dom->Gcc._ke-1)*s2b];
}

// pressure; bottom; Neumann
__global__ void BC_p_B_N(real *p, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b; 
  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    p[ti + tj*s1b + dom->Gcc._ksb*s2b] = p[ti + tj*s1b + dom->Gcc._ks*s2b];
}

// pressure; top; periodic
__global__ void BC_p_T_P(real *p, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    p[ti + tj*s1b + (dom->Gcc._keb-1)*s2b] = p[ti + tj*s1b + dom->Gcc._ks*s2b];
}

// pressure; top; Neumann
__global__ void BC_p_T_N(real *p, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gcc._s1b;
  int s2b = dom->Gcc._s2b;

  if((ti < dom->Gcc._inb) && (tj < dom->Gcc._jnb))
    p[ti + tj*s1b + (dom->Gcc._keb-1)*s2b] = p[ti
      + tj*s1b + (dom->Gcc._ke-1)*s2b];
}

// u-velocity; west; periodic
__global__ void BC_u_W_P(real *u, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    u[dom->Gfx._isb + tj*s1b + tk*s2b] = u[(dom->Gfx._ie-2) + tj*s1b + tk*s2b];
    //u[dom->Gfx._is + tj*s1b + tk*s2b] = u[(dom->Gfx._ie-1) + tj*s1b + tk*s2b];
  }
}

// u-velocity; west; Dirichlet
__global__ void BC_u_W_D(real *u, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    u[dom->Gfx._isb + tj*s1b + tk*s2b] = 2 * bc
      - u[(dom->Gfx._is+1) + tj*s1b + tk*s2b];
    u[dom->Gfx._is + tj*s1b + tk*s2b] = bc;
  }
}

// u-velocity; west; Neumann
__global__ void BC_u_W_N(real *u, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb))
    u[dom->Gfx._isb + tj*s1b + tk*s2b] = u[dom->Gfx._is + tj*s1b + tk*s2b];
}

// u-velocity; west; Turbulent precursor
__global__ void BC_u_W_T(real *u, dom_struct *dom, real* bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    u[dom->Gfx._isb + tj*s1b + tk*s2b] = 2 * bc[tj + tk*dom->Gfx.jnb]
      - u[(dom->Gfx._is+1) + tj*s1b + tk*s2b];
    u[dom->Gfx._is + tj*s1b + tk*s2b] = bc[tj + tk*dom->Gfx.jnb];
  }
}

// u-velocity; east; periodic
__global__ void BC_u_E_P(real *u, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    u[(dom->Gfx._ieb-1) + tj*s1b + tk*s2b] = u[(dom->Gfx._is+1) + tj*s1b + tk*s2b];
    //u[(dom->Gfx._ie-1) + tj*s1b + tk*s2b] = u[dom->Gfx._is + tj*s1b + tk*s2b];
  }
}

// u-velocity; east; Dirichlet
__global__ void BC_u_E_D(real *u, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    u[(dom->Gfx._ieb-1) + tj*s1b + tk*s2b] = 2 * bc - u[(dom->Gfx._ie-2)
      + tj*s1b + tk*s2b];
    u[(dom->Gfx._ie-1) + tj*s1b + tk*s2b] = bc;
  }
}

// u-velocity; east; Neumann
__global__ void BC_u_E_N(real *u, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb))
    u[(dom->Gfx._ieb-1) + tj*s1b + tk*s2b] = u[(dom->Gfx._ie-1)
      + tj*s1b + tk*s2b];
}

// u-velocity; east; Turbulent precursor
__global__ void BC_u_E_T(real *u, dom_struct *dom, real* bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    u[(dom->Gfx._ieb-1) + tj*s1b + tk*s2b] = 2 * bc[tj + tk*dom->Gfx.jnb]
      - u[(dom->Gfx._ie-2) + tj*s1b + tk*s2b];
    u[(dom->Gfx._ie-1) + tj*s1b + tk*s2b] = bc[tj + tk*dom->Gfx.jnb];
  }
}

// u-velocity; south; periodic
__global__ void BC_u_S_P(real *u, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb)) {
    u[ti + dom->Gfx._jsb*s1b + tk*s2b] = u[ti + (dom->Gfx._je-1)*s1b + tk*s2b];
  }
}

// u-velocity; south; Dirichlet
__global__ void BC_u_S_D(real *u, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb))
    u[ti + dom->Gfx._jsb*s1b + tk*s2b] = 2 * bc
      - u[ti + dom->Gfx._js*s1b + tk*s2b];
}

// u-velocity; south; Neumann
__global__ void BC_u_S_N(real *u, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb))
    u[ti + dom->Gfx._jsb*s1b + tk*s2b] = u[ti + dom->Gfx._js*s1b + tk*s2b];
}

// u-velocity; south; Turbulent precursor
__global__ void BC_u_S_T(real *u, dom_struct *dom, real* bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb)) {
    u[ti + dom->Gfx._jsb*s1b + tk*s2b] = 2 * bc[tk + ti*dom->Gfx.knb]
      - u[ti + dom->Gfx._js*s1b + tk*s2b];
  }
}

// u-velocity; north; periodic
__global__ void BC_u_N_P(real *u, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb))
    u[ti + (dom->Gfx._jeb-1)*s1b + tk*s2b] = u[ti
      + dom->Gfx._js*s1b + tk*s2b];
}

// u-velocity; north; Dirichlet
__global__ void BC_u_N_D(real *u, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb))
    u[ti + (dom->Gfx._jeb-1)*s1b + tk*s2b] = 2 * bc - u[ti
      + (dom->Gfx._je-1)*s1b + tk*s2b];
}

// u-velocity; north; Neumann
__global__ void BC_u_N_N(real *u, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb))
    u[ti + (dom->Gfx._jeb-1)*s1b + tk*s2b] = u[ti
      + (dom->Gfx._je-1)*s1b + tk*s2b];
}

// u-velocity; north; Turbulent precursor
__global__ void BC_u_N_T(real *u, dom_struct *dom, real* bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tk < dom->Gfx._knb)) {
    u[ti + (dom->Gfx._jeb-1)*s1b + tk*s2b] = 2 * bc[tk + ti*dom->Gfx.knb]
      - u[ti + (dom->Gfx._je-1)*s1b + tk*s2b];
  }
}

// u-velocity; bottom; periodic
__global__ void BC_u_B_P(real *u, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb))
    u[ti + tj*s1b + dom->Gfx._ksb*s2b] = u[ti + tj*s1b + (dom->Gfx._ke-1)*s2b];
}

// u-velocity; bottom; Dirichlet
__global__ void BC_u_B_D(real *u, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb))
    u[ti + tj*s1b + dom->Gfx._ksb*s2b] = 2 * bc
      - u[ti + tj*s1b + dom->Gfx._ks*s2b];
}

// u-velocity; bottom; Neumann
__global__ void BC_u_B_N(real *u, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb))
    u[ti + tj*s1b + dom->Gfx._ksb*s2b] = u[ti + tj*s1b + dom->Gfx._ks*s2b];
}

// u-velocity; bottom; Turbulent precursor
__global__ void BC_u_B_T(real *u, dom_struct *dom, real* bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb)) {
    u[ti + tj*s1b + dom->Gfx._ksb*s2b] = 2 * bc[ti + tj*dom->Gfx.inb]
      - u[ti + tj*s1b + dom->Gfx._ks*s2b];
  }
}

// u-velocity; top; periodic
__global__ void BC_u_T_P(real *u, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb))
    u[ti + tj*s1b + (dom->Gfx._keb-1)*s2b] = u[ti
      + tj*s1b + dom->Gfx._ks*s2b];
}

// u-velocity; top; Dirichlet
__global__ void BC_u_T_D(real *u, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb))
    u[ti + tj*s1b + (dom->Gfx._keb-1)*s2b] = 2 * bc - u[ti + tj*s1b +
      (dom->Gfx._ke-1)*s2b];
}

// u-velocity; top; Neumann
__global__ void BC_u_T_N(real *u, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb))
    u[ti + tj*s1b + (dom->Gfx._keb-1)*s2b] = u[ti
      + tj*s1b + (dom->Gfx._ke-1)*s2b];
}

// u-velocity; top; Turbulent precursor
__global__ void BC_u_T_T(real *u, dom_struct *dom, real* bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfx._s1b;
  int s2b = dom->Gfx._s2b;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb)) {
    u[ti + tj*s1b + (dom->Gfx._keb-1)*s2b] = 2 * bc[ti + tj*dom->Gfx.inb]
      - u[ti + tj*s1b + (dom->Gfx._ke-1)*s2b];
  }
}

// v-velocity; west; periodic
__global__ void BC_v_W_P(real *v, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb))
    v[dom->Gfy._isb + tj*s1b + tk*s2b] = v[(dom->Gfy._ie-1) + tj*s1b + tk*s2b];
}

// v-velocity; west; Dirichlet
__global__ void BC_v_W_D(real *v, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb))
    v[dom->Gfy._isb + tj*s1b + tk*s2b] = 2 * bc
      - v[dom->Gfy._is + tj*s1b + tk*s2b];
}

// v-velocity; west; Neumann
__global__ void BC_v_W_N(real *v, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb))
    v[dom->Gfy._isb + tj*s1b + tk*s2b] = v[dom->Gfy._is + tj*s1b + tk*s2b];
}

// v-velocity; west; Turbulent precursor
__global__ void BC_v_W_T(real *v, dom_struct *dom, real* bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb)) {
    v[dom->Gfy._isb + tj*s1b + tk*s2b] = 2 * bc[tj + tk*dom->Gfy.jnb]
      - v[(dom->Gfy._is) + tj*s1b + tk*s2b];
  }
}

// v-velocity; east; periodic
__global__ void BC_v_E_P(real *v, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb))
    v[(dom->Gfy._ieb-1) + tj*s1b + tk*s2b] = v[dom->Gfy._is
      + tj*s1b + tk*s2b];
}

// v-velocity; east; Dirichlet
__global__ void BC_v_E_D(real *v, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb))
    v[(dom->Gfy._ieb-1) + tj*s1b + tk*s2b] = 2 * bc - v[(dom->Gfy._ie-1)
      + tj*s1b + tk*s2b];
}

// v-velocity; east; Neumann
__global__ void BC_v_E_N(real *v, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb))
    v[(dom->Gfy._ieb-1) + tj*s1b + tk*s2b] = v[(dom->Gfy._ie-1)
      + tj*s1b + tk*s2b];
}

// v-velocity; east; Turbulent precursor
__global__ void BC_v_E_T(real *v, dom_struct *dom, real* bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb)) {
    v[(dom->Gfy._ieb-1) + tj*s1b + tk*s2b] = 2 * bc[tj + tk*dom->Gfy.jnb]
      - v[(dom->Gfy._ie-1) + tj*s1b + tk*s2b];
  }
}

// v-velocity; south; periodic
__global__ void BC_v_S_P(real *v, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    v[ti + dom->Gfy._jsb*s1b + tk*s2b] = v[ti + (dom->Gfy._je-2)*s1b + tk*s2b];
    //v[ti + dom->Gfy._js*s1b + tk*s2b] = v[ti + (dom->Gfy._je-1)*s1b + tk*s2b];
  }
}

// v-velocity; south; Dirichlet
__global__ void BC_v_S_D(real *v, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    v[ti + dom->Gfy._jsb*s1b + tk*s2b] = 2 * bc - v[ti
      + (dom->Gfy._js+1)*s1b + tk*s2b];
    v[ti + dom->Gfy._js*s1b + tk*s2b] = bc;
  }
}

// v-velocity; south; Neumann
__global__ void BC_v_S_N(real *v, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb))
    v[ti + dom->Gfy._jsb*s1b + tk*s2b] = v[ti + dom->Gfy._js*s1b + tk*s2b];
}

// v-velocity; south; Turbulent precursor
__global__ void BC_v_S_T(real *v, dom_struct *dom, real* bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    v[ti + dom->Gfy._jsb*s1b + tk*s2b] = 2 * bc[tk + ti*dom->Gfy.knb]
      - v[ti + (dom->Gfy._js+1)*s1b + tk*s2b];
    v[ti + dom->Gfy._js*s1b + tk*s2b] = bc[tk + ti*dom->Gfy.knb];
  }
}

// v-velocity; north; periodic
__global__ void BC_v_N_P(real *v, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    v[ti + (dom->Gfy._jeb-1)*s1b + tk*s2b] = v[ti
      + (dom->Gfy._js+1)*s1b + tk*s2b];
    //v[ti + (dom->Gfy._je-1)*s1b + tk*s2b] = v[ti
      //+ dom->Gfy._js*s1b + tk*s2b];
  }
}

// v-velocity; north; Dirichlet
__global__ void BC_v_N_D(real *v, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    v[ti + (dom->Gfy._jeb-1)*s1b + tk*s2b] = 2 * bc - v[ti +
      (dom->Gfy._je-2)*s1b + tk*s2b];
    v[ti + (dom->Gfy._je-1)*s1b + tk*s2b] = bc;
  }
}

// v-velocity; north; Neumann
__global__ void BC_v_N_N(real *v, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb))
    v[ti + (dom->Gfy._jeb-1)*s1b + tk*s2b] = v[ti
      + (dom->Gfy._je-1)*s1b + tk*s2b];
}

// v-velocity; north; Turbulent precursor
__global__ void BC_v_N_T(real *v, dom_struct *dom, real* bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    v[ti + (dom->Gfy._jeb-1)*s1b + tk*s2b] = 2 * bc[tk + ti*dom->Gfy.knb]
      - v[ti + (dom->Gfy._je-2)*s1b + tk*s2b];
    v[ti + (dom->Gfy._je-1)*s1b + tk*s2b] = bc[tk + ti*dom->Gfy.knb];
  }
}

// v-velocity; bottom; periodic
__global__ void BC_v_B_P(real *v, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb))
    v[ti + tj*s1b + dom->Gfy._ksb*s2b] = v[ti
      + tj*s1b + (dom->Gfy._ke-1)*s2b];
}

// v-velocity; bottom; Dirichlet
__global__ void BC_v_B_D(real *v, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb))
    v[ti + tj*s1b + dom->Gfy._ksb*s2b] = 2 * bc
      - v[ti + tj*s1b + dom->Gfy._ks*s2b];
}

// v-velocity; bottom; Neumann
__global__ void BC_v_B_N(real *v, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb))
    v[ti + tj*s1b + dom->Gfy._ksb*s2b] = v[ti + tj*s1b + dom->Gfy._ks*s2b];
}

// v-velocity; bottom; Turbulent precursor
__global__ void BC_v_B_T(real *v, dom_struct *dom, real* bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb)) {
    v[ti + tj*s1b + dom->Gfy._ksb*s2b] = 2 * bc[ti + tj*dom->Gfy.inb]
      - v[ti + tj*s1b + dom->Gfy._ks*s2b];
  }
}

// v-velocity; top; periodic
__global__ void BC_v_T_P(real *v, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb))
    v[ti + tj*s1b + (dom->Gfy._keb-1)*s2b] = v[ti
      + tj*s1b + dom->Gfy._ks*s2b];
}

// v-velocity; top; Dirichlet
__global__ void BC_v_T_D(real *v, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb))
    v[ti + tj*s1b + (dom->Gfy._keb-1)*s2b] = 2 * bc - v[ti + tj*s1b +
      (dom->Gfy._ke-1)*s2b];
}

// v-velocity; top; Neumann
__global__ void BC_v_T_N(real *v, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb))
    v[ti + tj*s1b + (dom->Gfy._keb-1)*s2b] = v[ti
      + tj*s1b + (dom->Gfy._ke-1)*s2b];
}

// v-velocity; top; Turbulent precursor
__global__ void BC_v_T_T(real *v, dom_struct *dom, real* bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfy._s1b;
  int s2b = dom->Gfy._s2b;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb)) {
    v[ti + tj*s1b + (dom->Gfy._keb-1)*s2b] = 2 * bc[ti + tj*dom->Gfy.inb]
      - v[ti + tj*s1b + (dom->Gfy._ke-1)*s2b];
  }
}

// w-velocity; west; periodic
__global__ void BC_w_W_P(real *w, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb))
    w[dom->Gfz._isb + tj*s1b + tk*s2b] = w[(dom->Gfz._ie-1) + tj*s1b + tk*s2b];
}

// w-velocity; west; Dirichlet
__global__ void BC_w_W_D(real *w, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb))
    w[dom->Gfz._isb + tj*s1b + tk*s2b] = 2 * bc
      - w[dom->Gfz._is + tj*s1b + tk*s2b];
}

// w-velocity; west; Neumann
__global__ void BC_w_W_N(real *w, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb))
    w[dom->Gfz._isb + tj*s1b + tk*s2b] = w[dom->Gfz._is + tj*s1b + tk*s2b];
}

// w-velocity; west; Turbulent precursor
__global__ void BC_w_W_T(real *w, dom_struct *dom, real* bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb)) {
    w[dom->Gfz._isb + tj*s1b + tk*s2b] = 2 * bc[tj + tk*dom->Gfz.jnb]
      - w[(dom->Gfz._is) + tj*s1b + tk*s2b];
  }
}

// w-velocity; east; periodic
__global__ void BC_w_E_P(real *w, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb))
    w[(dom->Gfz._ieb-1) + tj*s1b + tk*s2b] = w[dom->Gfz._is
      + tj*s1b + tk*s2b];
}

// w-velocity; east; Dirichlet
__global__ void BC_w_E_D(real *w, dom_struct *dom, real bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb))
    w[(dom->Gfz._ieb-1) + tj*s1b + tk*s2b] = 2 * bc - w[(dom->Gfz._ie-1)
      + tj*s1b + tk*s2b];
}

// w-velocity; east; Neumann
__global__ void BC_w_E_N(real *w, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb))
    w[(dom->Gfz._ieb-1) + tj*s1b + tk*s2b] = w[(dom->Gfz._ie-1)
      + tj*s1b + tk*s2b];
}

// w-velocity; east; Turbulent precursor
__global__ void BC_w_E_T(real *w, dom_struct *dom, real* bc)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb)) {
    w[(dom->Gfz._ieb-1) + tj*s1b + tk*s2b] = 2 * bc[tj + tk*dom->Gfz.jnb]
      - w[(dom->Gfz._ie-1) + tj*s1b + tk*s2b];
  }
}

// w-velocity; south; periodic
__global__ void BC_w_S_P(real *w, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb)) {
    w[ti + dom->Gfz._jsb*s1b + tk*s2b] = w[ti + (dom->Gfz._je-1)*s1b + tk*s2b];
  }
}

// w-velocity; south; Dirichlet
__global__ void BC_w_S_D(real *w, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb))
    w[ti + dom->Gfz._jsb*s1b + tk*s2b] = 2 * bc
      - w[ti + dom->Gfz._js*s1b + tk*s2b];
}

// w-velocity; south; Neumann
__global__ void BC_w_S_N(real *w, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb))
    w[ti + dom->Gfz._jsb*s1b + tk*s2b] = w[ti + dom->Gfz._js*s1b + tk*s2b];
}

// w-velocity; south; Turbulent precursor
__global__ void BC_w_S_T(real *w, dom_struct *dom, real* bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb)) {
    w[ti + dom->Gfz._jsb*s1b + tk*s2b] = 2 * bc[tk + ti*dom->Gfz.knb]
      - w[ti + dom->Gfz._js*s1b + tk*s2b];
  }
}

// w-velocity; north; periodic
__global__ void BC_w_N_P(real *w, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb))
    w[ti + (dom->Gfz._jeb-1)*s1b + tk*s2b] = w[ti
      + dom->Gfz._js*s1b + tk*s2b];
}

// w-velocity; north; Dirichlet
__global__ void BC_w_N_D(real *w, dom_struct *dom, real bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb))
    w[ti + (dom->Gfz._jeb-1)*s1b + tk*s2b] = 2 * bc - w[ti +
      (dom->Gfz._je-1)*s1b + tk*s2b];
}

// w-velocity; north; Neumann
__global__ void BC_w_N_N(real *w, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb))
    w[ti + (dom->Gfz._jeb-1)*s1b + tk*s2b] = w[ti
      + (dom->Gfz._je-1)*s1b + tk*s2b];
}

// w-velocity; north; Turbulent precursor
__global__ void BC_w_N_T(real *w, dom_struct *dom, real* bc)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb)) {
    w[ti + (dom->Gfz._jeb-1)*s1b + tk*s2b] = 2 * bc[tk + ti*dom->Gfz.knb]
      - w[ti + (dom->Gfz._je-1)*s1b + tk*s2b];
  }
}

// w-velocity; bottom; periodic
__global__ void BC_w_B_P(real *w, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    w[ti + tj*s1b + dom->Gfz._ksb*s2b] = w[ti + tj*s1b + (dom->Gfz._ke-2)*s2b];
    //w[ti + tj*s1b + dom->Gfz._ks*s2b] = w[ti + tj*s1b + (dom->Gfz._ke-1)*s2b];
  }
}

// w-velocity; bottom; Dirichlet
__global__ void BC_w_B_D(real *w, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    w[ti + tj*s1b + dom->Gfz._ksb*s2b] = 2 * bc - w[ti
      + tj*s1b + (dom->Gfz._ks+1)*s2b];
    w[ti + tj*s1b + dom->Gfz._ks*s2b] = bc;
  }
}

// w-velocity; bottom; Neumann
__global__ void BC_w_B_N(real *w, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb))
    w[ti + tj*s1b + dom->Gfz._ksb*s2b] = w[ti + tj*s1b + dom->Gfz._ks*s2b];
}

// w-velocity; bottom; Turbulent precursor
__global__ void BC_w_B_T(real *w, dom_struct *dom, real* bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    w[ti + tj*s1b + dom->Gfz._ksb*s2b] = 2 * bc[ti + tj*dom->Gfz.inb]
      - w[ti + tj*s1b + (dom->Gfz._ks+1)*s2b];
    w[ti + tj*s1b + dom->Gfz._ks*s2b] = bc[ti + tj*dom->Gfz.inb];
  }
}

// w-velocity; top; periodic
__global__ void BC_w_T_P(real *w, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    w[ti + tj*s1b + (dom->Gfz._keb-1)*s2b] = w[ti
      + tj*s1b + (dom->Gfz._ks+1)*s2b];
    //w[ti + tj*s1b + (dom->Gfz._ke-1)*s2b] = w[ti
      //+ tj*s1b + dom->Gfz._ks*s2b];
  }
}

// w-velocity; top; Dirichlet
__global__ void BC_w_T_D(real *w, dom_struct *dom, real bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    w[ti + tj*s1b + (dom->Gfz._keb-1)*s2b] = 2 * bc - w[ti + tj*s1b +
      (dom->Gfz._ke-2)*s2b];
    w[ti + tj*s1b + (dom->Gfz._ke-1)*s2b] = bc;
  }
}

// w-velocity; top; Neumann
__global__ void BC_w_T_N(real *w, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb))
    w[ti + tj*s1b + (dom->Gfz._keb-1)*s2b] = w[ti
      + tj*s1b + (dom->Gfz._ke-1)*s2b];
}

// w-velocity; top; Turbulent precursor
__global__ void BC_w_T_T(real *w, dom_struct *dom, real* bc)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  int s1b = dom->Gfz._s1b;
  int s2b = dom->Gfz._s2b;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    w[ti + tj*s1b + (dom->Gfz._keb-1)*s2b] = 2 * bc[ti + tj*dom->Gfz.inb]
      - w[ti + tj*s1b + (dom->Gfz._ke-2)*s2b];
    w[ti + tj*s1b + (dom->Gfz._ke-1)*s2b] = bc[ti + tj*dom->Gfz.inb];
  }
}

__global__ void project_u(real *u_star, real *p, real rho_f, real dt,
  real *u, dom_struct *dom, real ddx, int *flag_u)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tk = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if(tj < dom->Gfx._je && tk < dom->Gfx._ke) {
    for(int i = dom->Gfx._is; i < dom->Gfx._ie; i++) {
      real gradPhi = abs(flag_u[i + tj*dom->Gfx._s1b
        + tk*dom->Gfx._s2b])
        * ddx * (p[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p[(i-1) + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]);
      u[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = u_star[i + tj*dom->Gfx._s1b
        + tk*dom->Gfx._s2b] - dt / rho_f * gradPhi;
//if(i==16&&tj==8&&tk==8) printf("u_star grad_P %f\n",gradPhi);
    }
  }
}

__global__ void project_v(real *v_star, real *p, real rho_f, real dt,
  real *v, dom_struct *dom, real ddy, int *flag_v)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int ti = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if(tk < dom->Gfy._ke && ti < dom->Gfy._ie) {
    for(int j = dom->Gfy._js; j < dom->Gfy._je; j++) {
      real gradPhi = abs(flag_v[ti + j*dom->Gfy._s1b
        + tk*dom->Gfy._s2b])
        * ddy * (p[ti + j*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p[ti + (j-1)*dom->Gcc._s1b + tk*dom->Gcc._s2b]);
      v[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b] = v_star[ti + j*dom->Gfy._s1b
        + tk*dom->Gfy._s2b] - dt / rho_f * gradPhi;
//if(ti==16&j==8&tk==8) printf("gradPhi %f\n",gradPhi);

  }
  }
}

__global__ void project_w(real *w_star, real *p, real rho_f, real dt,
  real *w, dom_struct *dom, real ddz, int *flag_w)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x + DOM_BUF;
  int tj = blockIdx.y * blockDim.y + threadIdx.y + DOM_BUF;

  if(ti < dom->Gfz._ie && tj < dom->Gfz._je) {
    for(int k = dom->Gfz._ks; k < dom->Gfz._ke; k++) {
      real gradPhi = abs(flag_w[ti + tj*dom->Gfz._s1b
        + k*dom->Gfz._s2b])
        * ddz * (p[ti + tj*dom->Gcc._s1b + k*dom->Gcc._s2b]
        - p[ti + tj*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b]);
      w[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b] = w_star[ti + tj*dom->Gfz._s1b
        + k*dom->Gfz._s2b] - dt / rho_f * gradPhi;

//if(fabs(dt / rho_f * gradPhi)>0.1f) printf("\nproject id %d %d %d %f %f %f\n",ti,tj,k,w_star[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b],gradPhi, dt / rho_f * gradPhi);

    }
  }
}

__global__ void copy_p_ghost(real *p, real *p_tmp, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gcc.je-DOM_BUF && tk < dom->Gcc.ke-DOM_BUF) {
    for(int i = dom->Gcc.is-DOM_BUF; i < dom->Gcc.ie-DOM_BUF; i++) {
      p[(i+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc.s1b
        + (tk+DOM_BUF)*dom->Gcc.s2b] = p_tmp[i + tj*dom->Gcc.s1
        + tk*dom->Gcc.s2];
    }
  }
}

__global__ void copy_p_noghost(real *p_noghost, real *p_ghost, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gcc.je-DOM_BUF && tk < dom->Gcc.ke-DOM_BUF) {
    for(int i = dom->Gcc.is-DOM_BUF; i < dom->Gcc.ie-DOM_BUF; i++) {
      p_noghost[i + tj*dom->Gcc._s1 + tk*dom->Gcc._s2]
        = p_ghost[(i+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc._s1b
        + (tk+DOM_BUF)*dom->Gcc._s2b];
    }
  }
}

__global__ void copy_u_noghost(real *u_noghost, real *u_ghost, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gfx.je-DOM_BUF && tk < dom->Gfx.ke-DOM_BUF) {
    for(int i = dom->Gfx.is-DOM_BUF; i < dom->Gfx.ie-DOM_BUF; i++) {
      u_noghost[i + tj*dom->Gfx._s1 + tk*dom->Gfx._s2]
        = u_ghost[(i+DOM_BUF) + (tj+DOM_BUF)*dom->Gfx._s1b
        + (tk+DOM_BUF)*dom->Gfx._s2b];
    }
  }
}

__global__ void copy_v_noghost(real *v_noghost, real *v_ghost, dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if(tk < dom->Gfy.ke-DOM_BUF && ti < dom->Gfy.ie-DOM_BUF) {
    for(int j = dom->Gfy.js-DOM_BUF; j < dom->Gfy.je-DOM_BUF; j++) {
      v_noghost[ti + j*dom->Gfy._s1 + tk*dom->Gfy._s2]
        = v_ghost[(ti+DOM_BUF) + (j+DOM_BUF)*dom->Gfy._s1b
        + (tk+DOM_BUF)*dom->Gfy._s2b];
    }
  }
}

__global__ void copy_w_noghost(real *w_noghost, real *w_ghost, dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if(ti < dom->Gfz.ie-DOM_BUF && tj < dom->Gfz.je-DOM_BUF) {
    for(int k = dom->Gfz.ks-DOM_BUF; k < dom->Gfz.ke-DOM_BUF; k++) {
      w_noghost[ti + tj*dom->Gfz._s1 + k*dom->Gfz._s2]
        = w_ghost[(ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz._s1b
        + (k+DOM_BUF)*dom->Gfz._s2b];
    }
  }
}

__global__ void copy_u_fluid(real *u_noghost, real *u_ghost, int *phase, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gfx.je-DOM_BUF && tk < dom->Gfx.ke-DOM_BUF) {
    for(int i = dom->Gfx.is-DOM_BUF; i < dom->Gfx.ie-DOM_BUF; i++) {
      int boo = 1;
      if(phase[(i+DOM_BUF-1) + (tj+DOM_BUF)*dom->Gcc._s1b + (tk+DOM_BUF)*dom->Gcc._s2b] > -1) boo = 0;
      else if(phase[(i+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc._s1b + (tk+DOM_BUF)*dom->Gcc._s2b] > -1) boo = 0;
      u_noghost[i + tj*dom->Gfx._s1 + tk*dom->Gfx._s2]
        = boo * u_ghost[(i+DOM_BUF) + (tj+DOM_BUF)*dom->Gfx._s1b
        + (tk+DOM_BUF)*dom->Gfx._s2b];
    }
  }
}

__global__ void copy_v_fluid(real *v_noghost, real *v_ghost, int *phase, dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if(tk < dom->Gfy.ke-DOM_BUF && ti < dom->Gfy.ie-DOM_BUF) {
    for(int j = dom->Gfy.js-DOM_BUF; j < dom->Gfy.je-DOM_BUF; j++) {
      int boo = 1;
      if(phase[(ti+DOM_BUF) + (j+DOM_BUF-1)*dom->Gcc._s1b + (tk+DOM_BUF)*dom->Gcc._s2b] > -1) boo = 0;
      else if(phase[(ti+DOM_BUF) + (j+DOM_BUF)*dom->Gcc._s1b + (tk+DOM_BUF)*dom->Gcc._s2b] > -1) boo = 0;
      v_noghost[ti + j*dom->Gfy._s1 + tk*dom->Gfy._s2]
        = boo * v_ghost[(ti+DOM_BUF) + (j+DOM_BUF)*dom->Gfy._s1b
        + (tk+DOM_BUF)*dom->Gfy._s2b];
    }
  }
}

__global__ void copy_w_fluid(real *w_noghost, real *w_ghost, int *phase, dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if(ti < dom->Gfz.ie-DOM_BUF && tj < dom->Gfz.je-DOM_BUF) {
    for(int k = dom->Gfz.ks-DOM_BUF; k < dom->Gfz.ke-DOM_BUF; k++) {
      int boo = 1;
      if(phase[(ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc._s1b + (k+DOM_BUF-1)*dom->Gcc._s2b] > -1) boo = 0;
      else if(phase[(ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gcc._s1b + (k+DOM_BUF)*dom->Gcc._s2b] > -1) boo = 0;
      w_noghost[ti + tj*dom->Gfz._s1 + k*dom->Gfz._s2]
        = boo * w_ghost[(ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz._s1b
        + (k+DOM_BUF)*dom->Gfz._s2b];
    }
  }
}

__global__ void u_star_2(real rho_f, real nu,
  real *u0, real *v0, real *w0, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *u_star,
  dom_struct *dom, real dt0, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_u0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u back
  __shared__ real s_u1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u center
  __shared__ real s_u2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u forward
  __shared__ real s_v01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v back
  __shared__ real s_v12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v forward
  __shared__ real s_w01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // w back
  __shared__ real s_w12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // w forward
  __shared__ real s_d0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // diff0
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // conv0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_u_star[MAX_THREADS_DIM * MAX_THREADS_DIM];  // solution
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // x-force

  // working constants
  real ab0 = 0.5 * dt / dt0;  // for Adams-Bashforth stepping
  real ab = 1 + ab0;          // for Adams-Bashforth stepping
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
      s_d0[tj + tk*blockDim.x] = diff0[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_c0[tj + tk*blockDim.x] = conv0[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u0[tj + tk*blockDim.x] = u0[(i-1) + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u1[tj + tk*blockDim.x] = u0[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u2[tj + tk*blockDim.x] = u0[(i+1) + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_f[tj + tk*blockDim.x] = f[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((k >= dom->Gfy._ksb && k < dom->Gfy._keb)
      && (j >= dom->Gfy._jsb && j < dom->Gfy._jeb)) {
      s_v01[tj + tk*blockDim.x] = v0[(i-1) + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v12[tj + tk*blockDim.x] = v0[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
    if((k >= dom->Gfz._ksb && k < dom->Gfz._keb)
      && (j >= dom->Gfz._jsb && j < dom->Gfz._jeb)) {
      s_w01[tj + tk*blockDim.x] = w0[(i-1) + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w12[tj + tk*blockDim.x] = w0[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
    }

    s_u_star[tj + tk*blockDim.x] = 0.0;

    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute right-hand side
    // if off the shared memory block boundary
    if((tj > 0 && tj < blockDim.x-1) && (tk > 0 && tk < blockDim.y-1)) {
      real u011 = s_u0[tj + tk*blockDim.x];
      real u111 = s_u1[tj + tk*blockDim.x];
      real u211 = s_u2[tj + tk*blockDim.x];

      real u101 = s_u1[(tj-1) + tk*blockDim.x];
      real u121 = s_u1[(tj+1) + tk*blockDim.x];
      real v011 = s_v01[tj + tk*blockDim.x];
      real v111 = s_v12[tj + tk*blockDim.x];
      real v021 = s_v01[(tj+1) + tk*blockDim.x];
      real v121 = s_v12[(tj+1) + tk*blockDim.x];

      real u110 = s_u1[tj + (tk-1)*blockDim.x];
      real u112 = s_u1[tj + (tk+1)*blockDim.x];
      real w011 = s_w01[tj + tk*blockDim.x];
      real w111 = s_w12[tj + tk*blockDim.x];
      real w012 = s_w01[tj + (tk+1)*blockDim.x];
      real w112 = s_w12[tj + (tk+1)*blockDim.x];

      // compute convection term (Adams-Bashforth stepping)
      real duudx = (u211 + u111)*(u211 + u111) - (u111 + u011)*(u111 + u011);
      duudx *= 0.25 * ddx;

      real duvdy = (u121 + u111)*(v121 + v021) - (u111 + u101)*(v111 + v011);
      duvdy *= 0.25 * ddy;

      real duwdz = (u112 + u111)*(w112 + w012) - (u111 + u110)*(w111 + w011);
      duwdz *= 0.25 * ddz;

      s_c[tj + tk*blockDim.x] = duudx + duvdy + duwdz;

      // convection term sums into right-hand side
      s_u_star[tj + tk*blockDim.x] = - ab * s_c[tj + tk*blockDim.x]
        + ab0 * s_c0[tj + tk*blockDim.x];

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

      // diffusive term sums into right-hand side
      s_u_star[tj + tk*blockDim.x] += ab * s_d[tj + tk*blockDim.x]
        - ab0 * s_d0[tj + tk*blockDim.x];

      // add on imposed pressure gradient
      s_u_star[tj + tk*blockDim.x] += s_f[tj + tk*blockDim.x];

      // multiply by dt
      s_u_star[tj + tk*blockDim.x] *= dt;

      // velocity term sums into right-hand side
      s_u_star[tj + tk*blockDim.x] += u111;
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((k >= dom->Gfx._ks && k < dom->Gfx._ke)
      && (j >= dom->Gfx._js && j < dom->Gfx._je)
      && (tj > 0 && tj < (blockDim.x-1))
      && (tk > 0 && tk < (blockDim.y-1))) {
      u_star[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b]
        = s_u_star[tj + tk*blockDim.x];
      conv[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] = s_c[tj + tk*blockDim.x];
      diff[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] = s_d[tj + tk*blockDim.x];
    }
//if(i==16&&j==8&&k==8) printf("u_star diffusion %f\n",diff[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b]);
  }
}

__global__ void v_star_2(real rho_f, real nu,
  real *u0, real *v0, real *w0, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *v_star,
  dom_struct *dom, real dt0, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_v0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v back
  __shared__ real s_v1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v center
  __shared__ real s_v2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v forward
  __shared__ real s_w01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // w back
  __shared__ real s_w12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // w forward
  __shared__ real s_u01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u back
  __shared__ real s_u12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u forward
  __shared__ real s_d0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // diff0
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // conv0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_v_star[MAX_THREADS_DIM * MAX_THREADS_DIM];  // solution
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // y-force

  // working constants
  real ab0 = 0.5 * dt / dt0;  // for Adams-Bashforth stepping
  real ab = 1 + ab0;          // for Adams-Bashforth stepping
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
      s_d0[tk + ti*blockDim.x] = diff0[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_c0[tk + ti*blockDim.x] = conv0[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v0[tk + ti*blockDim.x] = v0[i + (j-1)*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v1[tk + ti*blockDim.x] = v0[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v2[tk + ti*blockDim.x] = v0[i + (j+1)*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_f[tk + ti*blockDim.x] = f[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
    if((i >= dom->Gfz._isb && i < dom->Gfz._ieb)
      && (k >= dom->Gfz._ksb && k < dom->Gfz._keb)) {
      s_w01[tk + ti*blockDim.x] = w0[i + (j-1)*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w12[tk + ti*blockDim.x] = w0[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
    }
    if((i >= dom->Gfx._isb && i < dom->Gfx._ieb)
      && (k >= dom->Gfx._ksb && k < dom->Gfx._keb)) {
      s_u01[tk + ti*blockDim.x] = u0[i + (j-1)*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u12[tk + ti*blockDim.x] = u0[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }

    s_v_star[tk + ti*blockDim.x] = 0.0;

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
      real w101 = s_w01[tk + ti*blockDim.x];
      real w111 = s_w12[tk + ti*blockDim.x];
      real w102 = s_w01[(tk+1) + ti*blockDim.x];
      real w112 = s_w12[(tk+1) + ti*blockDim.x];

      real v011 = s_v1[tk + (ti-1)*blockDim.x];
      real v211 = s_v1[tk + (ti+1)*blockDim.x];
      real u101 = s_u01[tk + ti*blockDim.x];
      real u111 = s_u12[tk + ti*blockDim.x];
      real u201 = s_u01[tk + (ti+1)*blockDim.x];
      real u211 = s_u12[tk + (ti+1)*blockDim.x];

      // compute convection term (Adams-Bashforth stepping)
      real dvudx = (v211 + v111)*(u211 + u201) - (v111 + v011)*(u111 + u101);
      dvudx *= 0.25 * ddx;

      real dvvdy = (v121 + v111)*(v121 + v111) - (v111 + v101)*(v111 + v101);
      dvvdy *= 0.25 * ddy;

      real dvwdz = (v112 + v111)*(w112 + w102) - (v111 + v110)*(w111 + w101);
      dvwdz *= 0.25 * ddz;

      s_c[tk + ti*blockDim.x] = dvudx + dvvdy + dvwdz;

      // convection term sums into right-hand side
      s_v_star[tk + ti*blockDim.x] = - ab * s_c[tk + ti*blockDim.x]
        + ab0 * s_c0[tk + ti*blockDim.x];

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

      // diffusive term sums into right-hand side
      s_v_star[tk + ti*blockDim.x] += ab * s_d[tk + ti*blockDim.x]
        - ab0 * s_d0[tk + ti*blockDim.x];

      // add on imposed pressure gradient
      s_v_star[tk + ti*blockDim.x] += s_f[tk + ti*blockDim.x];

      // multiply by dt
      s_v_star[tk + ti*blockDim.x] *= dt;

      // velocity term sums into right-hand side
      s_v_star[tk + ti*blockDim.x] += v111;
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((i >= dom->Gfy._is && i < dom->Gfy._ie)
      && (k >= dom->Gfy._ks && k < dom->Gfy._ke)
      && (tk > 0 && tk < (blockDim.x-1))
      && (ti > 0 && ti < (blockDim.y-1))) {
      v_star[i+ j*dom->Gfy._s1b + k*dom->Gfy._s2b]
        = s_v_star[tk + ti*blockDim.x];
      diff[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b] = s_d[tk + ti*blockDim.x];
      conv[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b] = s_c[tk + ti*blockDim.x];
    }
  }
}

__global__ void w_star_2(real rho_f, real nu,
  real *u0, real *v0, real *w0, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *w_star,
  dom_struct *dom, real dt0, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_w2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward
  __shared__ real s_u01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u back
  __shared__ real s_u12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u forward
  __shared__ real s_v01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v back
  __shared__ real s_v12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v forward
  __shared__ real s_d0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // diff0
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff0
  __shared__ real s_c0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // conv0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv0
  __shared__ real s_w_star[MAX_THREADS_DIM * MAX_THREADS_DIM];  // solution
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];  // solution

  // working constants
  real ab0 = 0.5 * dt / dt0;  // for Adams-Bashforth stepping
  real ab = 1 + ab0;          // for Adams-Bashforth stepping
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
      s_d0[ti + tj*blockDim.x] = diff0[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_c0[ti + tj*blockDim.x] = conv0[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w0[ti + tj*blockDim.x] = w0[i + j*dom->Gfz._s1b + (k-1)*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w0[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w2[ti + tj*blockDim.x] = w0[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
      s_f[ti + tj*blockDim.x] = f[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
    }
    if((j >= dom->Gfx._jsb && j < dom->Gfx._jeb)
      && (i >= dom->Gfx._isb && i < dom->Gfx._ieb)) {
      s_u01[ti + tj*blockDim.x] = u0[i + j*dom->Gfx._s1b + (k-1)*dom->Gfx._s2b];
      s_u12[ti + tj*blockDim.x] = u0[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._jsb && j < dom->Gfy._jeb)
      && (i >= dom->Gfy._isb && i < dom->Gfy._ieb)) {
      s_v01[ti + tj*blockDim.x] = v0[i + j*dom->Gfy._s1b + (k-1)*dom->Gfy._s2b];
      s_v12[ti + tj*blockDim.x] = v0[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }

    s_w_star[ti + tj*blockDim.x] = 0.0;

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
      real u110 = s_u01[ti + tj*blockDim.x];
      real u111 = s_u12[ti + tj*blockDim.x];
      real u210 = s_u01[(ti+1) + tj*blockDim.x];
      real u211 = s_u12[(ti+1) + tj*blockDim.x];

      real w101 = s_w1[ti + (tj-1)*blockDim.x];
      real w121 = s_w1[ti + (tj+1)*blockDim.x];
      real v110 = s_v01[ti + tj*blockDim.x];
      real v111 = s_v12[ti + tj*blockDim.x];
      real v120 = s_v01[ti + (tj+1)*blockDim.x];
      real v121 = s_v12[ti + (tj+1)*blockDim.x];

      // compute convection term (Adams-Bashforth stepping)
      real dwudx = (w211 + w111)*(u211 + u210) - (w111 + w011)*(u111 + u110);
      dwudx *= 0.25 * ddx;

      real dwvdy = (w121 + w111)*(v121 + v120) - (w111 + w101)*(v111 + v110);
      dwvdy *= 0.25 * ddy;

      real dwwdz = (w112 + w111)*(w112 + w111) - (w111 + w110)*(w111 + w110);
      dwwdz *= 0.25 * ddz;

      s_c[ti + tj*blockDim.x] = dwudx + dwvdy + dwwdz;

      // convection term sums into right-hand side
      s_w_star[ti + tj*blockDim.x] = - ab * s_c[ti + tj*blockDim.x]
        + ab0 * s_c0[ti + tj*blockDim.x];

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

      // diffusive term sums into right-hand side
      s_w_star[ti + tj*blockDim.x] += ab * s_d[ti + tj*blockDim.x]
        - ab0 * s_d0[ti + tj*blockDim.x];

      // add on imposed pressure gradient
      s_w_star[ti + tj*blockDim.x] += s_f[ti + tj*blockDim.x];

      // multiply by dt
      s_w_star[ti + tj*blockDim.x] *= dt;

      // velocity term sums into right-hand side
      s_w_star[ti + tj*blockDim.x] += w111;
//if(fabs(s_f[ti + tj*blockDim.x])>0.1f) printf("\nw_star %d %d %d %f %f %f\n",i,j,k, s_f[ti + tj*blockDim.x],w111,s_w_star[ti + tj*blockDim.x]);
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((j >= dom->Gfz._js && j < dom->Gfz._je)
      && (i >= dom->Gfz._is && i < dom->Gfz._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      w_star[i+ j*dom->Gfz._s1b + k*dom->Gfz._s2b]
        = s_w_star[ti + tj*blockDim.x];
      diff[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b] = s_d[ti + tj*blockDim.x];
      conv[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b] = s_c[ti + tj*blockDim.x];
    }
  }
}

__global__ void u_star_2_init(real rho_f, real nu,
  real *u, real *v, real *w, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *u_star,
  dom_struct *dom, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_u0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u back
  __shared__ real s_u1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u center
  __shared__ real s_u2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u forward
  __shared__ real s_v01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v back
  __shared__ real s_v12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v forward
  __shared__ real s_w01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // w back
  __shared__ real s_w12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // w forward
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_u_star[MAX_THREADS_DIM * MAX_THREADS_DIM];  // solution
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // x-force

  // working constants
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
    if((k >= dom->Gfx._ksb && k < dom->Gfx._keb)
      && (j >= dom->Gfx._jsb && j < dom->Gfx._jeb)) {
      s_u0[tj + tk*blockDim.x] = u[(i-1) + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u1[tj + tk*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u2[tj + tk*blockDim.x] = u[(i+1) + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_f[tj + tk*blockDim.x] = f[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((k >= dom->Gfy._ksb && k < dom->Gfy._keb)
      && (j >= dom->Gfy._jsb && j < dom->Gfy._jeb)) {
      s_v01[tj + tk*blockDim.x] = v[(i-1) + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v12[tj + tk*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
    if((k >= dom->Gfz._ksb && k < dom->Gfz._keb)
      && (j >= dom->Gfz._jsb && j < dom->Gfz._jeb)) {
      s_w01[tj + tk*blockDim.x] = w[(i-1) + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w12[tj + tk*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
    }

    s_u_star[tj + tk*blockDim.x] = 0.0;

    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute right-hand side
    // if off the shared memory block boundary
    if((tj > 0 && tj < blockDim.x-1) && (tk > 0 && tk < blockDim.y-1)) {
      // grab the required data points for calculations
      real u011 = s_u0[tj + tk*blockDim.x];
      real u111 = s_u1[tj + tk*blockDim.x];
      real u211 = s_u2[tj + tk*blockDim.x];

      real u101 = s_u1[(tj-1) + tk*blockDim.x];
      real u121 = s_u1[(tj+1) + tk*blockDim.x];
      real v011 = s_v01[tj + tk*blockDim.x];
      real v111 = s_v12[tj + tk*blockDim.x];
      real v021 = s_v01[(tj+1) + tk*blockDim.x];
      real v121 = s_v12[(tj+1) + tk*blockDim.x];

      real u110 = s_u1[tj + (tk-1)*blockDim.x];
      real u112 = s_u1[tj + (tk+1)*blockDim.x];
      real w011 = s_w01[tj + tk*blockDim.x];
      real w111 = s_w12[tj + tk*blockDim.x];
      real w012 = s_w01[tj + (tk+1)*blockDim.x];
      real w112 = s_w12[tj + (tk+1)*blockDim.x];

      // compute convection term (Euler stepping)
      real duudx = (u211 + u111)*(u211 + u111) - (u111 + u011)*(u111 + u011);
      duudx *= 0.25 * ddx;

      real duvdy = (u121 + u111)*(v121 + v021) - (u111 + u101)*(v111 + v011);
      duvdy *= 0.25 * ddy;

      real duwdz = (u112 + u111)*(w112 + w012) - (u111 + u110)*(w111 + w011);
      duwdz *= 0.25 * ddz;

      s_c[tj + tk*blockDim.x] = duudx + duvdy + duwdz;

      // convection term sums into right-hand side
      s_u_star[tj + tk*blockDim.x] = - s_c[tj + tk*blockDim.x];

      // compute diffusive term
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

      // diffusive term sums into right-hand side
      s_u_star[tj + tk*blockDim.x] += s_d[tj + tk*blockDim.x];

      // add on imposed pressure gradient
      s_u_star[tj + tk*blockDim.x] += s_f[tj + tk*blockDim.x];

      // multiply by dt
      s_u_star[tj + tk*blockDim.x] *= dt;

      // velocity term sums into right-hand side
      s_u_star[tj + tk*blockDim.x] += u111;
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((k >= dom->Gfx._ks && k < dom->Gfx._ke)
      && (j >= dom->Gfx._js && j < dom->Gfx._je)
      && (tj > 0 && tj < (blockDim.x-1))
      && (tk > 0 && tk < (blockDim.y-1))) {
      u_star[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b]
        = s_u_star[tj + tk*blockDim.x];
      diff[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] = s_d[tj + tk*blockDim.x];
      conv[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] = s_c[tj + tk*blockDim.x];
    }
  }
}

__global__ void v_star_2_init(real rho_f, real nu,
  real *u, real *v, real *w, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *v_star,
  dom_struct *dom, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_v0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v back
  __shared__ real s_v1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v center
  __shared__ real s_v2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v forward
  __shared__ real s_w01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // w back
  __shared__ real s_w12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // w forward
  __shared__ real s_u01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u back
  __shared__ real s_u12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u forward
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_v_star[MAX_THREADS_DIM * MAX_THREADS_DIM];  // solution
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // y-force

  // working constants
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

  // loop over v-planes
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
      s_v0[tk + ti*blockDim.x] = v[i + (j-1)*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v1[tk + ti*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v2[tk + ti*blockDim.x] = v[i + (j+1)*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_f[tk + ti*blockDim.x] = f[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
    if((i >= dom->Gfz._isb && i < dom->Gfz._ieb)
      && (k >= dom->Gfz._ksb && k < dom->Gfz._keb)) {
      s_w01[tk + ti*blockDim.x] = w[i + (j-1)*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w12[tk + ti*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
    }
    if((i >= dom->Gfx._isb && i < dom->Gfx._ieb)
      && (k >= dom->Gfx._ksb && k < dom->Gfx._keb)) {
      s_u01[tk + ti*blockDim.x] = u[i + (j-1)*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u12[tk + ti*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }

    s_v_star[tk + ti*blockDim.x] = 0.0;

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
      real w101 = s_w01[tk + ti*blockDim.x];
      real w111 = s_w12[tk + ti*blockDim.x];
      real w102 = s_w01[(tk+1) + ti*blockDim.x];
      real w112 = s_w12[(tk+1) + ti*blockDim.x];

      real v011 = s_v1[tk + (ti-1)*blockDim.x];
      real v211 = s_v1[tk + (ti+1)*blockDim.x];
      real u101 = s_u01[tk + ti*blockDim.x];
      real u111 = s_u12[tk + ti*blockDim.x];
      real u201 = s_u01[tk + (ti+1)*blockDim.x];
      real u211 = s_u12[tk + (ti+1)*blockDim.x];

      // compute convection term (Euler stepping)
      real dvudx = (v211 + v111)*(u211 + u201) - (v111 + v011)*(u111 + u101);
      dvudx *= 0.25 * ddx;

      real dvvdy = (v121 + v111)*(v121 + v111) - (v111 + v101)*(v111 + v101);
      dvvdy *= 0.25 * ddy;

      real dvwdz = (v112 + v111)*(w112 + w102) - (v111 + v110)*(w111 + w101);
      dvwdz *= 0.25 * ddz;

      s_c[tk + ti*blockDim.x] = dvudx + dvvdy + dvwdz;

      // convection term sums into right-hand side
      s_v_star[tk + ti*blockDim.x] = - s_c[tk + ti*blockDim.x];

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

      // diffusive term sums into right-hand side
      s_v_star[tk + ti*blockDim.x] += s_d[tk + ti*blockDim.x];

      // add on imposed pressure gradient
      s_v_star[tk + ti*blockDim.x] += s_f[tk + ti*blockDim.x];

      // multiply by dt
      s_v_star[tk + ti*blockDim.x] *= dt;

      // velocity term sums into right-hand side
      s_v_star[tk + ti*blockDim.x] += v111;
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((i >= dom->Gfy._is && i < dom->Gfy._ie)
      && (k >= dom->Gfy._ks && k < dom->Gfy._ke)
      && (tk > 0 && tk < (blockDim.x-1))
      && (ti > 0 && ti < (blockDim.y-1))) {
      v_star[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b]
        = s_v_star[tk + ti*blockDim.x];
      diff[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b] = s_d[tk + ti*blockDim.x];
      conv[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b] = s_c[tk + ti*blockDim.x];
    }
  }
}

__global__ void w_star_2_init(real rho_f, real nu,
  real *u, real *v, real *w, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *w_star,
  dom_struct *dom, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_w2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward
  __shared__ real s_u01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u back
  __shared__ real s_u12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // u forward
  __shared__ real s_v01[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v back
  __shared__ real s_v12[MAX_THREADS_DIM * MAX_THREADS_DIM];     // v forward
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_w_star[MAX_THREADS_DIM * MAX_THREADS_DIM];  // solution
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // z-force

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
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k-1)*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w2[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
      s_f[ti + tj*blockDim.x] = f[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
    }
    if((j >= dom->Gfx._jsb && j < dom->Gfx._jeb)
      && (i >= dom->Gfx._isb && i < dom->Gfx._ieb)) {
      s_u01[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + (k-1)*dom->Gfx._s2b];
      s_u12[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._jsb && j < dom->Gfy._jeb)
      && (i >= dom->Gfy._isb && i < dom->Gfy._ieb)) {
      s_v01[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + (k-1)*dom->Gfy._s2b];
      s_v12[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }

    s_w_star[ti + tj*blockDim.x] = 0.0;

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
      real u110 = s_u01[ti + tj*blockDim.x];
      real u111 = s_u12[ti + tj*blockDim.x];
      real u210 = s_u01[(ti+1) + tj*blockDim.x];
      real u211 = s_u12[(ti+1) + tj*blockDim.x];

      real w101 = s_w1[ti + (tj-1)*blockDim.x];
      real w121 = s_w1[ti + (tj+1)*blockDim.x];
      real v110 = s_v01[ti + tj*blockDim.x];
      real v111 = s_v12[ti + tj*blockDim.x];
      real v120 = s_v01[ti + (tj+1)*blockDim.x];
      real v121 = s_v12[ti + (tj+1)*blockDim.x];

      // compute convection term (Euler stepping)
      real dwudx = (w211 + w111)*(u211 + u210) - (w111 + w011)*(u111 + u110);
      dwudx *= 0.25 * ddx;

      real dwvdy = (w121 + w111)*(v121 + v120) - (w111 + w101)*(v111 + v110);
      dwvdy *= 0.25 * ddy;

      real dwwdz = (w112 + w111)*(w112 + w111) - (w111 + w110)*(w111 + w110);
      dwwdz *= 0.25 * ddz;

      s_c[ti + tj*blockDim.x] = dwudx + dwvdy + dwwdz;

      // convection term sums into right-hand side
      s_w_star[ti + tj*blockDim.x] = - s_c[ti + tj*blockDim.x];

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

      // diffusive term sums into right-hand side
      s_w_star[ti + tj*blockDim.x] += s_d[ti + tj*blockDim.x];

      // add on imposed pressure gradient
      s_w_star[ti + tj*blockDim.x] += s_f[ti + tj*blockDim.x];

      // multiply by dt
      s_w_star[ti + tj*blockDim.x] *= dt;

      // velocity term sums into right-hand side
      s_w_star[ti + tj*blockDim.x] += w111;
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((j >= dom->Gfz._js && j < dom->Gfz._je)
      && (i >= dom->Gfz._is && i < dom->Gfz._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      w_star[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b]
        = s_w_star[ti + tj*blockDim.x];
      diff[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b] = s_d[ti + tj*blockDim.x];
      conv[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b] = s_c[ti + tj*blockDim.x];
    }
  }
}


//u=cos(qx)cos(qy)sin(rt); v=sin(qx)sin(qy)sin(rt)
__global__ void forcing_test_taylor(real *fx, real *fy, dom_struct *dom,real t,real DIFF,point_struct *points, real n2,real n3)
{

//real n3=2;
//real n2=1;
 
real q=1;
real dia=2*points[0].r;
real L=2*PI;
real T=n3*dia*dia/DIFF;
real r=2*PI/T;
real U0=DIFF*L/dia/dia/n2;
 
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  for(int i = dom->Gfx._isb; i < dom->Gfx._ieb; i++) {
    if(tj < dom->Gfx._jnb && tk < dom->Gfx._knb) {
  real x = ((i-1) * dom->dx) + dom->xs;
  real y = ((tj-0.5) * dom->dy) + dom->ys;
      fx[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = U0*
cos(q*x)*
(r*cos(r*t)*cos(q*y)+
q*sin(r*t)*(2*DIFF*q*cos(q*y)-U0*sin(r*t)*sin(q*x)));
    }
  }

__syncthreads();

       tk = blockIdx.x * blockDim.x + threadIdx.x;
   int ti = blockIdx.y * blockDim.y + threadIdx.y;

//if(tk==dom->Gfy._knb-1&&ti== dom->Gfy._inb-1) printf("\ntest2 %d %d \n",ti,tk);


  for(int j = dom->Gfy._jsb; j < dom->Gfy._jeb; j++) {
    if(tk < dom->Gfy._knb && ti < dom->Gfy._inb) {
  real x = ((ti-0.5) * dom->dx) + dom->xs;
  real y = ((j-1) * dom->dy) + dom->ys;
      fy[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b] = U0*
sin(q*y)*(q*U0*cos(q*x)*cos(q*x)*sin(r*t)*sin(r*t)*cos(q*y)
+sin(q*x)*(r*cos(r*t)+q*sin(r*t)*(2*DIFF*q+U0*cos(q*y)*sin(r*t)*sin(q*x))));
    }
  }
}


//u=cos(qx)cos(qy)sin(rt); v=sin(qx)sin(qy)sin(rt)
__global__ void forcing_test_taylor_sc_ConvDiff(real *fx, real *fy, dom_struct *dom,real t,real nu,point_struct *points)
{
/*
 real q=6;
 real r=10*2*PI;
 real U0=20;
*/
 real q=2;
 real r=10*2*PI;
 real U0=10;
// real U0=0.f;

  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  for(int i = dom->Gfx._isb; i < dom->Gfx._ieb; i++) {
    if(tj < dom->Gfx._jnb && tk < dom->Gfx._knb) {
  real x = ((i-1) * dom->dx) + dom->xs;
  real y = ((tj-0.5) * dom->dy) + dom->ys;
      fx[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = U0*
cos(q*x)*
(r*cos(r*t)*cos(q*y)+
q*sin(r*t)*(2*nu*q*cos(q*y)-U0*sin(r*t)*sin(q*x)));
    }
  }

__syncthreads();

       tk = blockIdx.x * blockDim.x + threadIdx.x;
   int ti = blockIdx.y * blockDim.y + threadIdx.y;

//if(tk==dom->Gfy._knb-1&&ti== dom->Gfy._inb-1) printf("\ntest2 %d %d \n",ti,tk);


  for(int j = dom->Gfy._jsb; j < dom->Gfy._jeb; j++) {
    if(tk < dom->Gfy._knb && ti < dom->Gfy._inb) {
  real x = ((ti-0.5) * dom->dx) + dom->xs;
  real y = ((j-1) * dom->dy) + dom->ys;
      fy[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b] = U0*
sin(q*y)*(q*U0*cos(q*x)*cos(q*x)*sin(r*t)*sin(r*t)*cos(q*y)
+sin(q*x)*(r*cos(r*t)+q*sin(r*t)*(2*nu*q+U0*cos(q*y)*sin(r*t)*sin(q*x))));
    }
  }
}



__global__ void forcing_reset_x(real *fx, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  for(int i = dom->Gfx._isb; i < dom->Gfx._ieb; i++) {
    if(tj < dom->Gfx._jnb && tk < dom->Gfx._knb) {
      fx[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = 0.;
    }
  }
}

__global__ void forcing_reset_y(real *fy, dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  for(int j = dom->Gfy._jsb; j < dom->Gfy._jeb; j++) {
    if(tk < dom->Gfy._knb && ti < dom->Gfy._inb) {
      fy[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b] = 0.;
    }
  }
}

__global__ void forcing_reset_z(real *fz, dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  for(int k = dom->Gfz._ksb; k < dom->Gfz._keb; k++) {
    if(ti < dom->Gfz._inb && tj < dom->Gfz._jnb) {
      fz[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b] = 0.;
    }
  }
}

__global__ void forcing_add_x_const(real val, real *fx, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  for(int i = dom->Gfx._isb; i < dom->Gfx._ieb; i++) {
    if(tj < dom->Gfx._jnb && tk < dom->Gfx._knb) {
      fx[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] += val;
    }
  }
}

__global__ void forcing_add_y_const(real val, real *fy, dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  for(int j = dom->Gfy._jsb; j < dom->Gfy._jeb; j++) {
    if(tk < dom->Gfy._knb && ti < dom->Gfy._inb) {
      fy[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b] += val;
    }
  }
}

__global__ void forcing_add_z_const(real val, real *fz, dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  for(int k = dom->Gfz._ksb; k < dom->Gfz._keb; k++) {
    if(ti < dom->Gfz._inb && tj < dom->Gfz._jnb) {
      fz[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b] += val;
    }
  }
}


__global__ void forcing_add_x_field(real scale, real *val, real *fx,
  dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  for(int i = dom->Gfx._isb; i < dom->Gfx._ieb; i++) {
    if(tj < dom->Gfx._jnb && tk < dom->Gfx._knb) {
      fx[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b]
        += scale * val[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b];
    }
  }
}

__global__ void forcing_add_y_field(real scale, real *val, real *fy,
  dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  for(int j = dom->Gfy._jsb; j < dom->Gfy._jeb; j++) {
    if(tk < dom->Gfy._knb && ti < dom->Gfy._inb) {
      fy[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b]
        += scale * val[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b];
    }
  }
}

__global__ void forcing_add_z_field(real scale, real *val, real *fz,
  dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  for(int k = dom->Gfz._ksb; k < dom->Gfz._keb; k++) {
    if(ti < dom->Gfz._inb && tj < dom->Gfz._jnb) {
      fz[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b]
        += scale * val[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b];
    }
  }
}

__global__ void surf_int_x_copy(real *u_star, real *u_star_tmp,
  dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  if(tj < dom->Gfx.jn && tk < dom->Gfx.kn) {
    int C = dom->Gfx.is + (tj+DOM_BUF)*dom->Gfx._s1b + (tk+DOM_BUF)*dom->Gfx._s2b;
    int CC = 0 + tj + tk*dom->Gfx.jn;
    u_star_tmp[CC] = -u_star[C];
    C = dom->Gfx.ie-1 + (tj+DOM_BUF)*dom->Gfx._s1b + (tk+DOM_BUF)*dom->Gfx._s2b;
    CC = dom->Gfx.jn*dom->Gfx.kn + tj + tk*dom->Gfx.jn;
    u_star_tmp[CC] = u_star[C];
  }
}

__global__ void surf_int_y_copy(real *v_star, real *v_star_tmp,
  dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if(tk < dom->Gfy.kn && ti < dom->Gfy.in) {
    int C = (ti+DOM_BUF) + dom->Gfy.js*dom->Gfy._s1b + (tk+DOM_BUF)*dom->Gfy._s2b;
    int CC = ti + 0 + tk*dom->Gfy.in;
    v_star_tmp[CC] = -v_star[C];
    C = (ti+DOM_BUF) + (dom->Gfy.je-1)*dom->Gfy._s1b + (tk+DOM_BUF)*dom->Gfy._s2b;
    CC = ti + dom->Gfy.in*dom->Gfy.kn + tk*dom->Gfy.in;
    v_star_tmp[CC] = v_star[C];
  }
}

__global__ void surf_int_z_copy(real *w_star, real *w_star_tmp,
  dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if(ti < dom->Gfz.in && tj < dom->Gfz.jn) {
    int C = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz._s1b + dom->Gfz.ks*dom->Gfz._s2b;
    int CC = ti + tj*dom->Gfz.in + 0;
    w_star_tmp[CC] = -w_star[C];
    C = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz._s1b + (dom->Gfz.ke-1)*dom->Gfz._s2b;
    CC = ti + tj*dom->Gfz.in + dom->Gfz.in*dom->Gfz.jn;
    w_star_tmp[CC] = w_star[C];
  }
}

__global__ void plane_eps_x_W(real eps, real *u_star, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y; 

  if(tj < dom->Gfx.jn && tk < dom->Gfx.kn) {
    int C = dom->Gfx.is + (tj+DOM_BUF)*dom->Gfx._s1b + (tk+DOM_BUF)*dom->Gfx._s2b;
    u_star[C] = u_star[C] + eps;
  }
}

__global__ void plane_eps_x_E(real eps, real *u_star, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y; 

  if(tj < dom->Gfx.jn && tk < dom->Gfx.kn) {
    int C = dom->Gfx.ie-1 + (tj+DOM_BUF)*dom->Gfx._s1b + (tk+DOM_BUF)*dom->Gfx._s2b;
    u_star[C] = u_star[C] - eps;
  }
}

__global__ void plane_eps_y_S(real eps, real *v_star, dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if(tk < dom->Gfy.kn && ti < dom->Gfy.in) {
    int C = (ti+DOM_BUF) + (dom->Gfy.js)*dom->Gfy._s1b + (tk+DOM_BUF)*dom->Gfy._s2b;
    v_star[C] = v_star[C] + eps;
  }
}

__global__ void plane_eps_y_N(real eps, real *v_star, dom_struct *dom)
{
  int tk = blockIdx.x * blockDim.x + threadIdx.x;
  int ti = blockIdx.y * blockDim.y + threadIdx.y;

  if(tk < dom->Gfy.kn && ti < dom->Gfy.in) {
    int C = (ti+DOM_BUF) + (dom->Gfy.je-1)*dom->Gfy._s1b + (tk+DOM_BUF)*dom->Gfy._s2b;
    v_star[C] = v_star[C] - eps;
  }
}

__global__ void plane_eps_z_B(real eps, real *w_star, dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if(ti < dom->Gfz.in && tj < dom->Gfz.jn) {
    int C = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz._s1b + (dom->Gfz.ks)*dom->Gfz._s2b;
    w_star[C] = w_star[C] + eps;
  }
}

__global__ void plane_eps_z_T(real eps, real *w_star, dom_struct *dom)
{
  int ti = blockIdx.x * blockDim.x + threadIdx.x;
  int tj = blockIdx.y * blockDim.y + threadIdx.y;

  if(ti < dom->Gfz.in && tj < dom->Gfz.jn) {
    int C = (ti+DOM_BUF) + (tj+DOM_BUF)*dom->Gfz._s1b + (dom->Gfz.ke-1)*dom->Gfz._s2b;
    w_star[C] = w_star[C] - eps;
  }
}

/*
__global__ void move_parts_a(dom_struct *dom, part_struct *parts, int nparts,
  real dt, real dt0, g_struct g, real rho_f, real ttime)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x; // particle number
  real m = 4./3. * PI * parts[pp].rho * parts[pp].r*parts[pp].r*parts[pp].r;

  if(pp < nparts) {
    if(parts[pp].translating) {
      // update linear accelerations
      parts[pp].udot = (parts[pp].Fx + parts[pp].kFx + parts[pp].iFx
        + parts[pp].aFx) / m
        + (parts[pp].rho - rho_f) / rho_f * g.x;
      parts[pp].vdot = (parts[pp].Fy + parts[pp].kFy + parts[pp].iFy
        + parts[pp].aFy) / m
        + (parts[pp].rho - rho_f) / rho_f * g.y;
      parts[pp].wdot = (parts[pp].Fz + parts[pp].kFz + parts[pp].iFz
        + parts[pp].aFz) / m
        + (parts[pp].rho - rho_f) / rho_f * g.z;

      // update linear velocities
      parts[pp].u = parts[pp].u0 + parts[pp].udot * dt;
      parts[pp].v = parts[pp].v0 + parts[pp].vdot * dt;
      parts[pp].w = parts[pp].w0 + parts[pp].wdot * dt;

      // do not update position
    }
  }
}

__global__ void move_parts_b(dom_struct *dom, part_struct *parts, int nparts,
  real dt, real dt0, g_struct g, real rho_f, real ttime)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x; // particle number
  real m = 4./3. * PI * parts[pp].rho * parts[pp].r*parts[pp].r*parts[pp].r;
  //real dT = dt / dt0;

  if(pp < nparts) {
    if(parts[pp].translating) {
      // update linear accelerations
      parts[pp].udot = (parts[pp].Fx + parts[pp].kFx + parts[pp].iFx
        + parts[pp].aFx) / m
        + (parts[pp].rho - rho_f) / parts[pp].rho * g.x;
      parts[pp].vdot = (parts[pp].Fy + parts[pp].kFy + parts[pp].iFy
        + parts[pp].aFy) / m
        + (parts[pp].rho - rho_f) / parts[pp].rho * g.y;
      parts[pp].wdot = (parts[pp].Fz + parts[pp].kFz + parts[pp].iFz
        + parts[pp].aFz) / m
        + (parts[pp].rho - rho_f) / parts[pp].rho * g.z;

      // update linear velocities
      parts[pp].u = parts[pp].u0 + parts[pp].udot * dt;
      parts[pp].v = parts[pp].v0 + parts[pp].vdot * dt;
      parts[pp].w = parts[pp].w0 + parts[pp].wdot * dt;

      parts[pp].u0 = parts[pp].u;
      parts[pp].v0 = parts[pp].v;
      parts[pp].w0 = parts[pp].w;

      // update position
      parts[pp].x = parts[pp].x + parts[pp].u * dt;
      if(parts[pp].x < dom->xs) parts[pp].x = parts[pp].x + dom->xl;
      else if(parts[pp].x > dom->xe) parts[pp].x = parts[pp].x - dom->xl;
      parts[pp].y = parts[pp].y + parts[pp].v * dt;
      if(parts[pp].y < dom->ys) parts[pp].y = parts[pp].y + dom->yl;
      else if(parts[pp].y > dom->ye) parts[pp].y = parts[pp].y - dom->yl;
      parts[pp].z = parts[pp].z + parts[pp].w * dt;
      if(parts[pp].z < dom->zs) parts[pp].z = parts[pp].z + dom->zl;
      else if(parts[pp].z > dom->ze) parts[pp].z = parts[pp].z - dom->zl;
    }
    if(parts[pp].rotating) {
      // update angular accelerations
      real I = 0.4 * m * parts[pp].r*parts[pp].r;
      parts[pp].oxdot = (parts[pp].Lx + parts[pp].aLx) / I;
      parts[pp].oydot = (parts[pp].Ly + parts[pp].aLy) / I;
      parts[pp].ozdot = (parts[pp].Lz + parts[pp].aLz) / I;

      // update angular velocities
      parts[pp].ox = parts[pp].ox + parts[pp].oxdot * dt;
      parts[pp].oy = parts[pp].oy + parts[pp].oydot * dt;
      parts[pp].oz = parts[pp].oz + parts[pp].ozdot * dt;

      /// update basis vectors 
      // calculate rotation magnitude
      real mag = sqrt(parts[pp].ox*parts[pp].ox + parts[pp].oy*parts[pp].oy
        + parts[pp].oz*parts[pp].oz);
      // calculate normalized rotation axis
      real X = 0;
      real Y = 0;
      real Z = 0;
      if(mag > 0) {
        X = parts[pp].ox / mag;
        Y = parts[pp].oy / mag;
        Z = parts[pp].oz / mag;
      }
      // calculate rotation quaternion
      real theta = mag * dt;
      real qr = cos(0.5*theta);
      real qi = X * sin(0.5*theta);
      real qj = Y * sin(0.5*theta);
      real qk = Z * sin(0.5*theta);
      // compute quaternion conjugation to apply rotation to basis vectors
      rotate(qr, qi, qj, qk, &parts[pp].axx, &parts[pp].axy, &parts[pp].axz);
      rotate(qr, qi, qj, qk, &parts[pp].ayx, &parts[pp].ayy, &parts[pp].ayz);
      rotate(qr, qi, qj, qk, &parts[pp].azx, &parts[pp].azy, &parts[pp].azz);
    }
  }
}

__device__ void rotate(real qr, real qi, real qj, real qk,
  real *pi, real *pj, real *pk)
{
  real Pr = *pi*qi + *pj*qj + *pk*qk;
  real Pi = *pi*qr - *pj*qk + *pk*qj;
  real Pj = *pi*qk + *pj*qr - *pk*qi;
  real Pk = -*pi*qj + *pj*qi + *pk*qr;

  *pi = qr*Pi + qi*Pr + qj*Pk - qk*Pj;
  *pj = qr*Pj - qi*Pk + qj*Pr + qk*Pi;
  *pk = qr*Pk + qi*Pj - qj*Pi + qk*Pr;
}

__global__ void collision_init(part_struct *parts, int nparts)
{
  int j = threadIdx.x + blockIdx.x*blockDim.x;
  if(j < nparts) {
    parts[j].iFx = 0.;
    parts[j].iFy = 0.;
    parts[j].iFz = 0.;
    parts[j].iLx = 0.;
    parts[j].iLy = 0.;
    parts[j].iLz = 0.;
  }
}

__global__ void collision_parts(part_struct *parts, int i,
  dom_struct *dom, real eps, real *forces, real *moments, int nparts, real mu,
  BC bc)
{
  int j = threadIdx.x + blockIdx.x*blockDim.x;

  if(j < nparts) {
    if(i != j) {
      real ai = parts[i].r;
      real aj = parts[j].r;
      real B = aj / ai;
      real hN = 2.*(parts[i].rs - parts[i].r + parts[j].rs - parts[j].r);

      real ux, uy, uz;
      real rx, rx1, rx2, ry, ry1, ry2, rz, rz1, rz2, r;
      real h, ah, lnah;
      real nx, ny, nz, udotn;
      real unx, uny, unz, utx, uty, utz, ut;
      real tx, ty, tz, t, ox, oy, oz, o;
      real omegax, omegay, omegaz, omega, odoto;
      real omx, omy, omz;
      real opB;
      real Fnx, Fny, Fnz, Ftx, Fty, Ftz, Lox, Loy, Loz;

      real xi = parts[i].x;
      real xj = parts[j].x;
      // check for neighbors across the domain when using periodic boundaries
      rx = xi - xj;
      rx1 = xi - (xj + dom->xl);
      rx2 = xi - (xj - dom->xl);
      if(rx1*rx1 < rx*rx) rx = rx1;
      if(rx2*rx2 < rx*rx) rx = rx2;
      rx = (bc.uW == PERIODIC) * rx + (bc.uW != PERIODIC) * (xi - xj);

      real yi = parts[i].y;
      real yj = parts[j].y;
      // check for neighbors across the domain when using periodic boundaries
      ry = yi - yj;
      ry1 = yi - (yj + dom->yl);
      ry2 = yi - (yj - dom->yl);
      if(ry1*ry1 < ry*ry) ry = ry1;
      if(ry2*ry2 < ry*ry) ry = ry2;
      ry = (bc.vS == PERIODIC) * ry + (bc.vS != PERIODIC) * (yi - yj);

      real zi = parts[i].z;
      real zj = parts[j].z;
      // check for neighbors across the domain when using periodic boundaries
      rz = zi - zj;
      rz1 = zi - (zj + dom->zl);
      rz2 = zi - (zj - dom->zl);
      if(rz1*rz1 < rz*rz) rz = rz1;
      if(rz2*rz2 < rz*rz) rz = rz2;
      rz = (bc.wB == PERIODIC) * rz + (bc.wB != PERIODIC) * (zi - zj);

      ux = parts[i].u - parts[j].u;
      uy = parts[i].v - parts[j].v;
      uz = parts[i].w - parts[j].w;

      r = sqrt(rx*rx + ry*ry + rz*rz);

      omegax = parts[i].ox - parts[j].ox;
      omegay = parts[i].oy - parts[j].oy;
      omegaz = parts[i].oz - parts[j].oz;

      omega = sqrt(omegax*omegax + omegay*omegay + omegaz*omegaz);

      h = r - ai - aj;

      nx = rx / r;
      ny = ry / r;
      nz = rz / r;

      udotn = ux * nx + uy * ny + uz * nz;

      unx = udotn * nx;
      uny = udotn * ny;
      unz = udotn * nz;

      utx = ux - unx;
      uty = uy - uny;
      utz = uz - unz;

      ut = sqrt(utx*utx + uty*uty + utz*utz);

      if(ut > 0) {
        tx = utx / ut;
        ty = uty / ut;
        tz = utz / ut;

        ox = ny*tz - nz*ty;
        oy = -(nx*tz - nz*tx);
        oz = nx*ty - ny*tx;

        o = sqrt(ox*ox + oy*oy + oz*oz);

        ox = ox / o;
        oy = oy / o;
        oz = oz / o;

        odoto = omegax * ox + omegay * oy + omegaz * oz;

        omx = odoto * ox;
        omy = odoto * oy;
        omz = odoto * oz;
      } else if(omega > 0) {
        ox = omegax / omega;
        oy = omegay / omega;
        oz = omegaz / omega;

        omx = omegax;
        omy = omegay;
        omz = omegaz;

        tx = oy*nz - oz*ny;
        ty = -(ox*nz - oz*nx);
        tz = ox*ny - oy*nx;

        t = sqrt(tx*tx + ty*ty + tz*tz);
        tx = tx / t;
        ty = ty / t;
        tz = tz / t;
      } else {
        tx = 1.;
        ty = 0.;
        tz = 0.;

        ox = ny*tz - nz*ty;
        oy = -(nx*tz - nz*tx);
        oz = nx*ty - ny*tx;

        o = sqrt(ox*ox + oy*oy + oz*oz);

        ox = ox / o;
        oy = oy / o;
        oz = oz / o;

        omx = 0.;
        omy = 0.;
        omz = 0.;
      }

      ah = ai/h - ai/hN;
      lnah = log(ai/h);

      opB = 1 + B;

      if(h < hN) {
        Fnx = - B*B / (opB*opB) * ah - (1.+7.*B+B*B)/(5.*opB*opB*opB)*lnah;
        Fny = Fnx;
        Fnz = Fnx;
        Fnx *= 6.*PI*mu*ai*unx;
        Fny *= 6.*PI*mu*ai*uny;
        Fnz *= 6.*PI*mu*ai*unz;

        Lox = 0;
        Loy = 0;
        Loz = 0;
        Ftx = -6.*PI*mu*ai*utx*4.*B*(2.+B+2.*B*B)/(15.*opB*opB*opB)*lnah;
        Fty = -6.*PI*mu*ai*uty*4.*B*(2.+B+2.*B*B)/(15.*opB*opB*opB)*lnah;
        Ftz = -6.*PI*mu*ai*utz*4.*B*(2.+B+2.*B*B)/(15.*opB*opB*opB)*lnah;
        Ftx += 8.*PI*mu*ai*ai*omy*B*(4.+B)/(10.*opB*opB)*lnah;
        Ftx += 8.*PI*mu*ai*ai*omz*B*(4.+B)/(10.*opB*opB)*lnah;
        Fty += 8.*PI*mu*ai*ai*omx*B*(4.+B)/(10.*opB*opB)*lnah;
        Fty += 8.*PI*mu*ai*ai*omz*B*(4.+B)/(10.*opB*opB)*lnah;
        Ftz += 8.*PI*mu*ai*ai*omx*B*(4.+B)/(10.*opB*opB)*lnah;
        Ftz += 8.*PI*mu*ai*ai*omy*B*(4.+B)/(10.*opB*opB)*lnah;

        Lox = 8.*PI*mu*ai*ai*uty*B*(4.+B)/(10.*opB*opB)*lnah;
        Lox += 8.*PI*mu*ai*ai*utz*B*(4.+B)/(10.*opB*opB)*lnah;
        Loy = 8.*PI*mu*ai*ai*utx*B*(4.+B)/(10.*opB*opB)*lnah;
        Loy += 8.*PI*mu*ai*ai*utz*B*(4.+B)/(10.*opB*opB)*lnah;
        Loz = 8.*PI*mu*ai*ai*utx*B*(4.+B)/(10.*opB*opB)*lnah;
        Loz += 8.*PI*mu*ai*ai*uty*B*(4.+B)/(10.*opB*opB)*lnah;
        Lox += -8.*PI*mu*ai*ai*ai*omx*2.*B/(5.*opB)*lnah;
        Loy += -8.*PI*mu*ai*ai*ai*omy*2.*B/(5.*opB)*lnah;
        Loz += -8.*PI*mu*ai*ai*ai*omz*2.*B/(5.*opB)*lnah;
      } else {
        Fnx = 0;
        Fny = 0;
        Fnz = 0;
        Ftx = 0;
        Fty = 0;
        Ftz = 0;
        Lox = 0;
        Loy = 0;
        Loz = 0;
      }

      forces[    j*3] = Fnx + Ftx;
      forces[1 + j*3] = Fny + Fty;
      forces[2 + j*3] = Fnz + Ftz;
      moments[    j*3] = Lox;
      moments[1 + j*3] = Loy;
      moments[2 + j*3] = Loz;

      // extra short-range repulsive force
      // Make this force completely take the place of the lubrication force
      ah = 1./h - 1./eps;
      if(h < eps) {
        forces[    j*3] = parts[i].rho*4./3.*PI*ai*ai*ai*nx*ah + Ftx;
        forces[1 + j*3] = parts[i].rho*4./3.*PI*ai*ai*ai*ny*ah + Fty;
        forces[2 + j*3] = parts[i].rho*4./3.*PI*ai*ai*ai*nz*ah + Ftz;
      }
    } else {
      forces[    j*3] = 0;
      forces[1 + j*3] = 0;
      forces[2 + j*3] = 0;
      moments[    j*3] = 0;
      moments[1 + j*3] = 0;
      moments[2 + j*3] = 0;
    }

    __syncthreads();

    if(j == 0) {
      for(int k = 0; k < nparts; k++) {
        parts[i].iFx += forces[    k*3];
        parts[i].iFy += forces[1 + k*3];
        parts[i].iFz += forces[2 + k*3];
        parts[i].iLx += moments[    k*3];
        parts[i].iLy += moments[1 + k*3];
        parts[i].iLz += moments[2 + k*3];
      }
    }
  }
}

__global__ void collision_walls(dom_struct *dom, part_struct *parts,
  int nparts, BC bc, real eps, real mu)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  //// parallelize this further by using a CUDA block for each wall 

  if(i < nparts) {
    real dx = 0;
    real dy = 0;
    real dz = 0;
    real Un, Utx, Uty, Utz;
    real omx, omy, omz;

    real ai = parts[i].r;
    real h = 0;
    real hN = 2.*(parts[i].rs - parts[i].r);
    real ah, lnah;

    real Fnx, Fny, Fnz, Ftx, Fty, Ftz;
    real Lox, Loy, Loz;

    // west wall
    dx = fabs(parts[i].x - dom->xs);
    h = dx - ai;
    if(h < hN) {
      ah = ai/h - ai/hN;
      lnah = log(ai/h);

      Un = parts[i].u - bc.uWD;
      Utx = 0.;
      Uty = parts[i].v - bc.vWD;
      Utz = parts[i].w - bc.wWD;
      omx = parts[i].ox;
      omy = parts[i].oy;
      omz = parts[i].oz;

      Fnx = -6.*PI*mu*ai*Un*ah;
      Fny = 0.;
      Fnz = 0.;

      /// Need to reformulate these tangential forces using wall shear velocity  
      Ftx = 0.;
      Fty = -6.*PI*mu*ai*Uty*8./15.*lnah;
      Ftz = -6.*PI*mu*ai*Utz*8./15.*lnah;
      Ftx += 0.;
      Fty += 8.*PI*mu*ai*ai*omz*1./10.*lnah;
      Ftz += -8.*PI*mu*ai*ai*omy*1./10.*lnah;

      Lox = 0.;
      Loy = -8.*PI*mu*ai*ai*Utz*1./10.*lnah;
      Loz = 8.*PI*mu*ai*ai*Uty*1./10.*lnah;
      Lox += 0.;
      Loy += -8.*PI*mu*ai*ai*ai*omy*2./5.*lnah;
      Loz += -8.*PI*mu*ai*ai*ai*omz*2./5.*lnah;

      ah = 1./h - 1./eps;
      parts[i].iFx += (bc.uW == DIRICHLET) * (Fnx + Ftx);
      parts[i].iFy += (bc.uW == DIRICHLET) * (Fny + Fty);
      parts[i].iFz += (bc.uW == DIRICHLET) * (Fnz + Ftz);
      parts[i].iLx += (bc.uW == DIRICHLET) * Lox;
      parts[i].iLy += (bc.uW == DIRICHLET) * Loy;
      parts[i].iLz += (bc.uW == DIRICHLET) * Loz;
      if(h < eps) { // for short-range repulsive force, replace normal component
        parts[i].iFx -= (bc.uW == DIRICHLET) * Fnx;
        parts[i].iFx += (bc.uW == DIRICHLET) * parts[i].rho*4./3.*PI*ai*ai*ai*ah;
      }
    }

    // east wall
    dx = fabs(parts[i].x - dom->xe);
    h = dx - ai;
    if(h < hN) {
      ah = ai/h - ai/hN;
      lnah = log(ai/h);

      Un = parts[i].u - bc.uED;
      Utx = 0.;
      Uty = parts[i].v - bc.vED;
      Utz = parts[i].w - bc.wED;
      omx = parts[i].ox;
      omy = parts[i].oy;
      omz = parts[i].oz;

      Fnx = -6.*PI*mu*ai*Un*ah;
      Fny = 0.;
      Fnz = 0.;

      Ftx = 0.;
      Fty = -6.*PI*mu*ai*Uty*8./15.*lnah;
      Ftz = -6.*PI*mu*ai*Utz*8./15.*lnah;
      Ftx += 0.;
      Fty += -8.*PI*mu*ai*ai*omz*1./10.*lnah;
      Ftz += 8.*PI*mu*ai*ai*omy*1./10.*lnah;

      Lox = 0.;
      Loy = 8.*PI*mu*ai*ai*Utz*1./10.*lnah;
      Loz = -8.*PI*mu*ai*ai*Uty*1./10.*lnah;
      Lox += 0.;
      Loy += -8.*PI*mu*ai*ai*ai*omy*2./5.*lnah;
      Loz += -8.*PI*mu*ai*ai*ai*omz*2./5.*lnah;

      ah = 1./h - 1./eps;
      parts[i].iFx += (bc.uE == DIRICHLET) * (Fnx + Ftx);
      parts[i].iFy += (bc.uE == DIRICHLET) * (Fny + Fty);
      parts[i].iFz += (bc.uE == DIRICHLET) * (Fnz + Ftz);
      parts[i].iLx += (bc.uE == DIRICHLET) * Lox;
      parts[i].iLy += (bc.uE == DIRICHLET) * Loy;
      parts[i].iLz += (bc.uE == DIRICHLET) * Loz;
      if(h < eps) { // for short-range repulsive force, replace normal component
        parts[i].iFx -= (bc.uE == DIRICHLET) * Fnx;
        parts[i].iFx -= (bc.uE == DIRICHLET) * parts[i].rho*4./3.*PI*ai*ai*ai*ah;
      }
    } 

    // south wall
    dy = fabs(parts[i].y - dom->ys);
    h = dy - ai;
    if(h < hN) {
      ah = ai/h - ai/hN;
      lnah = log(ai/h);

      Un = parts[i].v - bc.vSD;
      Utx = parts[i].u - bc.uSD;
      Uty = 0.;
      Utz = parts[i].w - bc.wSD;
      omx = parts[i].ox;
      omy = parts[i].oy;
      omz = parts[i].oz;

      Fnx = 0.;
      Fny = -6.*PI*mu*ai*Un*ah;
      Fnz = 0.;

      Ftx = -6.*PI*mu*ai*Utx*8./15.*lnah;
      Fty = 0.;
      Ftz = -6.*PI*mu*ai*Utz*8./15.*lnah;
      Ftx += -8.*PI*mu*ai*ai*omz*1./10.*lnah;
      Fty += 0.;
      Ftz += 8.*PI*mu*ai*ai*omx*1./10.*lnah;

      Lox = 8.*PI*mu*ai*ai*Utz*1./10.*lnah;
      Loy = 0.;
      Loz = -8.*PI*mu*ai*ai*Utx*1./10.*lnah;
      Lox += -8.*PI*mu*ai*ai*ai*omx*2./5.*lnah;
      Loy += 0.;
      Loz += -8.*PI*mu*ai*ai*ai*omz*2./5.*lnah;

      ah = 1./h - 1./eps;
      parts[i].iFx += (bc.vS == DIRICHLET) * (Fnx + Ftx);
      parts[i].iFy += (bc.vS == DIRICHLET) * (Fny + Fty);
      parts[i].iFz += (bc.vS == DIRICHLET) * (Fnz + Ftz);
      parts[i].iLx += (bc.vS == DIRICHLET) * Lox;
      parts[i].iLy += (bc.vS == DIRICHLET) * Loy;
      parts[i].iLz += (bc.vS == DIRICHLET) * Loz;
      if(h < eps) { // for short-range repulsive force, replace normal component
        parts[i].iFy -= (bc.vS == DIRICHLET) * Fny;
        parts[i].iFy += (bc.vS == DIRICHLET) * parts[i].rho*4./3.*PI*ai*ai*ai*ah;
      }
    }

    // north wall
    dy = fabs(parts[i].y - dom->ye);
    h = dy - ai;
    if(h < hN) {
      ah = ai/h - ai/hN;
      lnah = log(ai/h);

      Un = parts[i].v - bc.vND;
      Utx = parts[i].u - bc.uND;
      Uty = 0.;
      Utz = parts[i].w - bc.wND;
      omx = parts[i].ox;
      omy = parts[i].oy;
      omz = parts[i].oz;

      Fnx = 0.;
      Fny = -6.*PI*mu*ai*Un*ah;
      Fnz = 0.;

      Ftx = -6.*PI*mu*ai*Utx*8./15.*lnah;
      Fty = 0.;
      Ftz = -6.*PI*mu*ai*Utz*8./15.*lnah;
      Ftx += 8.*PI*mu*ai*ai*omz*1./10.*lnah;
      Fty += 0.;
      Ftz += -8.*PI*mu*ai*ai*omx*1./10.*lnah;

      Lox = -8.*PI*mu*ai*ai*Utz*1./10.*lnah;
      Loy = 0.;
      Loz = 8.*PI*mu*ai*ai*Utx*1./10.*lnah;
      Lox += -8.*PI*mu*ai*ai*ai*omx*2./5.*lnah;
      Loy += 0.;
      Loz += -8.*PI*mu*ai*ai*ai*omz*2./5.*lnah;

      ah = 1./h - 1./eps;
      parts[i].iFx += (bc.vN == DIRICHLET) * (Fnx + Ftx);
      parts[i].iFy += (bc.vN == DIRICHLET) * (Fny + Fty);
      parts[i].iFz += (bc.vN == DIRICHLET) * (Fnz + Ftz);
      parts[i].iLx += (bc.vN == DIRICHLET) * Lox;
      parts[i].iLy += (bc.vN == DIRICHLET) * Loy;
      parts[i].iLz += (bc.vN == DIRICHLET) * Loz;
      if(h < eps) { // for short-range repulsive force, replace normal component
        parts[i].iFy -= (bc.vN == DIRICHLET) * Fny;
        parts[i].iFy -= (bc.vN == DIRICHLET) * parts[i].rho*4./3.*PI*ai*ai*ai*ah;
      }
    }

    // bottom wall
    dz = fabs(parts[i].z - dom->zs);
    h = dz - ai;
    if(h < hN) {
      ah = ai/h - ai/hN;
      lnah = log(ai/h);

      Un = parts[i].w - bc.wBD;
      Utx = parts[i].u - bc.uBD;
      Uty = parts[i].v - bc.vBD;
      Utz = 0.;
      omx = parts[i].ox;
      omy = parts[i].oy;
      omz = parts[i].oz;

      Fnx = 0.;
      Fny = 0.;
      Fnz = -6.*PI*mu*ai*Un*ah;

      Ftx = -6.*PI*mu*ai*Utx*8./15.*lnah;
      Fty = -6.*PI*mu*ai*Uty*8./15.*lnah;
      Ftz = 0.;
      Ftx += 8.*PI*mu*ai*ai*omy*1./10.*lnah;
      Fty += -8.*PI*mu*ai*ai*omx*1./10.*lnah;
      Ftz += 0.;

      Lox = -8.*PI*mu*ai*ai*Uty*1./10.*lnah;
      Loy = 8.*PI*mu*ai*ai*Utx*1./10.*lnah;
      Loz = 0.;
      Lox += -8.*PI*mu*ai*ai*ai*omx*2./5.*lnah;
      Loy += -8.*PI*mu*ai*ai*ai*omy*2./5.*lnah;
      Loz += 0.;

      ah = 1./h - 1./eps;
      parts[i].iFx += (bc.wB == DIRICHLET) * (Fnx + Ftx);
      parts[i].iFy += (bc.wB == DIRICHLET) * (Fny + Fty);
      parts[i].iFz += (bc.wB == DIRICHLET) * (Fnz + Ftz);
      parts[i].iLx += (bc.wB == DIRICHLET) * Lox;
      parts[i].iLy += (bc.wB == DIRICHLET) * Loy;
      parts[i].iLz += (bc.wB == DIRICHLET) * Loz;
      if(h < eps) { // for short-range repulsive force, replace normal component
        parts[i].iFz -= (bc.wB == DIRICHLET) * Fnz;
        parts[i].iFz += (bc.wB == DIRICHLET) * parts[i].rho*4./3.*PI*ai*ai*ai*ah;
      }
    }

    // top wall
    dz = fabs(parts[i].z - dom->ze);
    h = dz - ai;
    if(h < hN) {
      ah = ai/h - ai/hN;
      lnah = log(ai/h);

      Un = parts[i].w - bc.wTD;
      Utx = parts[i].u - bc.uTD;
      Uty = parts[i].v - bc.vTD;
      Utz = 0.;
      omx = parts[i].ox;
      omy = parts[i].oy;
      omz = parts[i].oz;

      Fnx = 0.;
      Fny = 0.;
      Fnz = -6.*PI*mu*ai*Un*ah;

      Ftx = -6.*PI*mu*ai*Utx*8./15.*lnah;
      Fty = -6.*PI*mu*ai*Uty*8./15.*lnah;
      Ftz = 0.;
      Ftx += -8.*PI*mu*ai*ai*omy*1./10.*lnah;
      Fty += 8.*PI*mu*ai*ai*omx*1./10.*lnah;
      Ftz += 0.;

      Lox = 8.*PI*mu*ai*ai*Uty*1./10.*lnah;
      Loy = -8.*PI*mu*ai*ai*Utx*1./10.*lnah;
      Loz = 0.;
      Lox += -8.*PI*mu*ai*ai*ai*omx*2./5.*lnah;
      Loy += -8.*PI*mu*ai*ai*ai*omy*2./5.*lnah;
      Loz += 0.;

      ah = 1./h - 1./eps;
      parts[i].iFx += (bc.wT == DIRICHLET) * (Fnx + Ftx);
      parts[i].iFy += (bc.wT == DIRICHLET) * (Fny + Fty);
      parts[i].iFz += (bc.wT == DIRICHLET) * (Fnz + Ftz);
      parts[i].iLx += (bc.wT == DIRICHLET) * Lox;
      parts[i].iLy += (bc.wT == DIRICHLET) * Loy;
      parts[i].iLz += (bc.wT == DIRICHLET) * Loz;
      if(h < eps) { // for short-range repulsive force, replace normal component
        parts[i].iFz -= (bc.wT == DIRICHLET) * Fnz;
        parts[i].iFz -= (bc.wT == DIRICHLET) * parts[i].rho*4./3.*PI*ai*ai*ai*ah;
      }
    }
   }
}

__global__ void spring_parts(part_struct *parts, int nparts)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;

  if(i < nparts) {
    real dx = parts[i].x-parts[i].spring_x;
    real dy = parts[i].y-parts[i].spring_y;
    real dz = parts[i].z-parts[i].spring_z;

    parts[i].kFx = - parts[i].spring_k * dx;
    parts[i].kFy = - parts[i].spring_k * dy;
    parts[i].kFz = - parts[i].spring_k * dz;
  }
}
*/
__global__ void yank_u_WE(real *u, dom_struct *dom, real *plane, real xpos,
  real vel)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  real ddx = 1. / dom->dx;

  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int i = floor((xpos - dom->xs) * ddx) + DOM_BUF;
    if(i < dom->Gfx.is) i += dom->Gfx.inb;
    if(i > dom->Gfx.ie-1) i -= dom->Gfx.inb;
    real xx = (i-DOM_BUF) * dom->dx + dom->xs;
    int W = i + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;
    int E = (i+1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b;
    real dudx = (u[E] - u[W]) * ddx;

    plane[tj + tk*dom->Gfx.jnb] = u[W] + dudx * (xpos - xx) + vel;
  }
}

__global__ void yank_v_WE(real *v, dom_struct *dom, real *plane, real xpos,
  real vel)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  real ddx = 1. / dom->dx;

  if((tj < dom->Gfy._jnb) && (tk < dom->Gfy._knb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int i = floor((xpos - dom->xs) * ddx - 0.5) + DOM_BUF;
    if(i < dom->Gfy.is) i += dom->Gfy.inb;
    if(i > dom->Gfy.ie-1) i -= dom->Gfy.inb;
    real xx = (i-DOM_BUF+0.5) * dom->dx + dom->xs;
    int W = i + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b;
    int E = (i+1) + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b;
    real dvdx = (v[E] - v[W]) * ddx;

    plane[tj + tk*dom->Gfy.jnb] = v[W] + dvdx * (xpos - xx);
  }
}

__global__ void yank_w_WE(real *w, dom_struct *dom, real *plane, real xpos,
  real vel)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  real ddx = 1. / dom->dx;

  if((tj < dom->Gfz._jnb) && (tk < dom->Gfz._knb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int i = floor((xpos - dom->xs) * ddx - 0.5) + DOM_BUF;
    if(i < dom->Gfz.is) i += dom->Gfz.inb;
    if(i > dom->Gfz.ie-1) i -= dom->Gfz.inb;
    real xx = (i-DOM_BUF + 0.5) * dom->dx + dom->xs;
    int W = i + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b;
    int E = (i+1) + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b;
    real dwdx = (w[E] - w[W]) * ddx;

    plane[tj + tk*dom->Gfz.jnb] = w[W] + dwdx * (xpos - xx);
  }
}

__global__ void yank_u_SN(real *u, dom_struct *dom, real *plane, real ypos,
  real vel)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  real ddy = 1. / dom->dy;

  if((tk < dom->Gfx._inb) && (ti < dom->Gfx._inb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int j = floor((ypos - dom->ys) * ddy - 0.5) + DOM_BUF;
    if(j < dom->Gfx.js) j += dom->Gfx.jnb;
    if(j > dom->Gfx.je-1) j -= dom->Gfx.jnb;
    real yy = (j-DOM_BUF + 0.5) * dom->dy + dom->ys;
    int S = ti + j*dom->Gfx.s1b + tk*dom->Gfx.s2b;
    int N = ti + (j+1)*dom->Gfx.s1b + tk*dom->Gfx.s2b;
    real dudy = (u[N] - u[S]) * ddy;

    plane[tk + ti*dom->Gfx.knb] = u[S] + dudy * (ypos - yy);
  }
}

__global__ void yank_v_SN(real *v, dom_struct *dom, real *plane, real ypos,
  real vel)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  real ddy = 1. / dom->dy;

  if((ti < dom->Gfy._inb) && (tk < dom->Gfy._knb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int j = floor((ypos - dom->ys) * ddy) + DOM_BUF;
    if(j < dom->Gfy.js) j += dom->Gfy.jnb;
    if(j > dom->Gfy.je-1) j -= dom->Gfy.jnb;
    real yy = (j-DOM_BUF) * dom->dy + dom->ys;
    int S = ti + j*dom->Gfy.s1b + tk*dom->Gfy.s2b;
    int N = ti + (j+1)*dom->Gfy.s1b + tk*dom->Gfy.s2b;
    real dvdy = (v[N] - v[S]) * ddy;

    plane[tk + ti*dom->Gfy.knb] = v[S] + dvdy * (ypos - yy) + vel;
  }
}

__global__ void yank_w_SN(real *w, dom_struct *dom, real *plane, real ypos,
  real vel)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  real ddy = 1. / dom->dy;

  if((ti < dom->Gfz._inb) && (tk < dom->Gfz._knb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int j = floor((ypos - dom->ys) * ddy - 0.5) + DOM_BUF;
    if(j < dom->Gfz.js) j += dom->Gfz.jnb;
    if(j > dom->Gfz.je-1) j -= dom->Gfz.jnb;
    real yy = (j-DOM_BUF + 0.5) * dom->dy + dom->ys;
    int S = ti + j*dom->Gfz.s1b + tk*dom->Gfz.s2b;
    int N = ti + (j+1)*dom->Gfz.s1b + tk*dom->Gfz.s2b;
    real dwdy = (w[N] - w[S]) * ddy;

    plane[tk + ti*dom->Gfz.knb] = w[S] + dwdy * (ypos - yy);
  }
}

__global__ void yank_u_BT(real *u, dom_struct *dom, real *plane, real zpos,
  real vel)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  real ddz = 1. / dom->dz;

  if((ti < dom->Gfx._inb) && (tj < dom->Gfx._jnb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int k = floor((zpos - dom->zs) * ddz - 0.5) + DOM_BUF;
    if(k < dom->Gfx.ks) k += dom->Gfx.knb;
    if(k > dom->Gfx.ke-1) k -= dom->Gfx.knb;
    real zz = (k-DOM_BUF + 0.5) * dom->dz + dom->zs;
    int B = ti + tj*dom->Gfx.s1b + k*dom->Gfx.s2b;
    int T = ti + tj*dom->Gfx.s1b + (k+1)*dom->Gfx.s2b;
    real dudz = (u[T] - u[B]) * ddz;

    plane[ti + tj*dom->Gfx.inb] = u[B] + dudz * (zpos - zz);
  }
}

__global__ void yank_v_BT(real *v, dom_struct *dom, real *plane, real zpos,
  real vel)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  real ddz = 1. / dom->dz;

  if((ti < dom->Gfy._inb) && (tj < dom->Gfy._jnb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int k = floor((zpos - dom->zs) * ddz - 0.5) + DOM_BUF;
    if(k < dom->Gfy.ks) k += dom->Gfy.knb;
    if(k > dom->Gfy.ke-1) k -= dom->Gfy.knb;
    real zz = (k-DOM_BUF + 0.5) * dom->dz + dom->zs;
    int B = ti + tj*dom->Gfy.s1b + k*dom->Gfy.s2b;
    int T = ti + tj*dom->Gfy.s1b + (k+1)*dom->Gfy.s2b;
    real dvdz = (v[T] - v[B]) * ddz;

    plane[ti + tj*dom->Gfy.inb] = v[B] + dvdz * (zpos - zz);
  }
}

__global__ void yank_w_BT(real *w, dom_struct *dom, real *plane, real zpos,
  real vel)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  real ddz = 1. / dom->dx;

  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    // find index of node
    // for now, ignore motion tangential to plane
    int k = floor((zpos - dom->zs) * ddz) + DOM_BUF;
    if(k < dom->Gfz.ks) k += dom->Gfz.knb;
    if(k > dom->Gfz.ke-1) k -= dom->Gfz.knb;
    real zz = (k-DOM_BUF) * dom->dz + dom->zs;
    int B = ti + tj*dom->Gfz.s1b + k*dom->Gfz.s2b;
    int T = ti + tj*dom->Gfz.s1b + (k+1)*dom->Gfz.s2b;
    real dwdz = (w[T] - w[B]) * ddz;

    plane[ti + tj*dom->Gfz.inb] = w[B] + dwdz * (zpos - zz) + vel;
  }
}

__global__ void colocate_Gfx(real *u, real *u_co, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x + DOM_BUF;
  int tk = blockDim.y*blockIdx.y + threadIdx.y + DOM_BUF;

  if((tj < dom->Gfx.jnb-1) && (tk < dom->Gfx.knb-1)) {
    for(int i = dom->Gfx.is; i < dom->Gfx.ie-1; i++) {
      u_co[(i-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2] =
        0.5 * (u[i + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b]
        + u[(i+1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b]);
    }
  }
}

__global__ void colocate_Gfy(real *v, real *v_co, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x + DOM_BUF;
  int ti = blockDim.y*blockIdx.y + threadIdx.y + DOM_BUF;

  if((tk < dom->Gfy.knb-1) && (ti < dom->Gfy.inb-1)) {
    for(int j = dom->Gfy.js; j < dom->Gfy.je-1; j++) {
      v_co[(ti-DOM_BUF) + (j-DOM_BUF)*dom->Gcc.s1 + (tk-DOM_BUF)*dom->Gcc.s2] =
        0.5 * (v[ti + j*dom->Gfy.s1b + tk*dom->Gfy.s2b]
        + v[ti + (j+1)*dom->Gfy.s1b + tk*dom->Gfy.s2b]);
    }
  }
}

__global__ void colocate_Gfz(real *w, real *w_co, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x + DOM_BUF;
  int tj = blockDim.y*blockIdx.y + threadIdx.y + DOM_BUF;

  if((ti < dom->Gfz.inb-1) && (tj < dom->Gfz.jnb-1)) {
    for(int k = dom->Gfz.ks; k < dom->Gfz.ke-1; k++) {
      w_co[(ti-DOM_BUF) + (tj-DOM_BUF)*dom->Gcc.s1 + (k-DOM_BUF)*dom->Gcc.s2] =
        0.5 * (w[ti + tj*dom->Gfz.s1b + k*dom->Gfz.s2b]
        + w[ti + tj*dom->Gfz.s1b + (k+1)*dom->Gfz.s2b]);
    }
  }
}

__global__ void energy_multiply(real *u_co, real *v_co, real *w_co, real *co,
  dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  int C;  // memory location

  if((tj < dom->Gcc.jn) && (tk < dom->Gcc.kn)) {
    for(int i = dom->Gcc.is-DOM_BUF; i < dom->Gcc.ie-DOM_BUF; i++) {
      C = i + tj*dom->Gcc.s1 + tk*dom->Gcc.s2;
      co[C] = u_co[C]*u_co[C] + v_co[C]*v_co[C] + w_co[C]*w_co[C];
    }
  }
}
