/****h* Bluebottle/cuda_bluebottle_kernel
 * NAME
 *  cuda_bluebottle_kernel
 * FUNCTION
 *  Bluebottle CUDA kernel functions.
 ******
 */

#ifndef _CUDA_SCALAR_H
#define _CUDA_SCALAR_H

extern "C"
{
#include "bluebottle.h"
#include "point.h"
}

extern texture<float,1,cudaReadModeElementType> texRefGaussian;
extern texture<int,1,cudaReadModeElementType> texRefDomInfo;



__global__ void copy_sc_noghost(real *sc_noghost, real *sc_ghost, dom_struct *dom);
__global__ void scalar_coeffs_init(dom_struct *dom, int pitch, real *values);
//__global__ void scalar_coeffs(real DIFF, real dt, dom_struct *dom, int pitch,real *values);
__global__ void scalar_coeffs(real DIFF, real dt, dom_struct *dom, int pitch,
  real *values,  int *flag_u, int *flag_v, int *flag_w,real *epsp);
__global__ void copy_sc_ghost(real *sc_ghost, real *sc_noghost, dom_struct *dom);

__global__ void scalar_rhs_FTCS(real rho_f, real DIFF, real *u, real *v, real *w,
 real *epsp, real *f, real *conv0,real *conv,real *diff, real *sc0,real *sc_rhs, dom_struct *dom, real dt, real dt0);

__global__ void scalar_rhs_upwind_1st(real rho_f, real DIFF, 
real *u, real *v, real *w, 
 real *epsp,real *epsp0,
 real *f,real *f0,
 real *conv0,real *conv,real *diff,
 real *sc0, real *sc_rhs, dom_struct *dom,
 real dt, real dt0);







__global__ void diffScalar_coeffs_init(dom_struct *dom, int pitch, real *values,int coordiSys);
//This kernel only works for diffusion of a scalar kernel.
//This coeff matrix depend on the flow boundary condition rather than the diffScalar BC!!!
//Diff_dt is the length scale square; and it is equall to (delta_f^2-dx^2)/16/ln2
 
__global__ void copy_diffSc_noghost(real *sc_noghost, real *sc_ghost, dom_struct *dom,int coordiSys);

//__global__ void diffScalar_rhs_CN(real DIFF_dt, real *sc0, real *sc_rhs, dom_struct *dom, real theta);
__global__ void diffScalar_rhs_CN(real DIFF_dt, real *sc0, real *sc_rhs, dom_struct *dom);
//Diff_dt is the length scale square; and it is equall to (delta_f^2-dx^2)/16/ln2
__global__ void diffScalar_coeffs_CN(real thetaD, dom_struct *dom, int pitch,
  real *values,  int *flag_u, int *flag_v, int *flag_w,int coordiSys);



//Using central difference scheme in space, and adam-bashforth scheme in time.
__global__ void diffScalar_explicitD(real *sc,real *sc0,
				dom_struct *dom, real DIFF_dt, int coordiSys);


__global__ void diff_Gcc_explicitD(real *sc,real *sc0,
dom_struct *dom, real DIFF_dt);
//Using central difference scheme in space, and adam-bashforth scheme in time.
__global__ void diff_Gfx_explicitD(real *sc,real *sc0,
dom_struct *dom, real DIFF_dt);
__global__ void diff_Gfy_explicitD(real *sc,real *sc0,
dom_struct *dom, real DIFF_dt);
__global__ void diff_Gfz_explicitD(real *sc,real *sc0,
dom_struct *dom, real DIFF_dt);



__global__ void scalar_coeffs_periodic_W(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values, real *epsp,int *flag);
__global__ void scalar_coeffs_periodic_E(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values, real *epsp,int *flag);
__global__ void scalar_coeffs_periodic_S(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values, real *epsp,int *flag);
__global__ void scalar_coeffs_periodic_N(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values, real *epsp,int *flag);
__global__ void scalar_coeffs_periodic_B(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values, real *epsp,int *flag);
__global__ void scalar_coeffs_periodic_T(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values, real *epsp,int *flag);


//For Diffusion of source after mollification
__global__ void diffScalar_coeffs_periodic_W(real DIFF_dt, dom_struct *dom,
  int pitch, real *values, int coordiSys);
__global__ void diffScalar_coeffs_periodic_E(real DIFF_dt, dom_struct *dom,
  int pitch, real *values, int coordiSys);
__global__ void diffScalar_coeffs_periodic_S(real DIFF_dt, dom_struct *dom,
  int pitch, real *values, int coordiSys);
__global__ void diffScalar_coeffs_periodic_N(real DIFF_dt, dom_struct *dom,
  int pitch, real *values, int coordiSys);
__global__ void diffScalar_coeffs_periodic_B(real DIFF_dt, dom_struct *dom,
  int pitch, real *values, int coordiSys);
__global__ void diffScalar_coeffs_periodic_T(real DIFF_dt, dom_struct *dom,
  int pitch, real *values, int coordiSys);



/*
__global__ void scalar_coeffs_periodic_W(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values);
__global__ void scalar_coeffs_periodic_E(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values);
__global__ void scalar_coeffs_periodic_S(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values);
__global__ void scalar_coeffs_periodic_N(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values);
__global__ void scalar_coeffs_periodic_B(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values);
__global__ void scalar_coeffs_periodic_T(real DIFF, real dt, dom_struct *dom,
  int pitch, real *values);
*/





__global__ void BC_scSrc_W_P(int coordiSys, real *scSrc, dom_struct *dom);
__global__ void BC_scSrc_W_N(int coordiSys, real *scSrc, dom_struct *dom);
__global__ void BC_scSrc_E_P(int coordiSys, real *scSrc, dom_struct *dom);
__global__ void BC_scSrc_E_N(int coordiSys, real *scSrc, dom_struct *dom);
__global__ void BC_scSrc_S_P(int coordiSys, real *scSrc, dom_struct *dom);
__global__ void BC_scSrc_S_N(int coordiSys, real *scSrc, dom_struct *dom);
__global__ void BC_scSrc_N_P(int coordiSys, real *scSrc, dom_struct *dom);
__global__ void BC_scSrc_N_N(int coordiSys, real *scSrc, dom_struct *dom);
__global__ void BC_scSrc_B_P(int coordiSys, real *scSrc, dom_struct *dom);
__global__ void BC_scSrc_B_N(int coordiSys, real *scSrc, dom_struct *dom);
__global__ void BC_scSrc_T_P(int coordiSys, real *scSrc, dom_struct *dom);
__global__ void BC_scSrc_T_N(int coordiSys, real *scSrc, dom_struct *dom);


__global__ void BC_scSrc_W_D(int coordiSys, real *scSrc, dom_struct *dom, real bc);
__global__ void BC_scSrc_E_D(int coordiSys, real *scSrc, dom_struct *dom, real bc);
__global__ void BC_scSrc_S_D(int coordiSys, real *scSrc, dom_struct *dom, real bc);
__global__ void BC_scSrc_N_D(int coordiSys, real *scSrc, dom_struct *dom, real bc);
__global__ void BC_scSrc_B_D(int coordiSys, real *scSrc, dom_struct *dom, real bc);
__global__ void BC_scSrc_T_D(int coordiSys, real *scSrc, dom_struct *dom, real bc);



/****f* cuda_bluebottle_kernel/BC_sc_W_P<<<>>>()
 * NAME
 *  BC_sc_W_P<<<>>>()
 * USAGE
 */
__global__ void BC_sc_W_P(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the west face scalar field.
 * ARGUMENTS
 *  * p -- the device scalar field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_sc_W_N<<<>>>()
 * NAME
 *  BC_sc_W_N<<<>>>()
 * USAGE
 */
__global__ void BC_sc_W_N(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the west face scalar field.
 * ARGUMENTS
 *  * p -- the device scalar field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_sc_E_P<<<>>>()
 * NAME
 *  BC_sc_E_P<<<>>>()
 * USAGE
 */
__global__ void BC_sc_E_P(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the east face scalar field.
 * ARGUMENTS
 *  * p -- the device scalar field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_sc_E_N<<<>>>()
 * NAME
 *  BC_sc_E_N<<<>>>()
 * USAGE
 */
__global__ void BC_sc_E_N(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the east face scalar field.
 * ARGUMENTS
 *  * p -- the device scalar field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_sc_S_P<<<>>>()
 * NAME
 *  BC_sc_S_P<<<>>>()
 * USAGE
 */
__global__ void BC_sc_S_P(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the south face scalar field.
 * ARGUMENTS
 *  * p -- the device scalar field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_sc_S_N<<<>>>()
 * NAME
 *  BC_sc_S_N<<<>>>()
 * USAGE
 */
__global__ void BC_sc_S_N(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the south face scalar field.
 * ARGUMENTS
 *  * p -- the device scalar field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_sc_N_P<<<>>>()
 * NAME
 *  BC_sc_N_P<<<>>>()
 * USAGE
 */
__global__ void BC_sc_N_P(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the north face scalar field.
 * ARGUMENTS
 *  * p -- the device scalar field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_sc_N_N<<<>>>()
 * NAME
 *  BC_sc_N_N<<<>>>()
 * USAGE
 */
__global__ void BC_sc_N_N(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the north face scalar field.
 * ARGUMENTS
 *  * p -- the device scalar field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_sc_B_P<<<>>>()
 * NAME
 *  BC_sc_B_P<<<>>>()
 * USAGE
 */
__global__ void BC_sc_B_P(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the bottom face scalar field.
 * ARGUMENTS
 *  * p -- the device scalar field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_sc_B_N<<<>>>()
 * NAME
 *  BC_sc_B_N<<<>>>()
 * USAGE
 */
__global__ void BC_sc_B_N(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the bottom face scalar field.
 * ARGUMENTS
 *  * p -- the device scalar field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_sc_T_P<<<>>>()
 * NAME
 *  BC_sc_T_P<<<>>>()
 * USAGE
 */
__global__ void BC_sc_T_P(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the top face scalar field.
 * ARGUMENTS
 *  * p -- the device scalar field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_sc_T_N<<<>>>()
 * NAME
 *  BC_sc_T_N<<<>>>()
 * USAGE
 */
__global__ void BC_sc_T_N(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the top face scalar field.
 * ARGUMENTS
 *  * p -- the device scalar field subdomain
 *  * dom -- the device subdomain
 ******
 */


/****f* cuda_bluebottle_kernel/BC_sc_W_D<<<>>>()
 * NAME
 *  BC_sc_W_D<<<>>>()
 * USAGE
 */
__global__ void BC_sc_W_D(real *sc, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the west face sc-scalar field.
 * ARGUMENTS
 *  * sc -- the device sc-scalar field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */



/****f* cuda_bluebottle_kernel/BC_sc_E_D<<<>>>()
 * NAME
 *  BC_sc_E_D<<<>>>()
 * USAGE
 */
__global__ void BC_sc_E_D(real *sc, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the east face sc-scalar field.
 * ARGUMENTS
 *  * sc -- the device sc-scalar field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */


/****f* cuda_bluebottle_kernel/BC_sc_S_D<<<>>>()
 * NAME
 *  BC_sc_S_D<<<>>>()
 * USAGE
 */
__global__ void BC_sc_S_D(real *sc, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the south face sc-scalar field.
 * ARGUMENTS
 *  * sc -- the device sc-scalar field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */


/****f* cuda_bluebottle_kernel/BC_sc_N_D<<<>>>()
 * NAME
 *  BC_sc_N_D<<<>>>()
 * USAGE
 */
__global__ void BC_sc_N_D(real *sc, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the north face sc-scalar field.
 * ARGUMENTS
 *  * sc -- the device sc-scalar field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */



/****f* cuda_bluebottle_kernel/BC_sc_T_D<<<>>>()
 * NAME
 *  BC_sc_T_D<<<>>>()
 * USAGE
 */
__global__ void BC_sc_T_D(real *sc, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the top face sc-scalar field.
 * ARGUMENTS
 *  * sc -- the device sc-scalar field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */


/****f* cuda_bluebottle_kernel/BC_sc_B_D<<<>>>()
 * NAME
 *  BC_sc_B_D<<<>>>()
 * USAGE
 */
__global__ void BC_sc_B_D(real *sc, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the bottom face sc-scalar field.
 * ARGUMENTS
 *  * sc -- the device sc-scalar field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */





__global__ void advance_sc_upwind_1st_init(real DIFF, real *u, real *v, real *w,real *f,real *epsp,
  real *diff0, real *conv0, real *diff, real *conv, real *sc,real *sc0,
  dom_struct *dom, real dt0, real dt);


//first order upwind
__global__ void advance_sc_upwind_1st(real DIFF, real *u, real *v, real *w,real *f,real *epsp,
  real *diff0, real *conv0, real *diff, real *conv, real *sc,real *sc0,
  dom_struct *dom, real dt0, real dt);






/*
Advance the scalar based on fluid velocity and particle volume fraction, using Adam-bashforth method in time and
central difference in space
	DIFF is scalar diffusion coeffcient,
	u~v are the face-center velocity of flow, f is the body force on flow
	epsp is the  particle volume fraction at each grid cell
	diff,conv are the diffusive and convective term of the scalar.  diff0,conv0 are for old step
	sc is the scalar field, sc0 for old step. 
	dt,dt0 are time step of scalar diffusion equation
*/
__global__ void advance_sc(real DIFF, real *u, real *v, real *w,real *f,real *epsp,
  real *diff0, real *conv0, real *diff, real *conv,real *sc,real *sc0,
  dom_struct *dom, real dt0, real dt);

//at the very 1st step to advance scalar without Adam-bashforth interpolation
__global__ void advance_sc_init(real DIFF, real *u, real *v, real *w,real *f,real *epsp,
  real *diff0, real *conv0, real *diff, real *conv,real *sc,real *sc0,
  dom_struct *dom, real dt0, real dt);

//Advance scalar using MacCormack scheme.
__global__ void advance_sc_macCormack(	real DIFF, 
				real *u, real *v, real *w,
				real *f,real *epsp,
				real *diff, real *conv, 
				real *sc,real *sc0,
				dom_struct *dom, real dt);

__global__ void advance_sc_QUICK(real DIFF, 
				real *u, real *v, real *w,
				real *f,real *epsp,
				real *diff, real *conv, 
				real *sc,real *sc0,
				dom_struct *dom, real dt);

//initialize flow field array to be 0 on device, include scalar source and particle volume fraction, array length is dom->s3b
__global__ void lpt_scalar_source_init(real *A, real *B,dom_struct *dom);

//initialize particle property array on device, include volume of point particle and source contributed by point particle, array length is equal to particle number
__global__ void   lpt_point_source_init(int npoints,point_struct *points,dom_struct *dom,real *volPoint,real *scSrc);

//clip volume fraction of particle in each grid cell, if it's bigger than a threshold(default:0.8), then clip it
__global__ void lpt_epsp_clip(real *epsp,dom_struct *dom);


//lpt_source_test
__global__ void lpt_scalar_source_init_test(real *src,dom_struct *dom,real t, real DIFF);

__global__ void lpt_scalar_source_convDiff_test(real *src,dom_struct *dom,real t, real DIFF);

#endif
