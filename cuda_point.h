/****h* Bluebottle/cuda_point_point_particle_kernel
 * NAME
 *  cuda_bluebottle_kernel
 * FUNCTION
 *  Bluebottle CUDA point_point_particle kernel functions.
 ******
 */

//#include "bluebottle.h"
#ifndef _CUDA_POINT_H
#define _CUDA_POINT_H

extern "C"
{
#include "bluebottle.h"
//#include "point.h"
}


//initialize array A with length n, initialize the array to be a constant C 
__global__ void array_init(real *A,dom_struct *dom, int n, real C);

/*
The following two subroutines are for mollifying point source with a Gaussian kernel
	A is located at cell center, Ksi is the mollify coefficients, 
	dirc indicate whether mollify on cell center or face center. dirc=0 is cell center, dirc =1~3 are x-face to z-face center.  When on face center, the subroutine mollifies momentum source
	dirc2 is for cell center only, indicate what property to mollify.  dirc2=0 is for mollification of scalar source, dirc2=1 is  for particle volume fraction 'epsp '
*/
__global__ void lpt_mollify_sc(int npoints, real *A, point_struct *points, dom_struct *dom,real *Ksi, int dirc, int dirc2);
__global__ void lpt_sum_ksi(int npoints, real *A, point_struct *points, dom_struct *dom,real *Ksi, int dirc,int dirc2);


//Make the position (x,y,z) periodic
__device__ void periodic_grid_position(real &x,real &y,real &z,dom_struct *dom);


//Make the index (i,j,k) periodic
__device__ void periodic_grid_index(int &ic,int &jc,int &kc,dom_struct *dom, int dirc);

/*
Gaussian kernel with object position index (ic,jc,kc) and Gaussian center position (xp,yp,zp), basically calculate their distance substitute into Gaussian kernel. 
	dx~dz are grid cell size, xs~zs are start position of the whole fluid domain
	dirc =0 indicate cell-center; dirc =1,2,3 indicate x,y,z-face-center;
*/
__device__ real lpt_integrate_mol(int ic,int jc,int kc,real xp,real yp,real zp, real dx,real dy,real dz,real xs,real ys,real zs, int dirc);


//initialize the the velocity ug~wg and stress of fluid lpt_stress_u~w, which will be interpreted at the point particle position. All arrays will be initialized to 0
__global__ void point_interp_init(int npoints,point_struct *points,real *ug,real *vg,real *wg,real *lpt_stress_u,real *lpt_stress_v,real *lpt_stress_w,real *scg);


/*
Interpolate flow field property at the particle position, it's for cell-center system
	A is cell-centerred flow field data 
	Ag is the property after interpolation, located at each particle
*/
__global__ void interpolate_point_scalar_Lag2(int npoints,real *A,real *Ag, point_struct *points, dom_struct *dom);


/*
Interpolate flow field property at the particle position, it's for face-center system
	u,v,w are face-centerred values, can be flow velocity field	
	ug,vg,wg are values after the interpolation, located at the particle position
	bc is the boundary condition for particles, USE PERIODIC BC FOR THE TIME BEING
*/
__global__ void interpolate_point_vel_Lag2(real *u, real *v, real *w,int npoints,real rho_f, real nu,real *ug,real *vg,real *wg, point_struct *points, dom_struct *dom, BC bc);


/*
calculate fluid stress  \nabla \cdot \sigma =-\nabla P +\mu * nabla^2 u
	u,v,w are velocity field at face center
	flag_u~w are tags which indicate interior or boundary of the domain.  flag=Dirichlet（inside particle）, Neuman（domain boundary）, otherwise at interiol of the domain
	dt and dt0 are for Adam-bashforth method in time, not used yet
*/
__global__ void stress_u(real rho_f, real nu, real *u0,real *p,real *p0,real *stress, dom_struct *dom,int *flag_u,real dt,real dt0);

__global__ void stress_v(real rho_f, real nu, real *v0,real *p,real *p0,real *stress, dom_struct *dom,int *flag_v,real dt,real dt0);

__global__ void stress_w(real rho_f, real nu, real *w0,real *p,real *p0,real *stress, dom_struct *dom,int *flag_w,real dt,real dt0);


/*
Calculate fluid force acted on particles, include drag force,fluid stress, added mass effect. Also calculate the exchange rate of soluble scalar mass between particle and flow field.
	ug,vg,wg are the fluid velocity at the particle location
	lpt_stress_u~w are the fluid stresses at the particle location
	gradP is the body force act on fluid
	g is gravity act on particle and fluid.
	C_add,C_stress,C_drag are coefficients of added mass,fluid stress, and drag force. Default values are 0.5,1,1 correspondingly
	sc_eq is the saturation concentration of the scalar at the particle-fluid interface
	DIFF is the diffusivity of the scalar
*/
__global__ void drag_points(point_struct *points, int npoints,
real *ug,real *vg,real *wg,
real *lpt_stress_u,real *lpt_stress_v,real *lpt_stress_w,real *scg,
real rho_f,real mu, g_struct g,gradP_struct gradP,
real C_add,real C_stress,real C_drag,
real sc_eq,real DIFF);

//advance particle position based on Euler prediction correction, it can be divided into 2 steps. See "Predictor–corrector method" in wiki for ref
__global__ void move_points_a(point_struct *points, int npoints, real dt);
__global__ void move_points_b(dom_struct *dom,point_struct *points, int npoints, real dt);


//Locate the grid cell index (i,j,k) of each particle
__global__ void lpt_localize(int npoints, point_struct *points, dom_struct *dom, BC bc);







/****f* cuda_particle_kernel/reset_flag_u<<<>>>()
 * NAME
 *  reset_flag_u<<<>>>()
 * USAGE
 */
__global__ void reset_flag_u(int *flag_u, dom_struct *dom, BC bc);
/*
 * FUNCTION
 *  Set all flag_u nodes to fluid (= 1).
 * ARGUMENTS
 *  * flag_u -- the device x-direction velocity flag array subdomain
 *  * dom -- the device subdomain
 *  * bc -- the boundary conditions
 ******
 */

/****f* cuda_particle_kernel/reset_flag_v<<<>>>()
 * NAME
 *  reset_flag_v<<<>>>()
 * USAGE
 */
__global__ void reset_flag_v(int *flag_v, dom_struct *dom, BC bc);
/*
 * FUNCTION
 *  Set all flag_v nodes to fluid (= 1).
 * ARGUMENTS
 *  * flag_v -- the device y-direction velocity flag array subdomain
 *  * dom -- the device subdomain
 *  * bc -- the boundary conditions
 ******
 */

/****f* cuda_particle_kernel/reset_flag_w<<<>>>()
 * NAME
 *  reset_flag_w<<<>>>()
 * USAGE
 */
__global__ void reset_flag_w(int *flag_w, dom_struct *dom, BC bc);
/*
 * FUNCTION
 *  Set all flag_w nodes to fluid (= 1).
 * ARGUMENTS
 *  * flag_w -- the device z-direction velocity flag array subdomain
 *  * dom -- the device subdomain
 *  * bc -- the boundary conditions
 ******
 */


/****f* cuda_particle_kernel/cage_flag_u_periodic_W<<<>>>()
 * NAME
 *  cage_flag_u_periodic_W<<<>>>()
 * USAGE
 */
__global__ void cage_flag_u_periodic_W(int *flag_u, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain W boundary flag_u for periodic conditions.
 * ARGUMENTS
 *  * flag_u -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_u_periodic_E<<<>>>()
 * NAME
 *  cage_flag_u_periodic_E<<<>>>()
 * USAGE
 */
__global__ void cage_flag_u_periodic_E(int *flag_u, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain E boundary flag_u for periodic conditions.
 * ARGUMENTS
 *  * flag_u -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_u_periodic_S<<<>>>()
 * NAME
 *  cage_flag_u_periodic_S<<<>>>()
 * USAGE
 */
__global__ void cage_flag_u_periodic_S(int *flag_u, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain S boundary flag_u for periodic conditions.
 * ARGUMENTS
 *  * flag_u -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_u_periodic_N<<<>>>()
 * NAME
 *  cage_flag_u_periodic_N<<<>>>()
 * USAGE
 */
__global__ void cage_flag_u_periodic_N(int *flag_u, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain N boundary flag_u for periodic conditions.
 * ARGUMENTS
 *  * flag_u -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_u_periodic_B<<<>>>()
 * NAME
 *  cage_flag_u_periodic_B<<<>>>()
 * USAGE
 */
__global__ void cage_flag_u_periodic_B(int *flag_u, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain B boundary flag_u for periodic conditions.
 * ARGUMENTS
 *  * flag_u -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_u_periodic_T<<<>>>()
 * NAME
 *  cage_flag_u_periodic_T<<<>>>()
 * USAGE
 */
__global__ void cage_flag_u_periodic_T(int *flag_u, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain T boundary flag_u for periodic conditions.
 * ARGUMENTS
 *  * flag_u -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_v_periodic_W<<<>>>()
 * NAME
 *  cage_flag_v_periodic_W<<<>>>()
 * USAGE
 */
__global__ void cage_flag_v_periodic_W(int *flag_v, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain W boundary flag_v for periodic conditions.
 * ARGUMENTS
 *  * flag_v -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_v_periodic_E<<<>>>()
 * NAME
 *  cage_flag_v_periodic_E<<<>>>()
 * USAGE
 */
__global__ void cage_flag_v_periodic_E(int *flag_v, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain E boundary flag_v for periodic conditions.
 * ARGUMENTS
 *  * flag_v -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_v_periodic_S<<<>>>()
 * NAME
 *  cage_flag_v_periodic_S<<<>>>()
 * USAGE
 */
__global__ void cage_flag_v_periodic_S(int *flag_v, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain S boundary flag_v for periodic conditions.
 * ARGUMENTS
 *  * flag_v -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_v_periodic_N<<<>>>()
 * NAME
 *  cage_flag_v_periodic_N<<<>>>()
 * USAGE
 */
__global__ void cage_flag_v_periodic_N(int *flag_v, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain N boundary flag_v for periodic conditions.
 * ARGUMENTS
 *  * flag_v -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_v_periodic_B<<<>>>()
 * NAME
 *  cage_flag_v_periodic_B<<<>>>()
 * USAGE
 */
__global__ void cage_flag_v_periodic_B(int *flag_v, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain B boundary flag_v for periodic conditions.
 * ARGUMENTS
 *  * flag_v -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_v_periodic_T<<<>>>()
 * NAME
 *  cage_flag_v_periodic_T<<<>>>()
 * USAGE
 */
__global__ void cage_flag_v_periodic_T(int *flag_v, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain T boundary flag_v for periodic conditions.
 * ARGUMENTS
 *  * flag_v -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_w_periodic_W<<<>>>()
 * NAME
 *  cage_flag_w_periodic_W<<<>>>()
 * USAGE
 */
__global__ void cage_flag_w_periodic_W(int *flag_w, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain W boundary flag_w for periodic conditions.
 * ARGUMENTS
 *  * flag_w -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_w_periodic_E<<<>>>()
 * NAME
 *  cage_flag_w_periodic_E<<<>>>()
 * USAGE
 */
__global__ void cage_flag_w_periodic_E(int *flag_w, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain E boundary flag_w for periodic conditions.
 * ARGUMENTS
 *  * flag_w -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_w_periodic_S<<<>>>()
 * NAME
 *  cage_flag_w_periodic_S<<<>>>()
 * USAGE
 */
__global__ void cage_flag_w_periodic_S(int *flag_w, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain S boundary flag_w for periodic conditions.
 * ARGUMENTS
 *  * flag_w -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_w_periodic_N<<<>>>()
 * NAME
 *  cage_flag_w_periodic_N<<<>>>()
 * USAGE
 */
__global__ void cage_flag_w_periodic_N(int *flag_w, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain N boundary flag_w for periodic conditions.
 * ARGUMENTS
 *  * flag_w -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_w_periodic_B<<<>>>()
 * NAME
 *  cage_flag_w_periodic_B<<<>>>()
 * USAGE
 */
__global__ void cage_flag_w_periodic_B(int *flag_w, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain B boundary flag_w for periodic conditions.
 * ARGUMENTS
 *  * flag_w -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_particle_kernel/cage_flag_w_periodic_T<<<>>>()
 * NAME
 *  cage_flag_w_periodic_T<<<>>>()
 * USAGE
 */
__global__ void cage_flag_w_periodic_T(int *flag_w, dom_struct *dom);
/*
 * FUNCTION
 *  Update domain T boundary flag_w for periodic conditions.
 * ARGUMENTS
 *  * flag_w -- the device flag array subdomain
 *  * dom -- the device subdomain
 ******
 */














#endif
