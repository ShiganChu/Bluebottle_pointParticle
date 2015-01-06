/****h* Bluebottle/cuda_point_point_particle_kernel
 * NAME
 *  cuda_bluebottle_kernel
 * FUNCTION
 *  Bluebottle CUDA point_point_particle kernel functions.
 ******
 */

#ifndef _CUDA_point_particle_H
#define _CUDA_point_particle_H

extern "C"
{
#include "bluebottle.h"
#include "point.h"
}

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
