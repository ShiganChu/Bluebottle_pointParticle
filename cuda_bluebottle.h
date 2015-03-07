/****h* Bluebottle/cuda_bluebottle_kernel
 * NAME
 *  cuda_bluebottle_kernel
 * FUNCTION
 *  Bluebottle CUDA kernel functions.
 ******
 */

#ifndef _CUDA_BLUEBOTTLE_H
#define _CUDA_BLUEBOTTLE_H

extern "C"
{
#include "bluebottle.h"
#include "point.h"
}

/****f* cuda_bluebottle_kernel/BC_p_W_P<<<>>>()
 * NAME
 *  BC_p_W_P<<<>>>()
 * USAGE
 */
__global__ void BC_p_W_P(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the west face pressure field.
 * ARGUMENTS
 *  * p -- the device pressure field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_p_W_N<<<>>>()
 * NAME
 *  BC_p_W_N<<<>>>()
 * USAGE
 */
__global__ void BC_p_W_N(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the west face pressure field.
 * ARGUMENTS
 *  * p -- the device pressure field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_p_E_P<<<>>>()
 * NAME
 *  BC_p_E_P<<<>>>()
 * USAGE
 */
__global__ void BC_p_E_P(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the east face pressure field.
 * ARGUMENTS
 *  * p -- the device pressure field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_p_E_N<<<>>>()
 * NAME
 *  BC_p_E_N<<<>>>()
 * USAGE
 */
__global__ void BC_p_E_N(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the east face pressure field.
 * ARGUMENTS
 *  * p -- the device pressure field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_p_S_P<<<>>>()
 * NAME
 *  BC_p_S_P<<<>>>()
 * USAGE
 */
__global__ void BC_p_S_P(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the south face pressure field.
 * ARGUMENTS
 *  * p -- the device pressure field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_p_S_N<<<>>>()
 * NAME
 *  BC_p_S_N<<<>>>()
 * USAGE
 */
__global__ void BC_p_S_N(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the south face pressure field.
 * ARGUMENTS
 *  * p -- the device pressure field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_p_N_P<<<>>>()
 * NAME
 *  BC_p_N_P<<<>>>()
 * USAGE
 */
__global__ void BC_p_N_P(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the north face pressure field.
 * ARGUMENTS
 *  * p -- the device pressure field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_p_N_N<<<>>>()
 * NAME
 *  BC_p_N_N<<<>>>()
 * USAGE
 */
__global__ void BC_p_N_N(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the north face pressure field.
 * ARGUMENTS
 *  * p -- the device pressure field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_p_B_P<<<>>>()
 * NAME
 *  BC_p_B_P<<<>>>()
 * USAGE
 */
__global__ void BC_p_B_P(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the bottom face pressure field.
 * ARGUMENTS
 *  * p -- the device pressure field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_p_B_N<<<>>>()
 * NAME
 *  BC_p_B_N<<<>>>()
 * USAGE
 */
__global__ void BC_p_B_N(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the bottom face pressure field.
 * ARGUMENTS
 *  * p -- the device pressure field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_p_T_P<<<>>>()
 * NAME
 *  BC_p_T_P<<<>>>()
 * USAGE
 */
__global__ void BC_p_T_P(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the top face pressure field.
 * ARGUMENTS
 *  * p -- the device pressure field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_p_T_N<<<>>>()
 * NAME
 *  BC_p_T_N<<<>>>()
 * USAGE
 */
__global__ void BC_p_T_N(real *p, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the top face pressure field.
 * ARGUMENTS
 *  * p -- the device pressure field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_W_P<<<>>>()
 * NAME
 *  BC_u_W_P<<<>>>()
 * USAGE
 */
__global__ void BC_u_W_P(real *u, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the west face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_W_D<<<>>>()
 * NAME
 *  BC_u_W_D<<<>>>()
 * USAGE
 */
__global__ void BC_u_W_D(real *u, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the west face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_W_N<<<>>>()
 * NAME
 *  BC_u_W_N<<<>>>()
 * USAGE
 */
__global__ void BC_u_W_N(real *u, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the west face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_W_T<<<>>>()
 * NAME
 *  BC_u_W_T<<<>>>()
 * USAGE
 */
__global__ void BC_u_W_T(real *u, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the west face u-velocity
 *  field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device u-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_E_P<<<>>>()
 * NAME
 *  BC_u_E_P<<<>>>()
 * USAGE
 */
__global__ void BC_u_E_P(real *u, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the east face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_E_D<<<>>>()
 * NAME
 *  BC_u_E_D<<<>>>()
 * USAGE
 */
__global__ void BC_u_E_D(real *u, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the east face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_E_N<<<>>>()
 * NAME
 *  BC_u_E_N<<<>>>()
 * USAGE
 */
__global__ void BC_u_E_N(real *u, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the east face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_E_T<<<>>>()
 * NAME
 *  BC_u_E_T<<<>>>()
 * USAGE
 */
__global__ void BC_u_E_T(real *u, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the east face u-velocity
 *  field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device u-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_S_P<<<>>>()
 * NAME
 *  BC_u_S_P<<<>>>()
 * USAGE
 */
__global__ void BC_u_S_P(real *u, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the south face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_S_D<<<>>>()
 * NAME
 *  BC_u_S_D<<<>>>()
 * USAGE
 */
__global__ void BC_u_S_D(real *u, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the south face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_S_N<<<>>>()
 * NAME
 *  BC_u_S_N<<<>>>()
 * USAGE
 */
__global__ void BC_u_S_N(real *u, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the south face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_S_T<<<>>>()
 * NAME
 *  BC_u_S_T<<<>>>()
 * USAGE
 */
__global__ void BC_u_S_T(real *u, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the south face u-velocity
 *  field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device u-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_N_P<<<>>>()
 * NAME
 *  BC_u_N_P<<<>>>()
 * USAGE
 */
__global__ void BC_u_N_P(real *u, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the north face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_N_D<<<>>>()
 * NAME
 *  BC_u_N_D<<<>>>()
 * USAGE
 */
__global__ void BC_u_N_D(real *u, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the north face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_N_N<<<>>>()
 * NAME
 *  BC_u_N_N<<<>>>()
 * USAGE
 */
__global__ void BC_u_N_N(real *u, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the north face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_N_T<<<>>>()
 * NAME
 *  BC_u_N_T<<<>>>()
 * USAGE
 */
__global__ void BC_u_N_T(real *u, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the north face u-velocity
 *  field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device u-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_B_P<<<>>>()
 * NAME
 *  BC_u_B_P<<<>>>()
 * USAGE
 */
__global__ void BC_u_B_P(real *u, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the bottom face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_B_D<<<>>>()
 * NAME
 *  BC_u_B_D<<<>>>()
 * USAGE
 */
__global__ void BC_u_B_D(real *u, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the bottom face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_B_N<<<>>>()
 * NAME
 *  BC_u_B_N<<<>>>()
 * USAGE
 */
__global__ void BC_u_B_N(real *u, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the bottom face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_B_T<<<>>>()
 * NAME
 *  BC_u_B_T<<<>>>()
 * USAGE
 */
__global__ void BC_u_B_T(real *u, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the bottom face u-velocity
 *  field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device u-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_T_P<<<>>>()
 * NAME
 *  BC_u_T_P<<<>>>()
 * USAGE
 */
__global__ void BC_u_T_P(real *u, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the top face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_T_D<<<>>>()
 * NAME
 *  BC_u_T_D<<<>>>()
 * USAGE
 */
__global__ void BC_u_T_D(real *u, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the top face u-velocity field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_T_N<<<>>>()
 * NAME
 *  BC_u_T_N<<<>>>()
 * USAGE
 */
__global__ void BC_u_T_N(real *u, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the top face u-velocity field.
 * ARGUMENTS
 *  * u -- the device pressure field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_u_T_T<<<>>>()
 * NAME
 *  BC_u_T_T<<<>>>()
 * USAGE
 */
__global__ void BC_u_T_T(real *u, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the top face u-velocity
 *  field.
 * ARGUMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device u-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_W_P<<<>>>()
 * NAME
 *  BC_v_W_P<<<>>>()
 * USAGE
 */
__global__ void BC_v_W_P(real *v, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the west face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_W_D<<<>>>()
 * NAME
 *  BC_v_W_D<<<>>>()
 * USAGE
 */
__global__ void BC_v_W_D(real *v, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the west face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_W_N<<<>>>()
 * NAME
 *  BC_v_W_N<<<>>>()
 * USAGE
 */
__global__ void BC_v_W_N(real *v, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the west face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_W_T<<<>>>()
 * NAME
 *  BC_v_W_T<<<>>>()
 * USAGE
 */
__global__ void BC_v_W_T(real *v, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the west face v-velocity
 *  field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device v-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_E_P<<<>>>()
 * NAME
 *  BC_v_E_P<<<>>>()
 * USAGE
 */
__global__ void BC_v_E_P(real *v, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the east face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_E_D<<<>>>()
 * NAME
 *  BC_v_E_D<<<>>>()
 * USAGE
 */
__global__ void BC_v_E_D(real *v, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the east face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_E_N<<<>>>()
 * NAME
 *  BC_v_E_N<<<>>>()
 * USAGE
 */
__global__ void BC_v_E_N(real *v, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the east face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_E_T<<<>>>()
 * NAME
 *  BC_v_E_T<<<>>>()
 * USAGE
 */
__global__ void BC_v_E_T(real *v, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the east face v-velocity
 *  field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device v-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_S_P<<<>>>()
 * NAME
 *  BC_v_S_P<<<>>>()
 * USAGE
 */
__global__ void BC_v_S_P(real *v, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the south face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_S_D<<<>>>()
 * NAME
 *  BC_v_S_D<<<>>>()
 * USAGE
 */
__global__ void BC_v_S_D(real *v, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the south face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_S_N<<<>>>()
 * NAME
 *  BC_v_S_N<<<>>>()
 * USAGE
 */
__global__ void BC_v_S_N(real *v, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the south face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_S_T<<<>>>()
 * NAME
 *  BC_v_S_T<<<>>>()
 * USAGE
 */
__global__ void BC_v_S_T(real *v, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the south face v-velocity
 *  field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device v-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_N_P<<<>>>()
 * NAME
 *  BC_v_N_P<<<>>>()
 * USAGE
 */
__global__ void BC_v_N_P(real *v, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the north face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_N_D<<<>>>()
 * NAME
 *  BC_v_N_D<<<>>>()
 * USAGE
 */
__global__ void BC_v_N_D(real *v, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the north face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_N_N<<<>>>()
 * NAME
 *  BC_v_N_N<<<>>>()
 * USAGE
 */
__global__ void BC_v_N_N(real *v, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the north face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_N_T<<<>>>()
 * NAME
 *  BC_v_N_T<<<>>>()
 * USAGE
 */
__global__ void BC_v_N_T(real *v, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the north face v-velocity
 *  field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device v-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_B_P<<<>>>()
 * NAME
 *  BC_v_B_P<<<>>>()
 * USAGE
 */
__global__ void BC_v_B_P(real *v, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the bottom face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_B_D<<<>>>()
 * NAME
 *  BC_v_B_D<<<>>>()
 * USAGE
 */
__global__ void BC_v_B_D(real *v, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the bottom face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_B_N<<<>>>()
 * NAME
 *  BC_v_B_N<<<>>>()
 * USAGE
 */
__global__ void BC_v_B_N(real *v, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the bottom face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_B_T<<<>>>()
 * NAME
 *  BC_v_B_T<<<>>>()
 * USAGE
 */
__global__ void BC_v_B_T(real *v, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the bottom face v-velocity
 *  field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device v-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_T_P<<<>>>()
 * NAME
 *  BC_v_T_P<<<>>>()
 * USAGE
 */
__global__ void BC_v_T_P(real *v, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the top face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_T_D<<<>>>()
 * NAME
 *  BC_v_T_D<<<>>>()
 * USAGE
 */
__global__ void BC_v_T_D(real *v, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the top face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_T_N<<<>>>()
 * NAME
 *  BC_v_T_N<<<>>>()
 * USAGE
 */
__global__ void BC_v_T_N(real *v, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the top face v-velocity field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_v_T_T<<<>>>()
 * NAME
 *  BC_v_T_T<<<>>>()
 * USAGE
 */
__global__ void BC_v_T_T(real *v, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the top face v-velocity
 *  field.
 * ARGUMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device v-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_W_P<<<>>>()
 * NAME
 *  BC_w_W_P<<<>>>()
 * USAGE
 */
__global__ void BC_w_W_P(real *w, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the west face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_W_D<<<>>>()
 * NAME
 *  BC_w_W_D<<<>>>()
 * USAGE
 */
__global__ void BC_w_W_D(real *w, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the west face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_W_N<<<>>>()
 * NAME
 *  BC_w_W_N<<<>>>()
 * USAGE
 */
__global__ void BC_w_W_N(real *w, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the west face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_W_T<<<>>>()
 * NAME
 *  BC_w_W_T<<<>>>()
 * USAGE
 */
__global__ void BC_w_W_T(real *w, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the west face w-velocity
 *  field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device u-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_E_P<<<>>>()
 * NAME
 *  BC_w_E_P<<<>>>()
 * USAGE
 */
__global__ void BC_w_E_P(real *w, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the east face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_E_D<<<>>>()
 * NAME
 *  BC_w_E_D<<<>>>()
 * USAGE
 */
__global__ void BC_w_E_D(real *w, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the east face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_E_N<<<>>>()
 * NAME
 *  BC_w_E_N<<<>>>()
 * USAGE
 */
__global__ void BC_w_E_N(real *w, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the east face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_E_T<<<>>>()
 * NAME
 *  BC_w_E_T<<<>>>()
 * USAGE
 */
__global__ void BC_w_E_T(real *w, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the east face w-velocity
 *  field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device u-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_S_P<<<>>>()
 * NAME
 *  BC_w_S_P<<<>>>()
 * USAGE
 */
__global__ void BC_w_S_P(real *w, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the south face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_S_D<<<>>>()
 * NAME
 *  BC_w_S_D<<<>>>()
 * USAGE
 */
__global__ void BC_w_S_D(real *w, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the south face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_S_N<<<>>>()
 * NAME
 *  BC_w_S_N<<<>>>()
 * USAGE
 */
__global__ void BC_w_S_N(real *w, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the south face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_S_T<<<>>>()
 * NAME
 *  BC_w_S_T<<<>>>()
 * USAGE
 */
__global__ void BC_w_S_T(real *w, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the south face w-velocity
 *  field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device u-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_N_P<<<>>>()
 * NAME
 *  BC_w_N_P<<<>>>()
 * USAGE
 */
__global__ void BC_w_N_P(real *w, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the north face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_N_D<<<>>>()
 * NAME
 *  BC_w_N_D<<<>>>()
 * USAGE
 */
__global__ void BC_w_N_D(real *w, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the north face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_N_N<<<>>>()
 * NAME
 *  BC_w_N_N<<<>>>()
 * USAGE
 */
__global__ void BC_w_N_N(real *w, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the north face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_N_T<<<>>>()
 * NAME
 *  BC_w_N_T<<<>>>()
 * USAGE
 */
__global__ void BC_w_N_T(real *w, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the north face w-velocity
 *  field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device u-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_B_P<<<>>>()
 * NAME
 *  BC_w_B_P<<<>>>()
 * USAGE
 */
__global__ void BC_w_B_P(real *w, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the bottom face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_B_D<<<>>>()
 * NAME
 *  BC_w_B_D<<<>>>()
 * USAGE
 */
__global__ void BC_w_B_D(real *w, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the bottom face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_B_N<<<>>>()
 * NAME
 *  BC_w_B_N<<<>>>()
 * USAGE
 */
__global__ void BC_w_B_N(real *w, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the bottom face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_B_T<<<>>>()
 * NAME
 *  BC_w_B_T<<<>>>()
 * USAGE
 */
__global__ void BC_w_B_T(real *w, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the bottom face w-velocity
 *  field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device u-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_T_P<<<>>>()
 * NAME
 *  BC_w_T_P<<<>>>()
 * USAGE
 */
__global__ void BC_w_T_P(real *w, dom_struct *dom);
/*
 * FUNCTION
 *  Apply periodic boundary conditions to the top face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_T_D<<<>>>()
 * NAME
 *  BC_w_T_D<<<>>>()
 * USAGE
 */
__global__ void BC_w_T_D(real *w, dom_struct *dom, real bc);
/*
 * FUNCTION
 *  Apply Dirichlet boundary conditions to the top face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the value to be set on the boundary
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_T_N<<<>>>()
 * NAME
 *  BC_w_T_N<<<>>>()
 * USAGE
 */
__global__ void BC_w_T_N(real *w, dom_struct *dom);
/*
 * FUNCTION
 *  Apply Neumann boundary conditions to the top face w-velocity field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/BC_w_T_T<<<>>>()
 * NAME
 *  BC_w_T_T<<<>>>()
 * USAGE
 */
__global__ void BC_w_T_T(real *w, dom_struct *dom, real *bc);
/*
 * FUNCTION
 *  Apply the precursor domain boundary conditions to the top face w-velocity
 *  field.
 * ARGUMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * bc -- the device u-velocity precursor plane subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/project_u<<<>>>()
 * NAME
 *  project_u<<<>>>()
 * USAGE
 */
__global__ void project_u(real *u_star, real *p, real rho_f, real dt,
  real *u, dom_struct *dom, real ddx, int *flag_u);
/*
 * FUNCTION
 *  Project the intermediate velocity u_star onto a divergence-free space via
 *  p.
 * ARGUMENTS
 *  * u_star -- the intermediate velocity field
 *  * p -- the projected pressure
 *  * rho_f -- the fluid density
 *  * dt -- the timestep
 *  * u -- the solution velocity field
 *  * dom -- the subdomain on which to operate
 *  * ddx -- 1 / dx
 *  * flag_u -- the x-direction boundary flag
 ******
 */

/****f* cuda_bluebottle_kernel/project_v<<<>>>()
 * NAME
 *  project_v<<<>>>()
 * USAGE
 */
__global__ void project_v(real *v_star, real *p, real rho_f, real dt,
  real *v, dom_struct *dom, real ddy, int *flag_v);
/*
 * FUNCTION
 *  Project the intermediate velocity v_star onto a divergence-free space via
 *  p.
 * ARGUMENTS
 *  * v_star -- the intermediate velocity field
 *  * p -- the projected pressure
 *  * rho_f -- the fluid density
 *  * dt -- the timestep
 *  * v -- the solution velocity field
 *  * dom -- the subdomain on which to operate
 *  * ddy -- 1 / dy
 *  * flag_v -- the y-direction boundary flag
 ******
 */

/****f* cuda_bluebottle_kernel/project_w<<<>>>()
 * NAME
 *  project_w<<<>>>()
 * USAGE
 */
__global__ void project_w(real *w_star, real *p, real rho_f, real dt,
  real *w, dom_struct *dom, real ddz, int *flag_w);
/*
 * FUNCTION
 *  Project the intermediate velocity w_star onto a divergence-free space via
 *  p.
 * ARGUMENTS
 *  * w_star -- the intermediate velocity field
 *  * p -- the projected pressure
 *  * rho_f -- the fluid density
 *  * dt -- the timestep
 *  * w -- the solution velocity field
 *  * dom -- the subdomain on which to operate
 *  * ddz -- 1 / dz
 *  * flag_w -- the z-direction boundary flag
 ******
 */

/****f* cuda_bluebottle_kernel/copy_p_ghost<<<>>>()
 * NAME
 *  copy_p_ghost<<<>>>()
 * USAGE
 */
__global__ void copy_p_ghost(real *p, real *p_tmp, dom_struct *dom);
/*
 * FUNCTION
 *  Copy the p solution of the pressure-Poisson problem (without ghost
 *  cells) to _p (with ghost cells).
 * ARGUMENTS
 *  * p -- the destination data structure with ghost cells
 *  * p_tmp -- the source data structure without ghost cells
 *  * dom -- the subdomain on which to operate
 ******
 */

/****f* cuda_bluebottle_kernel/copy_p_noghost<<<>>>()
 * NAME
 *  copy_p_noghost<<<>>>()
 * USAGE
 */
__global__ void copy_p_noghost(real *p_noghost, real *p_ghost, dom_struct *dom);
/*
 * FUNCTION
 *  Copy the pressure field containing ghost cells to a new array without
 *  ghost cells.
 * ARGUMENTS
 *  p_noghost -- the destination data structure without ghost cells
 *  p_ghost -- the source data structure with ghost cells
 *  dom -- the subdomain on which to operate
 ******
 */

/****f* cuda_bluebottle_kernel/copy_u_noghost<<<>>>()
 * NAME
 *  copy_u_noghost<<<>>>()
 * USAGE
 */
__global__ void copy_u_noghost(real *u_noghost, real *u_ghost, dom_struct *dom);
/*
 * FUNCTION
 *  Copy the u-velocity field containing ghost cells to a new array without
 *  ghost cells.
 * ARGUMENTS
 *  u_noghost -- the destination data structure without ghost cells
 *  u_ghost -- the source data structure with ghost cells
 *  dom -- the subdomain on which to operate
 ******
 */

/****f* cuda_bluebottle_kernel/copy_v_noghost<<<>>>()
 * NAME
 *  copy_v_noghost<<<>>>()
 * USAGE
 */
__global__ void copy_v_noghost(real *v_noghost, real *v_ghost, dom_struct *dom);
/*
 * FUNCTION
 *  Copy the v-velocity field containing ghost cells to a new array without
 *  ghost cells.
 * ARGUMENTS
 *  v_noghost -- the destination data structure without ghost cells
 *  v_ghost -- the source data structure with ghost cells
 *  dom -- the subdomain on which to operate
 ******
 */

/****f* cuda_bluebottle_kernel/copy_w_noghost<<<>>>()
 * NAME
 *  copy_w_noghost<<<>>>()
 * USAGE
 */
__global__ void copy_w_noghost(real *w_noghost, real *w_ghost, dom_struct *dom);
/*
 * FUNCTION
 *  Copy the w-velocity field containing ghost cells to a new array without
 *  ghost cells.
 * ARGUMENTS
 *  w_noghost -- the destination data structure without ghost cells
 *  w_ghost -- the source data structure with ghost cells
 *  dom -- the subdomain on which to operate
 ******
 */


__global__ void copy_u_fluid(real *u_noghost, real *u_ghost, int *phase, dom_struct *dom);
__global__ void copy_v_fluid(real *v_noghost, real *v_ghost, int *phase, dom_struct *dom);
__global__ void copy_w_fluid(real *w_noghost, real *w_ghost, int *phase, dom_struct *dom);


/****f* cuda_bluebottle_kernel/u_star_2<<<>>>()
 * NAME
 *  u_star_2<<<>>>()
 * USAGE
 */
__global__ void u_star_2(real rho_f, real nu,
  real *u0, real *v0, real *w0, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *u_star,
  dom_struct *dom, real dt0, real dt);
/*
 * FUNCTION
 *  Compute the intermediate velocity field u_star (2nd-order in time).
 * ARGUMENTS
 *  * rho_f -- fluid density
 *  * nu -- fluid kinematic viscosity
 *  * u0 -- device subdomain u-component velocity field from previous timestep
 *  * v0 -- device subdomain v-component velocity field from previous timestep
 *  * w0 -- device subdomain w-component velocity field from previous timestep
 *  * f -- the forcing array
 *  * diff0 -- device subdomain previous diffusion term
 *  * conv0 -- device subdomain previous convection term
 *  * u_star -- the intermediate velocity field
 *  * dom -- the subdomain in which this device is operating
 *  * dt0 -- the previous timestep
 *  * dt -- the current timestep
 ******
 */

/****f* cuda_bluebottle_kernel/v_star_2<<<>>>()
 * NAME
 *  v_star_2<<<>>>()
 * USAGE
 */
__global__ void v_star_2(real rho_f, real nu,
  real *u0, real *v0, real *w0, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *v_star,
  dom_struct *dom, real dt0, real dt);
/*
 * FUNCTION
 *  Compute the intermediate velocity field v_star (2nd-order in time).
 * ARGUMENTS
 *  * rho_f -- fluid density
 *  * nu -- fluid kinematic viscosity
 *  * u0 -- device subdomain u-component velocity field from previous timestep
 *  * v0 -- device subdomain v-component velocity field from previous timestep
 *  * w0 -- device subdomain w-component velocity field from previous timestep
 *  * f -- the forcing array
 *  * diff0 -- device subdomain previous diffusion term
 *  * conv0 -- device subdomain previous convection term
 *  * v_star -- the intermediate velocity field
 *  * dom -- the subdomain in which this device is operating
 *  * dt0 -- the previous timestep
 *  * dt -- the current timestep
 ******
 */

/****f* cuda_bluebottle_kernel/w_star_2<<<>>>()
 * NAME
 *  w_star_2<<<>>>()
 * USAGE
 */
__global__ void w_star_2(real rho_f, real nu,
  real *u0, real *v0, real *w0, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *w_star,
  dom_struct *dom, real dt0, real dt);
/*
 * FUNCTION
 *  Compute the intermediate velocity field w_star (2nd-order in time).
 * ARGUMENTS
 *  * rho_f -- fluid density
 *  * nu -- fluid kinematic viscosity
 *  * u0 -- device subdomain u-component velocity field from previous timestep
 *  * v0 -- device subdomain v-component velocity field from previous timestep
 *  * w0 -- device subdomain w-component velocity field from previous timestep
 *  * f -- the forcing array
 *  * diff0 -- device subdomain previous diffusion term
 *  * conv0 -- device subdomain previous convection term
 *  * w_star -- the intermediate velocity field
 *  * dom -- the subdomain in which this device is operating
 *  * dt0 -- the previous timestep
 *  * dt -- the current timestep
 ******
 */

/****f* cuda_bluebottle_kernel/u_star_2_init<<<>>>()
 * NAME
 *  u_star_2_init<<<>>>()
 * USAGE
 */
__global__ void u_star_2_init(real rho_f, real nu,
  real *u, real *v, real *w, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *u_star,
  dom_struct *dom, real dt);
/*
 * FUNCTION
 *  Compute the intermediate velocity field u_star (2nd-order in time).
 *  This version uses the Euler method for initializing the Adams-Bashforth
 *  stepping.
 * ARGUMENTS
 *  * rho_f -- fluid density
 *  * nu -- fluid kinematic viscosity
 *  * u -- device subdomain u-component velocity field
 *  * v -- device subdomain v-component velocity field
 *  * w -- device subdomain w-component velocity field
 *  * f -- the forcing array
 *  * diff0 -- device subdomain previous diffusion term
 *  * conv0 -- device subdomain previous convection term
 *  * u_star -- the intermediate velocity field
 *  * dom -- the subdomain in which this device is operating
 *  * dt -- the current timestep
 ******
 */

/****f* cuda_bluebottle_kernel/v_star_2_init<<<>>>()
 * NAME
 *  v_star_2_init<<<>>>()
 * USAGE
 */
__global__ void v_star_2_init(real rho_f, real nu,
  real *u, real *v, real *w, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *v_star,
  dom_struct *dom, real dt);
/*
 * FUNCTION
 *  Compute the intermediate velocity field v_star (2nd-order in time).
 *  This version uses the Euler method for initializing the Adams-Bashforth
 *  stepping.
 * ARGUMENTS
 *  * rho_f -- fluid density
 *  * nu -- fluid kinematic viscosity
 *  * u -- device subdomain u-component velocity field
 *  * v -- device subdomain v-component velocity field
 *  * w -- device subdomain w-component velocity field
 *  * f -- the forcing array
 *  * diff0 -- device subdomain previous diffusion term
 *  * conv0 -- device subdomain previous convection term
 *  * v_star -- the intermediate velocity field
 *  * dom -- the subdomain in which this device is operating
 *  * dt -- the current timestep
 ******
 */

/****f* cuda_bluebottle_kernel/w_star_2_init<<<>>>()
 * NAME
 *  w_star_2_init<<<>>>()
 * USAGE
 */
__global__ void w_star_2_init(real rho_f, real nu,
  real *u, real *v, real *w, real *f,
  real *diff0, real *conv0, real *diff, real *conv, real *w_star,
  dom_struct *dom, real dt);
/*
 * FUNCTION
 *  Compute the intermediate velocity field w_star (2nd-order in time).
 *  This version uses the Euler method for initializing the Adams-Bashforth
 *  stepping.
 * ARGUMENTS
 *  * rho_f -- fluid density
 *  * nu -- fluid kinematic viscosity
 *  * u -- device subdomain u-component velocity field
 *  * v -- device subdomain v-component velocity field
 *  * w -- device subdomain w-component velocity field
 *  * f -- the forcing array
 *  * diff0 -- device subdomain diffusion term
 *  * conv0 -- device subdomain convection term
 *  * w_star -- the intermediate velocity field
 *  * dom -- the subdomain in which this device is operating
 *  * dt -- the current timestep
 ******
 */


//Testing of the force subroutine   u=cos(qx)cos(qy)sin(rt); v=sin(qx)sin(qy)sin(rt)
__global__ void forcing_test_taylor(real *fx,real *fy, dom_struct *dom,real t,real DIFF,point_struct *points,real n2,real n3);

//combined with scalar diffusion, to test scalar schemes
__global__ void forcing_test_taylor_sc_ConvDiff(real *fx, real *fy, dom_struct *dom,real t,real nu,point_struct *points);






/****f* cuda_bluebottle_kernel/forcing_reset_x<<<>>>()
 * NAME
 *  forcing_reset_x<<<>>>()
 * USAGE
 */
__global__ void forcing_reset_x(real *fx, dom_struct *dom);
/*
 * FUNCTION 
 *  Reset the x-direction forcing array to zero.
 * ARGUMENTS
 *  * fx -- the forcing array to be reset
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/forcing_reset_y<<<>>>()
 * NAME
 *  forcing_reset_y<<<>>>()
 * USAGE
 */
__global__ void forcing_reset_y(real *fy, dom_struct *dom);
/*
 * FUNCTION 
 *  Reset the y-direction forcing array to zero.
 * ARGUMENTS
 *  * fy -- the forcing array to be reset
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/forcing_reset_z<<<>>>()
 * NAME
 *  forcing_reset_z<<<>>>()
 * USAGE
 */
__global__ void forcing_reset_z(real *fz, dom_struct *dom);
/*
 * FUNCTION 
 *  Reset the z-direction forcing array to zero.
 * ARGUMENTS
 *  * fz -- the forcing array to be reset
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/forcing_add_x_const<<<>>>()
 * NAME
 *  forcing_add_x_const<<<>>>()
 * USAGE
 */
__global__ void forcing_add_x_const(real val, real *fx, dom_struct *dom);
/*
 * FUNCTION 
 *  Add a constant force to the forcing array.
 * ARGUMENTS
 *  * val -- the value of the force to be added to the array
 *  * fx -- the forcing array to be reset
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/forcing_add_y_const<<<>>>()
 * NAME
 *  forcing_add_y_const<<<>>>()
 * USAGE
 */
__global__ void forcing_add_y_const(real val, real *fy, dom_struct *dom);
/*
 * FUNCTION 
 *  Add a constant force to the forcing array.
 * ARGUMENTS
 *  * val -- the value of the force to be added to the array
 *  * fy -- the forcing array to be reset
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/forcing_add_z_const<<<>>>()
 * NAME
 *  forcing_add_z_const<<<>>>()
 * USAGE
 */
__global__ void forcing_add_z_const(real val, real *fz, dom_struct *dom);
/*
 * FUNCTION 
 *  Add a constant force to the forcing array.
 * ARGUMENTS
 *  * val -- the value of the force to be added to the array
 *  * fz -- the forcing array to be reset
 *  * dom -- the current subdomain
 ******
 */

__global__ void forcing_add_z_gravity_pressure(real val, real *fz, dom_struct *dom);



/****f* cuda_bluebottle_kernel/forcing_add_x_field<<<>>>()
 * NAME
 *  forcing_add_x_field<<<>>>()
 * USAGE
 */
__global__ void forcing_add_x_field(real scale, real *val, real *fx,
  dom_struct *dom);
/*
 * FUNCTION 
 *  Add a field force to the forcing array.
 * ARGUMENTS
 *  * scale -- a constant scaling
 *  * val -- the value of the force to be added to the array
 *  * fx -- the forcing array to be reset
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/forcing_add_y_field<<<>>>()
 * NAME
 *  forcing_add_y_field<<<>>>()
 * USAGE
 */
__global__ void forcing_add_y_field(real scale, real *val, real *fy,
  dom_struct *dom);
/*
 * FUNCTION 
 *  Add a field force to the forcing array.
 * ARGUMENTS
 *  * scale -- a constant scaling
 *  * val -- the value of the force to be added to the array
 *  * fy -- the forcing array to be reset
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/forcing_add_z_field<<<>>>()
 * NAME
 *  forcing_add_z_field<<<>>>()
 * USAGE
 */
__global__ void forcing_add_z_field(real scale, real *val, real *fz,
  dom_struct *dom);
/*
 * FUNCTION 
 *  Add a field force to the forcing array.
 * ARGUMENTS
 *  * scale -- a constant scaling
 *  * val -- the value of the force to be added to the array
 *  * fz -- the forcing array to be reset
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/surf_int_x_copy<<<>>>()
 * NAME
 *  surf_int_x_copy<<<>>>()
 * USAGE
 */
__global__ void surf_int_x_copy(real *u_star, real *u_star_tmp,
  dom_struct *dom);
/*
 * FUNCTION
 *  Copy the West and East faces of u_star to a temporary vector to be
 *  summed in order to calculate the surface integral u*.n on these faces.
 * ARGUMENTS
 *  * u_star -- subdomain intermediate velocity in the x-direction on the GPU
 *  * u_star_tmp -- the location to which to copy point of u_star
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/surf_int_y_copy<<<>>>()
 * NAME
 *  surf_int_y_copy<<<>>>()
 * USAGE
 */
__global__ void surf_int_y_copy(real *v_star, real *v_star_tmp,
  dom_struct *dom);
/*
 * FUNCTION
 *  Copy the South and North faces of v_star to a temporary vector to be
 *  summed in order to calculate the surface integral u*.n on these faces.
 * ARGUMENTS
 *  * v_star -- subdomain intermediate velocity in the y-direction on the GPU
 *  * v_star_tmp -- the location to which to copy point of v_star
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/surf_int_z_copy<<<>>>()
 * NAME
 *  surf_int_z_copy<<<>>>()
 * USAGE
 */
__global__ void surf_int_z_copy(real *w_star, real *w_star_tmp,
  dom_struct *dom);
/*
 * FUNCTION
 *  Copy the Bottom and Top faces of w_star to a temporary vector to be
 *  summed in order to calculate the surface integral u*.n on these faces.
 * ARGUMENTS
 *  * w_star -- subdomain intermediate velocity in the w-direction on the GPU
 *  * w_star_tmp -- the location to which to copy point of w_star
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/plane_eps_x_W<<<>>>()
 * NAME
 *  plane_eps_x_W<<<>>>()
 * USAGE
 */
__global__ void plane_eps_x_W(real eps, real *u_star, dom_struct *dom);
/*
 * FUNCTION
 *  Subtract eps from each node of the WEST face to fudge the outflot plane
 *  just enough for the solvablility condition to hold to machine accuracy.
 * ARGUMENTS
 *  * eps -- the epsilon value to subtract from each node of the outflow plane
 *  * u_star -- the subdomain u_star velocity from which to subtract eps
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/plane_eps_x_E<<<>>>()
 * NAME
 *  plane_eps_x_E<<<>>>()
 * USAGE
 */
__global__ void plane_eps_x_E(real eps, real *u_star, dom_struct *dom);
/*
 * FUNCTION
 *  Subtract eps from each node of the EAST face to fudge the outflot plane
 *  just enough for the solvablility condition to hold to machine accuracy.
 * ARGUMENTS
 *  * eps -- the epsilon value to subtract from each node of the outflow plane
 *  * u_star -- the subdomain u_star velocity from which to subtract eps
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/plane_eps_y_S<<<>>>()
 * NAME
 *  plane_eps_y_S<<<>>>()
 * USAGE
 */
__global__ void plane_eps_y_S(real eps, real *v_star, dom_struct *dom);
/*
 * FUNCTION
 *  Subtract eps from each node of the SOUTH face to fudge the outflot plane
 *  just enough for the solvablility condition to hold to machine accuracy.
 * ARGUMENTS
 *  * eps -- the epsilon value to subtract from each node of the outflow plane
 *  * v_star -- the subdomain v_star velocity from which to subtract eps
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/plane_eps_y_N<<<>>>()
 * NAME
 *  plane_eps_y_N<<<>>>()
 * USAGE
 */
__global__ void plane_eps_y_N(real eps, real *v_star, dom_struct *dom);
/*
 * FUNCTION
 *  Subtract eps from each node of the North face to fudge the outflot plane
 *  just enough for the solvablility condition to hold to machine accuracy.
 * ARGUMENTS
 *  * eps -- the epsilon value to subtract from each node of the outflow plane
 *  * v_star -- the subdomain v_star velocity from which to subtract eps
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/plane_eps_z_B<<<>>>()
 * NAME
 *  plane_eps_z_B<<<>>>()
 * USAGE
 */
__global__ void plane_eps_z_B(real eps, real *w_star, dom_struct *dom);
/*
 * FUNCTION
 *  Subtract eps from each node of the BOTTOM face to fudge the outflot plane
 *  just enough for the solvablility condition to hold to machine accuracy.
 * ARGUMENTS
 *  * eps -- the epsilon value to subtract from each node of the outflow plane
 *  * w_star -- the subdomain w_star velocity from which to subtract eps
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/plane_eps_z_T<<<>>>()
 * NAME
 *  plane_eps_z_T<<<>>>()
 * USAGE
 */
__global__ void plane_eps_z_T(real eps, real *w_star, dom_struct *dom);
/*
 * FUNCTION
 *  Subtract eps from each node of the TOP face to fudge the outflot plane
 *  just enough for the solvablility condition to hold to machine accuracy.
 * ARGUMENTS
 *  * eps -- the epsilon value to subtract from each node of the outflow plane
 *  * w_star -- the subdomain w_star velocity from which to subtract eps
 *  * dom -- the current subdomain
 ******
 */

/****f* cuda_bluebottle_kernel/move_points_a<<<>>>()
 * NAME
 *  move_points_a<<<>>>()
 * USAGE
 */
//__global__ void move_points_a(dom_struct *dom, point_struct *points, int npoints,
//  real dt, real dt0, g_struct g, real rho_f, real ttime);
/*
 * FUNCTION
 *  Update the point_particle velocities and move the point_particles. point A: does
 *  everything but move the point_particles. This way we can use the velocity
 *  expected at the end of the timestep to update the lubrication forces
 *  before updating the point_particle position.
 * ARGUMENTS
 *  * dom -- the subdomain
 *  * points -- the point_particles
 *  * npoints -- the number of point_particles
 *  * dt -- timestep
 *  * dt0 -- previous timestep
 *  * g -- body force
 ******
 */

/****f* cuda_bluebottle_kernel/move_points_b<<<>>>()
 * NAME
 *  move_points_b<<<>>>()
 * USAGE
 */
//__global__ void move_points_b(dom_struct *dom, point_struct *points, int npoints,
//  real dt, real dt0, g_struct g, real rho_f, real ttime);
/*
 * FUNCTION
 *  Update the point_particle velocities and move the point_particles. point B: does
 *  everything includin move the point_particles. This way we can use the velocity
 *  expected at the end of the timestep to update the lubrication forces
 *  before updating the point_particle position.
 * ARGUMENTS
 *  * dom -- the subdomain
 *  * points -- the point_particles
 *  * npoints -- the number of point_particles
 *  * dt -- timestep
 *  * dt0 -- previous timestep
 *  * g -- body force
 ******
 */

/****f* cuda_bluebottle_kernel/detect_collision<<<>>>()
 * NAME
 *  detect_collision<<<>>>()
 * USAGE
 */
/*__global__ void detect_collision(dom_struct *dom, point_struct *points,
  real *dt, real *dt0, real *dt_sum, real *dt_tmp, real *udotnstore,
  real *nxstore, real *nystore, real *nzstore, g_struct g, int npoints,
  int bcuw, int bcvs, int bcwb);
*/
/*
 * FUNCTION
 *  Update the point_particle velocities and move the point_particles.
 * ARGUMENTS
 *  * dom -- the subdomain
 *  * points -- the point_particles
 *  * dt -- timestep
 *  * dt0 -- previous timestep
 *  * g -- body force
 *  * npoints -- number of point_particles
 ******
 */

//__global__ void update_collision(point_struct *points);

/****f* cuda_bluebottle_kernel/rotate()
 * NAME
 *  rotate
 * USAGE
 */
__device__ void rotate(real qr, real qi, real qj, real qk,
  real *pi, real *pj, real *pk);
/*
 * FUNCTION
 *  Apply quaternion conjugation p <-- q^-1 * p * q to rotate p according to
 *  rotation described by quaternion q.
 * ARGUMENTS
 *  * qr -- quaternion real component
 *  * qi -- quaternion imaginary component i
 *  * qj -- quaternion imaginary component j
 *  * qk -- quaternion imaginary component k
 *  * pi -- vector component i
 *  * pj -- vector component j
 *  * pk -- vector component k
 ******
 */

/****f* cuda_bluebottle_kernel/collision_init<<<>>>()
 * NAME
 *  collision_init<<<>>>()
 * USAGE
 */
__global__ void collision_init(point_struct *points, int npoints);
/*
 * FUNCTION
 *  Set interaction forces equal to zero to start.
 * ARGUMENTS
 *  * points -- the device point_particle array subdomain
 *  * npoints -- the number of point_particles
 ******
 */

/****f* cuda_bluebottle_kernel/collision_points<<<>>>()
 * NAME
 *  collision_points<<<>>>()
 * USAGE
 */
__global__ void collision_points(point_struct *points, int i,
  dom_struct *dom, real eps, real *forces, real *moments, int npoints, real mu,
  BC bc);
/*
 * FUNCTION
 *  Calculate collision forcing between point_particle i and all other point_particles.
 * ARGUMENTS
 *  * points -- the device point_particle array subdomain
 *  * i -- the point_particle number to calculate
 *  * eps -- magnitude of forcing
 *  * forces -- temporary working array
 *  * moments -- temporary working array
 *  * npoints -- number of point_particles
 *  * mu -- fluid viscosity
 ******
 */

/****f* cuda_bluebottle_kernel/collision_walls<<<>>>()
 * NAME
 *  collision_walls<<<>>>()
 * USAGE
 */
__global__ void collision_walls(dom_struct *dom, point_struct *points,
  int npoints, BC bc, real eps, real mu);
/*
 * FUNCTION
 *  Calculate collision forcing between point_particle i and all other point_particles.
 * ARGUMENTS
 *  * dom -- the device domain array
 *  * points -- the device point_particle array subdomain
 *  * npoints -- the number of point_particles
 *  * bc -- boundary condition data
 *  * eps -- magnitude of forcing
 *  * mu -- fluid viscosity
 ******
 */

/****f* cuda_bluebottle_kernel/spring_points<<<>>>()
 * NAME
 *  spring_points<<<>>>()
 * USAGE
 */
__global__ void spring_points(point_struct *points, int npoints);
/*
 * FUNCTION
 *  Calculate spring force pulling point_particle back to origin.
 * ARGUMENTS
 *  * points -- the device point_particle array subdomain
 *  * npoints -- the number of point_particles
 ******
 */

/****v* cuda_bluebottle_kernel/yank_u_WE<<<>>>()
 * NAME
 *  yank_u_WE<<<>>>()
 * USAGE
 */
__global__ void yank_u_WE(real *u, dom_struct *dom, real *plane, real xpos,
  real vel);
/*
 * FUNCTION
 *  Interpolate and copy the turbulent precursor u-velocity plane parallel to
 *  the W and E faces at x-position xpos.
 * ARUGMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * plane -- u-velocity plane parallel to the W and E faces to which to copy
 *    the interpolated values
 *  * xpos -- the x-position of plane
 *  * vel -- the velocity at which the plane is propagating (the inflow
 *    velocity)
 ******
 */

/****v* cuda_bluebottle_kernel/yank_v_WE<<<>>>()
 * NAME
 *  yank_v_WE<<<>>>()
 * USAGE
 */
__global__ void yank_v_WE(real *v, dom_struct *dom, real *plane, real xpos,
  real vel);
/*
 * FUNCTION
 *  Interpolate and copy the turbulent precursor v-velocity plane parallel to
 *  the W and E faces at x-position xpos.
 * ARUGMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * plane -- v-velocity plane parallel to the W and E faces to which to copy
 *    the interpolated values
 *  * xpos -- the x-position of plane
 *  * vel -- the velocity at which the plane is propagating (the inflow
 *    velocity)
 ******
 */

/****v* cuda_bluebottle_kernel/yank_w_WE<<<>>>()
 * NAME
 *  yank_w_WE<<<>>>()
 * USAGE
 */
__global__ void yank_w_WE(real *w, dom_struct *dom, real *plane, real xpos,
  real vel);
/*
 * FUNCTION
 *  Interpolate and copy the turbulent precursor w-velocity plane parallel to
 *  the W and E faces at x-position xpos.
 * ARUGMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * plane -- w-velocity plane parallel to the W and E faces to which to copy
 *    the interpolated values
 *  * xpos -- the x-position of plane
 *  * vel -- the velocity at which the plane is propagating (the inflow
 *    velocity)
 ******
 */

/****v* cuda_bluebottle_kernel/yank_u_SN<<<>>>()
 * NAME
 *  yank_u_SN<<<>>>()
 * USAGE
 */
__global__ void yank_u_SN(real *u, dom_struct *dom, real *plane, real ypos,
  real vel);
/*
 * FUNCTION
 *  Interpolate and copy the turbulent precursor u-velocity plane parallel to
 *  the S and N faces at y-position ypos.
 * ARUGMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * plane -- u-velocity plane parallel to the W and E faces to which to copy
 *    the interpolated values
 *  * ypos -- the y-position of plane
 *  * vel -- the velocity at which the plane is propagating (the inflow
 *    velocity)
 ******
 */

/****v* cuda_bluebottle_kernel/yank_v_SN<<<>>>()
 * NAME
 *  yank_v_SN<<<>>>()
 * USAGE
 */
__global__ void yank_v_SN(real *v, dom_struct *dom, real *plane, real ypos,
  real vel);
/*
 * FUNCTION
 *  Interpolate and copy the turbulent precursor v-velocity plane parallel to
 *  the S and N faces at y-position ypos.
 * ARUGMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * plane -- v-velocity plane parallel to the W and E faces to which to copy
 *    the interpolated values
 *  * ypos -- the y-position of plane
 *  * vel -- the velocity at which the plane is propagating (the inflow
 *    velocity)
 ******
 */

/****v* cuda_bluebottle_kernel/yank_w_SN<<<>>>()
 * NAME
 *  yank_w_SN<<<>>>()
 * USAGE
 */
__global__ void yank_w_SN(real *w, dom_struct *dom, real *plane, real ypos,
  real vel);
/*
 * FUNCTION
 *  Interpolate and copy the turbulent precursor w-velocity plane parallel to
 *  the S and N faces at y-position ypos.
 * ARUGMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * plane -- w-velocity plane parallel to the W and E faces to which to copy
 *    the interpolated values
 *  * ypos -- the y-position of plane
 *  * vel -- the velocity at which the plane is propagating (the inflow
 *    velocity)
 ******
 */

/****v* cuda_bluebottle_kernel/yank_u_BT<<<>>>()
 * NAME
 *  yank_u_BT<<<>>>()
 * USAGE
 */
__global__ void yank_u_BT(real *u, dom_struct *dom, real *plane, real zpos,
  real vel);
/*
 * FUNCTION
 *  Interpolate and copy the turbulent precursor u-velocity plane parallel to
 *  the B and T faces at z-position zpos.
 * ARUGMENTS
 *  * u -- the device u-velocity field subdomain
 *  * dom -- the device subdomain
 *  * plane -- u-velocity plane parallel to the W and E faces to which to copy
 *    the interpolated values
 *  * zpos -- the z-position of plane
 *  * vel -- the velocity at which the plane is propagating (the inflow
 *    velocity)
 ******
 */

/****v* cuda_bluebottle_kernel/yank_v_BT<<<>>>()
 * NAME
 *  yank_v_BT<<<>>>()
 * USAGE
 */
__global__ void yank_v_BT(real *v, dom_struct *dom, real *plane, real zpos,
  real vel);
/*
 * FUNCTION
 *  Interpolate and copy the turbulent precursor v-velocity plane parallel to
 *  the B and T faces at z-position zpos.
 * ARUGMENTS
 *  * v -- the device v-velocity field subdomain
 *  * dom -- the device subdomain
 *  * plane -- v-velocity plane parallel to the W and E faces to which to copy
 *    the interpolated values
 *  * zpos -- the z-position of plane
 *  * vel -- the velocity at which the plane is propagating (the inflow
 *    velocity)
 ******
 */

/****v* cuda_bluebottle_kernel/yank_w_BT<<<>>>()
 * NAME
 *  yank_w_BT<<<>>>()
 * USAGE
 */
__global__ void yank_w_BT(real *w, dom_struct *dom, real *plane, real zpos,
  real vel);
/*
 * FUNCTION
 *  Interpolate and copy the turbulent precursor w-velocity plane parallel to
 *  the B and T faces at z-position zpos.
 * ARUGMENTS
 *  * w -- the device w-velocity field subdomain
 *  * dom -- the device subdomain
 *  * plane -- w-velocity plane parallel to the W and E faces to which to copy
 *    the interpolated values
 *  * zpos -- the z-position of plane
 *  * vel -- the velocity at which the plane is propagating (the inflow
 *    velocity)
 ******
 */

/****f* bluebottle_kernel/colocate_Gfx<<<>>>()
 * NAME
 *  colocate_Gfx<<<>>>()
 * TYPE
 */
__global__ void colocate_Gfx(real *u, real *u_co, dom_struct *dom);
/*
 * PURPOSE
 *  Interpolate a Gfx grid (e.g. u-velocity) to a Gcc grid. This routine
 *  interpolates from a Gfx grid with ghost cell buffers to a Gcc grid with
 *  no ghost cell buffers.
 * ARGUMENTS
 *  * u -- device Gfx field with ghost cell buffers
 *  * u_co -- device colocated Gcc field with no ghost cell buffers
 *  * dom -- device domain data
 ******
 */

/****f* bluebottle_kernel/colocate_Gfy<<<>>>()
 * NAME
 *  colocate_Gfy<<<>>>()
 * TYPE
 */
__global__ void colocate_Gfy(real *v, real *v_co, dom_struct *dom);
/*
 * PURPOSE
 *  Interpolate a Gfy grid (e.g. v-velocity) to a Gcc grid. This routine
 *  interpolates from a Gfy grid with ghost cell buffers to a Gcc grid with
 *  no ghost cell buffers.
 * ARGUMENTS
 *  * v -- device Gfy field with ghost cell buffers
 *  * v_co -- device colocated Gcc field with no ghost cell buffers
 *  * dom -- device domain data
 ******
 */

/****f* bluebottle_kernel/colocate_Gfz<<<>>>()
 * NAME
 *  colocate_Gfz<<<>>>()
 * TYPE
 */
__global__ void colocate_Gfz(real *w, real *w_co, dom_struct *dom);
/*
 * PURPOSE
 *  Interpolate a Gfz grid (e.g. w-velocity) to a Gcc grid. This routine
 *  interpolates from a Gfz grid with ghost cell buffers to a Gcc grid with
 *  no ghost cell buffers.
 * ARGUMENTS
 *  * w -- device Gfz field with ghost cell buffers
 *  * w_co -- device colocated Gcc field with no ghost cell buffers
 *  * dom -- device domain data
 ******
 */

/****f* bluebottle_kernel/energy_multiply<<<>>>()
 * NAME
 *  energy_multiply<<<>>>()
 * TYPE
 */
__global__ void energy_multiply(real *u_co, real *v_co, real *w_co, real *co,
  dom_struct *dom);
/* PURPOSE
 *  CUDA kernel to multiply the three colocated Gcc velocity fields in place.
 * ARGUMENTS
 *  * u_co -- colocated u-velocity field
 *  * v_co -- colocated v-velocity field
 *  * w_co -- colocated w-velocity field
 *  * co -- colocated result field
 ******
 */

#endif
