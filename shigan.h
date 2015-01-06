/* Bluebottle/shigan.h */
#include "bluebottle.h"

#ifndef _SHIGAN_H
#define _SHIGAN_H

/*
extern "C"
{
#include "bluebottle.h"
#include "point.h"
}
*/

//extern "C"
//void cuda_move_points();





//extern "C"
void cgns_point_particles(real dtout);
//void cuda_flow_stress();

__global__ void interpolate_point_pre_shigan(real *p0, real *p,real *pg,int npoints,point_struct *points, dom_struct *dom, real dt0, real dt, BC bc);

__global__ void interpolate_point_vel_shigan(real *u, real *v, real *w,int npoints,real rho_f, real nu,real *ug,real *vg,real *wg, point_struct *points, dom_struct *dom, real dt0, real dt, BC bc);

__global__ void stress_u(real rho_f, real nu, real *u0,real *p,real *p0,real *stress, dom_struct *dom,int *flag_u,real dt,real dt0);

__global__ void stress_v(real rho_f, real nu, real *v0,real *p,real *p0,real *stress, dom_struct *dom,int *flag_v,real dt,real dt0);

__global__ void stress_w(real rho_f, real nu, real *w0,real *p,real *p0,real *stress, dom_struct *dom,int *flag_w,real dt,real dt0);


__global__ void interpolate_point_pre(real *p0, real *p,point_struct *points, dom_struct *dom, real dt0, real dt, BC bc);

__global__ void interpolate_point_vel(real *u, real *v, real *w,int npoints,real rho_f, real nu,real *ug,real *vg,real *wg, point_struct *points, dom_struct *dom, real dt0, real dt, BC bc);


__global__ void drag_points(point_struct *points, int npoints,
real *ug,real *vg,real *wg,
real *lpt_stress_u,real *lpt_stress_v,real *lpt_stress_w,
real rho_f,real mu, g_struct g,gradP_struct gradP,
real C_add,real C_stress,real C_drag);

__global__ void  stress_points(real *ug, real *vg, real *wg,real *ug0, real *vg0, real *wg0,real *conv_ug, real *conv_vg, real *conv_wg,point_struct *points, int npoints,real rho_f, real dt);

__global__ void move_points_b(dom_struct *dom, point_struct *points, int npoints, real dt, real dt0, g_struct g, real rho_f, real ttime,
real C_add,real C_stress,real C_drag);



#endif
