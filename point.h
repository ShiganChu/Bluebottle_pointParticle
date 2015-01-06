/****h* Bluebottle/point_particle
 * NAME
 *  point_particle
 * FUNCTION
 *  Bluebottle point_particle functions.
 ******
 */

#ifndef _point_particle_H
#define _point_particle_H

#include "bluebottle.h"

/****d* point_particle/NNODES
 * NAME
 *  NNODES
 * TYPE
 */
#define NNODES 26
/*
 * PURPOSE
 *  Define the number of nodes used for the Lebedev quadrature scheme.
 ******
 */



/****v* point_particle/
 * NAME
 *  flag_u
 * TYPE
 */
extern int *flag_u;
/* 
 * PURPOSE
 *  Flag x-direction components of velocity field that are set as boundaries.
 ******
 */

/****v* point_particle/_flag_u
 * NAME
 *  _flag_u
 * TYPE
 */
extern int **_flag_u;
/*
 * PURPOSE
 *  CUDA device analog for flag_u.  It contains pointers to arrays containing
 *  the subdomain fields on which each device operates.
 ******
 */

/****v* point_particle/flag_v
 * NAME
 *  flag_v
 * TYPE
 */
extern int *flag_v;
/* 
 * PURPOSE
 *  Flag y-direction components of velocity field that are set as boundaries.
 ******
 */

/****v* point_particle/_flag_v
 * NAME
 *  _flag_v
 * TYPE
 */
extern int **_flag_v;
/*
 * PURPOSE
 *  CUDA device analog for flag_v.  It contains pointers to arrays containing
 *  the subdomain fields on which each device operates.
 ******
 */

/****v* point_particle/flag_w
 * NAME
 *  flag_w
 * TYPE
 */
extern int *flag_w;
/* 
 * PURPOSE
 *  Flag z-direction components of velocity field that are set as boundaries.
 ******
 */

/****v* point_particle/_flag_w
 * NAME
 *  _flag_w
 * TYPE
 */
extern int **_flag_w;
/*
 * PURPOSE
 *  CUDA device analog for flag_w.  It contains pointers to arrays containing
 *  the subdomain fields on which each device operates.
 ******
 */


/****s* point_particle/point_struct
 * NAME
 *  point_struct
 * TYPE
 */
typedef struct point_struct {
  real r;
  real x;
  real y;
  real z;
  real u;
  real v;
  real w;
  real u0;
  real v0;
  real w0;
  real udot;
  real vdot;
  real wdot;

  real ox;
  real oy;
  real oz;
  real ox0;
  real oy0;
  real oz0;
  real oxdot;
  real oydot;
  real ozdot;
  real Fx;
  real Fy;
  real Fz;
  real Lx;
  real Ly;
  real Lz;
/* 
  real aFx;
  real aFy;
  real aFz;
  real aLx;
  real aLy;
  real aLz;
  real kFx;
  real kFy;
  real kFz;
*/
  real iFx;
  real iFy;
  real iFz;
  real iLx;
  real iLy;
  real iLz;
 
  real rho;
 // real mass; //add by shigan 8/9/2014   maybe can change rho to mass 

  real rs;
/*
  real spring_k;
  real spring_x;
  real spring_y;
  real spring_z;
*/
  int rotating;
/*
int i;
int j;
int k;
long int id;
real dt;
*/
} point_struct;
/*
 * PURPOSE
 *  Carry physical information regarding a point_particle.
 * MEMBERS
 *  * r -- the point_particle radius size
 *  * x -- the point_particle location component
 *  * y -- the point_particle location component
 *  * z -- the point_particle location component
 *  * u -- linear velocity in x-direction
 *  * v -- linear velocity in y-direction
 *  * w -- linear velocity in z-direction
 *  * u0 -- initial linear velocity in x-direction
 *  * v0 -- initial linear velocity in y-direction
 *  * w0 -- initial linear velocity in z-direction
 *  * udot -- linear acceleration in x-direction
 *  * vdot -- linear acceleration in y-direction
 *  * wdot -- linear acceleration in z-direction
 *  * axx -- x-component of basis vector initially coincident with x-axis
 *  * axy -- y-component of basis vector initially coincident with x-axis
 *  * axz -- z-component of basis vector initially coincident with x-axis
 *  * ayx -- x-component of basis vector initially coincident with y-axis
 *  * ayy -- y-component of basis vector initially coincident with y-axis
 *  * ayz -- z-component of basis vector initially coincident with y-axis
 *  * azx -- x-component of basis vector initially coincident with z-axis
 *  * azy -- y-component of basis vector initially coincident with z-axis
 *  * azz -- z-component of basis vector initially coincident with z-axis
 *  * ox -- angular velocity in x-direction
 *  * oy -- angular velocity in y-direction
 *  * oz -- angular velocity in z-direction
 *  * oxdot -- angular acceleration in x-direction
 *  * oydot -- angular acceleration in y-direction
 *  * ozdot -- angular acceleration in z-direction
 *  * Fx -- hydrodynamic force in the x-direction
 *  * Fy -- hydrodynamic force in the y-direction
 *  * Fz -- hydrodynamic force in the z-direction
 *  * Lx -- hydrodynamic moment in the x-direction
 *  * Ly -- hydrodynamic moment in the y-direction
 *  * Lz -- hydrodynamic moment in the z-direction
 *  * aFx -- applied force in the x-direction
 *  * aFy -- applied force in the y-direction
 *  * aFz -- applied force in the z-direction
 *  * aLx -- applied moment in the x-direction
 *  * aLy -- applied moment in the y-direction
 *  * aLz -- applied moment in the z-direction
 *  * kFx -- applied spring force in the x-direction
 *  * kFy -- applied spring force in the y-direction
 *  * kFz -- applied spring force in the z-direction
 *  * iFx -- interaction force in the x-direction
 *  * iFy -- interaction force in the y-direction
 *  * iFz -- interaction force in the z-direction
 *  * iLx -- interaction moment in the x-direction
 *  * iLy -- interaction moment in the y-direction
 *  * iLz -- interaction moment in the z-direction
 *  * nodes -- the status of the nodes
 *  * rho -- point_particle density
 *  * cage -- the cage_struct that defines the point_particle within the domain
 *  * order -- the order above which to truncate the Lamb's series solution
 *  * rs -- the radius of integration for scalar products
 *  * ncoeff -- the number of Lamb's coefficients required order truncation
 *  * spring_k -- strength of spring pulling point_particle back to origin
 *  * spring_x -- x location of spring connection
 *  * spring_y -- y location of spring connection
 *  * spring_z -- z location of spring connection
 *  * translating -- 1 if allowed to translate, 0 if not
 *  * rotating -- 1 if allowed to rotate, 0 if not
 ******
 */

/****v* point_particle/npoints
 * NAME
 *  npoints
 * TYPE
 */
extern int npoints;
/*
 * PURPOSE
 *  Define the total number of point_particles in the domain.
 ******
 */

/****v* point_particle/points
 * NAME
 *  points
 * TYPE
 */
extern point_struct *points;
/*
 * PURPOSE
 *  A list of all point_particles.
 ******
 */



/****v* point_particle/_points
 * NAME
 *  _points
 * TYPE
 */
extern point_struct **_points;
/*
 * PURPOSE
 *  CUDA device analog for points.  It contains pointers to arrays containing
 *  the point_particles in domain on which each device operates.  NOTE: for now,
 *  point_particles are designed to function on only one device.
 ******
 */

/****f* point_particle/points_read_input()
 * NAME
 *  points_read_input()
 * USAGE
 */
void points_read_input(int turb);
/*
 * FUNCTION
 *  Read point_particle specfications from point.input.
 * ARGUMENTS
 *  * turb -- boolean operator for removing point_particles from turbulence precursor
 ******
 */

/****f* point_particle/points_show_config()
 * NAME
 *  points_show_config()
 * USAGE
 */
void points_show_config(void);
/*
 * FUNCTION
 *  Write point_particle specifications to screen.
 ******
 */

/****f* point_particle/points_init()
 * NAME
 *  points_init()
 * USAGE
 *
 */
int points_init(void);
/*
 * FUNCTION
 *  Initialize the point_particles (build initial cages) and phase.
 * RESULT
 *  EXIT_SUCCESS if successful, EXIT_FAILURE otherwise.
 ******
 */



/****f* point_particle/flags_reset()
 * NAME
 *  flags_reset()
 * USAGE
 */
void flags_reset(void);
/*
 * FUNCTION
 *  Reinitializes flag arrays to no boundaries (1).
 ******
 */

/****f* point_particle/points_clean()
 * NAME
 *  points_clean()
 * USAGE
 */
void points_clean(void);
/*
 * FUNCTION
 *  Clean up.  Free any allocated host memory.
 ******
 */


#endif
