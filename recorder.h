/****h* Bluebottle/recorder
 * NAME
 *  recorder
 * FUNCTION
 *  A utility for recording simulation metrics.
 ******
 */

#ifndef _RECORDER_H
#define _RECORDER_H

#include "bluebottle.h"
#include "scalar.h"
#include "point.h"
#include <cgnslib.h>

/****f* recorder/recorder_read_config()
 * NAME
 *  recorder_read_config()
 * TYPE
 */
void recorder_read_config(void);
/*
 * FUNCTION
 *  Read the record.config file to determine what output to write.
 ******
 */

/****f* recorder/cgns_grid()
 * NAME
 *  cgns_grid()
 * TYPE
 */
void cgns_grid(void);
/*
 * FUNCTION
 *  Write the CGNS grid output file.
 ******
 */

/****f* recorder/cgns_flow_field()
 * NAME
 *  cgns_flow_field()
 * TYPE
 */
void cgns_flow_field(real dtout);
/*
 * FUNCTION
 *  Write the CGNS flow_field output file.
 * ARGUMENTS
 *  * dtout -- the output timestep size
 ******
 */

/****f* recorder/cgns_particles()
 * NAME
 *  cgns_particles()
 * TYPE
 */
void cgns_point_particles(real dtout);
void cgns_scalar_field(real dtout);
/*
 * FUNCTION
 *  Write the CGNS particles output file.
 * ARGUMENTS
 *  * dtout -- the output timestep size
 ******
 */

/****f* recorder/recorder_bicgstab_init()
 * NAME
 *  recorder_bicgstab_init()
 * TYPE
 */
void recorder_bicgstab_init(char *name);
/*
 * FUNCTION
 *  Create the file name for writing and summarize fields to be written for the
 *  BICGSTAB solver.
 * ARGUMENTS
 *  * name -- the name of the file to be written
 ******
 */

/****f* recorder/recorder_lamb_init()
 * NAME
 *  recorder_lamb_init()
 * TYPE
 */
void recorder_lamb_init(char *name);
/*
 * FUNCTION
 *  Create the file name for writing Lamb's coefficients.
 * ARGUMENTS
 *  * name -- the name of the file to be written
 ******
 */

/****f* recorder/recorder_cfl_init()
 * NAME
 *  recorder_cfl_init()
 * TYPE
 */
void recorder_cfl_init(char *name);
/*
 * FUNCTION
 *  Create the file name for writing and summarize fields to be written for the
 *  CFL condition
 * ARGUMENTS
 *  * name -- the name of the file to be written
 ******
 */

/****f* recorder/recorder_bicgstab()
 * NAME
 *  recorder_bicgstab()
 * TYPE
 */
void recorder_bicgstab(char *name, int niter, real resid);
/*
 * FUNCTION 
 *  Write out BICGSTAB solver information to file name.
 * ARGUMENTS
 *  * name -- the name of the file to which to write
 *  * niter -- the number of iterations to convergence
 *  * resid -- the residual at convergence
 ******
 */

/****f* recorder/recorder_lamb()
 * NAME
 *  recorder_lamb()
 * TYPE
 */
void recorder_lamb(char *name, int iter);
/*
 * FUNCTION 
 *  Write out Lamb's coefficients to file.
 * ARGUMENTS
 *  * name -- the name of the file to which to write
 *  * iter -- the number of the iteration
 ******
 */

/****f* recorder/recorder_cfl()
 * NAME
 *  recorder_cfl()
 * TYPE
 */
void recorder_cfl(char *name, real cfl);
/*
 * FUNCTION 
 *  Write out BICGSTAB solver information to file name.
 * ARGUMENTS
 *  * name -- the name of the file to which to write
 *  * cfl -- the cfl number
 ******
 */

#endif
