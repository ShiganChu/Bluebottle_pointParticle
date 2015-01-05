/****h* Bluebottle/precursor
 * NAME
 *  precursor
 * FUNCTION
 *  A utility for performing the various data manipulations needed for
 *  transferring precursor flow data to the inflow plane of the experimental
 *  domain.
 ******
 */

#ifndef _PRECURSOR_H
#define _PRECURSOR_H

#include <mpi.h>
#include "bluebottle.h"
#include "scalar.h"

/****f* precursor/expd_init_BC()
 * NAME
 *  expd_init_BC()
 * TYPE
 */
void expd_init_BC(int np);
/*
 * FUNCTION
 *  Communiciate the boundary condition configuration information with the
 *  precursor.
 * ARGUMENTS
 *  * np -- number of processors
 ******
 */

/****f* precursor/expd_update_BC()
 * NAME
 *  expd_update_BC()
 * TYPE
 */
void expd_update_BC(int np, MPI_Status status);
/*
 * FUNCTION
 *  Communicate the boundary condition configuration information with the
 *  precursor.
 * ARGUMENTS
 *  * np -- number of processors
 *  * status -- the MPI_Status for MPI_Recv
 ******
 */

/****f* precursor/expd_compare_dt()
 * NAME
 *  expd_compare_dt()
 * TYPE
 */
void expd_compare_dt(int np, MPI_Status status);
/*
 * FUNCTION
 *  Negotiate with the precursor domain on the appropriate timestep size to
 *  take.
 * ARGUMENTS
 *  * np -- number of processors
 *  * status -- the MPI_Status for MPI_Recv
 ******
 */

/****f* precursor/prec_init_BC()
 * NAME
 *  prec_init_BC()
 * TYPE
 */
void prec_init_BC(int np, int rank, MPI_Status status);
/*
 * FUNCTION
 *  Communicate the boundary condition configuration information with the
 *  experimental domain.
 * ARGUMENTS
 *  * np -- number of processors
 *  * rank -- the MPI rank of this process
 *  * status -- the MPI_Status for MPI_Recv
 ******
 */

/****f* precursor/prec_update_BC()
 * NAME
 *  prec_update_BC()
 * TYPE
 */
void prec_update_BC(int np, int rank, MPI_Status status);
/*
 * FUNCTION
 *  Communicate the boundary condition configuration information with the
 *  experimental domain.
 * ARGUMENTS
 *  * np -- number of processors
 *  * rank -- the MPI rank of this process
 ******
 */

/****f* precursor/prec_send_BC()
 * NAME
 *  prec_send_BC()
 * TYPE
 */
void prec_send_BC(int np, int rank, MPI_Status status);
/*
 * FUNCTION
 *  Communicate the boundary condition to the experimental domain.
 * ARGUMENTS
 *  * np -- number of processors
 *  * rank -- the MPI rank of this process
 *  * status -- the MPI_Status for MPI_Recv
 ******
 */

/****f* precursor/prec_compare_dt()
 * NAME
 *  prec_compare_dt()
 * TYPE
 */
void prec_compare_dt(int np, int rank, MPI_Status status);
/*
 * FUNCTION
 *  Negotiate with the precursor domain on the appropriate timestep size to
 *  take.
 * ARGUMENTS
 *  * np -- number of processors
 *  * rank -- the MPI rank of this process
 *  * status -- the MPI_Status for MPI_Recv
 ******
 */

#endif
