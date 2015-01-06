/****h* Bluebottle/vtk
 * NAME
 *  vtk
 * FUNCTION
 *  Write Paraview VTK files.
 ******
 */

#ifndef _VTK_H
#define _VTK_H

/****f* vtk/init_VTK()
 * NAME
 *  init_VTK()
 * USAGE
 */
void init_VTK(void);
/*
 * FUNCTION
 *  Set up the Paraview VTK output PVD file for writing timestep information
 *  about each output file.
 ******
 */

/****f* vtk/out_VTK()
 * NAME
 *  out_VTK()
 * USAGE
 */
void out_VTK(void);
/*
 * FUNCTION
 *  Output ParaView VTK files for the timestep (flow domain and point_particle data).
 ******
 */

/****f* vtk/dom_out_VTK()
 * NAME
 *  dom_out_VTK()
 * USAGE
 */
void dom_out_VTK(void);
/*
 * FUNCTION
 *  Output ParaView VTK file for viewing velocity and pressure fields.  The
 *  velocity field is interpolated from the faces and assigned to a cell.
 ******
 */

/****f* vtk/init_VTK_ghost()
 * NAME
 *  init_VTK_ghost()
 * USAGE
 */
void init_VTK_ghost(void);
/*
 * FUNCTION
 *  Set up the Paraview VTK output PVD file for writing timestep information
 *  about each output file. This version outputs ghost cells.
 ******
 */

/****f* vtk/out_VTK_ghost()
 * NAME
 *  out_VTK_ghost()
 * USAGE
 */
void out_VTK_ghost(void);
/*
 * FUNCTION
 *  Output ParaView VTK files for the timestep (flow domain and point_particle data).
 *  This version outputs ghost cells.
 ******
 */

/****f* vtk/dom_out_VTK_ghost()
 * NAME
 *  dom_out_VTK_ghost()
 * USAGE
 */
void dom_out_VTK_ghost(void);
/*
 * FUNCTION
 *  Output ParaView VTK file for viewing velocity and pressure fields.  The
 *  velocity field is interpolated from the faces and assigned to a cell. This
 *  version outputs ghost cells.
 ******
 */

/****f* vtk/init_VTK_turb()
 * NAME
 *  init_VTK_turb()
 * USAGE
 */
void init_VTK_turb(void);
/*
 * FUNCTION
 *  Set up the Paraview VTK output PVD file for writing timestep information
 *  about each output file.
 ******
 */

/****f* vtk/out_VTK_turb()
 * NAME
 *  out_VTK_turb()
 * USAGE
 */
void out_VTK_turb(void);
/*
 * FUNCTION
 *  Output ParaView VTK files for the timestep (flow domain and point_particle data).
 ******
 */

/****f* vtk/dom_out_VTK_turb()
 * NAME
 *  dom_out_VTK_turb()
 * USAGE
 */
void dom_out_VTK_turb(void);
/*
 * FUNCTION
 *  Output ParaView VTK file for viewing velocity and pressure fields.  The
 *  velocity field is interpolated from the faces and assigned to a cell.
 ******
 */

/****f* vtk/point_out_VTK()
 * NAME
 *  point_out_VTK()
 * USAGE
 */
void point_out_VTK();
/*
 * FUNCTION
 *  Output ParaView VTK file for viewing point_particle data.
 ******
 */

/****f* vtk/quadnodes_out_VTK()
 * NAME
 *  quadnodes_out_VTK()
 * USAGE
 */
void quadnodes_out_VTK();
/*
 * FUNCTION
 *  Output ParaView VTK file for viewing point_particle quadrature node data.
 ******
 */

#endif
