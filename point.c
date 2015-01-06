#include "point.h"

int *flag_u;
int **_flag_u;
int *flag_v;
int **_flag_v;
int *flag_w;
int **_flag_w;

int npoints;
point_struct *points;
point_struct **_points;

void points_read_input(int turb)
{
  int i;  // iterator

  int fret = 0;
  fret = fret; // prevent compiler warning

  // open configuration file for reading
  char fname[FILE_NAME_SIZE];
  sprintf(fname, "%s/input/point.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // read point_point_particle list
  fret = fscanf(infile, "n %d\n", &npoints);
//  if(turb) npoints = 0; // remove point_point_particles from turbulence precursor simulation
  // allocate point_point_particle list
  points = (point_struct*) malloc(npoints * sizeof(point_struct));
  cpumem += npoints * sizeof(point_struct);

  // read npoints point_point_particles
  for(i = 0; i < npoints; i++) {
    fret = fscanf(infile, "\n");
#ifdef DOUBLE
    fret = fscanf(infile, "r %lf\n", &points[i].r);
    fret = fscanf(infile, "(x, y, z) %lf %lf %lf\n",
      &points[i].x, &points[i].y, &points[i].z);
/*
    fret = fscanf(infile, "(aFx, aFy, aFz) %lf %lf %lf\n",
      &points[i].aFx, &points[i].aFy, &points[i].aFz);
    fret = fscanf(infile, "(aLx, aLy, aLz) %lf %lf %lf\n",
      &points[i].aLx, &points[i].aLy, &points[i].aLz);
*/
    fret = fscanf(infile, "rho %lf\n", &points[i].rho);
#else // single precision
    fret = fscanf(infile, "r %f\n", &points[i].r);
    fret = fscanf(infile, "(x, y, z) %f %f %f\n",
      &points[i].x, &points[i].y, &points[i].z);
/*
    fret = fscanf(infile, "(aFx, aFy, aFz) %f %f %f\n",
      &points[i].aFx, &points[i].aFy, &points[i].aFz);
    fret = fscanf(infile, "(aLx, aLy, aLz) %f %f %f\n",
      &points[i].aLx, &points[i].aLy, &points[i].aLz);
*/
    fret = fscanf(infile, "rho %f\n", &points[i].rho);
#endif
    fret = fscanf(infile, "rotating %d\n", &points[i].rotating);
  }

  fclose(infile);
}

void points_show_config(void)
{
  int i;  // iterator

  printf("point_particles:\n");
  for(i = 0; i < npoints; i++) {
    printf("  point_particle %d:\n", i);
    printf("    r = %e\n", points[i].r);
    printf("    (x, y, z) = (%e, %e, %e)\n",
      points[i].x, points[i].y, points[i].z);
    printf("    (u, v, w) = (%e, %e, %e)\n",
      points[i].u, points[i].v, points[i].w);
    printf("    (udot, vdot, wdot) = (%e, %e, %e)\n",
      points[i].udot, points[i].vdot, points[i].wdot);
    printf("    (ox, oy, oz) = (%e, %e, %e)\n",
      points[i].ox, points[i].oy, points[i].oz);
    printf("    (oxdot, oydot, ozdot) = (%e, %e, %e)\n",
      points[i].oxdot, points[i].oydot, points[i].ozdot);
    printf("    (Fx, Fy, Fz) = (%e %e %e)\n",
      points[i].Fx, points[i].Fy, points[i].Fz);
    printf("    (Lx, Ly, Lz) = (%e %e %e)\n",
      points[i].Lx, points[i].Ly, points[i].Lz);
/*
    printf("    (aFx, aFy, aFz) = (%e %e %e)\n",
      points[i].aFx, points[i].aFy, points[i].aFz);
    printf("    (aLx, aLy, aLz) = (%e %e %e)\n",
      points[i].aLx, points[i].aLy, points[i].aLz);
*/
    printf("    rho = %f\n", points[i].rho);
    printf("    rotating = %d\n", points[i].rotating);
  }
}

int points_init(void)
{
  int i;  // iterators
  // allocate and initialize phase array
  flag_u = (int*) malloc(Dom.Gfx.s3b * sizeof(int));
  cpumem += Dom.Gfx.s3b * sizeof(int);
  flag_v = (int*) malloc(Dom.Gfy.s3b * sizeof(int));
  cpumem += Dom.Gfy.s3b * sizeof(int);
  flag_w = (int*) malloc(Dom.Gfz.s3b * sizeof(int));
  cpumem += Dom.Gfz.s3b * sizeof(int);
  flags_reset();

  // TODO: check that the point_point_particles are in valid locations (i.e. exist
  // completely inside the domain and do not overlap)


  for(i = 0; i < npoints; i++) {

    // calculate the number of coefficients needed

    // initialize velocity and acceleration to zero
    points[i].u = 0.;
    points[i].v = 0.;
    points[i].w = 0.;

    points[i].u0 = 0.;
    points[i].v0 = 0.;
    points[i].w0 = 0.;

    points[i].udot = 0.;
    points[i].vdot = 0.;
    points[i].wdot = 0.;
    /* set initial position of point_point_particle reference basis to match the global
     * domain basis */
    points[i].ox = 0.;
    points[i].oy = 0.;
    points[i].oz = 0.;
    points[i].ox0 = 0.;
    points[i].oy0 = 0.;
    points[i].oz0 = 0.;
    points[i].oxdot = 0.;
    points[i].oydot = 0.;
    points[i].ozdot = 0.;

    // initialize the hydrodynamic forces and moments to zero
    points[i].Fx = 0.;
    points[i].Fy = 0.;
    points[i].Fz = 0.;
    points[i].Lx = 0.;
    points[i].Ly = 0.;
    points[i].Lz = 0.;


    // initialize the point_point_particle interaction force to zero
    points[i].iFx = 0.;
    points[i].iFy = 0.;
    points[i].iFz = 0.;
    points[i].iLx = 0.;
    points[i].iLy = 0.;
    points[i].iLz = 0.;

//TODO  change ms initialization in the future
 //   points[i].ms = 4*PI/3.0f *points[i].r*points[i].r*points[i].r*points[i].rho;

    real m=4*PI/3.0f *points[i].r*points[i].r*points[i].r*points[i].rho;
    points[i].ms = sc_init_percent*m;

//printf("\npoint %d %f %f %f %f\n",i,sc_init_percent,m,points[i].r,points[i].rho);

    points[i].x0 =points[i].x;
    points[i].y0 =points[i].y;
    points[i].z0 =points[i].z;
    points[i].ms0 =points[i].ms ;
    points[i].msdot =0;

    //point id
    points[i].id = i+1;
 
   //index of grid that the particle reside in
    points[i].i = 0;
    points[i].j = 0;
    points[i].k = 0;
 
    //point iteration sub-timestep
    points[i].dt = 0.;

  }


  return EXIT_SUCCESS;
}

void flags_reset(void)
{
  int i, j, k;  // iterators
  int C;        // cell location
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        flag_u[C] = 1;  // not a boundary
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        flag_v[C] = 1;  // not a boundary
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        flag_w[C] = 1;  // not a boundary
      }
    }
  }
}



void points_clean(void)
{

  free(points);
  free(flag_u);
  free(flag_v);
  free(flag_w);
}

//TODO remember to change the corresponding part in domain.c when modify point
//point_out_restart is in out_restart of domain.c
//point_in_restart is in in_restart of domain.c
//The read and cgns write of scalar are in recorder.c!!!
