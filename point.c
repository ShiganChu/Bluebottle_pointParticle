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

void points_read_input(void)
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
  if(npoints<=0) return;
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

//This subroutine delete the old scalar, and inject new scalar into the flow field based on point.config&&scalar.config
void scalar_inject(void)
{

//free scalar on device and host
      cuda_scalar_free();
      scalar_clean();

//read and initialize scalar
      scalar_read_input();
    // initialize the scalar 
      int scalar_init_flag = scalar_init();
    /*
  if(scalar_init_flag == EXIT_FAILURE) {
        printf("\nThe initial scalar configuration is not allowed.\n");
        return EXIT_FAILURE;
      }

*/
      point_ms_init();
      cuda_point_pull();

//malloc device memory of scalar and point, and push host data to device
      cuda_scalar_malloc();
      cuda_scalar_push();
      
//Initialize time again
	ttime=0.f;
	dt0=0.f;	
	dt=0.f;
//write initial field 
          cuda_dom_pull();
          if(rec_flow_field_dt > 0) {
            cgns_grid();
            cgns_flow_field(rec_flow_field_dt);
            rec_flow_field_stepnum_out++;
//printf("\nrec_flow %d\n", rec_flow_field_stepnum_out);
          }
          if(rec_point_particle_dt > 0) {
            cgns_point_particles(rec_point_particle_dt);
            rec_point_particle_stepnum_out++;
          }
      
         if(rec_scalar_field_dt > 0) {
            cgns_scalar_field(rec_scalar_field_dt);
            rec_scalar_stepnum_out++;
          }
        

}





//This subroutine delete the old particle, and inject new particle into the flow field based on point.config&&scalar.config
void points_inject(void)
{

//free points on device and host
      cuda_point_free();
      points_clean();


//read and initialize points	
      points_read_input();
      int points_init_flag = points_init();
      /*
	if(points_init_flag == EXIT_FAILURE) {
        printf("\nThe initial point_particle configuration is not allowed.\n");
        return EXIT_FAILURE;
      }
      */
      cuda_point_malloc();
      cuda_point_push();

if(npoints>0) cuda_flow_stress();

//The domain velocity has already been pushed to device
//Match device point velocity with flow field based on point position, which is copied from host
      match_point_vel_with_flow();
//pull the new point infomation to host
       cuda_point_pull();
//Initialize time again
	ttime=0.f;
	dt0=0.f;
	dt=0.f;
//write initial field 
          cuda_dom_pull();
          if(rec_flow_field_dt > 0) {
            cgns_grid();
            cgns_flow_field(rec_flow_field_dt);
            rec_flow_field_stepnum_out++;
//printf("\nrec_flow %d\n", rec_flow_field_stepnum_out);
          }
          if(rec_point_particle_dt > 0) {
            cgns_point_particles(rec_point_particle_dt);
            rec_point_particle_stepnum_out++;
          }
      
         if(rec_scalar_field_dt > 0) {
            cgns_scalar_field(rec_scalar_field_dt);
            rec_scalar_stepnum_out++;
          }
        

}



//This subroutine delete the old particle&scalar, and inject new particle and scalar into the flow field based on point.config&&scalar.config
void points_scalar_inject(void)
{

//free points on device and host
      cuda_point_free();
      points_clean();
//free scalar on device and host
      cuda_scalar_free();
      scalar_clean();

//read and initialize points	
      points_read_input();
      int points_init_flag = points_init();
      fflush(stdout);
/*
      if(points_init_flag == EXIT_FAILURE) {
        printf("\nThe initial point_particle configuration is not allowed.\n");
        return EXIT_FAILURE;
      }
*/
//read and initialize scalar
      scalar_read_input();
    // initialize the scalar 
      int scalar_init_flag = scalar_init();
/*
      fflush(stdout);
      if(scalar_init_flag == EXIT_FAILURE) {
        printf("\nThe initial scalar configuration is not allowed.\n");
        return EXIT_FAILURE;
      }
*/
//malloc device memory of scalar and point, and push host data to device
      cuda_scalar_malloc();
      cuda_scalar_push();
      
      cuda_point_malloc();
      cuda_point_push();

if(npoints>0) cuda_flow_stress();

//The domain velocity has already been pushed to device
//Match device point velocity with flow field based on point position, which is copied from host
      match_point_vel_with_flow();
//pull the new point infomation to host
          cuda_point_pull();
//Initialize time again
	ttime=0.f;
	dt0=0.f;
	dt=0.f;
//write initial field 
          cuda_dom_pull();
          if(rec_flow_field_dt > 0) {
            cgns_grid();
            cgns_flow_field(rec_flow_field_dt);
            rec_flow_field_stepnum_out++;
//printf("\nrec_flow %d\n", rec_flow_field_stepnum_out);
          }
          if(rec_point_particle_dt > 0) {
            cgns_point_particles(rec_point_particle_dt);
            rec_point_particle_stepnum_out++;
          }
      
         if(rec_scalar_field_dt > 0) {
            cgns_scalar_field(rec_scalar_field_dt);
            rec_scalar_stepnum_out++;
          }
        

}


//This subroutine delete the old particle&scalar, and inject new particle and scalar into the flow field based on point.config&&scalar.config
void bubble_scalar_inject(void)
{

//free points on device and host
      cuda_point_free();
      points_clean();
//free scalar on device and host
      cuda_scalar_free();
      scalar_clean();

//read and initialize points	
      points_read_input();
      int points_init_flag = points_init();
      fflush(stdout);
 
//read and initialize scalar
      scalar_read_input();
    // initialize the scalar 
      int scalar_init_flag = scalar_init();
 
//malloc device memory of scalar and point, and push host data to device
      cuda_scalar_malloc();
      cuda_scalar_push();
      
      cuda_point_malloc();
      cuda_point_push();

if(npoints>0) cuda_flow_stress();

//The domain velocity has already been pushed to device
//Match device point velocity with flow field based on point position, which is copied from host
      match_bubble_vel_with_flow();
//pull the new point infomation to host
          cuda_point_pull();
//Initialize time again
	ttime=0.f;
	dt0=0.f;
	dt=0.f;
//write initial field 
          cuda_dom_pull();
          if(rec_flow_field_dt > 0) {
            cgns_grid();
            cgns_flow_field(rec_flow_field_dt);
            rec_flow_field_stepnum_out++;
//printf("\nrec_flow %d\n", rec_flow_field_stepnum_out);
          }
          if(rec_point_particle_dt > 0) {
            cgns_point_particles(rec_point_particle_dt);
            rec_point_particle_stepnum_out++;
          }
      
         if(rec_scalar_field_dt > 0) {
            cgns_scalar_field(rec_scalar_field_dt);
            rec_scalar_stepnum_out++;
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

//To make sure that the insoluble and soluble together make the real oil density, rho is the density of both initialy. As time goes, the soluble component will decrease in density.  mp=ms+m_is|t=0  =rho*Vp*(1+sc_init_percent)=rho*Vp'
//We conclude rp'=rp*(1.f+sc_init_percent)^(1/3)
//There are two ways, the 1st one is decrease insoluble density; the other is increase total mass
//    points[i].r =points[i].r*powf((1.f+sc_init_percent),1/3.f);
//    points[i].rho =points[i].rho*(1.f-sc_init_percent);

//Keep the density to be a constant

//printf("\npoint %d %f %f %f %f\n",i,sc_init_percent,m,points[i].r,points[i].rho);

    points[i].x0 =points[i].x;
    points[i].y0 =points[i].y;
    points[i].z0 =points[i].z;
    points[i].ms0 =points[i].ms ;
    points[i].msdot =0;
    points[i].Nu =0;

    //point id
    points[i].id = i+1;
 
   //index of grid that the particle reside in
    points[i].i = 0;
    points[i].j = 0;
    points[i].k = 0;
 
    //point iteration sub-timestep
    points[i].dt = points[i].rho *2.f*points[i].r*points[i].r/(9.0f*mu);

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
