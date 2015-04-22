#include "recorder.h"

void recorder_read_config(void)
{
  int i, j;  // iterator

  char *fret = malloc(sizeof(char) * CHAR_BUF_SIZE);
  fret = fret;  // prevent compiler warning

  // initialize recorder values
  // if rec_X_dt < 0, the output is off
  // otherwise, the fields are all boolean values with 0 = false and 1 = true
  rec_flow_field_dt = -1;
  rec_flow_field_vel = 0;
  rec_flow_field_p = 0;
  rec_flow_field_phase = 0;

  rec_paraview_dt = -1;

  rec_point_particle_dt = -1;
  rec_point_particle_pos = 0;
  rec_point_particle_a = 0;
  rec_point_particle_vel = 0;
  rec_point_particle_omega = 0;
  rec_point_particle_force = 0;
  rec_point_particle_moment = 0;

  rec_restart_dt = -1;
  rec_precursor_dt = -1;
  rec_restart_stop = 0;

//  rec_scalar_field=0;
  rec_scalar_field_dt=-1;

  // read the config file
  char fname[FILE_NAME_SIZE];
  sprintf(fname, "%s/input/record.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  char buf[CHAR_BUF_SIZE];  // character read buffer

  /** list of recognized output configurations **/
  char ***configs;
//  int nconfigs = 6;
  int nconfigs = 7;
  int *nconfigsn = (int*) malloc(nconfigs * sizeof(int));
  // cpumem += nconfigs * sizeof(int);
  configs = (char***) malloc(nconfigs * sizeof(char**));
  // cpumum += nconfigs * sizeof(char**);
  int N;
  // flow field info
  N = 0;
  nconfigsn[N] = 3;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  // cpumem += nconfigsn[N] * sizeof(char*);
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    // cpumem += CHAR_BUF_SIZE * sizeof(char);
  }
  sprintf(configs[N][0], "FLOW_FIELD");
  sprintf(configs[N][1], "velocity");
  sprintf(configs[N][2], "pressure");
  // ParaView info
  N = 1;
  nconfigsn[N] = 1;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  // cpumem += nconfigsn[N] * sizeof(char*);
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    // cpumem += CHAR_BUF_SIZE * sizeof(char);
  }
  sprintf(configs[N][0], "PARAVIEW");
  // point_particle info
  N = 2;
  nconfigsn[N] = 2;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  // cpumem += nconfigs[N] * sizeof(char*);
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    // cpumem += CHAR_BUF_SIZE * sizeof(char);
  }
  sprintf(configs[N][0], "POINT_PARTICLE");
  sprintf(configs[N][1], "position");
  // restart info
  N = 3;
  nconfigsn[N] = 1;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  // cpumem += nconfigs[N] * sizeof(char*);
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    // cpumem += CHAR_BUF_SIZE * sizeof(char);
  }
  sprintf(configs[N][0], "RESTART");
  // precursor info
  N = 4;
  nconfigsn[N] = 1;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
  }
  sprintf(configs[N][0], "PRECURSOR");

 // restart_stop info
  N = 5;
  nconfigsn[N] = 1;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  // cpumem += nconfigs[N] * sizeof(char*);
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    // cpumem += CHAR_BUF_SIZE * sizeof(char);
  }
  sprintf(configs[N][0], "RESTART_STOP");

 // scalar info
  N = 6;
  nconfigsn[N] = 1;
  configs[N] = (char**) malloc(nconfigsn[N] * sizeof(char*));
  // cpumem += nconfigs[N] * sizeof(char*);
  for(i = 0; i < nconfigsn[N]; i++) {
    configs[N][i] = (char*) malloc(CHAR_BUF_SIZE * sizeof(char));
    // cpumem += CHAR_BUF_SIZE * sizeof(char);
  }
  sprintf(configs[N][0], "SCALAR");


  // read config file
  char bufconfig[CHAR_BUF_SIZE];
  while(fgets(buf, CHAR_BUF_SIZE, infile) != NULL) {
    // compare first line with predefined configurations
    int foundi = 0;
    
  // remove newline character
      real dtout = 0;
      sscanf(buf, "%s %lf\n", bufconfig, &dtout);

//      printf("\ni %s %f\n",bufconfig, dtout);
//      fflush(stdout);

  for(i = 0; i < nconfigs; i++) {	

 //     printf("\ni %s %f %d\n",bufconfig, dtout,i);
 //     fflush(stdout);
   
   // check config list
      if(strcmp(configs[i][0], bufconfig) == 0) {
        // found: continue reading config
        foundi = 1;
        switch(i) {
          case 0: // FLOW_FIELD
            rec_flow_field_dt = dtout;
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
            while(strcmp("\n", buf) != 0) {
              int foundj = 0;
              // remove trailing newline
              int ln = strlen(buf) - 1;
              if(buf[ln] == '\n') buf[ln] = '\0';
              // check options
              for(j = 1; j < nconfigsn[i]; j++) {
                if(strcmp(configs[i][j], buf) == 0) {
                  foundj = 1;
                  switch(j) {
                    case 1: // velocity
                      rec_flow_field_vel = 1;
                      break;
                    case 2: // pressure
                      rec_flow_field_p = 1;
                      break;
                    case 3: // phase
                      rec_flow_field_phase = 1;
                      break;
                    default:
                      printf("UNRECOGNIZED OPTION\n");
                  }
                }
              }
              if(!foundj) {
                fprintf(stderr, "Unrecognized record.config ");
                fprintf(stderr, "%s option %s\n", bufconfig, buf);

                // clean up
                fclose(infile);
                for(i = 0; i < nconfigs; i++) {
                  for(j = 0; j < nconfigsn[i]; j++) {
                    free(configs[i][j]);
                  }
                  free(configs[i]);
                }
                free(configs);
                free(nconfigsn);

                exit(EXIT_FAILURE);
              }
              // read next line
              fret = fgets(buf, CHAR_BUF_SIZE, infile);
            }
            break;
          case 1: // PARAVIEW
            rec_paraview_dt = dtout;
            // read next line
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
            break;
          case 2: // point_particle
            rec_point_particle_dt = dtout;
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
            while(strcmp("\n", buf) != 0) {
              int foundj = 0;
              // remove trailing newline
              int ln = strlen(buf) - 1;
              if(buf[ln] == '\n') buf[ln] = '\0';
              // check options
              for(j = 1; j < nconfigsn[i]; j++) {
                if(strcmp(configs[i][j], buf) == 0) {
                  foundj = 1;
                  switch(j) {
                    case 1: // position
                      rec_point_particle_pos = 1;
                      break;
                    case 2: // radius
                      rec_point_particle_a = 1;
                      break;
                    case 3: // velocity
                      rec_point_particle_vel = 1;
                      break;
                    case 4: // angular velocity
                      rec_point_particle_omega = 1;
                      break;
                    case 5: // hydrodynamic force
                      rec_point_particle_force = 1;
                      break;
                    case 6: // hydrodynamic moment
                      rec_point_particle_moment = 1;
                      break;
                    default:
                      printf("UNRECOGNIZED OPTION\n");
                  }
                }
              }
              if(!foundj) {
                fprintf(stderr, "Unrecognized record.config ");
                fprintf(stderr, "%s option %s\n", bufconfig, buf);

                // clean up
                fclose(infile);
                for(i = 0; i < nconfigs; i++) {
                  for(j = 0; j < nconfigsn[i]; j++) {
                    free(configs[i][j]);
                  }
                  free(configs[i]);
                }
                free(configs);
                free(nconfigsn);

                exit(EXIT_FAILURE);
              }
              // read next line
              fret = fgets(buf, CHAR_BUF_SIZE, infile);
            }
            break;
          case 3: // RESTART
            rec_restart_dt = dtout;
            // read next line
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
            break;
          case 4: // PRECURSOR
            rec_precursor_dt = dtout;
            // read next line
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
            break;
	  case 5: // RESTART_STOP
            rec_restart_stop = (int)dtout;
            // read next line
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
            break;
          case 6: // SCALAR_FIELD
// printf("\nscalar dtout case6 %f %f\n",rec_scalar_field_dt,dtout);
//      fflush(stdout);
            rec_scalar_field_dt = dtout;
            // read next line
            fret = fgets(buf, CHAR_BUF_SIZE, infile);
	    break;
	  default:
            printf("UNRECOGNIZED TYPE\n");
        }
      }
    }
    // if we've parsed the entire list but don't find anything: error
    if(!foundi) {
      fprintf(stderr, "Unrecognized record.config type %s\n", bufconfig);

      // clean up
      fclose(infile);
      for(i = 0; i < nconfigs; i++) {
        for(j = 0; j < nconfigsn[i]; j++) {
          free(configs[i][j]);
        }
        free(configs[i]);
      }
      free(configs);
      free(nconfigsn);

      exit(EXIT_FAILURE);
    }
  }

  // clean up
  fclose(infile);
  for(i = 0; i < nconfigs; i++) {
    for(j = 0; j < nconfigsn[i]; j++) {
      free(configs[i][j]);
    }
    free(configs[i]);
  }
  free(configs);
  free(nconfigsn);
}

void cgns_grid(void)
{
  // create the file
  char fname[FILE_NAME_SIZE];
  sprintf(fname, "%s/output/%s", ROOT_DIR, "grid.cgns");
  int fn;
  int bn;
  int zn;
  int gn;
  int cn;
  cg_open(fname, CG_MODE_WRITE, &fn);
  cg_base_write(fn, "Base", 3, 3, &bn);
  cgsize_t size[9];
  size[0] = Dom.xn+1; // cells -> vertices
  size[1] = Dom.yn+1;
  size[2] = Dom.zn+1;
  size[3] = Dom.xn;
  size[4] = Dom.yn;
  size[5] = Dom.zn;
  size[6] = 0;
  size[7] = 0;
  size[8] = 0;
  cg_zone_write(fn, bn, "Zone0", size, Structured, &zn);
  cg_grid_write(fn, bn, zn, "GridCoordinates", &gn);

  real *x = malloc((Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real));
  // cpumem = (Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real);
  real *y = malloc((Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real));
  // cpumem = (Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real);
  real *z = malloc((Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real));
  // cpumem = (Dom.xn+1)*(Dom.yn+1)*(Dom.zn+1) * sizeof(real);
  for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke+1; k++) {
    for(int j = Dom.Gcc.js; j < Dom.Gcc.je+1; j++) {
      for(int i = Dom.Gcc.is; i < Dom.Gcc.ie+1; i++) {
        int C = (i-1) + (j-1)*(Dom.xn+1) + (k-1)*(Dom.xn+1)*(Dom.yn+1);
        x[C] = Dom.xs + (i-1)*Dom.dx;
        y[C] = Dom.ys + (j-1)*Dom.dy;
        z[C] = Dom.zs + (k-1)*Dom.dz;
      }
    }
  }

  cg_coord_write(fn, bn, zn, RealDouble, "CoordinateX", x, &cn);
  cg_coord_write(fn, bn, zn, RealDouble, "CoordinateY", y, &cn);
  cg_coord_write(fn, bn, zn, RealDouble, "CoordinateZ", z, &cn);

  free(x);
  free(y);
  free(z);

  cg_close(fn);
}





void cgns_flow_field(real dtout)
{
  // create the solution file
  char fname[FILE_NAME_SIZE]="";
  char fname2[FILE_NAME_SIZE]="";
  char fnameall[FILE_NAME_SIZE]="";
  char fnameall2[FILE_NAME_SIZE]="";
  char gname[FILE_NAME_SIZE]="";
  char gnameall[FILE_NAME_SIZE]="";
  real tout = ttime; //  = rec_flow_field_stepnum_out * dtout;
  char format[CHAR_BUF_SIZE]="";
  char snodename[CHAR_BUF_SIZE]="";
  char snodenameall[CHAR_BUF_SIZE]="";
  int sigfigs = ceil(log10(1. / dtout));
  if(sigfigs < 1) sigfigs = 1;
  sprintf(format, "%%.%df", sigfigs);
  sprintf(fname2, "flow-%s.cgns", format);
  sprintf(fnameall2, "%s/output/flow-%s.cgns", ROOT_DIR, format);
  sprintf(snodename, "Solution-");
  sprintf(snodenameall, "/Base/Zone0/Solution-");
  sprintf(snodename, "%s%s", snodename, format);//Some problem is with this line??
  sprintf(snodenameall, "%s%s", snodenameall, format);
  sprintf(fname, fname2, tout);
  sprintf(fnameall, fnameall2, tout);
  sprintf(snodename, snodename, tout);
  sprintf(snodenameall, snodenameall, tout);
  sprintf(gname, "grid.cgns");
  sprintf(gnameall, "%s/output/%s", ROOT_DIR, "grid.cgns");
  int fn;
  int bn;
  int zn;
  int sn;
  int fnpress;
  int fnu;
  int fnv;
  int fnw;
  cg_open(fnameall, CG_MODE_WRITE, &fn);
  cg_base_write(fn, "Base", 3, 3, &bn);
  cgsize_t size[9];
  size[0] = Dom.xn+1; // cells -> vertices
  size[1] = Dom.yn+1;
  size[2] = Dom.zn+1;
  size[3] = Dom.xn;
  size[4] = Dom.yn;
  size[5] = Dom.zn;
  size[6] = 0;
  size[7] = 0;
  size[8] = 0;
  cg_zone_write(fn, bn, "Zone0", size, Structured, &zn);
  cg_goto(fn, bn, "Zone_t", zn, "end");
  // check that grid.cgns exists
  /*int fng;
  if(cg_open(gnameall, CG_MODE_READ, &fng) != 0) {
    fprintf(stderr, "CGNS flow field write failure: no grid.cgns\n");
    exit(EXIT_FAILURE);
  } else {
    cg_close(fng);
  }
    cg_close(fng);
*/
  
  cg_link_write("GridCoordinates", gname, "Base/Zone0/GridCoordinates");

  cg_sol_write(fn, bn, zn, "Solution", CellCenter, &sn);
  real *pout = malloc(Dom.Gcc.s3 * sizeof(real));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
    for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
      for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        pout[C] = p[CC];

      }
    }
  }
  cg_field_write(fn, bn, zn, sn, RealDouble, "Pressure", pout, &fnpress);

  real *uout = malloc(Dom.Gcc.s3 * sizeof(real));

  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(int j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(int i = Dom.Gfx.is; i < Dom.Gfx.ie-1; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC0 = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        int CC1 = (i+1) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        uout[C] = 0.5*(u[CC1] + u[CC0]);
//        fxout[C] = 0.5*(f_x[CC1] + f_x[CC0]);
      }
    }
  }
  cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityX", uout, &fnu);

  real *vout = malloc(Dom.Gcc.s3 * sizeof(real));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(int j = Dom.Gfy.js; j < Dom.Gfy.je-1; j++) {
      for(int i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC0 = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        int CC1 = i + (j+1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        vout[C] = 0.5*(v[CC1] + v[CC0]);
 //       fyout[C] = 0.5*(f_y[CC1] + f_y[CC0]);

     }
    }
  }
  cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityY", vout, &fnv);

  real *wout = malloc(Dom.Gcc.s3 * sizeof(real));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfz.ks; k < Dom.Gfz.ke-1; k++) {
    for(int j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(int i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC0 = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        int CC1 = i + j*Dom.Gfz.s1b + (k+1)*Dom.Gfz.s2b;
        wout[C] = 0.5*(w[CC1] + w[CC0]);
  //      fzout[C] = 0.5*(f_z[CC1] + f_z[CC0]);
      }
    }
  }
  cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityZ", wout, &fnw);
/*
  real *phaseout = malloc(Dom.Gcc.s3 * sizeof(real));
  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
    for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
      for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        phaseout[C] = phase[CC];
      }
    }
  }
  cg_field_write(fn, bn, zn, sn, Integer, "Phase", phaseout, &fnpress);
*/

  real *fxout = malloc(Dom.Gcc.s3 * sizeof(real));
  real *fyout = malloc(Dom.Gcc.s3 * sizeof(real));
  real *fzout = malloc(Dom.Gcc.s3 * sizeof(real));

  real *fx_x = malloc(Dom.Gfx.s3 * sizeof(real));
  real *fy_y = malloc(Dom.Gfy.s3 * sizeof(real));
  real *fz_z = malloc(Dom.Gfz.s3 * sizeof(real));

  real *u_x = malloc(Dom.Gfx.s3 * sizeof(real));
  real *v_y = malloc(Dom.Gfy.s3 * sizeof(real));
  real *w_z = malloc(Dom.Gfz.s3 * sizeof(real));


if(OUT_SIMPLE<=0)
{
// cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(int j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(int i = Dom.Gfx.is; i < Dom.Gfx.ie-1; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC0 = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        int CC1 = (i+1) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        fxout[C] = 0.5*(f_x[CC1] + f_x[CC0]);
      }
    }
  }
 
   // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(int j = Dom.Gfy.js; j < Dom.Gfy.je-1; j++) {
      for(int i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC0 = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        int CC1 = i + (j+1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        fyout[C] = 0.5*(f_y[CC1] + f_y[CC0]);

      }
    }
  }
  
   // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfz.ks; k < Dom.Gfz.ke-1; k++) {
    for(int j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(int i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC0 = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        int CC1 = i + j*Dom.Gfz.s1b + (k+1)*Dom.Gfz.s2b;
         fzout[C] = 0.5*(f_z[CC1] + f_z[CC0]);
/*
if(fabs(f_z[CC1])>EPSILON||fabs(f_z[CC0])>EPSILON) 
{printf("\nwrite %f %f %f\n",fzout[C],f_z[CC1],f_z[CC0]);
fflush(stdout);}
*/
      }
    }
  }

 real tag_fz=0.f;

 real tag_w=0.f;


//Write in different coordinate system
 for(int k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(int j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(int i = Dom.Gfx.is; i < Dom.Gfx.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        fx_x[C] = f_x[CC];
        u_x[C] = u[CC];
      }
    }
  }
 
   // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(int j = Dom.Gfy.js; j < Dom.Gfy.je; j++) {
      for(int i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        fy_y[C] = f_y[CC];
        v_y[C] = v[CC];
      }
    }
  }
  
   // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gfz.ks; k < Dom.Gfz.ke; k++) {
    for(int j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(int i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        fz_z[C] = f_z[CC];
        w_z[C] = w[CC];
tag_fz +=f_z[CC];
tag_w +=w[CC];
/*
int Ct=10 + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
if(fabs(f_z[CC])>fabs(f_z[Ct])+EPSILON) printf("\nForce and velocity: %f %f %f %f %d %d %d \n",f_z[CC]-f_z[Ct],w[CC],f_z[Ct],w[Ct],i,j,k);
fflush(stdout);
*/
      }
    }
  }
/*
real Vcell=powf(Dom.dx,3);
printf("\nSum of force and velocity: %f %f %f\n",tag_fz*Vcell,tag_w*Vcell,Vcell);
fflush(stdout);
*/
 cg_field_write(fn, bn, zn, sn, RealDouble, "f_x", fxout, &fnw);
 cg_field_write(fn, bn, zn, sn, RealDouble, "f_y", fyout, &fnw);
 cg_field_write(fn, bn, zn, sn, RealDouble, "f_z", fzout, &fnw);

 cg_field_write(fn, bn, zn, sn, RealDouble, "fx_x", fx_x, &fnw);
 cg_field_write(fn, bn, zn, sn, RealDouble, "fy_y", fy_y, &fnw);
 cg_field_write(fn, bn, zn, sn, RealDouble, "fz_z", fz_z, &fnw);

 cg_field_write(fn, bn, zn, sn, RealDouble, "u_x", u_x, &fnw);
 cg_field_write(fn, bn, zn, sn, RealDouble, "v_y", v_y, &fnw);
 cg_field_write(fn, bn, zn, sn, RealDouble, "w_z", w_z, &fnw);

}

  cg_user_data_write("Etc");
  cg_goto(fn, bn, "Zone_t", zn, "Etc", 0, "end");
  cgsize_t *N = malloc(sizeof(cgsize_t));
  N[0] = 1;
  cg_array_write("Time", RealDouble, 1, N, &ttime);
  free(N);

  cg_close(fn);
  free(pout);
  free(uout);
  free(vout);
  free(wout);

  free(fxout);
  free(fyout);
  free(fzout);
//  free(phaseout);

  free(fx_x);
  free(fy_y);
  free(fz_z);
  	
  free(u_x);
  free(v_y);
  free(w_z);

}



#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
//add By shigan
void check_output_dir(void)
{

char fnameOut[FILE_NAME_SIZE]="";
char fnameRec[FILE_NAME_SIZE]="";

sprintf(fnameOut, "%s/output", ROOT_DIR);
sprintf(fnameRec, "%s/record", ROOT_DIR);


//int errOut=mkdir(fnameOut)
struct stat s1={0};
if(stat(fnameOut, &s1)==-1) {mkdir(fnameOut,0700);}
struct stat s2={0};
if(stat(fnameRec, &s2)==-1) {mkdir(fnameRec,0700);}


/*
int err = stat("/path/to/possible_dir", &s);
if(-1 == err) {
    if(ENOENT == errno) {
//         does not exist 
    } else {
        perror("stat");
        exit(1);
    }
} else {
    if(S_ISDIR(s.st_mode)) {
        // it's a dir 
    } else {
        // exists but is no dir 
    }
  }
}
*/
}




//add by shigan
//TODO make sure the read and write of scalar are correct!!
void cgns_scalar_field(real dtout)
{
  // create the solution file
  char fname[FILE_NAME_SIZE]="";
  char fname2[FILE_NAME_SIZE]="";
  char fnameall[FILE_NAME_SIZE]="";
  char fnameall2[FILE_NAME_SIZE]="";
  char gname[FILE_NAME_SIZE]="";
  char gnameall[FILE_NAME_SIZE]="";
  real tout = ttime; //  = rec_scalar_field_stepnum_out * dtout;
  char format[CHAR_BUF_SIZE]="";
  char snodename[CHAR_BUF_SIZE]="";
  char snodenameall[CHAR_BUF_SIZE]="";
  int sigfigs = ceil(log10(1. / dtout));
  if(sigfigs < 1) sigfigs = 1;
  sprintf(format, "%%.%df", sigfigs);
  sprintf(fname2, "scalar-%s.cgns", format);
  sprintf(fnameall2, "%s/output/scalar-%s.cgns", ROOT_DIR, format);
  sprintf(snodename, "Solution-");
  sprintf(snodenameall, "/Base/Zone0/Solution-");
  sprintf(snodename, "%s%s", snodename, format);
  sprintf(snodenameall, "%s%s", snodenameall, format);
  
  sprintf(fname, fname2, tout);
  sprintf(fnameall, fnameall2, tout);
  sprintf(snodename, snodename, tout);
  sprintf(snodenameall, snodenameall, tout);
  sprintf(gname, "grid.cgns");
  sprintf(gnameall, "%s/output/%s", ROOT_DIR, "grid.cgns");
  int fn;
  int bn;
  int zn;
  int sn;
  int fnpress;
  cg_open(fnameall, CG_MODE_WRITE, &fn);
  cg_base_write(fn, "Base", 3, 3, &bn);
  cgsize_t size[9];
  size[0] = Dom.xn+1; // cells -> vertices
  size[1] = Dom.yn+1;
  size[2] = Dom.zn+1;
  size[3] = Dom.xn;
  size[4] = Dom.yn;
  size[5] = Dom.zn;
  size[6] = 0;
  size[7] = 0;
  size[8] = 0;
  cg_zone_write(fn, bn, "Zone0", size, Structured, &zn);
  cg_goto(fn, bn, "Zone_t", zn, "end");
  
  cg_link_write("GridCoordinates", gname, "Base/Zone0/GridCoordinates");

  cg_sol_write(fn, bn, zn, "Solution", CellCenter, &sn);
  real *sc_out = malloc(Dom.Gcc.s3 * sizeof(real));
  real *diff_out = malloc(Dom.Gcc.s3 * sizeof(real));
  real *conv_out = malloc(Dom.Gcc.s3 * sizeof(real));
  real *scSrc_out = malloc(Dom.Gcc.s3 * sizeof(real));
  real *epsp_out = malloc(Dom.Gcc.s3 * sizeof(real));

//Vorticity
  real *omegaX_out = malloc(Dom.Gcc.s3 * sizeof(real));
  real *omegaY_out = malloc(Dom.Gcc.s3 * sizeof(real));
  real *omegaZ_out = malloc(Dom.Gcc.s3 * sizeof(real));

  // cpumem += Dom.Gcc.s3 * sizeof(real);
  for(int k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
    for(int j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
      for(int i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
        int C = (i-DOM_BUF) + (j-DOM_BUF)*Dom.Gcc.s1 + (k-DOM_BUF)*Dom.Gcc.s2;
        int CC = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        sc_out[C] = sc[CC];
        conv_out[C] = conv_sc[CC];
        diff_out[C] = diff_sc[CC];
        scSrc_out[C] = scSrc[CC];
	epsp_out[C]=epsp[CC];

	omegaX_out[C]=omega_x[CC];
	omegaY_out[C]=omega_y[CC];
	omegaZ_out[C]=omega_z[CC];

      }
    }
  }

  cg_field_write(fn, bn, zn, sn, RealDouble, "SourceScalar", scSrc_out, &fnpress);
  cg_field_write(fn, bn, zn, sn, RealDouble, "Scalar", sc_out, &fnpress);
  cg_field_write(fn, bn, zn, sn, RealDouble, "ConvScalar", conv_out, &fnpress);
  cg_field_write(fn, bn, zn, sn, RealDouble, "DiffScalar", diff_out, &fnpress);
  cg_field_write(fn, bn, zn, sn, RealDouble, "EPSP", epsp_out, &fnpress);

  cg_field_write(fn, bn, zn, sn, RealDouble, "OmegaX", omegaX_out, &fnpress);
  cg_field_write(fn, bn, zn, sn, RealDouble, "OmegaY", omegaY_out, &fnpress);
  cg_field_write(fn, bn, zn, sn, RealDouble, "OmegaZ", omegaZ_out, &fnpress);

  cg_user_data_write("Etc");
  cg_goto(fn, bn, "Zone_t", zn, "Etc", 0, "end");
  cgsize_t *N = malloc(sizeof(cgsize_t));
  N[0] = 1;
  cg_array_write("Time", RealDouble, 1, N, &ttime);
  free(N);

  cg_close(fn);
  free(sc_out);
  free(conv_out);
  free(diff_out);
  free(scSrc_out);
  free(epsp_out);

  free(omegaX_out);
  free(omegaY_out);
  free(omegaZ_out);

}




void cgns_point_particles(real dtout)
{
  if(npoints > 0) {
    // create the solution file
    char fname[FILE_NAME_SIZE]="";
    char fname2[FILE_NAME_SIZE]="";
    char fnameall[FILE_NAME_SIZE]="";
    char fnameall2[FILE_NAME_SIZE]="";
    real tout = ttime; // = rec_point_particle_stepnum_out * dtout;
    char format[CHAR_BUF_SIZE]="";
    int sigfigs = ceil(log10(1. / dtout));
    if(sigfigs < 1) sigfigs = 1;
    sprintf(format, "%%.%df", sigfigs);
    sprintf(fname2, "point-%s.cgns", format);
    sprintf(fnameall2, "%s/output/point-%s.cgns", ROOT_DIR, format);
    sprintf(fname, fname2, tout);
    sprintf(fnameall, fnameall2, tout);
    int fn;
    int bn;
    int zn;
    int en;
    int sn;
    int Xn;
    int Yn;
    int Zn;
    int fnr;
    cg_open(fnameall, CG_MODE_WRITE, &fn);
    cg_base_write(fn, "Base", 3, 3, &bn);
    cgsize_t size[3][1];
    size[0][0] = npoints;
    size[1][0] = 0;
    size[2][0] = 0;
    cg_zone_write(fn, bn, "Zone0", size[0], Unstructured, &zn);

    // write point_particle locations
    real *x = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *y = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *z = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    cgsize_t *conn = malloc(npoints * sizeof(cgsize_t));
    // cpumem += npoints * sizeof(int);
    real *a = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(int);
 //   int *order = malloc(npoints * sizeof(int));
    // cpumem += npoints * sizeof(real);
    real *u = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *v = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *w = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *udot = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *vdot = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *wdot = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *ox = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *oy = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *oz = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *Fx = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *Fy = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *Fz = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *Lx = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *Ly = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *Lz = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real); 
    real *iFx = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *iFy = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *iFz = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *iLx = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *iLy = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *iLz = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real); 
/* 
    real *kFx = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *kFy = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *kFz = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
  */
    real *hFx = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *hFy = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *hFz = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *hLx = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *hLy = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *hLz = malloc(npoints * sizeof(real));


    real *msdot = malloc(npoints * sizeof(real));
    real *Nu = malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *ms = malloc(npoints * sizeof(real));
    real *den = malloc(npoints * sizeof(real));
    real *pdt = malloc(npoints * sizeof(real));

   int *ii = malloc(npoints * sizeof(int));
   int *jj = malloc(npoints * sizeof(int));
   int *kk = malloc(npoints * sizeof(int));
  
   long int *id = malloc(npoints * sizeof(long int));




    // cpumem += npoints * sizeof(real); 
    for(int i = 0; i < npoints; i++) {
      real mass = 4./3.*PI*(points[i].rho-rho_f)*points[i].r*points[i].r*points[i].r;
      x[i] = points[i].x;
      y[i] = points[i].y;
      z[i] = points[i].z;
      conn[i] = npoints-i;
      a[i] = points[i].r;
//      order[i] = points[i].order;
      u[i] = points[i].u;
      v[i] = points[i].v;
      w[i] = points[i].w;
      udot[i] = points[i].udot;
      vdot[i] = points[i].vdot;
      wdot[i] = points[i].wdot;
  /*
      kFx[i] = points[i].kFx;
      kFy[i] = points[i].kFy;
      kFz[i] = points[i].kFz;
  */
      iFx[i] = points[i].iFx;
      iFy[i] = points[i].iFy;
      iFz[i] = points[i].iFz;
      iLx[i] = points[i].iLx;
      iLy[i] = points[i].iLy;
      iLz[i] = points[i].iLz;
      hFx[i] = points[i].Fx;
      hFy[i] = points[i].Fy;
      hFz[i] = points[i].Fz;
      hLx[i] = points[i].Lx;
      hLy[i] = points[i].Ly;
      hLz[i] = points[i].Lz;
      Fx[i] = iFx[i] + hFx[i] + mass*g.x;
      Fy[i] = iFy[i] + hFy[i] + mass*g.y;
      Fz[i] = iFz[i] + hFz[i] + mass*g.z;
      Lx[i] = iLx[i] + hLx[i];
      Ly[i] = iLy[i] + hLy[i];
      Lz[i] = iLz[i] + hLz[i];
      ox[i] = points[i].ox;
      oy[i] = points[i].oy;
      oz[i] = points[i].oz;

	msdot[i]=points[i].msdot;
	Nu[i]=points[i].Nu;
	ms[i]=points[i].ms;
	den[i]=points[i].rho;

	ii[i]=points[i].i;
	jj[i]=points[i].j;
	kk[i]=points[i].k;

	id[i]=points[i].id;
	pdt[i]=points[i].dt;


}
    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateX", x, &Xn);
    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateY", y, &Yn);
    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateZ", z, &Zn);

    cg_section_write(fn, bn, zn, "Elements", NODE, 0, npoints-1, 0, conn, &en);

    cg_sol_write(fn, bn, zn, "Solution", Vertex, &sn);
    cg_field_write(fn, bn, zn, sn, RealDouble, "Radius", a, &fnr);
//    cg_field_write(fn, bn, zn, sn, Integer, "LambOrder", order, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityX", u, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityY", v, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "VelocityZ", w, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AccelerationX", udot, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AccelerationY", vdot, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AccelerationZ", wdot, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularVelocityX", ox, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularVelocityY", oy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "AngularVelocityZ", oz, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "HydroForceX", hFx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "HydroForceY", hFy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "HydroForceZ", hFz, &fnr);
/*
    cg_field_write(fn, bn, zn, sn, RealDouble, "SpringForceX", kFx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "SpringForceY", kFy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "SpringForceZ", kFz, &fnr);
*/
    cg_field_write(fn, bn, zn, sn, RealDouble, "InteractionForceX", iFx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "InteractionForceY", iFy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "InteractionForceZ", iFz, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "TotalForceX", Fx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "TotalForceY", Fy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "TotalForceZ", Fz, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "MomentX", Lx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "MomentY", Ly, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "MomentZ", Lz, &fnr);
    
//add by shigan 10_30_2014
    cg_field_write(fn, bn, zn, sn, RealDouble, "SolvableMassAcc", msdot, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "SolvableMassNu", Nu, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "SolvableMass", ms, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "Density", den, &fnr);
 
    cg_field_write(fn, bn, zn, sn, Integer, "gridIDi", ii, &fnr);
    cg_field_write(fn, bn, zn, sn, Integer, "gridIDj", jj, &fnr);
    cg_field_write(fn, bn, zn, sn, Integer, "gridIDk", kk, &fnr);
//    cg_field_write(fn, bn, zn, sn, Integer, "ID", id, &fnr);
    cg_field_write(fn, bn, zn, sn, Integer, "ID", id, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "dt", pdt, &fnr);

    cg_goto(fn, bn, "Zone_t", zn, "end");
    cg_user_data_write("Etc");
    cg_goto(fn, bn, "Zone_t", zn, "Etc", 0, "end");
    cgsize_t *N = malloc(sizeof(cgsize_t));
    N[0] = 1;
    cg_array_write("Time", RealDouble, 1, N, &ttime);
    free(N);

    cg_close(fn);
    free(x);
    free(y);
    free(z);
    free(conn);
    free(a);
  //  free(order);
    free(u);
    free(v);
    free(w);
    free(udot);
    free(vdot);
    free(wdot);
/*
    free(kFx);
    free(kFy);
    free(kFz);
*/
    free(iFx);
    free(iFy);
    free(iFz);
    free(iLx);
    free(iLy);
    free(iLz);
    free(hFx);
    free(hFy);
    free(hFz);
    free(hLx);
    free(hLy);
    free(hLz);
    free(ox);
    free(oy);
    free(oz);
    free(Fx);
    free(Fy);
    free(Fz);
    free(Lx);
    free(Ly);
    free(Lz);
    
    free(ms);
    free(msdot);
    free(Nu);
    free(den);
    free(pdt);
    free(ii);
    free(jj);
    free(kk);
    free(id);

  }
}

void recorder_bicgstab_init(char *name)
{
  // create the file
  char path[FILE_NAME_SIZE];
  sprintf(path, "%s/record/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "w");
  if(rec == NULL) {
    fprintf(stderr, "Could not open file %s\n", name);
    exit(EXIT_FAILURE);
  }

  fprintf(rec, "%-12s", "stepnum");
  fprintf(rec, "%-15s", "ttime");
  fprintf(rec, "%-15s", "dt");
  fprintf(rec, "%-8s", "niter");
  fprintf(rec, "%-15s", "resid");

  // close the file
  fclose(rec);
}

void recorder_lamb_init(char *name)
{
  // create the file
  char path[FILE_NAME_SIZE];
  sprintf(path, "%s/record/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "w");
  if(rec == NULL) {
    fprintf(stderr, "Could not open file %s\n", name);
    exit(EXIT_FAILURE);
  }

  // close the file
  fclose(rec);
}

void recorder_bicgstab(char *name, int niter, real resid)
{
  // open the file
  char path[FILE_NAME_SIZE];
  sprintf(path, "%s/record/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "r+");

  // move to the end of the file
  fseek(rec, 0, SEEK_END);

  fprintf(rec, "\n");
  fprintf(rec, "%-12d", stepnum);
  fprintf(rec, "%-15e", ttime);
  fprintf(rec, "%-15e", dt);
  fprintf(rec, "%-8d", niter);
  fprintf(rec, "%-15e", resid);

  // close the file
  fclose(rec);
}
/*
void recorder_lamb(char *name, int iter)
{
  // open the file
  char path[FILE_NAME_SIZE];
  sprintf(path, "%s/record/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "r+");

  int n = 0;
  int m = 0;
  int c = 0;

  // move to the end of the file
  fseek(rec, 0, SEEK_END);

  fprintf(rec, "ttime = %e, iter = %d\n", ttime, iter);
  for(int i = 0; i < npoints; i++) {
    fprintf(rec, "point_particle %d\n", i);
    fprintf(rec, "%4s%4s%12s%12s", "n", "m", "pnm_re", "pnm_im");
    fprintf(rec, "%12s%12s", "phinm_re", "phinm_im");
    fprintf(rec, "%12s%12s\n", "chinm_re", "chinm_im");
    n = 0;
    m = 0;
    c = 0;
    while(c < points[i].ncoeff) {
      fprintf(rec, "%4d%4d", n, m);
      fprintf(rec, "%12.3e%12.3e",
        pnm_re[coeff_stride*i+c], pnm_im[coeff_stride*i+c]);
      fprintf(rec, "%12.3e%12.3e",
        phinm_re[coeff_stride*i+c], phinm_im[coeff_stride*i+c]);
      fprintf(rec, "%12.3e%12.3e\n",
        chinm_re[coeff_stride*i+c], chinm_im[coeff_stride*i+c]);
      m++;
      c++;
      if(m > n) {
        n++;
        m = 0;
      }
    }
  }

  // close the file
  fclose(rec);
}
*/
void recorder_cfl_init(char *name)
{
  // create the file
  char path[FILE_NAME_SIZE];
  sprintf(path, "%s/record/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "w");
  if(rec == NULL) {
    fprintf(stderr, "Could not open file %s\n", name);
    exit(EXIT_FAILURE);
  }

  fprintf(rec, "%-12s", "stepnum");
  fprintf(rec, "%-15s", "ttime");
  fprintf(rec, "%-15s", "dt");
  fprintf(rec, "%-15s", "cfl");

  // close the file
  fclose(rec);
}

void recorder_cfl(char *name, real cfl)
{
  // open the file
  char path[FILE_NAME_SIZE];
  sprintf(path, "%s/record/%s", ROOT_DIR, name);
  FILE *rec = fopen(path, "r+");

  // move to the end of the file
  fseek(rec, 0, SEEK_END);

  fprintf(rec, "\n");
  fprintf(rec, "%-12d", stepnum);
  fprintf(rec, "%-15e", ttime);
  fprintf(rec, "%-15e", dt);
  fprintf(rec, "%-15e", cfl);

  // close the file
  fclose(rec);
}
