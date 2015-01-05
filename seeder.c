#include <time.h>

#include "bluebottle.h"
#include "scalar.h"
#include "domain.h"
#include "point.h"

void seeder(int N, real a, real rho, int o, int t, int r) {
  printf("Running bluebottle seeder for %d point_particles...\n\n", N);
  fflush(stdout);
  real xx, yy, zz;
  int fits = 1;
  int attempts = 1;
  int fail = 0;
  int redo = 1;

  // seed the random number generator
  srand(time(NULL));

  // read domain input
 // domain_read_input();
 // domain_init();

  turb_read_input();
  domain_init_turb();

  npoints = N;

  // allocate pointicle list
  points = (point_struct*) malloc(npoints * sizeof(point_struct));
  cpumem += npoints * sizeof(point_struct);

  // place the first pointicle
  points[0].r = a;
  redo = 1;
  while(redo == 1) {
    redo = 0;
    points[0].x = rand() / (real)RAND_MAX;
    points[0].x *= Dom.xl;
    points[0].x += Dom.xs;
    if((bc.uW != PERIODIC) && (points[0].x < (Dom.xs + 1.05*points[0].r)))
      redo = 1;
    if((bc.uE != PERIODIC) && (points[0].x > (Dom.xe - 1.05*points[0].r)))
      redo = 1;
  }
  redo = 1;
  while(redo == 1) {
    redo = 0;
    points[0].y = rand() / (real)RAND_MAX;
    points[0].y *= Dom.yl;
    points[0].y += Dom.ys;
    if((bc.vS != PERIODIC) && (points[0].y < (Dom.ys + 1.05*points[0].r)))
      redo = 1;
    if((bc.vN != PERIODIC) && (points[0].y > (Dom.ye - 1.05*points[0].r)))
      redo = 1;
  }
  redo = 1;
  while(redo == 1) {
    redo = 0;
    points[0].z = rand() / (real)RAND_MAX;
    //points[0].z = acos(2.*(points[0].z-0.5))/PI;
    points[0].z *= Dom.zl;
    points[0].z += Dom.zs;
    if((bc.wB != PERIODIC) && (points[0].z < (Dom.zs + 1.05*points[0].r)))
      redo = 1;
    if((bc.wT != PERIODIC) && (points[0].z > (Dom.ze - 1.05*points[0].r)))
      redo = 1;
  }

  points[0].u = 0;
  points[0].v = 0;
  points[0].w = 0;
 
  points[0].rho = rho;
/* 
  points[0].aLx = 0;
  points[0].aLy = 0;
  points[0].aLz = 0;
  points[0].order = o;
  points[0].rs = 1.05;
  points[0].ncoeff = 0;
  points[0].translating = t;
  */
  points[0].rotating = r;

  // place the rest of the point_particles
  int i = 0;
  for(i = 1; i < npoints; i++) {
    fits = !fits;
    if(fail) break;
    while(!fits) {
      attempts++;
      // place pointicle
      points[i].r = a;
      redo = 1;
      while(redo == 1) {
        redo = 0;
        points[i].x = rand() / (real)RAND_MAX;
        points[i].x *= Dom.xl;
        points[i].x += Dom.xs;
        if((bc.uW != PERIODIC) && (points[i].x < (Dom.xs + 1.05*points[i].r)))
          redo = 1;
        if((bc.uE != PERIODIC) && (points[i].x > (Dom.xe - 1.05*points[i].r)))
          redo = 1;
      }
      redo = 1;
      while(redo == 1) {
        redo = 0;
        points[i].y = rand() / (real)RAND_MAX;
        points[i].y *= Dom.yl;
        points[i].y += Dom.ys;
        if((bc.vS != PERIODIC) && (points[i].y < (Dom.ys + 1.05*points[i].r)))
          redo = 1;
        if((bc.vN != PERIODIC) && (points[i].y > (Dom.ye - 1.05*points[i].r)))
          redo = 1;
      }
      redo = 1;
      while(redo == 1) {
        redo = 0;
        points[i].z = rand() / (real)RAND_MAX;
        points[i].z *= Dom.zl;
        points[i].z += Dom.zs;
        if((bc.wB != PERIODIC) && (points[i].z < (Dom.zs + 1.05*points[i].r)))
          redo = 1;
        if((bc.wT != PERIODIC) && (points[i].z > (Dom.ze - 1.05*points[i].r)))
          redo = 1;
      }

      points[i].u = 0;
      points[i].v = 0;
      points[i].w = 0;
  /*    points[i].aFx = 0;
      points[i].aFy = 0;
      points[i].aFz = 0;
      points[i].aLx = 0;
      points[i].aLy = 0;
      points[i].aLz = 0;
      points[i].order = o;
      points[i].rs = 1.05;
      points[i].ncoeff = 0;
      points[i].translating = t;
    */
      points[i].rho = rho;
      points[i].rotating = r;

      // check that this pointicle does not intersect any other pointicle
      fits = !fits;
      for(int j = 0; j < i; j++) {
        xx = points[i].x - points[j].x;
        xx = xx * xx;
        yy = points[i].y - points[j].y;
        yy = yy * yy;
        zz = points[i].z - points[j].z;
        zz = zz * zz;
        if(sqrt(xx + yy + zz) < (1.05*points[i].r + 1.05*points[j].r)) {
          fits = !fits;
          break;
        }

        // also use virtual pointicle to check if a pointicle is too close in
        // a periodic direction
        // x only
        if(bc.uW == PERIODIC) {
          xx = points[i].x - points[j].x;
          xx = xx * xx;
          yy = points[i].y - points[j].y;
          yy = yy * yy;
          zz = points[i].z - points[j].z;
          zz = zz * zz;
          if(points[i].x < (Dom.xs + points[i].r))
            xx = points[i].x + Dom.xl - points[j].x;
          if(points[i].x > (Dom.xe - points[i].r))
            xx = points[i].x - Dom.xl - points[j].x;
          xx = xx * xx;
          yy = points[i].y - points[j].y;
          yy = yy * yy;
          zz = points[i].z - points[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[i].r + 1.05*points[j].r)) {
            fits = !fits;
            break;
          }
        }

        // y only
        if(bc.vS == PERIODIC) {
          xx = points[i].x - points[j].x;
          xx = xx * xx;
          yy = points[i].y - points[j].y;
          yy = yy * yy;
          zz = points[i].z - points[j].z;
          zz = zz * zz;
          xx = points[i].x - points[j].x;
          xx = xx * xx;
          if(points[i].y < (Dom.ys + points[i].r))
            yy = points[i].y + Dom.yl - points[j].y;
          if(points[i].y > (Dom.ye - points[i].r))
            yy = points[i].y - Dom.yl - points[j].y;
          yy = yy * yy;
          zz = points[i].z - points[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[i].r + 1.05*points[j].r)) {
            fits = !fits;
            break;
          }
        }

        // z only
        if(bc.wB == PERIODIC) {
          xx = points[i].x - points[j].x;
          xx = xx * xx;
          yy = points[i].y - points[j].y;
          yy = yy * yy;
          zz = points[i].z - points[j].z;
          zz = zz * zz;
          xx = points[i].x - points[j].x;
          xx = xx * xx;
          yy = points[i].y - points[j].y;
          yy = yy * yy;
          if(points[i].z < (Dom.zs + points[i].r))
            zz = points[i].z + Dom.zl - points[j].z;
          if(points[i].z > (Dom.ze - points[i].r))
            zz = points[i].z - Dom.zl - points[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[i].r + 1.05*points[j].r)) {
            fits = !fits;
            break;
          }
        }

        // x and y
        if(bc.uW == PERIODIC && bc.vS == PERIODIC) {
          xx = points[i].x - points[j].x;
          xx = xx * xx;
          yy = points[i].y - points[j].y;
          yy = yy * yy;
          zz = points[i].z - points[j].z;
          zz = zz * zz;
          if(points[i].x < (Dom.xs + points[i].r))
            xx = points[i].x + Dom.xl - points[j].x;
          if(points[i].x > (Dom.xe - points[i].r))
            xx = points[i].x - Dom.xl - points[j].x;
          xx = xx * xx;
          if(points[i].y < (Dom.ys + points[i].r))
            yy = points[i].y + Dom.yl - points[j].y;
          if(points[i].y > (Dom.ye - points[i].r))
            yy = points[i].y - Dom.yl - points[j].y;
          yy = yy * yy;
          zz = points[i].z - points[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[i].r + 1.05*points[j].r)) {
            fits = !fits;
            break;
          }
        }

        // y and z
        if(bc.vS == PERIODIC && bc.wB == PERIODIC) {
          xx = points[i].x - points[j].x;
          xx = xx * xx;
          yy = points[i].y - points[j].y;
          yy = yy * yy;
          zz = points[i].z - points[j].z;
          zz = zz * zz;
          xx = points[i].x - points[j].x;
          xx = xx * xx;
          if(points[i].y < (Dom.ys + points[i].r))
            yy = points[i].y + Dom.yl - points[j].y;
          if(points[i].y > (Dom.ye - points[i].r))
            yy = points[i].y - Dom.yl - points[j].y;
          yy = yy * yy;
          if(points[i].z < (Dom.zs + points[i].r))
            zz = points[i].z + Dom.zl - points[j].z;
          if(points[i].z > (Dom.ze - points[i].r))
            zz = points[i].z - Dom.zl - points[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[i].r + 1.05*points[j].r)) {
            fits = !fits;
            break;
          }
        }

        // z and x
        if(bc.uW == PERIODIC && bc.wB == PERIODIC) {
          xx = points[i].x - points[j].x;
          xx = xx * xx;
          yy = points[i].y - points[j].y;
          yy = yy * yy;
          zz = points[i].z - points[j].z;
          zz = zz * zz;
          if(points[i].x < (Dom.xs + points[i].r))
            xx = points[i].x + Dom.xl - points[j].x;
          if(points[i].x > (Dom.xe - points[i].r))
            xx = points[i].x - Dom.xl - points[j].x;
          xx = xx * xx;
          yy = points[i].y - points[j].y;
          yy = yy * yy;
          if(points[i].z < (Dom.zs + points[i].r))
            zz = points[i].z + Dom.zl - points[j].z;
          if(points[i].z > (Dom.ze - points[i].r))
            zz = points[i].z - Dom.zl - points[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[i].r + 1.05*points[j].r)) {
            fits = !fits;
            break;
          }
        }

        // x, y, and z
        if(bc.uW == PERIODIC && bc.vS == PERIODIC && bc.wB == PERIODIC) {
          xx = points[i].x - points[j].x;
          xx = xx * xx;
          yy = points[i].y - points[j].y;
          yy = yy * yy;
          zz = points[i].z - points[j].z;
          zz = zz * zz;
          if(points[i].x < (Dom.xs + points[i].r))
            xx = points[i].x + Dom.xl - points[j].x;
          if(points[i].x > (Dom.xe - points[i].r))
            xx = points[i].x - Dom.xl - points[j].x;
          xx = xx * xx;
          if(points[i].y < (Dom.ys + points[i].r))
            yy = points[i].y + Dom.yl - points[j].y;
          if(points[i].y > (Dom.ye - points[i].r))
            yy = points[i].y - Dom.yl - points[j].y;
          yy = yy * yy;
          if(points[i].z < (Dom.zs + points[i].r))
            zz = points[i].z + Dom.zl - points[j].z;
          if(points[i].z > (Dom.ze - points[i].r))
            zz = points[i].z - Dom.zl - points[j].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[i].r + 1.05*points[j].r)) {
            fits = !fits;
            break;
          }
        }

        // check both ways
        // x only
        if(bc.uW == PERIODIC) {
          xx = points[j].x - points[i].x;
          xx = xx * xx;
          yy = points[j].y - points[i].y;
          yy = yy * yy;
          zz = points[j].z - points[i].z;
          zz = zz * zz;
          if(points[j].x < (Dom.xs + points[j].r))
            xx = points[j].x + Dom.xl - points[i].x;
          if(points[j].x > (Dom.xe - points[j].r))
            xx = points[j].x - Dom.xl - points[i].x;
          xx = xx * xx;
          yy = points[j].y - points[i].y;
          yy = yy * yy;
          zz = points[j].z - points[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[j].r + 1.05*points[i].r)) {
            fits = !fits;
            break;
          }
        }

        // y only
        if(bc.vS == PERIODIC) {
          xx = points[j].x - points[i].x;
          xx = xx * xx;
          yy = points[j].y - points[i].y;
          yy = yy * yy;
          zz = points[j].z - points[i].z;
          zz = zz * zz;
          xx = points[j].x - points[i].x;
          xx = xx * xx;
          if(points[j].y < (Dom.ys + points[j].r))
            yy = points[j].y + Dom.yl - points[i].y;
          if(points[j].y > (Dom.ye - points[j].r))
            yy = points[j].y - Dom.yl - points[i].y;
          yy = yy * yy;
          zz = points[j].z - points[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[j].r + 1.05*points[i].r)) {
            fits = !fits;
            break;
          }
        }

        // z only
        if(bc.wB == PERIODIC) {
          xx = points[j].x - points[i].x;
          xx = xx * xx;
          yy = points[j].y - points[i].y;
          yy = yy * yy;
          zz = points[j].z - points[i].z;
          zz = zz * zz;
          xx = points[j].x - points[i].x;
          xx = xx * xx;
          yy = points[j].y - points[i].y;
          yy = yy * yy;
          if(points[j].z < (Dom.zs + points[j].r))
            zz = points[j].z + Dom.zl - points[i].z;
          if(points[j].z > (Dom.ze - points[j].r))
            zz = points[j].z - Dom.zl - points[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[j].r + 1.05*points[i].r)) {
            fits = !fits;
            break;
          }
        }

        // x and y
        if(bc.uW == PERIODIC && bc.vS == PERIODIC) {
          xx = points[j].x - points[i].x;
          xx = xx * xx;
          yy = points[j].y - points[i].y;
          yy = yy * yy;
          zz = points[j].z - points[i].z;
          zz = zz * zz;
          if(points[j].x < (Dom.xs + points[j].r))
            xx = points[j].x + Dom.xl - points[i].x;
          if(points[j].x > (Dom.xe - points[j].r))
            xx = points[j].x - Dom.xl - points[i].x;
          xx = xx * xx;
          if(points[j].y < (Dom.ys + points[j].r))
            yy = points[j].y + Dom.yl - points[i].y;
          if(points[j].y > (Dom.ye - points[j].r))
            yy = points[j].y - Dom.yl - points[i].y;
          yy = yy * yy;
          zz = points[j].z - points[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[j].r + 1.05*points[i].r)) {
            fits = !fits;
            break;
          }
        }

        // y and z
        if(bc.vS == PERIODIC && bc.wB == PERIODIC) {
          xx = points[j].x - points[i].x;
          xx = xx * xx;
          yy = points[j].y - points[i].y;
          yy = yy * yy;
          zz = points[j].z - points[i].z;
          zz = zz * zz;
          xx = points[j].x - points[i].x;
          xx = xx * xx;
          if(points[j].y < (Dom.ys + points[j].r))
            yy = points[j].y + Dom.yl - points[i].y;
          if(points[j].y > (Dom.ye - points[j].r))
            yy = points[j].y - Dom.yl - points[i].y;
          yy = yy * yy;
          if(points[j].z < (Dom.zs + points[j].r))
            zz = points[j].z + Dom.zl - points[i].z;
          if(points[j].z > (Dom.ze - points[j].r))
            zz = points[j].z - Dom.zl - points[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[j].r + 1.05*points[i].r)) {
            fits = !fits;
            break;
          }
        }

        // z and x
        if(bc.uW == PERIODIC && bc.wB == PERIODIC) {
          xx = points[j].x - points[i].x;
          xx = xx * xx;
          yy = points[j].y - points[i].y;
          yy = yy * yy;
          zz = points[j].z - points[i].z;
          zz = zz * zz;
          if(points[j].x < (Dom.xs + points[j].r))
            xx = points[j].x + Dom.xl - points[i].x;
          if(points[j].x > (Dom.xe - points[j].r))
            xx = points[j].x - Dom.xl - points[i].x;
          xx = xx * xx;
          yy = points[j].y - points[i].y;
          yy = yy * yy;
          if(points[j].z < (Dom.zs + points[j].r))
            zz = points[j].z + Dom.zl - points[i].z;
          if(points[j].z > (Dom.ze - points[j].r))
            zz = points[j].z - Dom.zl - points[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[j].r + 1.05*points[i].r)) {
            fits = !fits;
            break;
          }
        }

        // x, y, and z
        if(bc.uW == PERIODIC && bc.vS == PERIODIC && bc.wB == PERIODIC) {
          xx = points[j].x - points[i].x;
          xx = xx * xx;
          yy = points[j].y - points[i].y;
          yy = yy * yy;
          zz = points[j].z - points[i].z;
          zz = zz * zz;
          if(points[j].x < (Dom.xs + points[j].r))
            xx = points[j].x + Dom.xl - points[i].x;
          if(points[j].x > (Dom.xe - points[j].r))
            xx = points[j].x - Dom.xl - points[i].x;
          xx = xx * xx;
          if(points[j].y < (Dom.ys + points[j].r))
            yy = points[j].y + Dom.yl - points[i].y;
          if(points[j].y > (Dom.ye - points[j].r))
            yy = points[j].y - Dom.yl - points[i].y;
          yy = yy * yy;
          if(points[j].z < (Dom.zs + points[j].r))
            zz = points[j].z + Dom.zl - points[i].z;
          if(points[j].z > (Dom.ze - points[j].r))
            zz = points[j].z - Dom.zl - points[i].z;
          zz = zz * zz;
          if(sqrt(xx + yy + zz) < (1.05*points[j].r + 1.05*points[i].r)) {
            fits = !fits;
            break;
          }
        }

      }
      if(attempts == 1e5*npoints) {
        fail = !fail;
        break;
      }
    }
  }

  if(fail) {
    printf("After %d attempts, the seeder has placed", attempts);
    printf(" %d of %d point_particles (a = %f).\n\n", i-1, npoints, a);
    printf("...bluebottle seeder done.\n\n");
    exit(EXIT_FAILURE);
  }

  printf("It took %d attempts to place %d", attempts, npoints);
  printf(" point_particles (a = %f) with no intersections.\n\n", a);
  fflush(stdout);

  printf("Writing point_seeder.input...");
  fflush(stdout);
  // write pointicle configuration to file
  char fname[FILE_NAME_SIZE];
  // open file for writing
  sprintf(fname, "%spoint_seeder.input", INPUT_DIR);
  FILE *ofile = fopen(fname, "w");
  if(ofile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write the number of point_particles
  fprintf(ofile, "n %d\n", npoints);

  // write each pointicle configuration
  for(int i = 0; i < npoints; i++) {
    fprintf(ofile, "\n");
    fprintf(ofile, "r %f\n", points[i].r);
    fprintf(ofile, "(x, y, z) %f %f %f\n", points[i].x, points[i].y, points[i].z);
 /*
   fprintf(ofile, "(aFx, aFy, aFz) %f %f %f\n", points[i].aFx, points[i].aFy,
      points[i].aFz);
    fprintf(ofile, "(aLx, aLy, aLz) %f %f %f\n", points[i].aLx, points[i].aLy,
      points[i].aLz);
    fprintf(ofile, "order %d\n", points[i].order);
    fprintf(ofile, "spring_k %f\n", 0.);//points[i].spring_k);
    fprintf(ofile, "spring (x, y, z) %f %f %f\n", 0., 0., 0.);
    fprintf(ofile, "translating %d\n", points[i].translating);
*/ 
    fprintf(ofile, "rho %f\n", points[i].rho);
    fprintf(ofile, "rotating %d\n", points[i].rotating);
  }

  // close the file
  fclose(ofile);
  printf("done.\n");
  printf("\n...bluebottle seeder done.\n\n");
  fflush(stdout);

  // clean up
  domain_clean();
  points_clean();
}
