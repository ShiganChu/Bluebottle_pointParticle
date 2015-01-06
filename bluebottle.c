#include <mpi.h>
#include "bluebottle.h"
//#include "shigan.h"
#include "point.h"
#include "precursor.h"
#include "recorder.h"

//#include "cuPrintf.cu"
// define global variables that were declared in header file
int dev_start;
int dev_end;
dom_struct *dom;
dom_struct **_dom;
dom_struct Dom;

real **_stress_u;
real **_stress_v;
real **_stress_w;

real C_add;
real C_stress;
real C_drag;

real **_omega_x;
real **_omega_y;
real **_omega_z;

real *p0;
real *p;
real *divU;
real **_p0;
real **_p;
real **_divU;
real *u;
real *u0;
real **_u;
real **_u0;
real *v;
real *v0;
real **_v;
real **_v0;
real *w;
real *w0;
real **_w;
real **_w0;
real *f_x;
real *f_y;
real *f_z;
real **_f_x;
real **_f_y;
real **_f_z;
real *diff0_u;
real *diff0_v;
real *diff0_w;
real **_diff0_u;
real **_diff0_v;
real **_diff0_w;
real *conv0_u;
real *conv0_v;
real *conv0_w;
real **_conv0_u;
real **_conv0_v;
real **_conv0_w;
real *diff_u;
real *diff_v;
real *diff_w;
real **_diff_u;
real **_diff_v;
real **_diff_w;
real *conv_u;
real *conv_v;
real *conv_w;
real **_conv_u;
real **_conv_v;
real **_conv_w;
real **_u_star;
real **_v_star;
real **_w_star;
real *u_star;
real *v_star;
real *w_star;
real *u_WE;
real *u_SN;
real *u_BT;
real *v_WE;
real *v_SN;
real *v_BT;
real *w_WE;
real *w_SN;
real *w_BT;
real **_u_WE;
real **_u_SN;
real **_u_BT;
real **_v_WE;
real **_v_SN;
real **_v_BT;
real **_w_WE;
real **_w_SN;
real **_w_BT;
real **_rhs_p;
real duration;
real ttime;
real dt;
real dt0;
real CFL;
int pp_max_iter;
real pp_residual;
int lamb_max_iter;
real lamb_residual;
real lamb_relax;
int out_plane;
int stepnum;
int rec_flow_field_stepnum_out;
int rec_paraview_stepnum_out;
int rec_point_particle_stepnum_out;
int rec_precursor_stepnum_out;
real rec_flow_field_dt;
int rec_flow_field_vel;
int rec_flow_field_p;
int rec_flow_field_phase;
real rec_paraview_dt;
real rec_point_particle_dt;
real rec_restart_dt;
real rec_precursor_dt;
int rec_point_particle_pos;
int rec_point_particle_a;
int rec_point_particle_vel;
int rec_point_particle_omega;
int rec_point_particle_force;
int rec_point_particle_moment;
g_struct g;
real rho_f;
real mu;
real nu;
BC bc;
int init_cond;
gradP_struct gradP;
real turbA;
real turbl;
long int cpumem;
long int gpumem;
int bc_flow_configs[18];
real bc_flow_vels[18];
real bc_plane_pos[9];

int main(int argc, char *argv[]) {

  int np = 0;
  int rank = 0;

  // set up MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank == MASTER) {
    int turb = 0;
    MPI_Status status; // for MPI communication

    // parse command-line arguments
    // if none given, run normally
    // if -s given, run seeder program only
    // if anything else given, exit
    // this code comes from K&R, p. 117

    int argin;
    int runseeder = 0;
    int runrestart = 0;
    int NP = 0;
    real radius = -1.;
    real density = -1.;
    int order = -1.;
    int translating = -1;
    int rotating = -1;
    while(--argc > 0 && (*++argv)[0] == '-') {
      while((argin = *++argv[0])) {
        switch(argin) {
          case 's':
            runseeder = 1;
            break;
          case 'r':
            runrestart = 1;
            break;
          default:
            runseeder = 2;
            runrestart = 2;
            //printf("bluebottle: illegal option %c\n", argin);
            argc = 0;
            break;
        }
      }
    }
    
    NP = 0;
    radius = -1.;
    order = -1.;
    if(runseeder == 1) {
      if(argc != 6) {
        //printf("Usage: bluebottle -s N a d o t r\n");
        //printf("       N is the number of point_particles\n");
        //printf("       a is the point_particle radius\n");
        //printf("       d is the point_particle density (rho)\n");
        //printf("       o is the Lamb's solution truncation order\n");
        //printf("       t = 1 if point_particles can translate (0 if not)\n");
        //printf("       r = 1 if point_particles can rotate (0 if not)\n");
        return EXIT_FAILURE;
      }
      NP = atoi(argv[0]);
      radius = atof(argv[1]);
      density = atof(argv[2]);
      order = atoi(argv[3]);
      translating = atoi(argv[4]);
      rotating = atoi(argv[5]);
      if(NP < 1 || radius < 0 || density < 0 || order < 0
        || translating > 1 || rotating > 1) {
        //printf("Usage: bluebottle -s N a d o t r\n");
        //printf("       N is the number of point_particles\n");
        //printf("       a is the point_particle radius\n");
        //printf("       d is the point_particle density (rho)\n");
        //printf("       o is the Lamb's solution truncation order\n");
        //printf("       t = 1 if point_particles can translate (0 if not)\n");
        //printf("       r = 1 if point_particles can rotate (0 if not)\n");
        return EXIT_FAILURE;
      }
      seeder(NP, radius, density, order, translating, rotating);
      return EXIT_SUCCESS;
    } else if(runrestart == 1 && argc > 0) {
      //printf("Usage restart simulation: bluebottle -r\n");
      return EXIT_FAILURE;
    } else if(runseeder == 2) {
      return EXIT_FAILURE;
    } else if(runrestart == 2) {
      return EXIT_FAILURE;
    } else {
      // read recorder config file
      recorder_read_config();

      // read simulation input configuration file
      //printf("\nRunning bluebottle_0.1...\n\n");
      //printf("Reading the domain and point_particle input files...\n\n");
      domain_read_input();
      points_read_input(turb);
      fflush(stdout);
      //printf("FLOW: Using devices %d through %d.\n\n", dev_start, dev_end);
      fflush(stdout);

      if(runrestart != 1) {
        // start BICGSTAB recorder
        recorder_bicgstab_init("solver_flow.rec");
        // start Lamb's coefficient recorder
        //recorder_lamb_init("lamb.rec");
      }

      // initialize the domain
      //printf("Initializing domain variables...");
      fflush(stdout);
      int domain_init_flag = domain_init();
      //printf("done.\n");
      fflush(stdout);
      if(domain_init_flag == EXIT_FAILURE) {
        //printf("\nThe number of devices in DEV RANGE is insufficient\n");
        //printf("for the given domain decomposition.  Exiting now.\n");
        return EXIT_FAILURE;
      }

      // set up the boundary condition config info to send to precursor
      expd_init_BC(np);

      // initialize the point_particles
      //printf("Initializing point_particle variables...");
      fflush(stdout);
      int points_init_flag = points_init();
      //printf("done.\n");
      fflush(stdout);
      if(points_init_flag == EXIT_FAILURE) {
        //printf("\nThe initial point_particle configuration is not allowed.\n");
        return EXIT_FAILURE;
      }

      // allocate device memory
      //printf("Allocating domain CUDA device memory...");
      fflush(stdout);
      cuda_dom_malloc();
      //printf("...done.\n");
      fflush(stdout);
      //printf("Allocating point_particle CUDA device memory...");
      fflush(stdout);
      cuda_point_malloc();
      //printf("...done.\n");
      fflush(stdout);

      // copy host data to devices
      //printf("Copying host domain data to devices...");
      fflush(stdout);
      cuda_dom_push();
      //printf("done.\n");
      fflush(stdout);
      //printf("Copying host point_particle data to devices...");
      fflush(stdout);
      cuda_point_push();
      //printf("done.\n");
      fflush(stdout);

      count_mem();

      // initialize ParaView VTK output PVD file
      if(runrestart != 1) {
//        #ifdef DEBUG
//          init_VTK_ghost();
//        #else
          if(rec_paraview_dt > 0) {
            init_VTK();
          }
//        #endif
      }

      real rec_flow_field_ttime_out = 0.;
      real rec_paraview_ttime_out = 0.;
      real rec_point_particle_ttime_out = 0.;
      real rec_restart_ttime_out = 0.;

      // set up point_particles
      cuda_build_cages();
      cuda_point_pull();

      // run restart if requested
      if(runrestart == 1) {
        //printf("\nRestart requested.\n\n");
        //printf("Reading restart file...");
        fflush(stdout);
        in_restart();
        //printf("done.\n");
        fflush(stdout);
        //printf("Copying host domain data to devices...");
        fflush(stdout);
        cuda_dom_push();
        //printf("done.\n");
        fflush(stdout);
        //printf("Copying host point_particle data to devices...");
        fflush(stdout);
        cuda_point_push();
        //printf("done.\n");
        fflush(stdout);
      }

      #ifdef DEBUG
        // write config to screen
        //printf("\n=====DEBUG");
        //printf("================================");
        //printf("======================================\n");
        fflush(stdout);
        cuda_dom_pull();
        cuda_point_pull();
        domain_show_config();
        points_show_config();
        //printf("========================================");
        //printf("========================================\n\n");
        fflush(stdout);
      #endif

      #ifdef TEST // run test code
        // ** note that some of these work only for DEV RANGE 0 0 **
        // test CUDA functionality
        //printf("\n=====TEST");
        //printf("=================================");
        //printf("======================================\n");
        fflush(stdout);
        rec_flow_field_stepnum_out = -1;
        rec_paraview_stepnum_out = -1;
        rec_point_particle_stepnum_out = -1;
        //rec_restart_stepnum_out = -1;
        rec_precursor_stepnum_out = -1;
        cuda_point_pull();
        //cuda_BC_test();
        //cuda_U_star_test_exp();
        //cuda_U_star_test_cos();
        //cuda_project_test();
        //cuda_quad_interp_test();
        cuda_lamb_test();
        //printf("========================================");
        //printf("========================================\n\n");
        fflush(stdout);

      #else // run simulation
        // begin simulation
        //printf("\n=====BLUEBOTTLE");
        //printf("===========================");
        //printf("======================================\n");
        fflush(stdout);

//Test1();

        // get initial dt; this is an extra check for the SHEAR initialization
        dt = cuda_find_dt();

        // share this with the precursor domain
        expd_compare_dt(np, status);

        // update the boundary condition config info to share with precursor
        expd_update_BC(np, status);
               
        // apply boundary conditions to field variables
        cuda_dom_BC();
       
        // write initial fields
        if(runrestart != 1) {
          cuda_dom_pull();
          cuda_point_pull();

//        #ifdef DEBUG
//          out_VTK_ghost();
//        #else
          if(rec_flow_field_dt > 0) {
            //printf("Writing flow field file t = %e...", ttime);
            fflush(stdout);
            cgns_grid();
            cgns_flow_field(rec_flow_field_dt);
            rec_flow_field_stepnum_out++;
            //printf("done.               \n");
            fflush(stdout);
          }
          if(rec_point_particle_dt > 0) {
            //printf("Writing point_particle file t = %e...", ttime);
            fflush(stdout);
            cgns_point_particles(rec_point_particle_dt);
//            recorder_lamb("lamb.rec", 0);
            rec_point_particle_stepnum_out++;
            //printf("done.               \n");
            fflush(stdout);
          }
          if(rec_paraview_dt > 0) {
            //printf("Writing ParaView file %d (t = %e)...", rec_paraview_stepnum_out, ttime);
            fflush(stdout);
            out_VTK();
            rec_paraview_stepnum_out++;
            //printf("done.               \n");
            fflush(stdout);
          }

//        #endif
        }

        /******************************************************************/
        /** Begin the main timestepping loop in the experimental domain. **/
        /******************************************************************/
        while(ttime <= duration) {
          ttime += dt;
          rec_flow_field_ttime_out += dt;
          rec_paraview_ttime_out += dt;
          rec_point_particle_ttime_out += dt;
          rec_restart_ttime_out += dt;
          stepnum++;
          //printf("FLOW: Time = %e of %e (dt = %e).\n", ttime, duration, dt);
          fflush(stdout);

          cuda_compute_forcing();
          compute_vel_BC();
          // update the boundary condition config info and share with precursor
          expd_update_BC(np, status);

          // TODO: save work by rebuilding only the cages that need to be rebuilt
          cuda_build_cages();

          int iter = 0;
   //       real iter_err = FLT_MAX;
  //        while(iter_err > lamb_residual) {  // iterate for Lamb's coefficients
            #ifndef BATCHRUN
              //printf("  Iteration %d: ", iter);
              fflush(stdout);
            #endif

            // solve for U_star
            cuda_U_star_2();
            // apply boundary conditions to U_star
  
            cuda_dom_BC_star();
            // enforce solvability condition
            cuda_solvability();
  
            cuda_dom_BC_star();
            // solve for pressure
           

//Test();

	    cuda_PP_bicgstab(rank);
            cuda_dom_BC_p();
            // solve for U
            cuda_project();
  
            cuda_dom_BC();

            // update Lamb's coefficients
            /*
            cuda_move_points_sub();
            cuda_Lamb();

            #ifdef STEPS // force no sub-timestep iteration
              iter_err = -1;
            #else
              // check error between this set of coefficients and previous set
              // of coefficients
              iter_err = cuda_lamb_err();
              // TODO: write error to lamb.rec
            #endif
          
            
            #ifndef BATCHRUN
              //printf("Error = %f\r", iter_err);
            #endif
            */
            
     //       iter++;
            // check iteration limit
  /*          if(iter == lamb_max_iter) {
              lambflag = !lambflag;
              //printf("Reached the maximum number of Lamb's");
              //printf(" coefficient iterations.");
              //printf(" Ending simulation.\n");
              break;
            }
   
         }

          //printf("  The Lamb's coefficients converged in");
          //printf(" %d iterations.\n", iter);
 */
       //   if(!lambflag) {
          
            // store u, conv, and coeffs for use in next timestep
            cuda_store_u();
        
       /*     if(npoints > 0)
              cuda_store_coeffs();
*/
            // compute div(U)
            //cuda_div_U();

            // update point_particle position
            if(npoints>0)
		{
		    cuda_flow_stress();
	            cuda_move_points();
		}
            // compute next timestep size
            dt0 = dt;
            dt = cuda_find_dt();

            // compare this timestep size to that in the precursor and
            // and synchronize the result
            expd_compare_dt(np, status);

   /*
       } else {
            return EXIT_FAILURE;
          }
*/
          if(rec_restart_dt > 0) {
            if(rec_restart_ttime_out >= rec_restart_dt) {
              //printf("  Writing restart file (t = %e)...", ttime);
              fflush(stdout);
              cuda_dom_pull();
              cuda_point_pull();
              out_restart();
              //printf("done.               \n");
              fflush(stdout);
              rec_restart_ttime_out = rec_restart_ttime_out - rec_restart_dt;
            }
          }

          if(rec_flow_field_dt > 0) {
            if(rec_flow_field_ttime_out >= rec_flow_field_dt) {
              // pull back data and write fields
              cuda_dom_pull();
              #ifndef BATCHRUN
                //printf("  Writing flow field file t = %e...                  \r", ttime);
                fflush(stdout);
              #endif
              cgns_flow_field(rec_flow_field_dt);
              //printf("  Writing flow field file t = %e...done.\n", ttime);
              fflush(stdout);
              rec_flow_field_ttime_out = rec_flow_field_ttime_out
                - rec_flow_field_dt;
              rec_flow_field_stepnum_out++;
            }
          }
          if(rec_paraview_dt > 0) {
            if(rec_paraview_ttime_out >= rec_paraview_dt) {
              // pull back data and write fields
              cuda_dom_pull();
              cuda_point_pull();
              #ifndef BATCHRUN
                //printf("  Writing ParaView output file");
                //printf(" %d (t = %e)...                  \r",rec_paraview_stepnum_out, ttime);
                fflush(stdout);
              #endif
//              #ifdef DEBUG
//                out_VTK_ghost();
//              #else
                out_VTK();
//              #endif
              //printf("  Writing ParaView file %d (t = %e)...done.\n", rec_paraview_stepnum_out, ttime);
              rec_paraview_stepnum_out++;
              fflush(stdout);
              rec_paraview_ttime_out = rec_paraview_ttime_out - rec_paraview_dt;
            }
          }
          if(rec_point_particle_dt > 0) {
            if(rec_point_particle_ttime_out >= rec_point_particle_dt) {
              // pull back data and write fields
              cuda_point_pull();
              #ifndef BATCHRUN
                //printf("  Writing point_particle file t = %e...                  \r", ttime);
                fflush(stdout);
              #endif
               
              cgns_point_particles(rec_point_particle_dt);
/*
              #ifdef DEBUG
      //          recorder_lamb("lamb.rec", iter);
              #else
                cgns_point_particles(rec_point_particle_dt);
      //         recorder_lamb("lamb.rec", iter);
              #endif
  */

            //printf("  Writing point_particle file t = %e...done.\n", ttime);
              fflush(stdout);
              rec_point_particle_ttime_out = rec_point_particle_ttime_out - rec_point_particle_dt;
              rec_point_particle_stepnum_out++;
            }
          }

          // check for blow-up condition
          if(dt < 1e-20) {
            //printf("The solution has diverged.  Ending simulation.              \n");
            break;
          }
        }

        //printf("========================================");
        //printf("========================================\n\n");
        fflush(stdout);
      #endif

      // clean up devices
      //printf("Cleaning up domain data on devices...");
      fflush(stdout);
      cuda_dom_free();
      //printf("done.     \n");
      fflush(stdout);
      //printf("Cleaning up point_particle data on devices...");
      fflush(stdout);
      cuda_point_free();
      //printf("done.\n");
      fflush(stdout);

      // clean up host
      //printf("Cleaning up point_particles...");
      fflush(stdout);
      points_clean();
      //printf("done.\n");
      fflush(stdout);
      //printf("Cleaning up domain...");
      fflush(stdout);
      domain_clean();
      //printf("done.\n");
      fflush(stdout);

      //printf("\n...bluebottle_0.1 done.\n\n");
    }
  } else {
    int turb = 1;   // boolean
    MPI_Status status; // for MPI communication

    // parse command-line arguments
    // if none given, run normally
    // if -s given, run seeder program only
    // if anything else given, exit
    // this code comes from K&R, p. 117
    int argin;
    int runrestart = 0;
    while(--argc > 0 && (*++argv)[0] == '-') {
      while((argin = *++argv[0])) {
        switch(argin) {
          case 'r':
            runrestart = 1;
            break;
          default:
            runrestart = 2;
            //printf("bluebottle: illegal option %c\n", argin);
            argc = 0;
            break;
        }
      }
    }

    // read simulation input configuration file
    recorder_read_config();
    turb_read_input();
    points_read_input(turb);

    //printf("TURB: Using devices %d through %d.\n\n", dev_start, dev_end);
    fflush(stdout);

    if(runrestart != 1) {
      // start BICGSTAB recorder
      recorder_bicgstab_init("solver_turb.rec");
    }

    // initialize the domain
    int domain_init_flag = domain_init_turb();
    if(domain_init_flag == EXIT_FAILURE) {
      //printf("\nThe number of devices in DEV RANGE is insufficient\n");
      //printf("for the given turbulence domain decomposition.  Exiting now.\n");
      return EXIT_FAILURE;
    }

    // receive the boundary condition config info from MASTER
    prec_init_BC(np, rank, status);

    // initialize the point_particles
    int points_init_flag = points_init();
    if(points_init_flag == EXIT_FAILURE) {
      //printf("\nThe initial point_particle configuration is not allowed.\n");
      return EXIT_FAILURE;
    }

    // allocate device memory
    cuda_dom_malloc();
    cuda_point_malloc();

    // copy host data to devices
    cuda_dom_push();
    cuda_dom_turb_planes_push(bc_flow_configs);
    cuda_point_push();

    //count_mem();

    // initialize ParaView VTK output PVD file
    if(rec_precursor_dt > 0) {
      init_VTK_turb();
    }

    real rec_precursor_ttime_out = 0.;
    real rec_restart_ttime_out = 0.;

    cuda_build_cages();
    cuda_point_pull();

    // run restart if requested
    if(runrestart == 1) {
      //printf("\nRestart requested.\n\n");
      //printf("Reading restart file...");
      fflush(stdout);
      in_restart_turb();
      //printf("done.\n");
      fflush(stdout);
      //printf("Copying host domain data to devices...");
      fflush(stdout);
      cuda_dom_push();
      //printf("done.\n");
      fflush(stdout);
      //printf("Copying host point_particle data to devices...");
      fflush(stdout);
      cuda_point_push();
      //printf("done.\n");
      fflush(stdout);
    }

    // initialize timestep size since turbulent velocity is nonzero
    dt = cuda_find_dt();

    // share this dt with the experimental domain
    prec_compare_dt(np, rank, status);

    // update the boundary conditions according to the experimental domain
    prec_update_BC(np, rank, status);

    // begin simulation
    // apply boundary conditions to field variables
    cuda_dom_BC();

    // write initial fields
    if(rec_precursor_dt > 0) {
      cuda_dom_pull();
      //printf("Writing precursor file %d (t = %e)...", rec_precursor_stepnum_out, ttime);
      fflush(stdout);
      out_VTK_turb();
      rec_precursor_stepnum_out++;
      //printf("done.               \n");
      fflush(stdout);
    }

    /***************************************************************/
    /** Begin the main timestepping loop in the precursor domain. **/
    /***************************************************************/
    while(ttime <= duration) {
      ttime += dt;
      rec_precursor_ttime_out += dt;
      rec_restart_ttime_out += dt;
      stepnum++;
      //printf("TURB: Time = %e of %e (dt = %e).\n", ttime, duration, dt);

      cuda_compute_forcing();
      cuda_compute_turb_forcing();
      compute_vel_BC();

      // solve for U_star
      cuda_U_star_2();
      // apply boundary conditions to U_star
      cuda_dom_BC_star();
      // force solvability condition
      cuda_solvability();
      cuda_dom_BC_star();
      // solve for pressure
      cuda_PP_bicgstab(rank);
      cuda_dom_BC_p();
      // solve for U
      cuda_project();
      // apply boundary conditions to field variables
      cuda_dom_BC();

//move velocity of points as a beganing
/*
if(npoints>0){
cuda_move_points_sub();
}
*/
      cuda_store_u();

//move points finally
/*
if(npoints>0){
cuda_move_points();
}
*/
      // compute next timestep size
      dt0 = dt;
      dt = cuda_find_dt();
     
      // check for blow-up condition
      if(dt < 1e-20) {
        //printf("The solution has diverged.  Ending simulation.              \n");
        return EXIT_FAILURE;
      }

      // communicate the boundary condition with the experimental domain
      prec_send_BC(np, rank, status);

      if(rec_precursor_dt > 0) {
        if(rec_precursor_ttime_out >= rec_precursor_dt) {
          // pull back data and write fields
          cuda_dom_pull();
          #ifndef BATCHRUN
            //printf("  Writing precursor output file");
            //printf(" %d (t = %e)...                  \r", rec_precursor_stepnum_out, ttime);
            fflush(stdout);
          #endif
          #ifdef DEBUG
            out_VTK_ghost();
          #else
            out_VTK_turb();
          #endif
          //printf("  Writing precursor file %d (t = %e)...done.\n", rec_precursor_stepnum_out, ttime);
          rec_precursor_stepnum_out++;
          fflush(stdout);
          rec_precursor_ttime_out = rec_precursor_ttime_out - rec_precursor_dt;
        }
      }
      if(rec_restart_dt > 0) {
        if(rec_restart_ttime_out >= rec_restart_dt) {
          //printf("  Writing precursor restart file (t = %e)...", ttime);
          fflush(stdout);
          cuda_dom_pull();
          out_restart_turb();
          //printf("done.               \n");
          fflush(stdout);
          rec_restart_ttime_out = rec_restart_ttime_out - rec_restart_dt;
        }
      }
    }

    // clean up devices
    cuda_dom_free();
    cuda_point_free();

    // clean up host
    points_clean();
    domain_clean();
  }

  // finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}
