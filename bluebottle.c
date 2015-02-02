#include <mpi.h>
#include "bluebottle.h"
#include "point.h"
#include "scalar.h"
#include "precursor.h"
#include "recorder.h"
//#include "cuda_point.h"
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

float *GaussianKernel; //Gaussian kernel weight contributed by each point particle 
int *_DomInfo;

//Temp array for particle integration
real **ug,**vg,**wg;//device pointer of the fluid velocity at the particle position
real **posX,**posY,**posZ;//device pointer of the particle position
real **posXold,**posYold,**posZold;//device pointer of the particle position
real **lptSourceVal; //Source from each point particle 
real **lptSourceValOld; //Source from each point particle 

real **lpt_stress_u,**lpt_stress_v,**lpt_stress_w;//device pointer of the fluid velocity at the particle position
real **scg;//device pointer of the fluid scalar at the particle position
real **Weight; //Gaussian kernel weight contributed by each point particle 
real **Ksi; //Gaussian kernel weight contributed by each point particle 
int  **cellStart;
int  **cellEnd;
int  **gridParticleIndex;
int  **gridParticleHash;



//Temp array for scalar integration

real C_add;
real C_stress;
real C_drag;

real lpt_twoway;
real sc_twoway;
/*
real **_omega_x;
real **_omega_y;
real **_omega_z;
*/

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

//scalar in bluebottle
real *sc;
real *sc0;
real **_sc;
real **_sc0;
real *diff0_sc;
real **_diff0_sc;
real *conv0_sc;
real **_conv0_sc;
real *diff_sc;
real **_diff_sc;
real *conv_sc;
real **_conv_sc;
real *sc_WE;
real *sc_SN;
real *sc_BT;
real **_sc_WE;
real **_sc_SN;
real **_sc_BT;

real *scSrc;
real **_scSrc;
real *epsp;
real **_epsp;

real dt_sc;
real dt_done;
real dt_try;
real dt0_try;
real ttime_done;


real n2;
real n3;


scBC sc_bc;
int sc_init_cond;
real sc_eq;
real sc_init_percent;

real **_lpt_mom_source_x;
real **_lpt_mom_source_y;
real **_lpt_mom_source_z;

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
int rec_scalar_stepnum_out;
int rec_precursor_stepnum_out;
real rec_flow_field_dt;
int rec_flow_field_vel;
int rec_flow_field_p;
int rec_flow_field_phase;
real rec_paraview_dt;
real rec_point_particle_dt;
real rec_scalar_field_dt;
//real rec_scalar_field;
real rec_restart_dt;

int rec_restart_stop;
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
//add by shigan
real DIFF;
real DIFF_eq;

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


real pid_int;
real pid_back;
real Kp;
real Ki;
real Kd;




int main(int argc, char *argv[]) {

  int np = 0;
  int rank = 0;
  int restart_stop=0;

  // set up MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  time_t startwalltime = time(NULL);
  time_t timestepwalltime = time(NULL);
  time_t diffwalltime = difftime(timestepwalltime, startwalltime);


  if(rank == MASTER) {
    MPI_Status status; // for MPI communication


    int argin;
    int runseeder = 0;
    int runrestart = 0;
    int runPointScalar_restart = 0;
    int runBubbleScalar_restart = 0;
    int runScalar_restart = 0;

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
//'p' represents inject particle and scalar simutaneously
          case 'p':
	    runPointScalar_restart = 1;
	    runrestart=1;
	    break;
//'b' represents inject bubble and scalar simutaneously
          case 'b':
	    runBubbleScalar_restart = 1;
	    runrestart=1;
	    break;
//'c' represents inject scalar only
          case 'c':
	    runScalar_restart = 1;
	    runrestart=1;
	    break;
          default:
            runseeder = 2;
            runrestart = 2;
	    runPointScalar_restart=2;
	    runBubbleScalar_restart = 2;
	    runScalar_restart=2;
            argc = 0;
            break;
        }
      }
    }
    
    NP = 0;
    radius = -1.;
    order = -1.;
     if(runseeder == 1) {
      NP = atoi(argv[0]);
      radius = atof(argv[1]);
      density = atof(argv[2]);
      order = atoi(argv[3]);
      translating = atoi(argv[4]);
      rotating = atoi(argv[5]);
      seeder(NP, radius, density, order, translating, rotating);
      return EXIT_SUCCESS;
    } else if(runrestart == 1 && argc > 0) {
      printf("Usage restart simulation: bluebottle -r\n");
      return EXIT_FAILURE;
    } else if(runseeder == 2) {
      return EXIT_FAILURE;
    } else if(runrestart == 2) {
      return EXIT_FAILURE;
    } else {
	//printf("\nReading config files\n");
      // read recorder config file
      recorder_read_config();

//printf("\nrec_scalar_field_dt %f\n",rec_scalar_field_dt);fflush(stdout);

     // domain_read_input();
      turb_read_input();
      points_read_input();
      scalar_read_input();

 /********* Messy hack for taking advantage of CUDA_VISIBLE_DEVICES
       ********* treatment in SLURM resource manager. */

      // read the environment variable
      char *cvdin;
      cvdin = getenv("CUDA_VISIBLE_DEVICES");
      if(cvdin != NULL) {
        // number of devices
        int n_CUDA_VISIBLE_DEVICES = 0.5*(strlen(cvdin)+1.);
        // list of devices
        int *CUDA_VISIBLE_DEVICES = malloc(n_CUDA_VISIBLE_DEVICES*sizeof(int));
        // fill list of devices assuming single-character separation
        int j = 0;
        for(int i = 0; i < 2*n_CUDA_VISIBLE_DEVICES-1; i+=2) {
          CUDA_VISIBLE_DEVICES[j] = cvdin[i] - '0';
          j++;
        }
        // use the first device available (devices are re-indexed by CUDA)
        if(n_CUDA_VISIBLE_DEVICES > 0) {
          dev_start = 0;
          dev_end = 0;
        } else { // exit if things aren't just right
          printf("Environment variable CUDA_VISIBLE_DEVICES is empty:\n");
          printf("  a. To use the config files to specify the device number,\n");
          printf("     type 'unset CUDA_VISIBLE_DEVICES'\n");
          printf("  b. To use CUDA_VISIBLE_DEVICES to specify the device number,\n");
          printf("     type 'export CUDA_VISIBLE_DEVICES=N1,N2,...',\n");
          printf("     where N1,N2 are comma-separated device numbers\n");
          exit(EXIT_FAILURE);
        }
      }


      if(runrestart != 1||runPointScalar_restart==1||runBubbleScalar_restart==1||runScalar_restart==1) {
        recorder_bicgstab_init("solver_flow.rec");
        recorder_bicgstab_init("solver_helmholtz_scalar.rec");
      }
/*
      int domain_init_flag = domain_init();
      if(domain_init_flag == EXIT_FAILURE) {
        printf("\nThe number of devices in DEV RANGE is insufficient\n");
        printf("for the given domain decomposition.  Exiting now.\n");
        return EXIT_FAILURE;
      }
*/

 // initialize the domain
    //printf("\nInitializing domain...\n");
    int domain_init_flag = domain_init_turb();
    if(domain_init_flag == EXIT_FAILURE) {
      printf("\nThe number of devices in DEV RANGE is insufficient\n");
      printf("for the given turbulence domain decomposition.  Exiting now.\n");
      return EXIT_FAILURE;
    }


      // set up the boundary condition config info to send to precursor
    //printf("\nInitializing particles...\n");
    
      // initialize the point_particles
      int points_init_flag = points_init();
      fflush(stdout);
      if(points_init_flag == EXIT_FAILURE) {
        printf("\nThe initial point_particle configuration is not allowed.\n");
        return EXIT_FAILURE;
      }

  // initialize the scalar 
      int scalar_init_flag = scalar_init();
      fflush(stdout);
      if(scalar_init_flag == EXIT_FAILURE) {
        printf("\nThe initial scalar configuration is not allowed.\n");
        return EXIT_FAILURE;
      }

      // allocate device memory
      cuda_dom_malloc();

      cuda_point_malloc();
      cuda_scalar_malloc();
fflush(stdout);

      // copy host data to devices
      cuda_dom_push();
      cuda_point_push();
      cuda_scalar_push();

      count_mem();

      // initialize ParaView VTK output PVD file
      if(runrestart != 1) {
          if(rec_paraview_dt > 0) {
            init_VTK();
          }
      }

      real rec_flow_field_ttime_out = 0.;
      real rec_paraview_ttime_out = 0.;
      real rec_point_particle_ttime_out = 0.;
      real rec_restart_ttime_out = 0.;
      real rec_scalar_ttime_out = 0.;

      // set up point_particles
      cuda_build_cages();
      ////cuda_point_pull();//why pull at the beginning???
      ////cuda_scalar_pull();

      // run restart if requested
      if(runrestart == 1) {
        printf("\nRestart requested.\n\n");
        printf("Reading restart file...");
        fflush(stdout);
        in_restart();
        cuda_dom_push();
        cuda_point_push();
        cuda_scalar_push();

	if(runPointScalar_restart==1) points_scalar_inject();
	if(runBubbleScalar_restart==1) bubble_scalar_inject();
	if(runScalar_restart==1) scalar_inject();			

 	if(ttime >= duration) {
	 printf("\n...simulation completed.\n");
	 restart_stop = 1;
	 }
      }

      #ifdef DEBUG
        cuda_dom_pull();
        cuda_point_pull();
	cuda_scalar_pull();
        domain_show_config();
        points_show_config();
        scalar_show_config();
      #endif


       // run simulation
        // get initial dt; this is an extra check for the SHEAR initialization
        dt = cuda_find_dt();
        dt_sc = cuda_find_dt_sc(dt);
	real dt_point;
        dt_point = cuda_find_dt_points(dt_sc);
	dt_sc=dt_point;

        // share this with the precursor domain
        //expd_compare_dt(np, status);

        // update the boundary condition config info to share with precursor
        //expd_update_BC(np, status);
               
        // apply boundary conditions to field variables
        compute_scalar_BC();
        compute_vel_BC();
        cuda_dom_BC();
//      cuda_scalar_BC();


        // write initial fields
        if(runrestart != 1) {
          cuda_dom_pull();
          cuda_point_pull();
          cuda_scalar_pull();
          if(rec_flow_field_dt > 0) {
            cgns_grid();
            cgns_flow_field(rec_flow_field_dt);
            rec_flow_field_stepnum_out++;
          }
          if(rec_point_particle_dt > 0) {
            cgns_point_particles(rec_point_particle_dt);
            rec_point_particle_stepnum_out++;
          }

	 if(rec_scalar_field_dt > 0) {
            cgns_scalar_field(rec_scalar_field_dt);
            rec_scalar_stepnum_out++;
          }

          if(rec_paraview_dt > 0) {
            out_VTK();
            rec_paraview_stepnum_out++;
          }
        }


       printf("\nBegan iteration...\n");
        /******************************************************************/
        /** Begin the main timestepping loop in the experimental domain. **/
        /******************************************************************/
        while(ttime <= duration) {
          ttime += dt;
          rec_flow_field_ttime_out += dt;
          rec_paraview_ttime_out += dt;
          rec_point_particle_ttime_out += dt;
          rec_restart_ttime_out += dt;
          rec_scalar_ttime_out += dt;
          stepnum++;
          printf("FLOW: Time = %e of %e (dt = %e).\n", ttime, duration, dt);
          fflush(stdout);

 // time_t  t1=time(NULL);

          cuda_compute_forcing();

if(turbA>EPSILON)          cuda_compute_turb_forcing();//add
if(npoints>0&&lpt_twoway>0)    	   lpt_point_twoway_forcing();
		
//  time_t  t2=time(NULL);
//printf("\n21 %f\n",difftime(t2,t1));

          compute_vel_BC();
          // update the boundary condition config info and share with precursor
          //expd_update_BC(np, status);

          // TODO: save work by rebuilding only the cages that need to be rebuilt
          cuda_build_cages();

          int iter = 0;
            #ifndef BATCHRUN
              printf("  Iteration %d: ", iter);
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
           
	    cuda_PP_bicgstab(rank);
            cuda_dom_BC_p();
            // solve for U
            cuda_project();
  
            cuda_dom_BC();

//For nvvp profiler
  cudaProfilerStart();

	//Calculate fluid stress on particles at time t_n
            if(npoints>0)  cuda_flow_stress();


      // update point_particle position in substep
	     if(dt0<0)     cuda_scalar_BC();
dt_done=0.0;
	printf("\n");
while(dt_done<dt)
{
	if(dt_sc<dt-dt_done) dt_try=dt_sc;
	else dt_try=dt-dt_done;
	ttime_done=ttime+dt_done;

          printf("SCLAR and Particle: Time = %e of %e (dt = %e).\n", ttime_done, duration, dt_try);
          fflush(stdout);

//printf("\nrank %d\n",rank);
//          fflush(stdout);

	//move points explicitly at every sub-timestep
            if(npoints>0)   cuda_move_points();
	//add particle back reaction in momentum to the flow during the substep, then add it to flow source term
            if(lpt_twoway>0&&npoints>0)   lpt_point_twoway_momentum();
            compute_scalar_BC();
            cuda_scalar_advance();
            cuda_scalar_BC();
    
	    cuda_store_scalar();

	//Store old time step
	dt0_try=dt_try;

	// Increment
        dt_done = dt_done+dt_try;

  }      

//For nvvp profiler
 cudaProfilerStop();
      // store u, conv, and coeffs for use in next timestep
            cuda_store_u();
        
            // compute next timestep size
            dt0 = dt;
            dt = cuda_find_dt();
            dt_sc = cuda_find_dt_sc(dt);
            dt_point = cuda_find_dt_points(dt_sc);
	    dt_sc=dt_point;
            // compare this timestep size to that in the precursor and
            // and synchronize the result
            //expd_compare_dt(np, status);


// write a restart file and exit when the time is appropriate
          timestepwalltime = time(NULL);
          diffwalltime = difftime(timestepwalltime, startwalltime);
          if(rec_restart_dt > 0) {
            if((real)diffwalltime/60. > rec_restart_dt) {
              printf("  Writing restart file (t = %e)...", ttime);
              fflush(stdout);
              cuda_dom_pull();
              cuda_point_pull();
              cuda_scalar_pull();

              out_restart();
              printf("done.               \n");
              fflush(stdout);
              rec_restart_ttime_out = rec_restart_ttime_out - rec_restart_dt;
              startwalltime = time(NULL);
              if(rec_restart_stop)
                break; // exit!
            }
          }

 if(rec_restart_dt > 0 && ttime >= duration && !restart_stop) {
 printf(" Writing final restart file (t = %e)...", ttime);
 fflush(stdout);
 cuda_dom_pull();
 cuda_point_pull();
 cuda_scalar_pull();

 out_restart();
 printf("done. \n");
 fflush(stdout);
 rec_restart_ttime_out = rec_restart_ttime_out - rec_restart_dt;
 startwalltime = time(NULL);
 }
          if(rec_flow_field_dt > 0) {
            if(rec_flow_field_ttime_out >= rec_flow_field_dt) {
              // pull back data and write fields
              cuda_dom_pull();
              #ifndef BATCHRUN
                printf("  Writing flow field file t = %e...                  \r", ttime);
                fflush(stdout);
              #endif
              cgns_flow_field(rec_flow_field_dt);
              printf("  Writing flow field file t = %e...done.\n", ttime);
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
              cuda_scalar_pull();

              #ifndef BATCHRUN
                printf("  Writing ParaView output file");
                printf(" %d (t = %e)...                  \r",rec_paraview_stepnum_out, ttime);
                fflush(stdout);
              #endif
                out_VTK();
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
                printf("  Writing point_particle file t = %e...                  \r", ttime);
                fflush(stdout);
              #endif
               
              cgns_point_particles(rec_point_particle_dt);
              rec_point_particle_ttime_out = rec_point_particle_ttime_out - rec_point_particle_dt;
              rec_point_particle_stepnum_out++;
            }
          }
//printf("\nscalar_dt %f %f\n",rec_scalar_field_dt,rec_scalar_ttime_out);
//                fflush(stdout);
          if(rec_scalar_field_dt > 0) {
            if(rec_scalar_ttime_out >= rec_scalar_field_dt) {
              // pull back data and write fields
              cuda_scalar_pull();
              #ifndef BATCHRUN
                printf("  Writing scalar file t = %e...                  \n", ttime);
                fflush(stdout);
              #endif               
              cgns_scalar_field(rec_scalar_field_dt);
              rec_scalar_ttime_out = rec_scalar_ttime_out - rec_scalar_field_dt;//???????
              rec_scalar_stepnum_out++;
            }
          }

          // check for blow-up condition
          if(dt < 1e-20) {
            printf("The solution has diverged.  Ending simulation.              \n");
            break;
          }
        }

           fflush(stdout);

      // clean up devices
      cuda_dom_free();
      cuda_point_free();
      cuda_scalar_free();

      // clean up host
      points_clean();
      domain_clean();
      scalar_clean();

    }
  }
  MPI_Finalize();


//For nvvp profiler to work
  cudaDeviceReset();

  if(restart_stop) return EXIT_FAILURE;
  else return EXIT_SUCCESS;


}
