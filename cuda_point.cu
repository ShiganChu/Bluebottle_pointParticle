//#include "cuda_point.h"
//#include "point.h"

#include <cuda.h>
#include <helper_cuda.h>

//#include "bluebottle.h"
#include "cuda_point.h"

extern "C"
void cuda_flow_stress()
{
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));
//printf("\ndev in cuda_flow_stress %d %d\n",dev,dev_start);
    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

   // u-component
    if(dom[dev].Gfx.jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfx.jnb + 2;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfx.knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfx.knb + 2;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) (threads_y-2));
    blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) (threads_z-2));



    dim3 dimBlocks_u(threads_y, threads_z);
    dim3 numBlocks_u(blocks_y, blocks_z);
    stress_u<<<numBlocks_u, dimBlocks_u>>>(rho_f, nu,_u[dev],_p[dev],_p0[dev], _stress_u[dev], _dom[dev],_flag_u[dev],dt,dt0);
fflush(stdout);
/*
printf("\nGfx %f\n",dom[dev].Gfx.jnb);
printf("\nthreads_y threads_z %d %d %d %d\n",threads_y, threads_z,blocks_y, blocks_z);
fflush(stdout);
*/

// v-component
    if(dom[dev].Gfy.knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfy.knb + 2;
    else
      threads_z = MAX_THREADS_DIM;

    if(dom[dev].Gfy.inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfy.inb + 2;
    else
      threads_x = MAX_THREADS_DIM;
  
    blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) (threads_z-2));
    blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) (threads_x-2));

    dim3 dimBlocks_v(threads_z, threads_x);
    dim3 numBlocks_v(blocks_z, blocks_x);

    stress_v<<<numBlocks_v, dimBlocks_v>>>(rho_f, nu,_v[dev],_p[dev],_p0[dev], _stress_v[dev], _dom[dev],_flag_v[dev],dt,dt0);
fflush(stdout);

// w-component
    if(dom[dev].Gfz.inb < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfz.inb + 2;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfz.jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfz.jnb + 2;
    else
      threads_y = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) (threads_x-2));
    blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) (threads_y-2));

    dim3 dimBlocks_w(threads_x, threads_y);
    dim3 numBlocks_w(blocks_x, blocks_y);

    stress_w<<<numBlocks_w, dimBlocks_w>>>(rho_f, nu,_w[dev],_p[dev],_p0[dev], _stress_w[dev], _dom[dev],_flag_w[dev],dt,dt0);
fflush(stdout);
 }

}


extern "C"
void cuda_move_points()
{
   // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

    int threads = MAX_THREADS_1D;
    int blocks = (int)ceil((real) npoints / (real) threads);

    dim3 dimBlocks(threads);
    dim3 numBlocks(blocks);


//printf("\nthreads blocks %d %d \n",threads,blocks);
//fflush(stdout);

    
    if(npoints > 0) {
      // do collision forcing
      /** if there are n point_particles in a close group, repeat this n times **/
    /*
      real *forces;
      checkCudaErrors(cudaMalloc((void**) &forces, 3*npoints*sizeof(real)));
      gpumem += 3 * npoints * sizeof(real);
      real *moments;
      checkCudaErrors(cudaMalloc((void**) &moments, 3*npoints*sizeof(real)));
      gpumem += 3 * npoints * sizeof(real);
      real eps = 0.1 * (Dom.dx + Dom.dy + Dom.dz) / 3.;
*/
real *ug,*vg,*wg;//device pointer of the fluid velocity at the particle position
real *lpt_stress_u,*lpt_stress_v,*lpt_stress_w;//device pointer of the fluid velocity at the particle position
real *scg;//device pointer of the fluid scalar at the particle position


real *ug0,*vg0,*wg0;//device pointer of the fluid velocity at the particle position for the previous time step
real *conv_ug,*conv_vg,*conv_wg;//device pointer of the convective part for fluid velocity at the particle position

      checkCudaErrors(cudaMalloc((void**) &ug, npoints*sizeof(real)));
      checkCudaErrors(cudaMalloc((void**) &vg, npoints*sizeof(real)));
      checkCudaErrors(cudaMalloc((void**) &wg, npoints*sizeof(real)));
      gpumem += 3 * npoints * sizeof(real);

      checkCudaErrors(cudaMalloc((void**) &lpt_stress_u, npoints*sizeof(real)));
      checkCudaErrors(cudaMalloc((void**) &lpt_stress_v, npoints*sizeof(real)));
      checkCudaErrors(cudaMalloc((void**) &lpt_stress_w, npoints*sizeof(real)));
      gpumem += 3 * npoints * sizeof(real);

      checkCudaErrors(cudaMalloc((void**) &scg, npoints*sizeof(real)));
      gpumem +=  npoints * sizeof(real);



      checkCudaErrors(cudaMalloc((void**) &ug0, npoints*sizeof(real)));
      checkCudaErrors(cudaMalloc((void**) &vg0, npoints*sizeof(real)));
      checkCudaErrors(cudaMalloc((void**) &wg0, npoints*sizeof(real)));
      gpumem += 3 * npoints * sizeof(real);
      checkCudaErrors(cudaMalloc((void**) &conv_ug, npoints*sizeof(real)));
      checkCudaErrors(cudaMalloc((void**) &conv_vg, npoints*sizeof(real)));
      checkCudaErrors(cudaMalloc((void**) &conv_wg, npoints*sizeof(real)));
      gpumem += 3 * npoints * sizeof(real);


//for(int l = 0; l < 10; l++) {
   //   collision_init<<<numBlocks, dimBlocks>>>(_points[dev], npoints);
    
 /*
       for(int i = 0; i < npoints; i++) {
          collision_points<<<numBlocks, dimBlocks>>>(_points[dev], i,
            _dom[dev], eps, forces, moments, npoints, mu, bc);
        }
        spring_points<<<numBlocks, dimBlocks>>>(_points[dev], npoints);
        collision_walls<<<numBlocks, dimBlocks>>>(_dom[dev], _points[dev],
          npoints, bc, eps, mu);
      }
   */

   
//bc is bc.uTD etc. Make sure which BC this is. 
      

point_interp_init<<<numBlocks, dimBlocks>>>(npoints,_points[dev],ug,vg,wg,lpt_stress_u,lpt_stress_v,lpt_stress_w,scg);

interpolate_point_vel_Lag2<<<numBlocks, dimBlocks>>>(_u[dev],_v[dev],_w[dev],npoints,rho_f,nu,ug,vg,wg,_points[dev],_dom[dev],bc);
fflush(stdout);

interpolate_point_scalar_Lag2<<<numBlocks, dimBlocks>>>(npoints,_sc[dev],scg,_points[dev],_dom[dev]);
fflush(stdout);



/*
C_add=0.5;
C_stress=1;
C_drag=1;
*/

//get lpt_stress
//TODO _stress_u is not available near the boundary(set to 0), while lpt_stress can be interpolated on BC
if(C_stress>0||C_add>0) interpolate_point_vel_Lag2<<<numBlocks, dimBlocks>>>(_stress_u[dev],_stress_v[dev],_stress_w[dev],npoints,rho_f,nu,lpt_stress_u,lpt_stress_v,lpt_stress_w,_points[dev],_dom[dev],bc);


drag_points<<<numBlocks, dimBlocks>>>(_points[dev],npoints,
ug,vg,wg,
lpt_stress_u,lpt_stress_v,lpt_stress_w,scg,
rho_f,mu,g,gradP,
C_add, C_stress,C_drag,
sc_eq,DIFF);
fflush(stdout);

      move_points_a<<<numBlocks, dimBlocks>>>(_points[dev], npoints,dt_try);

drag_points<<<numBlocks, dimBlocks>>>(_points[dev],npoints,
ug,vg,wg,
lpt_stress_u,lpt_stress_v,lpt_stress_w,scg,
rho_f,mu,g,gradP,
C_add, C_stress,C_drag,
sc_eq,DIFF);
fflush(stdout);

      move_points_b<<<numBlocks, dimBlocks>>>(_dom[dev], _points[dev], npoints,dt_try);

fflush(stdout);




 /*
      checkCudaErrors(cudaFree(forces));
      checkCudaErrors(cudaFree(moments));
 */

      checkCudaErrors(cudaFree(ug));
      checkCudaErrors(cudaFree(vg));
      checkCudaErrors(cudaFree(wg));

      checkCudaErrors(cudaFree(scg));


      checkCudaErrors(cudaFree(lpt_stress_u));
      checkCudaErrors(cudaFree(lpt_stress_v));
      checkCudaErrors(cudaFree(lpt_stress_w));


      checkCudaErrors(cudaFree(ug0));
      checkCudaErrors(cudaFree(vg0));
      checkCudaErrors(cudaFree(wg0));
      checkCudaErrors(cudaFree(conv_ug));
      checkCudaErrors(cudaFree(conv_vg));
      checkCudaErrors(cudaFree(conv_wg));
   		}
  	
    }
}



extern "C"
void lpt_point_twoway_forcing()
{
 // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));
    int threads = MAX_THREADS_1D;
    int blocks = (int)ceil((real) npoints / (real) threads);
    int blocks_st = blocks*STENCIL3;

    dim3 dimBlocks(threads);
    dim3 numBlocks(blocks);

//blocks for initialization
    dim3 numBlocks_st(blocks_st);

      real *Ksi; //Source of scalar contributed by each point particle 
      checkCudaErrors(cudaMalloc((void**) &Ksi, npoints*STENCIL3*sizeof(real)));

array_init<<<numBlocks_st, dimBlocks>>>(Ksi,_dom[dev],npoints*STENCIL3,0.);

lpt_mollify_sc<<<numBlocks, dimBlocks>>>(npoints, _f_x[dev],_points[dev],_dom[dev],Ksi,1,0);
fflush(stdout);
   lpt_sum_ksi<<<numBlocks, dimBlocks>>>(npoints, _f_x[dev],_points[dev],_dom[dev],Ksi,1,0);
fflush(stdout);

lpt_mollify_sc<<<numBlocks, dimBlocks>>>(npoints, _f_y[dev],_points[dev],_dom[dev],Ksi,2,0);
   lpt_sum_ksi<<<numBlocks, dimBlocks>>>(npoints, _f_y[dev],_points[dev],_dom[dev],Ksi,2,0);

lpt_mollify_sc<<<numBlocks, dimBlocks>>>(npoints, _f_z[dev],_points[dev],_dom[dev],Ksi,3,0);
   lpt_sum_ksi<<<numBlocks, dimBlocks>>>(npoints, _f_z[dev],_points[dev],_dom[dev],Ksi,3,0);
 checkCudaErrors(cudaFree(Ksi));	
     }
}





extern "C"
void cuda_point_malloc(void)
{
  // allocate device memory on host
  _points = (point_struct**) malloc(nsubdom * sizeof(point_struct*));
  cpumem += nsubdom * sizeof(point_struct*);

  _flag_u = (int**) malloc(nsubdom * sizeof(int*));
  cpumem += nsubdom * sizeof(int);
  _flag_v = (int**) malloc(nsubdom * sizeof(int*));
  cpumem += nsubdom * sizeof(int);
  _flag_w = (int**) malloc(nsubdom * sizeof(int*));
  cpumem += nsubdom * sizeof(int);

  // allocate device memory on device
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

    checkCudaErrors(cudaMalloc((void**) &(_points[dev]),
      sizeof(point_struct) * npoints));
    gpumem += sizeof(point_struct) * npoints;

   checkCudaErrors(cudaMalloc((void**) &(_flag_u[dev]),
      sizeof(int) * dom[dev].Gfx.s3b));
    gpumem += sizeof(int) * dom[dev].Gfx.s3b;
    checkCudaErrors(cudaMalloc((void**) &(_flag_v[dev]),
      sizeof(int) * dom[dev].Gfy.s3b));
    gpumem += sizeof(int) * dom[dev].Gfy.s3b;
    checkCudaErrors(cudaMalloc((void**) &(_flag_w[dev]),
      sizeof(int) * dom[dev].Gfz.s3b));
    gpumem += sizeof(int) * dom[dev].Gfz.s3b;

  }
}

extern "C"
void cuda_point_push(void)
{
  // copy host data to device
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

    checkCudaErrors(cudaMemcpy(_points[dev], points, sizeof(point_struct) * npoints,cudaMemcpyHostToDevice));
  //  checkCudaErrors(cudaMemcpy(points,_points[dev],sizeof(point_struct) * npoints,cudaMemcpyDeviceToHost));
 }
}

extern "C"
void cuda_point_pull(void)
{
  // all devices have the same point_particle data for now, so just copy one of them
  checkCudaErrors(cudaMemcpy(points, _points[0], sizeof(point_struct) * npoints,
    cudaMemcpyDeviceToHost));

}

extern "C"
void cuda_point_free(void)
{
  // free device memory on device
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

    checkCudaErrors(cudaFree(_points[dev]));
    checkCudaErrors(cudaFree(_flag_u[dev]));
    checkCudaErrors(cudaFree(_flag_v[dev]));
    checkCudaErrors(cudaFree(_flag_w[dev]));
  

}

  free(_points);
  free(_flag_u);
  free(_flag_v);
  free(_flag_w);


}


extern "C"
void cuda_build_cages(void)
{
  cuda_point_pull();

  // parallelize over domains
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;


    // reset flag_u
    if(dom[dev].Gfx.jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfx.jnb;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gfx.kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfx.knb;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);

    dim3 dimBlocks_u(threads_y, threads_z);
    dim3 numBlocks_u(blocks_y, blocks_z);

    reset_flag_u<<<numBlocks_u, dimBlocks_u>>>(_flag_u[dev], _dom[dev], bc);

    // reset flag_v
    if(dom[dev].Gfy.kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gfy.knb;
    else
      threads_z = MAX_THREADS_DIM;

    if(dom[dev].Gfy.in < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfy.inb;
    else
      threads_x = MAX_THREADS_DIM;

    blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);
    blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);

    dim3 dimBlocks_v(threads_z, threads_x);
    dim3 numBlocks_v(blocks_z, blocks_x);

    reset_flag_v<<<numBlocks_v, dimBlocks_v>>>(_flag_v[dev], _dom[dev], bc);

    // reset flag_w
    if(dom[dev].Gfz.in < MAX_THREADS_DIM)
      threads_x = dom[dev].Gfz.inb;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gfz.jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gfz.jnb;
    else
      threads_y = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);

    dim3 dimBlocks_w(threads_x, threads_y);
    dim3 numBlocks_w(blocks_x, blocks_y);

    reset_flag_w<<<numBlocks_w, dimBlocks_w>>>(_flag_w[dev], _dom[dev], bc);


 // push point_particle information to device
    checkCudaErrors(cudaMemcpy(_points[dev], points, sizeof(point_struct) * npoints,
      cudaMemcpyHostToDevice));



  
  
      // fill in ghost cells for periodic boundary conditions
      if(bc.uW == PERIODIC) {
        blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);
        blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);
        numBlocks_u.x = blocks_y;
        numBlocks_u.y = blocks_z;
        cage_flag_u_periodic_W<<<numBlocks_u, dimBlocks_u>>>(_flag_u[dev],
          _dom[dev]);
      }
      if(bc.uE == PERIODIC) {
        blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);
        blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);
        numBlocks_u.x = blocks_y;
        numBlocks_u.y = blocks_z;
        cage_flag_u_periodic_E<<<numBlocks_u, dimBlocks_u>>>(_flag_u[dev],
          _dom[dev]);
      }
      if(bc.uS == PERIODIC) {
        blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);
        blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);
        numBlocks_u.x = blocks_z;
        numBlocks_u.y = blocks_x;
        cage_flag_u_periodic_S<<<numBlocks_u, dimBlocks_u>>>(_flag_u[dev],
          _dom[dev]);
      }
      if(bc.uN == PERIODIC) {
        blocks_z = (int)ceil((real) dom[dev].Gfx.knb / (real) threads_z);
        blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);
        numBlocks_u.x = blocks_z;
        numBlocks_u.y = blocks_x;
        cage_flag_u_periodic_N<<<numBlocks_u, dimBlocks_u>>>(_flag_u[dev],
          _dom[dev]);
      }
      if(bc.uB == PERIODIC) {
        blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);
        blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);
        numBlocks_u.x = blocks_x;
        numBlocks_u.y = blocks_y;
        cage_flag_u_periodic_B<<<numBlocks_u, dimBlocks_u>>>(_flag_u[dev],
          _dom[dev]);
      }
      if(bc.uT == PERIODIC) {
        blocks_x = (int)ceil((real) dom[dev].Gfx.inb / (real) threads_x);
        blocks_y = (int)ceil((real) dom[dev].Gfx.jnb / (real) threads_y);
        numBlocks_u.x = blocks_x;
        numBlocks_u.y = blocks_y;
        cage_flag_u_periodic_T<<<numBlocks_u, dimBlocks_u>>>(_flag_u[dev],
          _dom[dev]);
      }



      if(bc.vW == PERIODIC) {
        blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);
        blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);
        numBlocks_v.x = blocks_y;
        numBlocks_v.y = blocks_z;
        cage_flag_v_periodic_W<<<numBlocks_v, dimBlocks_v>>>(_flag_v[dev],
          _dom[dev]);
      }
      if(bc.vE == PERIODIC) {
        blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);
        blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);
        numBlocks_v.x = blocks_y;
        numBlocks_v.y = blocks_z;
        cage_flag_v_periodic_E<<<numBlocks_v, dimBlocks_v>>>(_flag_v[dev],
          _dom[dev]);
      }
      if(bc.vS == PERIODIC) {
        blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);
        blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);
        numBlocks_v.x = blocks_z;
        numBlocks_v.y = blocks_x;
        cage_flag_v_periodic_S<<<numBlocks_v, dimBlocks_v>>>(_flag_v[dev],
          _dom[dev]);
      }
      if(bc.vN == PERIODIC) {
        blocks_z = (int)ceil((real) dom[dev].Gfy.knb / (real) threads_z);
        blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);
        numBlocks_v.x = blocks_z;
        numBlocks_v.y = blocks_x;
        cage_flag_v_periodic_N<<<numBlocks_v, dimBlocks_v>>>(_flag_v[dev],
          _dom[dev]);
      }
      if(bc.vB == PERIODIC) {
        blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);
        blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);
        numBlocks_v.x = blocks_x;
        numBlocks_v.y = blocks_y;
        cage_flag_v_periodic_B<<<numBlocks_v, dimBlocks_v>>>(_flag_v[dev],
          _dom[dev]);
      }
      if(bc.vT == PERIODIC) {
        blocks_x = (int)ceil((real) dom[dev].Gfy.inb / (real) threads_x);
        blocks_y = (int)ceil((real) dom[dev].Gfy.jnb / (real) threads_y);
        numBlocks_v.x = blocks_x;
        numBlocks_v.y = blocks_y;
        cage_flag_v_periodic_T<<<numBlocks_v, dimBlocks_v>>>(_flag_v[dev],
          _dom[dev]);
      }




      if(bc.wW == PERIODIC) {
        blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);
        blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);
        numBlocks_w.x = blocks_y;
        numBlocks_w.y = blocks_z;
        cage_flag_w_periodic_W<<<numBlocks_w, dimBlocks_w>>>(_flag_w[dev],
          _dom[dev]);
      }
      if(bc.wE == PERIODIC) {
        blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);
        blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);
        numBlocks_w.x = blocks_y;
        numBlocks_w.y = blocks_z;
        cage_flag_w_periodic_E<<<numBlocks_w, dimBlocks_w>>>(_flag_w[dev],
          _dom[dev]);
      }
      if(bc.wS == PERIODIC) {
        blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);
        blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);
        numBlocks_w.x = blocks_z;
        numBlocks_w.y = blocks_x;
        cage_flag_w_periodic_S<<<numBlocks_w, dimBlocks_w>>>(_flag_w[dev],
          _dom[dev]);
      }
      if(bc.wN == PERIODIC) {
        blocks_z = (int)ceil((real) dom[dev].Gfz.knb / (real) threads_z);
        blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);
        numBlocks_w.x = blocks_z;
        numBlocks_w.y = blocks_x;
        cage_flag_w_periodic_N<<<numBlocks_w, dimBlocks_w>>>(_flag_w[dev],
          _dom[dev]);
      }
      if(bc.wB == PERIODIC) {
        blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);
        blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);
        numBlocks_w.x = blocks_x;
        numBlocks_w.y = blocks_y;
        cage_flag_w_periodic_B<<<numBlocks_w, dimBlocks_w>>>(_flag_w[dev],
          _dom[dev]);
      }
      if(bc.wT == PERIODIC) {
        blocks_x = (int)ceil((real) dom[dev].Gfz.inb / (real) threads_x);
        blocks_y = (int)ceil((real) dom[dev].Gfz.jnb / (real) threads_y);
        numBlocks_w.x = blocks_x;
        numBlocks_w.y = blocks_y;
        cage_flag_w_periodic_T<<<numBlocks_w, dimBlocks_w>>>(_flag_w[dev],
          _dom[dev]);
                           }
   }
  
}


