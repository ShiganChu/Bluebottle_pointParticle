#include "cuda_point.h"

#include <cuda.h>
#include <helper_cuda.h>

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

//printf("test %lf",points[0].rho);
    checkCudaErrors(cudaMemcpy(_points[dev], points, sizeof(point_struct) * npoints,cudaMemcpyHostToDevice));
  //  checkCudaErrors(cudaMemcpy(points,_points[dev],sizeof(point_struct) * npoints,cudaMemcpyDeviceToHost));
//printf("test %lf %lf",_points[dev][0].rho,points[0].rho);
//printf("test %lf",points[0].rho);
int test=1;
 }
}

extern "C"
void cuda_point_pull(void)
{
  // all devices have the same point_particle data for now, so just copy one of them
  checkCudaErrors(cudaMemcpy(points, _points[0], sizeof(point_struct) * npoints,
    cudaMemcpyDeviceToHost));

}


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


