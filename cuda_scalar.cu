#include "cuda_bicgstab.h"
#include "cuda_bluebottle.h"
#include "cuda_point.h"
#include "cuda_scalar.h"
#include "entrySearch.h"

#include <cuda.h>
#include <helper_cuda.h>

#include <cusp/array1d.h>
#include <cusp/blas.h>
#include <cusp/dia_matrix.h>
#include <cusp/monitor.h>
#include <cusp/precond/diagonal.h>
#include <cusp/krylov/bicgstab.h>
#include <cusp/krylov/cg.h>
#include <cusp/print.h>



extern "C"
void cuda_diffScalar_sub_explicitH(int coordiSys,int dev, real *scSrc)
{
//Find the characteristic length scale for the diffusion equation
if(DIFF_dt<=0) return;
/*
printf("\ndiff s3 %d %d \n",s3,dom[dev].Gcc.s3);
fflush(stdout);
*/
int index=coordiSys*48 +2;
int in=DomInfo[index];
int jn=DomInfo[index+6];
int kn=DomInfo[index+12];

 index=coordiSys*48 +23;
int s3b=DomInfo[index];

// write coefficients using kernel
int threads_x = 0;
int threads_y = 0;
int threads_z = 0;
int blocks_x = 0;
int blocks_y = 0;
int blocks_z = 0;
if(in < MAX_THREADS_DIM)
threads_x = in;
else
threads_x = MAX_THREADS_DIM;
if(jn < MAX_THREADS_DIM)
threads_y = jn;
else
threads_y = MAX_THREADS_DIM;
if(kn < MAX_THREADS_DIM)
threads_z = kn;
else
threads_z = MAX_THREADS_DIM;
blocks_x = (int)ceil((real) in / (real) threads_x);
blocks_y = (int)ceil((real) jn / (real) threads_y);
blocks_z = (int)ceil((real) kn / (real) threads_z);
dim3 dimBlocks_x(threads_y, threads_z);
dim3 numBlocks_x(blocks_y, blocks_z);
dim3 dimBlocks_y(threads_z, threads_x);
dim3 numBlocks_y(blocks_z, blocks_x);
dim3 dimBlocks_z(threads_x, threads_y);
dim3 numBlocks_z(blocks_x, blocks_y);
dim3 dimBlocks,numBlocks;
switch(coordiSys)
{
case 0:
dimBlocks.x=dimBlocks_x.x;
dimBlocks.y=dimBlocks_x.y;
numBlocks.x=numBlocks_x.x;
numBlocks.y=numBlocks_x.y;
break;
case 1:
dimBlocks.x=dimBlocks_x.x;
dimBlocks.y=dimBlocks_x.y;
numBlocks.x=numBlocks_x.x;
numBlocks.y=numBlocks_x.y;
break;
case 2:
dimBlocks.x=dimBlocks_y.x;
dimBlocks.y=dimBlocks_y.y;
numBlocks.x=numBlocks_y.x;
numBlocks.y=numBlocks_y.y;
break;
case 3:
dimBlocks.x=dimBlocks_z.x;
dimBlocks.y=dimBlocks_z.y;
numBlocks.x=numBlocks_z.x;
numBlocks.y=numBlocks_z.y;
break;
default:break;
}
real dx=Dom.dx;
real dy=Dom.dx;
real dz=Dom.dz;
real ddx=1/dx/dx;
real ddy=1/dy/dy;
real ddz=1/dz/dz;
real DIFF_dt_done=0;
//real DIFF_dt_sub=dt_sc;
real DIFF_dt_diffScalar=min(0.5*CFL/(ddx+ddy+ddz),DIFF_dt);

//create buffer array
    real *scSrc_buf;
    checkCudaErrors(cudaMalloc((void**) &scSrc_buf, sizeof(real) * s3b));

//real niter=10.f;
//if( DIFF_dt_diffScalar*niter>DIFF_dt)  DIFF_dt_diffScalar=DIFF_dt/niter;

//int iter=0;
real DIFF_dt_sub;
while(DIFF_dt_done<DIFF_dt)
{
//DIFF_dt_sub=min(dt_sc,DIFF_dt-DIFF_dt_done);
DIFF_dt_sub=min(DIFF_dt_diffScalar,DIFF_dt-DIFF_dt_done);
//DIFF_dt_sub=DIFF_dt;
DIFF_dt_done +=DIFF_dt_sub;
/*
iter +=1;
printf("\niter %d %f %f %f %f\n",iter,DIFF_dt,DIFF_dt_done,DIFF_dt_diffScalar,DIFF_dt_sub);
fflush(stdout);
*/
BC_sc_W_P<<<numBlocks, dimBlocks>>>(scSrc, _dom[dev]);
BC_sc_E_P<<<numBlocks, dimBlocks>>>(scSrc, _dom[dev]);
BC_sc_T_P<<<numBlocks, dimBlocks>>>(scSrc, _dom[dev]);
BC_sc_S_P<<<numBlocks, dimBlocks>>>(scSrc, _dom[dev]);
BC_sc_B_P<<<numBlocks, dimBlocks>>>(scSrc, _dom[dev]);
BC_sc_N_P<<<numBlocks, dimBlocks>>>(scSrc, _dom[dev]);

checkCudaErrors(cudaMemcpy(scSrc_buf, scSrc, dom[dev].Gcc.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
diffScalar_explicitD<<<numBlocks, dimBlocks>>>(scSrc,scSrc_buf,_dom[dev],DIFF_dt_sub);
/*
BC_sc_W_P<<<numBlocks, dimBlocks>>>(scSrc, _dom[dev]);
BC_sc_E_P<<<numBlocks, dimBlocks>>>(scSrc, _dom[dev]);
BC_sc_T_P<<<numBlocks, dimBlocks>>>(scSrc, _dom[dev]);
BC_sc_S_P<<<numBlocks, dimBlocks>>>(scSrc, _dom[dev]);
BC_sc_B_P<<<numBlocks, dimBlocks>>>(scSrc, _dom[dev]);
BC_sc_N_P<<<numBlocks, dimBlocks>>>(scSrc, _dom[dev]);
*/
}

    checkCudaErrors(cudaFree(scSrc_buf));
}

extern "C"
void cuda_diffScalar_helmholtz_CN(int coordiSys,int dev, real *scSrc)
{

//Find the characteristic length scale for the diffusion equation
if(DIFF_dt<=0) return;

  // CPU thread for multi-GPU
int index=coordiSys*48        +20;
int s3=DomInfo[index];
int s2=DomInfo[index-1];
int s1=DomInfo[index-2];

int s3b=DomInfo[index+3];

index=coordiSys*48        +2;
int in=DomInfo[index];
int jn=DomInfo[index+6];
int kn=DomInfo[index+12];

real dx2=Dom.dx*Dom.dx;
real dy2=Dom.dy*Dom.dy;
real dz2=Dom.dz*Dom.dz;
real ddx=1/dx2;
real ddy=1/dy2;
real ddz=1/dz2;
 

//create buffer array
    real *scSrc_buf;
    checkCudaErrors(cudaMalloc((void**) &scSrc_buf, sizeof(real) * s3b));

    // create CUSP diagonal matrix for p
    cusp::dia_matrix<int, real, cusp::device_memory> *_A_diffScalar;
    _A_diffScalar = new cusp::dia_matrix<int, real, cusp::device_memory>
      (s3,s3, 0, 13);

    // create temporary diffScalar without ghost cells
    real *_diffScalar_noghost;
    checkCudaErrors(cudaMalloc((void**) &_diffScalar_noghost,
      sizeof(real) * s3));

//When theta=1/2 the scheme becomes Crank-Nicolson; theta=1 correspond to fully implicit diffusion term;theta=1 correspond to fully explicit diffusion term;
//The choice for this theta is to maitain monoticity while have higher temperal accuracy than fully implicit or explicit scheme!
//real theta=max(1-0.5f/DIFF_dt/(ddx+ddy+ddz),0.f);
real theta=0.5f;;
//real theta=1.f;;
//real thetaD=DIFF_dt_sub*theta; thetaR=DIFF_dt_sub*(1-theta);
real thetaD,thetaR;
real DIFF_dt_done=0;
real DIFF_dt_diffScalar=min(0.5f/(ddx+ddy+ddz)/(1-theta),DIFF_dt);
real DIFF_dt_sub;
int iter=0;
while(DIFF_dt_done<DIFF_dt)
{
DIFF_dt_sub=min(DIFF_dt_diffScalar,DIFF_dt-DIFF_dt_done);
DIFF_dt_done +=DIFF_dt_sub;
thetaD=DIFF_dt_sub*theta;
thetaR=DIFF_dt_sub*(1-theta);
iter +=1;
printf("\niter_implicit %d %f %f %f\n",iter,DIFF_dt,DIFF_dt_done,DIFF_dt_diffScalar);
fflush(stdout);
  // write right-hand side
    cuda_diffScalar_rhs_CN(dev,scSrc_buf,scSrc,thetaR);
    
 
    checkCudaErrors(cudaMemcpy(scSrc,scSrc_buf, s3b*sizeof(real), cudaMemcpyDeviceToDevice));

    // create temporary U* without ghost cells for bicgstab result
    cusp::array1d<real, cusp::device_memory> diffScalar_tmp(s3, 0.);


    // set up the coefficient matrix for Crank-Nicolson
    _A_diffScalar->diagonal_offsets[0]  = -s3 + s2;
    _A_diffScalar->diagonal_offsets[1]  = -s2;
    _A_diffScalar->diagonal_offsets[2]  = -s2 + s1;
    _A_diffScalar->diagonal_offsets[3]  = -s1;

if(coordiSys==0)
{
    _A_diffScalar->diagonal_offsets[4]  = -s1 + 1;//has to be 2 if use Gfx for staggered grid, 1 for Gcc
    _A_diffScalar->diagonal_offsets[8]  = s1 - 1;//has to be 2 if use Gfx for staggered grid,1 for Gcc
}
else
{
    _A_diffScalar->diagonal_offsets[4]  = -s1 + 2;//has to be 2 if use Gfx for staggered grid, 1 for Gcc
    _A_diffScalar->diagonal_offsets[8]  = s1 - 2;//has to be 2 if use Gfx for staggered grid,1 for Gcc
}
    _A_diffScalar->diagonal_offsets[5]  = -1;
    _A_diffScalar->diagonal_offsets[6]  = 0;
    _A_diffScalar->diagonal_offsets[7]  = 1;
    _A_diffScalar->diagonal_offsets[9]  = s1;
    _A_diffScalar->diagonal_offsets[10] = s2 - s1;
    _A_diffScalar->diagonal_offsets[11] = s2;
    _A_diffScalar->diagonal_offsets[12] = s3 - s2;

    // write coefficients using kernel
    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    if(in < MAX_THREADS_DIM)
      threads_x = in;
    else
      threads_x = MAX_THREADS_DIM;

    if(jn < MAX_THREADS_DIM)
      threads_y = jn;
    else
      threads_y = MAX_THREADS_DIM;

    if(kn < MAX_THREADS_DIM)
      threads_z = kn;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) in / (real) threads_x);
    blocks_y = (int)ceil((real) jn / (real) threads_y);
    blocks_z = (int)ceil((real) kn / (real) threads_z);

    dim3 dimBlocks_x(threads_y, threads_z);
    dim3 numBlocks_x(blocks_y, blocks_z);
    dim3 dimBlocks_y(threads_z, threads_x);
    dim3 numBlocks_y(blocks_z, blocks_x);
    dim3 dimBlocks_z(threads_x, threads_y);
    dim3 numBlocks_z(blocks_x, blocks_y);

    dim3 dimBlocks,numBlocks;

switch(coordiSys)
{
case 0:
dimBlocks.x=dimBlocks_x.x;
dimBlocks.y=dimBlocks_x.y;
numBlocks.x=numBlocks_x.x;
numBlocks.y=numBlocks_x.y;
break;
case 1:
dimBlocks.x=dimBlocks_x.x;
dimBlocks.y=dimBlocks_x.y;
numBlocks.x=numBlocks_x.x;
numBlocks.y=numBlocks_x.y;
break;
case 2:
dimBlocks.x=dimBlocks_y.x;
dimBlocks.y=dimBlocks_y.y;
numBlocks.x=numBlocks_y.x;
numBlocks.y=numBlocks_y.y;
break;
case 3:
dimBlocks.x=dimBlocks_z.x;
dimBlocks.y=dimBlocks_z.y;
numBlocks.x=numBlocks_z.x;
numBlocks.y=numBlocks_z.y;
break;
default:break;
}

    // copy u_star into noghost structure for Helmholtz right-hand side
    copy_diffSc_noghost<<<numBlocks, dimBlocks>>>(_diffScalar_noghost, scSrc,
      _dom[dev],coordiSys);



    // build pressure-Poisson coefficient matrix
    diffScalar_coeffs_init<<<numBlocks, dimBlocks>>>(_dom[dev],
      _A_diffScalar->values.pitch,
      thrust::raw_pointer_cast(&_A_diffScalar->values.values[0]),coordiSys);
  
    diffScalar_coeffs_CN<<<numBlocks, dimBlocks>>>(thetaD,_dom[dev],
      _A_diffScalar->values.pitch,
      thrust::raw_pointer_cast(&_A_diffScalar->values.values[0]),
      _flag_u[dev],_flag_v[dev],_flag_w[dev],coordiSys);

int bcW, bcE, bcS, bcN, bcB, bcT;  // cell locations


switch(coordiSys)
{
case 0:
bcW=sc_bc.scW;
bcE=sc_bc.scE;
bcS=sc_bc.scS;
bcN=sc_bc.scN;
bcB=sc_bc.scB;
bcT=sc_bc.scT;
break;
case 1:
bcW=bc.uW;
bcE=bc.uE;
bcS=bc.uS;
bcN=bc.uN;
bcB=bc.uB;
bcT=bc.uT;
break;
case 2:
bcW=bc.vW;
bcE=bc.vE;
bcS=bc.vS;
bcN=bc.vN;
bcB=bc.vB;
bcT=bc.vT;
break;
case 3:
bcW=bc.wW;
bcE=bc.wE;
bcS=bc.wS;
bcN=bc.wN;
bcB=bc.wB;
bcT=bc.wT;
break;
default:break; 
}


    // account for boundary conditions
    if(bcW == PERIODIC)
      diffScalar_coeffs_periodic_W<<<numBlocks_x, dimBlocks_x>>>(thetaD,_dom[dev],
        _A_diffScalar->values.pitch,
        thrust::raw_pointer_cast(&_A_diffScalar->values.values[0]),coordiSys);
    if(bcE == PERIODIC)
      diffScalar_coeffs_periodic_E<<<numBlocks_x, dimBlocks_x>>>(thetaD,_dom[dev],
        _A_diffScalar->values.pitch,
        thrust::raw_pointer_cast(&_A_diffScalar->values.values[0]),coordiSys);
    if(bcS == PERIODIC)
      diffScalar_coeffs_periodic_S<<<numBlocks_y, dimBlocks_y>>>(thetaD,_dom[dev],
        _A_diffScalar->values.pitch,
        thrust::raw_pointer_cast(&_A_diffScalar->values.values[0]),coordiSys);
    if(bcN == PERIODIC)
      diffScalar_coeffs_periodic_N<<<numBlocks_y, dimBlocks_y>>>(thetaD,_dom[dev],
        _A_diffScalar->values.pitch,
        thrust::raw_pointer_cast(&_A_diffScalar->values.values[0]),coordiSys);
    if(bcB == PERIODIC)
      diffScalar_coeffs_periodic_B<<<numBlocks_z, dimBlocks_z>>>(thetaD,_dom[dev],
        _A_diffScalar->values.pitch,
        thrust::raw_pointer_cast(&_A_diffScalar->values.values[0]),coordiSys);
    if(bcT == PERIODIC)
      diffScalar_coeffs_periodic_T<<<numBlocks_z, dimBlocks_z>>>(thetaD,_dom[dev],
        _A_diffScalar->values.pitch,
        thrust::raw_pointer_cast(&_A_diffScalar->values.values[0]),coordiSys);

  
    // create CUSP pointer to right-hand side
    thrust::device_ptr<real> _ptr_diffScalar(_diffScalar_noghost);
    cusp::array1d<real, cusp::device_memory> *_diffScalar_rhs;
    _diffScalar_rhs = new cusp::array1d<real, cusp::device_memory>(_ptr_diffScalar,
      _ptr_diffScalar + s3);

 
    // normalize the problem by the right-hand side before sending to CUSP
    real norm = cusp::blas::nrm2(*_diffScalar_rhs);
    if(norm == 0)       norm = 1.;
    cusp::blas::scal(*_diffScalar_rhs, 1. / norm);

    // call BiCGSTAB to solve for diffScalar_tmp
    cusp::convergence_monitor<real> monitor(*_diffScalar_rhs, pp_max_iter,
      pp_residual);
    cusp::precond::diagonal<real, cusp::device_memory> M(*_A_diffScalar);
    //cusp::krylov::bicgstab(*_A_diffScalar, diffScalar_tmp, *_diffScalar_rhs, monitor, M);
    cusp::krylov::cg(*_A_diffScalar, diffScalar_tmp, *_diffScalar_rhs, monitor, M);
    // write convergence data to file
      recorder_bicgstab("solver_helmholtz_diffScalar.rec", monitor.iteration_count(),monitor.residual_norm());
    if(!monitor.converged()) {
      printf("The diffScalar Helmholtz equation did not converge.              \n");
      exit(EXIT_FAILURE);
    }

    // unnormalize the solution
    cusp::blas::scal(diffScalar_tmp, norm);

    // copy solution back to scSrc
    copy_sc_ghost<<<numBlocks, dimBlocks>>>(scSrc,
      thrust::raw_pointer_cast(diffScalar_tmp.data()), _dom[dev]);

    delete(_diffScalar_rhs);
}
    // clean up
    checkCudaErrors(cudaFree(_diffScalar_noghost));
    checkCudaErrors(cudaFree(scSrc_buf));

}





extern "C"
void cuda_diffScalar_rhs_CN(int dev,real *scSrc,real *scSrc0,real thetaR)
{
  int threads_y = 0;
  int threads_z = 0;
  int blocks_y = 0;
  int blocks_z = 0;

  if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
    threads_y = dom[dev].Gcc.jnb + 2;
  else
    threads_y = MAX_THREADS_DIM;

  if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
    threads_z = dom[dev].Gcc.knb + 2;
  else
    threads_z = MAX_THREADS_DIM;

  blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) (threads_y-2));
  blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) (threads_z-2));

  dim3 dimBlocks(threads_y, threads_z);
  dim3 numBlocks(blocks_y, blocks_z);
 // diffScalar_rhs_CN<<<numBlocks, dimBlocks>>>(DIFF_dt,scSrc0 , scSrc,_dom[dev],theta);
  diffScalar_rhs_CN<<<numBlocks, dimBlocks>>>(thetaR,scSrc0 , scSrc,_dom[dev]);


//fflush(stdout);
getLastCudaError("Kernel execution failed.");
}

 

extern "C"
void cuda_scalar_malloc(void)
{
/*
  _omega_x = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _omega_y = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _omega_z = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
*/

//malloc device scalar on host
  _sc = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _sc0 = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _diff0_sc = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _diff_sc = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _conv0_sc = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _conv_sc = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _scSrc = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _scSrc0 = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);


  _epsp = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _epsp0 = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);



  // allocate device memory on device
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

//printf("\n dom[dev].Gfx.s3b,dom[dev].Gfy.s3b,dom[dev].Gfz.s3b %d %d %d \n",dom[dev].Gfx.s3b,dom[dev].Gfy.s3b,dom[dev].Gfz.s3b);
//printf("\ndom s3 %d %d %d %d \n",dom[dev].Gcc.s3,dom[dev].Gfx.s3,dom[dev].Gfy.s3,dom[dev].Gfz.s3);
fflush(stdout);

//allocate scalar on device
    checkCudaErrors(cudaMalloc((void**) &(_sc[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);
    checkCudaErrors(cudaMalloc((void**) &(_sc0[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);
    checkCudaErrors(cudaMalloc((void**) &(_diff0_sc[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);
    checkCudaErrors(cudaMalloc((void**) &(_diff_sc[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);
    checkCudaErrors(cudaMalloc((void**) &(_conv0_sc[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);
    checkCudaErrors(cudaMalloc((void**) &(_conv_sc[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);
    checkCudaErrors(cudaMalloc((void**) &(_scSrc[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);
    checkCudaErrors(cudaMalloc((void**) &(_scSrc0[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);


    checkCudaErrors(cudaMalloc((void**) &(_epsp[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);
    checkCudaErrors(cudaMalloc((void**) &(_epsp0[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);


    // TODO add CUSP solver data structures to memory usage count

    //printf("Device %d of %d using %f Mb global memory.\n", dev, nsubdom, mb);
  }



}


extern "C"
void cuda_scalar_free(void)
{
  // free device memory on device
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

    checkCudaErrors(cudaFree(_sc0[dev]));
    checkCudaErrors(cudaFree(_sc[dev]));
    checkCudaErrors(cudaFree(_diff0_sc[dev]));
    checkCudaErrors(cudaFree(_diff_sc[dev]));
    checkCudaErrors(cudaFree(_conv0_sc[dev]));
    checkCudaErrors(cudaFree(_conv_sc[dev]));
    checkCudaErrors(cudaFree(_scSrc[dev]));
    checkCudaErrors(cudaFree(_scSrc0[dev]));
    checkCudaErrors(cudaFree(_epsp[dev]));
    checkCudaErrors(cudaFree(_epsp0[dev]));


  }

  // free device memory on host

    free(_sc0);
    free(_sc);
    free(_diff0_sc);
    free(_diff_sc);
    free(_conv0_sc);
    free(_conv_sc);
    free(_scSrc);
    free(_scSrc0);
    free(_epsp);
    free(_epsp0);


}

extern "C"
void cuda_scalar_push(void)
{
  // copy host data to device
  #pragma omp parallel num_threads(nsubdom)
  {
    int i, j, k;          // iterators
    int ii, jj, kk;       // helper iterators
    int C, CC;            // cell references

    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

    // set up host working arrays for subdomain copy from host to device

/*
//add by shigan
    real *stress_uu = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    real *stress_vv = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    real *stress_ww = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));


    real *omega_xx = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    real *omega_yy = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    real *omega_zz = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));
*/
    real *scc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    real *scc0 = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    real *diff0_scc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    real *diff_scc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    real *conv0_scc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    real *conv_scc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    real *scSrcc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));

    real *epspc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    // select appropriate subdomain
    // add by shigan
/*
//stress_uu
    for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
      for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
        for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
          ii = i - dom[dev].Gfx.isb;
          jj = j - dom[dev].Gfx.jsb;
          kk = k - dom[dev].Gfx.ksb;
          C = i + j * Dom.Gfx.s1b + k * Dom.Gfx.s2b;
          CC = ii + jj * dom[dev].Gfx.s1b + kk * dom[dev].Gfx.s2b;
          stress_uu[CC] = 0;
        }
      }
    }

  //stress_vv
    for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
      for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
        for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
          ii = i - dom[dev].Gfy.isb;
          jj = j - dom[dev].Gfy.jsb;
          kk = k - dom[dev].Gfy.ksb;
          C = i + j * Dom.Gfy.s1b + k * Dom.Gfy.s2b;
          CC = ii + jj * dom[dev].Gfy.s1b + kk * dom[dev].Gfy.s2b;
          stress_vv[CC] =0;
        }
      }
    }

  // stress_ww
    for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
      for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
        for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
          ii = i - dom[dev].Gfz.isb;
          jj = j - dom[dev].Gfz.jsb;
          kk = k - dom[dev].Gfz.ksb;
          C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
          CC = ii + jj * dom[dev].Gfz.s1b + kk * dom[dev].Gfz.s2b;
          stress_ww[CC] =0;
        }
      }
    }

*/

/*
//omega_x
    for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
      for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
        for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
          ii = i - dom[dev].Gfx.isb;
          jj = j - dom[dev].Gfx.jsb;
          kk = k - dom[dev].Gfx.ksb;
          C = i + j * Dom.Gfx.s1b + k * Dom.Gfx.s2b;
          CC = ii + jj * dom[dev].Gfx.s1b + kk * dom[dev].Gfx.s2b;
          omega_xx[CC] = 0;
        }
      }
    }

  //omega_y
    for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
      for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
        for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
          ii = i - dom[dev].Gfy.isb;
          jj = j - dom[dev].Gfy.jsb;
          kk = k - dom[dev].Gfy.ksb;
          C = i + j * Dom.Gfy.s1b + k * Dom.Gfy.s2b;
          CC = ii + jj * dom[dev].Gfy.s1b + kk * dom[dev].Gfy.s2b;
          omega_yy[CC] =0;
        }
      }
    }

  // omega_z
    for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
      for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
        for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
          ii = i - dom[dev].Gfz.isb;
          jj = j - dom[dev].Gfz.jsb;
          kk = k - dom[dev].Gfz.ksb;
          C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
          CC = ii + jj * dom[dev].Gfz.s1b + kk * dom[dev].Gfz.s2b;
          omega_zz[CC] =0;
        }
      }
    }
*/

//scalar initialization
 for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
      for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
        for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
          ii = i - dom[dev].Gcc.isb;
          jj = j - dom[dev].Gcc.jsb;
          kk = k - dom[dev].Gcc.ksb;
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
          scc[CC] = sc[C];
          scc0[CC] = sc0[C];
          diff0_scc[CC] = diff0_sc[C];
          diff_scc[CC] = diff_sc[C];
          conv0_scc[CC] = conv0_sc[C];
          conv_scc[CC] = conv_sc[C];
	  scSrcc[CC] = scSrc[C];
	  epspc[CC] = epsp[C];
        }
      }
    }


/*
    // copy from host to device
    checkCudaErrors(cudaMemcpy(_stress_u[dev],stress_uu, sizeof(real) * dom[dev].Gfx.s3b,
      cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(_stress_v[dev],stress_vv, sizeof(real) * dom[dev].Gfy.s3b,
      cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(_stress_w[dev],stress_ww, sizeof(real) * dom[dev].Gfz.s3b,
      cudaMemcpyHostToDevice));


    checkCudaErrors(cudaMemcpy(_omega_x[dev],omega_xx, sizeof(real) * dom[dev].Gfx.s3b,
      cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(_omega_y[dev],omega_yy, sizeof(real) * dom[dev].Gfy.s3b,
      cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(_omega_z[dev],omega_zz, sizeof(real) * dom[dev].Gfz.s3b,
      cudaMemcpyHostToDevice));
*/
    checkCudaErrors(cudaMemcpy(_sc[dev],scc, sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(_sc0[dev],scc0, sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(_diff0_sc[dev],diff0_scc, sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(_diff_sc[dev],diff_scc, sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(_conv0_sc[dev],conv0_scc, sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(_conv_sc[dev],conv_scc, sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(_scSrc[dev],scSrcc, sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(_epsp[dev],epspc, sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyHostToDevice));

    // free host subdomain working arrays
/*
//add by shigan
    free(stress_uu);
    free(stress_vv);
    free(stress_ww);

    free(omega_xx);
    free(omega_yy);
    free(omega_zz);
*/
    free(scc);
    free(scc0);
    free(diff0_scc);
    free(diff_scc);  
    free(conv0_scc);
    free(conv_scc);
    free(scSrcc);

    free(epspc);

  }
}

//copy scalar from device to host
extern "C"
void cuda_scalar_pull(void)
{
  // copy device data to host
  #pragma omp parallel num_threads(nsubdom)
  {
    int i, j, k;          // iterators
    int ii, jj, kk;       // helper iterators
    int C, CC;            // cell references

    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

 real *scc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    cpumem += dom[dev].Gcc.s3b * sizeof(real);
 real *scc0 = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    cpumem += dom[dev].Gcc.s3b * sizeof(real);
 real *diff0_scc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    cpumem += dom[dev].Gcc.s3b * sizeof(real);
 real *diff_scc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    cpumem += dom[dev].Gcc.s3b * sizeof(real);
 real *conv0_scc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
    cpumem += dom[dev].Gcc.s3b * sizeof(real);
 real *conv_scc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
   cpumem += dom[dev].Gcc.s3b * sizeof(real);
 real *scSrcc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
   cpumem += dom[dev].Gcc.s3b * sizeof(real);
  
 real *epspc = (real*) malloc(dom[dev].Gcc.s3b * sizeof(real));
   cpumem += dom[dev].Gcc.s3b * sizeof(real);
//add by shigan
/*
    real *stress_uu = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    real *stress_vv = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    real *stress_ww = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));

    real *omega_xx = (real*) malloc(dom[dev].Gfx.s3b * sizeof(real));
    real *omega_yy = (real*) malloc(dom[dev].Gfy.s3b * sizeof(real));
    real *omega_zz = (real*) malloc(dom[dev].Gfz.s3b * sizeof(real));

 checkCudaErrors(cudaMemcpy(stress_uu, _stress_u[dev], sizeof(real) * dom[dev].Gfx.s3b,
      cudaMemcpyDeviceToHost)); 
 checkCudaErrors(cudaMemcpy(stress_vv, _stress_v[dev], sizeof(real) * dom[dev].Gfy.s3b,
      cudaMemcpyDeviceToHost)); 
 checkCudaErrors(cudaMemcpy(stress_ww, _stress_w[dev], sizeof(real) * dom[dev].Gfz.s3b,
      cudaMemcpyDeviceToHost)); 
  
 checkCudaErrors(cudaMemcpy(omega_xx, _omega_x[dev], sizeof(real) * dom[dev].Gfx.s3b,
      cudaMemcpyDeviceToHost)); 
 checkCudaErrors(cudaMemcpy(omega_yy, _omega_y[dev], sizeof(real) * dom[dev].Gfy.s3b,
      cudaMemcpyDeviceToHost)); 
 checkCudaErrors(cudaMemcpy(omega_zz, _omega_z[dev], sizeof(real) * dom[dev].Gfz.s3b,
      cudaMemcpyDeviceToHost)); 
 
*/
   // copy from device to host
     checkCudaErrors(cudaMemcpy(scc, _sc[dev], sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(scc0,_sc0[dev], sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(diff0_scc, _diff0_sc[dev], sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(diff_scc,_diff_sc[dev], sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(conv0_scc,_conv0_sc[dev], sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(conv_scc,_conv_sc[dev], sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaMemcpy(scSrcc,_scSrc[dev], sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyDeviceToHost));

    checkCudaErrors(cudaMemcpy(epspc,_epsp[dev], sizeof(real) * dom[dev].Gcc.s3b,
      cudaMemcpyDeviceToHost));


   
 // run simulation
    // fill in apropriate subdomain
    // scalar
    for(k = dom[dev].Gcc.ksb; k < dom[dev].Gcc.keb; k++) {
      for(j = dom[dev].Gcc.jsb; j < dom[dev].Gcc.jeb; j++) {
        for(i = dom[dev].Gcc.isb; i < dom[dev].Gcc.ieb; i++) {
          ii = i - dom[dev].Gcc.isb;
          jj = j - dom[dev].Gcc.jsb;
          kk = k - dom[dev].Gcc.ksb;
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          CC = ii + jj * dom[dev].Gcc.s1b + kk * dom[dev].Gcc.s2b;
          sc0[C] = scc0[CC];
          sc[C] = scc[CC];
          diff0_sc[C] = diff0_scc[CC];
          diff_sc[C] = diff_scc[CC];
          conv0_sc[C] = conv0_scc[CC];
          conv_sc[C] = conv_scc[CC];

          scSrc[C] = scSrcc[CC];
          epsp[C] = epspc[CC];
        }
      }
    }

/*
  //stress_uu
    for(k = dom[dev].Gfx.ksb; k < dom[dev].Gfx.keb; k++) {
      for(j = dom[dev].Gfx.jsb; j < dom[dev].Gfx.jeb; j++) {
        for(i = dom[dev].Gfx.isb; i < dom[dev].Gfx.ieb; i++) {
          ii = i - dom[dev].Gfx.isb;
          jj = j - dom[dev].Gfx.jsb;
          kk = k - dom[dev].Gfx.ksb;
          C = i + j * Dom.Gfx.s1b + k * Dom.Gfx.s2b;
          CC = ii + jj * dom[dev].Gfx.s1b + kk * dom[dev].Gfx.s2b;
          stress_u[CC]=stress_uu[CC] ;
	  omega_x[CC]=omega_xx[CC]

        }
      }
    }

  //stress_vv
    for(k = dom[dev].Gfy.ksb; k < dom[dev].Gfy.keb; k++) {
      for(j = dom[dev].Gfy.jsb; j < dom[dev].Gfy.jeb; j++) {
        for(i = dom[dev].Gfy.isb; i < dom[dev].Gfy.ieb; i++) {
          ii = i - dom[dev].Gfy.isb;
          jj = j - dom[dev].Gfy.jsb;
          kk = k - dom[dev].Gfy.ksb;
          C = i + j * Dom.Gfy.s1b + k * Dom.Gfy.s2b;
          CC = ii + jj * dom[dev].Gfy.s1b + kk * dom[dev].Gfy.s2b;
          stress_v[CC]=stress_vv[CC];
          omega_y[CC]=omega_yy[CC]
        }
      }
    }

  // stress_ww
    for(k = dom[dev].Gfz.ksb; k < dom[dev].Gfz.keb; k++) {
      for(j = dom[dev].Gfz.jsb; j < dom[dev].Gfz.jeb; j++) {
        for(i = dom[dev].Gfz.isb; i < dom[dev].Gfz.ieb; i++) {
          ii = i - dom[dev].Gfz.isb;
          jj = j - dom[dev].Gfz.jsb;
          kk = k - dom[dev].Gfz.ksb;
          C = i + j * Dom.Gfz.s1b + k * Dom.Gfz.s2b;
          CC = ii + jj * dom[dev].Gfz.s1b + kk * dom[dev].Gfz.s2b;
          stress_w[CC]=stress_ww[CC];
	  omega_z[CC]=omega_zz[CC]
        }
      }
    }

*/



    // free host subdomain working arrays
   
    free(scc);
    free(scc0);
    free(diff0_scc);
    free(diff_scc);
    free(conv0_scc);
    free(conv_scc);

    free(scSrcc);
    free(epspc);
/*
    free(stress_uu);
    free(stress_vv);
    free(stress_ww);
    free(omega_xx);
    free(omega_yy);
    free(omega_zz);
*/
  }
}





extern "C"
void cuda_store_scalar(void)
{
  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

    checkCudaErrors(cudaMemcpy(_conv0_sc[dev], _conv_sc[dev],
      dom[dev].Gcc.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
    checkCudaErrors(cudaMemcpy(_diff0_sc[dev], _diff_sc[dev],
      dom[dev].Gcc.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
    checkCudaErrors(cudaMemcpy(_sc0[dev], _sc[dev], 
      dom[dev].Gcc.s3b*sizeof(real), cudaMemcpyDeviceToDevice));

    checkCudaErrors(cudaMemcpy(_scSrc0[dev], _scSrc[dev],
      dom[dev].Gcc.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
    checkCudaErrors(cudaMemcpy(_epsp0[dev], _epsp[dev],
      dom[dev].Gcc.s3b*sizeof(real), cudaMemcpyDeviceToDevice));
  }
}

extern "C"
void cuda_scalar_advance(void)
{

// CPU threading for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

    int threads_y = 0;
    int threads_z = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    if(dom[dev].Gcc._jnb < MAX_THREADS_DIM)
      threads_y = dom[dev].Gcc._jnb;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gcc._knb < MAX_THREADS_DIM)
      threads_z = dom[dev].Gcc._knb;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_y = (int)ceil((real) dom[dev].Gcc._jnb / (real) (threads_y-2));
    blocks_z = (int)ceil((real) dom[dev].Gcc._knb / (real) (threads_z-2));

    dim3 dimBlocks_x(threads_y, threads_z);
    dim3 numBlocks_x(blocks_y, blocks_z);

 

// Add the point particle source to scalar equation
    int threads = MAX_THREADS_1D;
    int blocks = (int)ceil((real) npoints / (real) threads);
 //   int blocks_st = blocks*STENCIL3;

    dim3 dimBlocks(threads);
    dim3 numBlocks(blocks);
 //   dim3 numBlocks_st(blocks_st);
    dim3 dimBlocks_3d,numBlocks_3d;

int coordiSys=0;
int incGhost=1;
    block_thread_cell_3D(dimBlocks_3d,numBlocks_3d,dom[dev],coordiSys,incGhost);

//Locate the point particle in each grid cell, store the grid cell number in points.i~points.k
 // if(npoints>0)  lpt_localize<<<numBlocks, dimBlocks>>>(npoints,_points[dev], _dom[dev],bc);

//initialize flow field array to 0 on device, the array length is Nx*Ny*Nz
//include scalar source and particle volume fraction divided by cell volume
////    lpt_scalar_source_init<<<numBlocks_x, dimBlocks_x>>>(_scSrc[dev],_epsp[dev], _dom[dev]);
    lpt_scalar_source_init<<<numBlocks_3d, dimBlocks_3d>>>(_scSrc[dev],_epsp[dev], _dom[dev]);

//lpt_scalar_source_init_test<<<numBlocks_x, dimBlocks_x>>>(_scSrc[dev], _dom[dev],ttime_done,DIFF_eq);
//lpt_scalar_source_convDiff_test<<<numBlocks_x, dimBlocks_x>>>(_scSrc[dev], _dom[dev],ttime_done,DIFF_eq);


//Mollify volume fraction on device, don't need too much thread source
lpt_mollify_sc_optH(coordiSys,EPSP_TYPE,dev,_epsp[dev]);

//lpt_mollify_delta_scH(coordiSys,EPSP_TYPE,dev,_epsp[dev]);
getLastCudaError("Kernel execution failed.");

/*
    int lenSrc=dom[dev].Gcc.s3b;
    dim3 dimBlocks_s,numBlocks_s;
    block_thread_point(dimBlocks_s,numBlocks_s,lenSrc);
    print_kernel_array_real<<<numBlocks_s,dimBlocks_s>>>(_epsp[dev],lenSrc);
 */


lpt_epsp_clip<<<numBlocks_x, dimBlocks_x>>>(_epsp[dev],_dom[dev]);



//Mollify source of scalar on device
if(sc_twoway>0) 
  {   
lpt_mollify_sc_optH(coordiSys,SCALAR_TYPE,dev,_scSrc[dev]);

//lpt_mollify_delta_scH(coordiSys,SCALAR_TYPE,dev,_scSrc[dev]);
}
fflush(stdout);

/*
cudaEvent_t start, stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
float milliseconds = 0;
cudaEventRecord(start);
*/

/*
//advance scalar TODO add boundary condition to sc in the kernel!,  takes 2.6 ms compared to 5 ms by u_star_2
if(dt0 > 0.) {
advance_sc_upwind_1st<<<numBlocks_x, dimBlocks_x>>>(DIFF_eq, _u[dev], _v[dev], _w[dev], _scSrc[dev],_epsp[dev],  _diff0_sc[dev], _conv0_sc[dev], _diff_sc[dev], _conv_sc[dev], _sc[dev], _sc0[dev],_dom[dev],dt0_try,dt_try);
fflush(stdout);

//advance_sc<<<numBlocks_x, dimBlocks_x>>>(DIFF_eq, _u[dev], _v[dev], _w[dev], _scSrc[dev],_epsp[dev],  _diff0_sc[dev], _conv0_sc[dev], _diff_sc[dev], _conv_sc[dev], _sc[dev], _sc0[dev],_dom[dev],dt0_try,dt_try);
}
else
{
advance_sc_upwind_1st_init<<<numBlocks_x, dimBlocks_x>>>(DIFF_eq, _u[dev], _v[dev], _w[dev], _scSrc[dev],_epsp[dev], _diff0_sc[dev], _conv0_sc[dev], _diff_sc[dev], _conv_sc[dev], _sc[dev], _sc0[dev],_dom[dev],dt0_try,dt_try);

//advance_sc_init<<<numBlocks_x, dimBlocks_x>>>(DIFF_eq, _u[dev], _v[dev], _w[dev], _scSrc[dev],_epsp[dev], _diff0_sc[dev], _conv0_sc[dev], _diff_sc[dev], _conv_sc[dev], _sc[dev], _sc0[dev],_dom[dev],dt0_try,dt_try);
}
*/

/*
//Using MacCormack scheme to advance scalar
advance_sc_macCormack<<<numBlocks_x, dimBlocks_x>>>(DIFF_eq, _u[dev], _v[dev], _w[dev], _scSrc[dev],_epsp[dev], _diff_sc[dev], _conv_sc[dev], _sc[dev], _sc0[dev],_dom[dev],dt_try);
fflush(stdout);
*/

//advance_sc_QUICK<<<numBlocks_x, dimBlocks_x>>>(DIFF_eq, _u[dev], _v[dev], _w[dev], _scSrc[dev],_epsp[dev], _diff_sc[dev], _conv_sc[dev], _sc[dev], _sc0[dev],_dom[dev],dt_try);

 cuda_scalar_helmholtz();

getLastCudaError("Kernel execution failed.");



//boundary condition of scalar

 /*
cudaEventRecord(stop);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&milliseconds, start, stop);
printf("\ntime_sc %f\n",milliseconds);
fflush(stdout);
*/
 }
}

 

extern "C"
real cuda_find_dt_implicit(void)
{
  // results from all devices
  real *dts = (real*) malloc(nsubdom * sizeof(real));
    // cpumem += nsubdom * sizeof(real);

real max_dt=1.f;
  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));
    // search
    real u_max = find_max_mag(dom[dev].Gfx.s3, _u[dev]);
    real v_max = find_max_mag(dom[dev].Gfy.s3, _v[dev]);
    real w_max = find_max_mag(dom[dev].Gfz.s3, _w[dev]);

    dts[dev]= EPSILON/max_dt;
    dts[dev] +=u_max / dom[dev].dx;
    dts[dev] +=v_max / dom[dev].dy;
    dts[dev] +=w_max / dom[dev].dz;
    dts[dev] = CFL / dts[dev];

 //   dts[dev] = dts[dev]/100.f;
  }
  // find max of all devices
  real dt_buf = -1.;
  for(int i = 0; i < nsubdom; i++)
    if(dts[i] > dt_buf) dt_buf = dts[i];

  // clean up
  free(dts);

  dt_buf=min(dt_buf,max_dt);

  real alpha=0.7f;
  if(dt_buf>dt0&&dt0>0) dt_buf=alpha*dt_buf+(1-alpha)*dt0;

  return dt_buf;
}





extern "C"
real cuda_find_dt_sc(real dt)
{
  // results from all devices
  real *dts = (real*) malloc(nsubdom * sizeof(real));
    // cpumem += nsubdom * sizeof(real);

  // parallelize over CPU threads
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));


    // search
    real u_max = find_max_mag(dom[dev].Gfx.s3, _u[dev]);
    real v_max = find_max_mag(dom[dev].Gfy.s3, _v[dev]);
    real w_max = find_max_mag(dom[dev].Gfz.s3, _w[dev]);
/*
//FTCS scheme with Adam-Bashforth method
    dts[dev] = (u_max + 2 * DIFF_eq / dom[dev].dx) / dom[dev].dx  + u_max*u_max/2/DIFF_eq;
    dts[dev] += (v_max + 2 * DIFF_eq / dom[dev].dy) / dom[dev].dy + v_max*v_max/2/DIFF_eq;
    dts[dev] += (w_max + 2 * DIFF_eq / dom[dev].dz) / dom[dev].dz + w_max*w_max/2/DIFF_eq;
    dts[dev] = CFL / dts[dev];
*/


/*
//MacCormack scheme ||QUICK scheme
    dts[dev] =  (u_max + 2 * DIFF_eq / dom[dev].dx) / dom[dev].dx;
    dts[dev] += (v_max + 2 * DIFF_eq / dom[dev].dy) / dom[dev].dy;
    dts[dev] += (w_max + 2 * DIFF_eq / dom[dev].dz) / dom[dev].dz;
    dts[dev] = CFL / dts[dev];
*/

//1st upwind scheme
    dts[dev] =  (u_max + 2 * DIFF_eq / dom[dev].dx) / dom[dev].dx;
    dts[dev] += (v_max + 2 * DIFF_eq / dom[dev].dy) / dom[dev].dy;
    dts[dev] += (w_max + 2 * DIFF_eq / dom[dev].dz) / dom[dev].dz;
    dts[dev] = CFL / dts[dev];


  }

  // find max of all devices
  real max = -1.;
  for(int i = 0; i < nsubdom; i++)
    if(dts[i] > max) max = dts[i];

  // clean up
  free(dts);
//Make the 1st time step smaller
 // if(dt0<=0) dt=0.01f*dt;
//  if(dt>1e-5) dt=1e-5;
 
  if(max>dt) max=dt; 

  return max;
}


 
 

extern "C"
void cuda_scalar_helmholtz(void)
{
  // CPU thread for multi-GPU
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    cudaSetDevice(dev + dev_start);
    // write right-hand side
    cuda_scalar_rhs(dev);

    // create temporary U* without ghost cells for bicgstab result
    cusp::array1d<real, cusp::device_memory> scalar_tmp(dom[dev].Gcc.s3, 0.);

    // create CUSP diagonal matrix for p
    cusp::dia_matrix<int, real, cusp::device_memory> *_A_scalar;
    _A_scalar = new cusp::dia_matrix<int, real, cusp::device_memory>
      (dom[dev].Gcc._s3, dom[dev].Gcc._s3, 0, 13);

    // set up the coefficient matrix for Crank-Nicolson
    _A_scalar->diagonal_offsets[0]  = -dom[dev].Gcc._s3 + dom[dev].Gcc._s2;
    _A_scalar->diagonal_offsets[1]  = -dom[dev].Gcc._s2;
    _A_scalar->diagonal_offsets[2]  = -dom[dev].Gcc._s2 + dom[dev].Gcc._s1;
    _A_scalar->diagonal_offsets[3]  = -dom[dev].Gcc._s1;
    _A_scalar->diagonal_offsets[4]  = -dom[dev].Gcc._s1 + 1;//has to be 2 if use Gfx for staggered grid, 1 for Gcc
    _A_scalar->diagonal_offsets[5]  = -1;
    _A_scalar->diagonal_offsets[6]  = 0;
    _A_scalar->diagonal_offsets[7]  = 1;
    _A_scalar->diagonal_offsets[8]  = dom[dev].Gcc._s1 - 1;//has to be 2 if use Gfx for staggered grid,1 for Gcc
    _A_scalar->diagonal_offsets[9]  = dom[dev].Gcc._s1;
    _A_scalar->diagonal_offsets[10] = dom[dev].Gcc._s2 - dom[dev].Gcc._s1;
    _A_scalar->diagonal_offsets[11] = dom[dev].Gcc._s2;
    _A_scalar->diagonal_offsets[12] = dom[dev].Gcc._s3 - dom[dev].Gcc._s2;

    // write coefficients using kernel
    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    if(dom[dev].Gcc._in < MAX_THREADS_DIM)
      threads_x = dom[dev].Gcc._in;
    else
      threads_x = MAX_THREADS_DIM;

    if(dom[dev].Gcc._jn < MAX_THREADS_DIM)
      threads_y = dom[dev].Gcc._jn;
    else
      threads_y = MAX_THREADS_DIM;

    if(dom[dev].Gcc._kn < MAX_THREADS_DIM)
      threads_z = dom[dev].Gcc._kn;
    else
      threads_z = MAX_THREADS_DIM;

    blocks_x = (int)ceil((real) dom[dev].Gcc._in / (real) threads_x);
    blocks_y = (int)ceil((real) dom[dev].Gcc._jn / (real) threads_y);
    blocks_z = (int)ceil((real) dom[dev].Gcc._kn / (real) threads_z);

    dim3 dimBlocks_x(threads_y, threads_z);
    dim3 numBlocks_x(blocks_y, blocks_z);
    dim3 dimBlocks_y(threads_z, threads_x);
    dim3 numBlocks_y(blocks_z, blocks_x);
    dim3 dimBlocks_z(threads_x, threads_y);
    dim3 numBlocks_z(blocks_x, blocks_y);

    // create temporary scalar without ghost cells
    real *_scalar_noghost;
    checkCudaErrors(cudaMalloc((void**) &_scalar_noghost,
      sizeof(real) * dom[dev].Gcc.s3));
    // copy u_star into noghost structure for Helmholtz right-hand side
    copy_sc_noghost<<<numBlocks_x, dimBlocks_x>>>(_scalar_noghost, _sc[dev],
      _dom[dev]);

    // build pressure-Poisson coefficient matrix
    scalar_coeffs_init<<<numBlocks_x, dimBlocks_x>>>(_dom[dev],
      _A_scalar->values.pitch,
      thrust::raw_pointer_cast(&_A_scalar->values.values[0]));
    /*
    scalar_coeffs<<<numBlocks_x, dimBlocks_x>>>(DIFF_eq, dt_try, _dom[dev],
      _A_scalar->values.pitch,
      thrust::raw_pointer_cast(&_A_scalar->values.values[0]));
*/
    scalar_coeffs<<<numBlocks_x, dimBlocks_x>>>(DIFF_eq, dt_try, _dom[dev],
      _A_scalar->values.pitch,
      thrust::raw_pointer_cast(&_A_scalar->values.values[0]),
      _flag_u[dev],_flag_v[dev],_flag_w[dev],_epsp[dev]);

/*
    // account for boundary conditions
    if(sc_bc.scW == PERIODIC)
      scalar_coeffs_periodic_W<<<numBlocks_x, dimBlocks_x>>>(DIFF_eq, dt_try, _dom[dev],
        _A_scalar->values.pitch,
        thrust::raw_pointer_cast(&_A_scalar->values.values[0]));
    if(sc_bc.scE == PERIODIC)
      scalar_coeffs_periodic_E<<<numBlocks_x, dimBlocks_x>>>(DIFF_eq, dt_try, _dom[dev],
        _A_scalar->values.pitch,
        thrust::raw_pointer_cast(&_A_scalar->values.values[0]));
    if(sc_bc.scS == PERIODIC)
      scalar_coeffs_periodic_S<<<numBlocks_y, dimBlocks_y>>>(DIFF_eq, dt_try, _dom[dev],
        _A_scalar->values.pitch,
        thrust::raw_pointer_cast(&_A_scalar->values.values[0]));
    if(sc_bc.scN == PERIODIC)
      scalar_coeffs_periodic_N<<<numBlocks_y, dimBlocks_y>>>(DIFF_eq, dt_try, _dom[dev],
        _A_scalar->values.pitch,
        thrust::raw_pointer_cast(&_A_scalar->values.values[0]));
    if(sc_bc.scB == PERIODIC)
      scalar_coeffs_periodic_B<<<numBlocks_z, dimBlocks_z>>>(DIFF_eq, dt_try, _dom[dev],
        _A_scalar->values.pitch,
        thrust::raw_pointer_cast(&_A_scalar->values.values[0]));
    if(sc_bc.scT == PERIODIC)
      scalar_coeffs_periodic_T<<<numBlocks_z, dimBlocks_z>>>(DIFF_eq, dt_try, _dom[dev],
        _A_scalar->values.pitch,
        thrust::raw_pointer_cast(&_A_scalar->values.values[0]));
*/
   // account for boundary conditions
    if(sc_bc.scW == PERIODIC)
      scalar_coeffs_periodic_W<<<numBlocks_x, dimBlocks_x>>>(DIFF_eq, dt_try, _dom[dev],
        _A_scalar->values.pitch,
        thrust::raw_pointer_cast(&_A_scalar->values.values[0]), _epsp[dev],_flag_u[dev]);
    if(sc_bc.scE == PERIODIC)
      scalar_coeffs_periodic_E<<<numBlocks_x, dimBlocks_x>>>(DIFF_eq, dt_try, _dom[dev],
        _A_scalar->values.pitch,
        thrust::raw_pointer_cast(&_A_scalar->values.values[0]), _epsp[dev],_flag_u[dev]);
    if(sc_bc.scS == PERIODIC)
      scalar_coeffs_periodic_S<<<numBlocks_y, dimBlocks_y>>>(DIFF_eq, dt_try, _dom[dev],
        _A_scalar->values.pitch,
        thrust::raw_pointer_cast(&_A_scalar->values.values[0]), _epsp[dev],_flag_v[dev]);
    if(sc_bc.scN == PERIODIC)
      scalar_coeffs_periodic_N<<<numBlocks_y, dimBlocks_y>>>(DIFF_eq, dt_try, _dom[dev],
        _A_scalar->values.pitch,
        thrust::raw_pointer_cast(&_A_scalar->values.values[0]), _epsp[dev],_flag_v[dev]);
    if(sc_bc.scB == PERIODIC)
      scalar_coeffs_periodic_B<<<numBlocks_z, dimBlocks_z>>>(DIFF_eq, dt_try, _dom[dev],
        _A_scalar->values.pitch,
        thrust::raw_pointer_cast(&_A_scalar->values.values[0]), _epsp[dev],_flag_w[dev]);
    if(sc_bc.scT == PERIODIC)
      scalar_coeffs_periodic_T<<<numBlocks_z, dimBlocks_z>>>(DIFF_eq, dt_try, _dom[dev],
        _A_scalar->values.pitch,
        thrust::raw_pointer_cast(&_A_scalar->values.values[0]), _epsp[dev],_flag_w[dev]);
    // create CUSP pointer to right-hand side
    thrust::device_ptr<real> _ptr_scalar(_scalar_noghost);
    cusp::array1d<real, cusp::device_memory> *_scalar_rhs;
    _scalar_rhs = new cusp::array1d<real, cusp::device_memory>(_ptr_scalar,
      _ptr_scalar + dom[dev].Gcc._s3);

getLastCudaError("Kernel execution failed.");

    // normalize the problem by the right-hand side before sending to CUSP
    real norm = cusp::blas::nrm2(*_scalar_rhs);
    if(norm == 0)       norm = 1.;
    cusp::blas::scal(*_scalar_rhs, 1. / norm);

    // call BiCGSTAB to solve for scalar_tmp
    cusp::convergence_monitor<real> monitor(*_scalar_rhs, pp_max_iter,
      pp_residual);
    cusp::precond::diagonal<real, cusp::device_memory> M(*_A_scalar);
    //cusp::krylov::bicgstab(*_A_scalar, scalar_tmp, *_scalar_rhs, monitor, M);
    cusp::krylov::cg(*_A_scalar, scalar_tmp, *_scalar_rhs, monitor, M);
    // write convergence data to file
      recorder_bicgstab("solver_helmholtz_scalar.rec", monitor.iteration_count(),monitor.residual_norm());
    if(!monitor.converged()) {
      printf("The scalar Helmholtz equation did not converge.              \n");
      exit(EXIT_FAILURE);
    }

    // unnormalize the solution
    cusp::blas::scal(scalar_tmp, norm);

    // copy solution back to _sc
    copy_sc_ghost<<<numBlocks_x, dimBlocks_x>>>(_sc[dev],
      thrust::raw_pointer_cast(scalar_tmp.data()), _dom[dev]);

    // clean up
    delete(_scalar_rhs);
    delete(_A_scalar);
    checkCudaErrors(cudaFree(_scalar_noghost));

  }
}





extern "C"
void cuda_scalar_rhs(int dev)
{
  int threads_y = 0;
  int threads_z = 0;
  int blocks_y = 0;
  int blocks_z = 0;

  if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
    threads_y = dom[dev].Gcc.jnb + 2;
  else
    threads_y = MAX_THREADS_DIM;

  if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
    threads_z = dom[dev].Gcc.knb + 2;
  else
    threads_z = MAX_THREADS_DIM;

  blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) (threads_y-2));
  blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) (threads_z-2));

  dim3 dimBlocks(threads_y, threads_z);
  dim3 numBlocks(blocks_y, blocks_z);
/*
  scalar_rhs_FTCS<<<numBlocks, dimBlocks>>>(rho_f, DIFF_eq,
			 _u[dev], _v[dev], _w[dev],
			 _epsp[dev], _scSrc[dev],
			 _conv0_sc[dev],_conv_sc[dev],_diff_sc[dev],
			_sc0[dev],_sc[dev], 
			_dom[dev], dt_try, dt0_try);
*/
  scalar_rhs_upwind_1st<<<numBlocks, dimBlocks>>>(rho_f, DIFF_eq,
			 _u[dev], _v[dev], _w[dev],
			 _epsp[dev],_epsp0[dev],
			 _scSrc[dev],_scSrc0[dev],
			 _conv0_sc[dev],_conv_sc[dev],_diff_sc[dev],
			_sc0[dev],_sc[dev], 
			_dom[dev], dt_try, dt0_try);

//fflush(stdout);
getLastCudaError("Kernel execution failed.");
}




extern "C"
void cuda_scSrc_BC(int coordiSys,int valType, real *scSrc, int dev)
{
  // CPU threading for multi-GPU
    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

int   index=coordiSys*48        +5;
int inb=DomInfo[index];
int jnb=DomInfo[index+6];
int knb=DomInfo[index+12];


int scSrcW,scSrcE,scSrcS,scSrcN,scSrcT,scSrcB;
real scSrcWD,scSrcED,scSrcSD,scSrcND,scSrcTD,scSrcBD;
switch(coordiSys)
{
case 0:
	switch(valType)
	{
	case PRESSURE_TYPE:
	scSrcW=bc.pW;
	scSrcE=bc.pE;
	scSrcS=bc.pS;
	scSrcN=bc.pN;
	scSrcT=bc.pT;
	scSrcB=bc.pB;
	break;
	case SCALAR_TYPE:
	scSrcW=sc_bc.scW;
	scSrcE=sc_bc.scE;
	scSrcS=sc_bc.scS;
	scSrcN=sc_bc.scN;
	scSrcT=sc_bc.scT;
	scSrcB=sc_bc.scB;

	scSrcWD=sc_bc.scWD;	
	scSrcED=sc_bc.scED;
	scSrcTD=sc_bc.scTD;
	scSrcBD=sc_bc.scBD;
	scSrcSD=sc_bc.scSD;
	scSrcND=sc_bc.scND;

	break;
	default:break;
	}
break;
case 1:
	scSrcW=bc.uW;
	scSrcE=bc.uE;
	scSrcS=bc.uS;
	scSrcN=bc.uN;
	scSrcT=bc.uT;
	scSrcB=bc.uB;

	scSrcWD=bc.uWD;	
	scSrcED=bc.uED;
	scSrcTD=bc.uTD;
	scSrcBD=bc.uBD;
	scSrcSD=bc.uSD;
	scSrcND=bc.uND;

	break;
case 2:
	scSrcW=bc.vW;
	scSrcE=bc.vE;
	scSrcS=bc.vS;
	scSrcN=bc.vN;
	scSrcT=bc.vT;
	scSrcB=bc.vB;

	scSrcWD=bc.vWD;	
	scSrcED=bc.vED;
	scSrcTD=bc.vTD;
	scSrcBD=bc.vBD;
	scSrcSD=bc.vSD;
	scSrcND=bc.vND;

	break;
case 3:
	scSrcW=bc.wW;
	scSrcE=bc.wE;
	scSrcS=bc.wS;
	scSrcN=bc.wN;
	scSrcT=bc.wT;
	scSrcB=bc.wB;

	scSrcWD=bc.wWD;	
	scSrcED=bc.wED;
	scSrcTD=bc.wTD;
	scSrcBD=bc.wBD;
	scSrcSD=bc.wSD;
	scSrcND=bc.wND;

	break;
default: break;
}




    // check whether each subdomain boundary (E, W, N, S, T, B) is
    // an external boundary
    if(dom[dev].W == -1) {
      // set up kernel call
      // pressure
      if(jnb < MAX_THREADS_DIM)
        threads_y = jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(knb < MAX_THREADS_DIM)
        threads_z = knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) jnb / (real) threads_y);
      blocks_z = (int)ceil((real) knb / (real) threads_z);

      dim3 dimBlocks_p(threads_y, threads_z);
      dim3 numBlocks_p(blocks_y, blocks_z);

      // apply BC to all fields for this face
      switch(scSrcW) {
        case PERIODIC:
          BC_scSrc_W_P<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev]);
          break;
        case NEUMANN:
          BC_scSrc_W_N<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev]);
          break;
         case DIRICHLET:
          BC_scSrc_W_D<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev],scSrcWD);
          break;
  
    }
 }
   if(dom[dev].E == -1) {
      // set up kernel call
      // pressure
      if(jnb < MAX_THREADS_DIM)
        threads_y = jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(knb < MAX_THREADS_DIM)
        threads_z = knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) jnb / (real) threads_y);
      blocks_z = (int)ceil((real) knb / (real) threads_z);

      dim3 dimBlocks_p(threads_y, threads_z);
      dim3 numBlocks_p(blocks_y, blocks_z);

      // apply BC to all fields for this face
       switch(scSrcE) {
        case PERIODIC:
          BC_scSrc_E_P<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev]);
          break;
        case NEUMANN:
          BC_scSrc_E_N<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev]);
          break;
        case DIRICHLET:
          BC_scSrc_E_D<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev],scSrcED);
          break;

    }

}
    if(dom[dev].S == -1) {
      // set up kernel call
      // pressure
      if(knb < MAX_THREADS_DIM)
        threads_z = knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(inb < MAX_THREADS_DIM)
        threads_x = inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) knb / (real) threads_z);
      blocks_x = (int)ceil((real) inb / (real) threads_x);

      dim3 dimBlocks_p(threads_z, threads_x);
      dim3 numBlocks_p(blocks_z, blocks_x);

     
 // apply BC to all fields for this face
          switch(scSrcS) {
        case PERIODIC:
          BC_scSrc_S_P<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev]);
          break;
        case NEUMANN:
          BC_scSrc_S_N<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev]);
          break;
        case DIRICHLET:
          BC_scSrc_S_D<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev],scSrcSD);
          break;
    }
}
    if(dom[dev].N == -1) {
      // set up kernel call
      // pressure
      if(knb < MAX_THREADS_DIM)
        threads_z = knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(inb < MAX_THREADS_DIM)
        threads_x = inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) knb / (real) threads_z);
      blocks_x = (int)ceil((real) inb / (real) threads_x);

      dim3 dimBlocks_p(threads_z, threads_x);
      dim3 numBlocks_p(blocks_z, blocks_x);



      // apply BC to all fields for this face
      switch(scSrcN) {
        case PERIODIC:
          BC_scSrc_N_P<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev]);
          break;
        case NEUMANN:
          BC_scSrc_N_N<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev]);
          break;
        case DIRICHLET:
          BC_scSrc_N_D<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev],scSrcND);
          break;

    }
}
    if(dom[dev].B == -1) {
      // set up kernel call
      // pressure
      if(inb < MAX_THREADS_DIM)
        threads_x = inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(jnb < MAX_THREADS_DIM)
        threads_y = jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) inb / (real) threads_x);
      blocks_y = (int)ceil((real) jnb / (real) threads_y);

      dim3 dimBlocks_p(threads_x, threads_y);
      dim3 numBlocks_p(blocks_x, blocks_y);



      // apply BC to all fields for this face
           switch(scSrcB) {
        case PERIODIC:
          BC_scSrc_B_P<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev]);
          break;
        case NEUMANN:
          BC_scSrc_B_N<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev]);
          break;
        case DIRICHLET:
          BC_scSrc_B_D<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev],scSrcBD);
          break;

    }
}
    if(dom[dev].T == -1) {
      // set up kernel call
      // pressure
      if(inb < MAX_THREADS_DIM)
        threads_x = inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(jnb < MAX_THREADS_DIM)
        threads_y = jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) inb / (real) threads_x);
      blocks_y = (int)ceil((real) jnb / (real) threads_y);

      dim3 dimBlocks_p(threads_x, threads_y);
      dim3 numBlocks_p(blocks_x, blocks_y);


      // apply BC to all fields for this face
            switch(scSrcT) {
        case PERIODIC:
          BC_scSrc_T_P<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev]);
          break;
        case NEUMANN:
          BC_scSrc_T_N<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev]);
          break;
        case DIRICHLET:
          BC_scSrc_T_D<<<numBlocks_p, dimBlocks_p>>>(coordiSys, scSrc, _dom[dev],scSrcTD);
          break;

   	 	}

    	}
getLastCudaError("Kernel execution failed.");

}

extern "C"
void cuda_scalar_BC(void)
{
  // CPU threading for multi-GPU
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

    // check whether each subdomain boundary (E, W, N, S, T, B) is
    // an external boundary
    if(dom[dev].W == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);

      dim3 dimBlocks_p(threads_y, threads_z);
      dim3 numBlocks_p(blocks_y, blocks_z);

      // apply BC to all fields for this face
      switch(sc_bc.scW) {
        case PERIODIC:
          BC_sc_W_P<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_sc_W_N<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev]);
          break;
         case DIRICHLET:
          BC_sc_W_D<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev],sc_bc.scWD);
          break;
  
    }
 }
   if(dom[dev].E == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);
      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);

      dim3 dimBlocks_p(threads_y, threads_z);
      dim3 numBlocks_p(blocks_y, blocks_z);

      // apply BC to all fields for this face
       switch(sc_bc.scE) {
        case PERIODIC:
          BC_sc_E_P<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_sc_E_N<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_sc_E_D<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev],sc_bc.scED);
          break;

    }

}
    if(dom[dev].S == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);

      dim3 dimBlocks_p(threads_z, threads_x);
      dim3 numBlocks_p(blocks_z, blocks_x);

     
 // apply BC to all fields for this face
          switch(sc_bc.scS) {
        case PERIODIC:
          BC_sc_S_P<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_sc_S_N<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_sc_S_D<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev],sc_bc.scSD);
          break;

    }
}
    if(dom[dev].N == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.knb < MAX_THREADS_DIM)
        threads_z = dom[dev].Gcc.knb;
      else
        threads_z = MAX_THREADS_DIM;

      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      blocks_z = (int)ceil((real) dom[dev].Gcc.knb / (real) threads_z);
      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);

      dim3 dimBlocks_p(threads_z, threads_x);
      dim3 numBlocks_p(blocks_z, blocks_x);



      // apply BC to all fields for this face
      switch(sc_bc.scN) {
        case PERIODIC:
          BC_sc_N_P<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_sc_N_N<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_sc_N_D<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev],sc_bc.scND);
          break;

    }
}
    if(dom[dev].B == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);

      dim3 dimBlocks_p(threads_x, threads_y);
      dim3 numBlocks_p(blocks_x, blocks_y);



      // apply BC to all fields for this face
           switch(sc_bc.scB) {
        case PERIODIC:
          BC_sc_B_P<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_sc_B_N<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_sc_B_D<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev],sc_bc.scBD);
          break;

    }
}
    if(dom[dev].T == -1) {
      // set up kernel call
      // pressure
      if(dom[dev].Gcc.inb < MAX_THREADS_DIM)
        threads_x = dom[dev].Gcc.inb;
      else
        threads_x = MAX_THREADS_DIM;

      if(dom[dev].Gcc.jnb < MAX_THREADS_DIM)
        threads_y = dom[dev].Gcc.jnb;
      else
        threads_y = MAX_THREADS_DIM;

      blocks_x = (int)ceil((real) dom[dev].Gcc.inb / (real) threads_x);
      blocks_y = (int)ceil((real) dom[dev].Gcc.jnb / (real) threads_y);

      dim3 dimBlocks_p(threads_x, threads_y);
      dim3 numBlocks_p(blocks_x, blocks_y);


      // apply BC to all fields for this face
            switch(sc_bc.scT) {
        case PERIODIC:
          BC_sc_T_P<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev]);
          break;
        case NEUMANN:
          BC_sc_T_N<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev]);
          break;
        case DIRICHLET:
          BC_sc_T_D<<<numBlocks_p, dimBlocks_p>>>(_sc[dev], _dom[dev],sc_bc.scTD);
          break;

   	 	}

    	}
getLastCudaError("Kernel execution failed.");

 }
}

