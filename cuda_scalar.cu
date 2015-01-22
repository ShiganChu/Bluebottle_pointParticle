#include "cuda_bicgstab.h"
#include "cuda_bluebottle.h"
#include "cuda_point.h"
#include "cuda_scalar.h"
#include "entrySearch.h"

#include <cuda.h>
#include <helper_cuda.h>

extern "C"
void cuda_scalar_malloc(void)
{
  // allocate device memory on host
  //add by shigan_9_22_2014, fluid stress on face center
  _stress_u = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _stress_v = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);
  _stress_w = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);

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

  _epsp = (real**) malloc(nsubdom * sizeof(real*));
  cpumem += nsubdom * sizeof(real*);



  // allocate device memory on device
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

printf("\n dom[dev].Gfx.s3b,dom[dev].Gfy.s3b,dom[dev].Gfz.s3b %d %d %d \n",dom[dev].Gfx.s3b,dom[dev].Gfy.s3b,dom[dev].Gfz.s3b);
//add by shigan
    checkCudaErrors(cudaMalloc((void**) &(_stress_u[dev]),
      sizeof(real) * dom[dev].Gfx.s3b));
    checkCudaErrors(cudaMalloc((void**) &(_stress_v[dev]),
      sizeof(real) * dom[dev].Gfy.s3b));
    checkCudaErrors(cudaMalloc((void**) &(_stress_w[dev]),
      sizeof(real) * dom[dev].Gfz.s3b));

/*
    checkCudaErrors(cudaMalloc((void**) &(_omega_x[dev]),
      sizeof(real) * dom[dev].Gfx.s3b));
    checkCudaErrors(cudaMalloc((void**) &(_omega_y[dev]),
      sizeof(real) * dom[dev].Gfy.s3b));
    checkCudaErrors(cudaMalloc((void**) &(_omega_z[dev]),
      sizeof(real) * dom[dev].Gfz.s3b));
*/

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

    checkCudaErrors(cudaMalloc((void**) &(_epsp[dev]),
      sizeof(real) * dom[dev].Gcc.s3b));
    gpumem += dom[dev].Gcc.s3b * sizeof(real);


    // TODO add CUSP solver data structures to memory usage count

    //printf("Device %d of %d using %f Mb global memory.\n", dev, nsubdom, mb);
  }



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


#ifdef DEBUG // run test code
   
#else // run simulation
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


#endif

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
    checkCudaErrors(cudaFree(_epsp[dev]));


    checkCudaErrors(cudaFree(_stress_u[dev]));
    checkCudaErrors(cudaFree(_stress_v[dev]));
    checkCudaErrors(cudaFree(_stress_w[dev]));
/*
    checkCudaErrors(cudaFree(_omega_x[dev]));
    checkCudaErrors(cudaFree(_omega_y[dev]));
    checkCudaErrors(cudaFree(_omega_z[dev]));
*/
  }

  // free device memory on host
  free(_stress_u);
  free(_stress_v);
  free(_stress_w);
/*  
  free(_omega_x);
  free(_omega_y);
  free(_omega_z);
*/
    free(_sc0);
    free(_sc);
    free(_diff0_sc);
    free(_diff_sc);
    free(_conv0_sc);
    free(_conv_sc);
    free(_scSrc);
    free(_epsp);


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
    block_thread_cell_3D(dimBlocks_3d,numBlocks_3d,dom[dev],0);

//Locate the point particle in each grid cell, store the grid cell number in points.i~points.k
    lpt_localize<<<numBlocks, dimBlocks>>>(npoints,_points[dev], _dom[dev],bc);

//initialize flow field array to 0 on device, the array length is Nx*Ny*Nz
//include scalar source and particle volume fraction divided by cell volume
////    lpt_scalar_source_init<<<numBlocks_x, dimBlocks_x>>>(_scSrc[dev],_epsp[dev], _dom[dev]);
 //   lpt_scalar_source_init<<<numBlocks_3d, dimBlocks_3d>>>(_scSrc[dev],_epsp[dev], _dom[dev]);

//lpt_scalar_source_init_test<<<numBlocks_x, dimBlocks_x>>>(_scSrc[dev], _dom[dev],ttime_done,DIFF);
lpt_scalar_source_convDiff_test<<<numBlocks_x, dimBlocks_x>>>(_scSrc[dev], _dom[dev],ttime_done,DIFF);


int coordiSys=0;
int valType=1;
//Mollify volume fraction on device, don't need too much thread source

//lpt_mollify_sc_optH(coordiSys,valType,dev,_epsp[dev]);

lpt_mollify_delta_scH(coordiSys,valType,dev,_epsp[dev]);

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
valType=0;
//lpt_mollify_sc_optH(coordiSys,valType,dev,_scSrc[dev]);

lpt_mollify_delta_scH(coordiSys,valType,dev,_scSrc[dev]);
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

//advance_sc_upwind_1st<<<numBlocks_x, dimBlocks_x>>>(DIFF, _u[dev], _v[dev], _w[dev], _scSrc[dev],_epsp[dev],  _diff0_sc[dev], _conv0_sc[dev], _diff_sc[dev], _conv_sc[dev], _sc[dev], _sc0[dev],_dom[dev],dt0_try,dt_try);
fflush(stdout);

advance_sc<<<numBlocks_x, dimBlocks_x>>>(DIFF, _u[dev], _v[dev], _w[dev], _scSrc[dev],_epsp[dev],  _diff0_sc[dev], _conv0_sc[dev], _diff_sc[dev], _conv_sc[dev], _sc[dev], _sc0[dev],_dom[dev],dt0_try,dt_try);
}
else
{
//advance_sc_upwind_1st_init<<<numBlocks_x, dimBlocks_x>>>(DIFF, _u[dev], _v[dev], _w[dev], _scSrc[dev],_epsp[dev], _diff0_sc[dev], _conv0_sc[dev], _diff_sc[dev], _conv_sc[dev], _sc[dev], _sc0[dev],_dom[dev],dt0_try,dt_try);

advance_sc_init<<<numBlocks_x, dimBlocks_x>>>(DIFF, _u[dev], _v[dev], _w[dev], _scSrc[dev],_epsp[dev], _diff0_sc[dev], _conv0_sc[dev], _diff_sc[dev], _conv_sc[dev], _sc[dev], _sc0[dev],_dom[dev],dt0_try,dt_try);
fflush(stdout);
}
*/

/*
//Using MacCormack scheme to advance scalar
advance_sc_macCormack<<<numBlocks_x, dimBlocks_x>>>(DIFF, _u[dev], _v[dev], _w[dev], _scSrc[dev],_epsp[dev], _diff_sc[dev], _conv_sc[dev], _sc[dev], _sc0[dev],_dom[dev],dt_try);
fflush(stdout);
*/

advance_sc_QUICK<<<numBlocks_x, dimBlocks_x>>>(DIFF, _u[dev], _v[dev], _w[dev], _scSrc[dev],_epsp[dev], _diff_sc[dev], _conv_sc[dev], _sc[dev], _sc0[dev],_dom[dev],dt_try);
fflush(stdout);




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
    dts[dev] = (u_max + 2 * DIFF / dom[dev].dx) / dom[dev].dx  + u_max*u_max/2/DIFF;
    dts[dev] += (v_max + 2 * DIFF / dom[dev].dy) / dom[dev].dy + v_max*v_max/2/DIFF;
    dts[dev] += (w_max + 2 * DIFF / dom[dev].dz) / dom[dev].dz + w_max*w_max/2/DIFF;
    dts[dev] = CFL / dts[dev];


*/

//MacCormack scheme ||QUICK scheme
    dts[dev] =  (u_max + 2 * DIFF / dom[dev].dx) / dom[dev].dx;
    dts[dev] += (v_max + 2 * DIFF / dom[dev].dy) / dom[dev].dy;
    dts[dev] += (w_max + 2 * DIFF / dom[dev].dz) / dom[dev].dz;
    dts[dev] = CFL / dts[dev];
/*
//1st upwind scheme
    dts[dev] =  (u_max + 2 * DIFF / dom[dev].dx) / dom[dev].dx;
    dts[dev] += (v_max + 2 * DIFF / dom[dev].dy) / dom[dev].dy;
    dts[dev] += (w_max + 2 * DIFF / dom[dev].dz) / dom[dev].dz;
    dts[dev] = CFL / dts[dev];
*/

  }

  // find max of all devices
  real max = -1.;
  for(int i = 0; i < nsubdom; i++)
    if(dts[i] > max) max = dts[i];

  // clean up
  free(dts);

  if(max>dt) max=dt; 

  return max;
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
  }
}

