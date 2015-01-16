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
/*
    int threads = MAX_THREADS_1D;
    int blocks = (int)ceil((real) npoints / (real) threads);

    dim3 dimBlocks(threads);
    dim3 numBlocks(blocks);
*/

    dim3 dimBlocks_p;
    dim3 numBlocks_p;
    block_thread_point(dimBlocks_p,numBlocks_p,npoints);

    
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


//for(int l = 0; l < 10; l++) {
   //   collision_init<<<numBlocks_p, dimBlocks_p>>>(_points[dev], npoints);
    
 /*
       for(int i = 0; i < npoints; i++) {
          collision_points<<<numBlocks_p, dimBlocks_p>>>(_points[dev], i,
            _dom[dev], eps, forces, moments, npoints, mu, bc);
        }
        spring_points<<<numBlocks_p, dimBlocks_p>>>(_points[dev], npoints);
        collision_walls<<<numBlocks_p, dimBlocks_p>>>(_dom[dev], _points[dev],
          npoints, bc, eps, mu);
      }
   */



//bc is bc.uTD etc. Make sure which BC this is. 
point_interp_init<<<numBlocks_p, dimBlocks_p>>>(npoints,_points[dev],
						ug[dev],vg[dev],wg[dev],
						lpt_stress_u[dev],lpt_stress_v[dev],lpt_stress_w[dev],scg[dev]);
fflush(stdout);

interpolate_point_vel_Lag2<<<numBlocks_p, dimBlocks_p>>>(_u[dev],_v[dev],_w[dev],
							 npoints,rho_f,nu,
							 ug[dev],vg[dev],wg[dev],
							_points[dev],_dom[dev],bc);
fflush(stdout);

interpolate_point_scalar_Lag2<<<numBlocks_p, dimBlocks_p>>>(npoints,_sc[dev],scg[dev],_points[dev],_dom[dev]);
fflush(stdout);



/*
C_add=0.5;
C_stress=1;
C_drag=1;
*/

//get lpt_stress
//TODO _stress_u is not available near the boundary(set to 0), while lpt_stress can be interpolated on BC
if(C_stress>0||C_add>0) 
interpolate_point_vel_Lag2<<<numBlocks_p, dimBlocks_p>>>(_stress_u[dev],_stress_v[dev],_stress_w[dev],
							 npoints, rho_f, nu,
							 lpt_stress_u[dev],lpt_stress_v[dev],lpt_stress_w[dev],
							 _points[dev],_dom[dev],bc);


drag_points<<<numBlocks_p, dimBlocks_p>>>(_points[dev],npoints,
ug[dev],vg[dev],wg[dev],
lpt_stress_u[dev],lpt_stress_v[dev],lpt_stress_w[dev],scg[dev],
rho_f,mu,g,gradP,
C_add, C_stress,C_drag,
sc_eq,DIFF);
fflush(stdout);

      move_points_a<<<numBlocks_p, dimBlocks_p>>>(_points[dev], npoints,dt_try);

drag_points<<<numBlocks_p, dimBlocks_p>>>(_points[dev],npoints,
ug[dev],vg[dev],wg[dev],
lpt_stress_u[dev],lpt_stress_v[dev],lpt_stress_w[dev],scg[dev],
rho_f,mu,g,gradP,
C_add, C_stress,C_drag,
sc_eq,DIFF);
fflush(stdout);

      move_points_b<<<numBlocks_p, dimBlocks_p>>>(_dom[dev], _points[dev], npoints,dt_try);

fflush(stdout);




 /*
      checkCudaErrors(cudaFree(forces));
      checkCudaErrors(cudaFree(moments));
 */
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

int valType=0;
/* 
lpt_mollify_sc_optH(1,valType,dev,_f_x[dev]);
lpt_mollify_sc_optH(2,valType,dev,_f_y[dev]);
lpt_mollify_sc_optH(3,valType,dev,_f_z[dev]);
*/
lpt_mollify_delta_scH(1,valType,dev,_f_x[dev]);
lpt_mollify_delta_scH(2,valType,dev,_f_y[dev]);
lpt_mollify_delta_scH(3,valType,dev,_f_z[dev]);


     }
}


   













//extern "C"
void sortParticles(int *dgridParticleHash, int *dgridParticleIndex, int numParticles)
{
//sort dgridParticleIndex based on dgridParticleHash, and they both change order
        thrust::sort_by_key(thrust::device_ptr<int>(dgridParticleHash),
                            thrust::device_ptr<int>(dgridParticleHash + numParticles),
                            thrust::device_ptr<int>(dgridParticleIndex));
}

extern "C"
void lpt_mollify_delta_scH(int coordiSys,int valType,int dev,real *scSrc)
{
//    dim3 dimBlocks_3d,numBlocks_3d;
    dim3 dimBlocks_w,numBlocks_w;
    dim3 dimBlocks_p,numBlocks_p;
//    dim3 numBlocks_st;

//    dim3 numBlocks_print,dimBlocks_print;

//int coordiSys=0;//coordinate systerm, cell-center or face center
int planeDirc=3;//parallel x-y or y-z or x-z plane

block_thread_cell(dimBlocks_w,numBlocks_w,dom[dev],coordiSys,planeDirc);
//block_thread_cell_3D(dimBlocks_3d,numBlocks_3d,dom[dev],coordiSys);
block_thread_point(dimBlocks_p,numBlocks_p,npoints);
//block_thread_point(dimBlocks_p,numBlocks_st,npoints*STENCIL3);


int lenCell;
switch(coordiSys)
{
case 0:lenCell=dom[dev].Gcc.s3b;break;
case 1:lenCell=dom[dev].Gfx.s3b;break;
case 2:lenCell=dom[dev].Gfy.s3b;break;
case 3:lenCell=dom[dev].Gfz.s3b;break;
default: break;
}

//block_thread_point(dimBlocks_print,numBlocks_print,lenCell);

checkCudaErrors(cudaMemset(gridParticleHash[dev],-1,npoints*sizeof(int)));
checkCudaErrors(cudaMemset(gridParticleIndex[dev],-1,npoints*sizeof(int)));

checkCudaErrors(cudaMemset(cellStart[dev],-1,lenCell*sizeof(int)));
checkCudaErrors(cudaMemset(cellEnd[dev],-1,lenCell*sizeof(int)));


cudaEvent_t start, stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
float milliseconds = 0;

//printf("\ngridDim_p %d %d %d\n",numBlocks_p.x,dimBlocks_p.x,npoints);
//printf("\ngridDim_3d %d %d %d %d\n",numBlocks_3d.x,numBlocks_3d.y,numBlocks_3d.z,lenCell);

lpt_delta_point_position<<<numBlocks_p,dimBlocks_p>>>(_points[dev],_dom[dev],
posXold[dev],posYold[dev],posZold[dev],
lptSourceValOld[dev],npoints,coordiSys,valType);
 

calcHash_optD<<<numBlocks_p,dimBlocks_p>>>(gridParticleHash[dev],
					gridParticleIndex[dev],
					_dom[dev],
					posXold[dev],posYold[dev],posZold[dev],
					npoints,coordiSys);


sortParticles(gridParticleHash[dev],gridParticleIndex[dev],npoints);

milliseconds = 0;
cudaEventRecord(start);
findCellStart_deltaD<<<numBlocks_p,dimBlocks_p>>>(cellStart[dev],
						cellEnd[dev],
						gridParticleHash[dev],
						gridParticleIndex[dev],
						lptSourceVal[dev],
						lptSourceValOld[dev],
						npoints);

cudaEventRecord(stop);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&milliseconds, start, stop);
printf("\ntime_reoder %f\n",milliseconds);
fflush(stdout);
milliseconds = 0;
cudaEventRecord(start);
//particle volume fraction 1 or other cell-centerred parameter 0
//lpt_point_weight<<<numBlocks_p,dimBlocks_p>>>(_points[dev],_dom[dev],posX[dev],posY[dev],posZ[dev],Weight[dev],gridParticleIndex[dev],npoints,coordiSys,valType);

//lpt_point_ksi_opt<<<numBlocks_p,dimBlocks_p>>>(_points[dev],_dom[dev],posX[dev],posY[dev],posZ[dev],Ksi[dev],gridParticleIndex[dev],npoints,coordiSys,valType);
fflush(stdout);

cudaEventRecord(stop);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&milliseconds, start, stop);
printf("\ntime_weight %f\n",milliseconds);
fflush(stdout);

milliseconds = 0;
cudaEventRecord(start);


lpt_mollify_delta_scD<<<numBlocks_w,dimBlocks_w>>>(_dom[dev],scSrc,lptSourceVal[dev],cellStart[dev],cellEnd[dev],npoints,coordiSys);


fflush(stdout);
cudaEventRecord(stop);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&milliseconds, start, stop);
printf("\ntime_mollify %f\n",milliseconds);
fflush(stdout);

//print_kernel_array_int<<<numBlocks_print,dimBlocks_print>>>(cellEnd[dev],lenCell);
}


extern "C"
void lpt_mollify_sc_optH(int coordiSys,int valType,int dev,real *scSrc)
{
    dim3 dimBlocks_3d,numBlocks_3d;
    dim3 dimBlocks_w,numBlocks_w;
    dim3 dimBlocks_p,numBlocks_p;
//    dim3 numBlocks_st;

//    dim3 numBlocks_print,dimBlocks_print;

//int coordiSys=0;//coordinate systerm, cell-center or face center
int planeDirc=3;//parallel x-y or y-z or x-z plane

block_thread_cell(dimBlocks_w,numBlocks_w,dom[dev],coordiSys,planeDirc);
block_thread_cell_3D(dimBlocks_3d,numBlocks_3d,dom[dev],coordiSys);
block_thread_point(dimBlocks_p,numBlocks_p,npoints);
//block_thread_point(dimBlocks_p,numBlocks_st,npoints*STENCIL3);


int lenCell;
switch(coordiSys)
{
case 0:lenCell=dom[dev].Gcc.s3b;break;
case 1:lenCell=dom[dev].Gfx.s3b;break;
case 2:lenCell=dom[dev].Gfy.s3b;break;
case 3:lenCell=dom[dev].Gfz.s3b;break;
default: break;
}

//block_thread_point(dimBlocks_print,numBlocks_print,lenCell);

checkCudaErrors(cudaMemset(gridParticleHash[dev],-1,npoints*sizeof(int)));
checkCudaErrors(cudaMemset(gridParticleIndex[dev],-1,npoints*sizeof(int)));

checkCudaErrors(cudaMemset(cellStart[dev],-1,lenCell*sizeof(int)));
checkCudaErrors(cudaMemset(cellEnd[dev],-1,lenCell*sizeof(int)));


cudaEvent_t start, stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
float milliseconds = 0;

//printf("\ngridDim_p %d %d %d\n",numBlocks_p.x,dimBlocks_p.x,npoints);
//printf("\ngridDim_3d %d %d %d %d\n",numBlocks_3d.x,numBlocks_3d.y,numBlocks_3d.z,lenCell);

//array_init<<<numBlocks_st, dimBlocks_p>>>(Ksi[dev],_dom[dev],npoints*STENCIL3, 0.);
/*
array_init<<<numBlocks_p, dimBlocks_p>>>(posX[dev],_dom[dev],npoints, -10.f);
array_init<<<numBlocks_p, dimBlocks_p>>>(posY[dev],_dom[dev],npoints, -10.f);
array_init<<<numBlocks_p, dimBlocks_p>>>(posZ[dev],_dom[dev],npoints, -10.f);
array_init<<<numBlocks_p, dimBlocks_p>>>(posXold[dev],_dom[dev],npoints, -10.f);
array_init<<<numBlocks_p, dimBlocks_p>>>(posYold[dev],_dom[dev],npoints, -10.f);
array_init<<<numBlocks_p, dimBlocks_p>>>(posZold[dev],_dom[dev],npoints, -10.f);
array_init<<<numBlocks_p, dimBlocks_p>>>(Weight[dev],_dom[dev],npoints, 0.f);
*/

lpt_point_position<<<numBlocks_p,dimBlocks_p>>>(_points[dev],posXold[dev],posYold[dev],posZold[dev],npoints);

calcHash_optD<<<numBlocks_p,dimBlocks_p>>>(gridParticleHash[dev],
					gridParticleIndex[dev],
					_dom[dev],
					posXold[dev],posYold[dev],posZold[dev],
					npoints,coordiSys);


sortParticles(gridParticleHash[dev],gridParticleIndex[dev],npoints);

milliseconds = 0;
cudaEventRecord(start);
findCellStart_optD<<<numBlocks_p,dimBlocks_p>>>(cellStart[dev],
						cellEnd[dev],
						gridParticleHash[dev],
						gridParticleIndex[dev],
						posX[dev],posY[dev],posZ[dev],
						posXold[dev],posYold[dev],posZold[dev],
						npoints);

cudaEventRecord(stop);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&milliseconds, start, stop);
printf("\ntime_reoder %f\n",milliseconds);
fflush(stdout);
milliseconds = 0;
cudaEventRecord(start);
//particle volume fraction 1 or other cell-centerred parameter 0
//lpt_point_weight<<<numBlocks_p,dimBlocks_p>>>(_points[dev],_dom[dev],posX[dev],posY[dev],posZ[dev],Weight[dev],gridParticleIndex[dev],npoints,coordiSys,valType);

lpt_point_ksi_opt<<<numBlocks_p,dimBlocks_p>>>(_points[dev],_dom[dev],posX[dev],posY[dev],posZ[dev],Ksi[dev],gridParticleIndex[dev],npoints,coordiSys,valType);
fflush(stdout);

cudaEventRecord(stop);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&milliseconds, start, stop);
printf("\ntime_weight %f\n",milliseconds);
fflush(stdout);

milliseconds = 0;
cudaEventRecord(start);

//lpt_mollify_sc_optD<<<numBlocks_w,dimBlocks_w>>>(_dom[dev],scSrc,posX[dev],posY[dev],posZ[dev],Weight[dev],cellStart[dev],cellEnd[dev],gridParticleIndex[dev],npoints,coordiSys,valType);

//lpt_mollify_sc_ksi_optD<<<numBlocks_3d,dimBlocks_3d>>>(_dom[dev],scSrc,posX[dev],posY[dev],posZ[dev],Ksi[dev],cellStart[dev],cellEnd[dev],gridParticleIndex[dev],npoints,coordiSys,valType);

lpt_mollify_sc_ksi_optD<<<numBlocks_w,dimBlocks_w>>>(_dom[dev],scSrc,posX[dev],posY[dev],posZ[dev],Ksi[dev],cellStart[dev],cellEnd[dev],gridParticleIndex[dev],npoints,coordiSys,valType);


fflush(stdout);
cudaEventRecord(stop);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&milliseconds, start, stop);
printf("\ntime_mollify %f\n",milliseconds);
fflush(stdout);

//print_kernel_array_int<<<numBlocks_print,dimBlocks_print>>>(cellEnd[dev],lenCell);
}


//About Swap, and reference, dirc is the system direction, dirc2 is the plane direction
void block_thread_cell(dim3 &dimBlocks,dim3 &numBlocks,dom_struct dom,int dirc,int dirc2)
{

    int threads_y = 0;
    int threads_x = 0;
    int blocks_y = 0;
    int blocks_x = 0;

    int lenX=0;
    int lenY=0;

grid_info G;
switch(dirc)
{
case 0:G=dom.Gcc;break;
case 1:G=dom.Gfx;break;
case 2:G=dom.Gfy;break;
case 3:G=dom.Gfz;break;
default: break;
}

	switch(dirc2)
	{
	case 1:
		lenX=G._jnb;
		lenY=G._knb;
		break;
	case 2:
		lenX=G._knb;
		lenY=G._inb;
		break;
	case 3:
		lenX=G._inb;
		lenY=G._jnb;
		break;
	default: break;	
	}


    if(lenX < MAX_THREADS_DIM)
      threads_x = lenX+2;
    else
      threads_x = MAX_THREADS_DIM;

    if(lenY < MAX_THREADS_DIM)
      threads_y = lenY+2;
    else
      threads_y = MAX_THREADS_DIM;


    blocks_x = (int)ceil((real) lenX / (real) (threads_x-2));
    blocks_y = (int)ceil((real) lenY / (real) (threads_y-2));

    dimBlocks.x=threads_x;
    dimBlocks.y=threads_y;
    numBlocks.x=blocks_x;
    numBlocks.y=blocks_y;

}

/*
//About Swap, and reference, dirc is the system direction, get 3d blocks and threads
void block_thread_cell_3D(dim3 &dimBlocks,dim3 &numBlocks,dom_struct dom,int dirc)
{

    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    int lenX=0;
    int lenY=0;
    int lenZ=0;

grid_info G;
switch(dirc)
{
case 0:G=dom.Gcc;break;
case 1:G=dom.Gfx;break;
case 2:G=dom.Gfy;break;
case 3:G=dom.Gfz;break;
default: break;
}

		lenX=G._inb;
		lenY=G._jnb;
		lenZ=G._knb;


    if(lenX < MAX_THREADS_DIM)
      threads_x = lenX;
    else
      threads_x = MAX_THREADS_DIM;

    if(lenY < MAX_THREADS_DIM)
      threads_y = lenY;
    else
      threads_y = MAX_THREADS_DIM;

    int MAX_THREADS_Z=(int) 1024/(1.0f*MAX_THREADS_DIM*MAX_THREADS_DIM);
    if(lenZ < MAX_THREADS_Z)
      threads_z = lenZ;
    else
      threads_z = MAX_THREADS_Z;

    blocks_x = (int)ceil((real) lenX / (real) (threads_x));
    blocks_y = (int)ceil((real) lenY / (real) (threads_y));
    blocks_z = (int)ceil((real) lenZ / (real) (threads_z));

    dimBlocks.x=threads_x;
    dimBlocks.y=threads_y;
    dimBlocks.z=threads_z;

    numBlocks.x=blocks_x;
    numBlocks.y=blocks_y;
    numBlocks.z=blocks_z;
}

*/


//About Swap, and reference, dirc is the system direction, get 3d blocks and threads
//This method is faster than blockDim=16*16*4 by 10%
void block_thread_cell_3D(dim3 &dimBlocks,dim3 &numBlocks,dom_struct dom,int dirc)
{

    int threads_x = 0;
    int threads_y = 0;
    int threads_z = 0;
    int blocks_x = 0;
    int blocks_y = 0;
    int blocks_z = 0;

    int lenX=0;
    int lenY=0;
    int lenZ=0;

grid_info G;
switch(dirc)
{
case 0:G=dom.Gcc;break;
case 1:G=dom.Gfx;break;
case 2:G=dom.Gfy;break;
case 3:G=dom.Gfz;break;
default: break;
}

		lenX=G._inb;
		lenY=G._jnb;
		lenZ=G._knb;


    if(lenX < MAX_THREADS_DIM3)
      threads_x = lenX;
    else
      threads_x = MAX_THREADS_DIM3;

    if(lenY < MAX_THREADS_DIM3)
      threads_y = lenY;
    else
      threads_y = MAX_THREADS_DIM3;

    if(lenZ < MAX_THREADS_DIM3)
      threads_z = lenZ;
    else
      threads_z = MAX_THREADS_DIM3;

    blocks_x = (int)ceil((real) lenX / (real) (threads_x));
    blocks_y = (int)ceil((real) lenY / (real) (threads_y));
    blocks_z = (int)ceil((real) lenZ / (real) (threads_z));

    dimBlocks.x=threads_x;
    dimBlocks.y=threads_y;
    dimBlocks.z=threads_z;

    numBlocks.x=blocks_x;
    numBlocks.y=blocks_y;
    numBlocks.z=blocks_z;
}

void block_thread_point(dim3 &dimBlocks,dim3 &numBlocks,int npoints)
{
    int threads = MAX_THREADS_1D;
    int blocks = (int)ceil((real) npoints / (real) threads);

    dimBlocks.x=threads;
    numBlocks.x=blocks;

}


void cuda_malloc_array_int(int** &A,int lenArray)
{
A= (int**) malloc(nsubdom * sizeof(int*));
          cpumem += nsubdom * sizeof(int*);

  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

   checkCudaErrors(cudaMalloc((void**) &(A[dev]), sizeof(int) * lenArray));
    gpumem += sizeof(int) * lenArray;
  }
}

void cuda_malloc_array_real(real**& A,int lenArray)
{   
A= (real**) malloc(nsubdom * sizeof(real*));
          cpumem += nsubdom * sizeof(real*);
      
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));
    

   checkCudaErrors(cudaMalloc((void**) &(A[dev]), sizeof(real)*lenArray));
    gpumem += sizeof(real) * lenArray;
  } 
}   


void cuda_free_array_real(real**& A)
{
  // free device memory on device
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));
    checkCudaErrors(cudaFree(A[dev]));
}

  free(A);

}

void cuda_free_array_int(int** &A)
{
  // free device memory on device
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

    checkCudaErrors(cudaFree(A[dev]));
}   
    
  free(A);

}

extern "C"
void cuda_point_malloc(void)
{
  // allocate device memory on host
  _points = (point_struct**) malloc(nsubdom * sizeof(point_struct*));
  cpumem += nsubdom * sizeof(point_struct*);

  _flag_u = (int**) malloc(nsubdom * sizeof(int*));
  cpumem += nsubdom * sizeof(int*);
  _flag_v = (int**) malloc(nsubdom * sizeof(int*));
  cpumem += nsubdom * sizeof(int*);
  _flag_w = (int**) malloc(nsubdom * sizeof(int*));
  cpumem += nsubdom * sizeof(int*);

  lpt_point_source_mollify_init();


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

  lpt_point_source_mollify_final();

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



//extern "C"
void lpt_point_source_mollify_init()
{
cuda_malloc_array_real(ug,npoints);
cuda_malloc_array_real(vg,npoints);
cuda_malloc_array_real(wg,npoints);

cuda_malloc_array_real(posX,npoints);
cuda_malloc_array_real(posY,npoints);
cuda_malloc_array_real(posZ,npoints);
cuda_malloc_array_real(posXold,npoints);
cuda_malloc_array_real(posYold,npoints);
cuda_malloc_array_real(posZold,npoints);

cuda_malloc_array_real(lptSourceVal,npoints);
cuda_malloc_array_real(lptSourceValOld,npoints);

cuda_malloc_array_real(lpt_stress_u,npoints);
cuda_malloc_array_real(lpt_stress_v,npoints);
cuda_malloc_array_real(lpt_stress_w,npoints);

cuda_malloc_array_real(scg,npoints);
cuda_malloc_array_real(Weight,npoints);
cuda_malloc_array_real(Ksi,npoints*STENCIL3);

cuda_malloc_array_int(gridParticleIndex,npoints);
cuda_malloc_array_int(gridParticleHash,npoints);

//calculate the maximum length of coordinate system
int lenCell=Dom.Gcc.s3b;
int len1=Dom.Gfx.s3b;
int len2=Dom.Gfy.s3b;
int len3=Dom.Gfz.s3b;
if(lenCell<len1) lenCell=len1;
if(len2<len3) len2=len3;
if(lenCell<len2) lenCell=len2;
//lenCell=max(max(len1,len2),max(len3,lenCell));
cuda_malloc_array_int(cellStart,lenCell);
cuda_malloc_array_int(cellEnd,lenCell);


  // allocate device memory on device
  #pragma omp parallel num_threads(nsubdom)
  {
    int dev = omp_get_thread_num();
    checkCudaErrors(cudaSetDevice(dev + dev_start));

    //Initialize those two arrays since they will be initialized with smaller length later in lpt_mollify_scH
    checkCudaErrors(cudaMemset(cellStart[dev],-1,lenCell*sizeof(int)));
    checkCudaErrors(cudaMemset(cellEnd[dev],-1,lenCell*sizeof(int)));
  }

gaussian_array_initH();
domInfo_array_initH();
}

texture<float,1,cudaReadModeElementType> texRefGaussian;
texture<int,1,cudaReadModeElementType> texRefDomInfo;

void domInfo_array_initH()
{
//ref the structure of grid_info in domain.h
int indexTypeLen=3;//0~2 correspond to is,ie,in  if indexDir<3; 0~2 correspond to s1,s2,s3 if indexDir=3
int incGhostLen=2;//0~1 correspond to include ghost cell or not
int indexDirLen=4;//0~2 correspond to i,j,k
int hostDeviceLen=2; //0~1 correspond to host,device
int coordiSysLen=4; //0~3 correspond to Gcc,Gfx,Gfy,Gfz.
int domInfo_array_len=coordiSysLen*hostDeviceLen*indexDirLen*incGhostLen*indexTypeLen;
 
int *DomInfo;

DomInfo=(int *) malloc(sizeof(int) * domInfo_array_len);
checkCudaErrors(cudaMalloc((void**) &(_DomInfo),sizeof(int) * domInfo_array_len));

grid_info G;
int indexType=3;//0~2 correspond to is,ie,in  if indexDir<3; 0~2 correspond to s1,s2,s3 if indexDir=3
int incGhost=2;//0~1 correspond to include ghost cell or not
int indexDir=4;//0~2 correspond to i,j,k
int hostDevice=2; //0~1 correspond to host,device
int coordiSys=4; //0~3 correspond to Gcc,Gfx,Gfy,Gfz.
int index; 
//coordiSys= 0~3 correspond to Gcc,Gfx,Gfy,Gfz.
for(coordiSys=0;coordiSys<coordiSysLen;coordiSys++)
{
switch(coordiSys)
{
case 0:G=Dom.Gcc;break;
case 1:G=Dom.Gfx;break;
case 2:G=Dom.Gfy;break;
case 3:G=Dom.Gfz;break;
default: break;
}

hostDevice=0;
indexDir=0;
incGhost=0;
indexType=0;
index=indexType+incGhost*3+indexDir*6+ hostDevice*24+ coordiSys*48;

DomInfo[index]  =G.is;
DomInfo[index+1]=G.ie;
DomInfo[index+2]=G.in;
DomInfo[index+3]=G.isb;
DomInfo[index+4]=G.ieb;
DomInfo[index+5]=G.inb;

DomInfo[index+6] =G.js;
DomInfo[index+7] =G.je;
DomInfo[index+8] =G.jn;
DomInfo[index+9] =G.jsb;
DomInfo[index+10]=G.jeb;
DomInfo[index+11]=G.jnb;

DomInfo[index+12] =G.ks;
DomInfo[index+13] =G.ke;
DomInfo[index+14] =G.kn;
DomInfo[index+15] =G.ksb;
DomInfo[index+16] =G.keb;
DomInfo[index+17] =G.knb;

DomInfo[index+18] =G.s1;
DomInfo[index+19] =G.s2;
DomInfo[index+20] =G.s3;
DomInfo[index+21] =G.s1b;
DomInfo[index+22] =G.s2b;
DomInfo[index+23] =G.s3b;
//printf("\nDomInfo: %d %d %d %d\n",G.s1b,G.s2b,G.s3b,coordiSys);
DomInfo[index+24]=G._is;
DomInfo[index+25]=G._ie;
DomInfo[index+26]=G._in;
DomInfo[index+27]=G._isb;
DomInfo[index+28]=G._ieb;
DomInfo[index+29]=G._inb;

DomInfo[index+30]=G._js;
DomInfo[index+31]=G._je;
DomInfo[index+32]=G._jn;
DomInfo[index+33]=G._jsb;
DomInfo[index+34]=G._jeb;
DomInfo[index+35]=G._jnb;

DomInfo[index+36] =G._ks;
DomInfo[index+37] =G._ke;
DomInfo[index+38] =G._kn;
DomInfo[index+39] =G._ksb;
DomInfo[index+40] =G._keb;
DomInfo[index+41] =G._knb;

DomInfo[index+42] =G._s1;
DomInfo[index+43] =G._s2;
DomInfo[index+44] =G._s3;
DomInfo[index+45] =G._s1b;
DomInfo[index+46] =G._s2b;
DomInfo[index+47] =G._s3b;
//printf("\ndevice DomInfo: %d %d %d %d\n",G._s1b,G._s2b,G._s3b,coordiSys); //_s1b~_s3b =0 at this time
}


checkCudaErrors(cudaMemcpy(_DomInfo, DomInfo, sizeof(int) * domInfo_array_len,
      cudaMemcpyHostToDevice));
cudaBindTexture(0,texRefDomInfo,_DomInfo,sizeof(int)*domInfo_array_len);

free(DomInfo);
}

void gaussian_array_initH()
{
//put into malloc and free
//    int LEN_GAUSSIAN_ARRAY=5000;
    //texture<float,1,cudaReadModeElementType> texRefGaussian;

    checkCudaErrors(cudaMalloc((void**) &(GaussianKernel),sizeof(float) * LEN_GAUSSIAN_ARRAY));

    dim3 dimBlocks_p,numBlocks_p;
    block_thread_point(dimBlocks_p,numBlocks_p,LEN_GAUSSIAN_ARRAY);


real dx=Dom.dx;real dy=Dom.dy;real dz=Dom.dz;    
real min_meshsize=min(min(dx,dy),dz);
real max_meshsize=max(max(dx,dy),dz);

//real min_meshsize=min(min(dx,dy),dz);
//2.0f*sqrt(2.0f*log(2.0f))=2.3548f;
real sig= min_meshsize/2.3548f;
real norm=pow(sqrt(2*PI)*sig,3);
real maxLenR=2.6f*max_meshsize;
real maxLenR2=maxLenR*maxLenR;
//real dg=maxLenR/(float)(LEN_GAUSSIAN_ARRAY);

real dg2=maxLenR2/(float)(LEN_GAUSSIAN_ARRAY);

//real dg2_sig2=dg*dg/2.f/sig/sig;
real dg2_sig2=dg2/2.f/sig/sig;

    gaussian_array_initD<<<numBlocks_p,dimBlocks_p>>>(GaussianKernel,dg2_sig2,dg2,norm);
fflush(stdout);
    cudaBindTexture(0,texRefGaussian,GaussianKernel,sizeof(float)*LEN_GAUSSIAN_ARRAY);
}



//extern "C"
void lpt_point_source_mollify_final()
{
cudaFree(GaussianKernel);
cudaUnbindTexture(texRefGaussian);

cudaFree(_DomInfo);
cudaUnbindTexture(texRefDomInfo);

cuda_free_array_real(posX);
cuda_free_array_real(posY);
cuda_free_array_real(posZ);

cuda_free_array_real(posXold);
cuda_free_array_real(posYold);
cuda_free_array_real(posZold);

cuda_free_array_real(lptSourceVal);
cuda_free_array_real(lptSourceValOld);

cuda_free_array_real(ug);
cuda_free_array_real(vg);
cuda_free_array_real(wg);

cuda_free_array_real(scg);
cuda_free_array_real(Ksi);
cuda_free_array_real(Weight);

cuda_free_array_real(lpt_stress_u);
cuda_free_array_real(lpt_stress_v);
cuda_free_array_real(lpt_stress_w);

cuda_free_array_int(cellStart);
cuda_free_array_int(cellEnd);
cuda_free_array_int(gridParticleIndex);
cuda_free_array_int(gridParticleHash);
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


