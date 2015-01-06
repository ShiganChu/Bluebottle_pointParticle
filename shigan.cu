#include "cuda_bicgstab.h"
#include "bluebottle.h" //contain _u,_v,_w, rho_f,nu,dt0,dt
#include "shigan.h" //contain _u,_v,_w, rho_f,nu,dt0,dt

#include "entrySearch.h" //contain _u,_v,_w, rho_f,nu,dt0,dt

#include <cusp/dia_matrix.h>

#include <cuda.h>
#include <helper_cuda.h>


void shigan_max(int size, real *_stress)
{
//int size=dom[dev].Gfx.s3b;
real *_max_u;
real *max_u=(real*) malloc(size * sizeof(real));
checkCudaErrors(cudaMalloc((void**) &_max_u, sizeof(real) * size));

checkCudaErrors(cudaMemcpy(_max_u,_stress, size*sizeof(real),cudaMemcpyDeviceToDevice));
checkCudaErrors(cudaMemcpy(max_u,_max_u, size*sizeof(real),cudaMemcpyDeviceToHost));
 
//for(int ii=0;ii<size;ii++) printf(" %f  ",max_u[ii]);

//printf("Array_previous %f %f %f %f\n",max_u[0],max_u[1],max_u[2],max_u[3]);
//fflush(stdout);
 

//could  change this line to find max or min!!
real maxValue=find_max(size,_max_u);

printf("max %f \n",maxValue);
fflush(stdout);

free(max_u);
cudaFree(_max_u);
}

void shigan_min(int size, real *_stress)
{
//int size=dom[dev].Gfx.s3b;
real *_max_u;
real *max_u=(real*) malloc(size * sizeof(real));
checkCudaErrors(cudaMalloc((void**) &_max_u, sizeof(real) * size));

checkCudaErrors(cudaMemcpy(_max_u,_stress, size*sizeof(real),cudaMemcpyDeviceToDevice));
checkCudaErrors(cudaMemcpy(max_u,_max_u, size*sizeof(real),cudaMemcpyDeviceToHost));
/*
printf("Array_previous %f %f %f %f\n",max_u[0],max_u[1],max_u[2],max_u[3]);
fflush(stdout);
*/
//could  change this line to find max or min!!
real maxValue=find_min(size,_max_u);

printf("min %f \n",maxValue);
fflush(stdout);

free(max_u);
cudaFree(_max_u);
}



void cuda_flow_stress()
{
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

//printf("max and min\n");

//printf("p\n");
//shigan_max(dom[dev].Gcc.s3b, _p[dev]);
//shigan_min(dom[dev].Gcc.s3b, _p[dev]);
//
//printf("u\n");
//shigan_max(dom[dev].Gfx.s3b, _u[dev]);
//shigan_min(dom[dev].Gfx.s3b, _u[dev]);
//
//printf("w\n");
//shigan_max(dom[dev].Gfz.s3b, _w[dev]);
//shigan_min(dom[dev].Gfz.s3b, _w[dev]);

    stress_u<<<numBlocks_u, dimBlocks_u>>>(rho_f, nu,_u[dev],_p[dev],_p0[dev], _stress_u[dev], _dom[dev],_flag_u[dev],dt,dt0);
fflush(stdout);

//printf("max and min stress_u\n");
//shigan_max(dom[dev].Gfx.s3b, _stress_u[dev]);
//shigan_min(dom[dev].Gfx.s3b, _stress_u[dev]);




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

//extern "C"
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
      if(npoints > 0) {
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
      
point_struct *test=(point_struct *)malloc(sizeof(point_struct));

//interpolate_point_vel<<<numBlocks, dimBlocks>>>(_u[dev],_v[dev],_w[dev],npoints,rho_f,nu,ug,vg,wg,_points[dev],_dom[dev],dt0,ttime,bc);
fflush(stdout);

interpolate_point_vel_shigan<<<numBlocks, dimBlocks>>>(_u[dev],_v[dev],_w[dev],npoints,rho_f,nu,ug,vg,wg,_points[dev],_dom[dev],dt0,ttime,bc);
fflush(stdout);



 //need further modification for calculating fluid stress!!
//interpolate_point_vel<<<numBlocks, dimBlocks>>>(_u0[dev],_v0[dev],_w0[dev],npoints,rho_f,nu,ug0,vg0,wg0,_points[dev],_dom[dev],dt0,dt,bc);
//interpolate_point_vel<<<numBlocks, dimBlocks>>>(_conv_u[dev],_conv_v[dev],_conv_w[dev],npoints,rho_f,nu,conv_ug,conv_vg,conv_wg,_points[dev],_dom[dev],dt0,dt,bc);


//  checkCudaErrors(cudaMemcpy(test,_points[dev], sizeof(point_struct), cudaMemcpyDeviceToHost));

//real C_add=0.5;
//real C_add=0;
//real C_stress=1;
//real C_stress=0;
//real C_drag=0;
//real C_drag=1;

//get lpt_stress
//TODO _stress_u is not available near the boundary(set to 0), while lpt_stress can be interpolated on BC
if(C_stress>0||C_add>0) interpolate_point_vel_shigan<<<numBlocks, dimBlocks>>>(_stress_u[dev],_stress_v[dev],_stress_w[dev],npoints,rho_f,nu,lpt_stress_u,lpt_stress_v,lpt_stress_w,_points[dev],_dom[dev],dt0,dt,bc);

//if(C_stress>0||C_add>0) interpolate_point_vel<<<numBlocks, dimBlocks>>>(_stress_u[dev],_stress_v[dev],_stress_w[dev],npoints,rho_f,nu,lpt_stress_u,lpt_stress_v,lpt_stress_w,_points[dev],_dom[dev],dt0,dt,bc);

drag_points<<<numBlocks, dimBlocks>>>(_points[dev],npoints,
ug,vg,wg,
lpt_stress_u,lpt_stress_v,lpt_stress_w,
rho_f,mu,g,gradP,
C_add, C_stress,C_drag);
fflush(stdout);


      move_points_b<<<numBlocks, dimBlocks>>>(_dom[dev], _points[dev], npoints,
        dt, dt0, g, rho_f, ttime,
C_add, C_stress,C_drag);

fflush(stdout);
//printf("w is %f wdot is %f\n",test->w,test->wdot);
//printf("w is %f wdot is %f\n",test->x,test->y);
//fflush(stdout);




 /*
      checkCudaErrors(cudaFree(forces));
      checkCudaErrors(cudaFree(moments));
 */
      checkCudaErrors(cudaFree(ug));
      checkCudaErrors(cudaFree(vg));
      checkCudaErrors(cudaFree(wg));

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
}

//extern "C"
/*
void cgns_point_particles(real dtout)
{
  if(npoints > 0) {
    // create the solution file
    char fname[FILE_NAME_SIZE];
    char fname2[FILE_NAME_SIZE];
    char fnameall[FILE_NAME_SIZE];
    char fnameall2[FILE_NAME_SIZE];
    real tout = ttime; // = rec_point_particle_stepnum_out * dtout;
    char format[CHAR_BUF_SIZE];
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
    real *x = (real *)malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real);
    real *y = (real *)malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    real *z = (real *)malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    cgsize_t *conn = (cgsize_t *) malloc(npoints * sizeof(cgsize_t));
    // cpumem (real *)+= npoints * sizeof(int);
    real *a = (real *)malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(int);
    real *u = (real *)malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    real *v = (real *)malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    real *w = (real *)malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    real *udot =(real *) malloc(npoints * sizeof(real));
    // cpumem +=(real *) npoints * sizeof(real);
    real *vdot =(real *) malloc(npoints * sizeof(real));
    // cpumem +=(real *) npoints * sizeof(real);
    real *wdot =(real *) malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    real *ox =(real *) malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    real *oy =(real *) malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    real *oz =(real *) malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    real *Fx =(real *) malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    real *Fy =(real *) malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    real *Fz =(real *) malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    real *Lx =(real *) malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    real *Ly =(real *) malloc(npoints * sizeof(real));
    // cpumem (real *)+= npoints * sizeof(real);
    real *Lz =(real *) malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real); 
    real *iFx =(real *) malloc(npoints * sizeof(real));
    // cpumem +(real *)= npoints * sizeof(real);
    real *iFy =(real *) malloc(npoints * sizeof(real));
    // cpumem +(real *)= npoints * sizeof(real);
    real *iFz =(real *) malloc(npoints * sizeof(real));
    // cpumem +(real *)= npoints * sizeof(real);
    real *iLx =(real *) malloc(npoints * sizeof(real));
    // cpumem +(real *)= npoints * sizeof(real);
    real *iLy =(real *) malloc(npoints * sizeof(real));
    // cpumem +(real *)= npoints * sizeof(real);
    real *iLz =(real *) malloc(npoints * sizeof(real));
    // cpumem +(real *)= npoints * sizeof(real); 
    real *hFx =(real *) malloc(npoints * sizeof(real));
    // cpumem +(real *)= npoints * sizeof(real);
    real *hFy =(real *) malloc(npoints * sizeof(real));
    // cpumem +(real *)= npoints * sizeof(real);
    real *hFz =(real *) malloc(npoints * sizeof(real));
    // cpumem +(real *)= npoints * sizeof(real);
    real *hLx =(real *) malloc(npoints * sizeof(real));
    // cpumem +(real *)= npoints * sizeof(real);
    real *hLy =(real *) malloc(npoints * sizeof(real));
    // cpumem +(real *)= npoints * sizeof(real);
    real *hLz =(real *) malloc(npoints * sizeof(real));
    // cpumem += npoints * sizeof(real); 

    for(int i = 0; i < npoints; i++) {
      real mass = 4./3.*PI*(points[i].rho-rho_f)*points[i].r*points[i].r*points[i].r;
      x[i] = points[i].x;
      y[i] = points[i].y;
      z[i] = points[i].z;
      conn[i] = npoints-i;
      a[i] = points[i].r;
      u[i] = points[i].u;
      v[i] = points[i].v;
      w[i] = points[i].w;
      udot[i] = points[i].udot;
      vdot[i] = points[i].vdot;
      wdot[i] = points[i].wdot;
    
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
    }

    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateX", x, &Xn);
    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateY", y, &Yn);
    cg_coord_write(fn, bn, zn, RealDouble, "CoordinateZ", z, &Zn);

    cg_section_write(fn, bn, zn, "Elements", NODE, 0, npoints-1, 0, conn, &en);

    cg_sol_write(fn, bn, zn, "Solution", Vertex, &sn);
    cg_field_write(fn, bn, zn, sn, RealDouble, "Radius", a, &fnr);
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
    cg_field_write(fn, bn, zn, sn, RealDouble, "InteractionForceX", iFx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "InteractionForceY", iFy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "InteractionForceZ", iFz, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "TotalForceX", Fx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "TotalForceY", Fy, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "TotalForceZ", Fz, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "MomentX", Lx, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "MomentY", Ly, &fnr);
    cg_field_write(fn, bn, zn, sn, RealDouble, "MomentZ", Lz, &fnr);


    cg_goto(fn, bn, "Zone_t", zn, "end");
    cg_user_data_write("Etc");
    cg_goto(fn, bn, "Zone_t", zn, "Etc", 0, "end");
    cgsize_t *N =(cgsize_t *) malloc(sizeof(cgsize_t));
    N[0] = 1;
    cg_array_write("Time", RealDouble, 1, N, &ttime);
    free(N);

    cg_close(fn);
    free(x);
    free(y);
    free(z);
    free(conn);
    free(a);
  
    free(u);
    free(v);
    free(w);
    free(udot);
    free(vdot);
    free(wdot);
   
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
  }
}


*/





//mask(i1,i2), 如何表示ww[:]?
