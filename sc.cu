extern "C"
void cuda_points_locate()
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

lpt_localize<<<numBlocks, dimBlocks>>>(npoints,_points[dev],_dom[dev],bc);	
	}
}


__global__ void lpt_localize(int npoints, point_struct *points, dom_struct *dom, BC bc)
{
 int pp =  threadIdx.x + blockIdx.x*blockDim.x;
if(pp<npoints)
{
// Cartesian location of node
real  xp =  points[pp].x;
real  yp =  points[pp].y;
real  zp =  points[pp].z;

//TODO whether periodic BC for point particle need to be determined in future
  if(xp < dom->xs && bc.uW == PERIODIC) xp = xp + dom->xl;
  else if(xp > dom->xe && bc.uE == PERIODIC) xp = xp - dom->xl;
  if(yp < dom->ys && bc.vS == PERIODIC) yp = yp + dom->yl;
  else if(yp > dom->ye && bc.vN == PERIODIC) yp = yp - dom->yl;
  if(zp < dom->zs && bc.wB == PERIODIC) zp = zp + dom->zl;
  else if(zp > dom->ze && bc.wT == PERIODIC) zp = zp - dom->zl;


int i =  points[pp].i;
int j =  points[pp].j;
int k =  points[pp].k;

//TODO
  if(i < dom->Gcc.is) i = dom->Gcc.ie-1;
  if(j < dom->Gcc.js) j = dom->Gcc.je-1;
  if(k < dom->Gcc.ks) k = dom->Gcc.ke-1;
  if(i > dom->Gcc.ie-1) i = dom->Gcc.is;
  if(j > dom->Gcc.je-1) j = dom->Gcc.js;
  if(k > dom->Gcc.ke-1) k = dom->Gcc.ks;

int ip=i;
int jp=j;
int kp=k;

real  x = (ip-DOM_BUF) * dom->dx + dom->xs;
real  y = (jp-DOM_BUF) * dom->dy + dom->ys;
real  z = (kp-DOM_BUF) * dom->dz + dom->zs;

while(xp <x)
{
ip=ip-1;
x = (ip-DOM_BUF) * dom->dx + dom->xs;
}

while(xp >=x)
{
ip=ip+1;
x = (ip-DOM_BUF) * dom->dx + dom->xs;
}

while(yp <y)
{
jp=jp-1;
y = (jp-DOM_BUF) * dom->dy + dom->ys;
}
while(yp >=y)
{
jp=jp+1;
y = (jp-DOM_BUF) * dom->dy + dom->ys;
}

while(zp <z)
{
kp=kp-1;
z = (kp-DOM_BUF) * dom->dz + dom->zs;
}

while(zp >=z)
{
kp=kp+1;
z = (kp-DOM_BUF) * dom->dz + dom->zs;
}

points[pp].i= ip;
points[pp].j= jp;
points[pp].k= kp;
}

}








__global__ void lpt_smooth(int npoints, point_struct *points, dom_struct *dom, real *fpx,real *fpy,real *fpz,real *f_x,real *f_y,real *f_z, BC bc)
{
 int pp =  threadIdx.x + blockIdx.x*blockDim.x;
  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;

if(pp<npoints)
{
// Cartesian location of node
real  xp =  points[pp].x;
real  yp =  points[pp].y;
real  zp =  points[pp].z;

//TODO whether periodic BC for point particle need to be determined in future
  if(xp < dom->xs && bc.uW == PERIODIC) xp = xp + dom->xl;
  else if(xp > dom->xe && bc.uE == PERIODIC) xp = xp - dom->xl;
  if(yp < dom->ys && bc.vS == PERIODIC) yp = yp + dom->yl;
  else if(yp > dom->ye && bc.vN == PERIODIC) yp = yp - dom->yl;
  if(zp < dom->zs && bc.wB == PERIODIC) zp = zp + dom->zl;
  else if(zp > dom->ze && bc.wT == PERIODIC) zp = zp - dom->zl;

int i =  points[pp].i;
int j =  points[pp].j;
int k =  points[pp].k;

  if(i < dom->Gcc.is) i = dom->Gcc.ie-1;
  if(j < dom->Gcc.js) j = dom->Gcc.je-1;
  if(k < dom->Gcc.ks) k = dom->Gcc.ke-1;
  if(i > dom->Gcc.ie-1) i = dom->Gcc.is;
  if(j > dom->Gcc.je-1) j = dom->Gcc.js;
  if(k > dom->Gcc.ke-1) k = dom->Gcc.ks;


//TODO how to smooth the func????
fx[pp]=
fy[pp]=
fz[pp]= 

 	
  }


}



__device__ void lpt_source_scalar(real xd,real yd,real zd,int id,int jd,int kd,real ud,real dud,real vd,real dvd,real wd,real dwd,real md,real dmd,int nparceld)
{
// Create the source term
//TODO need to change scalar source!!!  
real  scalar_src =1.f
//mollify each source term
lpt_mollify_sc(tmp1,xd,yd,zd,id,jd,kd,scalar_src);
}


__global__ void lpt_source_scalar(real *u, real *v, real *w,int npoints,real rho_f, real nu,real *ug,real *vg,real *wg, point_struct *points, dom_struct *dom, real dt0, real dt, BC bc)

real *scSrc,


{
  int pp =  threadIdx.x + blockIdx.x*blockDim.x;
  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;


if(pp<npoints)
{
// Cartesian location of node
real  xd =  points[pp].x;
real  yd =  points[pp].y;
real  zd =  points[pp].z;

int id =  points[pp].i;
int jd =  points[pp].j;
int kd =  points[pp].k;

// Create the source term
//TODO need to change scalar source!!!  
real  scalar_src_points =1.f
//mollify each source term
lpt_mollify_sc(scSrc,xd,yd,zd,id,jd,kd,scalar_src_points);

//
entrySearch_avg_entries_kernel(real *iarr, real *maxarr, int size);
  }
}




/*==================================== !
! Spray -> gas phase momentum exchange !
! Careful, still need to divide by vol !
! ==================================== !*/

__device__ void lpt_source_momentum(real xd,real yd,real zd,int id,int jd,int kd,real ud,real dud,real vd,real dvd,real wd,real dwd,real md,real dmd,int nparceld)
{
  // Create the source term
real  mom_src_x = -(md*dud+dmd*ud)*nparceld
real  mom_src_y = -(md*dvd+dmd*vd)*nparceld
real  mom_src_z = -(md*dwd+dmd*wd)*nparceld
lpt_mollify_sc(tmp1,xd,yd,zd,id,jd,kd,mom_src_x);
lpt_mollify_sc(tmp2,xd,yd,zd,id,jd,kd,mom_src_y);
lpt_mollify_sc(tmp3,xd,yd,zd,id,jd,kd,mom_src_z);
}











//Case periodic in y and z directions
// xp~zp and ip~jp are the particle position and the cell index of the particle correspondingly!  
//A is the value on grid cell center(which is need to be interpretated), and Ap is the particle contribution strength
__device__ void lpt_mollify_sc(real *A,real xp,real yp,real zp,int ip,int jp,int kp,real Ap)
{
  if(ip < dom->Gcc.is||jp < dom->Gcc.js||kp < dom->Gcc.ks||ip > dom->Gcc.ie-1||jp > dom->Gcc.je-1||kp > dom->Gcc.ke-1)
{
printf("\nip,jp,kp,xp,yp,zp %d %d %d %d %d %d\n",ip,jp,kp,xp,yp,zp);
fprintf(stderr,"\nParticle has left the domain\n");
exit(EXIT_FAILURE);
}

real ksi[3][3][3];
real buf=0;
for(int dk=-1;dk<2;dk++)
for(int dj=-1;dj<2;dj++)
for(int di=-1;di<2;di++)
{{{
//TODO what's ksi????
ksi[di+1][dj+1][dk+1]=lpt_integrate_mol(ip+di,jp+dj,kp+dk,xp,yp,zp);
buf+=ksi[di+1][dj+1][dk+1];
}}}

//TODO add mask infomation as in lpt_mollify_sc in lpt_interpolator.f90

// Normalize  ksi = ksi/buf
if (buf>0.f){
for(int dk=-1;dk<2;dk++)
for(int dj=-1;dj<2;dj++)
for(int di=-1;di<2;di++)
{{{ksi[di+1][dj+1][dk+1]=ksi[di+1][dj+1][dk+1]/buf;}}}
}

// Perform the actual extrapolation on A
for(int dk=-1;dk<2;dk++)
for(int dj=-1;dj<2;dj++)
for(int di=-1;di<2;di++)
{{{
A[di+ip][dj+jp][dk+kp]+=ksi[di+1][dj+1][dk+1]*Ap;
}}}

}


__device__ real lpt_integrate_mol(int ic,int jc,int kc,real xp,real yp,real zp, real dx,real dy,real dz,real xs,real ys,real zs)
{
/*
real xm = (ic-DOM_BUF+0.5) * dom->dx + dom->xs;
real ym = (jc-DOM_BUF+0.5) * dom->dy + dom->ys;
real zm = (kc-DOM_BUF+0.5) * dom->dz + dom->zs;
real r = sqrt((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
//TODO make this as defined value avaible from host and device!!!
real min_meshsize=min(min(dom->dx,dom->dy),dom->dz);
real cellVol=dom->dx *dom->dy *dom->dz;
*/
real xm = (ic-DOM_BUF+0.5) * dx + xs;
real ym = (jc-DOM_BUF+0.5) * dy + ys;
real zm = (kc-DOM_BUF+0.5) * dz + zs;
real r = sqrt((xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm));
//TODO make this as defined value avaible from host and device!!!
real min_meshsize=min(min(dx,dy),dz);
real cellVol=dx *dy *dz;

real sig= min_meshsize/(2.0f*sqrt(2.0f*log(2.0f)));
real val = exp(-r*r/(2.0f*sig*sig));
real fs= cellVol*val; 
return fs;
}



__global__ void forcing_add_sc_const(real val, real *sc, dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  for(int i = dom->Gcc._isb; i < dom->Gcc._ieb; i++) {
    if(tj < dom->Gcc._jnb && tk < dom->Gcc._knb) {
      sc[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b] += val;
    }
  }
}


__global__ void forcing_add_sc_field(real scale, real *val, real *sc,
  dom_struct *dom)
{
  int tj = blockIdx.x * blockDim.x + threadIdx.x;
  int tk = blockIdx.y * blockDim.y + threadIdx.y;

  for(int i = dom->Gcc._isb; i < dom->Gcc._ieb; i++) {
    if(tj < dom->Gcc._jnb && tk < dom->Gcc._knb) {
      sc[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        += scale * val[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b];
    }
  }
}

