//A is located at cell center, Ap is specified property of each particle
//Mollify the particle property Ap, and add it to flow property A
//dirc=0 coresponds to cell-center; dirc=1 x-face-center;dirc=2 y-face-center;dirc=3 z-face-center

__global__ void lpt_mollify_sc(int npoints, real *A, point_struct *points, dom_struct *dom,real *Ap,real *Ksi, int dirc, int dirc2)
{

int indx =  threadIdx.x + blockIdx.x*blockDim.x;
int indy =  threadIdx.y + blockIdx.y*blockDim.y;
int pp=indx+indy*gridDim.x*blockDim.x;

if(pp<npoints)
{
//Define the STENCIL of the gausian filter
real ksi[STENCIL][STENCIL][STENCIL];
real xs=dom->xs;real ys=dom->ys;real zs=dom->zs;
real dx=dom->dx;real dy=dom->dy;real dz=dom->dz;
real cellVol=dx*dy*dz;

  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;
//real cellVol=dx*dy*dz;

real  xp =  points[pp].x;
real  yp =  points[pp].y;
real  zp =  points[pp].z;

real rad=points[pp].r;
real Vp=PI*4/3*rad*rad*rad;

switch(dirc)
{
case 0:
  switch(dirc2)
   { case 0:   
     Ap[pp]=-points[pp].msdot/cellVol;break;
     case 1:
     real rad=  points[pp].r;
     Ap[pp]=PI*4/3*rad*rad*rad/cellVol;break;
     default:printf("Wrong direction in cell center mollify");
	}
break;
case 1:
Ap[pp]=-points[pp].Fx/cellVol;break;
case 2:
Ap[pp]=-points[pp].Fy/cellVol;break;
case 3:
Ap[pp]=-points[pp].Fz/cellVol;break;
default:
  printf("Wrong direction in lpt_mollify");
}

// i,j,k should be the start number of cell center, so that xm(i)<=xp<xm(i+1)
//xm(i)= (i-0.5)*dom->dx +dom->xs
//i=floor((xm - dom->xs) * ddx+ 0.5);
//i=floor((xm - dom->xs) * ddx- 0.5) + DOM_BUF;
  int ip,jp,kp;
  int ic,jc,kc;
  int is,js,ks;

switch(dirc)
{
case 0:
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for xm[i]
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;//for ym[j]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;//for zm[k]
  break;
case 1:
  ip = floor((xp - dom->xs) * ddx + 0.5) + DOM_BUF;//for xm[i], ip should be i+1
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;	   //for y[i]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;	   //for z[i]
  break;
case 2:
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for x[i] 
  jp = floor((yp - dom->ys) * ddy + 0.5) + DOM_BUF;	   //for ym[i]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;	   //for z[i]
  break;
case 3:
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for x[i] 
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;	   //for y[i]
  kp = floor((zp - dom->zs) * ddz+ 0.5) + DOM_BUF;  //for zm[i]
  break;
default:
  printf("Wrong direction in lpt_mollify");
}


real buf=0;
int began=-(STENCIL-1)/2;//equal to -1 if STENCIL=3
int end=1+(STENCIL-1)/2;//equal to 2 if STENCIL=3

for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
{
ksi[di-began][dj-began][dk-began]=lpt_integrate_mol(ip+di,jp+dj,kp+dk,xp,yp,zp,dx,dy,dz,xs,ys,zs,dirc);
buf+=ksi[di-began][dj-began][dk-began];
}

//TODO add mask infomation as in lpt_mollify_sc in lpt_interpolator.f90

// Normalize  ksi = ksi/buf
//if(fabs(buf-1)>EPSILON)
if(buf>0)
{
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
{
ic=ip+di;jc=jp+dj;kc=kp+dk;
is=di-began;js=dj-began;ks=dk-began;
ksi[is][js][ks]=ksi[is][js][ks]/buf;
Ksi[is+js*STENCIL+ks*STENCIL2+pp*npoints*STENCIL3]=ksi[is][js][ks];
}
}


}

}


__global__ void lpt_sum_ksi(int npoints, real *A, point_struct *points, dom_struct *dom,real *Ap,real *Ksi, int dirc)
{
int pp =  threadIdx.x + blockIdx.x*blockDim.x;
if(pp<1)
{

int began=-(STENCIL-1)/2;//equal to -1 if STENCIL=3
int end=1+(STENCIL-1)/2;//equal to 2 if STENCIL=3
real contribution;
for(int np=0;np<npoints;np++)
{
// Perform the actual extrapolation on A
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
{//TODO make sure the scale is correct
//A[(di+ip)+(dj+jp)*dom->Gcc._s1b+(dk+kp)*dom->Gcc._s2b]+=ksi[di-began][dj-began][dk-began]*Ap*pointPartVol/epsp/domVol;
  ic=di+ip;jc=dj+jp;kc=dk+kp;
  is=di-began;js=dj-began;ks=dk-began;
  contribution=Ksi[is+js*STENCIL+ks*STENCIL2+np*npoints*STENCIL3]*Ap[np];
switch(dirc)
{
case 0:
  if(ic < dom->Gcc.is) ic = ic +dom->Gcc.in;
  if(jc < dom->Gcc.js) jc = jc +dom->Gcc.jn;
  if(kc < dom->Gcc.ks) kc = kc +dom->Gcc.kn;
  if(ic > dom->Gcc.ie-1) ic = ic -dom->Gcc.in;
  if(jc > dom->Gcc.je-1) jc = jc -dom->Gcc.jn;
  if(kc > dom->Gcc.ke-1) kc = kc -dom->Gcc.kn;
A[ic+jc*dom->Gcc._s1b+kc*dom->Gcc._s2b]+=contribution;
break;
case 1:
  if(ic < dom->Gfx.is) ic = ic+(dom->Gfx.in-1);
  if(jc < dom->Gfx.js) jc = jc+dom->Gfx.jn;
  if(kc < dom->Gfx.ks) kc = kc+dom->Gfx.kn;
  if(ic > dom->Gfx.ie-1) ic = ic-(dom->Gfx.in-1);
  if(jc > dom->Gfx.je-1) jc = jc-dom->Gfx.jn;
  if(kc > dom->Gfx.ke-1) kc = kc-dom->Gfx.kn;
A[ic+jc*dom->Gfx._s1b+kc*dom->Gfx._s2b]+=contribution;
break;
case 2:
  if(ic < dom->Gfy.is) ic = ic+dom->Gfy.in;
  if(jc < dom->Gfy.js) jc = jc+(dom->Gfy.jn-1);
  if(kc < dom->Gfy.ks) kc = kc+dom->Gfy.kn;
  if(ic > dom->Gfy.ie-1) ic = ic-dom->Gfy.in;
  if(jc > dom->Gfy.je-1) jc = jc-(dom->Gfy.jn-1);
  if(kc > dom->Gfy.ke-1) kc = kc-dom->Gfy.kn;
A[ic+jc*dom->Gfy._s1b+kc*dom->Gfy._s2b]+=contribution;
break;
case 3
  if(ic < dom->Gfz.is) ic = ic+dom->Gfz.in;
  if(jc < dom->Gfz.js) jc = jc+dom->Gfz.jn;
  if(kc < dom->Gfz.ks) kc = kc+(dom->Gfz.kn-1);
  if(ic > dom->Gfz.ie-1) ic = ic-dom->Gfz.in;
  if(jc > dom->Gfz.je-1) jc = jc-dom->Gfz.jn;
  if(kc > dom->Gfz.ke-1) kc = kc-(dom->Gfz.kn-1);
A[ic+jc*dom->Gfz._s1b+kc*dom->Gfz._s2b]+=contribution;
break;
default:
break;
			}		

		}
	}
  }
}



__device__ periodic_grid_index(int ic,int jc,int kc,real *dom, int dirc)
{
switch(dirc)
{
case 0:
  if(ic < dom->Gcc.is) ic = ic +dom->Gcc.in;
  if(jc < dom->Gcc.js) jc = jc +dom->Gcc.jn;
  if(kc < dom->Gcc.ks) kc = kc +dom->Gcc.kn;
  if(ic > dom->Gcc.ie-1) ic = ic -dom->Gcc.in;
  if(jc > dom->Gcc.je-1) jc = jc -dom->Gcc.jn;
  if(kc > dom->Gcc.ke-1) kc = kc -dom->Gcc.kn;
break;
case 1:
  if(ic < dom->Gfx.is) ic = ic+(dom->Gfx.in-1);
  if(jc < dom->Gfx.js) jc = jc+dom->Gfx.jn;
  if(kc < dom->Gfx.ks) kc = kc+dom->Gfx.kn;
  if(ic > dom->Gfx.ie-1) ic = ic-(dom->Gfx.in-1);
  if(jc > dom->Gfx.je-1) jc = jc-dom->Gfx.jn;
  if(kc > dom->Gfx.ke-1) kc = kc-dom->Gfx.kn;
break;
case 2:
  if(ic < dom->Gfy.is) ic = ic+dom->Gfy.in;
  if(jc < dom->Gfy.js) jc = jc+(dom->Gfy.jn-1);
  if(kc < dom->Gfy.ks) kc = kc+dom->Gfy.kn;
  if(ic > dom->Gfy.ie-1) ic = ic-dom->Gfy.in;
  if(jc > dom->Gfy.je-1) jc = jc-(dom->Gfy.jn-1);
  if(kc > dom->Gfy.ke-1) kc = kc-dom->Gfy.kn;
break;
case 3:
  if(ic < dom->Gfz.is) ic = ic+dom->Gfz.in;
  if(jc < dom->Gfz.js) jc = jc+dom->Gfz.jn;
  if(kc < dom->Gfz.ks) kc = kc+(dom->Gfz.kn-1);
  if(ic > dom->Gfz.ie-1) ic = ic-dom->Gfz.in;
  if(jc > dom->Gfz.je-1) jc = jc-dom->Gfz.jn;
  if(kc > dom->Gfz.ke-1) kc = kc-(dom->Gfz.kn-1);
break;
default:
break;
}

}


__device__ real lpt_integrate_mol(int ic,int jc,int kc,real xp,real yp,real zp, real dx,real dy,real dz,real xs,real ys,real zs, int dirc)
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

switch(dirc)
{
case 0:
real xm = (ic-DOM_BUF+0.5) * dx + xs;
real ym = (jc-DOM_BUF+0.5) * dy + ys;
real zm = (kc-DOM_BUF+0.5) * dz + zs;
real r2 = (xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm);
break;
case 1:
real x =  (ic-DOM_BUF) * dx + xs;
real ym = (jc-DOM_BUF+0.5) * dy + ys;
real zm = (kc-DOM_BUF+0.5) * dz + zs;
real r2 = (xp-x)*(xp-x)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm);
break;
case 2:
real xm = (ic-DOM_BUF+0.5) * dx + xs;
real y = (jc-DOM_BUF) * dy + ys;
real zm = (kc-DOM_BUF+0.5) * dz + zs;
real r2 = (xp-xm)*(xp-xm)+(yp-y)*(yp-y)+(zp-zm)*(zp-zm);
break;
case 3:
real xm = (ic-DOM_BUF+0.5) * dx + xs;
real ym = (jc-DOM_BUF+0.5) * dy + ys;
real z = (kc-DOM_BUF) * dz + zs;
real r2 = (xp-xm)*(xp-xm)+(yp-ym)*(yp-ym)+(zp-z)*(zp-z);
break;
default:
printf("Wrong dirc in lpt_integrate_mol");
}

//real r2 = (xp-xm)*(xp-xm)+(yp-ym)*(yp-ym);
//TODO make this as defined value avaible from host and device!!!
real min_meshsize=min(min(dx,dy),dz);
real cellVol=dx *dy *dz;
real sig= KERNEL_WIDTH *min_meshsize/(2.0f*sqrt(2.0f*log(2.0f)));

real val = exp(-r2/(2.0f*sig*sig));
return val;
}

