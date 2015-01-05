/*
A is fluid property located at cell center or face center depending on dirc 
dirc=0 coresponds to cell-center; dirc=1 x-face-center;dirc=2 y-face-center;dirc=3 z-face-center
dir
*/
__global__ void lpt_mollify_sc(int npoints, real *A, point_struct *points, dom_struct *dom,real *Ksi, int dirc, int dirc2)
{

int indx =  threadIdx.x + blockIdx.x*blockDim.x;
int indy =  threadIdx.y + blockIdx.y*blockDim.y;
int pp=indx+indy*gridDim.x*blockDim.x;

//if(pp<npoints*STENCIL3)  Ksi[pp]=0.0;
//__sycthreads();

if(pp<npoints)
{

////printf("\npp indx %d %d %d \n",pp, indx, indy);
////printf("\nblockDim.x,gridDim.x %d %d %d %d\n",blockDim.x,gridDim.x,blockDim.y,gridDim.y );

//Define the STENCIL of the Gausian filter, the spread length of Gaussian kernel
real ksi[STENCIL][STENCIL][STENCIL];
real xs=dom->xs;real ys=dom->ys;real zs=dom->zs;
real dx=dom->dx;real dy=dom->dy;real dz=dom->dz;
  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;
//real cellVol=dx*dy*dz;

real  xp =  points[pp].x;
real  yp =  points[pp].y;
real  zp =  points[pp].z;

// i,j,k should be the start number of cell center, so that xm(i)<=xp<xm(i+1)
//xm(i)= (i-0.5)*dom->dx +dom->xs
//i=floor((xm - dom->xs) * ddx+ 0.5);
//i=floor((xm - dom->xs) * ddx- 0.5) + DOM_BUF;
  int ip,jp,kp;
  int is,js,ks;

//The cell-center or face-center index to locate the particle 
switch(dirc)
{
case 0:
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for xm[i]
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;//for ym[j]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;//for zm[k]
  break;
case 1:
  ip = floor((xp - dom->xs) * ddx + 0.5) + DOM_BUF;//for xm[i+1], ip should be i+1
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;	   //for y[i]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;	   //for z[i]
  break;
case 2:
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for x[i] 
  jp = floor((yp - dom->ys) * ddy + 0.5) + DOM_BUF;	   //for ym[i+1]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;	   //for z[i]
  break;
case 3:
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for x[i] 
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;	   //for y[i]
  kp = floor((zp - dom->zs) * ddz+ 0.5) + DOM_BUF;  //for zm[i+1]
  break;
default:
  printf("Wrong direction in lpt_mollify");
}

points[pp].i=ip;
points[pp].j=jp;
points[pp].k=kp;

real buf=0;
int began=-(STENCIL-1)/2;//equal to -1 if STENCIL=3
int end=1+(STENCIL-1)/2;//equal to 2 if STENCIL=3


//Calculate the filter strength from Gaussian kernel
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
{
is=di-began;js=dj-began;ks=dk-began;
ksi[is][js][ks]=lpt_integrate_mol(ip+di,jp+dj,kp+dk,xp,yp,zp,dx,dy,dz,xs,ys,zs,dirc);
buf+=ksi[is][js][ks];
}

//TODO add mask infomation as in lpt_mollify_sc in lpt_interpolator.f90

real Ap;
switch(dirc)
{
case 0:
    switch(dirc2)
   { case 0:   
//Exchange rate of soluble mass from particles. All the following source should divide cellVol to ensure conservation law
     Ap=-points[pp].msdot/cellVol;break;
     case 1:
//Particle volume filter
     Ap=PI*4/3*points[pp].r*points[pp].r*points[pp].r/cellVol;break;
     default: Ap=0;break;
	}
break;
case 1:
//Particle x-momentum reaction to fluid
Ap=-points[pp].Fx/cellVol;
break;
case 2:
Ap=-points[pp].Fy/cellVol;
break;
case 3:
//Particle z-momentum reaction to fluid
Ap=-points[pp].Fz/cellVol;
break;
default:
break;
	}		


// Normalize  ksi = ksi/buf
if(buf>0)
	{
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
							{
is=di-began;js=dj-began;ks=dk-began;
ksi[is][js][ks]=ksi[is][js][ks]/buf;
Ksi[is+js*STENCIL+ks*STENCIL2+pp*STENCIL3]=ksi[is][js][ks]*Ap; }
		}


	}

}


//Add the point source to fluid source weighted by the Gaussian filter coefficient Ksi
__global__ void lpt_sum_ksi(int npoints, real *A, point_struct *points, dom_struct *dom,real *Ksi, int dirc,int dirc2)
{
int indx =  threadIdx.x + blockIdx.x*blockDim.x;
int indy =  threadIdx.y + blockIdx.y*blockDim.y;
int np=indx+indy*gridDim.x*blockDim.x;
if(np<1)
{

int began=-(STENCIL-1)/2;//equal to -1 if STENCIL=3
int end=1+(STENCIL-1)/2;//equal to 2 if STENCIL=3
  int ip,jp,kp;
  int is,js,ks;
  int ic,jc,kc;
real Ap,ksi;
real cellVol=dom->dx*dom->dy*dom->dz;

for(int pp=0;pp<npoints;pp++)
{
ip =  points[pp].i;
jp =  points[pp].j;
kp =  points[pp].k;

// Perform the actual extrapolation on A
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
{
  ic=di+ip;jc=dj+jp;kc=dk+kp;
  is=di-began;js=dj-began;ks=dk-began;
  ksi=Ksi[is+js*STENCIL+ks*STENCIL2+pp*STENCIL3];
switch(dirc)
{
case 0:
//Make the index inside the domain before we change the property at this location
periodic_grid_index(ic,jc,kc,dom,0);
A[ic+jc*dom->Gcc._s1b+kc*dom->Gcc._s2b]+=ksi;
break;
case 1:
periodic_grid_index(ic,jc,kc,dom,1);
A[ic+jc*dom->Gfx._s1b+kc*dom->Gfx._s2b]+=ksi;
break;
case 2:
periodic_grid_index(ic,jc,kc,dom,2);
A[ic+jc*dom->Gfy._s1b+kc*dom->Gfy._s2b]+=ksi;
break;
case 3:
periodic_grid_index(ic,jc,kc,dom,3);
A[ic+jc*dom->Gfz._s1b+kc*dom->Gfz._s2b]+=ksi;
break;
default:
break;
			}		

		}
	}
  }
}


