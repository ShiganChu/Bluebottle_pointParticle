__global__ void lpt_mollify_momentum_z(int npoints,real *epsp, real *fz, point_struct *points, dom_struct *dom, real dt)
{

int index =  threadIdx.x + blockIdx.x*blockDim.x;

if(index<1)
{
//Define the stencil of the gausian filter
real ksi[STENCIL][STENCIL][STENCIL];
real xs=dom->xs;real ys=dom->ys;real zs=dom->zs;
real dx=dom->dx;real dy=dom->dy;real dz=dom->dz;
  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;
//real cellVol=dx*dy*dz;
for(int pp=0;pp<npoints;pp++)
{

real  xp =  points[pp].x;
real  yp =  points[pp].y;
real  zp =  points[pp].z;

real rad=points[pp].r;
real Vp=PI*4/3*rad*rad*rad;
real mp=Vp*points[pp].rho + points[pp].ms;
//ms*w+m*du*dt
real srcW=-(points[pp].msdot*points[pp].w+points[pp].wdot*dt*mp);

  int ip,jp,kp;
  // interpolate w-momentum source
  ip = round((x - dom->xs) * ddx) + DOM_BUF;
  jp = floor((y - dom->ys) * ddy) + DOM_BUF;
  kp = floor((z - dom->zs) * ddz - 0.5) + DOM_BUF;

real buf=0;
int began=-(STENCIL-1)/2;//equal to -1 if STENCIL=3
int end=1+(STENCIL-1)/2;//equal to 2 if STENCIL=3

for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
{
ksi[di-began][dj-began][dk-began]=lpt_integrate_mol_w(ip+di,jp+dj,kp+dk,xp,yp,zp,dx,dy,dz,xs,ys,zs);
buf+=ksi[di-began][dj-began][dk-began];
}


// Normalize  ksi = ksi/buf
//if(fabs(buf-1)>EPSILON)
if(buf>0)
{
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
{ksi[di-began][dj-began][dk-began]=ksi[di-began][dj-began][dk-began]/buf;}
}


int ic,jc,kc;
// Perform the actual extrapolation on A
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
{
//TODO make sure the scale is correct
//A[(di+ip)+(dj+jp)*dom->Gfz._s1b+(dk+kp)*dom->Gfz._s2b]+=ksi[di-began][dj-began][dk-began]*Ap*pointPartVol/epsp/domVol;
  ic=di+ip;jc=dj+jp;kc=dk+kp;
  if(ic < dom->Gfz.is) ic = ic+(dom->Gfz.in-1);
  if(jc < dom->Gfz.js) jc = jc+dom->Gfz.jn;
  if(kc < dom->Gfz.ks) kc = kc+dom->Gfz.kn;
  if(ic > dom->Gfz.ie-1) ic = ic-(dom->Gfz.in-1);
  if(jc > dom->Gfz.je-1) jc = jc-dom->Gfz.jn;
  if(kc > dom->Gfz.ke-1) kc = kc-dom->Gfz.kn;
fz[ic+jc*dom->Gfz._s1b+kc*dom->Gfz._s2b]+=ksi[di-began][dj-began][dk-began]*srcW;
}

}
}
}

