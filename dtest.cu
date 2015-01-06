__device__ real lpt_integrate_mol_u(int ic,int jc,int kc,real xp,real yp,real zp, real dx,real dy,real dz,real xs,real ys,real zs)
{

real x = (ic-DOM_BUF) * dx + xs;
real ym = (jc-DOM_BUF+0.5) * dy + ys;
real zm = (kc-DOM_BUF+0.5) * dz + zs;
real r2 = (xp-x)*(xp-x)+(yp-ym)*(yp-ym)+(zp-zm)*(zp-zm);
//TODO make this as defined value avaible from host and device!!!
real min_meshsize=min(min(dx,dy),dz);
//real cellVol=dx *dy *dz;
real sig= min_meshsize/(2.0f*sqrt(2.0f*log(2.0f)));

real val = exp(-r2/(2.0f*sig*sig));
//real fs= cellVol*val; 
//return fs;
return val;
}

__device__ real lpt_integrate_mol_v(int ic,int jc,int kc,real xp,real yp,real zp, real dx,real dy,real dz,real xs,real ys,real zs)
{

real xm = (ic-DOM_BUF+0.5) * dx + xs;
real y = (jc-DOM_BUF) * dy + ys;
real zm = (kc-DOM_BUF+0.5) * dz + zs;
real r2 = (xp-xm)*(xp-xm)+(yp-y)*(yp-y)+(zp-zm)*(zp-zm);
//TODO make this as defined value avaible from host and device!!!
real min_meshsize=min(min(dx,dy),dz);
//real cellVol=dx *dy *dz;
real sig= min_meshsize/(2.0f*sqrt(2.0f*log(2.0f)));

real val = exp(-r2/(2.0f*sig*sig));
//real fs= cellVol*val; 
//return fs;
return val;
}

__device__ real lpt_integrate_mol_w(int ic,int jc,int kc,real xp,real yp,real zp, real dx,real dy,real dz,real xs,real ys,real zs)
{

real xm = (ic-DOM_BUF+0.5) * dx + xs;
real y = (jc-DOM_BUF+0.5) * dy + ys;
real z = (kc-DOM_BUF) * dz + zs;
real r2 = (xp-xm)*(xp-xm)+(yp-y)*(yp-y)+(zp-z)*(zp-z);
//TODO make this as defined value avaible from host and device!!!
real min_meshsize=min(min(dx,dy),dz);
//real cellVol=dx *dy *dz;
real sig= min_meshsize/(2.0f*sqrt(2.0f*log(2.0f)));

real val = exp(-r2/(2.0f*sig*sig));
//real fs= cellVol*val; 
//return fs;
return val;
}

