
__global__ void advance_sc_upwind_2nd(real DIFF, real *u, real *v, real *w,real *f,real *epsp,
  real *diff0, real *conv0, real *diff, real *conv, real *sc,real *sc0,
  dom_struct *dom, real dt0, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward
  __shared__ real s_wi[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back back
  __shared__ real s_w2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward 2
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
//store scalar related values such as convective and diffusive term.
  __shared__ real s_d0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // diff0
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // conv0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // source term


//store scalar value of difference direction at cell center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top	 center
  __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at  	 center

//store point volume percentage in each grid cell
  __shared__ real p_epsp[MAX_THREADS_DIM * MAX_THREADS_DIM];       
 

  // working constants
  real ab0 = 0.5 * dt / dt0;  // for Adams-Bashforth stepping
  real ab = 1 + ab0;          // for Adams-Bashforth stepping

  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

 
//To calculate the peclet number
  real ui= dom->dx/DIFF;
  real vi= dom->dy/DIFF;
  real wi= dom->dz/DIFF;
 

  int C;

  // loop over z-planes
//  for(int k = dom->Gcc._ksb; k < dom->Gcc._keb; k++) {
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int i = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;

    // shared memory indices
    int ti = threadIdx.x;
    int tj = threadIdx.y;

    // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations

    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];

s_c0[ti + tj*blockDim.x]=conv0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
s_d0[ti + tj*blockDim.x]=diff0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
 s_f[ti + tj*blockDim.x]=    f[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];

p_epsp[ti + tj*blockDim.x]= epsp[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
}

   // make sure all threads complete shared memory copy
    __syncthreads();

 
    if((j >= dom->Gfz._jsb && j < dom->Gfz._jeb)
      && (i >= dom->Gfz._isb && i < dom->Gfz._ieb)) {
      s_wi[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k-1)*dom->Gfz._s2b];
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
      s_w2[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
   }
    if((j >= dom->Gfx._jsb && j < dom->Gfx._jeb)
      && (i >= dom->Gfx._isb && i < dom->Gfx._ieb)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._jsb && j < dom->Gfy._jeb)
      && (i >= dom->Gfy._isb && i < dom->Gfy._ieb)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
 
 

    //right hand side for scalar term
 //   s_rhs[ti + tj*blockDim.x] = 0.0;
    //convective term
 //   s_c[ti + tj*blockDim.x] = 0.0;
    //diffusion term
 //   s_d[ti + tj*blockDim.x] = 0.0;

    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute convective term
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) 
 
//if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) 
{
 
//scalar on the west and east u face
real u_c=(s_u[ti + tj*blockDim.x]+s_u[(ti+1) + tj*blockDim.x])/2;// at xm(i)
real v_c=(s_v[ti + tj*blockDim.x]+s_v[ti + (tj+1)*blockDim.x])/2;
real w_c=(s_w1[ti + tj*blockDim.x]+s_w0[ti + tj*blockDim.x])/2;

real u_ce=(s_u[(ti+2) + tj*blockDim.x]+s_u[(ti+1) + tj*blockDim.x])/2;//at xm(i+1)
real u_cw=(s_u[(ti-1) + tj*blockDim.x]+s_u[ti + tj*blockDim.x])/2;

real v_cn=(s_v[ti + (tj+2)*blockDim.x]+s_v[ti + (tj+1)*blockDim.x])/2;
real v_cs=(s_v[ti + (tj-1)*blockDim.x]+s_v[ti + tj*blockDim.x])/2;

real w_ct=(s_w1[ti + tj*blockDim.x]+s_w2[ti + tj*blockDim.x])/2;
real w_cb=(s_wi[ti + tj*blockDim.x]+s_w0[ti + tj*blockDim.x])/2;

real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0;//at x(i)
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;

real Pe_u= u_c*ui;
real Pe_v= v_c*vi;
real Pe_w= w_c*wi;

real dsudx_up1=( 3*sc_c[ti + tj*blockDim.x]*u_c-4*sc_uw*s_u[ti + tj*blockDim.x]+ sc_c[ti-1 + tj*blockDim.x]*u_cw)*ddx;//backward
real dsudx_up2=(-3*sc_c[ti + tj*blockDim.x]*u_c+4*sc_ue*s_u[ti +1 + tj*blockDim.x] - sc_c[ti+1 + tj*blockDim.x]*u_ce)*ddx;//forward

real dsvdy_up1=( 3*sc_c[ti + tj*blockDim.x]*v_c-4*sc_vs*s_v[ti + tj*blockDim.x]+ sc_c[ti + (tj-1)*blockDim.x]*v_cs)*ddy;//backward
real dsvdy_up2=(-3*sc_c[ti + tj*blockDim.x]*v_c+4*sc_vn*s_v[ti + (tj+1)*blockDim.x]-sc_c[ti+(tj+1)*blockDim.x]*v_cn)*ddy;//forward

real dswdz_up1=( 3*sc_c[ti + tj*blockDim.x]*w_c-4*sc_wb*s_w0[ti + tj*blockDim.x]+sc_b[ti +tj*blockDim.x]*w_cb)*ddz;//backward
real dswdz_up2=(-3*sc_c[ti + tj*blockDim.x]*w_c+4*sc_wt*s_w1[ti + tj*blockDim.x]-sc_t[ti+tj*blockDim.x]*w_ct)*ddz;//forward


//convective term of scalar at cell center -- central difference
real dsudx_ce=(s_u[(ti+1) + tj*blockDim.x]*sc_ue - s_u[ti + tj*blockDim.x]*sc_uw) * ddx;//at xm(i)
real dsvdy_ce=(s_v[ti + (tj+1)*blockDim.x]*sc_vn - s_v[ti + tj*blockDim.x]*sc_vs) * ddy;
real dswdz_ce=(s_w1[ti + tj*blockDim.x]   *sc_wt - s_w0[ti + tj*blockDim.x]*sc_wb) * ddz;

if(Pe_u>2.0) dsudx=dsudx_up1;
else if(Pe_u<-2.0) dsudx=dsudx_up2;
else dsudx=dsudx_ce;

if(Pe_v>2.0) dsvdy=dsvdy_up1;
else if(Pe_v<-2.0) dsvdy=dsvdy_up2;
else dsvdy=dsvdy_ce;

if(Pe_w>2.0) dswdz=dswdz_up1;
else if(Pe_w<-2.0) dswdz=dswdz_up2;
else dswdz=dswdz_ce;


//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;

     //move convective term to right hand side 
//real sc_conv=-(dsudx+dsvdy+dswdz);
      s_c[ti + tj*blockDim.x]=-(dsudx+dsvdy+dswdz);
     //Diffusion term
//real sc_diff=(ddsdxx+ddsdyy+ddsdzz)*DIFF;
      s_d[ti + tj*blockDim.x]= (ddsdxx+ddsdyy+ddsdzz)*DIFF;

//adam-bashforth at t+1/2 step
real sc_conv= ab * s_c[ti + tj*blockDim.x]-ab0*s_c0[ti + tj*blockDim.x];
real sc_diff= ab * s_d[ti + tj*blockDim.x]-ab0*s_d0[ti + tj*blockDim.x];

//TODO should we have s_f at t+1/2 step as well???
//right hand side of scalar equation
//s_rhs[ti + tj*blockDim.x]=sc_conv+sc_diff+s_f[ti + tj*blockDim.x];
//advance scalar
//Take the particle volume fraction into consideration
//real rhs=(sc_diff+s_f[ti + tj*blockDim.x])/(1-p_epsp[ti + tj*blockDim.x]);
real rhs=sc_diff+s_f[ti + tj*blockDim.x];
sc_c[ti + tj*blockDim.x] +=(sc_conv+rhs)*dt;
}

    // copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
     && (i >= dom->Gcc._is && i < dom->Gcc._ie)
  //  if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)
  //    && (i >= dom->Gcc._isb && i < dom->Gcc._ieb)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      sc[C] = sc_c[ti + tj*blockDim.x];
      conv[C] = s_c[ti + tj*blockDim.x];
      diff[C] = s_d[ti + tj*blockDim.x];
    }
  }
}


//Using central difference scheme in space, and adam-bashforth scheme in time.
__global__ void advance_sc_upwind_1st(	real DIFF, 
				real *u, real *v, real *w,
				real *f,real *epsp,
				real *diff0, real *conv0, 
				real *diff, real *conv, 
				real *sc,real *sc0,
				dom_struct *dom, real dt0, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
//store scalar related values such as convective and diffusive term.
  __shared__ real s_d0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // diff0
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // conv0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // source term


//store scalar value of difference direction at cell center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top	 center
  __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at  	 center

//store point volume percentage in each grid cell
  __shared__ real p_epsp[MAX_THREADS_DIM * MAX_THREADS_DIM];       
 

  // working constants
  real ab0 = 0.5 * dt / dt0;  // for Adams-Bashforth stepping
  real ab = 1 + ab0;          // for Adams-Bashforth stepping

  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

 
  int C;

  // loop over z-planes
//  for(int k = dom->Gcc._ksb; k < dom->Gcc._keb; k++) {
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int i = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;

    // shared memory indices
    int ti = threadIdx.x;
    int tj = threadIdx.y;

    // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations

    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];

s_c0[ti + tj*blockDim.x]=conv0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
s_d0[ti + tj*blockDim.x]=diff0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
 s_f[ti + tj*blockDim.x]=    f[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];

p_epsp[ti + tj*blockDim.x]= epsp[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];

//if(isnan(sc_c[ti + tj*blockDim.x])) printf("\nglobal_var %f %f %f %f %d %d %d %d %d\n",sc_b[ti + tj*blockDim.x],sc_t[ti + tj*blockDim.x],s_f[ti + tj*blockDim.x],sc_c[ti + tj*blockDim.x],ti,tj,i,j,k);

}

   // make sure all threads complete shared memory copy
    __syncthreads();

    if((j >= dom->Gfz._js && j < dom->Gfz._je)
      && (i >= dom->Gfz._is && i < dom->Gfz._ie)) {
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
    }
    if((j >= dom->Gfx._js && j < dom->Gfx._je)
      && (i >= dom->Gfx._is && i < dom->Gfx._ie)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._js && j < dom->Gfy._je)
      && (i >= dom->Gfy._is && i < dom->Gfy._ie)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }

  
    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute convective term
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) {
//scalar at face center
real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0f;
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0f;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0f;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0f;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0f;

//convective term of scalar at cell center
real dsudx=(s_u[(ti+1) + tj*blockDim.x]*sc_ue - s_u[ti + tj*blockDim.x]*sc_uw) * ddx;
real dsvdy=(s_v[ti + (tj+1)*blockDim.x]*sc_vn - s_v[ti + tj*blockDim.x]*sc_vs) * ddy;
real dswdz=(s_w1[ti + tj*blockDim.x]   *sc_wt - s_w0[ti + tj*blockDim.x]*sc_wb) * ddz;

//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;

     //move convective term to right hand side 
      s_c[ti + tj*blockDim.x]=-(dsudx+dsvdy+dswdz);
     //Diffusion term
      s_d[ti + tj*blockDim.x]= (ddsdxx+ddsdyy+ddsdzz)*DIFF;

//if(isnan(sc_c[ti + tj*blockDim.x])) printf("\nconv_diff %f %f %f %f %f %d %d %d %d %d\n",dsudx,ddsdxx,ddx,ab,sc_c[ti + tj*blockDim.x],ti,tj,i,j,k);
//Adam-bashforth at t+1/2 step
real sc_conv= ab * s_c[ti + tj*blockDim.x]-ab0*s_c0[ti + tj*blockDim.x];
real sc_diff= ab * s_d[ti + tj*blockDim.x]-ab0*s_d0[ti + tj*blockDim.x];

//if(isnan(s_f[ti + tj*blockDim.x])) printf("\nsource %f %f %f %f %d %d %d %d %d\n",s_d[ti + tj*blockDim.x],s_c[ti + tj*blockDim.x],s_f[ti + tj*blockDim.x],sc_c[ti + tj*blockDim.x],ti,tj,i,j,k);

//TODO should we have s_f at t+1/2 step as well? 
//advance scalar. Take the particle volume fraction into consideration
//real rhs=(sc_diff+s_f[ti + tj*blockDim.x])/(1-p_epsp[ti + tj*blockDim.x]);
real rhs=sc_diff+s_f[ti + tj*blockDim.x];
//real rhs=s_f[ti + tj*blockDim.x];
sc_c[ti + tj*blockDim.x] +=(sc_conv+rhs)*dt;

//if(i==36&&j==38&&k==35) printf("\nsc_diff %f %f %f %f\n",sc_diff,sc_conv,s_f[ti + tj*blockDim.x],sc_c[ti + tj*blockDim.x]);
}

// copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
     && (i >= dom->Gcc._is && i < dom->Gcc._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
//if(isnan(sc_c[ti + tj*blockDim.x])) printf("\ntest_diff %f %f %f %f %d %d %d %d %d\n",s_d[ti + tj*blockDim.x],s_c[ti + tj*blockDim.x],s_f[ti + tj*blockDim.x],sc_c[ti + tj*blockDim.x],ti,tj,i,j,k);
      sc[C] = sc_c[ti + tj*blockDim.x];
      conv[C] = s_c[ti + tj*blockDim.x];
      diff[C] = s_d[ti + tj*blockDim.x];
    }
  }
}







//first order upwind
__global__ void advance_sc_upwind_1st(real DIFF, real *u, real *v, real *w,real *f,real *epsp,
  real *diff0, real *conv0, real *diff, real *conv, real *sc,real *sc0,
  dom_struct *dom, real dt0, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward
  __shared__ real s_wi[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back back
  __shared__ real s_w2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward 2
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
//store scalar related values such as convective and diffusive term.
  __shared__ real s_d0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // diff0
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // conv0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // source term


//store scalar value of difference direction at cell center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top	 center
  __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at  	 center

//store point volume percentage in each grid cell
  __shared__ real p_epsp[MAX_THREADS_DIM * MAX_THREADS_DIM];       
 

  // working constants
  real ab0 = 0.5 * dt / dt0;  // for Adams-Bashforth stepping
  real ab = 1 + ab0;          // for Adams-Bashforth stepping

  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

 
//To calculate the peclet number
  real ui= dom->dx/DIFF;
  real vi= dom->dy/DIFF;
  real wi= dom->dz/DIFF;
 

  int C;

  // loop over z-planes
//  for(int k = dom->Gcc._ksb; k < dom->Gcc._keb; k++) {
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int i = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;

    // shared memory indices
    int ti = threadIdx.x;
    int tj = threadIdx.y;

    // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations

    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];

s_c0[ti + tj*blockDim.x]=conv0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
s_d0[ti + tj*blockDim.x]=diff0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
 s_f[ti + tj*blockDim.x]=    f[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];

p_epsp[ti + tj*blockDim.x]= epsp[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
}

   // make sure all threads complete shared memory copy
    __syncthreads();
 
    if((j >= dom->Gfz._jsb && j < dom->Gfz._jeb)
      && (i >= dom->Gfz._isb && i < dom->Gfz._ieb)) {
      s_wi[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k-1)*dom->Gfz._s2b];
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
      s_w2[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
   }
    if((j >= dom->Gfx._jsb && j < dom->Gfx._jeb)
      && (i >= dom->Gfx._isb && i < dom->Gfx._ieb)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._jsb && j < dom->Gfy._jeb)
      && (i >= dom->Gfy._isb && i < dom->Gfy._ieb)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
 
 
    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute convective term
    // if off the shared memory block boundary
    if((ti > 1 && ti < blockDim.x-2) && (tj > 1 && tj < blockDim.y-2)) 
 
//if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) 
{
//scalar on the west and east u face
real u_c=(s_u[ti + tj*blockDim.x]+s_u[(ti+1) + tj*blockDim.x])/2;// at xm(i)
real v_c=(s_v[ti + tj*blockDim.x]+s_v[ti + (tj+1)*blockDim.x])/2;
real w_c=(s_w1[ti + tj*blockDim.x]+s_w0[ti + tj*blockDim.x])/2;
/*
real u_ce=(s_u[(ti+2) + tj*blockDim.x]+s_u[(ti+1) + tj*blockDim.x])/2;//at xm(i+1)
real u_cw=(s_u[(ti-1) + tj*blockDim.x]+s_u[ti + tj*blockDim.x])/2;

real v_cn=(s_v[ti + (tj+2)*blockDim.x]+s_v[ti + (tj+1)*blockDim.x])/2;
real v_cs=(s_v[ti + (tj-1)*blockDim.x]+s_v[ti + tj*blockDim.x])/2;

real w_ct=(s_w1[ti + tj*blockDim.x]+s_w2[ti + tj*blockDim.x])/2;
real w_cb=(s_wi[ti + tj*blockDim.x]+s_w0[ti + tj*blockDim.x])/2;
*/
real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0;//at x(i)
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;

//1st order upwind scheme
real dsudx,dsvdy,dswdz;
if(u_c>0)	dsudx= (sc_c[ti + tj*blockDim.x]*u_c-sc_uw*s_u[ti + tj*blockDim.x]) *2*ddx;//backward
else		dsudx=-(sc_c[ti + tj*blockDim.x]*u_c-sc_ue*s_u[ti+1 + tj*blockDim.x]) *2* ddx;//forward

if(v_c>0)	dsvdy= (sc_c[ti + tj*blockDim.x]*v_c-sc_vs*s_v[ti + tj*blockDim.x])*2*ddy;//backward
else		dsvdy=-(sc_c[ti + tj*blockDim.x]*v_c-sc_vn*s_v[ti + (tj+1)*blockDim.x])*2*ddy;//forward

if(w_c>0)	dswdz= (sc_c[ti + tj*blockDim.x]*w_c-sc_wb*s_w0[ti + tj*blockDim.x])*2*ddz;//backward
else		dswdz=-(sc_c[ti + tj*blockDim.x]*w_c-sc_wt*s_w1[ti + tj*blockDim.x])*2*ddz;//forward


//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;

     //move convective term to right hand side 
//real sc_conv=-(dsudx+dsvdy+dswdz);
      s_c[ti + tj*blockDim.x]=-(dsudx+dsvdy+dswdz);
     //Diffusion term
//real sc_diff=(ddsdxx+ddsdyy+ddsdzz)*DIFF;
      s_d[ti + tj*blockDim.x]= (ddsdxx+ddsdyy+ddsdzz)*DIFF;

//adam-bashforth at t+1/2 step
real sc_conv= ab * s_c[ti + tj*blockDim.x]-ab0*s_c0[ti + tj*blockDim.x];
real sc_diff= ab * s_d[ti + tj*blockDim.x]-ab0*s_d0[ti + tj*blockDim.x];

//TODO should we have s_f at t+1/2 step as well???
//right hand side of scalar equation
//s_rhs[ti + tj*blockDim.x]=sc_conv+sc_diff+s_f[ti + tj*blockDim.x];
//advance scalar
//Take the particle volume fraction into consideration
//real rhs=(sc_diff+s_f[ti + tj*blockDim.x])/(1-p_epsp[ti + tj*blockDim.x]);
real rhs=sc_diff+s_f[ti + tj*blockDim.x];
sc_c[ti + tj*blockDim.x] +=(sc_conv+rhs)*dt;
}

    // copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
     && (i >= dom->Gcc._is && i < dom->Gcc._ie)
  //  if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)
  //    && (i >= dom->Gcc._isb && i < dom->Gcc._ieb)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      sc[C] = sc_c[ti + tj*blockDim.x];
      conv[C] = s_c[ti + tj*blockDim.x];
      diff[C] = s_d[ti + tj*blockDim.x];
    }
  }
}


//first order upwind
__global__ void advance_sc_upwind_1st_init(real DIFF, real *u, real *v, real *w,real *f,real *epsp,
  real *diff0, real *conv0, real *diff, real *conv, real *sc,real *sc0,
  dom_struct *dom, real dt0, real dt)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward
  __shared__ real s_wi[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back back
  __shared__ real s_w2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward 2
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
//store scalar related values such as convective and diffusive term.
  __shared__ real s_d0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // diff0
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // conv0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv
  __shared__ real s_f[MAX_THREADS_DIM * MAX_THREADS_DIM];       // source term


//store scalar value of difference direction at cell center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top	 center
  __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at  	 center

//store point volume percentage in each grid cell
  __shared__ real p_epsp[MAX_THREADS_DIM * MAX_THREADS_DIM];       
 

  // working constants
  real ab0 = 0.5 * dt / dt0;  // for Adams-Bashforth stepping
  real ab = 1 + ab0;          // for Adams-Bashforth stepping

  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

 
//To calculate the peclet number
  real ui= dom->dx/DIFF;
  real vi= dom->dy/DIFF;
  real wi= dom->dz/DIFF;
 

  int C;

  // loop over z-planes
//  for(int k = dom->Gcc._ksb; k < dom->Gcc._keb; k++) {
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int i = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;

    // shared memory indices
    int ti = threadIdx.x;
    int tj = threadIdx.y;

    // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations

    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];

/*
s_c0[ti + tj*blockDim.x]=conv0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
s_d0[ti + tj*blockDim.x]=diff0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
*/
s_f[ti + tj*blockDim.x]=    f[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];

p_epsp[ti + tj*blockDim.x]= epsp[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
}

   // make sure all threads complete shared memory copy
    __syncthreads();
 
    if((j >= dom->Gfz._jsb && j < dom->Gfz._jeb)
      && (i >= dom->Gfz._isb && i < dom->Gfz._ieb)) {
      s_wi[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k-1)*dom->Gfz._s2b];
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
      s_w2[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
   }
    if((j >= dom->Gfx._jsb && j < dom->Gfx._jeb)
      && (i >= dom->Gfx._isb && i < dom->Gfx._ieb)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._jsb && j < dom->Gfy._jeb)
      && (i >= dom->Gfy._isb && i < dom->Gfy._ieb)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
 
 
    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute convective term
    // if off the shared memory block boundary
    if((ti > 1 && ti < blockDim.x-2) && (tj > 1 && tj < blockDim.y-2)) 
 
//if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) 
{
//scalar on the west and east u face
real u_c=(s_u[ti + tj*blockDim.x]+s_u[(ti+1) + tj*blockDim.x])/2;// at xm(i)
real v_c=(s_v[ti + tj*blockDim.x]+s_v[ti + (tj+1)*blockDim.x])/2;
real w_c=(s_w1[ti + tj*blockDim.x]+s_w0[ti + tj*blockDim.x])/2;
 
real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0;//at x(i)
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;

//1st order upwind scheme
real dsudx,dsvdy,dswdz;
if(u_c>0)	dsudx= (sc_c[ti + tj*blockDim.x]*u_c-sc_uw*s_u[ti + tj*blockDim.x]) *2*ddx;//backward
else		dsudx=-(sc_c[ti + tj*blockDim.x]*u_c-sc_ue*s_u[ti+1 + tj*blockDim.x]) *2* ddx;//forward

if(v_c>0)	dsvdy= (sc_c[ti + tj*blockDim.x]*v_c-sc_vs*s_v[ti + tj*blockDim.x])*2*ddy;//backward
else		dsvdy=-(sc_c[ti + tj*blockDim.x]*v_c-sc_vn*s_v[ti + (tj+1)*blockDim.x])*2*ddy;//forward

if(w_c>0)	dswdz= (sc_c[ti + tj*blockDim.x]*w_c-sc_wb*s_w0[ti + tj*blockDim.x])*2*ddz;//backward
else		dswdz=-(sc_c[ti + tj*blockDim.x]*w_c-sc_wt*s_w1[ti + tj*blockDim.x])*2*ddz;//forward


//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;

     //move convective term to right hand side 
//real sc_conv=-(dsudx+dsvdy+dswdz);
      s_c[ti + tj*blockDim.x]=-(dsudx+dsvdy+dswdz);
     //Diffusion term
//real sc_diff=(ddsdxx+ddsdyy+ddsdzz)*DIFF;
      s_d[ti + tj*blockDim.x]= (ddsdxx+ddsdyy+ddsdzz)*DIFF;

//adam-bashforth at t+1/2 step
real sc_conv=  s_c[ti + tj*blockDim.x] ;
real sc_diff=  s_d[ti + tj*blockDim.x] ;

//TODO should we have s_f at t+1/2 step as well???
//right hand side of scalar equation
//s_rhs[ti + tj*blockDim.x]=sc_conv+sc_diff+s_f[ti + tj*blockDim.x];
//advance scalar
//Take the particle volume fraction into consideration
//real rhs=(sc_diff+s_f[ti + tj*blockDim.x])/(1-p_epsp[ti + tj*blockDim.x]);
real rhs=sc_diff+s_f[ti + tj*blockDim.x];
sc_c[ti + tj*blockDim.x] +=(sc_conv+rhs)*dt;
}

    // copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
     && (i >= dom->Gcc._is && i < dom->Gcc._ie)
  //  if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)
  //    && (i >= dom->Gcc._isb && i < dom->Gcc._ieb)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      sc[C] = sc_c[ti + tj*blockDim.x];
      conv[C] = s_c[ti + tj*blockDim.x];
      diff[C] = s_d[ti + tj*blockDim.x];
    }
  }
}






__global__ void sc_conv_diff(real DIFF, real *u, real *v, real *w,
  real *diff0, real *conv0, real *diff, real *conv, real *sc0,
  dom_struct *dom)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward
  __shared__ real s_wi[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back back
  __shared__ real s_w2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward 2
  __shared__ real s_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       // u
  __shared__ real s_v[MAX_THREADS_DIM * MAX_THREADS_DIM];       // v
//store scalar related values such as convective and diffusive term.
  __shared__ real s_d0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // diff0
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  __shared__ real s_c0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // conv0
  __shared__ real s_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // conv



//store scalar value of difference direction at cell center
  __shared__ real sc_b[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at bottom center
  __shared__ real sc_t[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at top	 center
  __shared__ real sc_c[MAX_THREADS_DIM * MAX_THREADS_DIM];       // sc at  	 center
 
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

 
//To calculate the peclet number
  real ui= dom->dx/DIFF;
  real vi= dom->dy/DIFF;
  real wi= dom->dz/DIFF;
 

  int C;

  // loop over z-planes
//  for(int k = dom->Gcc._ksb; k < dom->Gcc._keb; k++) {
  for(int k = dom->Gcc._ks; k < dom->Gcc._ke; k++) {
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int i = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;

    // shared memory indices
    int ti = threadIdx.x;
    int tj = threadIdx.y;

    // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations

    if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)&& (i >= dom->Gcc._isb && i < dom->Gcc._ieb)) {
sc_c[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + k*dom->Gcc._s2b];
sc_b[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b];
sc_t[ti + tj*blockDim.x]=sc0[i+j*dom->Gcc._s1b + (k+1)*dom->Gcc._s2b];
}

   // make sure all threads complete shared memory copy
    __syncthreads();

 
    if((j >= dom->Gfz._jsb && j < dom->Gfz._jeb)
      && (i >= dom->Gfz._isb && i < dom->Gfz._ieb)) {
      s_wi[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k-1)*dom->Gfz._s2b];
      s_w0[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
      s_w2[ti + tj*blockDim.x] = w[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
   }
    if((j >= dom->Gfx._jsb && j < dom->Gfx._jeb)
      && (i >= dom->Gfx._isb && i < dom->Gfx._ieb)) {
      s_u[ti + tj*blockDim.x] = u[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
    if((j >= dom->Gfy._jsb && j < dom->Gfy._jeb)
      && (i >= dom->Gfy._isb && i < dom->Gfy._ieb)) {
      s_v[ti + tj*blockDim.x] = v[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
    }
     // make sure all threads complete shared memory copy
    __syncthreads();

    // compute convective term
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) {
 
//scalar on the west and east u face
real u_c=(s_u[ti + tj*blockDim.x]+s_u[(ti+1) + tj*blockDim.x])/2;// at xm(i)
real v_c=(s_v[ti + tj*blockDim.x]+s_v[ti + (tj+1)*blockDim.x])/2;
real w_c=(s_w1[ti + tj*blockDim.x]+s_w0[ti + tj*blockDim.x])/2;

real u_ce=(s_u[(ti+2) + tj*blockDim.x]+s_u[(ti+1) + tj*blockDim.x])/2;//at xm(i+1)
real u_cw=(s_u[(ti-1) + tj*blockDim.x]+s_u[ti + tj*blockDim.x])/2;


real v_cn=(s_v[ti + (tj+2)*blockDim.x]+s_v[ti + (tj+1)*blockDim.x])/2;
real v_cs=(s_v[ti + (tj-1)*blockDim.x]+s_v[ti + tj*blockDim.x])/2;

real w_ct=(s_w1[ti + tj*blockDim.x]+s_w2[ti + tj*blockDim.x])/2;
real w_cb=(s_wi[ti + tj*blockDim.x]+s_w0[ti + tj*blockDim.x])/2;

real sc_uw=(sc_c[ti + tj*blockDim.x]+sc_c[ti-1 + tj*blockDim.x])/2.0;//at x(i)
real sc_ue=(sc_c[ti + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x])/2.0;

real sc_vs=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj-1)*blockDim.x])/2.0;
real sc_vn=(sc_c[ti + tj*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x])/2.0;

real sc_wb=(sc_b[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;
real sc_wt=(sc_t[ti + tj*blockDim.x]+sc_c[ti + tj*blockDim.x])/2.0;

real Pe_u= u_c*ui;
real Pe_v= v_c*vi;
real Pe_w= w_c*wi;

real dsudx_up1=( 3*sc_c[ti + tj*blockDim.x]*u_c-4*sc_uw*s_u[ti + tj*blockDim.x]+ sc_c[ti-1 + tj*blockDim.x]*u_cw)*ddx;//backward
real dsudx_up2=(-3*sc_c[ti + tj*blockDim.x]*u_c+4*sc_ue*s_u[ti +1 + tj*blockDim.x] - sc_c[ti+1 + tj*blockDim.x]*u_ce)*ddx;//forward

real dsvdy_up1=( 3*sc_c[ti + tj*blockDim.x]*v_c-4*sc_vs*s_v[ti + tj*blockDim.x]+ sc_c[ti + (tj-1)*blockDim.x]*v_cs)*ddy;//backward
real dsvdy_up2=(-3*sc_c[ti + tj*blockDim.x]*v_c+4*sc_vn*s_v[ti + (tj+1)*blockDim.x]-sc_c[ti+(tj+1)*blockDim.x]*v_cn)*ddy;//forward

real dswdz_up1=( 3*sc_c[ti + tj*blockDim.x]*w_c-4*sc_wb*s_w0[ti + tj*blockDim.x]+sc_b[ti +tj*blockDim.x]*w_cb)*ddz;//backward
real dswdz_up2=(-3*sc_c[ti + tj*blockDim.x]*w_c+4*sc_wt*s_w1[ti + tj*blockDim.x]-sc_t[ti+tj*blockDim.x]*w_ct)*ddz;//forward


//convective term of scalar at cell center -- central difference
real dsudx_ce=(s_u[(ti+1) + tj*blockDim.x]*sc_ue - s_u[ti + tj*blockDim.x]*sc_uw) * ddx;//at xm(i)
real dsvdy_ce=(s_v[ti + (tj+1)*blockDim.x]*sc_vn - s_v[ti + tj*blockDim.x]*sc_vs) * ddy;
real dswdz_ce=(s_w1[ti + tj*blockDim.x]   *sc_wt - s_w0[ti + tj*blockDim.x]*sc_wb) * ddz;

if(Pe_u>2.0) dsudx=dsudx_up1;
else if(Pe_u<-2.0) dsudx=dsudx_up2;
else dsudx=dsudx_ce;

if(Pe_v>2.0) dsvdy=dsvdy_up1;
else if(Pe_v<-2.0) dsvdy=dsvdy_up2;
else dsvdy=dsvdy_ce;

if(Pe_w>2.0) dswdz=dswdz_up1;
else if(Pe_w<-2.0) dswdz=dswdz_up2;
else dswdz=dswdz_ce;


//diffusive term of scalar at cell center
real ddsdxx=(sc_c[ti-1 + tj*blockDim.x]+sc_c[ti+1 + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddx* ddx;
real ddsdyy=(sc_c[ti + (tj-1)*blockDim.x]+sc_c[ti + (tj+1)*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddy* ddy;
real ddsdzz=(sc_b[ti + tj*blockDim.x]+sc_t[ti + tj*blockDim.x]-2*sc_c[ti + tj*blockDim.x]) *ddz* ddz;

     //move convective term to right hand side 
//real sc_conv=-(dsudx+dsvdy+dswdz);
      s_c[ti + tj*blockDim.x]=-(dsudx+dsvdy+dswdz);
     //Diffusion term
//real sc_diff=(ddsdxx+ddsdyy+ddsdzz)*DIFF;
      s_d[ti + tj*blockDim.x]= (ddsdxx+ddsdyy+ddsdzz)*DIFF;
}

    // copy shared memory back to global, without copying boundary ghost values
//sc stores n+1 timestep, conv and diff are n timestep!!
    if((j >= dom->Gcc._js && j < dom->Gcc._je)
     && (i >= dom->Gcc._is && i < dom->Gcc._ie)
  //  if((j >= dom->Gcc._jsb && j < dom->Gcc._jeb)
  //    && (i >= dom->Gcc._isb && i < dom->Gcc._ieb)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
      C = i + j*dom->Gcc._s1b + k*dom->Gcc._s2b;
      conv[C] = s_c[ti + tj*blockDim.x];
      diff[C] = s_d[ti + tj*blockDim.x];

      conv0[C] = conv[C];
      diff0[C] = diff[C];
    }
  }
}


