extern "C"
{
#include "bluebottle.h"
}

#include <cuda.h>
#include <helper_cuda.h>

//#include "cuda_quadrature.h"

extern "C"
void cuda_U_star_test_exp(void)
{
  int i, j, k;  // iterators
  int C;    // cell locations

  real *u_a = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_c = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_e = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *v_a = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_c = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_e = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *w_a = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_c = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_e = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfz.s3b * sizeof(real);

  // min and max error search
  real u_err_min = FLT_MAX;
  real u_err_max = FLT_MIN;
  real v_err_min = FLT_MAX;
  real v_err_max = FLT_MIN;
  real w_err_min = FLT_MAX;
  real w_err_max = FLT_MIN;

  printf("\nIntermediate velocity calculation validation:\n\n");
  printf("  u = exp(x), v = exp(y), w = exp(z)\n\n");

  // set up expected solution
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_a[C] = nu;
        u_a[C] -= 2 * exp((i-1.5)*Dom.dx);
        u_a[C] -= exp((j-1.0)*Dom.dy);
        u_a[C] -= exp((k-1.0)*Dom.dz);
        u_a[C] *= dt * exp((i-1.5)*Dom.dx);
        u_a[C] += exp((i-1.5)*Dom.dx);
        u[C] = u_a[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_a[C] = nu;
        v_a[C] -= 2 * exp((j-1.5)*Dom.dy);
        v_a[C] -= exp((k-1.0)*Dom.dz);
        v_a[C] -= exp((i-1.0)*Dom.dx);
        v_a[C] *= dt * exp((j-1.5)*Dom.dy);
        v_a[C] += exp((j-1.5)*Dom.dy);
        v[C] = v_a[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_a[C] = nu;
        w_a[C] -= 2 * exp((k-1.5)*Dom.dz);
        w_a[C] -= exp((i-1.0)*Dom.dx);
        w_a[C] -= exp((j-1.0)*Dom.dy);
        w_a[C] *= dt * exp((k-1.5)*Dom.dz);
        w_a[C] += exp((k-1.5)*Dom.dz);
        w[C] = w_a[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("  Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // initialize input velocity fields
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u[C] = exp((i-1.5)*Dom.dx);
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v[C] = exp((j-1.5)*Dom.dy);
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w[C] = exp((k-1.5)*Dom.dz);
      }
    }
  }

  // write initial fields
  rec_paraview_stepnum_out++;
  printf("  Writing initial fields to:    out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // push fields to device
  printf("\n  Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // call code to test
  printf("  Running cuda_U_star_2()...");
  cuda_U_star_2();
  printf("done.\n");

  // pull fields back to host
  printf("  Pulling fields back to host...");
  cuda_dom_pull();
  printf("done.\n");

  // write computed solution
  rec_paraview_stepnum_out++;
  printf("\n  Writing computed solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // copy results and compute error
  for(k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(i = Dom.Gfx.is; i < Dom.Gfx.ie; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_c[C] = u[C];
        u_e[C] = (u_c[C] - u_a[C]) / u_a[C];
        if(fabs(u_e[C]) > u_err_max) u_err_max = fabs(u_e[C]);
        if(fabs(u_e[C]) < u_err_min) u_err_min = fabs(u_e[C]);
        u[C] = u_e[C];
      }
    }
  }
  for(k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(j = Dom.Gfy.js; j < Dom.Gfy.je; j++) {
      for(i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_c[C] = v[C];
        v_e[C] = (v_c[C] - v_a[C]) / v_a[C];
        if(fabs(v_e[C]) > v_err_max) v_err_max = fabs(v_e[C]);
        if(fabs(v_e[C]) < v_err_min) v_err_min = fabs(v_e[C]);
        v[C] = v_e[C];
      }
    }
  }
  for(k = Dom.Gfz.ks; k < Dom.Gfz.ke; k++) {
    for(j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_c[C] = w[C];
        w_e[C] = (w_c[C] - w_a[C]) / w_a[C];
        if(fabs(w_e[C]) > w_err_max) w_err_max = fabs(w_e[C]);
        if(fabs(w_e[C]) < w_err_min) w_err_min = fabs(w_e[C]);
        w[C] = w_e[C];
      }
    }
  }

  // write error difference
  rec_paraview_stepnum_out++;
  printf("  Writing error difference to:  out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  printf("\n  Error summary:\n");
  printf("  Velocity component:     minimum error:     maximum error:\n");
  printf("          u              %12.3e       %12.3e\n",
    u_err_min, u_err_max);
  printf("          v              %12.3e       %12.3e\n",
    v_err_min, v_err_max);
  printf("          w              %12.3e       %12.3e\n\n",
    w_err_min, w_err_max);

  // clean up
  free(u_a);
  free(u_c);
  free(u_e);
  free(v_a);
  free(v_c);
  free(v_e);
  free(w_a);
  free(w_c);
  free(w_e);
}

extern "C"
void cuda_U_star_test_cos(void)
{
  int i, j, k;  // iterators
  int C;    // cell locations

  real *u_a = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_c = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_e = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *v_a = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_c = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_e = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *w_a = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_c = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_e = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfz.s3b * sizeof(real);

  // min and max error search
  real u_err_min = FLT_MAX;
  real u_err_max = FLT_MIN;
  real v_err_min = FLT_MAX;
  real v_err_max = FLT_MIN;
  real w_err_min = FLT_MAX;
  real w_err_max = FLT_MIN;

  printf("Intermediate velocity calculation validation:\n\n");
  printf("  u = cos(x), v = cos(y), w = cos(z)\n\n");

  // set up expected solution
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_a[C] = - nu * cos((i-1.5)*Dom.dx);
        u_a[C] += 2 * sin((i-1.5)*Dom.dx) * cos((i-1.5)*Dom.dx);
        u_a[C] += cos((i-1.5)*Dom.dx) * sin((j-1.0)*Dom.dy);
        u_a[C] += cos((i-1.5)*Dom.dx) * sin((k-1.0)*Dom.dz);
        u_a[C] *= dt;
        u_a[C] += cos((i-1.5)*Dom.dx);
        u[C] = u_a[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_a[C] = - nu * cos((j-1.5)*Dom.dy);
        v_a[C] += 2 * sin((j-1.5)*Dom.dy) * cos((j-1.5)*Dom.dy);
        v_a[C] += cos((j-1.5)*Dom.dy) * sin((k-1.0)*Dom.dz);
        v_a[C] += cos((j-1.5)*Dom.dy) * sin((i-1.0)*Dom.dx);
        v_a[C] *= dt;
        v_a[C] += cos((j-1.5)*Dom.dy);
        v[C] = v_a[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_a[C] = - nu * cos((k-1.5)*Dom.dz);
        w_a[C] += 2 * sin((k-1.5)*Dom.dz) * cos((k-1.5)*Dom.dz);
        w_a[C] += cos((k-1.5)*Dom.dz) * sin((i-1.0)*Dom.dx);
        w_a[C] += cos((k-1.5)*Dom.dz) * sin((j-1.0)*Dom.dy);
        w_a[C] *= dt;
        w_a[C] += cos((k-1.5)*Dom.dz);
        w[C] = w_a[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("  Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // initialize input pressure and velocity fields
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u[C] = cos((i-1.5)*Dom.dx);
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v[C] = cos((j-1.5)*Dom.dy);
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w[C] = cos((k-1.5)*Dom.dz);
      }
    }
  }

  // write initial fields
  rec_paraview_stepnum_out++;
  printf("  Writing initial fields to:    out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // push fields to device
  printf("\n  Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // call code to test
  printf("  Running cuda_CN_rhs()...");
  cuda_U_star_2();
  printf("done.\n");

  // pull fields back to host
  printf("  Pulling fields back to host...");
  cuda_dom_pull();
  printf("done.\n");

  // write computed solution
  rec_paraview_stepnum_out++;
  printf("\n  Writing computed solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // copy results and compute error
  for(k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(i = Dom.Gfx.is; i < Dom.Gfx.ie; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_c[C] = u[C];
        u_e[C] = (u_c[C] - u_a[C]);// / u_a[C];
        if(fabs(u_e[C]) > u_err_max) u_err_max = fabs(u_e[C]);
        if(fabs(u_e[C]) < u_err_min) u_err_min = fabs(u_e[C]);
        u[C] = u_e[C];
      }
    }
  }
  for(k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(j = Dom.Gfy.js; j < Dom.Gfy.je; j++) {
      for(i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_c[C] = v[C];
        v_e[C] = (v_c[C] - v_a[C]);// / v_a[C];
        if(fabs(v_e[C]) > v_err_max) v_err_max = fabs(v_e[C]);
        if(fabs(v_e[C]) < v_err_min) v_err_min = fabs(v_e[C]);
        v[C] = v_e[C];
      }
    }
  }
  for(k = Dom.Gfz.ks; k < Dom.Gfz.ke; k++) {
    for(j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_c[C] = w[C];
        w_e[C] = (w_c[C] - w_a[C]);// / w_a[C];
        if(fabs(w_e[C]) > w_err_max) w_err_max = fabs(w_e[C]);
        if(fabs(w_e[C]) < w_err_min) w_err_min = fabs(w_e[C]);
        w[C] = w_e[C];
      }
    }
  }

  // write error difference
  rec_paraview_stepnum_out++;
  printf("  Writing error difference to:  out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  printf("\n  Error summary:\n");
  printf("  Velocity component:     minimum error:     maximum error:\n");
  printf("          u              %12.3e       %12.3e\n",
    u_err_min, u_err_max);
  printf("          v              %12.3e       %12.3e\n",
    v_err_min, v_err_max);
  printf("          w              %12.3e       %12.3e\n\n",
    w_err_min, w_err_max);

  // clean up
  free(u_a);
  free(u_c);
  free(u_e);
  free(v_a);
  free(v_c);
  free(v_e);
  free(w_a);
  free(w_c);
  free(w_e);
}

__global__ void memcpy_u_star_test(real *dst, real *src, dom_struct *dom)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  for(int i = dom->Gfx._is; i < dom->Gfx._ie; i++) {
    dst[i + j * dom->Gfx._s1b + k * dom->Gfx._s2b] = src[i + j * dom->Gfx._s1b
      + k * dom->Gfx._s2b];
  }
}

__global__ void memcpy_v_star_test(real *dst, real *src, dom_struct *dom)
{
  int k = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  for(int j = dom->Gfy._js; j < dom->Gfy._je; j++) {
    dst[i + j * dom->Gfy._s1b + k * dom->Gfy._s2b] = src[i + j * dom->Gfy._s1b
      + k * dom->Gfy._s2b];
  }
}

__global__ void memcpy_w_star_test(real *dst, real *src, dom_struct *dom)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  for(int k = dom->Gfz._ks; k < dom->Gfz._ke; k++) {
    dst[i + j * dom->Gfz._s1b + k * dom->Gfz._s2b] = src[i + j * dom->Gfz._s1b
      + k * dom->Gfz._s2b];
  }
}

__global__ void PP_memcpy_p_test(real *dst, real *src, dom_struct *dom)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int k = blockIdx.y * blockDim.y + threadIdx.y;
  for(int i = dom->Gcc._is-DOM_BUF; i < dom->Gcc._ie-DOM_BUF; i++) {
    dst[(i+DOM_BUF) + (j+DOM_BUF) * dom->Gfx._s1b
      + (k+DOM_BUF) * dom->Gcc._s2b] = src[i + j * dom->Gcc._s1b
      + k * dom->Gcc._s2b];
  }
}

extern "C"
void cuda_BC_test(void)
{
  int i, j, k;  // iterators
  int C;    // cell locations

  real *p_p_i = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // periodic input
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_p_o = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // periodic output
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_p_e = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // periodic error
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_d_i = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // Dirichlet input
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_d_o = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // Dirichlet output
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_d_e = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // Dirichlet error
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_n_i = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // Neumann input
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_n_o = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // Neumann output
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *p_n_e = (real*) malloc(Dom.Gcc.s3b * sizeof(real)); // Neumann error
  // cpumem += Dom.Gcc.s3b * sizeof(real);

  real *u_p_i = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // periodic input
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_p_o = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // periodic output
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_p_e = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // periodic error
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_d_i = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // Dirichlet input
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_d_o = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // Dirichlet output
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_d_e = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // Dirichlet error
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_n_i = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // Neumann input
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_n_o = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // Neumann output
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_n_e = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // Neumann error
  // cpumem += Dom.Gfx.s3b * sizeof(real);

  real *v_p_i = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // periodic input
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_p_o = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // periodic output
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_p_e = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // periodic error
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_d_i = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // Dirichlet input
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_d_o = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // Dirichlet output
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_d_e = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // Dirichlet error
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_n_i = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // Neumann input
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_n_o = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // Neumann output
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_n_e = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // Neumann error
  // cpumem += Dom.Gfy.s3b * sizeof(real);

  real *w_p_i = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // periodic input
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_p_o = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // periodic output
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_p_e = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // periodic error
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_d_i = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // Dirichlet input
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_d_o = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // Dirichlet output
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_d_e = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // Dirichlet error
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_n_i = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // Neumann input
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_n_o = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // Neumann output
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_n_e = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // Neumann error
  // cpumem += Dom.Gfz.s3b * sizeof(real);

  // min and max error search
  real p_p_err_min = FLT_MAX;
  real p_p_err_max = FLT_MIN;
  real p_d_err_min = FLT_MAX;
  real p_d_err_max = FLT_MIN;
  real p_n_err_min = FLT_MAX;
  real p_n_err_max = FLT_MIN;
  real u_p_err_min = FLT_MAX;
  real u_p_err_max = FLT_MIN;
  real u_d_err_min = FLT_MAX;
  real u_d_err_max = FLT_MIN;
  real u_n_err_min = FLT_MAX;
  real u_n_err_max = FLT_MIN;
  real v_p_err_min = FLT_MAX;
  real v_p_err_max = FLT_MIN;
  real v_d_err_min = FLT_MAX;
  real v_d_err_max = FLT_MIN;
  real v_n_err_min = FLT_MAX;
  real v_n_err_max = FLT_MIN;
  real w_p_err_min = FLT_MAX;
  real w_p_err_max = FLT_MIN;
  real w_d_err_min = FLT_MAX;
  real w_d_err_max = FLT_MIN;
  real w_n_err_min = FLT_MAX;
  real w_n_err_max = FLT_MIN;

  printf("\nBoundary condition application validation:\n");

  // periodic field (on -1 <= x <= 1, -1 <= y <= 1, -1 <= z <= 1)
  printf("\n  Periodic boundary conditions:\n");
  printf("    p = cos(pi*x) + cos(pi*y) + cos(pi*z)\n");
  printf("    u = cos(pi*x) + cos(pi*y) + cos(pi*z)\n");
  printf("    v = cos(pi*x) + cos(pi*y) + cos(pi*z)\n");
  printf("    w = cos(pi*x) + cos(pi*y) + cos(pi*z)\n\n");

  // write input fields
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p_p_i[C] = cos(PI * (i-0.5)*Dom.dx);
        p_p_i[C] += cos(PI * (j-0.5)*Dom.dy);
        p_p_i[C] += cos(PI * (k-0.5)*Dom.dz);
        p[C] = p_p_i[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
      C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
      u_p_i[C] = cos(PI * (i-1.0)*Dom.dx);
      u_p_i[C] += cos(PI * (j-0.5)*Dom.dy);
      u_p_i[C] += cos(PI * (k-0.5)*Dom.dz);
      u[C] = u_p_i[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
      C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
      v_p_i[C] = cos(PI * (i-0.5)*Dom.dx);
      v_p_i[C] += cos(PI * (j-1.0)*Dom.dy);
      v_p_i[C] += cos(PI * (k-0.5)*Dom.dz);
      v[C] = v_p_i[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
      C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
      w_p_i[C] = cos(PI * (i-0.5)*Dom.dx);
      w_p_i[C] += cos(PI * (j-0.5)*Dom.dy);
      w_p_i[C] += cos(PI * (k-1.0)*Dom.dz);
      w[C] = w_p_i[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("    Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // push to device
  printf("\n    Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // set BC to PERIODIC (overwrite input file for testing)
  printf("    Overwriting input file boundary conditions to ensure ");
  printf("PERIODIC...");
  bc.pW = PERIODIC;
  bc.pE = PERIODIC;
  bc.pS = PERIODIC;
  bc.pN = PERIODIC;
  bc.pB = PERIODIC;
  bc.pT = PERIODIC;
  bc.uW = PERIODIC;
  bc.uE = PERIODIC;
  bc.uS = PERIODIC;
  bc.uN = PERIODIC;
  bc.uB = PERIODIC;
  bc.uT = PERIODIC;
  bc.vW = PERIODIC;
  bc.vE = PERIODIC;
  bc.vS = PERIODIC;
  bc.vN = PERIODIC;
  bc.vB = PERIODIC;
  bc.vT = PERIODIC;
  bc.wW = PERIODIC;
  bc.wE = PERIODIC;
  bc.wS = PERIODIC;
  bc.wN = PERIODIC;
  bc.wB = PERIODIC;
  bc.wT = PERIODIC;
  printf("done.\n");

  // apply BC
  printf("    Running cuda_BC()...");
  cuda_dom_BC();
  printf("done.\n");

  // pull to host
  printf("    Pulling fields back to host...");
  cuda_dom_pull();
  printf("done.\n");

  // write computed solution
  rec_paraview_stepnum_out++;
  printf("\n    Writing computed solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // copy results and compute error
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p_p_o[C] = p[C];
        p_p_e[C] = p_p_o[C] - p_p_i[C];
        if(fabs(p_p_e[C]) > p_p_err_max) p_p_err_max = fabs(p_p_e[C]);
        if(fabs(p_p_e[C]) < p_p_err_min) p_p_err_min = fabs(p_p_e[C]);
        p[C] = p_p_e[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_p_o[C] = u[C];
        u_p_e[C] = u_p_o[C] - u_p_i[C];
        if(fabs(u_p_e[C]) > u_p_err_max) u_p_err_max = fabs(u_p_e[C]);
        if(fabs(u_p_e[C]) < u_p_err_min) u_p_err_min = fabs(u_p_e[C]);
        u[C] = u_p_e[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_p_o[C] = v[C];
        v_p_e[C] = v_p_o[C] - v_p_i[C];
        if(fabs(v_p_e[C]) > v_p_err_max) v_p_err_max = fabs(v_p_e[C]);
        if(fabs(v_p_e[C]) < v_p_err_min) v_p_err_min = fabs(v_p_e[C]);
        v[C] = v_p_e[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_p_o[C] = w[C];
        w_p_e[C] = w_p_o[C] - w_p_i[C];
        if(fabs(w_p_e[C]) > w_p_err_max) w_p_err_max = fabs(w_p_e[C]);
        if(fabs(w_p_e[C]) < w_p_err_min) w_p_err_min = fabs(w_p_e[C]);
        w[C] = w_p_e[C];
      }
    }
  }

  // write error difference
  rec_paraview_stepnum_out++;
  printf("    Writing error difference to:  out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  printf("\n    Error summary:\n");
  printf("    Field component:      minimum error:     maximum error:\n");
  printf("          p              %12.3e       %12.3e\n",
    p_p_err_min, p_p_err_max);
  printf("          u              %12.3e       %12.3e\n",
    u_p_err_min, u_p_err_max);
  printf("          v              %12.3e       %12.3e\n",
    v_p_err_min, v_p_err_max);
  printf("          w              %12.3e       %12.3e\n",
    w_p_err_min, w_p_err_max);

  // Dirichlet field (on -1 <= x <= 1, -1 <= y <= 1, -1 <= z <= 1)
  printf("\n  Dirichlet boundary conditions:\n");
  printf("    p = sin(pi*x) * sin(pi*y) * sin(pi*z)\n");
  printf("    u = sin(pi*x) * sin(pi*y) * sin(pi*z)\n");
  printf("    v = sin(pi*x) * sin(pi*y) * sin(pi*z)\n");
  printf("    w = sin(pi*x) * sin(pi*y) * sin(pi*z)\n\n");

  // write input field
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p_d_i[C] = sin(PI * (i-0.5)*Dom.dx);
        p_d_i[C] *= sin(PI * (j-0.5)*Dom.dy);
        p_d_i[C] *= sin(PI * (k-0.5)*Dom.dz);
        p[C] = p_d_i[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_d_i[C] = sin(PI * (i-1.0)*Dom.dx);
        u_d_i[C] *= sin(PI * (j-0.5)*Dom.dy);
        u_d_i[C] *= sin(PI * (k-0.5)*Dom.dz);
        u[C] = u_d_i[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_d_i[C] = sin(PI * (i-0.5)*Dom.dx);
        v_d_i[C] *= sin(PI * (j-1.0)*Dom.dy);
        v_d_i[C] *= sin(PI * (k-0.5)*Dom.dz);
        v[C] = v_d_i[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_d_i[C] = sin(PI * (i-0.5)*Dom.dx);
        w_d_i[C] *= sin(PI * (j-0.5)*Dom.dy);
        w_d_i[C] *= sin(PI * (k-1.0)*Dom.dz);
        w[C] = w_d_i[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("    Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // push to device
  printf("\n    Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // set BC to DIRICHLET (overwrite input file for testing)
  printf("    Overwriting input file boundary conditions to ensure ");
  printf("DIRICHLET...");
  bc.pW = DIRICHLET;
  bc.pWD = 0.;
  bc.pE = DIRICHLET;
  bc.pED = 0.;
  bc.pS = DIRICHLET;
  bc.pSD = 0.;
  bc.pN = DIRICHLET;
  bc.pND = 0.;
  bc.pB = DIRICHLET;
  bc.pBD = 0.;
  bc.pT = DIRICHLET;
  bc.pTD = 0.;

  bc.uW = DIRICHLET;
  bc.uWD = 0.;
  bc.uE = DIRICHLET;
  bc.uED = 0.;
  bc.uS = DIRICHLET;
  bc.uSD = 0.;
  bc.uN = DIRICHLET;
  bc.uND = 0.;
  bc.uB = DIRICHLET;
  bc.uBD = 0.;
  bc.uT = DIRICHLET;
  bc.uTD = 0.;

  bc.vW = DIRICHLET;
  bc.vWD = 0.;
  bc.vE = DIRICHLET;
  bc.vED = 0.;
  bc.vS = DIRICHLET;
  bc.vSD = 0.;
  bc.vN = DIRICHLET;
  bc.vND = 0.;
  bc.vB = DIRICHLET;
  bc.vBD = 0.;
  bc.vT = DIRICHLET;
  bc.vTD = 0.;

  bc.wW = DIRICHLET;
  bc.wWD = 0.;
  bc.wE = DIRICHLET;
  bc.wED = 0.;
  bc.wS = DIRICHLET;
  bc.wSD = 0.;
  bc.wN = DIRICHLET;
  bc.wND = 0.;
  bc.wB = DIRICHLET;
  bc.wBD = 0.;
  bc.wT = DIRICHLET;
  bc.wTD = 0.;
  printf("done.\n");

  // apply BC
  printf("    Running cuda_BC()...");
  cuda_dom_BC();
  printf("done.\n");

  // pull to host
  printf("    Pulling fields back to host...");
  cuda_dom_pull();
  printf("done.\n");

  // write computed solution
  rec_paraview_stepnum_out++;
  printf("\n    Writing computed solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // copy results and compute error
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p_d_o[C] = p[C];
        p_d_e[C] = p_d_o[C] - p_d_i[C];
        if(fabs(p_d_e[C]) > p_d_err_max) p_d_err_max = fabs(p_d_e[C]);
        if(fabs(p_d_e[C]) < p_d_err_min) p_d_err_min = fabs(p_d_e[C]);
        p[C] = p_d_e[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_d_o[C] = u[C];
        u_d_e[C] = u_d_o[C] - u_d_i[C];
        if(fabs(u_d_e[C]) > u_d_err_max) u_d_err_max = fabs(u_d_e[C]);
        if(fabs(u_d_e[C]) < u_d_err_min) u_d_err_min = fabs(u_d_e[C]);
        u[C] = u_d_e[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_d_o[C] = v[C];
        v_d_e[C] = v_d_o[C] - v_d_i[C];
        if(fabs(v_d_e[C]) > v_d_err_max) v_d_err_max = fabs(v_d_e[C]);
        if(fabs(v_d_e[C]) < v_d_err_min) v_d_err_min = fabs(v_d_e[C]);
        v[C] = v_d_e[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_d_o[C] = w[C];
        w_d_e[C] = w_d_o[C] - w_d_i[C];
        if(fabs(w_d_e[C]) > w_d_err_max) w_d_err_max = fabs(w_d_e[C]);
        if(fabs(w_d_e[C]) < w_d_err_min) w_d_err_min = fabs(w_d_e[C]);
        w[C] = w_d_e[C];
      }
    }
  }

  // write error difference
  rec_paraview_stepnum_out++;
  printf("    Writing error difference to:  out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  printf("\n    Error summary:\n");
  printf("    Field component:      minimum error:     maximum error:\n");
  printf("          p              %12.3e       %12.3e\n",
    p_d_err_min, p_d_err_max);
  printf("          u              %12.3e       %12.3e\n",
    u_d_err_min, u_d_err_max);
  printf("          v              %12.3e       %12.3e\n",
    v_d_err_min, v_d_err_max);
  printf("          w              %12.3e       %12.3e\n",
    w_d_err_min, w_d_err_max);

  // Neumann field (on -1 <= x <= 1, -1 <= y <= 1, -1 <= z <= 1)
  printf("\n  Neumann boundary conditions:\n");
  printf("    p = cos(pi*x) + cos(pi*y) + cos(pi*z)\n");
  printf("    u = cos(pi*x) + cos(pi*y) + cos(pi*z)\n");
  printf("    v = cos(pi*x) + cos(pi*y) + cos(pi*z)\n");
  printf("    w = cos(pi*x) + cos(pi*y) + cos(pi*z)\n\n");

  // write input field
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p_n_i[C] = cos(PI * (i-0.5)*Dom.dx);
        p_n_i[C] += cos(PI * (j-0.5)*Dom.dy);
        p_n_i[C] += cos(PI * (k-0.5)*Dom.dz);
        p[C] = p_n_i[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_n_i[C] = cos(PI * (i-1.0)*Dom.dx);
        u_n_i[C] += cos(PI * (j-0.5)*Dom.dy);
        u_n_i[C] += cos(PI * (k-0.5)*Dom.dz);
        u[C] = u_n_i[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_n_i[C] = cos(PI * (i-0.5)*Dom.dx);
        v_n_i[C] += cos(PI * (j-1.0)*Dom.dy);
        v_n_i[C] += cos(PI * (k-0.5)*Dom.dz);
        v[C] = v_n_i[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_n_i[C] = cos(PI * (i-0.5)*Dom.dx);
        w_n_i[C] += cos(PI * (j-0.5)*Dom.dy);
        w_n_i[C] += cos(PI * (k-1.0)*Dom.dz);
        w[C] = w_n_i[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("    Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // push to device
  printf("\n    Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // set BC to NEUMANN (overwrite input file for testing)
  printf("    Overwriting input file boundary conditions to ensure ");
  printf("NEUMANN...");
  bc.pW = NEUMANN;
  bc.pE = NEUMANN;
  bc.pS = NEUMANN;
  bc.pN = NEUMANN;
  bc.pB = NEUMANN;
  bc.pT = NEUMANN;
  bc.uW = NEUMANN;
  bc.uE = NEUMANN;
  bc.uS = NEUMANN;
  bc.uN = NEUMANN;
  bc.uB = NEUMANN;
  bc.uT = NEUMANN;
  bc.vW = NEUMANN;
  bc.vE = NEUMANN;
  bc.vS = NEUMANN;
  bc.vN = NEUMANN;
  bc.vB = NEUMANN;
  bc.vT = NEUMANN;
  bc.wW = NEUMANN;
  bc.wE = NEUMANN;
  bc.wS = NEUMANN;
  bc.wN = NEUMANN;
  bc.wB = NEUMANN;
  bc.wT = NEUMANN;
  printf("done.\n");

  // apply BC
  printf("    Running cuda_BC()...");
  cuda_dom_BC();
  printf("done.\n");

  // pull to host
  printf("    Pulling fields back to host...");
  cuda_dom_pull();
  printf("done.\n");

  // write computed solution
  rec_paraview_stepnum_out++;
  printf("\n    Writing computed solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // copy results and compute error
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p_n_o[C] = p[C];
        p_n_e[C] = p_n_o[C] - p_n_i[C];
        if(fabs(p_n_e[C]) > p_n_err_max) p_n_err_max = fabs(p_n_e[C]);
        if(fabs(p_n_e[C]) < p_n_err_min) p_n_err_min = fabs(p_n_e[C]);
        p[C] = p_n_e[C];
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_n_o[C] = u[C];
        u_n_e[C] = u_n_o[C] - u_n_i[C];
        if(fabs(u_n_e[C]) > u_n_err_max) u_n_err_max = fabs(u_n_e[C]);
        if(fabs(u_n_e[C]) < u_n_err_min) u_n_err_min = fabs(u_n_e[C]);
        u[C] = u_n_e[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_n_o[C] = v[C];
        v_n_e[C] = v_n_o[C] - v_n_i[C];
        if(fabs(v_n_e[C]) > v_n_err_max) v_n_err_max = fabs(v_n_e[C]);
        if(fabs(v_n_e[C]) < v_n_err_min) v_n_err_min = fabs(v_n_e[C]);
        v[C] = v_n_e[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_n_o[C] = w[C];
        w_n_e[C] = w_n_o[C] - w_n_i[C];
        if(fabs(w_n_e[C]) > w_n_err_max) w_n_err_max = fabs(w_n_e[C]);
        if(fabs(w_n_e[C]) < w_n_err_min) w_n_err_min = fabs(w_n_e[C]);
        w[C] = w_n_e[C];
      }
    }
  }

  // write error difference
  rec_paraview_stepnum_out++;
  printf("    Writing error difference to:  out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  printf("\n    Error summary:\n");
  printf("    Field component:      minimum error:     maximum error:\n");
  printf("          p              %12.3e       %12.3e\n",
    p_n_err_min, p_n_err_max);
  printf("          u              %12.3e       %12.3e\n",
    u_n_err_min, u_n_err_max);
  printf("          v              %12.3e       %12.3e\n",
    v_n_err_min, v_n_err_max);
  printf("          w              %12.3e       %12.3e\n",
    w_n_err_min, w_n_err_max);

  // clean up
  free(p_p_i);
  free(p_p_o);
  free(p_p_e);
  free(p_d_i);
  free(p_d_o);
  free(p_d_e);
  free(p_n_i);
  free(p_n_o);
  free(p_n_e);
  free(u_p_i);
  free(u_p_o);
  free(u_p_e);
  free(u_d_i);
  free(u_d_o);
  free(u_d_e);
  free(u_n_i);
  free(u_n_o);
  free(u_n_e);
  free(v_p_i);
  free(v_p_o);
  free(v_p_e);
  free(v_d_i);
  free(v_d_o);
  free(v_d_e);
  free(v_n_i);
  free(v_n_o);
  free(v_n_e);
  free(w_p_i);
  free(w_p_o);
  free(w_p_e);
  free(w_d_i);
  free(w_d_o);
  free(w_d_e);
  free(w_n_i);
  free(w_n_o);
  free(w_n_e);
}

void cuda_project_test(void)
{
  int i, j, k;  // iterators
  int C;    // cell locations

  real *u_a = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_c = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *u_e = (real*) malloc(Dom.Gfx.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfx.s3b * sizeof(real);
  real *v_a = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_c = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfy.s3b * sizeof(real);
  real *v_e = (real*) malloc(Dom.Gfy.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gyz.s3b * sizeof(real);
  real *w_a = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // expected solution
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_c = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // computed solution
  // cpumem += Dom.Gfz.s3b * sizeof(real);
  real *w_e = (real*) malloc(Dom.Gfz.s3b * sizeof(real)); // error difference
  // cpumem += Dom.Gfz.s3b * sizeof(real);

  // min and max error search
  real u_err_min = FLT_MAX;
  real u_err_max = FLT_MIN;
  real v_err_min = FLT_MAX;
  real v_err_max = FLT_MIN;
  real w_err_min = FLT_MAX;
  real w_err_max = FLT_MIN;

  printf("\nVelocity projection calculation validation:\n\n");
  printf("  u = exp(x), v = exp(y), w = exp(z), ");
  printf("p = exp(x) + exp(y) + exp(z)\n\n");

  // set up expected solution
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_a[C] = nu;
        u_a[C] -= 1. / rho_f;
        u_a[C] -= 2 * exp((i-1.5)*Dom.dx);
        u_a[C] -= exp((j-1.0)*Dom.dy);
        u_a[C] -= exp((k-1.0)*Dom.dz);
        u_a[C] *= dt * exp((i-1.5)*Dom.dx);
        u_a[C] += exp((i-1.5)*Dom.dx);
        u[C] = u_a[C];
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_a[C] = nu;
        v_a[C] -= 1. / rho_f;
        v_a[C] -= 2 * exp((j-1.5)*Dom.dy);
        v_a[C] -= exp((k-1.0)*Dom.dz);
        v_a[C] -= exp((i-1.0)*Dom.dx);
        v_a[C] *= dt * exp((j-1.5)*Dom.dy);
        v_a[C] += exp((j-1.5)*Dom.dy);
        v[C] = v_a[C];
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_a[C] = nu;
        w_a[C] -= 1. / rho_f;
        w_a[C] -= 2 * exp((k-1.5)*Dom.dz);
        w_a[C] -= exp((i-1.0)*Dom.dx);
        w_a[C] -= exp((j-1.0)*Dom.dy);
        w_a[C] *= dt * exp((k-1.5)*Dom.dz);
        w_a[C] += exp((k-1.5)*Dom.dz);
        w[C] = w_a[C];
      }
    }
  }

  // write expected solution
  rec_paraview_stepnum_out++;
  printf("  Writing expected solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // initialize input pressure and velocity fields
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p[C] = exp((i-1.0)*Dom.dx) + exp((j-1.0)*Dom.dy) + exp((k-1.0)*Dom.dz);
      }
    }
  }
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u[C] = exp((i-1.5)*Dom.dx);
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v[C] = exp((j-1.5)*Dom.dy);
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w[C] = exp((k-1.5)*Dom.dz);
      }
    }
  }

  // write initial fields
  rec_paraview_stepnum_out++;
  printf("  Writing initial fields to:    out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // push fields to device
  printf("\n  Pushing fields to devices...");
  cuda_dom_push();
  printf("done.\n");

  // call code to test
  printf("  Running cuda_U_star_2()...");
  cuda_U_star_2();
  cuda_project();
  printf("done.\n");

  // pull fields back to host
  printf("  Pulling fields back to host...");
  cuda_dom_pull();
  printf("done.\n");

  // write computed solution
  rec_paraview_stepnum_out++;
  printf("\n  Writing computed solution to: out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  // copy results and compute error
  for(k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(i = Dom.Gfx.is; i < Dom.Gfx.ie; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u_c[C] = u[C];
        u_e[C] = (u_c[C] - u_a[C]) / u_a[C];
        if(fabs(u_e[C]) > u_err_max) u_err_max = fabs(u_e[C]);
        if(fabs(u_e[C]) < u_err_min) u_err_min = fabs(u_e[C]);
        u[C] = u_e[C];
      }
    }
  }
  for(k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(j = Dom.Gfy.js; j < Dom.Gfy.je; j++) {
      for(i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v_c[C] = v[C];
        v_e[C] = (v_c[C] - v_a[C]) / v_a[C];
        if(fabs(v_e[C]) > v_err_max) v_err_max = fabs(v_e[C]);
        if(fabs(v_e[C]) < v_err_min) v_err_min = fabs(v_e[C]);
        v[C] = v_e[C];
      }
    }
  }
  for(k = Dom.Gfz.ks; k < Dom.Gfz.ke; k++) {
    for(j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w_c[C] = w[C];
        w_e[C] = (w_c[C] - w_a[C]) / w_a[C];
        if(fabs(w_e[C]) > w_err_max) w_err_max = fabs(w_e[C]);
        if(fabs(w_e[C]) < w_err_min) w_err_min = fabs(w_e[C]);
        w[C] = w_e[C];
      }
    }
  }

  // write error difference
  rec_paraview_stepnum_out++;
  printf("  Writing error difference to:  out_%d.pvtr...", rec_paraview_stepnum_out);
  out_VTK_ghost();
  printf("done.\n");

  printf("\n  Error summary:\n");
  printf("  Velocity component:     minimum error:     maximum error:\n");
  printf("          u              %12.3e       %12.3e\n",
    u_err_min, u_err_max);
  printf("          v              %12.3e       %12.3e\n",
    v_err_min, v_err_max);
  printf("          w              %12.3e       %12.3e\n\n",
    w_err_min, w_err_max);

  // clean up
  free(u_a);
  free(u_c);
  free(u_e);
  free(v_a);
  free(v_c);
  free(v_e);
  free(w_a);
  free(w_c);
  free(w_e);
}


