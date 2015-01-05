#include "bluebottle.h"
#include "scalar.h"
#include "point.h"

void init_VTK(void)
{
  char fname[FILE_NAME_SIZE];

  // open PVD file for writing
  sprintf(fname, "%sout.pvd", OUTPUT_DIR);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write PVD file header and footer
  fprintf(outfile, "<VTKFile type=\"Collection\">\n");
  fprintf(outfile, "<Collection>\n");
  fprintf(outfile, "</Collection>\n");
  fprintf(outfile, "</VTKFile>");

  // close the file
  fclose(outfile);
}

void out_VTK(void)
{
  char fname_pvd[FILE_NAME_SIZE]; // pvd filename
  char fname_pvtr[FILE_NAME_SIZE]; // pvtr filename
  char fname_vtp[FILE_NAME_SIZE]; // vtp filename
  char fnamenodes_vtp[FILE_NAME_SIZE]; // vtp filename

  sprintf(fname_pvd, "%sout.pvd", OUTPUT_DIR);
  sprintf(fname_pvtr, "out_%d.pvtr", rec_paraview_stepnum_out);
  sprintf(fname_vtp, "out_%d.vtp", rec_paraview_stepnum_out);
  sprintf(fnamenodes_vtp, "out_nodes_%d.vtp", rec_paraview_stepnum_out);

  FILE *pvdfile= fopen(fname_pvd, "r+");
  // moves back 2 lines from the end of the file (above the footer)
  fseek(pvdfile, -24, SEEK_END);

  fprintf(pvdfile, "<DataSet timestep=\"%e\" point=\"0\" file=\"%s\"/>\n",
    ttime, fname_pvtr);
  fprintf(pvdfile, "<DataSet timestep=\"%e\" point=\"1\" file=\"%s\"/>\n",
    ttime, fname_vtp);
  fprintf(pvdfile, "<DataSet timestep=\"%e\" point=\"2\" file=\"%s\"/>\n",
    ttime, fnamenodes_vtp);
  fprintf(pvdfile, "</Collection>\n");
  fprintf(pvdfile, "</VTKFile>");
  fclose(pvdfile);

  dom_out_VTK();
  point_out_VTK();
 // quadnodes_out_VTK();
}

void dom_out_VTK(void)
{
  int i, j, k, l; // iterators
  char fname[FILE_NAME_SIZE]; // output filename
  char fname_dom[FILE_NAME_SIZE]; // subdomain filename
  int C;  // cell center index
  int Cx;  // cell center index for interpolation
  int Cy;  // cell center index for interpolation
  int Cz;  // cell center index for interpolation

  // number of cells in a subdomain (length, start, end)
  int ncx_l, ncx_s, ncx_e;  
  int ncy_l, ncy_s, ncy_e;  
  int ncz_l, ncz_s, ncz_e;  

  sprintf(fname, "%sout_%d.pvtr", OUTPUT_DIR, rec_paraview_stepnum_out);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write Paraview pvtr file
  fprintf(outfile, "<VTKFile type=\"PRectilinearGrid\">\n");
  fprintf(outfile, "<PRectilinearGrid WholeExtent=");
  fprintf(outfile, "\"0 %d 0 %d 0 %d\" GhostLevel=\"0\">\n",
    Dom.xn, Dom.yn, Dom.zn);
  //fprintf(outfile, "<PCellData Scalars=\"p divU phase phase_shell\" Vectors=\"vel flag\">\n");
  fprintf(outfile, "<PCellData Scalars=\"p phase\" Vectors=\"vel\">\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"p\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"divU\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Int32\" Name=\"phase\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Int32\" Name=\"phase_shell\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"vel\"");
  fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"flag\"");
  //fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  fprintf(outfile, "</PCellData>\n");
  fprintf(outfile, "<PCoordinates>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"x\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"y\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"z\"/>\n");
  fprintf(outfile, "</PCoordinates>\n");
  for(l = 0; l < 6 * nsubdom; l += 6) {
    ncx_l = (dom[l].xe - dom[l].xs) / Dom.dx;
    ncx_s = (dom[l].xs - Dom.xs) / Dom.dx;
    ncx_e = ncx_s + ncx_l;
    ncy_l = (dom[l].ye - dom[l].ys) / Dom.dy;
    ncy_s = (dom[l].ys - Dom.ys) / Dom.dy;
    ncy_e = ncy_s + ncy_l;
    ncz_l = (dom[l].ze - dom[l].zs) / Dom.dz;
    ncz_s = (dom[l].zs - Dom.zs) / Dom.dz;
    ncz_e = ncz_s + ncz_l;
    fprintf(outfile, "<Piece Extent=\"");
    fprintf(outfile, "%d %d ", ncx_s, ncx_e);
    fprintf(outfile, "%d %d ", ncy_s, ncy_e);
    fprintf(outfile, "%d %d\" ", ncz_s, ncz_e);
    sprintf(fname_dom, "out_%d_%d.vtr", rec_paraview_stepnum_out, l/6);
    fprintf(outfile, "Source=\"%s\"/>\n", fname_dom);
  }
  fprintf(outfile, "</PRectilinearGrid>\n");
  fprintf(outfile, "</VTKFile>\n");
  fclose(outfile);

  // interpolate velocities to cell centers
  // cell-center working arrays
  real *uu = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *vv = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *ww = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  //real *flag_uu = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  //real *flag_vv = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  //real *flag_ww = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  for(k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
    for(j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
      for(i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
        // interpolate velocity
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        Cx = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        Cy = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        Cz = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        uu[C] = 0.5 * (u[Cx] + u[Cx+1]);
        vv[C] = 0.5 * (v[Cy] + v[Cy+Dom.Gfy.s1b]);
        ww[C] = 0.5 * (w[Cz] + w[Cz+Dom.Gfz.s2b]);
        // interpolate flags
        //flag_uu[C] = 0.5*(flag_u[Cx] + flag_u[Cx+1]);
        //flag_vv[C] = 0.5*(flag_v[Cy] + flag_v[Cy+Dom.Gfy.s1b]);
        //flag_ww[C] = 0.5*(flag_w[Cz] + flag_w[Cz+Dom.Gfz.s2b]);
      }
    }
  }

  // write each subdomain file
  for(l = 0; l < 6 * nsubdom; l += 6) {
    // number of cells in the subdomain
    ncx_l = (dom[l].xe - dom[l].xs) / Dom.dx;
    ncx_s = (dom[l].xs - Dom.xs) / Dom.dx;
    ncx_e = ncx_s + ncx_l;
    ncy_l = (dom[l].ye - dom[l].ys) / Dom.dy;
    ncy_s = (dom[l].ys - Dom.ys) / Dom.dy;
    ncy_e = ncy_s + ncy_l;
    ncz_l = (dom[l].ze - dom[l].zs) / Dom.dz;
    ncz_s = (dom[l].zs - Dom.zs) / Dom.dz;
    ncz_e = ncz_s + ncz_l;

    // open file for writing
    sprintf(fname, "%s/out_%d_%d.vtr", OUTPUT_DIR, rec_paraview_stepnum_out, l/6);
    FILE *outfile = fopen(fname, "w");
    if(outfile == NULL) {
      fprintf(stderr, "Could not open file %s\n", fname);
      exit(EXIT_FAILURE);
    }

    fprintf(outfile, "<VTKFile type=\"RectilinearGrid\">\n");
    fprintf(outfile, "<RectilinearGrid WholeExtent=");
    fprintf(outfile, "\"0 %d 0 %d 0 %d\" GhostLevel=\"0\">\n",
      Dom.xn, Dom.yn, Dom.zn);
    fprintf(outfile, "<Piece Extent=\"");
    fprintf(outfile, "%d %d ", ncx_s, ncx_e);
    fprintf(outfile, "%d %d ", ncy_s, ncy_e);
    fprintf(outfile, "%d %d\">\n", ncz_s, ncz_e);
    //fprintf(outfile, "<CellData Scalars=\"p divU phase phase_shell\" Vectors=\"vel flag\">\n");
    fprintf(outfile, "<CellData Scalars=\"p phase\" Vectors=\"vel\">\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"p\">\n");
    // write pressure for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e ", p[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    /*fprintf(outfile, "<DataArray type=\"Float32\" Name=\"divU\">\n");
    // write pressure for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e ", divU[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
*/
/*
    fprintf(outfile, "<DataArray type=\"Int32\" Name=\"phase\">\n");
    // write phase for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%d ", phase[C]);
        }
      }
    }
    fprintf(outfile, "\n");
*/
/*
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Int32\" Name=\"phase_shell\">\n");
    // write phase for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%d ", phase_shell[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
*/

    // write velocity vector
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"vel\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e %e %e ", uu[C], vv[C], ww[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");
    // write flag vector
    /*fprintf(outfile, "<DataArray type=\"Float32\" Name=\"flag\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
          fprintf(outfile, "%lf %lf %lf ", flag_uu[C], flag_vv[C], flag_ww[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");
*/
    fprintf(outfile, "</CellData>\n");

    fprintf(outfile, "<Coordinates>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"x\">\n");
    for(i = 0; i <= ncx_l; i++) {
      fprintf(outfile, "%lf ", dom[l].xs + Dom.dx * i);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"y\">\n");
    for(j = 0; j <= ncy_l; j++) {
      fprintf(outfile, "%lf ", dom[l].ys + Dom.dy * j);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"z\">\n");
    for(k = 0; k <= ncz_l; k++) {
      fprintf(outfile, "%lf ", dom[l].zs + Dom.dz * k);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "</Coordinates>\n");
    fprintf(outfile, "</Piece>\n");
    fprintf(outfile, "</RectilinearGrid>\n");
    fprintf(outfile, "</VTKFile>\n");
    fclose(outfile);
  }

  // clean up interpolated fields
  free(uu);
  free(vv);
  free(ww);
  //free(flag_uu);
  //free(flag_vv);
  //free(flag_ww);
}

void init_VTK_turb(void)
{
  char fname[FILE_NAME_SIZE];

  // open PVD file for writing
  sprintf(fname, "%sout_turb.pvd", OUTPUT_DIR);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write PVD file header and footer
  fprintf(outfile, "<VTKFile type=\"Collection\">\n");
  fprintf(outfile, "<Collection>\n");
  fprintf(outfile, "</Collection>\n");
  fprintf(outfile, "</VTKFile>");

  // close the file
  fclose(outfile);
}

void out_VTK_turb(void)
{
  char fname_pvd[FILE_NAME_SIZE]; // pvd filename
  char fname_pvtr[FILE_NAME_SIZE]; // pvtr filename

  sprintf(fname_pvd, "%sout_turb.pvd", OUTPUT_DIR);
  sprintf(fname_pvtr, "out_turb_%d.pvtr", rec_precursor_stepnum_out);

  FILE *pvdfile= fopen(fname_pvd, "r+");
  // moves back 2 lines from the end of the file (above the footer)
  fseek(pvdfile, -24, SEEK_END);

  fprintf(pvdfile, "<DataSet timestep=\"%e\" point=\"0\" file=\"%s\"/>\n",
    ttime, fname_pvtr);
  fprintf(pvdfile, "</Collection>\n");
  fprintf(pvdfile, "</VTKFile>");
  fclose(pvdfile);

  dom_out_VTK_turb();
}

void dom_out_VTK_turb(void)
{
  int i, j, k, l; // iterators
  char fname[FILE_NAME_SIZE]; // output filename
  char fname_dom[FILE_NAME_SIZE]; // subdomain filename
  int C;  // cell center index
  int Cx;  // cell center index for interpolation
  int Cy;  // cell center index for interpolation
  int Cz;  // cell center index for interpolation

  // number of cells in a subdomain (length, start, end)
  int ncx_l, ncx_s, ncx_e;  
  int ncy_l, ncy_s, ncy_e;  
  int ncz_l, ncz_s, ncz_e;  

  sprintf(fname, "%sout_turb_%d.pvtr", OUTPUT_DIR, rec_precursor_stepnum_out);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write Paraview pvtr file
  fprintf(outfile, "<VTKFile type=\"PRectilinearGrid\">\n");
  fprintf(outfile, "<PRectilinearGrid WholeExtent=");
  fprintf(outfile, "\"0 %d 0 %d 0 %d\" GhostLevel=\"0\">\n",
    Dom.xn, Dom.yn, Dom.zn);
  //fprintf(outfile, "<PCellData Scalars=\"p divU phase phase_shell\" Vectors=\"vel flag\">\n");
  fprintf(outfile, "<PCellData Scalars=\"p\" Vectors=\"vel\">\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"p\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"divU\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Int32\" Name=\"phase\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Int32\" Name=\"phase_shell\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"vel\"");
  fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"flag\"");
  //fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  fprintf(outfile, "</PCellData>\n");
  fprintf(outfile, "<PCoordinates>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"x\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"y\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"z\"/>\n");
  fprintf(outfile, "</PCoordinates>\n");
  for(l = 0; l < 6 * nsubdom; l += 6) {
    ncx_l = (dom[l].xe - dom[l].xs) / Dom.dx;
    ncx_s = (dom[l].xs - Dom.xs) / Dom.dx;
    ncx_e = ncx_s + ncx_l;
    ncy_l = (dom[l].ye - dom[l].ys) / Dom.dy;
    ncy_s = (dom[l].ys - Dom.ys) / Dom.dy;
    ncy_e = ncy_s + ncy_l;
    ncz_l = (dom[l].ze - dom[l].zs) / Dom.dz;
    ncz_s = (dom[l].zs - Dom.zs) / Dom.dz;
    ncz_e = ncz_s + ncz_l;
    fprintf(outfile, "<Piece Extent=\"");
    fprintf(outfile, "%d %d ", ncx_s, ncx_e);
    fprintf(outfile, "%d %d ", ncy_s, ncy_e);
    fprintf(outfile, "%d %d\" ", ncz_s, ncz_e);
    sprintf(fname_dom, "out_turb_%d_%d.vtr", rec_precursor_stepnum_out, l/6);
    fprintf(outfile, "Source=\"%s\"/>\n", fname_dom);
  }
  fprintf(outfile, "</PRectilinearGrid>\n");
  fprintf(outfile, "</VTKFile>\n");
  fclose(outfile);

  // interpolate velocities to cell centers
  // cell-center working arrays
  real *uu = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *vv = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *ww = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  //real *flag_uu = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  //real *flag_vv = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  //real *flag_ww = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  for(k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
    for(j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
      for(i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
        // interpolate velocity
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        Cx = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        Cy = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        Cz = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        uu[C] = 0.5 * (u[Cx] + u[Cx+1]);
        vv[C] = 0.5 * (v[Cy] + v[Cy+Dom.Gfy.s1b]);
        ww[C] = 0.5 * (w[Cz] + w[Cz+Dom.Gfz.s2b]);
        // interpolate flags
        //flag_uu[C] = 0.5*(flag_u[Cx] + flag_u[Cx+1]);
        //flag_vv[C] = 0.5*(flag_v[Cy] + flag_v[Cy+Dom.Gfy.s1b]);
        //flag_ww[C] = 0.5*(flag_w[Cz] + flag_w[Cz+Dom.Gfz.s2b]);
      }
    }
  }

  // write each subdomain file
  for(l = 0; l < 6 * nsubdom; l += 6) {
    // number of cells in the subdomain
    ncx_l = (dom[l].xe - dom[l].xs) / Dom.dx;
    ncx_s = (dom[l].xs - Dom.xs) / Dom.dx;
    ncx_e = ncx_s + ncx_l;
    ncy_l = (dom[l].ye - dom[l].ys) / Dom.dy;
    ncy_s = (dom[l].ys - Dom.ys) / Dom.dy;
    ncy_e = ncy_s + ncy_l;
    ncz_l = (dom[l].ze - dom[l].zs) / Dom.dz;
    ncz_s = (dom[l].zs - Dom.zs) / Dom.dz;
    ncz_e = ncz_s + ncz_l;

    // open file for writing
    sprintf(fname, "%s/out_turb_%d_%d.vtr", OUTPUT_DIR, rec_precursor_stepnum_out, l/6);
    FILE *outfile = fopen(fname, "w");
    if(outfile == NULL) {
      fprintf(stderr, "Could not open file %s\n", fname);
      exit(EXIT_FAILURE);
    }

    fprintf(outfile, "<VTKFile type=\"RectilinearGrid\">\n");
    fprintf(outfile, "<RectilinearGrid WholeExtent=");
    fprintf(outfile, "\"0 %d 0 %d 0 %d\" GhostLevel=\"0\">\n",
      Dom.xn, Dom.yn, Dom.zn);
    fprintf(outfile, "<Piece Extent=\"");
    fprintf(outfile, "%d %d ", ncx_s, ncx_e);
    fprintf(outfile, "%d %d ", ncy_s, ncy_e);
    fprintf(outfile, "%d %d\">\n", ncz_s, ncz_e);
    //fprintf(outfile, "<CellData Scalars=\"p divU phase phase_shell\" Vectors=\"vel flag\">\n");
    fprintf(outfile, "<CellData Scalars=\"p\" Vectors=\"vel\">\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"p\">\n");
    // write pressure for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e ", p[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    //fprintf(outfile, "<DataArray type=\"Float32\" Name=\"divU\">\n");
    // write pressure for this subdomain
    /*for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e ", divU[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
*/
    /*fprintf(outfile, "<DataArray type=\"Int32\" Name=\"phase\">\n");
    // write phase for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%d ", phase[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Int32\" Name=\"phase_shell\">\n");
    // write phase for this subdomain
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%d ", phase_shell[C]);
        }
      }
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
*/

    // write velocity vector
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"vel\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e %e %e ", uu[C], vv[C], ww[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");
    // write flag vector
    /*fprintf(outfile, "<DataArray type=\"Float32\" Name=\"flag\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ks; k < ncz_e + Dom.Gcc.ks; k++) {
      for(j = ncy_s + Dom.Gcc.js; j < ncy_e + Dom.Gcc.js; j++) {
        for(i = ncx_s + Dom.Gcc.is; i < ncx_e + Dom.Gcc.is; i++) {
          C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
          fprintf(outfile, "%lf %lf %lf ", flag_uu[C], flag_vv[C], flag_ww[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");
*/
    fprintf(outfile, "</CellData>\n");

    fprintf(outfile, "<Coordinates>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"x\">\n");
    for(i = 0; i <= ncx_l; i++) {
      fprintf(outfile, "%lf ", dom[l].xs + Dom.dx * i);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"y\">\n");
    for(j = 0; j <= ncy_l; j++) {
      fprintf(outfile, "%lf ", dom[l].ys + Dom.dy * j);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"z\">\n");
    for(k = 0; k <= ncz_l; k++) {
      fprintf(outfile, "%lf ", dom[l].zs + Dom.dz * k);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "</Coordinates>\n");
    fprintf(outfile, "</Piece>\n");
    fprintf(outfile, "</RectilinearGrid>\n");
    fprintf(outfile, "</VTKFile>\n");
    fclose(outfile);
  }

  // clean up interpolated fields
  free(uu);
  free(vv);
  free(ww);
  //free(flag_uu);
  //free(flag_vv);
  //free(flag_ww);
}

void init_VTK_ghost(void)
{
  char fname[FILE_NAME_SIZE];

  // open PVD file for writing
  sprintf(fname, "%sout_ghost.pvd", OUTPUT_DIR);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write PVD file header and footer
  fprintf(outfile, "<VTKFile type=\"Collection\">\n");
  fprintf(outfile, "<Collection>\n");
  fprintf(outfile, "</Collection>\n");
  fprintf(outfile, "</VTKFile>");

  // close the file
  fclose(outfile);
}

void out_VTK_ghost(void)
{
  char fname_pvd[FILE_NAME_SIZE]; // pvd filename
  char fname_pvtr[FILE_NAME_SIZE]; // pvtr filename
  char fname_vtp[FILE_NAME_SIZE]; // vtp filename

  sprintf(fname_pvd, "%sout_ghost.pvd", OUTPUT_DIR);
  sprintf(fname_pvtr, "out_ghost_%d.pvtr", rec_paraview_stepnum_out);
  sprintf(fname_vtp, "out_ghost_%d.vtp", rec_paraview_stepnum_out);

  FILE *pvdfile= fopen(fname_pvd, "r+");
  // moves back 2 lines from the end of the file (above the footer)
  fseek(pvdfile, -24, SEEK_END);

  fprintf(pvdfile, "<DataSet timestep=\"%e\" point=\"0\" file=\"%s\"/>\n",
    ttime, fname_pvtr);
  fprintf(pvdfile, "<DataSet timestep=\"%e\" point=\"1\" file=\"%s\"/>\n",
    ttime, fname_vtp);
  fprintf(pvdfile, "</Collection>\n");
  fprintf(pvdfile, "</VTKFile>");
  fclose(pvdfile);

  dom_out_VTK_ghost();
  point_out_VTK();
}

void dom_out_VTK_ghost(void)
{
  int i, j, k, l; // iterators
  char fname[FILE_NAME_SIZE]; // output filename
  char fname_dom[FILE_NAME_SIZE]; // subdomain filename
  int C;  // cell center index
  int Cx;  // cell center index for interpolation
  int Cy;  // cell center index for interpolation
  int Cz;  // cell center index for interpolation

  // number of cells in a subdomain (length, start, end)
  int ncx_l, ncx_s, ncx_e;  
  int ncy_l, ncy_s, ncy_e;  
  int ncz_l, ncz_s, ncz_e;  

  sprintf(fname, "%sout_ghost_%d.pvtr", OUTPUT_DIR, rec_paraview_stepnum_out);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // write Paraview pvtr file
  fprintf(outfile, "<VTKFile type=\"PRectilinearGrid\">\n");
  fprintf(outfile, "<PRectilinearGrid WholeExtent=");
  fprintf(outfile, "\"0 %d 0 %d 0 %d\" GhostLevel=\"0\">\n",
    Dom.xn+2*DOM_BUF, Dom.yn+2*DOM_BUF, Dom.zn+2*DOM_BUF);
  //fprintf(outfile, "<PCellData Scalars=\"p divU phase phase_shell\" Vectors=\"vel vel_star flag\">\n");
  fprintf(outfile, "<PCellData Scalars=\"p phase phase_shell\" Vectors=\"vel vel_star flag\">\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"p\"/>\n");
  //fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"divU\"/>\n");
 /*
  fprintf(outfile, "<PDataArray type=\"Int32\" Name=\"phase\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Int32\" Name=\"phase_shell\"/>\n");
 */
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"vel\"");
  fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"vel_star\"");
  fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"flag\"");
  fprintf(outfile, " NumberOfComponents=\"3\"/>\n");
  fprintf(outfile, "</PCellData>\n");
  fprintf(outfile, "<PCoordinates>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"x\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"y\"/>\n");
  fprintf(outfile, "<PDataArray type=\"Float32\" Name=\"z\"/>\n");
  fprintf(outfile, "</PCoordinates>\n");
  for(l = 0; l < 6 * nsubdom; l += 6) {
    ncx_l = (dom[l].xe - dom[l].xs) / Dom.dx + 2*DOM_BUF;
    ncx_s = (dom[l].xs - Dom.xs) / Dom.dx;
    ncx_e = ncx_s + ncx_l;
    ncy_l = (dom[l].ye - dom[l].ys) / Dom.dy + 2*DOM_BUF;
    ncy_s = (dom[l].ys - Dom.ys) / Dom.dy;
    ncy_e = ncy_s + ncy_l;
    ncz_l = (dom[l].ze - dom[l].zs) / Dom.dz + 2*DOM_BUF;
    ncz_s = (dom[l].zs - Dom.zs) / Dom.dz;
    ncz_e = ncz_s + ncz_l;
    fprintf(outfile, "<Piece Extent=\"");
    fprintf(outfile, "%d %d ", ncx_s, ncx_e);
    fprintf(outfile, "%d %d ", ncy_s, ncy_e);
    fprintf(outfile, "%d %d\" ", ncz_s, ncz_e);
    sprintf(fname_dom, "out_ghost_%d_%d.vtr", rec_paraview_stepnum_out, l/6);
    fprintf(outfile, "Source=\"%s\"/>\n", fname_dom);
  }
  fprintf(outfile, "</PRectilinearGrid>\n");
  fprintf(outfile, "</VTKFile>\n");
  fclose(outfile);

  // interpolate velocities to cell centers
  // cell-center working arrays
  real *uu = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *vv = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *ww = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *uu_star = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *vv_star = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *ww_star = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *flag_uu = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *flag_vv = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  real *flag_ww = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  // cpumem += Dom.Gcc.s3b * sizeof(real);
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        // interpolate velocity
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        Cx = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        Cy = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        Cz = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        uu[C] = 0.5 * (u[Cx] + u[Cx+1]);
        vv[C] = 0.5 * (v[Cy] + v[Cy+Dom.Gfy.s1b]);
        ww[C] = 0.5 * (w[Cz] + w[Cz+Dom.Gfz.s2b]);
        uu_star[C] = 0.5 * (u_star[Cx] + u_star[Cx+1]);
        vv_star[C] = 0.5 * (v_star[Cy] + v_star[Cy+Dom.Gfy.s1b]);
        ww_star[C] = 0.5 * (w_star[Cz] + w_star[Cz+Dom.Gfz.s2b]);
        // interpolate flags
        flag_uu[C] = 0.5*(flag_u[Cx] + flag_u[Cx+1]);
        flag_vv[C] = 0.5*(flag_v[Cy] + flag_v[Cy+Dom.Gfy.s1b]);
        flag_ww[C] = 0.5*(flag_w[Cz] + flag_w[Cz+Dom.Gfz.s2b]);
      }
    }
  }
  for(k = Dom.Gcc.ks; k < Dom.Gcc.ke; k++) {
    for(j = Dom.Gcc.js; j < Dom.Gcc.je; j++) {
      for(i = Dom.Gcc.is; i < Dom.Gcc.ie; i++) {
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        Cx = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        Cy = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        Cz = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        flag_uu[C] = 0.5*(flag_u[Cx] + flag_u[Cx+1]);
        flag_vv[C] = 0.5*(flag_v[Cy] + flag_v[Cy+Dom.Gfy.s1b]);
        flag_ww[C] = 0.5*(flag_w[Cz] + flag_w[Cz+Dom.Gfz.s2b]);
      }
    }
  }

  // write each subdomain file
  for(l = 0; l < 6 * nsubdom; l += 6) {
    // number of cells in the subdomain
    ncx_l = (dom[l].xe - dom[l].xs) / Dom.dx + 2*DOM_BUF;
    ncx_s = (dom[l].xs - Dom.xs) / Dom.dx;
    ncx_e = ncx_s + ncx_l;
    ncy_l = (dom[l].ye - dom[l].ys) / Dom.dy + 2*DOM_BUF;
    ncy_s = (dom[l].ys - Dom.ys) / Dom.dy;
    ncy_e = ncy_s + ncy_l;
    ncz_l = (dom[l].ze - dom[l].zs) / Dom.dz + 2*DOM_BUF;
    ncz_s = (dom[l].zs - Dom.zs) / Dom.dz;
    ncz_e = ncz_s + ncz_l;

    // open file for writing
    sprintf(fname, "%s/out_ghost_%d_%d.vtr", OUTPUT_DIR, rec_paraview_stepnum_out, l/6);
    FILE *outfile = fopen(fname, "w");
    if(outfile == NULL) {
      fprintf(stderr, "Could not open file %s\n", fname);
      exit(EXIT_FAILURE);
    }

    fprintf(outfile, "<VTKFile type=\"RectilinearGrid\">\n");
    fprintf(outfile, "<RectilinearGrid WholeExtent=");
    fprintf(outfile, "\"0 %d 0 %d 0 %d\" GhostLevel=\"0\">\n",
      Dom.xn+2*DOM_BUF, Dom.yn+2*DOM_BUF, Dom.zn+2*DOM_BUF);
    fprintf(outfile, "<Piece Extent=\"");
    fprintf(outfile, "%d %d ", ncx_s, ncx_e);
    fprintf(outfile, "%d %d ", ncy_s, ncy_e);
    fprintf(outfile, "%d %d\">\n", ncz_s, ncz_e);
    //fprintf(outfile, "<CellData Scalars=\"p divU phase phase_shell\" Vectors=\"vel vel_star flag\">\n");
    fprintf(outfile, "<CellData Scalars=\"p phase phase_shell\" Vectors=\"vel vel_star flag\">\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"p\">\n");
    // write pressure for this subdomain
    for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e ", p[C]);
        }
      }
    }
    //fprintf(outfile, "\n");
    //fprintf(outfile, "</DataArray>\n");
    //fprintf(outfile, "<DataArray type=\"Float32\" Name=\"divU\">\n");
    // write pressure for this subdomain
    /*for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          if(i >= Dom.Gcc.is && i < Dom.Gcc.ie
            && j >= Dom.Gcc.js && j < Dom.Gcc.je
            && k >= Dom.Gcc.ks && k < Dom.Gcc.ke) {
            fprintf(outfile, "%1.16e ", divU[C]);
          } else {
            fprintf(outfile, "%1.16e ", 0.);
          }
        }
      }
    }
*/
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Int32\" Name=\"phase\">\n");
/*  
  // write phase for this subdomain
    for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%d ", phase[C]);
        }
      }
    }
    fprintf(outfile, "\n");
  */
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Int32\" Name=\"phase_shell\">\n");
/*  
  // write phase for this subdomain
    for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%d ", phase_shell[C]);
        }
      }
    }
    fprintf(outfile, "\n");
  */
    fprintf(outfile, "</DataArray>\n");
    // write velocity vector
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"vel\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e %e %e ", uu[C], vv[C], ww[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");
    // write velocity star vector
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"vel_star\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%e %e %e ", uu_star[C], vv_star[C], ww_star[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");

    // write flag vector
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"flag\"");
    fprintf(outfile, " NumberOfComponents=\"3\" format=\"ascii\">\n");
    for(k = ncz_s + Dom.Gcc.ksb; k < ncz_e + Dom.Gcc.ksb; k++) {
      for(j = ncy_s + Dom.Gcc.jsb; j < ncy_e + Dom.Gcc.jsb; j++) {
        for(i = ncx_s + Dom.Gcc.isb; i < ncx_e + Dom.Gcc.isb; i++) {
          C = i + j * Dom.Gcc.s1b + k * Dom.Gcc.s2b;
          fprintf(outfile, "%lf %lf %lf ", flag_uu[C], flag_vv[C], flag_ww[C]);
        }
      }
    }
    fprintf(outfile, "\n</DataArray>\n");
    fprintf(outfile, "</CellData>\n");

    fprintf(outfile, "<Coordinates>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"x\">\n");
    for(i = -1; i < ncx_l; i++) {
      fprintf(outfile, "%lf ", dom[l].xs + Dom.dx * i);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"y\">\n");
    for(j = -1; j < ncy_l; j++) {
      fprintf(outfile, "%lf ", dom[l].ys + Dom.dy * j);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "<DataArray type=\"Float32\" Name=\"z\">\n");
    for(k = -1; k < ncz_l; k++) {
      fprintf(outfile, "%lf ", dom[l].zs + Dom.dz * k);
    }
    fprintf(outfile, "\n");
    fprintf(outfile, "</DataArray>\n");
    fprintf(outfile, "</Coordinates>\n");
    fprintf(outfile, "</Piece>\n");
    fprintf(outfile, "</RectilinearGrid>\n");
    fprintf(outfile, "</VTKFile>\n");
    fclose(outfile);
  }

  // clean up interpolated fields
  free(uu);
  free(vv);
  free(ww);
  free(uu_star);
  free(vv_star);
  free(ww_star);
  free(flag_uu);
  free(flag_vv);
  free(flag_ww);
}

void point_out_VTK(void)
{
  int i; // iterator
  char fname[FILE_NAME_SIZE]; // output filename
  int npoints_plot = 0;    // the number of point_particles to plot; may be greater
                          // may be greater than npoints if point_particles straddle

  // open file for writing
  sprintf(fname, "%s/out_%d.vtp", OUTPUT_DIR, rec_paraview_stepnum_out);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // go through all point_particles and check if any straddle boundaries
  npoints_plot = npoints;
  for(i = 0; i < npoints; i++) {
    int straddle_x = 0;
    int straddle_y = 0;
    int straddle_z = 0;
    // figure out which boudaries this point_particle straddles
    if(points[i].x < (Dom.xs + points[i].r)) straddle_x = 1;
    if(points[i].x > (Dom.xe - points[i].r)) straddle_x = 1;
    if(points[i].y < (Dom.ys + points[i].r)) straddle_y = 1;
    if(points[i].y > (Dom.ye - points[i].r)) straddle_y = 1;
    if(points[i].z < (Dom.zs + points[i].r)) straddle_z = 1;
    if(points[i].z > (Dom.ze - points[i].r)) straddle_z = 1;

    // now add the appropriate number of point_particles to plot
    if(straddle_x) npoints_plot += 1;
    if(straddle_y) npoints_plot += 1;
    if(straddle_z) npoints_plot += 1;
    if(straddle_x && straddle_y) npoints_plot += 1;
    if(straddle_x && straddle_z) npoints_plot += 1;
    if(straddle_y && straddle_z) npoints_plot += 1;
    if(straddle_x && straddle_y && straddle_z) npoints_plot += 1;
  }

  // create a temporary points list containing virtual point_particles
  point_struct *points_virt = (point_struct*) malloc(npoints_plot
    * sizeof(point_struct));
  // cpumem += npoints_plot * sizeof(point_struct);
  int *ind = (int*) malloc(npoints_plot * sizeof(int));
  // cpumem += npoints_plot * sizeof(int);

  int j = 0;  // virtual point_particle counter
  for(i = 0; i < npoints; i++) {
    // these take care of all of the point_particles stradding boundaries
    if(points[i].x < (Dom.xs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].y < (Dom.ys + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
 /*
  
*/
      j++;
    } if(points[i].y > (Dom.ye - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    }
    // if a point_particle straddles two boundaries, a fourth point_particle is needed
    if(points[i].x < (Dom.xs + points[i].r)
      && points[i].y < (Dom.ys + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].x < (Dom.xs + points[i].r)
      && points[i].y > (Dom.ye - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].y < (Dom.ys + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].y > (Dom.ye - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].z < (Dom.zs + points[i].r)
      && points[i].y < (Dom.ys + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].z < (Dom.zs + points[i].r)
      && points[i].y > (Dom.ye - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].z > (Dom.ze - points[i].r)
      && points[i].y < (Dom.ys + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].z > (Dom.ze - points[i].r)
      && points[i].y > (Dom.ye - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].x < (Dom.xs + points[i].r)
      && points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].x < (Dom.xs + points[i].r)
      && points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
/*
      
      
*/
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
/*
   
*/
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
   
      j++;
    }

    // if a point_particle straddles all three boundaries, an eighth is needed
    if(points[i].x < (Dom.xs + points[i].r)
      && points[i].y < (Dom.ys + points[i].r)
      && points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
   
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].y < (Dom.ys + points[i].r)
      && points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
   
      j++;
    } if(points[i].x < (Dom.xs + points[i].r)
      && points[i].y > (Dom.ye - points[i].r)
      && points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
   
      j++;
    } if(points[i].x < (Dom.xs + points[i].r)
      && points[i].y < (Dom.ys + points[i].r)
      && points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
   
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].y > (Dom.ye - points[i].r)
      && points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
   
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].y < (Dom.ys + points[i].r)
      && points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
   
      j++;
    } if(points[i].x < (Dom.xs + points[i].r)
      && points[i].y > (Dom.ye - points[i].r)
      && points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
   
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].y > (Dom.ye - points[i].r)
      && points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
   
      j++;
    }

    ind[j] = i;
    points_virt[j].r = points[i].r;
      points_virt[j].rho = points[i].rho;
    points_virt[j].x = points[i].x;
    points_virt[j].y = points[i].y;
    points_virt[j].z = points[i].z;
    points_virt[j].u = points[i].u;
    points_virt[j].v = points[i].v;
    points_virt[j].w = points[i].w;
    points_virt[j].udot = points[i].udot;
    points_virt[j].vdot = points[i].vdot;
    points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
    points_virt[j].Fx = points[i].Fx;
    points_virt[j].Fy = points[i].Fy;
    points_virt[j].Fz = points[i].Fz;
    points_virt[j].Lx = points[i].Lx;
    points_virt[j].Ly = points[i].Ly;
    points_virt[j].Lz = points[i].Lz;
      points_virt[j].iFx = points[i].iFx;
      points_virt[j].iFy = points[i].iFy;
      points_virt[j].iFz = points[i].iFz;
      points_virt[j].iLx = points[i].iLx;
      points_virt[j].iLy = points[i].iLy;
      points_virt[j].iLz = points[i].iLz;
   
    j++;
  }

  fprintf(outfile, "<VTKFile type=\"PolyData\">\n");
  fprintf(outfile, "<PolyData>\n");
  fprintf(outfile, "<Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" ",
    npoints_plot);
  fprintf(outfile, "NumberOfLines=\"0\" NumberOfStrips=\"0\" ");
  fprintf(outfile, "NumberOfPolys=\"0\">\n");
  fprintf(outfile, "<Points>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, "format=\"ascii\">\n");

  // write the locations of the point_particle centers
  for(i = 0; i < npoints_plot; i++) {
    fprintf(outfile, "%f %f %f ", points_virt[i].x, points_virt[i].y,
      points_virt[i].z);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "</Points>\n");
  fprintf(outfile, "<PointData Scalars=\"n r\" Vectors=\"pvel pacc");
  fprintf(outfile, " pax pay paz pomega pomegadot F L Fh Fk Fi\">\n");
  fprintf(outfile, "<DataArray type=\"Int32\" Name=\"n\" ");
  fprintf(outfile, "format=\"ascii\">\n");

  // write index of each point_particle
  for(i = 0; i < npoints_plot; i++) {
    fprintf(outfile, "%d ", ind[i]);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" Name=\"r\" ");
  fprintf(outfile, "format=\"ascii\">\n");

  // write radius of each point_particle
  for(i = 0; i < npoints_plot; i++) {
    fprintf(outfile, "%f ", points_virt[i].r);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"pvel\" format=\"ascii\">\n");

  // write velocity of each point_particle
  for(i = 0; i < npoints_plot; i++) {
    fprintf(outfile, "%f %f %f ", points_virt[i].u, points_virt[i].v,
      points_virt[i].w);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"pacc\" format=\"ascii\">\n");

  // write acceleration of each point_particle
  for(i = 0; i < npoints_plot; i++) {
    fprintf(outfile, "%f %f %f ", points_virt[i].udot, points_virt[i].vdot,
      points_virt[i].wdot);
  }
/*
  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"pax\" format=\"ascii\">\n");

  // write angular position of each point_particle
  for(i = 0; i < npoints_plot; i++) {
    fprintf(outfile, "%f %f %f ", points_virt[i].axx, points_virt[i].axy,
      points_virt[i].axz);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"pay\" format=\"ascii\">\n");

  // write angular position of each point_particle
  for(i = 0; i < npoints_plot; i++) {
    fprintf(outfile, "%f %f %f ", points_virt[i].ayx, points_virt[i].ayy,
      points_virt[i].ayz);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"paz\" format=\"ascii\">\n");

  // write angular position of each point_particle
  for(i = 0; i < npoints_plot; i++) {
    fprintf(outfile, "%f %f %f ", points_virt[i].azx, points_virt[i].azy,
      points_virt[i].azz);
  }
*/
  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"pomega\" format=\"ascii\">\n");

  // write angular velocity of each point_particle
  for(i = 0; i < npoints_plot; i++) {
    fprintf(outfile, "%f %f %f ", points_virt[i].ox, points_virt[i].oy,
      points_virt[i].oz);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"pomegadot\" format=\"ascii\">\n");

  // write angular acceleration of each point_particle
  for(i = 0; i < npoints_plot; i++) {
    fprintf(outfile, "%f %f %f ", points_virt[i].oxdot, points_virt[i].oydot,
      points_virt[i].ozdot);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"F\" format=\"ascii\">\n");

  // write linear force on each point_particle
  for(i = 0; i < npoints_plot; i++) {
    real mass = 4./3.*PI*(points_virt[i].rho-rho_f)
      *points_virt[i].r*points_virt[i].r*points_virt[i].r;
    fprintf(outfile, "%f %f %f ",
      points_virt[i].Fx + points_virt[i].iFx + mass*g.x,
      points_virt[i].Fy + points_virt[i].iFy + mass*g.y,
      points_virt[i].Fz + points_virt[i].iFz + mass*g.z);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"L\" format=\"ascii\">\n");

  // write angular moment on each point_particle
  for(i = 0; i < npoints_plot; i++) {
    fprintf(outfile, "%f %f %f ",
      points_virt[i].Lx + points_virt[i].iLx,
      points_virt[i].Ly + points_virt[i].iLy,
      points_virt[i].Lz + points_virt[i].iLz);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"Fh\" format=\"ascii\">\n");

  // write linear force on each point_particle
  for(i = 0; i < npoints_plot; i++) {
    fprintf(outfile, "%f %f %f ",
      points_virt[i].Fx,
      points_virt[i].Fy,
      points_virt[i].Fz);
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"Fk\" format=\"ascii\">\n");


  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, " Name=\"Fi\" format=\"ascii\">\n");

  // write linear force on each point_particle
  for(i = 0; i < npoints_plot; i++) {
    fprintf(outfile, "%f %f %f ",
      points_virt[i].iFx,
      points_virt[i].iFy,
      points_virt[i].iFz);
  }

  fprintf(outfile, "\n</DataArray>\n");

  fprintf(outfile, "</PointData>\n");
  fprintf(outfile, "</Piece>\n");
  fprintf(outfile, "</PolyData>\n");
  fprintf(outfile, "</VTKFile>\n");

  // clean up
  free(points_virt);
  free(ind);

  fclose(outfile);
}
/*
void quadnodes_out_VTK(void)
{
  int i, j; // iterator
  char fname[FILE_NAME_SIZE]; // output filename
  int npoints_plot = 0;    // the number of point_particles to plot; may be greater
                          // may be greater than npoints if point_particles straddle

  // create quadratue node locations
  real PI14 = 0.25 * PI;
  real PI12 = 0.5 * PI;
  real PI34 = 0.75 * PI;
  real PI54 = 1.25 * PI;
  real PI32 = 1.5 * PI;
  real PI74 = 1.75 * PI;
  real alph1 = 0.955316618124509; //54.736
  real alph2 = 2.186276035465284; //125.264

  real a1_t[6] = {PI12, PI12, PI12, PI12, 0.+DIV_ST, PI-DIV_ST};
  real a1_p[6] = {0., PI12, PI, PI32, 0., 0.};
  real a2_t[12] = {PI12, PI12, PI12, PI12,
                   PI14, PI14, PI14, PI14,
                   PI34, PI34, PI34, PI34};
  real a2_p[12] = {PI14, PI34, PI54, PI74,
                   0., PI12, PI, PI32,
                   0., PI12, PI, PI32};
  real a3_t[8] = {alph1, alph1, alph1, alph1,
                  alph2, alph2, alph2, alph2};
  real a3_p[8] = {PI14, PI34, PI54, PI74,
                  PI14, PI34, PI54, PI74};

  // put all quadrature nodes together for interpolation
  real node_t[NNODES];
  real node_p[NNODES];
  for(i = 0; i < 6; i++) {
    node_t[i] = a1_t[i];
    node_p[i] = a1_p[i];
  }
  for(i = 0; i < 12; i++) {
    node_t[6+i] = a2_t[i];
    node_p[6+i] = a2_p[i];
  }
  for(i = 0; i < 8; i++) {
    node_t[18+i] = a3_t[i];
    node_p[18+i] = a3_p[i];
  }


  // open file for writing
  sprintf(fname, "%s/out_nodes_%d.vtp", OUTPUT_DIR, rec_paraview_stepnum_out);
  FILE *outfile = fopen(fname, "w");
  if(outfile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // go through all point_particles and check if any straddle boundaries
  npoints_plot = npoints;
  for(i = 0; i < npoints; i++) {
    int straddle_x = 0;
    int straddle_y = 0;
    int straddle_z = 0;
    // figure out which boudaries this point_particle straddles
    if(points[i].x < (Dom.xs + points[i].r)) straddle_x = 1;
    if(points[i].x > (Dom.xe - points[i].r)) straddle_x = 1;
    if(points[i].y < (Dom.ys + points[i].r)) straddle_y = 1;
    if(points[i].y > (Dom.ye - points[i].r)) straddle_y = 1;
    if(points[i].z < (Dom.zs + points[i].r)) straddle_z = 1;
    if(points[i].z > (Dom.ze - points[i].r)) straddle_z = 1;

    // now add the appropriate number of point_particles to plot
    if(straddle_x) npoints_plot += 1;
    if(straddle_y) npoints_plot += 1;
    if(straddle_z) npoints_plot += 1;
    if(straddle_x && straddle_y) npoints_plot += 1;
    if(straddle_x && straddle_z) npoints_plot += 1;
    if(straddle_y && straddle_z) npoints_plot += 1;
    if(straddle_x && straddle_y && straddle_z) npoints_plot += 1;
  }

  // create a temporary points list containing virtual point_particles
  point_struct *points_virt = (point_struct*) malloc(npoints_plot
    * sizeof(point_struct));
  // cpumem += npoints_plot * sizeof(point_struct);
  int *ind = (int*) malloc(npoints_plot * sizeof(int));
  // cpumem += npoints_plot * sizeof(int);

  j = 0;  // virtual point_particle counter
  for(i = 0; i < npoints; i++) {
    // these take care of all of the point_particles stradding boundaries
    if(points[i].x < (Dom.xs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].y < (Dom.ys + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].y > (Dom.ye - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    }
    // if a point_particle straddles two boundaries, a fourth point_particle is needed
    if(points[i].x < (Dom.xs + points[i].r)
      && points[i].y < (Dom.ys + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x < (Dom.xs + points[i].r)
      && points[i].y > (Dom.ye - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].y < (Dom.ys + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].y > (Dom.ye - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].z < (Dom.zs + points[i].r)
      && points[i].y < (Dom.ys + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].z < (Dom.zs + points[i].r)
      && points[i].y > (Dom.ye - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].z > (Dom.ze - points[i].r)
      && points[i].y < (Dom.ys + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].z > (Dom.ze - points[i].r)
      && points[i].y > (Dom.ye - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x < (Dom.xs + points[i].r)
      && points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x < (Dom.xs + points[i].r)
      && points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    }

    // if a point_particle straddles all three boundaries, an eighth is needed
    if(points[i].x < (Dom.xs + points[i].r)
      && points[i].y < (Dom.ys + points[i].r)
      && points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].y < (Dom.ys + points[i].r)
      && points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x < (Dom.xs + points[i].r)
      && points[i].y > (Dom.ye - points[i].r)
      && points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x < (Dom.xs + points[i].r)
      && points[i].y < (Dom.ys + points[i].r)
      && points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].y > (Dom.ye - points[i].r)
      && points[i].z < (Dom.zs + points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z + Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].y < (Dom.ys + points[i].r)
      && points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y + Dom.yl;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x < (Dom.xs + points[i].r)
      && points[i].y > (Dom.ye - points[i].r)
      && points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x + Dom.xl;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    } if(points[i].x > (Dom.xe - points[i].r)
      && points[i].y > (Dom.ye - points[i].r)
      && points[i].z > (Dom.ze - points[i].r)) {
      ind[j] = i;
      points_virt[j].r = points[i].r;
      points_virt[j].x = points[i].x - Dom.xl;
      points_virt[j].y = points[i].y - Dom.yl;
      points_virt[j].z = points[i].z - Dom.zl;
      points_virt[j].u = points[i].u;
      points_virt[j].v = points[i].v;
      points_virt[j].w = points[i].w;
      points_virt[j].udot = points[i].udot;
      points_virt[j].vdot = points[i].vdot;
      points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
      points_virt[j].Fx = points[i].Fx;
      points_virt[j].Fy = points[i].Fy;
      points_virt[j].Fz = points[i].Fz;
      points_virt[j].Lx = points[i].Lx;
      points_virt[j].Ly = points[i].Ly;
      points_virt[j].Lz = points[i].Lz;
      j++;
    }

    ind[j] = i;
    points_virt[j].r = points[i].r;
    points_virt[j].x = points[i].x;
    points_virt[j].y = points[i].y;
    points_virt[j].z = points[i].z;
    points_virt[j].u = points[i].u;
    points_virt[j].v = points[i].v;
    points_virt[j].w = points[i].w;
    points_virt[j].udot = points[i].udot;
    points_virt[j].vdot = points[i].vdot;
    points_virt[j].wdot = points[i].wdot;
      
      
      points_virt[j].ox = points[i].ox;
      points_virt[j].oy = points[i].oy;
      points_virt[j].oz = points[i].oz;
      points_virt[j].oxdot = points[i].oxdot;
      points_virt[j].oydot = points[i].oydot;
      points_virt[j].ozdot = points[i].ozdot;
    points_virt[j].Fx = points[i].Fx;
    points_virt[j].Fy = points[i].Fy;
    points_virt[j].Fz = points[i].Fz;
    points_virt[j].Lx = points[i].Lx;
    points_virt[j].Ly = points[i].Ly;
    points_virt[j].Lz = points[i].Lz;
    j++;
  }

  fprintf(outfile, "<VTKFile type=\"PolyData\">\n");
  fprintf(outfile, "<PolyData>\n");
  fprintf(outfile, "<Piece NumberOfPoints=\"%d\" NumberOfVerts=\"0\" ",
    NNODES*npoints_plot);
  fprintf(outfile, "NumberOfLines=\"0\" NumberOfStrips=\"0\" ");
  fprintf(outfile, "NumberOfPolys=\"0\">\n");
  fprintf(outfile, "<Points>\n");
  fprintf(outfile, "<DataArray type=\"Float32\" NumberOfComponents=\"3\" ");
  fprintf(outfile, "format=\"ascii\">\n");

  real r, X, Y, Z;

  // write the locations of the nodes
  for(i = 0; i < npoints_plot; i++) {
    for(j = 0; j < NNODES; j++) {
      r = points[i].rs;
      // convert node (r, theta, phi) to (x, y, z)
      X = r * sin(node_t[j]) * cos(node_p[j]);
      Y = r * sin(node_t[j]) * sin(node_p[j]);
      Z = r * cos(node_t[j]);

      // shift from point_particle center
      X = X + points[i].x;
      Y = Y + points[i].y;
      Z = Z + points[i].z;

      if(X < dom->xs) X = X + dom->xl;
      else if(X > dom->xe) X = X - dom->xl;
      if(Y < dom->ys) Y = Y + dom->yl;
      else if(Y > dom->ye) Y = Y - dom->yl;
      if(Z < dom->zs) Z = Z + dom->zl;
      else if(Z > dom->xe) Z = Z - dom->zl;
      
      // write to file
      fprintf(outfile, "%f %f %f ", X, Y, Z);
    }
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "</Points>\n");
  fprintf(outfile, "<PointData Scalars=\"n status\">\n");
  fprintf(outfile, "<DataArray type=\"Int32\" Name=\"n\" ");
  fprintf(outfile, "format=\"ascii\">\n");

  // write index of each point_particle to each of the nodes
  for(i = 0; i < npoints_plot; i++) {
    for(j = 0; j < NNODES; j++) {
      fprintf(outfile, "%d ", ind[i]);
    }
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "<DataArray type=\"Int32\" Name=\"status\" ");
  fprintf(outfile, "format=\"ascii\">\n");

  // write status of each node
  for(i = 0; i < npoints_plot; i++) {
    for(j = 0; j < NNODES; j++) {
      fprintf(outfile, "%d ", points[i].nodes[j]);
    }
  }

  fprintf(outfile, "\n</DataArray>\n");
  fprintf(outfile, "</PointData>\n");
  fprintf(outfile, "</Piece>\n");
  fprintf(outfile, "</PolyData>\n");
  fprintf(outfile, "</VTKFile>\n");

  // clean up
  free(points_virt);
  free(ind);

  fclose(outfile);
}
*/
