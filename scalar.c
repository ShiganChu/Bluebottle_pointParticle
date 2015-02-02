#include <time.h>

#include "bluebottle.h"
#include "scalar.h"
#include "domain.h"
//#include "point.h"

int nsubdom;
real sc_init;
void scalar_read_input(void)
{

  int fret = 0;
  fret = fret; // prevent compiler warning

  cpumem = 0;
  gpumem = 0;

  // open configuration file for reading
  char fname[FILE_NAME_SIZE];
  sprintf(fname, "%s/input/scalar.config", ROOT_DIR);
  FILE *infile = fopen(fname, "r");
  if(infile == NULL) {
    fprintf(stderr, "Could not open file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  char buf[CHAR_BUF_SIZE];  // character read buffer

  // read domain
  fret = fscanf(infile, "n %d\n", &nsubdom);
  fret = fscanf(infile, "PHYSICAL PARAMETERS\n");

//default value for DIFF_eq
DIFF_eq=1.f;
#ifdef DOUBLE
  fret = fscanf(infile, "DIFF %lf\n", &DIFF);
  fret = fscanf(infile, "DIFF_eq %lf\n", &DIFF_eq);
#else // single
  fret = fscanf(infile, "DIFF %f\n", &DIFF);
  fret = fscanf(infile, "DIFF_eq %f\n", &DIFF_eq);
#endif
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "BOUNDARY CONDITIONS\n");
  fret = fscanf(infile, "SCALAR\n");

#ifdef DOUBLE
  fret = fscanf(infile, "sc_bc.scW %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    sc_bc.scW = PERIODIC;
  else if(strcmp(buf, "NEUMANN") == 0)
    sc_bc.scW = NEUMANN;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    sc_bc.scW = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &sc_bc.scWDm, &sc_bc.scWDa);
    sc_bc.scWD = 0;
  }
  else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "sc_bc.scE %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    sc_bc.scE = PERIODIC;
  else if(strcmp(buf, "NEUMANN") == 0)
    sc_bc.scE = NEUMANN;
 else if(strcmp(buf, "DIRICHLET") == 0) {
    sc_bc.scE = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &sc_bc.scEDm, &sc_bc.scEDa);
    sc_bc.scED = 0;
  }

  else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "sc_bc.scS %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    sc_bc.scS = PERIODIC;
  else if(strcmp(buf, "NEUMANN") == 0)
    sc_bc.scS = NEUMANN;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    sc_bc.scS = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &sc_bc.scSDm, &sc_bc.scSDa);
    sc_bc.scSD = 0;
  }

  else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "sc_bc.scN %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    sc_bc.scN = PERIODIC;
  else if(strcmp(buf, "NEUMANN") == 0)
    sc_bc.scN = NEUMANN;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    sc_bc.scN = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &sc_bc.scNDm, &sc_bc.scNDa);
    sc_bc.scND = 0;
  }
  else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "sc_bc.scB %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    sc_bc.scB = PERIODIC;
  else if(strcmp(buf, "NEUMANN") == 0)
    sc_bc.scB = NEUMANN;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    sc_bc.scB = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &sc_bc.scBDm, &sc_bc.scBDa);
    sc_bc.scBD = 0;
  }
  else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "sc_bc.scT %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    sc_bc.scT = PERIODIC;
  else if(strcmp(buf, "NEUMANN") == 0)
    sc_bc.scT = NEUMANN;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    sc_bc.scT = DIRICHLET;
    fret = fscanf(infile, "%lf %lf", &sc_bc.scTDm, &sc_bc.scTDa);
    sc_bc.scTD = 0;
  }
#else
  fret = fscanf(infile, "sc_bc.scW %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    sc_bc.scW = PERIODIC;
  else if(strcmp(buf, "NEUMANN") == 0)
    sc_bc.scW = NEUMANN;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    sc_bc.scW = DIRICHLET;
    fret = fscanf(infile, "%f %f", &sc_bc.scWDm, &sc_bc.scWDa);
    sc_bc.scWD = 0;
  }
  else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "sc_bc.scE %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    sc_bc.scE = PERIODIC;
  else if(strcmp(buf, "NEUMANN") == 0)
    sc_bc.scE = NEUMANN;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    sc_bc.scE = DIRICHLET;
    fret = fscanf(infile, "%f %f", &sc_bc.scEDm, &sc_bc.scEDa);
    sc_bc.scED = 0;
  }
  else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "sc_bc.scS %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    sc_bc.scS = PERIODIC;
  else if(strcmp(buf, "NEUMANN") == 0)
    sc_bc.scS = NEUMANN;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    sc_bc.scS = DIRICHLET;
    fret = fscanf(infile, "%f %f", &sc_bc.scSDm, &sc_bc.scSDa);
    sc_bc.scSD = 0;
  }

  else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "sc_bc.scN %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    sc_bc.scN = PERIODIC;
  else if(strcmp(buf, "NEUMANN") == 0)
    sc_bc.scN = NEUMANN;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    sc_bc.scN = DIRICHLET;
    fret = fscanf(infile, "%f %f", &sc_bc.scNDm, &sc_bc.scNDa);
    sc_bc.scND = 0;
  }
  else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "sc_bc.scB %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    sc_bc.scB = PERIODIC;
  else if(strcmp(buf, "NEUMANN") == 0)
    sc_bc.scB = NEUMANN;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    sc_bc.scB = DIRICHLET;
    fret = fscanf(infile, "%f %f", &sc_bc.scBDm, &sc_bc.scBDa);
    sc_bc.scBD = 0;
  }
  else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "sc_bc.scT %s", buf);
  if(strcmp(buf, "PERIODIC") == 0)
    sc_bc.scT = PERIODIC;
  else if(strcmp(buf, "NEUMANN") == 0)
    sc_bc.scT = NEUMANN;
  else if(strcmp(buf, "DIRICHLET") == 0) {
    sc_bc.scT = DIRICHLET;
    fret = fscanf(infile, "%f %f", &sc_bc.scTDm, &sc_bc.scTDa);
    sc_bc.scTD = 0;
  }

#endif

  else {
    fprintf(stderr, "turb.config read error.\n");
    exit(EXIT_FAILURE);
  }
  fret = fscanf(infile, "\n");

  fret = fscanf(infile, "INITIAL CONDITION\n");
  fret = fscanf(infile, "sc_init_cond %s", buf);
  if(strcmp(buf, "QUIESCENT") == 0)
    sc_init_cond = QUIESCENT;
  else if(strcmp(buf, "SHEAR") == 0)
    sc_init_cond = SHEAR;
  else if(strcmp(buf, "POISEUILLE") == 0)
    sc_init_cond = POISEUILLE;
  else if(strcmp(buf, "OSCILLATORY") == 0)
    sc_init_cond = OSCILLATORY;
  else if(strcmp(buf, "BULK") == 0)
{    sc_init_cond = BULK;
  #ifdef DOUBLE
   fret = fscanf(infile, "%lf", &sc_init);
  #else
   fret = fscanf(infile, "%f", &sc_init);
  #endif
}
  else if(strcmp(buf, "TURB") == 0)
    sc_init_cond = TURB;
  else if(strcmp(buf, "TAYLOR") == 0)
    sc_init_cond = TAYLOR;
  else if(strcmp(buf, "DIFFUSION") == 0)
    sc_init_cond = DIFFUSION;
  else if(strcmp(buf, "EVEN") == 0)
    sc_init_cond = EVEN;
  else {
    fprintf(stderr, "scalar.config read error.\n");
    fprintf(stderr, "INITIAL CONDITION not found.\n");
    exit(EXIT_FAILURE);
  }

 fret = fscanf(infile, "\n");
#ifdef DOUBLE
  fret = fscanf(infile, "sc_soluble_proportion %lf\n",&sc_init_percent);
#else
  fret = fscanf(infile, "sc_soluble_proportion %f\n",&sc_init_percent);
#endif

  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "EQUILIBRIUM CONCENTRATION\n");

#ifdef DOUBLE
  fret = fscanf(infile, "sc_eq %lf\n",&sc_eq);
#else
  fret = fscanf(infile, "sc_eq %f\n",&sc_eq);
#endif

  fret = fscanf(infile, "\n");
  fret = fscanf(infile, "TWO WAY FORCING\n");
#ifdef DOUBLE
  fret = fscanf(infile, "lpt_twoway %lf\n", &lpt_twoway);
  fret = fscanf(infile, "sc_twoway %lf\n", &sc_twoway);
#else
  fret = fscanf(infile, "lpt_twoway %f\n", &lpt_twoway);
  fret = fscanf(infile, "sc_twoway %f\n", &sc_twoway);
#endif

#ifdef DOUBLE
  fret = fscanf(infile, "n2 %lf\n", &n2);
  fret = fscanf(infile, "n3 %lf\n", &n3);
#else
  fret = fscanf(infile, "n2 %f\n", &n2);
  fret = fscanf(infile, "n3 %f\n", &n3);
#endif


}

void compute_scalar_BC(void)
{
  // scWD
  if(sc_bc.scWDa == 0) sc_bc.scWD = sc_bc.scWDm;
  else if(fabs(ttime_done*sc_bc.scWDa) > fabs(sc_bc.scWDm)) sc_bc.scWD = sc_bc.scWDm;
  else sc_bc.scWD = ttime_done*sc_bc.scWDa;
  // scED
  if(sc_bc.scEDa == 0) sc_bc.scED = sc_bc.scEDm;
  else if(fabs(ttime_done*sc_bc.scEDa) > fabs(sc_bc.scEDm)) sc_bc.scED = sc_bc.scEDm;
  else sc_bc.scED = ttime_done*sc_bc.scEDa;
  // scSD
  if(sc_bc.scSDa == 0) sc_bc.scSD = sc_bc.scSDm;
  else if(fabs(ttime_done*sc_bc.scSDa) > fabs(sc_bc.scSDm)) sc_bc.scSD = sc_bc.scSDm;
  else sc_bc.scSD = ttime_done*sc_bc.scSDa;
  // scND
  if(sc_bc.scNDa == 0) sc_bc.scND = sc_bc.scNDm;
  else if(fabs(ttime_done*sc_bc.scNDa) > fabs(sc_bc.scNDm)) sc_bc.scND = sc_bc.scNDm;
  else sc_bc.scND = ttime_done*sc_bc.scNDa;
  // scBD
  if(sc_bc.scBDa == 0) sc_bc.scBD = sc_bc.scBDm;
  else if(fabs(ttime_done*sc_bc.scBDa) > fabs(sc_bc.scBDm)) sc_bc.scBD = sc_bc.scBDm;
  else sc_bc.scBD = ttime_done*sc_bc.scBDa;
  // scTD
  if(sc_bc.scTDa == 0) sc_bc.scTD = sc_bc.scTDm;
  else if(fabs(ttime_done*sc_bc.scTDa) > fabs(sc_bc.scTDm)) sc_bc.scTD = sc_bc.scTDm;
  else sc_bc.scTD = ttime_done*sc_bc.scTDa;
}


void scalar_show_config(void)
{
  int i;  // iterator

  printf("GPU Domain Decomposition:\n");
  printf("  nsubdom = %d\n", nsubdom);
  for(i = 0; i < nsubdom; i++) {
  printf("Boundary Conditions: (0 = PERIODIC, 1 = DIRICHLET, 2 = NEUMANN)\n");
  printf("  sc_bc.scW = %d", sc_bc.scW);
  if(sc_bc.scW == DIRICHLET) printf(" %f", sc_bc.scWD);
  printf(", sc_bc.scE = %d", sc_bc.scE);
  if(sc_bc.scE == DIRICHLET) printf(" %f", sc_bc.scED);
  printf("\n");
  printf("  sc_bc.scS = %d", sc_bc.scS);
  if(sc_bc.scS == DIRICHLET) printf(" %f", sc_bc.scSD);
  printf(", sc_bc.scN = %d", sc_bc.scN);
  if(sc_bc.scN == DIRICHLET) printf(" %f", sc_bc.scND);
  printf("\n");
  printf("  sc_bc.scB = %d", sc_bc.scB);
  if(sc_bc.scB == DIRICHLET) printf(" %f", sc_bc.scBD);
  printf(", sc_bc.scT = %d", sc_bc.scT);
  if(sc_bc.scT == DIRICHLET) printf(" %f", sc_bc.scTD);
  printf("\n");

printf("\nsc_init_percent %f\n", sc_init_percent);
printf("\nsc_eq %f\n", sc_eq);
printf("\nlpt_twoway %f\n", lpt_twoway);
printf("\nsc_twoway %f\n", sc_twoway);

  }
}

int scalar_init(void)
{
  int i;    // iterator
 // int i, j, k;    // iterator
 // int C, W, E, S, N, B, T;

  // make sure there are enough GPU devices in the given range
  if(nsubdom > dev_end - dev_start + 1) {
    return EXIT_FAILURE;
  }

  // set up grid index structs
  // allocate and initialize scalar vectors
  sc0 = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  sc = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  diff0_sc = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  conv0_sc = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  diff_sc = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  conv_sc = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  scSrc = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);

  epsp = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);

// initialize scalar as 0 (default) 
  for(i = 0; i < Dom.Gcc.s3b; i++) {
    sc0[i] = 0.;
    sc[i] = 0.;
    diff0_sc[i]=0.;
    diff_sc[i]=0.;
    conv0_sc[i]=0.;
    conv_sc[i]=0.;
    scSrc[i]=0.;
    epsp[i]=0.;
  }


if(sc_init_cond==DIFFUSION)
{
//printf("\ninitial cond DIFFUSION\n");

 real t0=0.1;
 int dim=3;

if(dim==2)
{
 for(int k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
      for(int j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
        for(int i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
          int C = i+j*Dom.Gcc.s1b+k*Dom.Gcc.s2b;
 real x = ((i-0.5) * Dom.dx) + Dom.xs;
 real y = ((j-0.5) * Dom.dy) + Dom.ys;
 real a=1/4.0/DIFF_eq/t0;

    sc[C] = a*exp(-(x*x+y*y)*a)/PI;
    sc0[C] = sc[C];

    diff_sc[C]=DIFF_eq*(4*a*a*(x*x+y*y)-4*a)*sc[C];
    diff0_sc[C]=diff_sc[C];
    conv0_sc[C]=0.;
    conv_sc[C]=0.;
    scSrc[C]=0.;
		 }
		}	
	       }
}
else if(dim==3)
{
 for(int k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
      for(int j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
        for(int i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
          int C = i+j*Dom.Gcc.s1b+k*Dom.Gcc.s2b;
 real x = ((i-0.5) * Dom.dx) + Dom.xs;
 real y = ((j-0.5) * Dom.dy) + Dom.ys;
 real z = ((k-0.5) * Dom.dz) + Dom.zs;
 real a=1/4.0/DIFF_eq/t0;
    sc[C] = pow(a/PI,1.5)*exp(-(x*x+y*y+z*z)*a);
    sc0[C] = sc[C];
		}	
	       }
}
}

}

if(sc_init_cond==BULK)
{
 for(int k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) 
      for(int j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) 
        for(int i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
          int C = i+j*Dom.Gcc.s1b+k*Dom.Gcc.s2b;
    sc[C] = sc_init;
//printf("\nsc_init %f \n",sc_init);
    sc0[C] = sc[C];
}
}



if(sc_init_cond==EVEN)
{
//printf("\ninitial cond EVEN\n");
int dim=3;

if(dim==2){
 for(int k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
      for(int j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
        for(int i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
          int C = i+j*Dom.Gcc.s1b+k*Dom.Gcc.s2b;
 real x = ((i-0.5) * Dom.dx) + Dom.xs;
 real y = ((j-0.5) * Dom.dy) + Dom.ys;

 real k=1;

//    sc[C] = 8*cos(k*x)*cos(k*y);
    sc[C] = cos(k*x)*cos(k*y);
    sc0[C] = sc[C];

    diff_sc[C]=0;
    diff0_sc[C]=0;
    conv0_sc[C]=0.;
    conv_sc[C]=0.;
    scSrc[C]=0.;
		 }
		}	
	       }
}
else if(dim==3)
{
 for(int k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
      for(int j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
        for(int i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
          int C = i+j*Dom.Gcc.s1b+k*Dom.Gcc.s2b;
 real x = ((i-0.5) * Dom.dx) + Dom.xs;
 real y = ((j-0.5) * Dom.dy) + Dom.ys;
 real z = ((k-0.5) * Dom.dz) + Dom.zs;

 real k=1;

    sc[C] = cos(k*x)*cos(k*y)*cos(k*z);
    sc0[C] = sc[C];

    diff_sc[C]=0;
    diff0_sc[C]=0;
    conv0_sc[C]=0.;
    conv_sc[C]=0.;
    scSrc[C]=0.;
		 }
		}	
	       }

}

}



  dt_sc = 2 * nu / (Dom.dx * Dom.dx);
  dt_sc += 2 * nu / (Dom.dy * Dom.dy);
  dt_sc += 2 * nu / (Dom.dz * Dom.dz);
  dt_sc = CFL / dt_sc;

  dt_try=dt_sc;
  dt0_try=-1.;


  // initialize some variables

  rec_scalar_stepnum_out = 0;
  return EXIT_SUCCESS;
}


/*
int scalar_test(void)
{
  int i;    // iterator
 // int i, j, k;    // iterator
 // int C, W, E, S, N, B, T;

  // make sure there are enough GPU devices in the given range
  if(nsubdom > dev_end - dev_start + 1) {
    return EXIT_FAILURE;
  }
int k =(int) ((Dom.Gcc.ksb+Dom.Gcc.keb)/2);
      for(int j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
        for(int i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
          int C = i+j*Dom.Gcc.s1b+k*Dom.Gcc.s2b;
int C0 = i+j*Dom.Gcc.s1b+(k-1)*Dom.Gcc.s2b;
int C1 = i+j*Dom.Gcc.s1b+(k+1)*Dom.Gcc.s2b;
real ddsdzz=(sc[C0]+sc[C1]-2*sc[C])/Dom.dz/Dom.dz;
//if(i==33) printf("\ntest ddsdzz sc j %f %f %d\n",ddsdzz,sc[C],j);
//if(i==42&&j==57&&k==65)  printf("\nsc test %f %f %f %f\n",sc[C0],sc[C1],sc[C],ddsdzz);
if(j==57&&k==65)  printf("\nsc test,i %d %f %f %f %f\n",i,sc[C0],sc[C1],sc[C],ddsdzz);

	}	
	 }

  return EXIT_SUCCESS;
}

*/

void scalar_clean(void)
{
 free(sc);
 free(sc0);
 free(diff0_sc);
 free(diff_sc);
 free(conv0_sc);
 free(conv_sc);
 free(scSrc);
 free(epsp);
}



//TODO remember to change the corresponding part in domain.c when modify scalar
//sc_out_restart is in out_restart of domain.c
//sc_in_restart is in in_restart of domain.c
//The read and cgns write of scalar are in recorder.c!!!
