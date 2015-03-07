#include "cuda_scalar.h"
#include "cuda_point.h"

__global__ void scSrc_value_init(dom_struct *dom,real *scSrc,real a,int coordiSys)
{
int i,j,k,C;
    i = blockIdx.x*blockDim.x + threadIdx.x;
    j = blockIdx.y*blockDim.y + threadIdx.y;
int   index=coordiSys*48	+21;
int s1b=tex1Dfetch(texRefDomInfo,index);
int s2b=tex1Dfetch(texRefDomInfo,index+1);

int is,js,ks,ie,je,ke;
//get domain start and end index
int incGhost=1;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);
for(k=ks;k<ke;k++)
{
if(i<ie&&j<je)
{
C=i+j*s1b+k*s2b;
scSrc[C] =a;
}
}


}

__global__ void scSrc_value_add(dom_struct *dom,real *scSrc,real *scSrc_buf,int coordiSys)
{
int i,j,k,C;

    i = blockIdx.x*blockDim.x + threadIdx.x;
    j = blockIdx.y*blockDim.y + threadIdx.y;
int   index=coordiSys*48	+21;
int s1b=tex1Dfetch(texRefDomInfo,index);
int s2b=tex1Dfetch(texRefDomInfo,index+1);

int is,js,ks,ie,je,ke;
//get domain start and end index
int incGhost=1;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);
for(k=ks;k<ke;k++)
{
if(i<ie&&j<je)
{
C=i+j*s1b+k*s2b;
scSrc[C] +=scSrc_buf[C];
}
}
}


__global__ void boundary_face_value_periodic_start(dom_struct *dom,real *scSrc,int coordiSys)
{
int i,j,k;
/*
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
*/
int   index=coordiSys*48	+21;
int s1b=tex1Dfetch(texRefDomInfo,index);
int s2b=tex1Dfetch(texRefDomInfo,index+1);

int is,js,ks,ie,je,ke;
//get domain start and end index
int incGhost=0;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);
int C1,C2;
switch(coordiSys)
{
case 1:
     k = blockIdx.x*blockDim.x + threadIdx.x;
     j = blockIdx.y*blockDim.y + threadIdx.y;
if(k>=ks&&k<ke &&j>=js&&j<je)
{
C1=is+j*s1b+k*s2b;
C2=ie-1+j*s1b+k*s2b;
scSrc[C1] +=scSrc[C2];
}
break;
case 2:
     k = blockIdx.x*blockDim.x + threadIdx.x;
     i = blockIdx.y*blockDim.y + threadIdx.y;
if(i>=is&&i<ie &&k>=ks&&k<ke)
{
C1=i+js*s1b+k*s2b;
C2=i+(je-1)*s1b+k*s2b;
scSrc[C1] +=scSrc[C2];
}
break;
case 3:
     i = blockIdx.x*blockDim.x + threadIdx.x;
     j = blockIdx.y*blockDim.y + threadIdx.y;
if(i>=is&&i<ie &&j>=js&&j<je)
{
C1=i+j*s1b+ks*s2b;
C2=i+j*s1b+(ke-1)*s2b;

//if(fabs(scSrc[C1])>EPSILON) printf("\nscSrc %f %f %d %d %d\n",scSrc[C1],scSrc[C2],i,j,ks);

scSrc[C1] +=scSrc[C2];
}
break;
default: break;
}
}


__global__ void boundary_face_value_periodic_end(dom_struct *dom,real *scSrc,int coordiSys)
{

int i,j,k;
/*
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
*/
int   index=coordiSys*48	+21;
int s1b=tex1Dfetch(texRefDomInfo,index);
int s2b=tex1Dfetch(texRefDomInfo,index+1);

int is,js,ks,ie,je,ke;
//get domain start and end index
int incGhost=0;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);
int C1,C2;
//int C3,C4,C5,C6;
switch(coordiSys)
{
case 1:
     k = blockIdx.x*blockDim.x + threadIdx.x;
     j = blockIdx.y*blockDim.y + threadIdx.y;
if(k>=ks&&k<ke &&j>=js&&j<je)
{
C1=is+j*s1b+k*s2b;
C2=ie-1+j*s1b+k*s2b;
scSrc[C2] =scSrc[C1];

/*
//Periodic BC
C3=is-1+j*s1b+k*s2b;
C4=ie+j*s1b+k*s2b;
C5=is+1+j*s1b+k*s2b;
C6=ie-2+j*s1b+k*s2b;
scSrc[C3]=scSrc[C6];
scSrc[C4]=scSrc[C5];
*/
}
break;
case 2:
     k = blockIdx.x*blockDim.x + threadIdx.x;
     i = blockIdx.y*blockDim.y + threadIdx.y;
if(i>=is&&i<ie &&k>=ks&&k<ke)
{
C1=i+js*s1b+k*s2b;
C2=i+(je-1)*s1b+k*s2b;
scSrc[C2] =scSrc[C1];

/*
//Periodic BC
C3=i+(js-1)*s1b+k*s2b;
C4=i+je*s1b+k*s2b;
C5=i+(js+1)*s1b+k*s2b;
C6=i+(je-2)*s1b+k*s2b;
scSrc[C3]=scSrc[C6];
scSrc[C4]=scSrc[C5];
*/
}
break;
case 3:
     i = blockIdx.x*blockDim.x + threadIdx.x;
     j = blockIdx.y*blockDim.y + threadIdx.y;
if(i>=is&&i<ie &&j>=js&&j<je)
{
C1=i+j*s1b+ks*s2b;
C2=i+j*s1b+(ke-1)*s2b;
scSrc[C2] =scSrc[C1];

/*
//Periodic BC
C3=i+j*s1b+(ks-1)*s2b;
C4=i+j*s1b+ke*s2b;
C5=i+j*s1b+(ks+1)*s2b;
C6=i+j*s1b+(ke-2)*s2b;

scSrc[C3]=scSrc[C6];
scSrc[C4]=scSrc[C5];
*/
}
break;
default: break;
}
}

__global__ void boundary_face_value_homo_end(dom_struct *dom,real *scSrc,int coordiSys)
{

int i,j,k;
int   index=coordiSys*48	+21;
int s1b=tex1Dfetch(texRefDomInfo,index);
int s2b=tex1Dfetch(texRefDomInfo,index+1);

int is,js,ks,ie,je,ke;
//get domain start and end index
int incGhost=0;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);
int C1,C2,C3;
if(coordiSys==3)
{
     i = blockIdx.x*blockDim.x + threadIdx.x;
     j = blockIdx.y*blockDim.y + threadIdx.y;
if(i>=is-1&&i<ie+1 &&j>=js-1&&j<je+1)
{
C1=i+j*s1b+ke*s2b;
C2=i+j*s1b+(ke-1)*s2b;
C3=i+j*s1b+(ks-1)*s2b;
scSrc[C1] =0;
scSrc[C2] =0;
scSrc[C3] =0;
}

     i = blockIdx.x*blockDim.x + threadIdx.x;
     k = blockIdx.y*blockDim.y + threadIdx.y;
if(i>=is-1&&i<ie+1 &&k>=ks-1&&k<ke+1)
{
C1=i+(js-1)*s1b+k*s2b;
C2=i+je*s1b+k*s2b;
scSrc[C1] =0;
scSrc[C2] =0;
}

     j = blockIdx.x*blockDim.x + threadIdx.x;
     k = blockIdx.y*blockDim.y + threadIdx.y;
if(j>=js-1&&j<je+1 &&k>=ks-1&&k<ke+1)
{
C1=is-1+j*s1b+k*s2b;
C2=ie+j*s1b+k*s2b;
scSrc[C1] =0;
scSrc[C2] =0;
}
}
}
__global__ void array_init(real *A,dom_struct *dom, int n, real a)
{
int indx =  threadIdx.x + blockIdx.x*blockDim.x;
int indy =  threadIdx.y + blockIdx.y*blockDim.y;
int pp=indx+indy*gridDim.x*blockDim.x;
if(pp<n) A[pp]=a;
}




__global__ void point_vel_specify(real *ug,real *vg,real *wg,point_struct *points,int npoints)
{
int pp = blockIdx.x*blockDim.x + threadIdx.x;
if(pp>=npoints) return;

points[pp].u=ug[pp];
points[pp].v=vg[pp];
points[pp].w=wg[pp];
}

__global__ void point_ms_initD(point_struct *points, int npoints,int percent)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;
  if(pp>=npoints) return;

real r=points[pp].r;
real rhod=points[pp].rho;
real volume=4*PI/3.f *r*r*r;

points[pp].ms=rhod*volume*percent;
points[pp].msdot=0;
}

//calculate particle diameter
__global__ void points_dp(point_struct *points, real *dp, int npoints)

{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;
if(pp<npoints) 
{
real rp=points[pp].r;
dp[pp]=2*rp;
}
}

//calculate particle diameter
__global__ void points_Fz(point_struct *points, real *dp, int npoints)

{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;
if(pp<npoints) 
{
dp[pp]=points[pp].Fz;
}
}
__global__
void copy_points_dt(real *pdt,point_struct *points,int npoints)
{
int pp = blockIdx.x*blockDim.x + threadIdx.x;
if(pp>=npoints) return;
pdt[pp]=points[pp].dt;
}

/* The GPU kernel that performs the power-of-two maximum value search
 * algorithm.
 */
__global__ void entrySearch_max_int_kernel(int *g_iarr, int *g_max_intarr,
  int size)
{
    // create shared memory
    extern __shared__ int sarr[];

    // load shared mem
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x * blockDim.x * 2 + threadIdx.x;
 
    if(i + blockDim.x < size) {
      if(g_iarr[i] > g_iarr[i + blockDim.x]) {
        sarr[tid] = g_iarr[i];
      } else {
        sarr[tid] = g_iarr[i + blockDim.x];
      }
    } else if (i < size) {
      sarr[tid] = g_iarr[i];
    } else {
#ifdef DOUBLE
        sarr[tid] = DBL_MIN;
#else
        sarr[tid] = FLT_MIN;
#endif
    }

    __syncthreads();

    // do comparison in shared mem
    for(unsigned int s=blockDim.x/2; s>0; s>>=1) {
      if(tid < s) {
        if(sarr[tid] < sarr[tid + s]) {
          sarr[tid] = sarr[tid + s];
        }
      }
      __syncthreads();
    }
  
    // write result for this block to global mem
    if(tid == 0) {
      g_max_intarr[blockIdx.x] = sarr[0];
    }
}




__global__
void gaussian_array_initD(float * GaussianKernel,real dg2_sig2,real dg2,real norm)
{
int k = blockIdx.x*blockDim.x + threadIdx.x;
if(k>=LEN_GAUSSIAN_ARRAY) return; 
/*
dg2_sig2=dg*dg/2.f/sig/sig
real min_meshsize=min(min(dx,dy),dz);
//2.0f*sqrt(2.0f*log(2.0f))=2.3548f;
real sig= min_meshsize/2.3548f;
real norm=pow(sqrt(2*PI)*sig,3);
real dg=3.f*sig/(float)LEN_GAUSSIAN_ARRAY;
*/
if(k<LEN_GAUSSIAN_ARRAY-1) GaussianKernel[k]=exp(-k*dg2_sig2)/norm;
if(k==LEN_GAUSSIAN_ARRAY-1) GaussianKernel[k]=1/dg2;

//if(GaussianKernel[k]<EPSILON) printf("\nk %d %f\n",k,GaussianKernel[k]);
}


//This method will work with higher efficiency for thousands of particles
__global__
void lpt_mollify_delta_scD(dom_struct *dom,
	      real *flowSource,
	      real *lptSourceVal,
              int   *cellStart,
              int   *cellEnd,
              int    npoints,
	      int    coordiSys)
{

int is,js,ks,ie,je,ke;

__shared__ int cell_start[MAX_THREADS_DIM * MAX_THREADS_DIM];
__shared__ int cell_end[MAX_THREADS_DIM * MAX_THREADS_DIM];

 
    // subdomain indices
    int i = blockIdx.x*blockDim.x + threadIdx.x ;
    int j = blockIdx.y*blockDim.y + threadIdx.y ;

    // shared memory indices
    int ti = threadIdx.x;
    int tj = threadIdx.y;

//get domain start and end index
int incGhost=0;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);
//if(i==2&&j==3) printf("\ncoordi %d %d %d %d %d\n",is,ie,je,ke,coordiSys);

   __syncthreads();


for(int k=ks;k<ke;k++)
{
int gridHash=calcGridHash(i,j,k,dom,coordiSys);
if(i<ie&&i>=is && j<je&&j>=js)
{
//TODO reduce Repeat calculation of gridHash
cell_start[ti+tj*blockDim.x]=cellStart[gridHash];
cell_end[ti+tj*blockDim.x]=cellEnd[gridHash];
}

   __syncthreads();

//Use Shared memory just read the particles around it 
//Determin int the future, if we need 3D blocks and threads


//	real force=0.f;
	real forceCell = 0.f;
	int start_index,end_index;

    if((j >= js && j < je)
     && (i >= is && i < ie))
		{
		// store grid hash and partisle pp
		start_index=cell_start[ti+tj*blockDim.x];
	
if(start_index>=0)
			{ 
	  	end_index=cell_end[ti+tj*blockDim.x];
				for (int index=start_index; index<end_index; index++)
{				forceCell += lptSourceVal[index];
//if(coordiSys==3||coordiSys==0)
//if(coordiSys==0)  printf("\nlptSourceVal %d %d %d %d %d %f\n",index,coordiSys,i,j,k,lptSourceVal[index]);
}
			}
		}

    __syncthreads();

    if((j >= js && j < je)
     && (i >= is && i < ie))
		{
		 // examine neighbouring cells
		// flowSource[gridHash] +=forceCell;
		 flowSource[gridHash] =forceCell;
/*
if(fabs(forceCell)>EPSILON) 
{
printf("\nforceCell %f %d %d %d\n",forceCell,i,j,k);
//printf("\nke %d %d\n",ke,coordiSys);
}
*/
//if(coordiSys==0) if(fabs(forceCell)>EPSILON)  printf("\nflowSource %d %d %d %d %d %f\n",gridHash,coordiSys,i,j,k,forceCell);
		}

	}
}



//Use periodic BC for the particle source; Use cache to read layer by layer
__global__
void lpt_mollify_sc_ksi_optD(dom_struct *dom,
              real *A,
	      real *posX,
	      real *posY,
	      real *posZ,
	      real *Ksi,
              int   *cellStart,
              int   *cellEnd,
              int   *gridFlowHash,
              int   *gridParticleIndex,
              int    npoints,int maxPointsPerCell,
              int coordiSys,int valType)
{

//get domain start and end index
int is,js,ks,ie,je,ke;
int incGhost=0;//using ghost will cause problem, that the array maximum will be violated
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);

//int isb,jsb,ksb,ieb,jeb,keb;
int isb,jsb,ieb,jeb;
isb=is-DOM_BUF;
jsb=js-DOM_BUF;
//ksb=ks-DOM_BUF;
ieb=ie+DOM_BUF;
jeb=je+DOM_BUF;
//keb=ke+DOM_BUF;

int   index=coordiSys*48	+21;
int s1b=tex1Dfetch(texRefDomInfo,index);
int s2b=tex1Dfetch(texRefDomInfo,index+1);

__shared__ int cell_start0[MAX_THREADS_DIM * MAX_THREADS_DIM];
__shared__ int cell_start1[MAX_THREADS_DIM * MAX_THREADS_DIM];
__shared__ int cell_start2[MAX_THREADS_DIM * MAX_THREADS_DIM];

__shared__ int cell_end0[MAX_THREADS_DIM * MAX_THREADS_DIM];
__shared__ int cell_end1[MAX_THREADS_DIM * MAX_THREADS_DIM];
__shared__ int cell_end2[MAX_THREADS_DIM * MAX_THREADS_DIM];

//int maxThreadsDim2=MAX_THREADS_DIM * MAX_THREADS_DIM;


int gridHash0,gridHash1,gridHash2;
    // the extra 2*blockIdx.X terms implement the necessary overlapping of shared memory blocks in the subdomain
    int i = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int ti = threadIdx.x;
    int tj = threadIdx.y;


//Only combine the restriction of ti&tj and i&j could we get unique i&j
//Should contain boundary here
if(i<ieb&&i>=isb && j<jeb&&j>=jsb)
{
//using cache to store gridHash is more efficient than using shared memory
 gridHash0=gridFlowHash[i+j*s1b+s2b*(ks-1)];
 gridHash1=gridFlowHash[i+j*s1b+s2b*ks];
 gridHash2=gridFlowHash[i+j*s1b+s2b*(ks+1)];

//representing three layers
cell_start0[ti+tj*blockDim.x]=cellStart[gridHash0];
cell_start1[ti+tj*blockDim.x]=cellStart[gridHash1];
cell_start2[ti+tj*blockDim.x]=cellStart[gridHash2];

cell_end0[ti+tj*blockDim.x]=cellEnd[gridHash0];
cell_end1[ti+tj*blockDim.x]=cellEnd[gridHash1];
cell_end2[ti+tj*blockDim.x]=cellEnd[gridHash2];

}


for(int k=ks;k<ke;k++)
{
//Should include ghost value
if(k>ks&&i<ieb&&i>=isb && j<jeb&&j>=jsb)
{
gridHash2=gridFlowHash[i+j*s1b+s2b*(k+1)];

//Read one layer each time
cell_start0[ti+tj*blockDim.x]=cell_start1[ti+tj*blockDim.x];
cell_start1[ti+tj*blockDim.x]=cell_start2[ti+tj*blockDim.x];
cell_start2[ti+tj*blockDim.x]=cellStart[gridHash2];

cell_end0[ti+tj*blockDim.x]=cell_end1[ti+tj*blockDim.x];
cell_end1[ti+tj*blockDim.x]=cell_end2[ti+tj*blockDim.x];
cell_end2[ti+tj*blockDim.x]=cellEnd[gridHash2];

}


   __syncthreads();
	int it,jt,kt;
	int tip,tjp;
	int ksi_ind;
	real forceCell = 0.f;
	int start_index[3];
	int end_index[3];
//should not contain boundary here
if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)&&(i<ie&&i>=is && j<je&&j>=js))
{
      //  int ip,jp,kp;//index of particle cell

//Make the sum_ksi_cell a device kernel will save half time of the calculation
for(int di=-1;di<2;di++)
for(int dj=-1;dj<2;dj++)
{
		// store grid hash and partisle pp
		tip=ti-di;tjp=tj-dj;
		start_index[0]=cell_start0[tip+tjp*blockDim.x];
		start_index[1]=cell_start1[tip+tjp*blockDim.x];
		start_index[2]=cell_start2[tip+tjp*blockDim.x];

                end_index[0]=cell_end0[tip+tjp*blockDim.x];
                end_index[1]=cell_end1[tip+tjp*blockDim.x];
                end_index[2]=cell_end2[tip+tjp*blockDim.x];


		for(int dk=-1;dk<2;dk++)
		{
			it=di+1;jt=dj+1;kt=dk+1;
		//	ip=i-di;jp=j-dj;kp=k-dk;
// if(start_index[kt]>=0)
//			{ 

ksi_ind=it+jt*STENCIL+kt*STENCIL2;

		//		for (int index=start_index[kt]; index<end_index[kt]; index++)
int inCellTag=(start_index[kt]>=0);
				for (int count=0; count<maxPointsPerCell; count++)
        			{

index=start_index[kt]+count;
inCellTag=(index<end_index[kt])&&inCellTag;
				//Read of Global memory, how to coalesce it???
				//int pp=gridParticleIndex[index]; 
				//Try to avoid bank conflicts here, if we have shared memory
				//forceCell += Ksi[it+jt*STENCIL+kt*STENCIL2+pp*STENCIL3];//calc index for Ksi takes no timea.
				//Reading the position takes little time, about 160 ct(clock_t)
			//	xp=posX[index];
			//	yp=posY[index];
			//	zp=posZ[index];
				
			//tag is the criterial for whether add or not
			//	ksi=inCellTag*Ksi[it+jt*STENCIL+kt*STENCIL2+index*inCellTag*STENCIL3];
			//	ksi=inCellTag*Ksi[ksi_ind+index*inCellTag*STENCIL3];
				forceCell +=inCellTag*Ksi[ksi_ind+index*inCellTag*STENCIL3];
			//	forceCell += ksi;
//				forceCell +=lpt_integrate_mol_opt(i,j,k,xp,yp,zp,dom,coordiSys)*Weight[index];
		//printf("\nstartEndIndex %d %d %d %f %f\n",start_index[kt],end_index[kt],index,test,Weight[index]);
		//printf("\nposition %f %f %f %d %d\n",xp,yp,zp,index,start_index[kt]);

	 
				}


//			}
		}
}

//The above are cheap, takes 1 ms or so.  The following has significant thread divergence
		// Any operation inside will be expensive
//		force+= forceCell;
//		force+= point_cell_ksi(gridHash,is,js,ks,Ksi,cellStart,cellEnd,gridParticleIndex);
}
	
//    int gridHash=calcGridHash(i,j,k,dom,coordiSys);
    int gridHash=gridFlowHash[i+j*s1b+k*s2b];

    __syncthreads();

    if((j >= js && j < je)
     && (i >= is && i < ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1)))
		{
		 // examine neighbouring cells
		 //A[gridHash] +=forceCell;
		 A[gridHash] =forceCell;
//if(forceCell>EPSILON) printf("\nforceCell %f %f %f %f %d %d %d\n",xp,yp,zp,forceCell,i,j,k,gridHash);
//	 A[gridHash]=sum_ksi_cell(i,j,k,points,dom,Ksi,cellStart,cellEnd,gridParticleIndex,coordiSys,valType);
		}

	}
}



// calculate grid hash value for each particle
__global__ void calcHash_optD(int   *gridParticleHash,  // output
               int   *gridParticleIndex, // output
               dom_struct *dom,               // input: dom info
	       real *posX,real *posY,real *posZ,
               int    npoints,
	       int coordiSys)
{
int pp =  threadIdx.x + blockIdx.x*blockDim.x;

if(pp>=npoints) return;
real xp=posX[pp];
real yp=posY[pp];
real zp=posZ[pp];

__syncthreads();

int ip,jp,kp;
//The cell-center or face-center index to locate the particle, i.e get points[pp].i,points[pp].j,points[pp].k 
calcGridPos_opt(ip,jp,kp,xp,yp,zp,dom,coordiSys);

    // calculate grid index where the particle lives in
int hash=calcGridHash(ip,jp,kp,dom,coordiSys);

//__syncthreads();
    // store grid hash and particle pp
    gridParticleHash[pp] = hash;
    gridParticleIndex[pp] = pp;
//printf("\ncalcHash %d %d %d %d %d\n",ip,jp,kp,hash,pp);
}



__global__ void lpt_point_position(point_struct *points,real *posX,real *posY,real *posZ, int npoints)
{
int pp =  threadIdx.x + blockIdx.x*blockDim.x;
if(pp>=npoints) return;
posX[pp]=points[pp].x;
posY[pp]=points[pp].y;
posZ[pp]=points[pp].z;
//printf("\npoint_position %f %f %d\n",posX[pp],points[pp].x,pp);
}


__global__ void lpt_delta_point_position(point_struct *points,
dom_struct *dom,
real *posX,real *posY,real *posZ,
real *lptSourceVal, 
int npoints, int coordiSys,int valType)
{
int pp =  threadIdx.x + blockIdx.x*blockDim.x;
if(pp>=npoints) return;

posX[pp]=points[pp].x;
posY[pp]=points[pp].y;
posZ[pp]=points[pp].z;

lptSourceVal[pp]=lpt_mol_typeVal(points,dom,pp,coordiSys,valType);
//if(isnan(lptSourceVal[pp])) printf("\npoint_position %f %f %d %d %d\n",posX[pp],points[pp].x,pp,coordiSys,valType);
}



/*
A is fluid property located at cell center or face center depending on coordiSys 
coordiSys=0 coresponds to cell-center; coordiSys=1 x-face-center;coordiSys=2 y-face-center;coordiSys=3 z-face-center
dir
*/
__global__ void lpt_point_ksi_opt(point_struct *points,
				 dom_struct *dom,
				 real *posX,real *posY,real *posZ, 
				 real *Ksi,
				 int   *gridParticleIndex,
				 int npoints,int coordiSys, int valType)
{

int index =  threadIdx.x + blockIdx.x*blockDim.x;
if(index>=npoints) return;
int pp=gridParticleIndex[index];
////printf("\nindex pp %d %d\n",index,pp);


//Define the STENCIL of the Gausian filter, the slpt_mopread length of Gaussian kernel
real ksi[3][3][3];

real  xp =  posX[index];
real  yp =  posY[index];
real  zp =  posZ[index];

__syncthreads();

// i,j,k should be the start number of cell centelpt_mor, so that xm(i)<=xp<xm(i+1)
//xm(i)= (i-0.5)*dom->dx +dom->xs
//i=floor((xm - dom->xs) * ddx+ 0.5);
//i=floor((xm - dom->xs) * ddx- 0.5) + DOM_BUF;
  int ip,jp,kp;//index of particle cell
  int ic,jc,kc;//index of object cell
  int is,js,ks;

//Calculate ip,jp,kp independently
calcGridPos_opt(ip,jp,kp,xp,yp,zp,dom,coordiSys);

real buf=0;
int began=-floor((STENCIL-1)/2.0f);//equal to -1,-1 if STENCIL=3,4 
int end=1+ceil((STENCIL-1)/2.0f);//equal to 2,3 if STENCIL=3,4

//printf("\ni~j %d %d %d %f %f %f\n",ip,jp,kp,xp,yp,zp);

//Calculate the filter strength from Gaussian kernel
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
{
is=di-began;js=dj-began;ks=dk-began;
ic=ip+di;jc=jp+dj;kc=kp+dk;
//ip~kp could be beyond domain boundary
ksi[is][js][ks]=lpt_integrate_mol_opt(ic,jc,kc,xp,yp,zp,dom,coordiSys);
//ksi[is][js][ks]=lpt_integrate_mol(ic,jc,kc,xp,yp,zp,coordiSys,xs,ys,zs,dx,dy,dz);
buf+=ksi[is][js][ks];
}

//TODO add mask infomation as in lpt_mollify_sc in lpt_interpolator.f90

//This calculation is cheap, takes about 1 ms 
real Ap=lpt_mol_typeVal(points,dom,pp,coordiSys,valType);

// Normalize  ksi = ksi/buf
if(buf>0)
		{
for(int dk=began;dk<end;dk++)
for(int dj=began;dj<end;dj++)
for(int di=began;di<end;di++)
			{
		is=di-began;js=dj-began;ks=dk-began;
		ksi[is][js][ks]=ksi[is][js][ks]/buf;
//printf("\nksi %d %d %d %f %f %f\n",di,dj,dk,ksi[is][js][ks],buf,Ap);
                //Add gausssian weight between cell center and partisle position to object
		Ksi[is+js*STENCIL+ks*STENCIL2+index*STENCIL3]=ksi[is][js][ks]*Ap;
			}
		}
}



/*
A is fluid property located at cell center or face center depending on coordiSys 
coordiSys=0 coresponds to cell-center; coordiSys=1 x-face-center;coordiSys=2 y-face-center;coordiSys=3 z-face-center
dir
*/

//calculate the maximum number of particles in a single grid cell
__global__ void calcMaxPointsPerCell_optD(dom_struct *dom, int *cellStart, int *cellEnd,int *pointsNum, int coordiSys)
{
  int i = threadIdx.x + blockIdx.x*blockDim.x;
  int j = threadIdx.x + blockIdx.x*blockDim.x;

int is,js,ks,ie,je,ke;
int incGhost=1;//using ghost will cause problem, that the array maximum will be violated
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);

int   index=coordiSys*48        +21;
int s1b=tex1Dfetch(texRefDomInfo,index);
int s2b=tex1Dfetch(texRefDomInfo,index+1);

int C;
for(int k=ks;k<ke;k++)
{
C=i+j*s1b+k*s2b;
pointsNum[C]=cellEnd[C]-cellStart[C];
}

}


// rearrange particle data into sorted order, and find the start of each cell in the sorted hash array
// cellStart[hash+1]=cellEnd[hash]; 
__global__ void findCellStart_deltaD(int   *cellStart,        // output: cell start pp
                                  int   *cellEnd,          // output: cell end pp
                                  int   *gridParticleHash, // input: sorted grid hashes
                                  int   *gridParticleIndex,   // input: sorted particle indices
				  real *lptSourceVal,
				  real *lptSourceValOld,
                                  int    npoints)
{
__shared__ int sharedHash[MAX_THREADS_1D+1]; // blockSize + 1 elements

int pp =  threadIdx.x + blockIdx.x*blockDim.x;
int ti =  threadIdx.x;

// printf("\n0wong ti pp %d %d \n",ti,pp);

    // handle case when no. of particles not multiple of block size
	if(pp>=npoints) return;

       int hash = gridParticleHash[pp];
        // Load hash data into shared memory so that we can look
        // at neighboring particle's hash value without loading
        // two hash values per thread
        sharedHash[ti+1] = hash;
        if (pp > 0 && ti == 0)
        {
            // first thread in block must load neighbor particle hash
            sharedHash[0] = gridParticleHash[pp-1];
        }

    __syncthreads();

        // If this particle has a different cell pp to the previous
        // particle then it must be the first particle in the cell,
        // so store the pp of this particle in the cell.
        // As it isn't the first particle, it must also be the cell end of
        // the previous particle's cell


        if (pp == 0 || hash != sharedHash[ti])
        {
            cellStart[hash] = pp;
	// printf("\nstart %d %d\n",hash,pp);
            if (pp > 0)
              {
        cellEnd[sharedHash[ti]] = pp;
	// printf("\nend %d %d\n",sharedHash[ti],pp);
		}
        }

        if (pp == npoints - 1)
        {
            cellEnd[hash] = pp + 1;
	 //  printf("\nend %d %d\n",hash,pp+1);
        }
   
         __syncthreads();

	int sortedIndex = gridParticleIndex[pp];

//TODO These three commands takes time,try to use shared memory to optimize it
/*
        posX[pp] = posXold[sortedIndex];
        posY[pp] = posYold[sortedIndex];
        posZ[pp] = posZold[sortedIndex];
*/
        lptSourceVal[pp] = lptSourceValOld[sortedIndex];
	//printf("\nsorterdIndex %d %d %f\n",pp,sortedIndex,posX[pp]);
}
 



// rearrange particle data into sorted order, and find the start of each cell in the sorted hash array
// cellStart[hash+1]=cellEnd[hash]; 
__global__ void findCellStart_optD(int   *cellStart,        // output: cell start pp
                                  int   *cellEnd,          // output: cell end pp
                                  int   *gridParticleHash, // input: sorted grid hashes
                                  int   *gridParticleIndex,   // input: sorted particle indices
				  real *posX,real *posY,real *posZ,
				  real *posXold,real *posYold,real *posZold,
                                  int    npoints)
{
__shared__ int sharedHash[MAX_THREADS_1D+1]; // blockSize + 1 elements

int pp =  threadIdx.x + blockIdx.x*blockDim.x;
int ti =  threadIdx.x;


    // handle case when no. of particles not multiple of block size
	if(pp>=npoints) return;

       int hash = gridParticleHash[pp];
        // Load hash data into shared memory so that we can look
        // at neighboring particle's hash value without loading
        // two hash values per thread
        sharedHash[ti+1] = hash;
        if (pp > 0 && ti == 0)
        {
            // first thread in block must load neighbor particle hash
            sharedHash[0] = gridParticleHash[pp-1];
        }

    __syncthreads();

        // If this particle has a different cell pp to the previous
        // particle then it must be the first particle in the cell,
        // so store the pp of this particle in the cell.
        // As it isn't the first particle, it must also be the cell end of
        // the previous particle's cell


        if (pp == 0 || hash != sharedHash[ti])
        {
            cellStart[hash] = pp;
	// printf("\nstart %d %d\n",hash,pp);
            if (pp > 0)
              {
        cellEnd[sharedHash[ti]] = pp;
	// printf("\nend %d %d\n",sharedHash[ti],pp);
		}
        }

        if (pp == npoints - 1)
        {
            cellEnd[hash] = pp + 1;
	//   printf("\nend %d %d\n",hash,pp+1);
        }
   
         __syncthreads();

	int sortedIndex = gridParticleIndex[pp];

//TODO These three commands takes time,try to use shared memory to optimize it
        posX[pp] = posXold[sortedIndex];
        posY[pp] = posYold[sortedIndex];
        posZ[pp] = posZold[sortedIndex];
	//printf("\nsorterdIndex %d %d %f\n",pp,sortedIndex,posX[pp]);
}
 

//calculate which grid the particle is in, won't make the grid index periodic
//Passing the dom and points pointer won't take time!!
__device__ void calcGridPos_opt(int &ip,int &jp,int &kp,real xp,real yp,real zp,dom_struct *dom,int coordiSys)
{
real xs=dom->xs;
real ys=dom->ys;
real zs=dom->zs;
real ddx=1.f/dom->dx;
real ddy=1.f/dom->dy;
real ddz=1.f/dom->dz;

int index=coordiSys*48 +coordiSys*6-4;
int Ng=tex1Dfetch(texRefDomInfo,index);
switch(coordiSys)
{
case 0:
  ip = floor((xp - xs) * ddx) + DOM_BUF;//for xm[i]
  jp = floor((yp - ys) * ddy) + DOM_BUF;//for ym[j]
  kp = floor((zp - zs) * ddz) + DOM_BUF;//for zm[k]
  break;
case 1:
  ip = floor((xp - xs) * ddx + 0.5f) + DOM_BUF;//for xm[i+1], ip should be i+1
  jp = floor((yp - ys) * ddy) + DOM_BUF;	   //for y[i]
  kp = floor((zp - zs) * ddz) + DOM_BUF;	   //for z[i]

//if(ip>=Ng) ip -=Ng-DOM_BUF;
  break;
case 2:
  ip = floor((xp - xs) * ddx) + DOM_BUF;//for x[i] 
  jp = floor((yp - ys) * ddy + 0.5f) + DOM_BUF;	   //for ym[i+1]
  kp = floor((zp - zs) * ddz) + DOM_BUF;	   //for z[i]

//if(jp>=Ng) jp -=Ng-DOM_BUF;
  break;
case 3:
  ip = floor((xp - xs) * ddx) + DOM_BUF;//for x[i] 
  jp = floor((yp - ys) * ddy) + DOM_BUF;	   //for y[i]
  kp = floor((zp - zs) * ddz+ 0.5f) + DOM_BUF;  //for zm[i+1]
//if(kp>=Ng) kp -=Ng-DOM_BUF;
  break;
default:break;
}

//printf("\nip xp xs dx  %d %f %f %f\n",ip,xp,xs,ddx);
}

//Use texture memory to fetch it!!!
//Gausian kernel to calculate weight coefficient of the filter
//ic~kc is the grid cell index, xp~kp is the particle position
//this subroutine calculates the particle contribution to the cell center of (ic,jc,kc)
__device__ real lpt_integrate_mol_opt(int ic,int jc,int kc,real xp,real yp,real zp, dom_struct *dom, int coordiSys)
{

//clock_t time11 = clock();	 
real xs=dom->xs;real ys=dom->ys;real zs=dom->zs;
real dx=dom->dx;real dy=dom->dy;real dz=dom->dz;
real xr,yr,zr;
real x,y,z,r2;
//real x[4];
//real r2;

//clock_t time12 = clock();	 

switch(coordiSys)
{
case 0:
x = (ic-DOM_BUF+0.5f) * dx + xs;
y = (jc-DOM_BUF+0.5f) * dy + ys;
z = (kc-DOM_BUF+0.5f) * dz + zs;
break;
case 1:
x =  (ic-DOM_BUF) * dx + xs;
y = (jc-DOM_BUF+0.5f) * dy + ys;
z = (kc-DOM_BUF+0.5f) * dz + zs;
break;
case 2:
x = (ic-DOM_BUF+0.5f) * dx + xs;
y = (jc-DOM_BUF) * dy + ys;
z = (kc-DOM_BUF+0.5f) * dz + zs;
break;
case 3:
x = (ic-DOM_BUF+0.5f) * dx + xs;
y = (jc-DOM_BUF+0.5f) * dy + ys;
z = (kc-DOM_BUF) * dz + zs;
break;
default:break;
}

xr=xp-x;yr=yp-y;zr=zp-z;
r2=xr*xr+yr*yr+zr*zr;

/*
//TODO make this as defined value avaible from host and device!!!
real min_meshsize=min(min(dx,dy),dz);
real max_meshsize=max(max(dx,dy),dz);
//2.0f*sqrt(2.0f*log(2.0f))=2.3548;
real sig= min_meshsize/2.3548f;
real val = exp(-r2/(2.0f*sig*sig));
//real dg= 2.6f*max_meshsize/LEN_GAUSSIAN_ARRAY;

//Ref gaussian_array_initD, the last value of the GaussianKernel is the discretization length dg.
real dg=(real) tex1Dfetch(texRefGaussian,LEN_GAUSSIAN_ARRAY-1);
real val2=(real) tex1Dfetch(texRefGaussian,sqrt(r2)/dg);

clock_t time13 = clock();	 
//1 fetch takes 200 clock cycles to finish
real idg2=(real) tex1Dfetch(texRefGaussian,LEN_GAUSSIAN_ARRAY-1);
clock_t time14 = clock();	 
real val3=(real) tex1Dfetch(texRefGaussian,r2*idg2);

clock_t time15 = clock();	 


real val1=val/pow(sqrt(2*PI)*sig,3);
real diff=val1-val3;
 printf("\nval1,val2,diff %f %f %f %f\n",val1,val3,diff,sqrt(r2)*idg2);
//if(fabs(diff*pow(sqrt(2*PI)*sig,3))>0.001) printf("\nval1,val2,diff %f %f %f %f\n",val1,val3,diff,sqrt(r2)*idg2);

int dt11=(int) (time12-time11);
int dt12=(int) (time13-time12);
int dt13=(int) (time14-time13);
int dt14=(int) (time15-time14);
// printf("\nlpt_integrate_mol: %d %d %d %d\n",dt11,dt12,dt13,dt14);
// printf("\nfsqrt: %f %f %f %f %f\n",__fsqrt_rd(3.f),__fsqrt_rn(3.f),__fsqrt_ru(3.f),__fsqrt_rz(3.f));
*/

/*
real min_meshsize=min(min(dx,dy),dz);
real max_meshsize=max(max(dx,dy),dz);
//2.0f*sqrt(2.0f*log(2.0f))=2.3548;
real sig= min_meshsize/2.3548f;
real val = exp(-r2/(2.0f*sig*sig));
*/
real idg2=(real) tex1Dfetch(texRefGaussian,LEN_GAUSSIAN_ARRAY-1);
real val=(real) tex1Dfetch(texRefGaussian,r2*idg2);
// printf("\nxr,yr,zr %f %f %f %f\n",xr,yr,zr,sqrt(r2));
return val;
}


__device__ void periodic_grid_position(real &x,real &y,real &z,dom_struct *dom)
{
real xl=dom->xl;real yl=dom->yl;real zl=dom->zl;
real xs=dom->xs;real ys=dom->ys;real zs=dom->zs;
real xe=dom->xe;real ye=dom->ye;real ze=dom->ze;

real xr=x-xs;real yr=y-ys;real zr=z-zs;

if(x<xs||x>xe)  xr=xr-floor(__fdiv_rd(xr,xl))*xl;
if(y<ys||y>ye)  yr=yr-floor(__fdiv_rd(yr,yl))*yl;
if(z<zs||z>ze)  zr=zr-floor(__fdiv_rd(zr,zl))*zl;

x=xs+xr;
y=ys+yr;
z=zs+zr;

}

//Eg: for Gfx system, if i> Gfx.ie-1, then we change it
__device__ void periodic_grid_index(int &ic,int &jc,int &kc,dom_struct *dom, int coordiSys)
{
/*
clock_t time10 = clock();	 
int is,js,ks;
int ie,je,ke;
int in,jn,kn;
clock_t time11 = clock();
*/
int index=coordiSys*48;

int is=tex1Dfetch(texRefDomInfo,index);
int js=tex1Dfetch(texRefDomInfo,index+6);
int ks=tex1Dfetch(texRefDomInfo,index+12);

int ie=tex1Dfetch(texRefDomInfo,index+1);
int je=tex1Dfetch(texRefDomInfo,index+7);
int ke=tex1Dfetch(texRefDomInfo,index+13);

int in=tex1Dfetch(texRefDomInfo,index+2);
int jn=tex1Dfetch(texRefDomInfo,index+8);
int kn=tex1Dfetch(texRefDomInfo,index+14);

//clock_t time10 = clock();
in -=(coordiSys==1);
jn -=(coordiSys==2);
kn -=(coordiSys==3);

/*
//This is not enough for ic~kc, they could be still out of bounds after operation
ic=ic+in*(ic<is);
ic=ic-in*(ic>ie-1);

jc=jc+jn*(jc<js);
jc=jc-jn*(jc>je-1);

kc=kc+kn*(kc<ks);
kc=kc-kn*(kc>ke-1);
*/


/*
if(ic<is||ic > ie-1)  ir= ir -floorf(__fdiv_rd(ir,in))*in;
if(jc<js||jc > je-1)  jr= jr -floorf(__fdiv_rd(jr,jn))*jn;
if(kc<ks||kc > ke-1)  kr= kr -floorf(__fdiv_rd(kr,kn))*kn;
*/

/*
clock_t time11 = clock();
switch(coordiSys)
{
case 1:in-=1;break;
case 2:jn-=1;break;
case 3:kn-=1;break;
default:break;}
clock_t time12 = clock();
int dt10=(int) (time11-time10);
int dt11=(int) (time12-time11);
printf("\ninner periodic_grid: %d %d\n",dt10,dt11);
*/
/* 
clock_t time12 = clock();
switch(coordiSys)
{
case 0:
is=dom->Gcc.is;
js=dom->Gcc.js;
ks=dom->Gcc.ks;

ie=dom->Gcc.ie;
je=dom->Gcc.je;
ke=dom->Gcc.ke;

in=dom->Gcc.in;
jn=dom->Gcc.jn;
kn=dom->Gcc.kn;
break;
case 1:
is=dom->Gfx.is;
js=dom->Gfx.js;
ks=dom->Gfx.ks;

ie=dom->Gfx.ie;
je=dom->Gfx.je;
ke=dom->Gfx.ke;

in=dom->Gfx.in-1;
jn=dom->Gfx.jn;
kn=dom->Gfx.kn;
break;
case 2:
is=dom->Gfy.is;
js=dom->Gfy.js;
ks=dom->Gfy.ks;

ie=dom->Gfy.ie;
je=dom->Gfy.je;
ke=dom->Gfy.ke;

in=dom->Gfy.in;
jn=dom->Gfy.jn-1;
kn=dom->Gfy.kn;
break;
case 3:
is=dom->Gfz.is;
js=dom->Gfz.js;
ks=dom->Gfz.ks;

ie=dom->Gfz.ie;
je=dom->Gfz.je;
ke=dom->Gfz.ke;

in=dom->Gfz.in;
jn=dom->Gfz.jn;
kn=dom->Gfz.kn-1;
break;
default:break;
}
clock_t time13 = clock();	 


int dt10=(int) (time11-time10);
int dt11=(int) (time12-time11);
int dt12=(int) (time13-time12);
printf("\ninner periodic_grid: %d %d %d %d %d %d\n",dt10,dt11,dt12,in1-in,jn1-jn,kn1-kn);
*/

//ic=-65;jc=129;
int ir=ic-is;
int jr=jc-js;
int kr=kc-ks;

//clock_t time13 = clock();	 
if(ic<is||ic > ie-1)  ir= ir -floorf(__fdiv_rd(ir,in))*in;
if(jc<js||jc > je-1)  jr= jr -floorf(__fdiv_rd(jr,jn))*jn;
if(kc<ks||kc > ke-1)  kr= kr -floorf(__fdiv_rd(kr,kn))*kn;
//clock_t time14 = clock();	 
/*
if(ic>jc) ;
clock_t time15 = clock();	 
ks=fminf(1,__exp10f(ic-jc)) ;
clock_t time16 = clock();	 
ks=floorf(3.2f);
clock_t time17 = clock();	 
if(ic==jc) ;
clock_t time18 = clock();	 
ks=fminf(0,__cosf(ic-jc)) ;
clock_t time19 = clock();	 

int dt10=(int) (time11-time10);
int dt11=(int) (time12-time11);
int dt12=(int) (time13-time12);
int dt13=(int) (time14-time13);
printf("\ninner periodic_grid: %d %d %d %d\n",dt10,dt11,dt12,dt13);

int dt14=(int) (time15-time14);
int dt15=(int) (time16-time15);
int dt16=(int) (time17-time16);
int dt17=(int) (time18-time17);
int dt18=(int) (time19-time18);
printf("\nbasic operations: %d %d %d %d %d\n",dt14,dt15,dt16,dt17,dt18);
*/
ic=ir+is;
jc=jr+js;
kc=kr+ks;
}


__device__ void dom_startEnd_index(int &is, int &js,int &ks,int &ie, int &je,int &ke,dom_struct *dom,int coordiSys,int incGhost)
{
/*
int indexDir=3;
int hostDevice=0;
int incGhost=1;
int index=coordiSys*48+hostDevice*24+indexDir*6+3*incGhost;
*/
int index=coordiSys*48+				3*incGhost;

is=tex1Dfetch(texRefDomInfo,index);
ie=tex1Dfetch(texRefDomInfo,index+1);
js=tex1Dfetch(texRefDomInfo,index+6);
je=tex1Dfetch(texRefDomInfo,index+7);
ks=tex1Dfetch(texRefDomInfo,index+12);
ke=tex1Dfetch(texRefDomInfo,index+13);

/*
switch(incGhost)
{
case 0:
switch(coordiSys)
{
case 0:
is=dom->Gcc._is;
ie=dom->Gcc._ie;
js=dom->Gcc._js;
je=dom->Gcc._je;
ks=dom->Gcc._ks;
ke=dom->Gcc._ke;
break;
case 1:
is=dom->Gfx._is;
ie=dom->Gfx._ie;
js=dom->Gfx._js;
je=dom->Gfx._je;
ks=dom->Gfx._ks;
ke=dom->Gfx._ke;
break;
case 2:
is=dom->Gfy._is;
ie=dom->Gfy._ie;
js=dom->Gfy._js;
je=dom->Gfy._je;
ks=dom->Gfy._ks;
ke=dom->Gfy._ke;
break;
case 3:
is=dom->Gfz._is;
ie=dom->Gfz._ie;
js=dom->Gfz._js;
je=dom->Gfz._je;
ks=dom->Gfz._ks;
ke=dom->Gfz._ke;
break;
default: break;
}
break;
case 1:
switch(coordiSys)
{
case 0:
is=dom->Gcc._isb;
ie=dom->Gcc._ieb;
js=dom->Gcc._jsb;
je=dom->Gcc._jeb;
ks=dom->Gcc._ksb;
ke=dom->Gcc._keb;
break;
case 1:
is=dom->Gfx._isb;
ie=dom->Gfx._ieb;
js=dom->Gfx._jsb;
je=dom->Gfx._jeb;
ks=dom->Gfx._ksb;
ke=dom->Gfx._keb;
break;
case 2:
is=dom->Gfy._isb;
ie=dom->Gfy._ieb;
js=dom->Gfy._jsb;
je=dom->Gfy._jeb;
ks=dom->Gfy._ksb;
ke=dom->Gfy._keb;
break;
case 3:
is=dom->Gfz._isb;
ie=dom->Gfz._ieb;
js=dom->Gfz._jsb;
je=dom->Gfz._jeb;
ks=dom->Gfz._ksb;
ke=dom->Gfz._keb;
break;
default: break;
}
break;
default: break;}
*/	
}

//Impose periodic BC, and calculate the grid index from (i,j,k) to i+j*s1b+s2b
//Include ghost boundary
__global__
void calcGridFlowHash_optD(dom_struct *dom,
              int *gridFlowHash,
              int coordiSys)
{

int is,js,ks,ie,je,ke;

//get domain start and end index
int incGhost=1;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);

int   index=coordiSys*48	+21;
int s1b=tex1Dfetch(texRefDomInfo,index);
int s2b=tex1Dfetch(texRefDomInfo,index+1);

for(int k=ks;k<ke;k++)
{
    // subdomain indices
    // shared memory blocks in the subdomain
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;

__syncthreads();

if(i<ie&&i>=is && j<je&&j>=js)
{
//Impose periodic BC, and calculate the grid index from (i,j,k) to i+j*s1b+s2b
//periodic eg: gridFlowHash[isb+j*s1b+k*s2b]=gridFlowHash[ie-2+j*s1b+k*s2b]
//if (i,j,k) is inside the domain,  then gridFlowHash[i+j*s1b+k*s2b]=i+j*s1b+k*s2b
gridFlowHash[i+j*s1b+k*s2b]=calcGridHash(i,j,k,dom,coordiSys);
}


}
}




//Impose periodic BC, and calculate the grid index from (i,j,k) to i+j*s1b+s2b
//calculate the grid index, won't change the ic~kc of the input;
__device__ int calcGridHash(int ic, int jc,int kc,dom_struct *dom,int coordiSys)
{
//int hash1;
//clock_t time10 = clock();	 
periodic_grid_index(ic,jc,kc,dom,coordiSys);
/*
int indexDir=3;
int hostDevice=0;
int incGhost=1;
*/
//clock_t time11 = clock();	 
//int index=coordiSys*48+hostDevice*24+indexDir*6+3*incGhost;
int   index=coordiSys*48	+21;
int s1b=tex1Dfetch(texRefDomInfo,index);
int s2b=tex1Dfetch(texRefDomInfo,index+1);

int hash=ic+jc*s1b+kc*s2b;
return hash;

/*
clock_t time12 = clock();	 

switch(coordiSys)
{
case 0:
hash=ic+jc*dom->Gcc._s1b+kc*dom->Gcc._s2b;break;
case 1:
hash=ic+jc*dom->Gfx._s1b+kc*dom->Gfx._s2b;break;
case 2:
hash=ic+jc*dom->Gfy._s1b+kc*dom->Gfy._s2b;break;
case 3:
hash=ic+jc*dom->Gfz._s1b+kc*dom->Gfz._s2b;break;
default: break;
}

clock_t time13 = clock();	 

int dt10=(int) (time11-time10);
int dt11=(int) (time12-time11);
int dt12=(int) (time13-time12);
printf("\ninner calcGridHash: %d %d %d %d\n",dt10,dt11,dt12,hash-hash1);
*/
}


__device__ real lpt_mol_typeVal(point_struct *points,dom_struct *dom,int pp,int coordiSys,int valType)
{
real Ap;
real cellVol=dom->dx*dom->dy*dom->dz;

switch(coordiSys)
{
case 0:
    switch(valType)
   { case 0:   
//Exchange rate of soluble mass from particles. All the following source should divide cellVol to ensure conservation law
    Ap=-points[pp].msdot;break;
 //    Ap=1;break;
     case 1:
{//Particle volume filter
     real rad=points[pp].r;
     Ap=PI*4/3*rad*rad*rad;
 //    Ap=1;
     break;}
     default:break;
	}
break;
case 1:
//Particle x-momentum reaction to fluid
Ap=-points[pp].Fx;
break;
case 2:
//Particle y-momentum reaction to fluid
Ap=-points[pp].Fy;
break;
case 3:
//Particle z-momentum reaction to fluid
Ap=-points[pp].Fz;
//printf("\nAp %d %f %f %f\n",pp,Ap,cellVol,Ap/cellVol);
break;
default:
break;
}
return Ap/cellVol;
}

__global__ void print_kernel_array_int(int *cell,int lenCell)
{

    int i = blockIdx.x*blockDim.x + threadIdx.x ;
    int j = blockIdx.y*blockDim.y + threadIdx.y ;
    int index=i+j*blockDim.x*gridDim.x;
	
	if(index>=lenCell) return;

        if(cell[index]>0) printf("\nindex cell %d %d\n",index,cell[index]);
        if(index==lenCell-1) printf("\nEnd index cell %d %d\n",index,cell[index]);
}

__global__ void print_kernel_array_real(real *cell,int lenCell)
{

    int i = blockIdx.x*blockDim.x + threadIdx.x ;
    int j = blockIdx.y*blockDim.y + threadIdx.y ;
    int index=i+j*blockDim.x*gridDim.x;

        if(index>=lenCell) return;

        if(index==lenCell-1)     printf("\nindex cell %d %f\n",index,cell[index]);
        //if(cell[index]>EPSILON||cell[index]<-EPSILON)     printf("\nindex cell %d %f\n",index,cell[index]);
	//printf("\nindex cell %d %f\n",index,cell[index]);
}


// xp~zp is the particle location, x(1~N) is the grid face position
// if x(ip)<=xp<x(ip+1) points[pp].i=ip, this is the grid face number rather than grid center number
__global__ void lpt_localize(int npoints, point_struct *points, dom_struct *dom, BC bc)
{
 int pp =  threadIdx.x + blockIdx.x*blockDim.x;

if(pp<npoints)
{
  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;

// Cartesian location of node
real  xp =  points[pp].x;
real  yp =  points[pp].y;
real  zp =  points[pp].z;

//TODO whether periodic BC for point particle need to be determined in future
periodic_grid_position(xp,yp,zp,dom);
 
  int ip,jp,kp;
//x(i)=(i-DOM_BUF)*dom->dx+dom->xs
//i=floor((x(i) - dom->xs) * ddx) + DOM_BUF;
  ip = floor((xp - dom->xs) * ddx) + DOM_BUF;//for x[i]
  jp = floor((yp - dom->ys) * ddy) + DOM_BUF;//for y[j]
  kp = floor((zp - dom->zs) * ddz) + DOM_BUF;//for z[k]

points[pp].i= ip;
points[pp].j= jp;
points[pp].k= kp;
}

}


//Should Add after we have perfect two-way coupling!!
//Note that this subroutine calculate dudt only 1st order in time
__global__ void DvelDt(real *u0,real *u,
real *conv_u, 
real *dudt,
dom_struct *dom, real dt, 
int coordiSys)
{
  // create shared memory
  __shared__ real s_u0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u back
  __shared__ real s_u1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u center

  __shared__ real s_dudt[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
//  __shared__ real s_conv0_u[MAX_THREADS_DIM * MAX_THREADS_DIM];       
  __shared__ real s_conv_u[MAX_THREADS_DIM * MAX_THREADS_DIM];        


// working constants
/*
  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??
  real ab0= 0.5 * a ;
  real ab=1. + ab0 ;

//Weight for Adam-Bashforth interpolation
  real a = dt0/dt;
  a = (a + 2.)/(a + 1.);
  real ab0= a-1.0 ;
  real ab=a;
  */
 
int is,js,ks,ie,je,ke;
//get domain start and end index
int incGhost=1;
dom_startEnd_index(is,js,ks,ie,je,ke,dom,coordiSys,incGhost);
int index=coordiSys*48	+21;
int s1b=tex1Dfetch(texRefDomInfo,index);
int s2b=tex1Dfetch(texRefDomInfo,index+1);

  // loop over u-planes
  for(int i =is; i < ie; i++) {

 // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int j = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int k = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int tj = threadIdx.x;
    int tk = threadIdx.y;

 // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations
    // TODO: THIS CAN BE FIXED BY PADDING ALL OF THESE ARRAYS WHEN COPYING FROM
    // HOST TO DEVICE
    if((k >= ks && k < ke)
      && (j >= js && j < je)) {
      s_u0[tj + tk*blockDim.x] = u0[i + j*s1b + k*s2b];
      s_u1[tj + tk*blockDim.x] = u[i + j*s1b + k*s2b];
//      s_conv0_u[tj + tk*blockDim.x] = conv0_u[i + j*s1b + k*s2b];
      s_conv_u[tj + tk*blockDim.x] = conv_u[i + j*s1b + k*s2b];
    }
   
    __syncthreads();

    real idt=__fdiv_rd(1,dt);
    // Calculate the value at time t_n-0.5, could be more close to t_n instead
    if((tj > 0 && tj < blockDim.x-1) && (tk > 0 && tk < blockDim.y-1)) {

//     s_dudt[tj + tk*blockDim.x]=(s_u1[tj + tk*blockDim.x]-s_u0[tj + tk*blockDim.x])*idt + (s_conv0_u[tj + tk*blockDim.x]+s_conv_u[tj + tk*blockDim.x])/2.f;
     s_dudt[tj + tk*blockDim.x]=(s_u1[tj + tk*blockDim.x]-s_u0[tj + tk*blockDim.x])*idt + s_conv_u[tj + tk*blockDim.x];

    }
   // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((k >= ks && k < ke)
      && (j >= js && j < je)
      && (tj > 0 && tj < (blockDim.x-1))
      && (tk > 0 && tk < (blockDim.y-1))) { 
 	 dudt[i + j*s1b + k*s2b] =s_dudt[tj + tk*blockDim.x];  }
  }
} 





__global__ void stress_u(real rho_f, real nu, real *u0,real *p,real *p0,real *stress, dom_struct *dom,int *flag_u, real dt, real dt0)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_u0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u back
  __shared__ real s_u1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u center
  __shared__ real s_u2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // u forward

  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
  
  __shared__ real grad_P[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
//  __shared__ real grad_P0[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
 
// working constants

/*
  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??
  real ab0= 0.5 * a ;
  real ab=1. + ab0 ;

//Weight for Adam-Bashforth interpolation
  real a = dt0/dt;
  a = (a + 2.)/(a + 1.);
  real ab0= a-1.0 ;
  real ab=a;
  */

  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed


  // loop over u-planes
  for(int i = dom->Gfx._is; i < dom->Gfx._ie; i++) {

 // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int j = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int k = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int tj = threadIdx.x;
    int tk = threadIdx.y;

 // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations
    // TODO: THIS CAN BE FIXED BY PADDING ALL OF THESE ARRAYS WHEN COPYING FROM
    // HOST TO DEVICE
    if((k >= dom->Gfx._ksb && k < dom->Gfx._keb)
      && (j >= dom->Gfx._jsb && j < dom->Gfx._jeb)) {
      s_u0[tj + tk*blockDim.x] = u0[(i-1) + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u1[tj + tk*blockDim.x] = u0[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
      s_u2[tj + tk*blockDim.x] = u0[(i+1) + j*dom->Gfx._s1b + k*dom->Gfx._s2b];
    }
   
    __syncthreads();


    // if off the shared memory block boundary
    if((tj > 0 && tj < blockDim.x-1) && (tk > 0 && tk < blockDim.y-1)) {
      real u011 = s_u0[tj + tk*blockDim.x];
      real u111 = s_u1[tj + tk*blockDim.x];
      real u211 = s_u2[tj + tk*blockDim.x];

      real u101 = s_u1[(tj-1) + tk*blockDim.x];
      real u121 = s_u1[(tj+1) + tk*blockDim.x];

      real u110 = s_u1[tj + (tk-1)*blockDim.x];
      real u112 = s_u1[tj + (tk+1)*blockDim.x];

 // compute diffusion term (Adams-Bashforth stepping)
      real dud1 = (u211 - u111) * ddx;
      real dud0 = (u111 - u011) * ddx;
      real ddudxx = (dud1 - dud0) * ddx;

      dud1 = (u121 - u111) * ddy;
      dud0 = (u111 - u101) * ddy;
      real ddudyy = (dud1 - dud0) * ddy;

      dud1 = (u112 - u111) * ddz;
      dud0 = (u111 - u110) * ddz;
      real ddudzz = (dud1 - dud0) * ddz;

      s_d[tj + tk*blockDim.x] = nu * (ddudxx + ddudyy + ddudzz);
      grad_P[tj + tk*blockDim.x]=abs(flag_u[i + tj*dom->Gfx._s1b
        + tk*dom->Gfx._s2b])
        * ddx * (p[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p[(i-1) + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]);
/*
      grad_P0[tj + tk*blockDim.x]=abs(flag_u[i + tj*dom->Gfx._s1b
        + tk*dom->Gfx._s2b])
        * ddx * (p0[i + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p0[(i-1) + tj*dom->Gcc._s1b + tk*dom->Gcc._s2b]);
  */    

    }
   // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((k >= dom->Gfx._ks && k < dom->Gfx._ke)
      && (j >= dom->Gfx._js && j < dom->Gfx._je)
      && (tj > 0 && tj < (blockDim.x-1))
      && (tk > 0 && tk < (blockDim.y-1))) { 
//     stress[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] =rho_f* s_d[tj + tk*blockDim.x]-ab*grad_P[tj + tk*blockDim.x]+ab0*grad_P0[tj + tk*blockDim.x]-gradP_x;
//     stress[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] =rho_f* s_d[tj + tk*blockDim.x]-ab*grad_P[tj + tk*blockDim.x]+ab0*grad_P0[tj + tk*blockDim.x];
	 stress[i + j*dom->Gfx._s1b + k*dom->Gfx._s2b] =rho_f* s_d[tj + tk*blockDim.x]-grad_P[tj + tk*blockDim.x];
    }
  }
} 


__global__ void stress_v(real rho_f, real nu, real *v0,real *p,real *p0,real *stress, dom_struct *dom,int *flag_v, real dt, real dt0)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_v0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v back
  __shared__ real s_v1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v center
  __shared__ real s_v2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // v forward

  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff

  __shared__ real grad_P[MAX_THREADS_DIM * MAX_THREADS_DIM];       // y-force
// __shared__ real grad_P0[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
 
// working constants
/*  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??
  real ab0= 0.5 * a ;
  real ab=1. + ab0 ;
 real a = dt0/dt;
  a = (a + 2.)/(a + 1.);
  real ab0= a-1.0 ;
  real ab=a;
  */

  // working constants
 
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

  // loop over u-planes
  for(int j = dom->Gfy._js; j < dom->Gfy._je; j++) {
    // subdomain indices
    // the extra 2*blockIdx.X terms implement the necessary overlapping of
    // shared memory blocks in the subdomain
    int k = blockIdx.x*blockDim.x + threadIdx.x - 2*blockIdx.x;
    int i = blockIdx.y*blockDim.y + threadIdx.y - 2*blockIdx.y;
    // shared memory indices
    int tk = threadIdx.x;
    int ti = threadIdx.y;

    // load shared memory
    // TODO: look into the effect of removing these if statements and simply
    // allowing memory overruns for threads that don't matter for particular
    // discretizations
    if((i >= dom->Gfy._isb && i < dom->Gfy._ieb)
      && (k >= dom->Gfy._ksb && k < dom->Gfy._keb)) {
     
      s_v0[tk + ti*blockDim.x] = v0[i + (j-1)*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v1[tk + ti*blockDim.x] = v0[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b];
      s_v2[tk + ti*blockDim.x] = v0[i + (j+1)*dom->Gfy._s1b + k*dom->Gfy._s2b];
     }
   
    // make sure all threads complete shared memory copy
    __syncthreads();

    // compute right-hand side
    // if off the shared memory block boundary
    if((tk > 0 && tk < blockDim.x-1) && (ti > 0 && ti < blockDim.y-1)) {
      // grab the required data points for calculations
      real v101 = s_v0[tk + ti*blockDim.x];
      real v111 = s_v1[tk + ti*blockDim.x];
      real v121 = s_v2[tk + ti*blockDim.x];

      real v110 = s_v1[(tk-1) + ti*blockDim.x];
      real v112 = s_v1[(tk+1) + ti*blockDim.x];
   
      real v011 = s_v1[tk + (ti-1)*blockDim.x];
      real v211 = s_v1[tk + (ti+1)*blockDim.x];

      // compute diffusive term
      real dvd1 = (v211 - v111) * ddx;
      real dvd0 = (v111 - v011) * ddx;
      real ddvdxx = (dvd1 - dvd0) * ddx;

      dvd1 = (v121 - v111) * ddy;
      dvd0 = (v111 - v101) * ddy;
      real ddvdyy = (dvd1 - dvd0) * ddy;

      dvd1 = (v112 - v111) * ddz;
      dvd0 = (v111 - v110) * ddz;
      real ddvdzz = (dvd1 - dvd0) * ddz;

      s_d[tk + ti*blockDim.x] = nu * (ddvdxx + ddvdyy + ddvdzz);
     grad_P[tk + ti*blockDim.x]=abs(flag_v[ti + j*dom->Gfy._s1b
        + tk*dom->Gfy._s2b])
        * ddy * (p[ti + j*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p[ti + (j-1)*dom->Gcc._s1b + tk*dom->Gcc._s2b]);
/* 
      grad_P0[tk + ti*blockDim.x]=abs(flag_v[ti + j*dom->Gfy._s1b
        + tk*dom->Gfy._s2b])
        * ddy * (p0[ti + j*dom->Gcc._s1b + tk*dom->Gcc._s2b]
        - p0[ti + (j-1)*dom->Gcc._s1b + tk*dom->Gcc._s2b]);

*/
     
	}

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((i >= dom->Gfy._is && i < dom->Gfy._ie)
      && (k >= dom->Gfy._ks && k < dom->Gfy._ke)
      && (tk > 0 && tk < (blockDim.x-1))
      && (ti > 0 && ti < (blockDim.y-1))) {
  //    stress[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b] = rho_f*s_d[tk + ti*blockDim.x]-ab*grad_P[tk + ti*blockDim.x]+ab0*grad_P0[tk + ti*blockDim.x]-gradP_y;
     stress[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b] = rho_f*s_d[tk + ti*blockDim.x]-grad_P[tk + ti*blockDim.x];
//    stress[i + j*dom->Gfy._s1b + k*dom->Gfy._s2b] = rho_f*s_d[tk + ti*blockDim.x]-ab*grad_P[tk + ti*blockDim.x]+ab0*grad_P0[tk + ti*blockDim.x];
      }
  }
}

__global__ void stress_w(real rho_f, real nu, real *w0,real *p,real *p0,real *stress, dom_struct *dom,int *flag_w, real dt, real dt0)
{
  // create shared memory
  // no reason to load pressure into shared memory, but leaving it in global
  // will require additional if statements, so keep it in shared
  __shared__ real s_w0[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w back
  __shared__ real s_w1[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w center
  __shared__ real s_w2[MAX_THREADS_DIM * MAX_THREADS_DIM];      // w forward
  
  __shared__ real s_d[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff0
  
  __shared__ real grad_P[MAX_THREADS_DIM * MAX_THREADS_DIM];  // pressure gradient
//  __shared__ real grad_P0[MAX_THREADS_DIM * MAX_THREADS_DIM];       // diff
 
// working constants
/*  real a = (0.5 * dt) / (0.5 * dt0 + 0.5 * dt); //why in time??
  real ab0= 0.5 * a ;
  real ab=1. + ab0 ;
  real a = dt0/dt;
  a = (a + 2.)/(a + 1.);
  real ab0= a-1.0 ;
  real ab=a;

  */
  // working constants
  real ddx = 1. / dom->dx;     // to limit the number of divisions needed
  real ddy = 1. / dom->dy;     // to limit the number of divisions needed
  real ddz = 1. / dom->dz;     // to limit the number of divisions needed

  // loop over w-planes
  for(int k = dom->Gfz._ks; k < dom->Gfz._ke; k++) {
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
    if((j >= dom->Gfz._jsb && j < dom->Gfz._jeb)
      && (i >= dom->Gfz._isb && i < dom->Gfz._ieb)) {
      
      s_w0[ti + tj*blockDim.x] = w0[i + j*dom->Gfz._s1b + (k-1)*dom->Gfz._s2b];
      s_w1[ti + tj*blockDim.x] = w0[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b];
      s_w2[ti + tj*blockDim.x] = w0[i + j*dom->Gfz._s1b + (k+1)*dom->Gfz._s2b];
      
    }
     // make sure all threads complete shared memory copy
    __syncthreads();

    // compute right-hand side
    // if off the shared memory block boundary
    if((ti > 0 && ti < blockDim.x-1) && (tj > 0 && tj < blockDim.y-1)) {
      // grab the required data points for calculations
      real w110 = s_w0[ti + tj*blockDim.x];
      real w111 = s_w1[ti + tj*blockDim.x];
      real w112 = s_w2[ti + tj*blockDim.x];

      real w011 = s_w1[(ti-1) + tj*blockDim.x];
      real w211 = s_w1[(ti+1) + tj*blockDim.x];
   
      real w101 = s_w1[ti + (tj-1)*blockDim.x];
      real w121 = s_w1[ti + (tj+1)*blockDim.x];

      // compute diffusive term
      real dwd1 = (w211 - w111) * ddx;
      real dwd0 = (w111 - w011) * ddx;
      real ddwdxx = (dwd1 - dwd0) * ddx;

      dwd1 = (w121 - w111) * ddy;
      dwd0 = (w111 - w101) * ddy;
      real ddwdyy = (dwd1 - dwd0) * ddy;

      dwd1 = (w112 - w111) * ddz;
      dwd0 = (w111 - w110) * ddz;
      real ddwdzz = (dwd1 - dwd0) * ddz;

      s_d[ti + tj*blockDim.x] = nu * (ddwdxx + ddwdyy + ddwdzz);
      grad_P[ti + tj*blockDim.x] = abs(flag_w[ti + tj*dom->Gfz._s1b
        + k*dom->Gfz._s2b])
        * ddz * (p[ti + tj*dom->Gcc._s1b + k*dom->Gcc._s2b]
        - p[ti + tj*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b]);
/*
      grad_P0[ti + tj*blockDim.x] = abs(flag_w[ti + tj*dom->Gfz._s1b
        + k*dom->Gfz._s2b])
        * ddz * (p0[ti + tj*dom->Gcc._s1b + k*dom->Gcc._s2b]
        - p0[ti + tj*dom->Gcc._s1b + (k-1)*dom->Gcc._s2b]);

*/	
 
    }

    // make sure all threads complete computations
    __syncthreads();

    // copy shared memory back to global
    if((j >= dom->Gfz._js && j < dom->Gfz._je)
      && (i >= dom->Gfz._is && i < dom->Gfz._ie)
      && (ti > 0 && ti < (blockDim.x-1))
      && (tj > 0 && tj < (blockDim.y-1))) {
//      stress[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b]= rho_f* s_d[ti + tj*blockDim.x]-ab*grad_P[ti + tj*blockDim.x]+ab0*grad_P0[ti + tj*blockDim.x]-gradP_z;
stress[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b]= rho_f* s_d[ti + tj*blockDim.x]-grad_P[ti + tj*blockDim.x];
      
//stress[i + j*dom->Gfz._s1b + k*dom->Gfz._s2b]= rho_f* s_d[ti + tj*blockDim.x]-ab*grad_P[ti + tj*blockDim.x]+ab0*grad_P0[ti + tj*blockDim.x];
    }
  }
}



//Interpolate cell-center flow field data A to the particle position, and get Ag at each particle position
__global__ void interpolate_point_scalar_Lag2(int npoints,real *A,real *Ag, point_struct *points, dom_struct *dom)
{
  int pp =  threadIdx.x + blockIdx.x*blockDim.x;
  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;
 
if(pp<npoints)
{

real  x =  points[pp].x;
real  y =  points[pp].y;
real  z =  points[pp].z;
//real Ap=points[pp].mdot;
 
//TODO threat BC in the future

 __syncthreads();
int i,j,k;   
real x1,x2,y1,y2,z1,z2;
int ic,jc,kc;
real wx[2],wy[2],wz[2];
real weight[2][2][2];

  i = floor((x - dom->xs) * ddx- 0.5f) + DOM_BUF;//for xm[i]
  j = floor((y - dom->ys) * ddy- 0.5f) + DOM_BUF;//for ym[j]
  k = floor((z - dom->zs) * ddz- 0.5f) + DOM_BUF;//for zm[k]


  x1 = (i-0.5f) * dom->dx + dom->xs;
  x2 = (i+0.5f) * dom->dx + dom->xs;
  y1=  (j-0.5f) * dom->dy + dom->ys;
  y2=  (j+0.5f) * dom->dy + dom->ys;
  z1=  (k-0.5f) * dom->dz + dom->zs;
  z2=  (k+0.5f) * dom->dz + dom->zs;

//2nd order lagragian interpolation
 wx[0]=(x2-x)/(x2-x1);
 wx[1]=(x-x1)/(x2-x1);
 wy[0]=(y2-y)/(y2-y1);
 wy[1]=(y-y1)/(y2-y1);
 wz[0]=(z2-z)/(z2-z1);
 wz[1]=(z-z1)/(z2-z1);

real buf=0.f;
//Normalize the weight
 for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
      { 
	weight[ii][jj][kk]=wx[ii]*wy[jj]*wz[kk];	
	buf+=weight[ii][jj][kk];
	}

//Normalize
if(fabs(buf-1)>EPSILON)
{
 for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
	weight[ii][jj][kk]=weight[ii][jj][kk]/buf;	
}

Ag[pp]=0;
//Add them to the scalar field source
 for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
      { 
  ic=i+ii;jc=j+jj;kc=k+kk;
periodic_grid_index(ic,jc,kc,dom,0);
  Ag[pp]+=weight[ii][jj][kk]*A[ic +jc*dom->Gcc.s1b + kc*dom->Gcc.s2b];

//if(pp==1) printf("\nweight %f %f %d %d %d\n",weight[ii][jj][kk],A[ic +jc*dom->Gcc.s1b + kc*dom->Gcc.s2b],ic,jc,kc);
 }

//TODO add mask infomation as in lpt_mollify_sc in lpt_interpolator.f90
  }
}

//Bilinear interpolation!  This skeme is accurate enough, which has been tested by oscllating flows.
//TODO treat boudary specification,such as mask information
__global__ void interpolate_point_vel_Lag2(real *u, real *v, real *w,int npoints,real rho_f, real nu,real *ug,real *vg,real *wg, point_struct *points, dom_struct *dom, BC bc)
{
  // int node = threadIdx.x;
  int pp =  threadIdx.x + blockIdx.x*blockDim.x;
  real ddx = 1. / dom->dx;
  real ddy = 1. / dom->dy;
  real ddz = 1. / dom->dz;


if(pp>=npoints) return;

// Cartesian location of node
real  x =  points[pp].x;
real  y =  points[pp].y;
real  z =  points[pp].z;

// periodic_grid_position(x,y,z,dom);

  __syncthreads();

int i,j,k;   
real x1,x2,y1,y2,z1,z2;
int ic,jc,kc;
real wx[2],wy[2],wz[2];
real weight[2][2][2];

  // interpolate velocities
  // interpolate u-velocity
  i = floor((x - dom->xs) * ddx) + DOM_BUF; 	//for x[i]
  j = floor((y - dom->ys) * ddy- 0.5f) + DOM_BUF;//for ym[j]
  k = floor((z - dom->zs) * ddz- 0.5f) + DOM_BUF;//for zm[k]

//periodic_grid_index(ic,jc,kc,dom,1);

  x1 = (i-DOM_BUF) * dom->dx + dom->xs;
  x2 = (i+1-DOM_BUF) * dom->dx + dom->xs;
  y1=  (j-0.5f) * dom->dy + dom->ys;
  y2=  (j+0.5f) * dom->dy + dom->ys;
  z1=  (k-0.5f) * dom->dz + dom->zs;
  z2=  (k+0.5f) * dom->dz + dom->zs;

 
 wx[0]=(x2-x)/(x2-x1);
 wx[1]=(x-x1)/(x2-x1);
 wy[0]=(y2-y)/(y2-y1);
 wy[1]=(y-y1)/(y2-y1);
 wz[0]=(z2-z)/(z2-z1);
 wz[1]=(z-z1)/(z2-z1);

ug[pp]=0;
  for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
      { 
	
weight[ii][jj][kk]=wx[ii]*wy[jj]*wz[kk];	
ic=i+ii;
jc=j+jj;
kc=k+kk;
periodic_grid_index(ic,jc,kc,dom,1);
	ug[pp]+=weight[ii][jj][kk]*u[ic +jc*dom->Gfx.s1b + kc*dom->Gfx.s2b];
//if(isnan(ug[pp])) printf("\ninterp_u %f %f %d %d %d\n",weight[ii][jj][kk],u[ic +jc*dom->Gfx.s1b + kc*dom->Gfx.s2b],ic,jc,kc);
	}

//if(isnan(ug[pp])) printf("\nposition %f %f %f %f %f %f %f %f %f\n",x2,x1,y2,y1,z2,z1,dom->dx,dom->dy,dom->dz);

  // interpolate V-velocity
  i = floor((x - dom->xs) * ddx- 0.5f) + DOM_BUF;
  j = floor((y - dom->ys) * ddy) + DOM_BUF;
  k = floor((z - dom->zs) * ddz- 0.5f) + DOM_BUF;

//periodic_grid_index(ic,jc,kc,dom,2);

  x1 = (i-0.5f) * dom->dx + dom->xs;
  x2 = (i+0.5f) * dom->dx + dom->xs;
  y1=  (j-DOM_BUF) * dom->dy + dom->ys;
  y2=  (j+1-DOM_BUF) * dom->dy + dom->ys;
  z1=  (k-0.5f) * dom->dz + dom->zs;
  z2=  (k+0.5f) * dom->dz + dom->zs;


 wx[0]=(x2-x)/(x2-x1);
 wx[1]=(x-x1)/(x2-x1);
 wy[0]=(y2-y)/(y2-y1);
 wy[1]=(y-y1)/(y2-y1);
 wz[0]=(z2-z)/(z2-z1);
 wz[1]=(z-z1)/(z2-z1);

vg[pp]=0;
  for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
      { 
	weight[ii][jj][kk]=wx[ii]*wy[jj]*wz[kk];
ic=i+ii;
jc=j+jj;
kc=k+kk;
periodic_grid_index(ic,jc,kc,dom,2);
	vg[pp]+=weight[ii][jj][kk]*v[ic +jc*dom->Gfy.s1b + kc*dom->Gfy.s2b];
	}


  // interpolate W-velocity
  i = floor((x - dom->xs) * ddx- 0.5f) + DOM_BUF;
  j = floor((y - dom->ys) * ddy- 0.5f) + DOM_BUF;
  k = floor((z - dom->zs) * ddz) + DOM_BUF;
//periodic_grid_index(ic,jc,kc,dom,3);

  x1 = (i-0.5f) * dom->dx + dom->xs;
  x2 = (i+0.5f) * dom->dx + dom->xs;
  y1=  (j-0.5f) * dom->dy + dom->ys;
  y2=  (j+0.5f) * dom->dy + dom->ys;
  z1=  (k-DOM_BUF) * dom->dz + dom->zs;
  z2=  (k+1-DOM_BUF) * dom->dz + dom->zs;


 wx[0]=(x2-x)/(x2-x1);
 wx[1]=(x-x1)/(x2-x1);
 wy[0]=(y2-y)/(y2-y1);
 wy[1]=(y-y1)/(y2-y1);
 wz[0]=(z2-z)/(z2-z1);
 wz[1]=(z-z1)/(z2-z1);

wg[pp]=0;

  for(int kk = 0; kk < 2; kk++) 
   for(int jj = 0; jj < 2; jj++) 
    for(int ii = 0; ii < 2; ii++) 
      { 
	weight[ii][jj][kk]=wx[ii]*wy[jj]*wz[kk];
ic=i+ii;
jc=j+jj;
kc=k+kk;
periodic_grid_index(ic,jc,kc,dom,3);
	wg[pp]+=weight[ii][jj][kk]*w[ic + jc*dom->Gfz.s1b + kc*dom->Gfz.s2b];
//if(fabs(w[ic + jc*dom->Gfz.s1b + kc*dom->Gfz.s2b])>0.1f) printf("\nw %d %d %d %f\n",ic,jc,kc,w[ic + jc*dom->Gfz.s1b + kc*dom->Gfz.s2b]);
	}
//if(pp==0) printf("\nwg %f %f %f %f\n",x,y,z,wg[pp]);
 }


__global__ void point_interp_init(int npoints,point_struct *points,real *ug,real *vg,real *wg)
{
 int pp = threadIdx.x + blockIdx.x*blockDim.x;
  if(pp >= npoints) return;
//real rad=  points[pp].r;
ug[pp]=0;
vg[pp]=0;
wg[pp]=0;


}

/*
__global__ void point_interp_init(int npoints,point_struct *points,real *ug,real *vg,real *wg,real *lpt_stress_u,real *lpt_stress_v,real *lpt_stress_w,real *scg)
{
 int pp = threadIdx.x + blockIdx.x*blockDim.x;
  if(pp >= npoints) return;
//real rad=  points[pp].r;
ug[pp]=0;
vg[pp]=0;
wg[pp]=0;

//if(pp==npoints-1) printf("\nug[pp] %f\n",ug[pp]); 

lpt_stress_u[pp]=0;
lpt_stress_v[pp]=0;
lpt_stress_w[pp]=0;

scg[pp]=0;
//volPoint[pp]=PI*4/3*rad*rad*rad;

}
*/

//Note for this method, C_drag has to be greater than 0!!
__global__ void drag_move_bubbles(point_struct *points,dom_struct *dom, int npoints,
real *ug,real *vg,real *wg,
real *lpt_stress_u,real *lpt_stress_v,real *lpt_stress_w,real *scg,
real rho_f,real mu, g_struct g,gradP_struct gradP,
real sc_eq,real DIFF,real dt)
//gradP serve as bodyforce for the time being
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;

  if(pp >= npoints) return;

real up=points[pp].u;
real vp=points[pp].v;
real wp=points[pp].w;

real xp=points[pp].x;
real yp=points[pp].y;
real zp=points[pp].z;

real dia=2*points[pp].r;
real idia=1/dia;
real rhod=points[pp].rho;//rhod is the particle density

//fluid velocity at the particle position
real uf=ug[pp];
real vf=vg[pp];
real wf=wg[pp];

//realative velocity between point_particle and fluid
real ur=sqrt((up-uf)*(up-uf)+(vp-vf)*(vp-vf)+(wp-wf)*(wp-wf));
real nu=mu/rho_f;

//Particle Reynolds number
real Rep=ur*dia/nu+EPSILON;

//Based on hp*dp/D=2+0.6Re_p^0.5 Sc^{1/3}, ref eq (22) in Oresta&&Prosperetti(2014)
real Nu =2+0.6*sqrt(Rep)*powf(nu/DIFF,1.0/3.0);
points[pp].Nu=Nu;

//Nu=2;
real hp =Nu*DIFF/dia;

//Modification to stokes theory when Rep>1
real F = 1.0f+Rep/(8+0.5f*(Rep+3.315f*sqrt(Rep)));
real taud = dia*dia/(24.0f*nu);
real itau=F/taud;

//if(pp==1) printf("\nitau %f %f %f %f\n",Rep,F,taud,itau);

real volume=1./6. * PI * dia*dia*dia;
//real ivolume=1/volume;

//Including soluble part:rhod*volume  and insoluble part: ms
real mp =  rhod *volume;
points[pp].ms=mp;

/*
Exchange rate of soluble mass into scalar field
dms/dt = pi *dp^2*hp*(rho_s-rho_{sat})   ref eq (20) in Oresta&&Prosperetti(2014)
*/
//real kappa=3.2e-5;
//real z0=-1615*1000;//TODO tricky to determin !!!
//TODO This is tricky, what's the altitute for initialization
//Record the accleration of bubble radius basically
//real msdot,sc_sat;
real msdot;

if(points[pp].r>0) 
{
//sc_sat=kappa*rho_f*g.z*(zp+z0);
msdot= PI*hp*dia*dia*(scg[pp]-sc_eq);
}
else   msdot= 0;


//fluid mass in particle volume
//real mf =  rho_f *volume;
real ddiadt=2*msdot/rhod/PI/dia/dia;


//drag acc on particle
real drag_x=(uf-up)*itau;
real drag_y=(vf-vp)*itau;
real drag_z=(wf-wp)*itau;

//Fluid stress at the particle position
real stress_x=lpt_stress_u[pp];
real stress_y=lpt_stress_v[pp];
real stress_z=lpt_stress_w[pp];

//idvdt=1/volume *d(volume)/dt 
real idvdt=3*ddiadt*idia;

//Virtual buoyancy term
real virt_buoy_x=(uf-up)*idvdt;
real virt_buoy_y=(vf-vp)*idvdt;
real virt_buoy_z=(wf-wp)*idvdt;


//Total add mass,  -gradP.x~z are the body force on fluid apart from gravity
/*
real dudt=(stress_x/rho_f-gradP.x);
real dvdt=(stress_y/rho_f-gradP.y);
real dwdt=(stress_z/rho_f-gradP.z);
*/
real dudt=stress_x;
real dvdt=stress_y;
real dwdt=stress_z;


//Store the fluid force on particle including add mass, drag, fluid force and gravity
//default C_add=0.5; C_stress=1; C_drag=1;
real Fx=(3*dudt+drag_x+virt_buoy_x)-2*g.x;
real Fy=(3*dvdt+drag_y+virt_buoy_y)-2*g.y;
real Fz=(3*dwdt+drag_z+virt_buoy_z)-2*g.z;

real itaub=itau+idvdt;
real taub=1/itaub;
real rhs_x=Fx+ up*itaub;
real rhs_y=Fy+ vp*itaub;
real rhs_z=Fz+ wp*itaub;

     //Treat taup as very small and integrate directly!!
      real ratio=exp(-dt*itaub);
	
      real termVel_x=rhs_x*taub;
      real termVel_y=rhs_y*taub;
      real termVel_z=rhs_z*taub;

     // update linear velocities
      real u1 = ratio*(up -termVel_x) + termVel_x;
      real v1 = ratio*(vp -termVel_y) + termVel_y;
      real w1 = ratio*(wp -termVel_z) + termVel_z;


      points[pp].u = u1;
      points[pp].v = v1;
      points[pp].w = w1;

    // update position 
      real x1 = xp + (up-termVel_x)*(1-ratio)/itaub+termVel_x*dt;
      real y1 = yp + (vp-termVel_y)*(1-ratio)/itaub+termVel_y*dt;
      real z1 = zp + (wp-termVel_z)*(1-ratio)/itaub+termVel_z*dt;

      points[pp].x = x1;
      points[pp].y = y1;
      points[pp].z = z1;
/*
Store the fluid Momentum force on particle including add mass, drag, fluid force and gravity
added mass effect on particle-fluid force interaction
gravity is not implemented on fluid, but we need to add -mf*g to reaction force
F=d(mp*u)/dt -mp*g, according to eq 6 from dropDiffusionForceImpleBlue.pdf
*/
//This is in fact momentum back reaction in particle sub-timestep, it will be accumulated during the sub-step, and then been added to fluid source force.
/*
points[pp].Fx=(rhs_x-mp*g.x)*dt-(C_add*(u1-up)*mf+C_drag*(x1-xp)*mp*itau);
points[pp].Fy=(rhs_y-mp*g.y)*dt-(C_add*(v1-vp)*mf+C_drag*(y1-yp)*mp*itau);
points[pp].Fz=(rhs_z-mp*g.z)*dt-(C_add*(w1-wp)*mf+C_drag*(z1-zp)*mp*itau);
*/
     // update soluble mass
      points[pp].ms = points[pp].ms0 + msdot * dt;
      points[pp].msdot=msdot;

     // update bubble radius
      points[pp].r = (dia + ddiadt * dt)/2.f;
  
     // update particle time step
      points[pp].dt =  dt;


//TODO periodic BC for particles, may need to change in future
periodic_grid_position(points[pp].x,points[pp].y,points[pp].z,dom);

//update old values
      points[pp].x0 = points[pp].x;
      points[pp].y0 = points[pp].y;
      points[pp].z0 = points[pp].z;

      points[pp].u0 = points[pp].u;
      points[pp].v0 = points[pp].v;
      points[pp].w0 = points[pp].w;

      points[pp].ms0 = points[pp].ms;

}





//TODO Needs to change rho_sat!!
__global__ void drag_points(point_struct *points, int npoints,
real *ug,real *vg,real *wg,
real *lpt_stress_u,real *lpt_stress_v,real *lpt_stress_w,real *scg,
real rho_f,real mu, g_struct g,gradP_struct gradP,
real C_add,real C_stress,real C_drag,
real sc_eq,real DIFF)
//gradP serve as bodyforce for the time being
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;

  if(pp >= npoints) return;
real up=points[pp].u;
real vp=points[pp].v;
real wp=points[pp].w;
real dia=2*points[pp].r;
real rhod=points[pp].rho;//rhod is the particle density

//particle interaction force
real iFx=points[pp].iFx;
real iFy=points[pp].iFy;
real iFz=points[pp].iFz;

//fluid velocity at the particle position
real uf=ug[pp];
real vf=vg[pp];
real wf=wg[pp];



//realative velocity between point_particle and fluid
real ur=sqrt((up-uf)*(up-uf)+(vp-vf)*(vp-vf)+(wp-wf)*(wp-wf));
real nu=mu/rho_f;

//Particle Reynolds number
real Rep=ur*dia/nu+EPSILON;


//Based on hp*dp/D=2+0.6Re_p^0.5 Sc^{1/3}, ref eq (22) in Oresta&&Prosperetti(2014)
real Nu =2+0.6*sqrt(Rep)*powf(nu/DIFF,1.0/3.0);

points[pp].Nu=Nu;

//Nu=2;
real hp =Nu*DIFF/dia;
/*
real hh =(0.6*sqrt(Rep)*powf(nu/DIFF,1.0/3.0));
int n3=5;
int n2=1;
real r=n3*DIFF/PI/2;
real strength=r/n2;
if(ur>0.5*strength)  printf("\nhh %f %f %f %f\n",hh,Rep,ur,DIFF);
*/

//Modification to stokes theory when Rep>1
real F = 1.0f+0.15f*powf(Rep,0.687f); 
//real F = 1.0f; 
real taud = rhod*dia*dia/(18.0f*mu);
real itau=F/taud;


real volume=1./6. * PI * dia*dia*dia;

//Including soluble part:rhod*volume  and insoluble part: ms
real mp =  rhod *volume + points[pp].ms;
//fluid mass in particle volume
real mf =  rho_f *volume;
real msdot=points[pp].msdot;
real gammar=mp/mf;

//drag force on particle
real drag_x=(uf-up)*itau*mp;
real drag_y=(vf-vp)*itau*mp;
real drag_z=(wf-wp)*itau*mp;

//Fluid stress at the particle position
real stress_x=lpt_stress_u[pp];
real stress_y=lpt_stress_v[pp];
real stress_z=lpt_stress_w[pp];

//Total add mass,  -gradP.x~z are the body force on fluid apart from gravity
real add_x=(stress_x/rho_f-gradP.x)*mf;
real add_y=(stress_y/rho_f-gradP.y)*mf;
real add_z=(stress_z/rho_f-gradP.z)*mf;


//Store the fluid force on particle including add mass, drag, fluid force and gravity
//default C_add=0.5; C_stress=1; C_drag=1;
real Fx=(C_add*add_x+C_stress*stress_x*volume+C_drag*drag_x)+(mp-mf)*g.x;
real Fy=(C_add*add_y+C_stress*stress_y*volume+C_drag*drag_y)+(mp-mf)*g.y;
real Fz=(C_add*add_z+C_stress*stress_z*volume+C_drag*drag_z)+(mp-mf)*g.z;


//Store the temp particle acceleration, also including the particle soluble mass change  in the last term!
real udot =(Fx+ iFx -up*msdot)/mp;
real vdot =(Fy+ iFy -vp*msdot)/mp;
real wdot =(Fz+ iFz -wp*msdot)/mp;

//printf("\nFx iFx up msdot,C_add,gammar %f %f %f %f %f %f\n",Fx,iFx,up,msdot,C_add,gammar);

//acount for added mass effect since it appears also on the left handside of particle governing equation
if(fabs(C_add-0.5f)<EPSILON)
{
      udot = udot/(1+C_add/gammar);
      vdot = vdot/(1+C_add/gammar);
      wdot = wdot/(1+C_add/gammar);
}

//particle acceleration
      points[pp].udot = udot;
      points[pp].vdot = vdot;
      points[pp].wdot = wdot;

//if(pp==1) printf("\nitau %f %f %f %f\n",Rep,F,taud,itau);
//if(pp==1) printf("\ndrag %f %f %f %f %f %f %f\n",drag_z,itau,mp-mf,wp,wf,wdot,Fz);
//if(pp==1) printf("\n%f %f %f %f %f %f %f %f %f\n",add_x,stress_x,drag_x,itau,mp-mf,up,uf,udot,Fx);

/*
Store the fluid force on particle including add mass, drag, fluid force and gravity
added mass effect on particle-fluid force interaction
gravity is not implemented on fluid, but we need to add -mf*g to reaction force
F=d(mp*u)/dt -mp*g, according to eq 6 from dropDiffusionForceImpleBlue.pdf
*/
points[pp].Fx=Fx-C_add*udot*mf-mp*g.x;
points[pp].Fy=Fy-C_add*vdot*mf-mp*g.y;
points[pp].Fz=Fz-C_add*wdot*mf-mp*g.z;


/*
Exchange rate of soluble mass into scalar field
dms/dt = pi *dp^2*hp*(rho_s-rho_{sat})   ref eq (20) in Oresta&&Prosperetti(2014)
*/
//printf("\nmsdot %f %f %f %f\n",points[pp].msdot,hp,scg[pp],sc_eq);

if(points[pp].ms>0)  points[pp].msdot= PI*dia*dia*hp*(scg[pp]- sc_eq);
else   points[pp].msdot= 0;
}



//update point velocity -> 1st step of Eulerian prediction
__global__ void move_points_a(point_struct *points, int npoints,
  real dt)
{

  int pp = threadIdx.x + blockIdx.x*blockDim.x; // point_point_particle number
  
  //real m = 4./3. * PI * points[pp].rho * points[pp].r*points[pp].r*points[pp].r;
  //real dT = dt / dt0;

  if(pp < npoints) {
    // update position
      points[pp].x = points[pp].x0 + 0.5*points[pp].u * dt;
      points[pp].y = points[pp].y0 + 0.5*points[pp].v * dt;
      points[pp].z = points[pp].z0 + 0.5*points[pp].w * dt;
//printf("\nmove_points_a %f %f %f %f\n",points[pp].x,points[pp].x0,points[pp].u,dt);
     // update linear velocities
      points[pp].u = points[pp].u0 + 0.5*points[pp].udot * dt;
      points[pp].v = points[pp].v0 + 0.5*points[pp].vdot * dt;
      points[pp].w = points[pp].w0 + 0.5*points[pp].wdot * dt;

     // update soluble mass
      points[pp].ms = points[pp].ms0 + 0.5*points[pp].msdot * dt;
  }
}

//update point velocity -> 2nd step of Eulerian prediction
__global__ void move_points_b(dom_struct *dom,point_struct *points, int npoints,
  real dt)
{

  int pp = threadIdx.x + blockIdx.x*blockDim.x; // point_point_particle number
  
  //real m = 4./3. * PI * points[pp].rho * points[pp].r*points[pp].r*points[pp].r;
  //real dT = dt / dt0;

  if(pp < npoints) {
    // update position
      points[pp].x = points[pp].x0 + points[pp].u * dt;
      points[pp].y = points[pp].y0 + points[pp].v * dt;
      points[pp].z = points[pp].z0 + points[pp].w * dt;

     // update linear velocities
      points[pp].u = points[pp].u0 + points[pp].udot * dt;
      points[pp].v = points[pp].v0 + points[pp].vdot * dt;
      points[pp].w = points[pp].w0 + points[pp].wdot * dt;

     // update soluble mass
      points[pp].ms = points[pp].ms0 + points[pp].msdot * dt;
  

     // update particle time step
      points[pp].dt =  dt;


//TODO periodic BC for particles, may need to change in future
periodic_grid_position(points[pp].x,points[pp].y,points[pp].z,dom);

//update old values
      points[pp].x0 = points[pp].x;
      points[pp].y0 = points[pp].y;
      points[pp].z0 = points[pp].z;

      points[pp].u0 = points[pp].u;
      points[pp].v0 = points[pp].v;
      points[pp].w0 = points[pp].w;

      points[pp].ms0 = points[pp].ms;
  }
}


__global__ void points_vel_square(point_struct *points, real *vel, int npoints)
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;
  if(pp>=npoints) return;

real up=points[pp].u;
real vp=points[pp].v;
real wp=points[pp].w;

  vel[pp]=up*up+vp*vp+wp*wp;

}




//The droplet mass is constant for momentum in this case
__global__ void drag_move_points(point_struct *points,dom_struct *dom, int npoints,
real *ug,real *vg,real *wg,
real *lpt_stress_u,real *lpt_stress_v,real *lpt_stress_w,real *scg,
real rho_f,real mu, g_struct g,gradP_struct gradP,
real C_add,real C_stress,real C_drag,
real sc_eq,real DIFF,real dt)
//gradP serve as bodyforce for the time being
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;

  if(pp >= npoints) return;

real up=points[pp].u;
real vp=points[pp].v;
real wp=points[pp].w;

real xp=points[pp].x;
real yp=points[pp].y;
real zp=points[pp].z;

real dia=2*points[pp].r;
real rhod=points[pp].rho;//rhod is the particle density


//fluid velocity at the particle position
real uf=ug[pp];
real vf=vg[pp];
real wf=wg[pp];

//realative velocity between point_particle and fluid
real ur=sqrt((up-uf)*(up-uf)+(vp-vf)*(vp-vf)+(wp-wf)*(wp-wf));
real nu=mu/rho_f;

//Particle Reynolds number
real Rep=ur*dia/nu+EPSILON;


//Modification to stokes theory when Rep>1
real F = 1.0f+0.15f*powf(Rep,0.687f); 
//real F = 1.0f; 
real taud = rhod*dia*dia/(18.0f*mu);
real itau=F/taud;

//if(pp==1) printf("\nitau %f %f %f %f\n",Rep,F,taud,itau);

real volume=1./6. * PI * dia*dia*dia;

//Including soluble part:rhod*volume  and insoluble part: ms
real ms =  points[pp].ms;
//real mp =  rhod *volume + ms;
real mp =  rhod *volume;
//fluid mass in particle volume
real mf =  rho_f *volume;
//real gammar=mp/mf;

//drag force on particle
real drag_x=(uf-up)*itau*mp;
real drag_y=(vf-vp)*itau*mp;
real drag_z=(wf-wp)*itau*mp;

//if(pp==0) printf("\ntau %f %f %f\n",1.f/itau,taud,Rep);
//if(pp==0) printf("\nflow vel %f %f %f\n",wf,wp,drag_z);
//Fluid stress at the particle position

real stress_x=lpt_stress_u[pp];
real stress_y=lpt_stress_v[pp];
real stress_z=lpt_stress_w[pp];


//Total add mass,  -gradP.x~z are the body force on fluid apart from gravity
real dudt=(stress_x/rho_f-gradP.x);
real dvdt=(stress_y/rho_f-gradP.y);
real dwdt=(stress_z/rho_f-gradP.z);

real add_x=dudt*mf;
real add_y=dvdt*mf;
real add_z=dwdt*mf;

//Store the fluid force on particle including add mass, drag, fluid force and gravity
//default C_add=0.5; C_stress=1; C_drag=1;
real Fx=(C_add*add_x+C_stress*stress_x*volume+C_drag*drag_x)+(mp-mf)*g.x;
real Fy=(C_add*add_y+C_stress*stress_y*volume+C_drag*drag_y)+(mp-mf)*g.y;
real Fz=(C_add*add_z+C_stress*stress_z*volume+C_drag*drag_z)+(mp-mf)*g.z;


real rhs_x=Fx+ C_drag*up*itau*mp;
real rhs_y=Fy+ C_drag*vp*itau*mp;
real rhs_z=Fz+ C_drag*wp*itau*mp;

//if(pp==0) printf("\n Fz %f %f %f %f %f %f\n",add_z,stress_z*volume,drag_z,(mp-mf)*g.z,Fz,rhs_z);



//Exchange rate of soluble mass into scalar field dms/dt = pi *dp^2*hp*(rho_s-rho_{sat})   ref eq (20) in Oresta&&Prosperetti(2014)


//Based on hp*dp/D=2+0.6Re_p^0.5 Sc^{1/3}, ref eq (22) in Oresta&&Prosperetti(2014)
real Nu =2+0.6*sqrt(Rep)*powf(nu/DIFF,1.0/3.0);
points[pp].Nu=Nu;
//Nu=2;
real hp =Nu*DIFF/dia;

//sc_eq=log(K_{ow}), and rho_{sat}= \kappa *ms/vp, where \kappa=rho_f/rhod*1/K_{ow}
real kappa=rho_f/rhod/__exp10f(sc_eq);
//real kappa=rho_f/rhod/__expf(sc_eq);
real rho_sat=kappa*ms/volume;
real msdot= PI*dia*dia*hp*(scg[pp]- rho_sat);


 
     //Treat taup as very small and integrate directly!!
     // real acc_m=points[pp].msdot+C_drag*mp*itau;
      real acc_m=C_drag*mp*itau;
      real total_m=mp+C_add*mf; 
      real idt=acc_m/total_m;
      real ratio=exp(-dt*idt);
	
      real x1,y1,z1;
      real u1,v1,w1;

 if(fabs(acc_m)>EPSILON)
{
      real termVel_x=rhs_x/acc_m;
      real termVel_y=rhs_y/acc_m;
      real termVel_z=rhs_z/acc_m;
//if(isnan(acc_m)||isnan(total_m)) printf("\nacc_m %f %f %f %f %d\n",points[pp].msdot,hp,scg[pp],sc_eq,pp);

 
     // update linear velocities
       u1 = ratio*(up -termVel_x) + termVel_x;
       v1 = ratio*(vp -termVel_y) + termVel_y;
       w1 = ratio*(wp -termVel_z) + termVel_z;


//if(fabs(points[pp].u)>100||fabs(points[pp].v)>100||fabs(points[pp].w)>100) printf("\nrhs %f %f %f\n",rhs_x,total_m,acc_m);
//if(fabs(points[pp].u)>100||fabs(points[pp].v)>100||fabs(points[pp].w)>100) printf("\nu0~w0 %f %f %f %f\n",rhs_x/total_m,rhs_x/acc_m,points[pp].u0,dt);

    // update position 
       x1 = xp + (up-termVel_x)*(1-ratio)/idt+termVel_x*dt;
       y1 = yp + (vp-termVel_y)*(1-ratio)/idt+termVel_y*dt;
       z1 = zp + (wp-termVel_z)*(1-ratio)/idt+termVel_z*dt;

}

else
{
real acc_x=rhs_x/total_m;
real acc_y=rhs_y/total_m;
real acc_z=rhs_z/total_m;

 u1=acc_x*dt+up;
 v1=acc_y*dt+vp;
 w1=acc_z*dt+wp;

 x1=xp+(0.5f *acc_x*dt+up)*dt;
 y1=yp+(0.5f *acc_y*dt+vp)*dt;
 z1=zp+(0.5f *acc_z*dt+wp)*dt;

}

      points[pp].u = u1;
      points[pp].v = v1;
      points[pp].w = w1;

      points[pp].x = x1;
      points[pp].y = y1;
      points[pp].z = z1;
 
//Store the fluid Momentum force on particle including add mass, drag, fluid force and gravity added mass effect on particle-fluid force interaction  gravity is not implemented on fluid, but we need to add -mf*g to reaction force; F=d(mp*u)/dt -mp*g, according to eq 6 from dropDiffusionForceImpleBlue.pdf


//This is in fact momentum back reaction in particle sub-timestep, it will be accumulated during the sub-step, and then been added to fluid source force.
//Usually we use this one

points[pp].Fx=rhs_x*dt-(C_add*(u1-up)*mf+C_drag*(x1-xp)*mp*itau);
points[pp].Fy=rhs_y*dt-(C_add*(v1-vp)*mf+C_drag*(y1-yp)*mp*itau);
points[pp].Fz=rhs_z*dt-(C_add*(w1-wp)*mf+C_drag*(z1-zp)*mp*itau);
/*
points[pp].Fx=(rhs_x-mp*g.x)*dt-(C_add*(u1-up)*mf+C_drag*(x1-xp)*mp*itau);
points[pp].Fy=(rhs_y-mp*g.y)*dt-(C_add*(v1-vp)*mf+C_drag*(y1-yp)*mp*itau);
points[pp].Fz=(rhs_z-mp*g.z)*dt-(C_add*(w1-wp)*mf+C_drag*(z1-zp)*mp*itau);
*/


///TODO periodic BC for particles, may need to change in future
periodic_grid_position(points[pp].x,points[pp].y,points[pp].z,dom);
//update old values
      points[pp].x0 = points[pp].x;
      points[pp].y0 = points[pp].y;
      points[pp].z0 = points[pp].z;

      points[pp].u0 = points[pp].u;
      points[pp].v0 = points[pp].v;
      points[pp].w0 = points[pp].w;

     // update soluble mass
      points[pp].ms = points[pp].ms0 + msdot * dt;
      points[pp].msdot = msdot;
  
     // update particle time step
      points[pp].dt =  dt;

      points[pp].ms0 = points[pp].ms;

}

/*
//Note for this method, C_drag has to be greater than 0!!
//The relaxed condition would be that C_drag*dms/dt can't be 0 together!
__global__ void drag_move_points(point_struct *points,dom_struct *dom, int npoints,
real *ug,real *vg,real *wg,
real *lpt_stress_u,real *lpt_stress_v,real *lpt_stress_w,real *scg,
real rho_f,real mu, g_struct g,gradP_struct gradP,
real C_add,real C_stress,real C_drag,
real sc_eq,real DIFF,real dt)
//gradP serve as bodyforce for the time being
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;

  if(pp >= npoints) return;

real up=points[pp].u;
real vp=points[pp].v;
real wp=points[pp].w;

real xp=points[pp].x;
real yp=points[pp].y;
real zp=points[pp].z;

real dia=2*points[pp].r;
real rhod=points[pp].rho;//rhod is the particle density


//fluid velocity at the particle position
real uf=ug[pp];
real vf=vg[pp];
real wf=wg[pp];

//realative velocity between point_particle and fluid
real ur=sqrt((up-uf)*(up-uf)+(vp-vf)*(vp-vf)+(wp-wf)*(wp-wf));
real nu=mu/rho_f;

//Particle Reynolds number
real Rep=ur*dia/nu+EPSILON;


//Modification to stokes theory when Rep>1
real F = 1.0f+0.15f*powf(Rep,0.687f); 
//real F = 1.0f; 
real taud = rhod*dia*dia/(18.0f*mu);
real itau=F/taud;

//if(pp==1) printf("\nitau %f %f %f %f\n",Rep,F,taud,itau);

real volume=1./6. * PI * dia*dia*dia;

//Including soluble part:rhod*volume  and insoluble part: ms
real ms =  points[pp].ms;
real mp =  rhod *volume + ms;
//fluid mass in particle volume
real mf =  rho_f *volume;
//real gammar=mp/mf;

//drag force on particle
real drag_x=(uf-up)*itau*mp;
real drag_y=(vf-vp)*itau*mp;
real drag_z=(wf-wp)*itau*mp;

//if(pp==0) printf("\ntau %f %f %f\n",1.f/itau,taud,Rep);
//if(pp==0) printf("\nflow vel %f %f %f\n",wf,wp,drag_z);
//Fluid stress at the particle position

real stress_x=lpt_stress_u[pp];
real stress_y=lpt_stress_v[pp];
real stress_z=lpt_stress_w[pp];


//Total add mass,  -gradP.x~z are the body force on fluid apart from gravity
real dudt=(stress_x/rho_f-gradP.x);
real dvdt=(stress_y/rho_f-gradP.y);
real dwdt=(stress_z/rho_f-gradP.z);

real add_x=dudt*mf;
real add_y=dvdt*mf;
real add_z=dwdt*mf;

//Store the fluid force on particle including add mass, drag, fluid force and gravity
//default C_add=0.5; C_stress=1; C_drag=1;
real Fx=(C_add*add_x+C_stress*stress_x*volume+C_drag*drag_x)+(mp-mf)*g.x;
real Fy=(C_add*add_y+C_stress*stress_y*volume+C_drag*drag_y)+(mp-mf)*g.y;
real Fz=(C_add*add_z+C_stress*stress_z*volume+C_drag*drag_z)+(mp-mf)*g.z;


real rhs_x=Fx+ C_drag*up*itau*mp;
real rhs_y=Fy+ C_drag*vp*itau*mp;
real rhs_z=Fz+ C_drag*wp*itau*mp;

//if(pp==0) printf("\n Fz %f %f %f %f %f %f\n",add_z,stress_z*volume,drag_z,(mp-mf)*g.z,Fz,rhs_z);



//Exchange rate of soluble mass into scalar field dms/dt = pi *dp^2*hp*(rho_s-rho_{sat})   ref eq (20) in Oresta&&Prosperetti(2014)


//Based on hp*dp/D=2+0.6Re_p^0.5 Sc^{1/3}, ref eq (22) in Oresta&&Prosperetti(2014)
real Nu =2+0.6*sqrt(Rep)*powf(nu/DIFF,1.0/3.0);
points[pp].Nu=Nu;
//Nu=2;
real hp =Nu*DIFF/dia;

//sc_eq=log(K_{ow}), and rho_{sat}= \kappa *ms/vp, where \kappa=rho_f/rhod*1/K_{ow}
real kappa=rho_f/rhod/__exp10f(sc_eq);
//real kappa=rho_f/rhod/__expf(sc_eq);
real rho_sat=kappa*ms/volume;
real msdot= PI*dia*dia*hp*(scg[pp]- rho_sat);



//if(pp==1) printf("\nhp %f %f \n",dia,hp);
//if(pp==1) printf("\nsc_eq %f %f %f %f %f %f\n",kappa,ms,volume,scg[pp], rho_sat, msdot);
//if(msdot>0) printf("\nmsdot %d %f %f %f %f %f\n",pp,ms,volume,scg[pp], rho_sat, msdot);
//if(scg[pp]<0) printf("\nscg %d %f %f %f %f %f\n",pp,ms,volume,scg[pp], rho_sat, msdot);

     //Treat taup as very small and integrate directly!!
      real acc_m=points[pp].msdot+C_drag*mp*itau;
      real total_m=mp+C_add*mf; 
      real idt=acc_m/total_m;
      real ratio=exp(-dt*idt);
	
      real x1,y1,z1;
      real u1,v1,w1;

 if(fabs(acc_m)>EPSILON)
{
      real termVel_x=rhs_x/acc_m;
      real termVel_y=rhs_y/acc_m;
      real termVel_z=rhs_z/acc_m;
//if(isnan(acc_m)||isnan(total_m)) printf("\nacc_m %f %f %f %f %d\n",points[pp].msdot,hp,scg[pp],sc_eq,pp);

 
     // update linear velocities
       u1 = ratio*(up -termVel_x) + termVel_x;
       v1 = ratio*(vp -termVel_y) + termVel_y;
       w1 = ratio*(wp -termVel_z) + termVel_z;


//if(fabs(points[pp].u)>100||fabs(points[pp].v)>100||fabs(points[pp].w)>100) printf("\nrhs %f %f %f\n",rhs_x,total_m,acc_m);
//if(fabs(points[pp].u)>100||fabs(points[pp].v)>100||fabs(points[pp].w)>100) printf("\nu0~w0 %f %f %f %f\n",rhs_x/total_m,rhs_x/acc_m,points[pp].u0,dt);

    // update position 
       x1 = xp + (up-termVel_x)*(1-ratio)/idt+termVel_x*dt;
       y1 = yp + (vp-termVel_y)*(1-ratio)/idt+termVel_y*dt;
       z1 = zp + (wp-termVel_z)*(1-ratio)/idt+termVel_z*dt;

}

else
{
real acc_x=rhs_x/total_m;
real acc_y=rhs_y/total_m;
real acc_z=rhs_z/total_m;

 u1=acc_x*dt+up;
 v1=acc_y*dt+vp;
 w1=acc_z*dt+wp;

 x1=xp+(0.5f *acc_x*dt+up)*dt;
 y1=yp+(0.5f *acc_y*dt+vp)*dt;
 z1=zp+(0.5f *acc_z*dt+wp)*dt;

}

      points[pp].u = u1;
      points[pp].v = v1;
      points[pp].w = w1;

      points[pp].x = x1;
      points[pp].y = y1;
      points[pp].z = z1;
 
//Store the fluid Momentum force on particle including add mass, drag, fluid force and gravity added mass effect on particle-fluid force interaction  gravity is not implemented on fluid, but we need to add -mf*g to reaction force; F=d(mp*u)/dt -mp*g, according to eq 6 from dropDiffusionForceImpleBlue.pdf


//This is in fact momentum back reaction in particle sub-timestep, it will be accumulated during the sub-step, and then been added to fluid source force.
//Usually we use this one

points[pp].Fx=(rhs_x-mp*g.x)*dt-(C_add*(u1-up)*mf+C_drag*(x1-xp)*mp*itau);
points[pp].Fy=(rhs_y-mp*g.y)*dt-(C_add*(v1-vp)*mf+C_drag*(y1-yp)*mp*itau);
points[pp].Fz=(rhs_z-mp*g.z)*dt-(C_add*(w1-wp)*mf+C_drag*(z1-zp)*mp*itau);


//TODO periodic BC for particles, may need to change in future
periodic_grid_position(points[pp].x,points[pp].y,points[pp].z,dom);
//update old values
      points[pp].x0 = points[pp].x;
      points[pp].y0 = points[pp].y;
      points[pp].z0 = points[pp].z;

      points[pp].u0 = points[pp].u;
      points[pp].v0 = points[pp].v;
      points[pp].w0 = points[pp].w;

     // update soluble mass
      points[pp].ms = points[pp].ms0 + msdot * dt;
      points[pp].msdot = msdot;
  
     // update particle time step
      points[pp].dt =  dt;

      points[pp].ms0 = points[pp].ms;

}

*/

//Twoway coupling of point motion
//The relaxed condition would be that C_drag*dms/dt can't be 0 together!
__global__ void drag_move_points_twoway(point_struct *points,dom_struct *dom, int npoints,
real *ug,real *vg,real *wg,
real *lpt_stress_u,real *lpt_stress_v,real *lpt_stress_w,
real *lpt_dudt,real *lpt_dvdt,real *lpt_dwdt,
real *scg,
real rho_f,real mu, g_struct g,gradP_struct gradP,
real C_add,real C_stress,real C_drag,
real sc_eq,real DIFF,real dt)
//gradP serve as bodyforce for the time being
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;

  if(pp >= npoints) return;

real up=points[pp].u;
real vp=points[pp].v;
real wp=points[pp].w;

real xp=points[pp].x;
real yp=points[pp].y;
real zp=points[pp].z;

real dia=2*points[pp].r;
real rhod=points[pp].rho;//rhod is the particle density


//fluid velocity at the particle position
real uf=ug[pp];
real vf=vg[pp];
real wf=wg[pp];

//realative velocity between point_particle and fluid
real ur=sqrt((up-uf)*(up-uf)+(vp-vf)*(vp-vf)+(wp-wf)*(wp-wf));
real nu=mu/rho_f;

//Particle Reynolds number
real Rep=ur*dia/nu+EPSILON;


//Modification to stokes theory when Rep>1
real F = 1.0f+0.15f*powf(Rep,0.687f); 
//real F = 1.0f; 
real taud = rhod*dia*dia/(18.0f*mu);
real itau=F/taud;

//if(pp==1) printf("\nitau %f %f %f %f\n",Rep,F,taud,itau);

real volume=1./6. * PI * dia*dia*dia;

//Including soluble part:rhod*volume  and insoluble part: ms
real ms =  points[pp].ms;

//TODO change:Let ms don't influence the momentum for the time being
//real mp =  rhod *volume + ms;
real mp =  rhod *volume ;
//fluid mass in particle volume
real mf =  rho_f *volume;
//real gammar=mp/mf;

//drag force on particle
real drag_x=(uf-up)*itau*mp;
real drag_y=(vf-vp)*itau*mp;
real drag_z=(wf-wp)*itau*mp;

//if(pp==0) printf("\ntau %f %f %f\n",1.f/itau,taud,Rep);
//if(pp==0) printf("\nflow vel %f %f %f\n",wf,wp,drag_z);
//Fluid stress at the particle position
 

//Fluid acceleration evaluated at the particle position
real stress_x=lpt_stress_u[pp];
real stress_y=lpt_stress_v[pp];
real stress_z=lpt_stress_w[pp];

real dudt=lpt_dudt[pp];
real dvdt=lpt_dvdt[pp];
real dwdt=lpt_dwdt[pp];

/*
//Total add mass,  -gradP.x~z are the body force on fluid apart from gravity
real dudt=(stress_x/rho_f-gradP.x);
real dvdt=(stress_y/rho_f-gradP.y);
real dwdt=(stress_z/rho_f-gradP.z);

*/
 
real add_x=dudt*mf;
real add_y=dvdt*mf;
real add_z=dwdt*mf;

//Store the fluid force on particle including add mass, drag, fluid force and gravity
//default C_add=0.5; C_stress=1; C_drag=1;
real Fx=(C_add*add_x+C_stress*stress_x*volume+C_drag*drag_x)+(mp-mf)*g.x;
real Fy=(C_add*add_y+C_stress*stress_y*volume+C_drag*drag_y)+(mp-mf)*g.y;
real Fz=(C_add*add_z+C_stress*stress_z*volume+C_drag*drag_z)+(mp-mf)*g.z;


real rhs_x=Fx+ C_drag*up*itau*mp;
real rhs_y=Fy+ C_drag*vp*itau*mp;
real rhs_z=Fz+ C_drag*wp*itau*mp;


/*
Exchange rate of soluble mass into scalar field
dms/dt = pi *dp^2*hp*(rho_s-rho_{sat})   ref eq (20) in Oresta&&Prosperetti(2014)
*/
//Based on hp*dp/D=2+0.6Re_p^0.5 Sc^{1/3}, ref eq (22) in Oresta&&Prosperetti(2014)
real Nu =2+0.6*sqrt(Rep)*powf(nu/DIFF,1.0/3.0);
points[pp].Nu=Nu;
//Nu=2;
real hp =Nu*DIFF/dia;

//sc_eq=log(K_{ow}), and rho_{sat}= \kappa *ms/vp, where \kappa=rho_f/rhod*1/K_{ow}
real kappa=rho_f/rhod/__exp10f(sc_eq);
//real kappa=rho_f/rhod/__expf(sc_eq);
real rho_sat=kappa*ms/volume;
real msdot= PI*dia*dia*hp*(scg[pp]- rho_sat);
 
     //Treat taup as very small and integrate directly!!
      real acc_m=points[pp].msdot+C_drag*mp*itau;
      real total_m=mp+C_add*mf; 
      real idtp=acc_m/total_m;
      real dtp=1/idtp;
      real ratio=exp(-dt*idtp);
	
      real x1,y1,z1;
      real u1,v1,w1;

 if(fabs(acc_m)>EPSILON||C_drag>0)
{
      real termVel_x=rhs_x/acc_m;
      real termVel_y=rhs_y/acc_m;
      real termVel_z=rhs_z/acc_m;

     // update linear velocities
       u1 = ratio*(up -termVel_x) + termVel_x;
       v1 = ratio*(vp -termVel_y) + termVel_y;
       w1 = ratio*(wp -termVel_z) + termVel_z;

//if(pp==0) printf("\ntermVz %f %f %f %f %f %f\n",termVel_z,w1,wp,wf,Rep,ratio);

//if(fabs(points[pp].u)>100||fabs(points[pp].v)>100||fabs(points[pp].w)>100) printf("\nrhs %f %f %f\n",rhs_x,total_m,acc_m);
//if(fabs(points[pp].u)>100||fabs(points[pp].v)>100||fabs(points[pp].w)>100) printf("\nu0~w0 %f %f %f %f\n",rhs_x/total_m,rhs_x/acc_m,points[pp].u0,dt);

    // update position 
       x1 = xp + (up-termVel_x)*(1-ratio)*dtp+termVel_x*dt;
       y1 = yp + (vp-termVel_y)*(1-ratio)*dtp+termVel_y*dt;
       z1 = zp + (wp-termVel_z)*(1-ratio)*dtp+termVel_z*dt;

}

else
{
real acc_x=rhs_x/total_m;
real acc_y=rhs_y/total_m;
real acc_z=rhs_z/total_m;

 u1=acc_x*dt+up;
 v1=acc_y*dt+vp;
 w1=acc_z*dt+wp;

 x1=xp+(0.5f *acc_x*dt+up)*dt;
 y1=yp+(0.5f *acc_y*dt+vp)*dt;
 z1=zp+(0.5f *acc_z*dt+wp)*dt;

}

      real idt=1/dt;
      points[pp].udot = (u1-up)*idt;
      points[pp].vdot = (v1-vp)*idt;
      points[pp].wdot = (w1-wp)*idt;

      points[pp].u = u1;
      points[pp].v = v1;
      points[pp].w = w1;

      points[pp].x = x1;
      points[pp].y = y1;
      points[pp].z = z1;

/*
      points[pp].w = 100;
      points[pp].z = points[pp].z0+points[pp].w*dt;
*/

// printf("\npartVel %f %f %f\n",points[pp].z,points[pp].z0,points[pp].w);

/*
Store the fluid Momentum force on particle including add mass, drag, fluid force and gravity
added mass effect on particle-fluid force interaction
gravity is not implemented on fluid, but we need to add -mf*g to reaction force
F=d(mp*u)/dt -mp*g, according to eq 6 from dropDiffusionForceImpleBlue.pdf
*/
//This is in fact momentum back reaction in particle sub-timestep, it will be accumulated during the sub-step, and then been added to fluid source force.
//Usually we use this one
/*
points[pp].Fx=(rhs_x-mp*g.x)*dt-(C_add*(u1-up)*mf+C_drag*(x1-xp)*mp*itau);
points[pp].Fy=(rhs_y-mp*g.y)*dt-(C_add*(v1-vp)*mf+C_drag*(y1-yp)*mp*itau);
points[pp].Fz=(rhs_z-mp*g.z)*dt-(C_add*(w1-wp)*mf+C_drag*(z1-zp)*mp*itau);

points[pp].Fx=rhs_x*dt-(C_add*(u1-up)*mf+C_drag*(x1-xp)*mp*itau);
points[pp].Fy=rhs_y*dt-(C_add*(v1-vp)*mf+C_drag*(y1-yp)*mp*itau);
points[pp].Fz=rhs_z*dt-(C_add*(w1-wp)*mf+C_drag*(z1-zp)*mp*itau);

*/

points[pp].Fx=(rhs_x-(mp-mf)*g.x)*dt-(C_add*(u1-up)*mf+C_drag*(x1-xp)*mp*itau);
points[pp].Fy=(rhs_y-(mp-mf)*g.y)*dt-(C_add*(v1-vp)*mf+C_drag*(y1-yp)*mp*itau);
points[pp].Fz=(rhs_z-(mp-mf)*g.z)*dt-(C_add*(w1-wp)*mf+C_drag*(z1-zp)*mp*itau);

/*
points[pp].ox +=points[pp].Fx;
points[pp].oy +=points[pp].Fy;
points[pp].oz +=points[pp].Fz;
real dx=dom->dx;
real Vg=dx*dx*dx;
//if(pp==0) printf("\nFz on point %f %f %f %f %f %f\n",add_z,stress_z*volume,drag_z,(mp-mf)*g.z,Fz,rhs_z);
//if(pp==0) printf("\nFz on point %f %f %f %f\n",drag_z,(mp-mf)*g.z,Fz,rhs_z);
if(pp==0) printf("\nFz on point %f %f %f %f\n",drag_z,(mp-mf)*g.z*dt,points[pp].Fz,rhs_z);
if(pp==0) printf("\nreaction %f %f %f %f %f\n",(rhs_z-mp*g.z)*dt,rhs_z,C_add*(u1-up)*mf,C_drag*(z1-zp)*mp*itau,points[pp].Fz);
if(pp==0) printf("\nvelPredic %f %f %f\n",points[pp].Fz/Vg,points[pp].Fz,Vg);
*/

//if(pp==0) printf("\nposition %f %f %f %f %f\n",x1-xp,x1,z1-zp,z1,zp);
//if(pp==0) printf("\npart vel %f %f %f %f %f %f\n",w1-wp,w1,u1,u1-up,dt,rhs_z);

//TODO periodic BC for particles, may need to change in future
periodic_grid_position(points[pp].x,points[pp].y,points[pp].z,dom);
//update old values
      points[pp].x0 = points[pp].x;
      points[pp].y0 = points[pp].y;
      points[pp].z0 = points[pp].z;

      points[pp].u0 = points[pp].u;
      points[pp].v0 = points[pp].v;
      points[pp].w0 = points[pp].w;

     // update soluble mass
      points[pp].ms = points[pp].ms0 + msdot * dt;
      points[pp].msdot = msdot;
  
     // update particle time step
      points[pp].dt =  dt;

      points[pp].ms0 = points[pp].ms;

}




//Note for this method, C_drag has to be greater than 0!!
__global__ void drag_move_points_init(point_struct *points,dom_struct *dom, int npoints,
real *ug,real *vg,real *wg,
real *lpt_stress_u,real *lpt_stress_v,real *lpt_stress_w,
real rho_f,real mu, g_struct g,gradP_struct gradP,
real C_add,real C_stress,real C_drag,
real sc_eq,real DIFF)
//gradP serve as bodyforce for the time being
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;

  if(pp >= npoints) return;
/*
real up=points[pp].u;
real vp=points[pp].v;
real wp=points[pp].w;
real dia=2*points[pp].r;
real rhod=points[pp].rho;//rhod is the particle density

//particle interaction force
real iFx=points[pp].iFx;
real iFy=points[pp].iFy;
real iFz=points[pp].iFz;

//fluid velocity at the particle position
real uf=ug[pp];
real vf=vg[pp];
real wf=wg[pp];

//realative velocity between point_particle and fluid
real ur=sqrt((up-uf)*(up-uf)+(vp-vf)*(vp-vf)+(wp-wf)*(wp-wf));
real nu=mu/rho_f;

//Particle Reynolds number
real Rep=ur*dia/nu+EPSILON;

//Based on hp*dp/D=2+0.6Re_p^0.5 Sc^{1/3}, ref eq (22) in Oresta&&Prosperetti(2014)
real Nu =2+0.6*sqrt(Rep)*powf(nu/DIFF,1.0/3.0);
points[pp].Nu=Nu;

 
//Modification to stokes theory when Rep>1
real F = 1.0f+0.15f*powf(Rep,0.687f); 
//real F = 1.0f; 
real taud = rhod*dia*dia/(18.0f*mu);
real itau=F/taud;

//if(pp==1) printf("\nitau %f %f %f %f\n",Rep,F,taud,itau);

real volume=1./6. * PI * dia*dia*dia;

//Including soluble part:rhod*volume  and insoluble part: ms
real mp =  rhod *volume + points[pp].ms;
//fluid mass in particle volume
real mf =  rho_f *volume;
real msdot=points[pp].msdot;
real gammar=mp/mf;

//drag force on particle
real drag_x=(uf-up)*itau*mp;
real drag_y=(vf-vp)*itau*mp;
real drag_z=(wf-wp)*itau*mp;

//Fluid stress at the particle position
real stress_x=lpt_stress_u[pp];
real stress_y=lpt_stress_v[pp];
real stress_z=lpt_stress_w[pp];

//Total add mass,  -gradP.x~z are the body force on fluid apart from gravity
real add_x=(stress_x/rho_f-gradP.x)*mf;
real add_y=(stress_y/rho_f-gradP.y)*mf;
real add_z=(stress_z/rho_f-gradP.z)*mf;

//Store the fluid force on particle including add mass, drag, fluid force and gravity
//default C_add=0.5; C_stress=1; C_drag=1;
real Fx=(C_add*add_x+C_stress*stress_x*volume+C_drag*drag_x)+(mp-mf)*g.x;
real Fy=(C_add*add_y+C_stress*stress_y*volume+C_drag*drag_y)+(mp-mf)*g.y;
real Fz=(C_add*add_z+C_stress*stress_z*volume+C_drag*drag_z)+(mp-mf)*g.z;



//Store the temp particle acceleration, also including the particle soluble mass change  in the last term!
real udot =(Fx+ iFx -up*msdot)/mp;
real vdot =(Fy+ iFy -vp*msdot)/mp;
real wdot =(Fz+ iFz -wp*msdot)/mp;

 
//acount for added mass effect since it appears also on the left handside of particle governing equation
if(fabs(C_add-0.5f)<EPSILON)
{
      udot = udot/(1+C_add/gammar);
      vdot = vdot/(1+C_add/gammar);
      wdot = wdot/(1+C_add/gammar);
}


 
points[pp].Fx=Fx-C_add*udot*mf-mp*g.x;
points[pp].Fy=Fy-C_add*vdot*mf-mp*g.y;
points[pp].Fz=Fz-C_add*wdot*mf-mp*g.z;

points[pp].udot=udot;
points[pp].vdot=vdot;
points[pp].wdot=wdot;

 

//TODO periodic BC for particles, may need to change in future
periodic_grid_position(points[pp].x,points[pp].y,points[pp].z,dom);
*/

//update old values
      points[pp].x0 = points[pp].x;
      points[pp].y0 = points[pp].y;
      points[pp].z0 = points[pp].z;

      points[pp].u0 = points[pp].u;
      points[pp].v0 = points[pp].v;
      points[pp].w0 = points[pp].w;

      points[pp].ms0 = points[pp].ms;

}


__global__ void store_pointsD(point_struct *points,dom_struct *dom, int npoints)
//gradP serve as bodyforce for the time being
{
  int pp = threadIdx.x + blockIdx.x*blockDim.x;

  if(pp >= npoints) return;
 
      //update old values
      points[pp].x0 = points[pp].x;
      points[pp].y0 = points[pp].y;
      points[pp].z0 = points[pp].z;

      points[pp].u0 = points[pp].u;
      points[pp].v0 = points[pp].v;
      points[pp].w0 = points[pp].w;

      points[pp].ms0 = points[pp].ms;

}




__global__ void reset_flag_u(int *flag_u, dom_struct *dom, BC bc)
{
  int i;    // iterator
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  // flag everything as fluid
  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    for(i = dom->Gfx._isb; i < dom->Gfx._ieb; i++) {
      flag_u[i + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = 1;
    }
  }

  __syncthreads();

  // flag external boundaries
  if((tj < dom->Gfx._jnb) && (tk < dom->Gfx._knb)) {
    if(bc.uW != PERIODIC)
      flag_u[dom->Gfx._is + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = 0;
    if(bc.uE != PERIODIC)
      flag_u[dom->Gfx._ie-1 + tj*dom->Gfx._s1b + tk*dom->Gfx._s2b] = 0;
  }
}

__global__ void reset_flag_v(int *flag_v, dom_struct *dom, BC bc)
{
  int j;    // iterator
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  // flag everything as fluid
  if((tk < dom->Gfy._knb) && (ti < dom->Gfy._inb)) {
    for(j = dom->Gfy._jsb; j < dom->Gfy._jeb; j++) {
      flag_v[ti + j*dom->Gfy._s1b + tk*dom->Gfy._s2b] = 1;
    }
  }

  __syncthreads();

  // flag external boundaries
  if((tk < dom->Gfy._knb) && (ti < dom->Gfy._inb)) {
    if(bc.vS != PERIODIC)
      flag_v[ti + dom->Gfy._js*dom->Gfy._s1b + tk*dom->Gfy._s2b] = 0;
    if(bc.vN != PERIODIC)
      flag_v[ti + (dom->Gfy._je-1)*dom->Gfy._s1b + tk*dom->Gfy._s2b] = 0;
  }
}

__global__ void reset_flag_w(int *flag_w, dom_struct *dom, BC bc)
{
  int k;    // iterator
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  // flag everything as fluid
  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    for(k = dom->Gfz._ksb; k < dom->Gfz._keb; k++) {
      flag_w[ti + tj*dom->Gfz._s1b + k*dom->Gfz._s2b] = 1;
    }
  }

  __syncthreads();

  // flag external boundaries
  if((ti < dom->Gfz._inb) && (tj < dom->Gfz._jnb)) {
    if(bc.wB != PERIODIC)
      flag_w[ti + tj*dom->Gfz._s1b + dom->Gfz._ks*dom->Gfz._s2b] = 0;
    if(bc.wT != PERIODIC)
      flag_w[ti + tj*dom->Gfz._s1b + (dom->Gfz._ke-1)*dom->Gfz._s2b] = 0;
  }
}


__global__ void cage_flag_u_periodic_W(int *flag_u, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfx.jnb && tk < dom->Gfx.knb) {
    flag_u[dom->Gfx.isb + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b]
      = flag_u[(dom->Gfx.ie-2) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_u_periodic_E(int *flag_u, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfx.jnb && tk < dom->Gfx.knb) {
    flag_u[(dom->Gfx.ieb-1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b]
      = flag_u[(dom->Gfx.is+1) + tj*dom->Gfx.s1b + tk*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_u_periodic_S(int *flag_u, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfx.knb && ti < dom->Gfx.inb) {
    flag_u[ti + dom->Gfx.jsb*dom->Gfx.s1b + tk*dom->Gfx.s2b]
      = flag_u[ti + (dom->Gfx.je-1)*dom->Gfx.s1b + tk*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_u_periodic_N(int *flag_u, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfx.knb && ti < dom->Gfx.inb) {
    flag_u[ti + (dom->Gfx.jeb-1)*dom->Gfx.s1b + tk*dom->Gfx.s2b]
      = flag_u[ti + dom->Gfx.js*dom->Gfx.s1b + tk*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_u_periodic_B(int *flag_u, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfx.inb && tj < dom->Gfx.jnb) {
    flag_u[ti + tj*dom->Gfx.s1b + dom->Gfx.ksb*dom->Gfx.s2b]
      = flag_u[ti + tj*dom->Gfx.s1b + (dom->Gfx.ke-1)*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_u_periodic_T(int *flag_u, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfx.inb && tj < dom->Gfx.jnb) {
    flag_u[ti + tj*dom->Gfx.s1b + (dom->Gfx.keb-1)*dom->Gfx.s2b]
      = flag_u[ti + tj*dom->Gfx.s1b + dom->Gfx.ks*dom->Gfx.s2b];
  }
}

__global__ void cage_flag_v_periodic_W(int *flag_v, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfy.jnb && tk < dom->Gfy.knb) {
    flag_v[dom->Gfy.isb + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b]
      = flag_v[(dom->Gfy.ie-1) + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_v_periodic_E(int *flag_v, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfy.jnb && tk < dom->Gfy.knb) {
    flag_v[(dom->Gfy.ieb-1) + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b]
      = flag_v[dom->Gfy.is + tj*dom->Gfy.s1b + tk*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_v_periodic_S(int *flag_v, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfy.knb && ti < dom->Gfy.inb) {
    flag_v[ti + dom->Gfy.jsb*dom->Gfy.s1b + tk*dom->Gfy.s2b]
      = flag_v[ti + (dom->Gfy.je-2)*dom->Gfy.s1b + tk*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_v_periodic_N(int *flag_v, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfy.knb && ti < dom->Gfy.inb) {
    flag_v[ti + (dom->Gfy.jeb-1)*dom->Gfy.s1b + tk*dom->Gfy.s2b]
      = flag_v[ti + (dom->Gfy.js+1)*dom->Gfy.s1b + tk*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_v_periodic_B(int *flag_v, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfy.inb && tj < dom->Gfy.jnb) {
    flag_v[ti + tj*dom->Gfy.s1b + dom->Gfy.ksb*dom->Gfy.s2b]
      = flag_v[ti + tj*dom->Gfy.s1b + (dom->Gfy.ke-1)*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_v_periodic_T(int *flag_v, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfy.inb && tj < dom->Gfy.jnb) {
    flag_v[ti + tj*dom->Gfy.s1b + (dom->Gfy.keb-1)*dom->Gfy.s2b]
      = flag_v[ti + tj*dom->Gfy.s1b + dom->Gfy.ks*dom->Gfy.s2b];
  }
}

__global__ void cage_flag_w_periodic_W(int *flag_w, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfz.jnb && tk < dom->Gfz.knb) {
    flag_w[dom->Gfz.isb + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b]
      = flag_w[(dom->Gfz.ie-1)+ tj*dom->Gfz.s1b + tk*dom->Gfz.s2b];
  }
}

__global__ void cage_flag_w_periodic_E(int *flag_w, dom_struct *dom)
{
  int tj = blockDim.x*blockIdx.x + threadIdx.x;
  int tk = blockDim.y*blockIdx.y + threadIdx.y;

  if(tj < dom->Gfz.jnb && tk < dom->Gfz.knb) {
    flag_w[(dom->Gfz.ieb-1) + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b]
      = flag_w[dom->Gfz.is + tj*dom->Gfz.s1b + tk*dom->Gfz.s2b];
  }
}

__global__ void cage_flag_w_periodic_S(int *flag_w, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfz.knb && ti < dom->Gfz.inb) {
    flag_w[ti + dom->Gfz.jsb*dom->Gfz.s1b + tk*dom->Gfz.s2b]
      = flag_w[ti + (dom->Gfz.je-1)*dom->Gfz.s1b + tk*dom->Gfz.s2b];
  }
}

__global__ void cage_flag_w_periodic_N(int *flag_w, dom_struct *dom)
{
  int tk = blockDim.x*blockIdx.x + threadIdx.x;
  int ti = blockDim.y*blockIdx.y + threadIdx.y;

  if(tk < dom->Gfz.knb && ti < dom->Gfz.inb) {
    flag_w[ti + (dom->Gfz.jeb-1)*dom->Gfz.s1b + tk*dom->Gfz.s2b]
      = flag_w[ti + dom->Gfz.js*dom->Gfz.s1b + tk*dom->Gfz.s2b];
  }
}

__global__ void cage_flag_w_periodic_B(int *flag_w, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfz.inb && tj < dom->Gfz.jnb) {
    flag_w[ti + tj*dom->Gfz.s1b + dom->Gfz.ksb*dom->Gfz.s2b]
      = flag_w[ti + tj*dom->Gfz.s1b + (dom->Gfz.ke-2)*dom->Gfz.s2b];
  }
}

__global__ void cage_flag_w_periodic_T(int *flag_w, dom_struct *dom)
{
  int ti = blockDim.x*blockIdx.x + threadIdx.x;
  int tj = blockDim.y*blockIdx.y + threadIdx.y;

  if(ti < dom->Gfz.inb && tj < dom->Gfz.jnb) {
    flag_w[ti + tj*dom->Gfz.s1b + (dom->Gfz.keb-1)*dom->Gfz.s2b]
      = flag_w[ti + tj*dom->Gfz.s1b + (dom->Gfz.ks+1)*dom->Gfz.s2b];
  }
}





