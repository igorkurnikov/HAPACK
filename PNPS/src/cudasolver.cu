#ifndef _CUDA_SOLVER_CU_
#define _CUDA_SOLVER_CU_

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>

#include <cuda.h>

#include "pnpstructs.h"
//#include "pnpdebug.h"
#define DefClock0 clock_t time0;timeval tvtime0,tvtime1;
#define StartClock0 time0=clock ();gettimeofday(&tvtime0,NULL);
#define StopClock0 gettimeofday(&tvtime1,NULL);DbgPrint0("Time : %g s(CPU Time) %g  s(Wall Time)\n",((double)(clock ()-time0))/CLOCKS_PER_SEC,double(tvtime1.tv_sec)+(double(tvtime1.tv_usec)/1000000.0)-double(tvtime0.tv_sec)-(double(tvtime0.tv_usec)/1000000.0));
#define StopClockWMes0(Massege) gettimeofday(&tvtime1,NULL);printf("Time for %s is %.5g s (CPU Time) %g  s(Wall Time)\n",(Massege),((double)(clock ()-time0))/CLOCKS_PER_SEC,double(tvtime1.tv_sec)+(double(tvtime1.tv_usec)/1000000.0)-double(tvtime0.tv_sec)-(double(tvtime0.tv_usec)/1000000.0));

#define FMUL __fmul_rn
#define FADD __fadd_rn
#define FMAF __fmaf_rn

__constant__ float  dc_om1;
__constant__ float  dc_om2d6;
__constant__ float* dc_P[8];
__constant__ int dc_Qnum[8];
__constant__ float* dc_Q[8];
__constant__ float* dc_Qmult[8];
__constant__ int* dc_Qpos[8];

__constant__ int dc_DBnum[8];
__constant__ float* dc_DielMult[48];
__constant__ int* dc_DBpos[8];
#define dc_lookupVirGridSIZE 1024
__constant__ int  dc_lookupVirGrid[dc_lookupVirGridSIZE];


__global__ void KFLaplaceB(float* d_P0,float* d_P1,float* d_P2,float* d_P3,
																float* d_P4,float* d_P5,float* d_P6,float* d_P7,
																const float om1,const float om2d6,
								const int pitchX,const int pitchXY,const int pitchXY_BS_Z,const int BS_X,const int BS_XY)
{
// 	int VirXblock = (blockIdx.x%VirGridX);
// 	int VirYblock = (blockIdx.x/VirGridX);
// 	int VirZblock = blockIdx.y;
// 	
// 	int tx = VirXblock*BS_X + threadIdx.x;
// 	int ty = VirYblock*BS_Y + threadIdx.y+1;
// 	int tz = VirZblock*BS_Z + threadIdx.z+1;
// 	
// 	int i=tx+ty*pitchX+tz*pitchXY;
	int t=threadIdx.x+threadIdx.y*pitchX+threadIdx.z*pitchXY;
	int i=dc_lookupVirGrid[blockIdx.x] + blockIdx.y*pitchXY_BS_Z + t;
	t=threadIdx.x+threadIdx.y*BS_X+threadIdx.z*BS_XY;
	
	float xP0,xP3,xP5,xP6;
	float yP0,yP3,yP5,yP6;
	float zP0,zP3,zP5,zP6;
	
	__shared__ float shP[256];
	
	
	//do over P1
	shP[t]=d_P1[i];
	__syncthreads();
	
	zP5=FADD(shP[t],d_P1[i+pitchXY]);
	yP3=FADD(shP[t],d_P1[i+pitchX]);
	if(threadIdx.x!=0)
		xP0=FADD(shP[t],shP[t-1]);
	else
		xP0=FADD(shP[t],d_P1[i-1]);
	__syncthreads();
	//do over P2
	shP[t]=d_P2[i];
	__syncthreads();
	
	zP6=FADD(shP[t],d_P2[i+pitchXY]);
	yP0=FADD(shP[t],d_P2[i-pitchX]);
	if(threadIdx.x!=blockDim.x-1)
		xP3=FADD(shP[t],shP[t+1]);
	else
		xP3=FADD(shP[t],d_P2[i+1]);
	__syncthreads();
	//do over P4
	shP[t]=d_P4[i];
	__syncthreads();
	zP0=FADD(shP[t],d_P4[i-pitchXY]);
	yP6=FADD(shP[t],d_P4[i+pitchX]);
	if(threadIdx.x!=blockDim.x-1)
		xP5=FADD(shP[t],shP[t+1]);
	else
		xP5=FADD(shP[t],d_P4[i+1]);
	__syncthreads();
	//do over P7
	shP[t]=d_P7[i];
	__syncthreads();
	zP3=FADD(shP[t],d_P7[i-pitchXY]);
	yP5=FADD(shP[t],d_P7[i-pitchX]);
	if(threadIdx.x!=0)
		xP6=FADD(shP[t],shP[t-1]);
	else
		xP6=FADD(shP[t],d_P7[i-1]);
	
	
	d_P0[i]=FADD(FMUL(om1,d_P0[i]),FMUL(om2d6,FADD(FADD(xP0,yP0),zP0)));
	d_P3[i]=FADD(FMUL(om1,d_P3[i]),FMUL(om2d6,FADD(FADD(xP3,yP3),zP3)));
	d_P5[i]=FADD(FMUL(om1,d_P5[i]),FMUL(om2d6,FADD(FADD(xP5,yP5),zP5)));
	d_P6[i]=FADD(FMUL(om1,d_P6[i]),FMUL(om2d6,FADD(FADD(xP6,yP6),zP6)));
}
__global__ void KFLaplaceW(float* d_P0,float* d_P1,float* d_P2,float* d_P3,
																float* d_P4,float* d_P5,float* d_P6,float* d_P7,
								const float om1,const float om2d6,
				const int pitchX,const int pitchXY,const int pitchXY_BS_Z,const int BS_X,const int BS_XY)
{
// 	int VirXblock = (blockIdx.x%VirGridX);
// 	int VirYblock = (blockIdx.x/VirGridX);
// 	int VirZblock = blockIdx.y;
// 	
// 	int tx = VirXblock*BS_X + threadIdx.x;
// 	int ty = VirYblock*BS_Y + threadIdx.y+1;
// 	int tz = VirZblock*BS_Z + threadIdx.z+1;
// 	
// 	int i=tx+ty*pitchX+tz*pitchXY;
	//int i=dc_lookupVirGrid[blockIdx.x] + blockIdx.y*pitchXY_BS_Z + threadIdx.x+threadIdx.y*pitchX+threadIdx.z*pitchXY;
	
	int t=threadIdx.x+threadIdx.y*pitchX+threadIdx.z*pitchXY;
	int i=dc_lookupVirGrid[blockIdx.x] + blockIdx.y*pitchXY_BS_Z + t;
	t=threadIdx.x+threadIdx.y*BS_X+threadIdx.z*BS_XY;
	
	float xP1,xP2,xP4,xP7;
	float yP1,yP2,yP4,yP7;
	float zP1,zP2,zP4,zP7;
	
	__shared__ float shP[256];
	
	//do dc_P[0]
	shP[t]=d_P0[i];
	__syncthreads();
	zP4=FADD(shP[t],d_P0[i+pitchXY]);
	yP2=FADD(shP[t],d_P0[i+pitchX]);
	if(threadIdx.x!=blockDim.x-1)
		xP1=FADD(shP[t],shP[t+1]);
	else
		xP1=FADD(shP[t],d_P0[i+1]);
	__syncthreads();
	//do d_P[3]
	shP[t]=d_P3[i];
	__syncthreads();
	zP7=FADD(shP[t],d_P3[i+pitchXY]);
	yP1=FADD(shP[t],d_P3[i-pitchX]);
	if(threadIdx.x!=0)
		xP2=FADD(shP[t],shP[t-1]);
	else
		xP2=FADD(shP[t],d_P3[i-1]);
	__syncthreads();
	//do d_P[5]
	shP[t]=d_P5[i];
	__syncthreads();
	zP1=FADD(shP[t],d_P5[i-pitchXY]);
	yP7=FADD(shP[t],d_P5[i+pitchX]);
	if(threadIdx.x!=0)
		xP4=FADD(shP[t],shP[t-1]);
	else
		xP4=FADD(shP[t],d_P5[i-1]);
	__syncthreads();
	//do d_P6
	shP[t]=d_P6[i];
	__syncthreads();
	zP2=FADD(shP[t],d_P6[i-pitchXY]);
	yP4=FADD(shP[t],d_P6[i-pitchX]);
	if(threadIdx.x!=blockDim.x-1)
		xP7=FADD(shP[t],shP[t+1]);
	else
		xP7=FADD(shP[t],d_P6[i+1]);
	
	//d_P6[i]=FADD(FMUL(om1,d_P6[i]),FMUL(om2d6,FADD(FADD(xP6,yP6),zP6)));
	
	d_P1[i]=FADD(FMUL(om1,d_P1[i]),FMUL(om2d6,FADD(FADD(xP1,yP1),zP1)));
	d_P2[i]=FADD(FMUL(om1,d_P2[i]),FMUL(om2d6,FADD(FADD(xP2,yP2),zP2)));
	d_P4[i]=FADD(FMUL(om1,d_P4[i]),FMUL(om2d6,FADD(FADD(xP4,yP4),zP4)));
	d_P7[i]=FADD(FMUL(om1,d_P7[i]),FMUL(om2d6,FADD(FADD(xP7,yP7),zP7)));
}
__global__ void KFPoissonQB()
{
	int t=threadIdx.x+blockIdx.x*blockDim.x;
	int i;
	if(t<dc_Qnum[0])
	{
		i=dc_Qpos[0][t];
		dc_P[0][i]=FADD(dc_P[0][i],dc_Q[0][t]);
	}
	if(t<dc_Qnum[3])
	{
		i=dc_Qpos[3][t];
		dc_P[3][i]=FADD(dc_P[3][i],dc_Q[3][t]);
	}
	if(t<dc_Qnum[5])
	{
		i=dc_Qpos[5][t];
		dc_P[5][i]=FADD(dc_P[5][i],dc_Q[5][t]);
	}
	if(t<dc_Qnum[6])
	{
		i=dc_Qpos[6][t];
		dc_P[6][i]=FADD(dc_P[6][i],dc_Q[6][t]);
	}
}
__global__ void KFPoissonQW()
{
	int t=threadIdx.x+blockIdx.x*blockDim.x;
	int i;
	if(t<dc_Qnum[1])
	{
		i=dc_Qpos[1][t];
		dc_P[1][i]=FADD(dc_P[1][i],dc_Q[1][t]);
	}
	if(t<dc_Qnum[2])
	{
		i=dc_Qpos[2][t];
		dc_P[2][i]=FADD(dc_P[2][i],dc_Q[2][t]);
	}
	if(t<dc_Qnum[4])
	{
		i=dc_Qpos[4][t];
		dc_P[4][i]=FADD(dc_P[4][i],dc_Q[4][t]);
	}
	if(t<dc_Qnum[7])
	{
		i=dc_Qpos[7][t];
		dc_P[7][i]=FADD(dc_P[7][i],dc_Q[7][t]);
	}
}
__global__ void KFPoissonDBB(float* d_P0,float* d_P1,float* d_P2,float* d_P3,
																float* d_P4,float* d_P5,float* d_P6,float* d_P7,
								const float om2d6,
				const int pitchX,const int pitchXY)
{
	int t=threadIdx.x+blockIdx.x*blockDim.x;
	int i;
	float xP,xM,yP,yM,zP,zM;
	if(t<dc_DBnum[0])
	{
		i=dc_DBpos[0][t];
		//P0 x
		xP=FMUL(dc_DielMult[PlusX][t],d_P1[i]);
		xM=FMUL(dc_DielMult[MinusX][t],d_P1[i-1]);
		//P0 y 
		yP=FMUL(dc_DielMult[PlusY][t],d_P2[i]);
		yM=FMUL(dc_DielMult[MinusY][t],d_P2[i-pitchX]);
		//P0 z
		zP=FMUL(dc_DielMult[PlusZ][t],d_P4[i]);
		zM=FMUL(dc_DielMult[MinusZ][t],d_P4[i-pitchXY]);
		
		xP=FADD(xM,xP);
		yP=FADD(yM,yP);
		zP=FADD(zM,zP);
		dc_P[0][i]=FADD(FADD(xP,yP),FADD(zP,dc_P[0][i]));
	}
	if(t<dc_DBnum[3])
	{
		i=dc_DBpos[3][t];
		//P3 x
		xM=FMUL(dc_DielMult[3*6+MinusX][t],d_P2[i]);
		xP=FMUL(dc_DielMult[3*6+PlusX][t],d_P2[i+1]);
		//P3 y
		yM=FMUL(dc_DielMult[3*6+MinusY][t],d_P1[i]);
		yP=FMUL(dc_DielMult[3*6+PlusY][t],d_P1[i+pitchX]);
		//P3 z
		zP=FMUL(dc_DielMult[3*6+PlusZ][t],d_P7[i]);
		zM=FMUL(dc_DielMult[3*6+MinusZ][t],d_P7[i-pitchXY]);
		
		xP=FADD(xM,xP);
		yP=FADD(yM,yP);
		zP=FADD(zM,zP);
		dc_P[3][i]=FADD(FADD(xP,yP),FADD(zP,dc_P[3][i]));
	}
	if(t<dc_DBnum[5])
	{
		i=dc_DBpos[5][t];
		//P5 x
		xM=FMUL(dc_DielMult[5*6+MinusX][t],d_P4[i]);
		xP=FMUL(dc_DielMult[5*6+PlusX][t],d_P4[i+1]);
		//P5 y
		yP=FMUL(dc_DielMult[5*6+PlusY][t],d_P7[i]);
		yM=FMUL(dc_DielMult[5*6+MinusY][t],d_P7[i-pitchX]);
		//P5 z
		zM=FMUL(dc_DielMult[5*6+MinusZ][t],d_P1[i]);
		zP=FMUL(dc_DielMult[5*6+PlusZ][t],d_P1[i+pitchXY]);
		
		xP=FADD(xM,xP);
		yP=FADD(yM,yP);
		zP=FADD(zM,zP);
		dc_P[5][i]=FADD(FADD(xP,yP),FADD(zP,dc_P[5][i]));
	}
	if(t<dc_DBnum[6])
	{
		i=dc_DBpos[6][t];
			//P6 x
		xP=FMUL(dc_DielMult[6*6+PlusX][t],d_P7[i]);
		xM=FMUL(dc_DielMult[6*6+MinusX][t],d_P7[i-1]);
		//P6 y
		yM=FMUL(dc_DielMult[6*6+MinusY][t],d_P4[i]);
		yP=FMUL(dc_DielMult[6*6+PlusY][t],d_P4[i+pitchX]);
		//P6 z
		zM=FMUL(dc_DielMult[6*6+MinusZ][t],d_P2[i]);
		zP=FMUL(dc_DielMult[6*6+PlusZ][t],d_P2[i+pitchXY]);
		
		xP=FADD(xM,xP);
		yP=FADD(yM,yP);
		zP=FADD(zM,zP);
		dc_P[6][i]=FADD(FADD(xP,yP),FADD(zP,dc_P[6][i]));
	}
}
__global__ void KFPoissonDBW(float* d_P0,float* d_P1,float* d_P2,float* d_P3,
																float* d_P4,float* d_P5,float* d_P6,float* d_P7,
								const float om2d6,
				const int pitchX,const int pitchXY)
{
	int t=threadIdx.x+blockIdx.x*blockDim.x;
	int i;
	float xP,xM,yP,yM,zP,zM;
	
	if(t<dc_DBnum[1])
	{
		i=dc_DBpos[1][t];
		//P1 x
		xM=FMUL(dc_DielMult[1*6+MinusX][t],d_P0[i]);
		xP=FMUL(dc_DielMult[1*6+PlusX][t],d_P0[i+1]);
		//P1 y
		yP=FMUL(dc_DielMult[1*6+PlusY][t],d_P3[i]);
		yM=FMUL(dc_DielMult[1*6+MinusY][t],d_P3[i-pitchX]);
		//P1 z
		zP=FMUL(dc_DielMult[1*6+PlusZ][t],d_P5[i]);
		zM=FMUL(dc_DielMult[1*6+MinusZ][t],d_P5[i-pitchXY]);
		
		xP=FADD(xM,xP);
		yP=FADD(yM,yP);
		zP=FADD(zM,zP);
		dc_P[1][i]=FADD(FADD(xP,yP),FADD(zP,dc_P[1][i]));
	}
	if(t<dc_DBnum[2])
	{
		i=dc_DBpos[2][t];
		//P2 x
		xP=FMUL(dc_DielMult[2*6+PlusX][t],d_P3[i]);
		xM=FMUL(dc_DielMult[2*6+MinusX][t],d_P3[i-1]);
		//P2 y
		yM=FMUL(dc_DielMult[2*6+MinusY][t],d_P0[i]);
		yP=FMUL(dc_DielMult[2*6+PlusY][t],d_P0[i+pitchX]);
		//P2 z
		zP=FMUL(dc_DielMult[2*6+PlusZ][t],d_P6[i]);
		zM=FMUL(dc_DielMult[2*6+MinusZ][t],d_P6[i-pitchXY]);
		
		xP=FADD(xM,xP);
		yP=FADD(yM,yP);
		zP=FADD(zM,zP);
		dc_P[2][i]=FADD(FADD(xP,yP),FADD(zP,dc_P[2][i]));
	}
	if(t<dc_DBnum[4])
	{
		i=dc_DBpos[4][t];
		//P4x
		xP=FMUL(dc_DielMult[4*6+PlusX][t],d_P5[i]);
		xM=FMUL(dc_DielMult[4*6+MinusX][t],d_P5[i-1]);
		//P4 y
		yP=FMUL(dc_DielMult[4*6+PlusY][t],d_P6[i]);
		yM=FMUL(dc_DielMult[4*6+MinusY][t],d_P6[i-pitchX]);
		//P4 z
		zM=FMUL(dc_DielMult[4*6+MinusZ][t],d_P0[i]);
		zP=FMUL(dc_DielMult[4*6+PlusZ][t],d_P0[i+pitchXY]);
		
		xP=FADD(xM,xP);
		yP=FADD(yM,yP);
		zP=FADD(zM,zP);
		dc_P[4][i]=FADD(FADD(xP,yP),FADD(zP,dc_P[4][i]));
	}
	if(t<dc_DBnum[7])
	{
		i=dc_DBpos[7][t];
		//P7 x
		xM=FMUL(dc_DielMult[7*6+MinusX][t],d_P6[i]);
		xP=FMUL(dc_DielMult[7*6+PlusX][t],d_P6[i+1]);
		//P7 y
		yM=FMUL(dc_DielMult[7*6+MinusY][t],d_P5[i]);
		yP=FMUL(dc_DielMult[7*6+PlusY][t],d_P5[i+pitchX]);
		//P7 z
		zM=FMUL(dc_DielMult[7*6+MinusZ][t],d_P3[i]);
		zP=FMUL(dc_DielMult[7*6+PlusZ][t],d_P3[i+pitchXY]);
		
		xP=FADD(xM,xP);
		yP=FADD(yM,yP);
		zP=FADD(zM,zP);
		dc_P[7][i]=FADD(FADD(xP,yP),FADD(zP,dc_P[7][i]));
	}
}

__global__ void KFPoissonQnCalcEB()
{
	int t=threadIdx.x+blockIdx.x*blockDim.x;
	int i;
	if(t<dc_Qnum[0])
	{
		i=dc_Qpos[0][t];
		dc_P[0][i]=FADD(dc_P[0][i],dc_Q[0][t]);
	}
	if(t<dc_Qnum[3])
	{
		i=dc_Qpos[3][t];
		dc_P[3][i]=FADD(dc_P[3][i],dc_Q[3][t]);
	}
	if(t<dc_Qnum[5])
	{
		i=dc_Qpos[5][t];
		dc_P[5][i]=FADD(dc_P[5][i],dc_Q[5][t]);
	}
	if(t<dc_Qnum[6])
	{
		i=dc_Qpos[6][t];
		dc_P[6][i]=FADD(dc_P[6][i],dc_Q[6][t]);
	}
}
__global__ void KFPoissonQnCalcW()
{
	int t=threadIdx.x+blockIdx.x*blockDim.x;
	int i;
	if(t<dc_Qnum[1])
	{
		i=dc_Qpos[1][t];
		dc_P[1][i]=FADD(dc_P[1][i],dc_Q[1][t]);
	}
	if(t<dc_Qnum[2])
	{
		i=dc_Qpos[2][t];
		dc_P[2][i]=FADD(dc_P[2][i],dc_Q[2][t]);
	}
	if(t<dc_Qnum[4])
	{
		i=dc_Qpos[4][t];
		dc_P[4][i]=FADD(dc_P[4][i],dc_Q[4][t]);
	}
	if(t<dc_Qnum[7])
	{
		i=dc_Qpos[7][t];
		dc_P[7][i]=FADD(dc_P[7][i],dc_Q[7][t]);
	}
}
int checkCUDAError(const char* msg);
int GetCUDADevStat();
extern "C" int DoPoissonSolverOnCudaFloat(PoissonSolverOnCudaStruct* PS, PoissonSolverOnCudaParamStruct CudaParm)
{
	printf("E1\n");
	GetCUDADevStat();
	DefClock0;
	int i,k;
	cudaEvent_t start, stop;
	float time;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	// Allocate vectors in device memory
	int ErrorCount=0;
	
	int GS_X=PS->GS[0];
	int GS_Y=PS->GS[1];
	int GS_Z=PS->GS[2];
	int GS_XY=GS_X*GS_Y;
	int GS_XYZ=GS_X*GS_Y*GS_Z;
	
	float om1 = 1.0-PS->Relaxation;
	float om2d6 = PS->Relaxation/6.0;
	cudaMemcpyToSymbol(dc_om1, &om1, sizeof(float), 0, cudaMemcpyHostToDevice );
	cudaMemcpyToSymbol(dc_om2d6, &om2d6, sizeof(float), 0, cudaMemcpyHostToDevice );
	
	printf("GS=[%d,%d,%d]=%d\n",GS_X,GS_Y,GS_Z,GS_XYZ);
	
	//Start Clock for GPU
	StartClock0;
	
	int BS_X=CudaParm.BS_X;
	int BS_Y=CudaParm.BS_Y;
	int BS_Z=CudaParm.BS_Z;
	int BS_XY=BS_X*BS_Y;
	int BS_XYZ=BS_X*BS_Y*BS_Z;
	
	dim3 dimBlock(BS_X,BS_Y,BS_Z);
	dim3 dimGridVirt(GS_X/BS_X/2, GS_Y/BS_Y/2, GS_Z/BS_Z/2);
	//d_P* is store in pitched array x has 16 for CUDA and Y/Z is +2 for BC
	dim3 spltGSWBC(PS->spltGSWBC[0],PS->spltGSWBC[1],PS->spltGSWBC[2]);
	dim3 dimGrid(dimGridVirt.x*dimGridVirt.y, dimGridVirt.z,1);
	
	
	int pitchX=spltGSWBC.x;
	int pitchXY=spltGSWBC.x*spltGSWBC.y;
	
	printf("dimBlock [%d,%d,%d]\n",dimBlock.x,dimBlock.y,dimBlock.z);
	printf("dimGrid [%d,%d,%d]\n",dimGrid.x,dimGrid.y,dimGrid.z);
	printf("dimGridVirt [%d,%d,%d]\n",dimGridVirt.x,dimGridVirt.y,dimGridVirt.z);
	printf("spltGSWBC [%d,%d,%d]\n",spltGSWBC.x,spltGSWBC.y,spltGSWBC.z);
	
	printf("Total number of threads %d\n",dimGrid.x*dimGrid.y*dimGrid.z*dimBlock.x*dimBlock.y*dimBlock.z);
	
	//fill dc_lookupVirGrid
	int  h_lookupVirGrid[dc_lookupVirGridSIZE];
	if(dc_lookupVirGridSIZE<dimGrid.x)
	{
		printf("ERROR: dc_lookupVirGridSIZE is smaller then dimGrid.x, make it at least %d\n",dimGrid.x);
		return 1;
	}
	int VirXblock,VirYblock;
	for(i=0;i<dimGrid.x;i++)
	{
		VirXblock = (i%dimGridVirt.x);
		VirYblock = (i/dimGridVirt.x);
		h_lookupVirGrid[i]=VirXblock*BS_X+(VirYblock*BS_Y + 1)*pitchX+pitchXY;
	}
	cudaMemcpyToSymbol(dc_lookupVirGrid, h_lookupVirGrid, dimGrid.x*sizeof(int), 0, cudaMemcpyHostToDevice );
	
	
	int GS_XYZsplit = spltGSWBC.x*spltGSWBC.y*spltGSWBC.z;
	int sizeGS_XYZsplit = GS_XYZsplit*sizeof(float);
	
	//allocate and copy to device lin-array
	float* d_P[8];
	for(i=0;i<8;i++)
		cudaMalloc((void**)&d_P[i], sizeGS_XYZsplit);
	for(i=0;i<8;i++)
		cudaMemcpy(d_P[i], PS->P[i], sizeGS_XYZsplit, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(dc_P, d_P, 8*sizeof(float*), 0, cudaMemcpyHostToDevice );
	//charges
	float* d_Q[8];
	int* d_Qpos[8];
	float* d_Qmult[8];
	int Qmax=0;
	for(i=0;i<8;i++)
	{
		d_Q[i]=NULL;
		d_Qpos[i]=NULL;
		d_Qmult[i]=NULL;
		if(PS->Qnum[i]>Qmax)Qmax=PS->Qnum[i];
		if(PS->Qnum[i]>0)
		{
			cudaMalloc((void**)&d_Q[i], PS->Qnum[i]*sizeof(float));
			cudaMalloc((void**)&d_Qpos[i], PS->Qnum[i]*sizeof(int));
			cudaMalloc((void**)&d_Qmult[i], PS->Qnum[i]*sizeof(int));
			
			cudaMemcpy(d_Q[i], PS->Q[i], PS->Qnum[i]*sizeof(float), cudaMemcpyHostToDevice);
			cudaMemcpy(d_Qpos[i], PS->Qpos[i], PS->Qnum[i]*sizeof(int), cudaMemcpyHostToDevice);
			cudaMemcpy(d_Qmult[i], PS->Qmult[i], PS->Qnum[i]*sizeof(int), cudaMemcpyHostToDevice);
		}
	}
	cudaMemcpyToSymbol(dc_Q, d_Q, 8*sizeof(float*), 0, cudaMemcpyHostToDevice );
	cudaMemcpyToSymbol(dc_Qnum, PS->Qnum, 8*sizeof(int), 0, cudaMemcpyHostToDevice );
	cudaMemcpyToSymbol(dc_Qpos, d_Qpos, 8*sizeof(int*), 0, cudaMemcpyHostToDevice );
	cudaMemcpyToSymbol(dc_Qmult, d_Qmult, 8*sizeof(int*), 0, cudaMemcpyHostToDevice );
	int Qblock=CudaParm.Qblock;
	int QGrid=Qmax/Qblock;
	if(Qmax%Qblock!=0)QGrid++;
	printf("Qmax=%d Qblock=%d QGrid=%d\n",Qmax,Qblock,QGrid);
	//Diel Border
	int d_DBNum[8];
	float* d_DielMult[48];
	int* d_DBPos[8];
	int DBmax=0;
	for(i=0;i<8;i++)
	{
		d_DBPos[i]=NULL;
		for(k=0;k<6;k++)
			d_DielMult[i*6+k]=NULL;
		if(PS->DielBordNum[i]>DBmax)DBmax=PS->DielBordNum[i];
		if(PS->DielBordNum[i]>0)
		{
			cudaMalloc((void**)&d_DBPos[i], PS->DielBordNum[i]*sizeof(int));
			cudaMemcpy(d_DBPos[i], PS->DielBordPos[i], PS->DielBordNum[i]*sizeof(int), cudaMemcpyHostToDevice);
			for(k=0;k<6;k++)
			{
				cudaMalloc((void**)&d_DielMult[i*6+k], PS->DielBordNum[i]*sizeof(float));
				cudaMemcpy(d_DielMult[i*6+k], PS->DielMult[i][k], PS->DielBordNum[i]*sizeof(float), cudaMemcpyHostToDevice);
			}
		}
	}
	cudaMemcpyToSymbol(dc_DBnum, PS->DielBordNum, 8*sizeof(int), 0, cudaMemcpyHostToDevice );
	cudaMemcpyToSymbol(dc_DBpos, d_DBPos, 8*sizeof(int*), 0, cudaMemcpyHostToDevice );
	cudaMemcpyToSymbol(dc_DielMult, d_DielMult, 48*sizeof(float*), 0, cudaMemcpyHostToDevice );
	
	
	int DBblock=CudaParm.DBblock;
	int DBGrid=DBmax/DBblock;
	if(DBmax%DBblock!=0)DBGrid++;
	printf("DBmax=%d DBblock=%d DBGrid=%d\n",DBmax,DBblock,DBGrid);
	//do loop
	cudaEventRecord( start, 0 );
	int j;
	double totalEnergy,dtmp1;
	double fpoh=4.0*M_PI*PS->GridScale;
	int AvrTECount=0;
	int CollectingTEforAvr=0;
	GetCUDADevStat();
	for(int iteration=1;iteration<=PS->MaxIterations;iteration++)
	{//pitchXY*BS_Z
		KFLaplaceB<<<dimGrid, dimBlock>>>(d_P[0],d_P[1],d_P[2],d_P[3], d_P[4],d_P[5],d_P[6],d_P[7], om1,om2d6, pitchX,pitchXY,pitchXY*BS_Z,BS_X,BS_XY);
		cudaThreadSynchronize();
		ErrorCount+=1-checkCUDAError("cuda kernel running: KFLaplaceB");
		
		if(QGrid>0)
		{
			KFPoissonQB<<<QGrid, Qblock>>>();
			cudaThreadSynchronize();
			ErrorCount+=1-checkCUDAError("cuda kernel running: KFPoissonQB");
		}
		if(DBGrid>0)
		{
			KFPoissonDBB<<<DBGrid, DBblock>>>(d_P[0],d_P[1],d_P[2],d_P[3], d_P[4],d_P[5],d_P[6],d_P[7], om2d6, pitchX,pitchXY);
			cudaThreadSynchronize();
			ErrorCount+=1-checkCUDAError("cuda kernel running: KFPoissonDBB");
		}
		
		KFLaplaceW<<<dimGrid, dimBlock>>>(d_P[0],d_P[1],d_P[2],d_P[3], d_P[4],d_P[5],d_P[6],d_P[7], om1,om2d6, pitchX,pitchXY,pitchXY*BS_Z,BS_X,BS_XY);
		cudaThreadSynchronize();
		ErrorCount+=1-checkCUDAError("cuda kernel running: KFLaplaceW");
		if(QGrid>0)
		{
			KFPoissonQW<<<QGrid, Qblock>>>();
			cudaThreadSynchronize();
			ErrorCount+=1-checkCUDAError("cuda kernel running: KFPoissonQW");
		}
		if(DBGrid>0)
		{
			KFPoissonDBW<<<DBGrid, DBblock>>>(d_P[0],d_P[1],d_P[2],d_P[3], d_P[4],d_P[5],d_P[6],d_P[7],om2d6, pitchX,pitchXY);
			cudaThreadSynchronize();
			ErrorCount+=1-checkCUDAError("cuda kernel running: KFPoissonDBW");
		}
		if(iteration%PS->ConvergenceCheck==0)
		{
			double OldTotalEnergy=totalEnergy;
			double totalChange;
			double relativeChange;
			double ConvFac;
			totalEnergy=0.0;
			for(i=0;i<8;i++)
			{
				cudaMemcpy(PS->P[i], d_P[i], sizeGS_XYZsplit, cudaMemcpyDeviceToHost);
				for(j=0;j<PS->Qnum[i];j++)
				{
					
					dtmp1=double(PS->P[i][PS->Qpos[i][j]])*double(PS->Q[i][j])/double(PS->Qmult[i][j]);
					
					totalEnergy+=dtmp1;
				}
			}
			totalEnergy=totalEnergy/(fpoh*2.0);
			
			totalChange=totalEnergy-OldTotalEnergy;
			relativeChange=totalChange/totalEnergy;
			
			printf("<PoissonIterations Nit=\"%8d\" E=\"%20.16e\" dE=\"%.4e\" rel.E=\"%.4e\" ConvFac=\"%.4e\"/>\n", iteration, totalEnergy, totalChange, relativeChange,ConvFac);
			if(PS->Tolerance!=0.0 && CollectingTEforAvr==0)
			{
				if(fabs(relativeChange)<=PS->Tolerance)
				{
					if(PS->AvrOverChecks>0)
					{
						printf("Solver has reached the requiered tolerance level\n");
						printf("Collecting total energies for final result\n");
						CollectingTEforAvr=1;
						PS->TotalEnergy=totalEnergy;
						PS->TEavr=0.0;
						PS->stdevTE=0.0;
						AvrTECount=0;
					}
					else
					{
						printf("Solver has reached the requiered tolerance level\n");
						break;
					}
				}
			}
			if(CollectingTEforAvr)
			{
				PS->TEavr+=totalEnergy;
				PS->stdevTE+=(totalEnergy-PS->TotalEnergy)*(totalEnergy-PS->TotalEnergy);
				AvrTECount++;
				if(AvrTECount==PS->AvrOverChecks)
				{
					printf("Finished collecting total energies for final result\n");
					PS->TEavr/=AvrTECount;
					PS->stdevTE/=AvrTECount;
					PS->stdevTE=PS->stdevTE-(PS->TEavr-PS->TotalEnergy)*(PS->TEavr-PS->TotalEnergy);
					//printf("TE=%.7f\tTEavr=%.7f\tstdevTE=%.3g\n",PS->TotalEnergy,PS->TEavr,PS->stdevTE);
					if(PS->stdevTE>0.0)
						PS->stdevTE=sqrt(PS->stdevTE);
					else
						PS->stdevTE=0.0;
					PS->TotalEnergy=totalEnergy;
					printf("TE=%.7f\tTEavr=%.7f\tstdevTE=%.3g\n",PS->TotalEnergy,PS->TEavr,PS->stdevTE);
					break;
				}
			}
		}
		if(ErrorCount)
			break;
	}
	GetCUDADevStat();
	ErrorCount+=1-checkCUDAError("cuda kernel running");
	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop );
	cudaEventElapsedTime( &time, start, stop );
	printf("Time on iterations: %e\n",time);
	time/=1000;//time in s
	printf("\tfor [%d,%d,%d] load store cycle: %f iter/s\n",GS_X,GS_Y,GS_Z,float(PS->MaxIterations)/time);
	//printf("block QGridQ time iter/s %d %d %g %g\n",Qblock,QGrid,time,float(PS->MaxIterations)/time);
	//copy from device lin-array		
	for(i=0;i<8;i++)
		cudaMemcpy(PS->P[i], d_P[i], sizeGS_XYZsplit, cudaMemcpyDeviceToHost);
	
	//free stuff
	for(i=0;i<8;i++)
	{
		if(PS->DielBordNum[i]>0)
		{
			for(k=0;k<6;k++)
			{
				cudaFree(d_DielMult[i*6+k]);
			}
			cudaFree(d_DBPos[i]);
		}
	}
	for(i=0;i<8;i++)
	{
		if(PS->Qnum[i]>0)
		{
			cudaFree(d_Qmult[i]);
			cudaFree(d_Qpos[i]);
			cudaFree(d_Q[i]);
		}
	}
	for(i=0;i<8;i++)
		cudaFree(d_P[i]);
	
	StopClockWMes0("GPU");
	GetCUDADevStat();
	if(ErrorCount)
	{
		printf("During GPU accelerated calculations found %d errors\n",ErrorCount);
		return 0;
	}
	return 1;
}
int checkCUDAError(const char *msg)
{
	cudaError_t err = cudaGetLastError();
	if( cudaSuccess != err) 
	{
		fprintf(stderr, "CUDA ERROR: %s: %s.\n", msg, 
						cudaGetErrorString( err) );
		return 0;
	}
	return 1;
}
int GetCUDADevStat()
{
	printf("########################################\n");
	int device;
	cudaGetDevice(&device);
			
	//cudaSetDevice(device);
	cudaDeviceProp properties;
	cudaGetDeviceProperties(&properties, device);
	int countDev;
	cudaGetDeviceCount(&countDev);
	printf("Totally %d CUDA devices present\n",countDev);
	printf("running on GPU #%d: %s\n", device,properties.name);
	printf("clockRate=%dkHz\n",properties.clockRate);
	printf("computeMode=%d\n",properties.computeMode);
	printf("Major compute capability: %d\n",properties.major);
	printf("Minor compute capability: %d\n",properties.minor);
	printf("Number of multiprocessors on device.: %d\n",properties.multiProcessorCount);
	printf("Maximum pitch in bytes allowed by memory copies: %d\n",properties.memPitch);
	printf("32-bit registers available per block: %d\n",properties.regsPerBlock);
	printf("Shared memory available per block in bytes: %d\n",properties.sharedMemPerBlock);
	printf("Constant memory available on device in bytes: %d\n",properties.totalConstMem);
	printf("Global memory available on device in bytes: %d\n",properties.totalGlobalMem);
	printf("Warp size in threads: %d\n",properties.warpSize);
	//CUdevice device;
	//cuDeviceGet(&device,  ordinal);
	//unsigned int totmem;
	//cuDeviceTotalMem(&totmem, device);
	//printf("memory available on device in bytes: %d\n",totmem);
	unsigned int free, total;
	cuMemGetInfo(&free, &total);
	if (properties.totalGlobalMem >= 1024*1024*1024) {
		printf("Total GPU Memory: %.4f GB\n", properties.totalGlobalMem/(1024.f*1024.f*1024.f) );
	} else {
		printf("Total GPU Memory: %.4f MB\n", properties.totalGlobalMem/(1024.f*1024.f) );
	}
	unsigned int free_mem,total_mem, used_mem;
	cuMemGetInfo( &free_mem, &total_mem );
	used_mem = total_mem-free_mem;
	printf("#CDS1 total mem: %0.3f MB, free: %0.3f MB, used : %0.3f MB\n",
		((double)total_mem)/1024.0/1024.0,
		((double)free_mem )/1024.0/1024.0,
		((double)used_mem )/1024.0/1024.0 );
	printf("#CDS2 total mem: %d, free: %d, used : %d\n",free_mem,total_mem, used_mem);
	printf("########################################\n");
	return 1;
}

#endif
