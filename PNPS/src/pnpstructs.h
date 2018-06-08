//
// C++ Interface: pnpstructs
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PNPSTRUCTS_H
#define PNPSTRUCTS_H

#define PlusX 0
#define MinusX 1
#define PlusY 2
#define MinusY 3
#define PlusZ 4
#define MinusZ 5

#define CUDAXTRAX 16



#if defined(WITH_CUDA)
#include <vector_types.h>
#else
struct float4
{
	float x, y, z, w;
};
struct float3
{
	float x, y, z;
};
struct int3
{
	int x, y, z;
};
struct int4
{
	int x, y, z, w;
};
#endif

typedef struct 
{
	int GS[3];//!<GridSize Without BC
	//!Size of splitted potential with bouders for CUDA optimal access
	//!PS.P* is store in pitched array x has 16 for CUDA and Y/Z is +2 for BC
	float GridScale;
	int spltGSWBC[3];
	int spltGSWBC_X;
	int spltGSWBC_XY;
	int spltGSWBC_XYZ;
	float* P[8];//!<splitted potential
	int MaxIterations;///<maximum number of iterations
	float Relaxation;///<relaxation paramter
	
	int Qnum[8];
	float* Q[8];
	float* Qmult[8];
	int* Qpos[8];
	
	int DielBordNum[8];
	float* DielMult[8][6];
	int* DielBordPos[8];
	
	int ConvergenceCheck;
	double Tolerance;
	
	double TotalEnergy;
	int AvrOverChecks;
	double TEavr;
	double stdevTE;
} PoissonSolverOnCudaStruct;

typedef struct
{
	int GS[3];//!<GridSize Without BC
	//!Size of splitted potential with bouders for CUDA optimal access
	//!PS.P* is store in pitched array x has 16 for CUDA and Y/Z is +2 for BC
	float GridScale;
	int spltGSWBC[3];
	int spltGSWBC_X;
	int spltGSWBC_XY;
	int spltGSWBC_XYZ;
	double* P[8];//!<splitted potential
	int MaxIterations;///<maximum number of iterations
	double Relaxation;///<relaxation paramter
	
	int Qnum[8];
	double* Q[8];
	double* Qmult[8];
	int* Qpos[8];
	
	int DielBordNum[8];
	double* DielMult[8][6];
	int* DielBordPos[8];
	
	int ConvergenceCheck;
	double Tolerance;
	
	double TotalEnergy;
	int AvrOverChecks;
	double TEavr;
	double stdevTE;
} PSolverOnCudaStructDouble;

//!PoissonSolverOnCuda4Struct
typedef struct 
{
	int3 GS;//!<GridSize Without BC
	//!Size of splitted potential with bouders for CUDA optimal access
	//!PS.P* is store in pitched array x has 16 for CUDA and Y/Z is +2 for BC
	int3 ps4GS;
	float GridScale;
	
	float* Pot;//!<regular potential
	//!splitted potential
	float  *P000 ,  *P100 ,  *P200 ,  *P300 ;
	float  *P010 ,  *P110 ,  *P210 ,  *P310 ;
	float  *P020 ,  *P120 ,  *P220 ,  *P320 ;
	float  *P030 ,  *P130 ,  *P230 ,  *P330 ;

	float  *P001 ,  *P101 ,  *P201 ,  *P301 ;
	float  *P011 ,  *P111 ,  *P211 ,  *P311 ;
	float  *P021 ,  *P121 ,  *P221 ,  *P321 ;
	float  *P031 ,  *P131 ,  *P231 ,  *P331 ;

	float  *P002 ,  *P102 ,  *P202 ,  *P302 ;
	float  *P012 ,  *P112 ,  *P212 ,  *P312 ;
	float  *P022 ,  *P122 ,  *P222 ,  *P322 ;
	float  *P032 ,  *P132 ,  *P232 ,  *P332 ;

	float  *P003 ,  *P103 ,  *P203 ,  *P303 ;
	float  *P013 ,  *P113 ,  *P213 ,  *P313 ;
	float  *P023 ,  *P123 ,  *P223 ,  *P323 ;
	float  *P033 ,  *P133 ,  *P233 ,  *P333 ;
	
	int3 spltGSWBC;
	int spltGSWBC_X;
	int spltGSWBC_XY;
	int spltGSWBC_XYZ;
	
	int MaxIter;///<maximum number of iterations
	float Rel;///<relaxation paramter
	double Tol;///<Tolerance
	
	int Qnum[8];
	float* Q[8];
	float* Qmult[8];
	int* Qpos[8];
	
	int DielBordNum[8];
	float* DielMult[8][6];
	int* DielBordPos[8];
	
	int ConvergenceCheck;
	
	double TotalEnergy;
	int AvrOverChecks;
	double TEavr;
	double stdevTE;
	float om1;
	float om2d6;
} PoissonSolverOnCuda4Struct;
typedef struct 
{
	int3 GS;//!<GridSize Without BC
	//!Size of splitted potential with bouders for CUDA optimal access
	//!PS.P* is store in pitched array x has 16 for CUDA and Y/Z is +2 for BC
	float GridScale;
	
	float* Pot;//!<regular potential
	//!splitted potential
	float  *PotCu;
	
	int3 spltGSWBC;
	int spltGSWBC_X;
	int spltGSWBC_XY;
	int spltGSWBC_XYZ;
	
	int MaxIter;///<maximum number of iterations
	float Rel;///<relaxation paramter
	double Tol;///<Tolerance
	
	float om1;
	float om2d6;
} PoissonSolverOnCuda1Struct;

struct PoissonSolverOnCudaParamStruct
{
	//!Block of Threads size for laplace
	int BS_X;
	int BS_Y;
	int BS_Z;
	int BS_XY;
	int BS_XYZ;
	
	int Qblock;
	int DBblock;
};
typedef struct
{
	int GS_X;
	int GS_Y;
	int GS_Z;
	int Natoms;
	float *r[3];
	float *R;
	float *q;
	
	float Rsmoth;
	
	int iDiel;
	int iDielBulk;
	
	float* Surf;
	int* iVtmp;
} GOAtomsStruct;
GOAtomsStruct* GOAtomsStruct_Create(int GS_X,int GS_Y,int GS_Z,int Natoms,float Rsmoth);
GOAtomsStruct* GOAtomsStruct_Delete(GOAtomsStruct* goatoms);
#endif
