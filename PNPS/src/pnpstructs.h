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

#include "pnpdebug.h"
#include <assert.h>


#define AligneSize 16

class FieldBW {
public:
    FieldBW()
    {
        B = nullptr;
        W = nullptr;
        Borg = nullptr;
        Worg = nullptr;
        GS_X = 0;
        GS_Y = 0;
        GS_Z = 0;
        StrideX = 0;
        StrideXY = 0;
        Bsize = 0;
        Wsize = 0;
    }
    ~FieldBW()
    {
        DeleteCArray(Borg);
        DeleteCArray(Worg);
    }

public:
    float* B;
    float* W;
    float* Borg;
    float* Worg;
    int GS_X;
    int GS_Y;
    int GS_Z;
    //int GS_XY;
    //int GS_XYZ;
    int StrideX;
    int StrideXY;

    int Bsize;
    int Wsize;

    void Init(int* GridSize)
    {
        GS_X = GridSize[0];
        GS_Y = GridSize[1];
        GS_Z = GridSize[2];
        //GS_XY = GS_X * GS_Y;
        //GS_XYZ = GS_X * GS_Y * GS_Z;

        int H_X = (GS_X + 1) / 2;

        StrideX = H_X;
        StrideXY = StrideX * GS_Y;

        // Only work with odd set
        //assert(GS_X % 2 == 1);
        //assert(GS_Y % 2 == 1);
        //assert(GS_Z % 2 == 1);

        Bsize = StrideX * GS_Y * GS_Z;
        Wsize = StrideX * GS_Y * GS_Z;

        DeleteCArray(B);
        DeleteCArray(W);

        B = new float[Bsize];
        W = new float[Wsize];
        for (int iu = 0; iu < Bsize; ++iu) {
            B[iu] = 0.0f;
            W[iu] = 0.0f;
        }
    }
    bool SameSize(int* GridSize)
    {
        return GS_X == GridSize[0] && GS_Y == GridSize[1] && GS_Z == GridSize[2];
    }
    void SetFromField(const float* F)
    {
        int GS_XY = GS_X * GS_Y;
        for (int iz = 0; iz < GS_Z; iz++) {
            for (int iy = 0; iy < GS_Y; iy++) {
                for (int ix = 0; ix < GS_X; ix++) {
                    int GrdPnt = ix + iy * GS_X + iz * GS_XY;
                    int iu = ix/2 + iy * StrideX + iz * StrideXY;
                    if((ix+iy+iz)%2==1){
                        B[iu] = F[GrdPnt];
                    }
                    else{
                        W[iu] = F[GrdPnt];
                    }
                }
            }
        }
    }
    void SetField(float* F) const
    {
        int GS_XY = GS_X * GS_Y;
        for (int iz = 0; iz < GS_Z; iz++) {
            for (int iy = 0; iy < GS_Y; iy++) {
                for (int ix = 0; ix < GS_X; ix++) {
                    int GrdPnt = ix + iy * GS_X + iz * GS_XY;
                    int iu = ix / 2 + iy * StrideX + iz * StrideXY;
                    if ((ix + iy + iz) % 2 == 1) {
                        F[GrdPnt] = B[iu];
                    }
                    else {
                        F[GrdPnt] = W[iu];
                    }
                }
            }
        }
    }
    void BorderExchange(const bool *PBC);
};


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
