//
// C++ Interface: pnpstructs
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2006-2018
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "pnpdebug.h"
#include "pnpstructs.h"
GOAtomsStruct* GOAtomsStruct_Create(int GS_X,int GS_Y,int GS_Z,int Natoms,float Rsmoth)
{
	GOAtomsStruct* atms=new GOAtomsStruct;
	atms->GS_X=GS_X;
	atms->GS_Y=GS_Y;
	atms->GS_Z=GS_Z;
	atms->Natoms=Natoms;
	atms->Rsmoth=Rsmoth;
	
	atms->r[0]=new float[Natoms];
	atms->r[1]=new float[Natoms];
	atms->r[2]=new float[Natoms];
	atms->R=new float[Natoms];
	atms->q=new float[Natoms];
	int j;
	for(j=0;j<Natoms;j++)
	{
		atms->r[0][j]=0.0f;
		atms->r[1][j]=0.0f;
		atms->r[2][j]=0.0f;
		atms->R[j]=0.0f;
		atms->q[j]=0.0f;
	}
	return atms;
}

GOAtomsStruct* GOAtomsStruct_Delete(GOAtomsStruct* goatoms)
{
	DeleteCArray(goatoms->r[0]);
	DeleteCArray(goatoms->r[1]);
	DeleteCArray(goatoms->r[2]);
	DeleteCArray(goatoms->R);
	DeleteCArray(goatoms->q);
	
	//DeleteCArray(goatoms->Surf);
	//DeleteCArray(goatoms->iVtmp);
	
	DeleteObjByPnt(goatoms);
	return goatoms;
}


void FieldBW::BorderExchange(const bool *PBC)
{
    int H_X = (GS_X + 1) / 2;
    if (PBC[0]){
        //along X, copy YZ bordering planes
        if (GS_X%2){
            // i.e. GS_X is Odd
            // if G=GS_X the to simulate periodicity
            // we need to copy YZ planes in following x-s
            //
            // phi[ix=0]   = phi[ix G-2]
            // phi[ix=G-1] = phi[ix 0]
            //
            // following scheme should help with setting indexes in 
            // separate black and white arrays
            //                         | padding
            //         B  W  B   W   B | W 
            //         W  B  W   B   W | B
            //         B  W  B   W   B | W 
            //         W  B  W   B   W | B
            // index   0  0     LL   L   L
            // ix      0  1     G-2 G-1
            //
            /*for (int iz = 0; iz < GS_Z; iz++) {
                for (int iy = 0; iy < GS_Y; iy++) {
                    int iu0 = iy * StrideX + iz * StrideXY;
                    int iuL = H_X - 1 + iy * StrideX + iz * StrideXY;
                    int iuLL = iuL - 1;
                    if ((iy + iz) % 2) {
                        //i.e. starts with black
                        B[iu0] = W[iuLL];
                        B[iuL] = W[iu0];
                    }
                    else {
                        //i.e. starts with white
                        W[iu0] = B[iuLL];
                        W[iuL] = B[iu0];
                    }
                }
            }*/
            for (int iz = 0; iz < GS_Z; iz ++) {
                //i.e. starts with black
                for (int iy = (iz + 1) % 2; iy < GS_Y; iy+=2) {
                    int iu0 = iy * StrideX + iz * StrideXY;
                    int iuL = H_X - 1 + iy * StrideX + iz * StrideXY;
                    int iuLL = iuL - 1;
                    B[iu0] = W[iuLL];
                    B[iuL] = W[iu0];
                }
                //i.e. starts with white
                for (int iy = iz % 2; iy < GS_Y; iy+=2) {
                    int iu0 = iy * StrideX + iz * StrideXY;
                    int iuL = H_X - 1 + iy * StrideX + iz * StrideXY;
                    int iuLL = iuL - 1;
                    W[iu0] = B[iuLL];
                    W[iuL] = B[iu0];
                }
            }
        }
        else{
            //i.e. GS_X is Even
            //                         
            //         B  W  B   W   B  W 
            //         W  B  W   B   W  B
            //         B  W  B   W   B  W 
            //         W  B  W   B   W  B
            // index   0  0          L   L
            // ix      0  1         G-2 G-1
            //
            for (int iz = 0; iz < GS_Z; iz++) {
                for (int iy = 0; iy < GS_Y; iy++) {
                    int iu0 = iy * StrideX + iz * StrideXY;
                    int iuL = H_X - 1 + iy * StrideX + iz * StrideXY;
                    if ((iy + iz) % 2) {
                        //i.e. starts with black
                        B[iu0] = B[iuL];
                        W[iuL] = W[iu0];
                    }
                    else {
                        //i.e. starts with white
                        B[iuL] = B[iu0];
                        W[iu0] = W[iuL];
                    }
                }
            }
        }
    }
    if (PBC[1]) {
        //along Y, copy XZ bordering planes
        if (GS_Y % 2) {
            // i.e. GS_Y is Odd
            // To simulate periodicity
            // we need to copy XZ planes in following y-s
            //
            // phi[iy=0]      = phi[iy=GS_Y-2]
            // phi[iy=GS_Y-1] = phi[iy=0]
            //
            // following scheme should help with setting indexes in 
            // separate black and white arrays
            //                         
            //    iy        iz=0
            // GS_Y-1  W  B  W   B   W  B
            // GS_Y-2  B  W  B   W   B  W 
            //         W  B  W   B   W  B
            //     1   B  W  B   W   B  W 
            //     0   W  B  W   B   W  B
            // index   0  0          L   L
            // ix      0  1         G-2 G-1
            //
            for (int iz = 0; iz < GS_Z; iz++) {
                for (int iu = iz * StrideXY; iu < H_X + iz * StrideXY; iu++) {
                    W[iu] = B[iu + (GS_Y - 2) * StrideX];
                    B[iu] = W[iu + (GS_Y - 2) * StrideX];
                    W[iu + (GS_Y - 1) * StrideX] = B[iu + StrideX];
                    B[iu + (GS_Y - 1) * StrideX] = W[iu + StrideX];
                }
            }
        }
        else {
            //i.e. GS_Y is Even
            //    iy        iz=0
            // GS_Y-1  B  W  B   W   B  W 
            // GS_Y-2  W  B  W   B   W  B
            //         B  W  B   W   B  W 
            //         W  B  W   B   W  B
            //     1   B  W  B   W   B  W 
            //     0   W  B  W   B   W  B
            // index   0  0          L   L
            // ix      0  1         G-2 G-1
            //
            for (int iz = 0; iz < GS_Z; iz++) {
                for (int iu = iz * StrideXY; iu < H_X + iz * StrideXY; iu++) {
                    W[iu] = W[iu + (GS_Y - 2) * StrideX];
                    B[iu] = B[iu + (GS_Y - 2) * StrideX];
                    W[iu + (GS_Y - 1) * StrideX] = W[iu + StrideX];
                    B[iu + (GS_Y - 1) * StrideX] = B[iu + StrideX];
                }
            }
        }
    }
    if (PBC[2]) {
        //along Z, copy XY bordering planes
        if (GS_Z % 2) {
            // i.e. GS_Z is Odd
            // To simulate periodicity
            // we need to copy XZ planes in following y-s
            //
            // phi[iz=0]      = phi[iz=GS_Z-2]
            // phi[iz=GS_Z-1] = phi[iz=0]
            //
            // following scheme should help with setting indexes in 
            // separate black and white arrays
            //                         
            //    iz        iy=0
            // GS_Z-1  W  B  W   B   W  B
            // GS_Z-2  B  W  B   W   B  W 
            //         W  B  W   B   W  B
            //     1   B  W  B   W   B  W 
            //     0   W  B  W   B   W  B
            // index   0  0          L   L
            // ix      0  1         G-2 G-1
            //
            for (int iy = 0; iy < GS_Y; iy++) {
                for (int iu = iy * StrideX; iu < H_X + iy * StrideX; iu++) {
                    W[iu] = B[iu + (GS_Z - 2) * StrideXY];
                    B[iu] = W[iu + (GS_Z - 2) * StrideXY];
                    W[iu + (GS_Z - 1) * StrideXY] = B[iu + StrideXY];
                    B[iu + (GS_Z - 1) * StrideXY] = W[iu + StrideXY];
                }
            }
        }
        else {
            //i.e. GS_Y is Even
            //    iy        iz=0
            // GS_Y-1  B  W  B   W   B  W 
            // GS_Y-2  W  B  W   B   W  B
            //         B  W  B   W   B  W 
            //         W  B  W   B   W  B
            //     1   B  W  B   W   B  W 
            //     0   W  B  W   B   W  B
            // index   0  0          L   L
            // ix      0  1         G-2 G-1
            //
            for (int iy = 0; iy < GS_Y; iy++) {
                for (int iu = iy * StrideX; iu < H_X + iy * StrideX; iu++) {
                    W[iu] = W[iu + (GS_Z - 2) * StrideXY];
                    B[iu] = B[iu + (GS_Z - 2) * StrideXY];
                    W[iu + (GS_Z - 1) * StrideXY] = W[iu + StrideXY];
                    B[iu + (GS_Z - 1) * StrideXY] = B[iu + StrideXY];
                }
            }
        }
    }
}
