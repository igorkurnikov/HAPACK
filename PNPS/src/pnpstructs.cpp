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
#ifndef PNPSTRUCTS_CPP
#define PNPSTRUCTS_CPP

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

#endif
