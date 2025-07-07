//
// C++ Interface: poissonsolver
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SISSGAUSPS_H
#define SISSGAUSPS_H

#include "pnpinterfaces.h"
#include "pnps.h"

class ContWorld;

typedef long int DielMatElt;
/**
	@author Nikolay Simakov <nsimakov@andrew.cmu.edu>
*/
class SISSGausPS : public GenericSolver,public PnpsObject
{
	public:
		SISSGausPS();
		~SISSGausPS();
		int InitZero();
		int Clear();
		
	public:
		//internal variables
		ContWorld* World;
		
	protected:
		//internal variables
		//!variables for solver, setuped in InitZero(CWorld *world, Data* Dt);
		int GS_X;
		int GS_Y;
		int GS_Z;
		int GS_XY;
		int GS_XYZ;
		int lGS_X;
		int lGS_Y;
		int lGS_Z;
		int lGS_XY;
		int lGS_XYZ;
		int DMepslen;
		float GridScale;
		
		int Qnum;
		int *Qpos;
		float *Qval;
		DielMatElt **QMult;
		DielMatElt **DielMat;
		
		DielMatElt **rDielMat;
		float *potential;//!< World->Potential
		
		int *PrimeNumbers;
		int NumOfPrimeNumbers;
		
	public:
		//methods
		int SetContWorld(ContWorld* _world);
		virtual int InitSolver();//!<initiate internal arrays
		virtual int Solve();//!< Solve problem
		
		virtual int ShowParameters();//!<show parameters of the solver
    //!show properies can be diffrent from parameters of the solver, valid after InitSolver()
		virtual int ShowProperties();
		
		int PrintDielMat();
		int PrintRDielMat();
		int PrintDielMatCurState();
		int PrintRDielMatCurState();
		int AddRowItoRowJGoDown(int I,DielMatElt multI,int J,DielMatElt multJ);//!< AddRowItoRowJGoDown store results in J,always dIJ=J-I>0and dIJ>=DMepslen
		int AddRowItoRowJGoUp(int I,DielMatElt multI,int J,DielMatElt multJ);
		bool IsItComMultForRowJ(int J,DielMatElt Mult);
		int DevRowJByComMult(int J,DielMatElt Mult);//!<just do it do not check
		int FindAndDevByComMult(int J);
		
		int RexchangeRows(int I, int J);
		int RAddRowItoRowJ(int I,DielMatElt multI,int J,DielMatElt multJ);
		bool RIsItComMultForRowJ(int J,DielMatElt Mult);
		int RDevRowJByComMult(int J,DielMatElt Mult);//!<just do it do not check
		int RFindAndDevByComMult(int J);
};

#endif
