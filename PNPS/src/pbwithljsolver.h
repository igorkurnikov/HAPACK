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
#ifndef PBWITHLJSOLVER_H
#define PBWITHLJSOLVER_H

#include "pnpinterfaces.h"
#include "pnps.h"

#include <vector>

#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif


class ContWorld;

const int PBwithLJConvFacMaxHistory=5;
/**
	@author Nikolay Simakov <nsimakov@andrew.cmu.edu>
*/
class PBwithLJSolver : public GenericSolver,public PnpsObject
{
	public:
		PBwithLJSolver();
		~PBwithLJSolver();
		friend class ElMod;
		int InitZero();
		int Clear();
	public:
		//external variables, parameters for solver
		enum {Auto=0,NodeIndexBased=1,ArrayDirect=2};
		int solver;///<Solver type {Auto=0,NodeIndexBased=1,ArrayDirect=2}
		int MaxIterations;///<maximum number of iterations
		int ConvergenceCheck;///<how often check energies and convergence
		float Convergence;///<if change less then convergence then stop
		float Relaxation;///<relaxation paramter, if MinIterations < 0 will estimate number of iterations and Relaxation,if Relaxation<0 will estimate it
		bool verbose;
		bool bPBExpAll;
		bool bAnalyseExplosion;
		//external variables, results
		double totalChange,relativeChange;
		double totalEnergy,totalEnergyInd;
		double ConvFac;
	public:
		//internal variables
		ContWorld* m_ContWorld;
		int NoSingularNum[3];
		int *IndexNoSingular;
		int SingularNum[3];
		int *IndexSingular;
		float *dielectricXS,*dielectricYS,*dielectricZS,*dielectricZSSUM;
		float *dielectricXmS,*dielectricYmS,*dielectricZmS;
		float *QstS;
		float *PhiSingular;
		int DielBoarderNum[3];
		int *IndexDielBoarder;
		float *dielectricXDB,*dielectricYDB,*dielectricZDB,*dielectricZDBSUM;
		float *dielectricXmDB,*dielectricYmDB,*dielectricZmDB;
		int ChargeNum[3];
		int *IndexCharge;
		float *dielectricCh;
		float *Qst;
		float *PhiCharge;
		
		int PBZoneNum[3];
		int *IndexPBZone;
		
		int PBLJZoneNum[3];
		int *IndexPBLJZone;
		
		int PBLJQstZoneNum[3];
		int *IndexPBLJQstZone;
		float *dielectricPBLJQst;
		float *QstPBLJ;
		
		int PBLJDBZoneNum[3];
		int *IndexPBLJDBZone;
		float *dielectricXPBLJDB,*dielectricYPBLJDB,*dielectricZPBLJDB;
		float *dielectricZPBLJDBSUM;
		float *dielectricXmPBLJDB,*dielectricYmPBLJDB,*dielectricZmPBLJDB;
		
		int PBLJDBQstZoneNum[3];
		int *IndexPBLJDBQstZone;
		float *dielectricXPBLJDBQst,*dielectricYPBLJDBQst,*dielectricZPBLJDBQst;
		float *dielectricZPBLJDBQstSUM;
		float *dielectricXmPBLJDBQst,*dielectricYmPBLJDBQst,*dielectricZmPBLJDBQst;
		float *QstPBLJDB;
		
		
		//!< NodeTypes for AD solver
		enum NodeTypeAD {NoSingular=0, Boarder=1, Charge=2, DielBoarder=3, ChargeAndDielBoarder=4, Singular=5};
		
	protected:
		//internal variables
		//!variables for solver, setuped in InitZero(CWorld *world, Data* Dt);
		int GS_X;
		int GS_Y;
		int GS_Z;
		int GS_XY;
		int GS_XYZ;
		float GridScale;
		float *potential;//!< World->Potential
		
    //!variables for solver, setuped in SetParameters(Data* Dt);
		float om2;//!< relaxation
		float om1;//!< 1.0-om2
		float om2d6;//!< om2/6.0
		int SetRelaxation(float _Relaxation);
		
		float ConvFacHistory[PBwithLJConvFacMaxHistory];
	public:
		//methods
		virtual int LoadParamFromPyDict(PyObject *dict);
		
		int SetContWorld(ContWorld* _world);
		
		
		virtual int InitSolver();//!<initiate internal arrays
		int InitSolverNIB();
		//int InitSolverAD();
		
		virtual int Solve();//!< Solve problem
		int AnalyseExplosion();
		
		virtual int CalcSystemEnergy(int iteration);
		virtual int CalcSystemEnergyStdDevPhi(int iteration);
		//virtual float CalculateEnergyPAD(float fpoh,float *Potential,float *StaticCharge,float *Epsilon,int *IndexCharge, int *IndexSingular,int ChargeNum,int SingularNum);
		virtual int ShowParameters();//!<show parameters of the solver
    //!show properies can be diffrent from parameters of the solver, valid after InitSolver()
		virtual int ShowProperties();
	private:
		std::vector<std::string> SolverStr;
		double* ChargesEnergy;

};
#ifdef SWIG
%pythoncode %{
#pnpmod.SolvePBSR(contworld, MaxInterations=100, Tolerance=0.0, Relaxation=1.6)
def SolvePBSR(contworld,**kwargs):
	pbsr=PBwithLJSolver()

	pbsr.MaxIterations = kwargs.get("MaxIterations",200)
	pbsr.Convergence = kwargs.get("Tolerance",0.0)
	pbsr.Relaxation = kwargs.get("Relaxation",1.6)
	pbsr.ConvergenceCheck = kwargs.get("ConvergenceCheck",20)
	pbsr.solver = kwargs.get("solver",0)
	pbsr.verbose = kwargs.get("verbose",True)
	pbsr.bPBExpAll = kwargs.get("bPBExpAll",False)
	pbsr.bAnalyseExplosion = kwargs.get("bAnalyseExplosion",False)

	pbsr.ShowParameters()
	#pbsr.SetRelaxation(pbsr.Relaxation)

	
	pbsr.SetContWorld(contworld)
	pbsr.InitSolver()
	pbsr.Solve()
	del pbsr
	
	

%}
#endif
#endif
