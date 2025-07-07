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
#ifndef POISSONSOLVER_H
#define POISSONSOLVER_H

#include <vector>

#include "pnpinterfaces.h"
#include "pnps.h"

#include "pnpstructs.h"

#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

class ContWorld;

//!definition for Warent code
#define POISSON_SOLVER_CONVERGENCE_CHECK 			1
#define POISSON_SOLVER_CHARGED_POINT				1
#define POISSON_SOLVER_UNCOMMON_DIELECTRIC_POINT		2
#define POISSON_SOLVER_CHARGED_UNCOMMON_DIELECTRIC_POINT	3
//!Conteiner for CPoisson to run Warent code
typedef struct _PoissonSolverData {
	/* grid points with charged points or uncommon dielectric points around
	* it
	*/
	long * borderPoints;
	/* inverse of the dielectric sum multiplied by the relaxation paramter
	* for borderPoints with uncommon dielectric points
	*/
	float * om2InverseDielectricSum;
	/* type of border point, uncommon dielectric, charged, or both */
	short * typeOfBorderPoint;
	/* maximum iterations for solving poisson */
	int maxIterations;
	/* stop when rmsChange < convergence */
	float convergence;
	/* relaxation paramter */
	float relaxation;
	/* do we have dynamic charges? */
	int hasDynamicCharges;
} PoissonSolverDataW;

/**
	@author Nikolay Simakov <nsimakov@andrew.cmu.edu>
*/
class PoissonSolver : public GenericSolver,public PnpsObject
{
	public:
		PoissonSolver();
		~PoissonSolver();
		friend class ElMod;
		int InitZero();
		int Clear();
	public:
		//external variables, parameters for solver
		enum {Auto=0,NodeIndexBased=1,ArrayDirect=2,PNPC=3};
		int solver;///<Solver type {Auto=0,NodeIndexBased=1,ArrayDirect=2,PNPC=3}
		int MaxIterations;///<maximum number of iterations
		int MinIterations;///<minimum number of iterations, if <0 will estimate number of iterations and Relaxation
		int ConvergenceCheck;///<how often check energies and convergence
		float Convergence;///<if change less then convergence then stop
		float Relaxation;///<relaxation paramter, if MinIterations < 0 will estimate number of iterations and Relaxation,if Relaxation<0 will estimate it
		int QmobMod;///<0-off, 1-on 2 Auto
		bool verbose;
		//external variables, results
		double totalChange,relativeChange;
		double totalEnergy,totalEnergyInd;
		double ConvFac;
		
		
		int WayToCalcSystemEnergy;//!<0-Auto,1-CalcSystemEnergyStdDevPhi
		
		int ConvFacMaxHistory;
		
	public:
		//internal variables
		ContWorld* World;
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
		
		float *ChargeSum;//!<SummaryCharge for AD solver
		//!< NodeTypes for AD solver
		enum NodeTypeAD {NoSingular=0, Boarder=1, Charge=2, DielBoarder=3, ChargeAndDielBoarder=4, Singular=5};
		int QmobNum[3];
		int *IndexQmob;
		float *Qmob;
		float *dielectricChMob;
		
		//bool *QmobFlag;
		bool *CalcVolume;//!<volume there run calculation currently not working
    
		int QmobDielBoarderNum[3];
		int *IndexQmobDielBoarder;
		float *QmobDielBoarder;
		
		float *dielectricXQmobDB;
		float *dielectricYQmobDB;
		float *dielectricZQmobDB;
		float *dielectricXmQmobDB;
		float *dielectricYmQmobDB;
		float *dielectricZmQmobDB;
		float *dielectricZQmobDBSUM;
		
		
		int QmobDielBoarderQstNum[3];
		int *IndexQmobDielBoarderQst;
		float *QmobDielBoarderQst;
		float *QstQmobDielBoarderQst;
		
		float *dielectricXQmobDBQst;
		float *dielectricYQmobDBQst;
		float *dielectricZQmobDBQst;
		float *dielectricXmQmobDBQst;
		float *dielectricYmQmobDBQst;
		float *dielectricZmQmobDBQst;
		float *dielectricZQmobDBSUMQst;
		///!to have Waren's code
		PoissonSolverDataW* poissonSolverData;
	protected:
		//internal variables
		//!variables for solver, setuped in InitZero(CWorld *world, Data* Dt);
		int GS_X;
		int GS_Y;
		int GS_Z;
		int GS_XY;
		int GS_XYZ;
		float GridScale;
#	ifndef PNPDOUBLE
		float *potential;//!< World->Potential
        FieldBW *PotBW;
#	else
		double *potential;//!< World->Potential
        FieldBW *PotBW;
#	endif
		
    //!variables for solver, setuped in SetParameters(Data* Dt);
		float om2;//!< relaxation
		float om1;//!< 1.0-om2
		float om2d6;//!< om2/6.0
		int GuessNumberOfIteration();
		int SetRelaxation(float _Relaxation);
	public:
		//methods
		virtual int LoadParamFromPyDict(PyObject *dict);
		int SetContWorld(ContWorld* _world);
		
		
		virtual int InitSolver();//!<initiate internal arrays
		int InitSolverNIB();
		int InitSolverW();
		int InitSolverAD();
		
		int SetQmobForPNP();//!<init
		int SetQmobFromConcentration();//!set particular charge
		int SetQmobFromConcentrationDouble();//!set particular charge
		
		virtual int Solve();//!< Solve problem
		int PoissonSolverNIB(bool ckenergy=true);
		int PoissonSolverAD();
		int PoissonSolverW();
		
		virtual int CalcSystemEnergy(int iteration);
		virtual int CalcSystemEnergyFloat(int iteration);
		virtual int CalcSystemEnergyDouble(int iteration);
		virtual int CalcSystemEnergyLongDouble(int iteration);
		virtual int CalcSystemEnergyAnalizer(int iteration);
		virtual int CalcSystemEnergy0(int iteration);
		virtual int CalcSystemEnergy1(int iteration);
		virtual int CalcSystemEnergy3(int iteration);
		//virtual int CalcSystemEnergy4PMF(int iteration);
		virtual int CalcSystemEnergyMaxPhiChange(int iteration);
		virtual int CalcSystemEnergyStdDevPhi(int iteration);
		//!Calculation of E for AD Solver
		float CalculateEnergyPAD(float fpoh,float *Potential,float *StaticCharge,float *Epsilon,int *IndexCharge, int *IndexSingular,int ChargeNum,int SingularNum);
		
		virtual int ShowParameters();//!<show parameters of the solver
    //!show properies can be diffrent from parameters of the solver, valid after InitSolver()
		virtual int ShowProperties();
	private:
		std::vector<std::string> SolverStr;
		double* ChargesEnergy;

};
#endif
