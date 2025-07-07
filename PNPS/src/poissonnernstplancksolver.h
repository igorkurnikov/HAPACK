//
// C++ Interface: poissonnernstplancksolver
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef POISSONNERNSTPLANCKSOLVER_H
#define POISSONNERNSTPLANCKSOLVER_H

#include "pnps.h"
#include "pnpinterfaces.h"

#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

class ContWorld;
class PoissonSolver;
class NernstPlankSolver;
class VectorField3D;

#include <vector>

/**
	@author Nikolay Simakov <nsimakov@andrew.cmu.edu>
*/
class PoissonNernstPlanckSolver : public PnpsObject, public GenericSolver
{
	public:
		friend class PoissonNernstPlanckMultiGridSolver;
		PoissonNernstPlanckSolver();
		~PoissonNernstPlanckSolver();
		int InitZero();//!<initialize zero values
		int Clear();//!< delete all internal arrays and objects
	public:
		//external variables, parameters for solver
		int MaxIterations;
		/* nernst planck tolerance used to determine pnp convergence,does not work */
		int ConvergenceCheck;
		float tolerance;
		/* dimension alint which to calculate current (x, y, or z) */
		char currentDimension;
		bool verbose;
		int PMFWeightMode;
		std::vector<std::string> PMFWeightModeStr;
		
		//external variables, results
		double Itot,Ipos,Ineg;
		double ItotErr,IposErr,InegErr;
		
		double * positiveCurrentProfile;//! Will be **I in the feture
		double * negativeCurrentProfile;
		
		bool bLimitCurrentCalc;
		float LimitCurrentCalcZ[2];
		bool SaveMemory;
		bool bDouble;
	public:
		//internal variables
		ContWorld *World;
		PoissonSolver *Poisson;
		NernstPlankSolver *NernstPlank;
	public:
		//internal variables
	public:
		//methods
		virtual int LoadParamFromPyDict(PyObject *dict);
		
		int SetContWorld(ContWorld* _world);

		virtual int ShowParameters();//!<show parameters of the solver
    //!show properies can be diffrent from parameters of the solver, valid after InitSolver()
		virtual int ShowProperties(){ return 0;}
    
		virtual int InitSolver();//!<initiate internal arrays
		virtual int Solve();//!< Solve problem
		int SolveSingle();//!< Solve problem
		int SolveDouble();//!< Solve problem
		
		VectorField3D* CalcCartI(int ion);//!< Calculate curent along axis
		VectorField3D* CalcIinout(VectorField3D* CartI);//!< Calculate curent using average from in and out
		VectorField3D* CalcAvrI(VectorField3D* I,int iavr);
	protected:
		int CombineAndPrintCurrents(int iteration);
		int CombineAndPrintCurrentsDouble(int iteration);
};

#ifdef SWIG
%pythoncode %{
def SolvePNPSR(contworld,**kwargs):
	pnps=PoissonNernstPlanckSolver()
	pnps.LoadParamFromPyDict(kwargs);
	pnps.SetContWorld(contworld)
	pnps.InitSolver()
	pnps.Solve()
%}
#endif

#endif
