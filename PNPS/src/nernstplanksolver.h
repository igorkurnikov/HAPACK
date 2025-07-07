//
// C++ Interface: nernstplanksolver
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef NERNSTPLANKSOLVER_H
#define NERNSTPLANKSOLVER_H

#include "pnps.h"
#include "pnpinterfaces.h"
#include "mapio.h"

#include <vector>

#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

class ContWorld;
class VectorField3D;
class IAVCalc;
//class HaPyDoubleVectorField3D;
//!definition for Warent code
#define NERNST_PLANCK_SOLVER_CONVERGENCE_CHECK 	5
#define NERNST_PLANCK_SOLVER_DELPHI_TOLERANCE		0.02
typedef struct _NernstPlanckSolverDataW {
	/* maximum number of nernst planck iterations */
	int maxIterations;
	/* done when rmsChange < convergence */
	float convergence;
	/* relaxation paramter */
	float relaxation;
	/* use non-linear interpolation scheme for conc at midpoints*/
	int expInterpolation;
	/* array to hold gridPoints to start solving NP */
	long * start;
	/* array to hold gridPoints to end solving NP */
	long * end;
	/* grid points where surrounding grid points are different */
	long * borderPoints;
	/* diffusion constants for borderPoints */
	float * dix[6];
	/* sum of diffusion constants for borderPoints */
	float * dixt;
} NernstPlanckSolverDataW;

/**
	@author Nikolay Simakov <nsimakov@andrew.cmu.edu>
*/
class NernstPlankSolver : public PnpsObject, public GenericSolver
{
	public:
		NernstPlankSolver();
		~NernstPlankSolver();
		
		friend class PoissonNernstPlanckMultiGridSolver;
		friend class PoissonNernstPlanckSolver;
		friend class IAVCalc;
		int InitZero();
		int Clear();
	public:
		//external variables, parameters for solver
		int MaxIterations;
		float Convergence;
		int ConvergenceCheck;///<how often check energies and convergence
		//int CurrentCalc;///<calculate current every CurrentCalc iterations
		float Relaxation;
		std::string QmobMaskFile;//!<QmobMaskFile
		int solver;
		bool verbose;
		float PMFWeight;
		//external variables, results
	public:
		//internal variables
		ContWorld *World;
		//! diffusion constants for borderPoints
		float * dix[2][6], *dphi[6];
		//! sum of diffusion constants for borderPoints
		float * dixt[2];
		float MaxChange;
    
		int NoSingularNum[4][3];
		int *IndexNoSingular[2];
		int SingularNum[4][3];
		int *IndexSingular[2];
	public:
		//internal variables
		NernstPlanckSolverDataW* nernstPlanckSolverData;
		//optional stuf
		int** QmobMask;
		bool *CalcVolume;//!<volume there run calculation
		bool writeStatQmob;
		bool bConcCorrection;
		//bool bDouble;
		int PntsWithNegativeConc[10];
	public:
		enum {Auto=0,NodeIndexBased=1,ArrayDirect=2,PNPC=3};
		std::vector<std::string> SolverStr;
	public:
		float *UTMPSingle;
		float *TMPSingle;
		double *UTMPDouble;
		double *TMPDouble;
	public:
		//methods
		virtual int SaveXML(TiXmlElement* Elt, HaContext* p_ctxt = NULL );
		virtual int LoadXML(const TiXmlElement* Elt, HaContext* p_ctxt = NULL );
		virtual int LoadParamFromPyDict(PyObject *dict);
		int SetContWorld(ContWorld* _world);
		
		int InitSolver();
		int InitSolverAD();
		int InitSolverD();
		int InitSolverDflex();
		int InitSolverW();
		
		int ShowParameters();
		int ShowProperties(){ return 0;}
	
		int Solve();
	public:
		int NernstPlanckSolverD();
		int NernstPlanckSolverDDouble();
	protected:
		int NernstPlanckSolverDpre(int IonType);
		int NernstPlanckSolverDiter(int IonType,bool calcchange);
		int NernstPlanckSolverDiterHalfStep(int IonType,bool calcchange,int BlackOrWhite);
		int NernstPlanckSolverDpost(int IonType);
		int NernstPlanckSolverDpreDouble(int IonType);
		int NernstPlanckSolverDiterDouble(int IonType,bool calcchange);
		int NernstPlanckSolverDpostDouble(int IonType);
	public:	
		int NernstPlanckSolverDOld();
		int NernstPlanckSolverN2tmp();
		int NernstPlanckSolverW();
		
		float nernstPlanckSolverTolerance();
		void nernstPlanckSolverCurrentProfile(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile);
		void nernstPlanckSolverCurrentProfileW(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile);
		void nernstPlanckSolverCurrentProfileN(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile);
		void nernstPlanckSolverCurrentProfileAD(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile);
		void nernstPlanckSolverCurrentProfileDflex(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile);
		VectorField3D* CalcJ();
		
		int CheckSystem();
};
#ifdef PNPDOUBLE
class NernstPlankSolverDouble : public HaObject, public GenericSolver
{
	public:
		NernstPlankSolverDouble();
		~NernstPlankSolverDouble();
		
		int InitZero();
		int Clear();
	public:
		//external variables, parameters for solver
		int MaxIterations;
		double Convergence;
		int ConvergenceCheck;///<how often check energies and convergence
		double Relaxation;
		std::string QmobMaskFile;//!<QmobMaskFile
		int solver;
		bool verbose;
		double PMFWeight;
		//external variables, results
	public:
		//internal variables
		ContWorld *World;
		//! diffusion constants for borderPoints
		double * dix[2][6], *dphi[2][6];
		//! sum of diffusion constants for borderPoints
		double * dixt[2];
		double MaxChange;
    
		int NoSingularNum[4][3];
		int *IndexNoSingular[2];
		int SingularNum[4][3];
		int *IndexSingular[2];
	public:
		//internal variables
		NernstPlanckSolverDataW* nernstPlanckSolverData;
		//optional stuf
		int* QmobMask;
	public:
		enum {Auto=0,NodeIndexBased=1,ArrayDirect=2,PNPC=3};
		std::vector<std::string> SolverStr;
	public:
		double *UTMP[2];
		double *TMP[2];
		double *locD;
	public:
		//methods
		virtual int SaveXML(TiXmlElement* Elt, HaContext* p_ctxt = NULL );
		virtual int LoadXML(const TiXmlElement* Elt, HaContext* p_ctxt = NULL );
		
		int SetContWorld(ContWorld* _world);
		
		int InitSolver();
		int InitSolverAD();
		int InitSolverD();
		int InitSolverDflex();
		
		int ShowParameters();
		int ShowProperties(){ return 0;}
		
		int Solve();
		int NernstPlanckSolverD();
		int NernstPlanckSolverDOld();
		int NernstPlanckSolverN2tmp();
		
		double nernstPlanckSolverTolerance();
		void nernstPlanckSolverCurrentProfile(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile);
		void nernstPlanckSolverCurrentProfileN(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile);
		void nernstPlanckSolverCurrentProfileAD(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile);
		void nernstPlanckSolverCurrentProfileDflex(int dim, double * positiveCurrentProfile, double * negativeCurrentProfile);
		HaPyDoubleVectorField3D* CalcJ();
		
		int CheckSystem();
};
#endif
#endif
