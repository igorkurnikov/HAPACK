//
// C++ Interface: pnputil
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PNPUTIL_H
#define PNPUTIL_H
#include <vector>
#include <string>
#ifdef HARLEM_MOD
#include <halinalg.h>
#endif
class ContWorld;
class TiXmlElement;
class VectorField3D;
class PBwithLJSolver;
class NernstPlankSolver;

class PNPUtil
{
	public:
		PNPUtil();
		~PNPUtil();
		static int RemoveLargePMFfromPNP(ContWorld *world,float LargePMF);
		static int RemoveSmallCfromPNP(ContWorld *world,float SmallC);
		static int RemoveSmallCfromPNPNearDiffBoarder(ContWorld *world,float SmallC);
		static int RemoveLargedPMFfromPNP(ContWorld *m_ContWorld,float dPMF);
		
		static int RemoveLargedPMFandLargePMFfromPNP(ContWorld *m_ContWorld,float dPMF,float PMF);
		static int RemoveLargedPMFandLargePMFfromPNPOld(ContWorld *m_ContWorld,float dPMF,float PMF);
		static int RemoveCavitiesAtDiffusionMap(ContWorld *m_ContWorld);
#		ifdef HARLEM_MOD
		static HaVec_float GetRlim(ContWorld *world,const TiXmlElement *RlimElt);
#		endif
		static int ScaleDiffusionInTheChannel(ContWorld *world,float x0,float y0,std::vector<float>* zRlim,std::vector<float>* Rlim,std::vector<float>* zDiffScale,std::vector<float>* DiffScale);
		static int ScaleDiffusionInTheChannel(ContWorld *world,float x0,float y0,float* Rlim,float *DiffScale);
		
		static int ConvertPBLJresultsToDynamicCharge(ContWorld *m_ContWorld);
		
		static int SetInternalCtoZero(ContWorld *m_ContWorld);
		static int SetCtoZeroWhereDZero(ContWorld *m_ContWorld);
		static int SetDzeroAtEps(ContWorld *m_ContWorld,int iEps);
		static int SetDZeroWithIndex(ContWorld* _ContWorld);
		static int RemoveQfromNI(ContWorld *m_ContWorld);
		
		//!will return number of removed nodes
		static int RemoveNodesFromPNPBasedOnNPCriteria(ContWorld *m_ContWorld,float Relaxation, int MaxCycles,float MaxdC, bool RemoveNegC);
		//static bool ProcessPNPUtilCmds(ContWorld *world,const TiXmlElement *Elt);
	protected:
		static bool ShouldRemoveLargedPMFandLargePMFfromPNP(ContWorld *m_ContWorld,float dPMF,float PMF,int ion, int GrdPnt);
		static int ShouldRemoveNodeBasedOnNPCriteria(ContWorld *m_ContWorld,int ion, int GrdPnt,float Relaxation,float MaxdC, bool RemoveNegC);
};
class NPMaskBuilder
{
	public:
		NPMaskBuilder(ContWorld *world);
		~NPMaskBuilder();
		int CmdNPMaskBuilder(const TiXmlElement *Elt);
		
		int InitNPMask();
		int SetToNIDiffusion();
		//!<All in Internal Units
		int RemoveTubeRegion(int ion,float X,float Y,float Z0, float Z1,float R0,float R2);
		int RemoveHighC(float HighC);
		
		int ReadNPMask(const char * filename);
		int WriteNPMask(const char * filename);
		
		ContWorld *m_ContWorld;
		VectorField3D *NPMask;
};

//Python helpers
#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

#ifndef SWIG


int haPyDict_GetItemAsBool(PyObject *dict, const char *key, bool *v);
int haPyDict_GetItemAsInt(PyObject *dict, const char *key, int *v);
int haPyDict_GetItemAsFloat(PyObject *dict, const char *key, float *v);
int haPyDict_GetItemAsString(PyObject *dict, const char *key, std::string *v);

bool haPyDict_GetItemValueAsBool(PyObject *dict, const char *key, bool vdefault);
int haPyDict_GetItemValueAsInt(PyObject *dict, const char *key, int vdefault);
float haPyDict_GetItemValueAsFloat(PyObject *dict, const char *key, float vdefault);
char* haPyDict_GetItemValueAsString(PyObject *dict, const char *key, const char *vdefault);

int haPy_SetCFloatArrFromListOfFloat(PyObject *vlist, float *carr);
#endif

//! Calculator of ion accessible volume
class IAVCalc
{
	public:
		IAVCalc(ContWorld *world);
		~IAVCalc();
		int CmdIAVCalcDict(PyObject *dict);
	protected:
		ContWorld *m_ContWorld;
		PBwithLJSolver *PBLJ;
		NernstPlankSolver *NP;
#ifdef SWIG
//	%pythoncode %{
//		%}
#endif
};

#ifdef SWIG
%pythoncode %{
def RefineIAV(contworld,**kwarg):
	"""
            MaxdC=0.05,
            MaxCycle=100000,
            PBSR_Param={"MaxIterations":100, "Tolerance":2.0e-06, "Relaxation":1.6},
            NP_Param={"Relaxation":1.0}
	"""
	iavcalc=IAVCalc(contworld)
	iavcalc.CmdIAVCalcDict(kwarg)
	del iavcalc
	print "RefineIAV"
%}
#endif
#endif
