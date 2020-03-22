/*!  \file qc_params.cpp

    Classes to define different parameters of Quantum Chemical Model and computations

    \author Igor Kurnikov 
    \date 2010-

*/

#include <mpi.h>

#include "haconst.h"
#include "hastl.h"
#include "haio.h"
#include "qc_params.h"

#define QC_PARAMS_CPP

IntStrMap SetQCIntEngineTypeLbls()
{
	IntStrMap lmap;
	lmap[QCIntEngineType::INT_ENGINE_GAUSS] = "Gaussian Integral Engine";
	lmap[QCIntEngineType::INT_ENGINE_IPACK] = "IPACK Integral Engine";
	return lmap;
}

IntStrMap QCIntEngineType::labels = SetQCIntEngineTypeLbls();  

QCIntEngineType::QCIntEngineType()
{
	v_= INT_ENGINE_IPACK;	
}
	
QCIntEngineType::~QCIntEngineType()
{

}

int QCIntEngineType::SetWithValue(int val)
{
	v_ = (Value) val;
	return TRUE;
}

int QCIntEngineType::SetWithLabel(const char* label_set)
{
	IntStrMap::iterator itr;
	for(itr = labels.begin(); itr != labels.end(); itr++)
	{
		if( (*itr).second == label_set ) 
		{
			v_ = (Value) (*itr).first;
			return TRUE;
		}
	}
	return FALSE;
}
