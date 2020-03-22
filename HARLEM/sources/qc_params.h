/*!  \file qc_params.h

    Classes to define different parameters of Quantum Chemical Model and computations

    \author Igor Kurnikov 
    \date 2010-

*/
#ifndef QC_PARAMS_H
#define QC_PARAMS_H

#include "hastring.h"
#include "hatypes.h"

class QCIntEngineType : public HaEnum1
//! Enum class for Quantum Chemical Integral Engine type 
{
public:
	QCIntEngineType();
	virtual ~QCIntEngineType();

	enum Value { INT_ENGINE_GAUSS = 0, INT_ENGINE_IPACK = 1};
    
	QCIntEngineType& operator=( int value) { SetWithValue(value); return (*this); }
	operator int() const { return v_; }
	bool operator==(const Value& val) const { return v_ == val; }
	bool operator!=(const Value& val) const { return v_ != val; } 

	virtual IntStrMap& GetLabelsMap() const { return labels; }
	virtual int& value() { return (int &) v_; }
	virtual const char* label() const { return labels[v_].c_str();}
	virtual int SetWithValue(int val);
	virtual int SetWithLabel(const char* label);

	QCIntEngineType(Value v): v_(v) {}
private:
	Value v_;
	static IntStrMap labels;
};

namespace swig {
	const QCIntEngineType INT_ENGINE_GAUSS   = QCIntEngineType::INT_ENGINE_GAUSS;
	const QCIntEngineType INT_ENGINE_IPACK   = QCIntEngineType::INT_ENGINE_IPACK;
}






#endif   // End #define QC_PARAMS_H