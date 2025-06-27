/*! \file hastl.h

    Machine Dependent header for STL Template library

    \date 1997-2002
    \author Igor Kurnikov 

*/
#ifndef HASTL_H
#define HASTL_H

#ifdef _MSC_VER
// disable warning about long names for MS C++ compiler
#pragma warning (disable:4786)
#endif

#include <queue>
#if !defined(linux)
#   include <limits>
#endif
#include <stack>
#include <list>
#include <map>
#include <vector>
#include <set>
#include <string>

#ifdef linux
#include <new>
#endif

#if !defined(__PRETTY_FUNCTION__) 
#if defined(_MSC_VER) 
#define __PRETTY_FUNCTION__ __FUNCSIG__
#else
#define __PRETTY_FUNCTION__ __func__
#endif
#endif


#if defined(SWIG)
 
class IntIntMap_parent
{
public:
	
	void clear();
	int count(int ikey);
	int size();
	bool empty();
};  


 
class IntPtrMap_parent
{
public:
		
	void clear();
	int count(int ikey);
	int size();
	bool empty();
};


class PtrIntMap_parent
{
public:
		
	void clear();
	int count(void* ptr);
	int size();
	bool empty();
};

class PtrPtrMap_parent
{
public:
		
	void clear();
	int count(void* ptr);
	int size();
	bool empty();
};

 
class PtrDoubleMap_parent
{
public:
		
	void clear();
	int count(void* ptr);
	int size();
	bool empty();
};


#else

typedef std::map<int, int>              IntIntMap_parent;
typedef std::map<int, int>::iterator    IntIntMap_itr_parent;
typedef std::map<int, void*>            IntPtrMap_parent;
typedef std::map<int, void*>::iterator  IntPtrMap_itr_parent;
typedef std::map<void*, int>            PtrIntMap_parent;
typedef std::map<void*, int>::iterator  PtrIntMap_itr_parent;
typedef std::map<void*, void*>           PtrPtrMap_parent;
typedef std::map<void*, void*>::iterator PtrPtrMap_itr_parent;
typedef std::map<void*, double>           PtrDoubleMap_parent;
typedef std::map<void*, double>::iterator PtrDoubleMap_itr_parent;

#endif

class IntIntMap_itr: public IntIntMap_itr_parent
//! Iterator class of a map of integers to integers 
{
public:
	IntIntMap_itr() {}
	IntIntMap_itr( const IntIntMap_itr_parent& ref);

	int GetKey();
	int GetVal();

    bool operator==(const IntIntMap_itr& ref) const;
	bool operator!=(const IntIntMap_itr& ref) const;	
	
	IntIntMap_itr& next();
};


class IntIntMap: public IntIntMap_parent
//! Map of integers to integers 
{
public:
    IntIntMap();

	typedef IntIntMap_itr iterator;

	IntIntMap_itr begin(); 
	IntIntMap_itr end();

	int GetVal(int ikey);
	void SetVal(int ikey, int val);
	
	int ierr; //!< indicate error
};

class IntPtrMap_itr: public IntPtrMap_itr_parent
//! Iterator class for a map of integers to void* pointers
{
public:

	IntPtrMap_itr() {}
	IntPtrMap_itr( const IntPtrMap_itr_parent& ref);

	int   GetKey();
	void* GetVal();

    bool operator==(const IntPtrMap_itr& ref) const;
	bool operator!=(const IntPtrMap_itr& ref) const;	

	
	IntPtrMap_itr& next();
};

class IntPtrMap: public IntPtrMap_parent
//! Map of integers to void* pointers
{
public:
	IntPtrMap();

	typedef IntPtrMap_itr iterator;

	IntPtrMap_itr begin(); 
	IntPtrMap_itr end();
	
	void* GetVal(int ikey);
	void SetVal(int ikey, void* val);
	
	int ierr; //!< indicate error
};


class PtrIntMap_itr: public PtrIntMap_itr_parent
//! Iterator class of a map of pointers to integers 
{
public:
	PtrIntMap_itr() {}
	PtrIntMap_itr( const PtrIntMap_itr_parent& ref);

	void*  GetKey();
	int    GetVal();

    bool operator==(const PtrIntMap_itr& ref) const;
	bool operator!=(const PtrIntMap_itr& ref) const;	
	
	PtrIntMap_itr& next();
};


class PtrIntMap: public PtrIntMap_parent
//! Map of void* pointers to integers 
{
public:
    PtrIntMap();

	typedef PtrIntMap_itr iterator;

	PtrIntMap_itr begin(); 
	PtrIntMap_itr end();

	int GetVal(void* ptr);
	void SetVal(void* ptr, int val);
	
	int ierr; //!< indicate error
};

class PtrPtrMap_itr: public PtrPtrMap_itr_parent
//! Iterator class of a map of pointers to pointers
{
public:

	PtrPtrMap_itr() {}
	PtrPtrMap_itr( const PtrPtrMap_itr_parent& ref);

	void* GetKey();
	void* GetVal();

    bool operator==(const PtrPtrMap_itr& ref) const;
	bool operator!=(const PtrPtrMap_itr& ref) const;	
	
	PtrPtrMap_itr& next();
};

class PtrPtrMap: public PtrPtrMap_parent
//! Map of pointers to pointers
{
public:
	PtrPtrMap();

	typedef PtrPtrMap_itr iterator;

	PtrPtrMap_itr begin(); 
	PtrPtrMap_itr end();
	
	void* GetVal(void* ptr);
	void SetVal(void* ptr, void* val);
	int ierr; //!< indicate error
};


class PtrDoubleMap_itr: public PtrDoubleMap_itr_parent
//! Iterator class of a map of pointers to real numbers
{
public:
	PtrDoubleMap_itr() {}
	PtrDoubleMap_itr( const PtrDoubleMap_itr_parent& ref);

	void*  GetKey();
	double GetVal();

    bool operator==(const PtrDoubleMap_itr& ref) const;
	bool operator!=(const PtrDoubleMap_itr& ref) const;	

	PtrDoubleMap_itr& next();
};

class PtrDoubleMap: public PtrDoubleMap_parent
//! Map of pointers to real numbers
{
public:
    PtrDoubleMap();

	typedef PtrDoubleMap_itr iterator;

	PtrDoubleMap_itr begin(); 
	PtrDoubleMap_itr end();

	double GetVal(void* ptr);
	void SetVal(void* ptr, double val);
	
	int ierr; //!< indicate error
};


class PtrSet: public std::set<void*>
{
//! Set of pointers
public:
	PtrSet();
	
	int HasAtom(void* const ptr) const;
};

class VecPtr: public std::vector<void*>
{
//! Vector of pointers
public:
	VecPtr();
	VecPtr( int n);

	int Delete(const PtrSet& ptr_set); //!< Delete from the vector pointers that are contained in ptr_set; return the number of pointers removed
};

#endif  /* !HASTL_H */ 

