/*! \file hastl.cpp

    Typres derived from STL Template library

    \date 2003
    \author Igor Kurnikov 

*/

#include "hastl.h"

IntIntMap_itr::IntIntMap_itr( const IntIntMap_itr_parent& ref):IntIntMap_itr_parent(ref)
{

}


int
IntIntMap_itr::GetKey()
{
	return (*(*this)).first;
}

int
IntIntMap_itr::GetVal()
{
	return (*(*this)).second;
}

bool
IntIntMap_itr::operator==(const IntIntMap_itr& ref) const
{
	return (const IntIntMap_itr_parent&) (*this) == (const IntIntMap_itr_parent&) ref;
}

bool
IntIntMap_itr::operator!=(const IntIntMap_itr& ref) const
{
	return (const IntIntMap_itr_parent&) (*this) != (const IntIntMap_itr_parent&) ref;
}



IntIntMap_itr& 
IntIntMap_itr::next()
{
   (*this)++;
   return( *this );
}


IntIntMap::IntIntMap()
{

}

IntIntMap_itr
IntIntMap::begin()
{
	return IntIntMap_parent::begin();
}

IntIntMap_itr
IntIntMap::end()
{
	return IntIntMap_parent::end();
}

int
IntIntMap::GetVal(int ikey)
{
	IntIntMap_parent::iterator itr;
	itr = this->find(ikey);
	if(itr == ((IntIntMap_parent*)this)->end() ) 
	{
		ierr = 1;
		return 0;
	}
	
	ierr = 0;
	return (*itr).second;	
}

void
IntIntMap::SetVal(int ikey, int ival)
{
	(*this)[ikey] = ival;
}


IntPtrMap_itr::IntPtrMap_itr( const IntPtrMap_itr_parent& ref):IntPtrMap_itr_parent(ref)
{
	
}


int
IntPtrMap_itr::GetKey()
{
	return (*(*this)).first;
}

void*
IntPtrMap_itr::GetVal()
{
	return (*(*this)).second;
}

bool
IntPtrMap_itr::operator==(const IntPtrMap_itr& ref) const
{
	return (const IntPtrMap_itr_parent&) (*this) == (const IntPtrMap_itr_parent&) ref;
}

bool
IntPtrMap_itr::operator!=(const IntPtrMap_itr& ref) const
{
	return (const IntPtrMap_itr_parent&) (*this) != (const IntPtrMap_itr_parent&) ref;
}


IntPtrMap_itr& 
IntPtrMap_itr::next()
{
   (*this)++;
   return( *this );
}


IntPtrMap::IntPtrMap()
{

}

IntPtrMap_itr
IntPtrMap::begin()
{
	return IntPtrMap_parent::begin();
}

IntPtrMap_itr
IntPtrMap::end()
{
	return IntPtrMap_parent::end();
}

void*
IntPtrMap::GetVal(int ikey)
{
	IntPtrMap_parent::iterator itr;
	itr = this->find(ikey);
	if(itr == ((IntPtrMap_parent*)this)->end()) 
	{
		ierr = 1;
		return NULL;
	}
	
	ierr = 0;
	return (*itr).second;	
}

void
IntPtrMap::SetVal(int ikey, void* ptr)
{
	(*this)[ikey] = ptr;
}



PtrIntMap_itr::PtrIntMap_itr( const PtrIntMap_itr_parent& ref):PtrIntMap_itr_parent(ref)
{
	
}


void*
PtrIntMap_itr::GetKey()
{
	return (*(*this)).first;
}

int
PtrIntMap_itr::GetVal()
{
	return (*(*this)).second;
}

bool
PtrIntMap_itr::operator==(const PtrIntMap_itr& ref) const
{
	return (const PtrIntMap_itr_parent&) (*this) == (const PtrIntMap_itr_parent&) ref;
}

bool
PtrIntMap_itr::operator!=(const PtrIntMap_itr& ref) const
{
	return (const PtrIntMap_itr_parent&) (*this) != (const PtrIntMap_itr_parent&) ref;
}

PtrIntMap_itr& 
PtrIntMap_itr::next()
{
   (*this)++;
   return( *this );
}


PtrIntMap::PtrIntMap()
{

}

PtrIntMap_itr
PtrIntMap::begin()
{
	return PtrIntMap_parent::begin();
}

PtrIntMap_itr
PtrIntMap::end()
{
	return PtrIntMap_parent::end();
}

int
PtrIntMap::GetVal(void* ptr)
{
	PtrIntMap_parent::iterator itr;
	itr = this->find(ptr);
	if(itr == ((PtrIntMap_parent*)this)->end() ) 
	{
		ierr = 1;
		return 0;
	}
	
	ierr = 0;
	return (*itr).second;	
}

void
PtrIntMap::SetVal(void* ptr, int ival)
{
	(*this)[ptr] = ival;
}


PtrDoubleMap_itr::PtrDoubleMap_itr( const PtrDoubleMap_itr_parent& ref):PtrDoubleMap_itr_parent(ref)
{

}


void*
PtrDoubleMap_itr::GetKey()
{
	return (*(*this)).first;
}

double
PtrDoubleMap_itr::GetVal()
{
	return (*(*this)).second;
}

bool
PtrDoubleMap_itr::operator==(const PtrDoubleMap_itr& ref) const
{
	return (const PtrDoubleMap_itr_parent&) (*this) == (const PtrDoubleMap_itr_parent&) ref;
}

bool
PtrDoubleMap_itr::operator!=(const PtrDoubleMap_itr& ref) const
{
	return (const PtrDoubleMap_itr_parent&) (*this) != (const PtrDoubleMap_itr_parent&) ref;
}

PtrDoubleMap_itr& 
PtrDoubleMap_itr::next()
{
   (*this)++;
   return( *this );
}


PtrDoubleMap::PtrDoubleMap()
{

}

PtrDoubleMap_itr
PtrDoubleMap::begin()
{
	return PtrDoubleMap_parent::begin();
}

PtrDoubleMap_itr
PtrDoubleMap::end()
{
	return PtrDoubleMap_parent::end();
}


double
PtrDoubleMap::GetVal(void* ptr)
{
	PtrDoubleMap_parent::iterator itr;
	itr = this->find(ptr);
	if(itr == ((PtrDoubleMap*)this)->end() ) 
	{
		ierr = 1;
		return 0.0;
	}
	
	ierr = 0;
	return (*itr).second;	
}

void
PtrDoubleMap::SetVal(void* ptr, double val)
{
	(*this)[ptr] = val;
}

PtrPtrMap_itr::PtrPtrMap_itr( const PtrPtrMap_itr_parent& ref):PtrPtrMap_itr_parent(ref)
{

}


void*
PtrPtrMap_itr::GetKey()
{
	return (*(*this)).first;
}

void*
PtrPtrMap_itr::GetVal()
{
	return (*(*this)).second;
}


bool
PtrPtrMap_itr::operator==(const PtrPtrMap_itr& ref) const
{
	return (const PtrPtrMap_itr_parent&) (*this) == (const PtrPtrMap_itr_parent&) ref;
}

bool
PtrPtrMap_itr::operator!=(const PtrPtrMap_itr& ref) const
{
	return (const PtrPtrMap_itr_parent&) (*this) != (const PtrPtrMap_itr_parent&) ref;
}


PtrPtrMap_itr& 
PtrPtrMap_itr::next()
{
   (*this)++;
   return( *this );
}


PtrPtrMap::PtrPtrMap()
{

}


PtrPtrMap_itr
PtrPtrMap::begin()
{
	return PtrPtrMap_parent::begin();
}

PtrPtrMap_itr
PtrPtrMap::end()
{
	return PtrPtrMap_parent::end();
}


void*
PtrPtrMap::GetVal(void* ptr)
{
	PtrPtrMap_parent::iterator itr;
	itr = this->find(ptr);
	if( itr == PtrPtrMap_parent::end() ) 
	{
		ierr = 1;
		return NULL;
	}
	
	ierr = 0;
	return (*itr).second;	
}

void
PtrPtrMap::SetVal(void* ptr, void* val)
{
	(*this)[ptr] = val;
}


PtrSet::PtrSet()
{

}

int
PtrSet::IsMember( void* const ptr) const
{
	if( this->find(ptr) != this->end()) return 1;

	return 0;
}

VecPtr::VecPtr()
{

}

VecPtr::VecPtr(int n): vector<void*>(n)
{

}


int 
VecPtr::Delete(const PtrSet& ptr_set)
{
	int n = this->size();
	VecPtr copy_vec;
	copy_vec.reserve(n);
	
	int i;
	for(i = 0; i < n; i++)
	{
		if( !ptr_set.IsMember( (*this)[i]) )
		{
			copy_vec.push_back( (*this)[i]);
		}
	}
	this->swap(copy_vec);

	return (n - this->size());
}

