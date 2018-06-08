/*! \file haobject.cpp
 *  
 *   Base Class for objects in HARLEM
 *
 *   \author Nikolay Simakov                                                       *
 *   nsimakov@andrew.cmu.edu                                               *
 *   Maria Kurnikova Research Group                                        *
 *   http://crete.chem.cmu.edu/                                            *
 *   Carnegie Mellon University, 
 *   \date 2006-                                      *
 ***************************************************************************/

#include "haobject.h"

HaObject::HaObject()
{
	name="";
}
    
HaObject::~HaObject()
{

}

int HaObject::Clear()
{ 
	return 0;
}

int HaObject::SetDefaultValues( HaContext* p_ctxt )
{ 
	return 0;
}

int HaObject::SaveXML(TiXmlElement* Elt, HaContext* p_ctxt )
{ 
	return 0;
}
    
int HaObject::LoadXML(const TiXmlElement* Elt, HaContext* p_ctxt )
{ 
	return 0;
}
    
void HaObject::SetName(const std::string& name_new)
{
	name = name_new;
}
	
const char* HaObject::GetName() const 
{ 
	return name.c_str(); 
}

const char* HaObject::GetCStrName() const
{
	return name.c_str();
}

std::string HaObject::GetStdStrName() const
{
	return name;
}

HaContext::HaContext()
{
	option_primary = 0;
}

HaContext::~HaContext()
{

}

int HaContext::GetPrimaryOption() const
{
	return option_primary;
}


void HaContext::SetPrimaryOption(int option_primary_new )
{
	option_primary = option_primary_new;
}
