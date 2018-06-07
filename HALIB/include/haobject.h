/*!  \file haobject.h

   Basic Classes to define objects in HARLEM
   
   \author Nikolay Simakov, Igor Kurnikov  
   \date 2006-

*/
#ifndef HAOBJECT_H
#define HAOBJECT_H

#include <string>
#include <vector>
#include <map>

class TiXmlElement;
class HaContext;

class HaObject
//! Basic object in HARLEM serves as a parent class for many objects in the program  
//! Provides functionality of saving to and reading from XML elements 
//! Returns its Name identificator 
//!
{
public:
	HaObject();
	virtual ~HaObject();
	virtual int Clear();
	virtual int SetDefaultValues(HaContext* p_ctxt = NULL);
	virtual int SaveXML(TiXmlElement* Elt, HaContext* p_ctxt = NULL );
	virtual int LoadXML(const TiXmlElement* Elt, HaContext* p_ctxt = NULL );

	void SetName(const std::string& name_new );
	const char* GetName() const; 
	const char* GetCStrName() const;
	std::string GetStdStrName() const;

protected:
	std::string name;
};

class HaContext
//!< Base Class to provide a context for creating and setting objects
{
public:
	HaContext();
	virtual ~HaContext();

	int  GetPrimaryOption() const;
	void SetPrimaryOption(int option_primary_new );

private:
	int option_primary;
};



#endif
