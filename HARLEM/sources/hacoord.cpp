/*! \file hacoord.cpp
 
    Basic Classes for coordinates in HARLEM 
    
    \author Igor Kurnikov  
    \date 2013- 
  
*/

#include <assert.h>
#include "hacoord.h"

harlem::HashID harlem::hash(const std::string& str)
{
	assert( sizeof(harlem::HashID)==8 );

	HashID sum=0;

	int len = std::min( (int)str.length(), 12);

	for(int i=len-1; i>=0; --i)
	{
		int c = (int)str[i];

		if( islower(c) )
		{
			c = c - 'a';
		}
		else if( isupper(c) )
		{
			c = tolower(c) - 'a';
		}
		else if( isdigit(c) )
		{
			c -= '0';
			c += 26;
		}
		else if( c=='_' )
		{
			c = 36;
		}
		else
		{
			continue;
		}

		sum = 40*sum + c;
	}

	return sum;
}

std::string harlem::unhash(const harlem::HashID& id)
{
	harlem::HashID code = id;
	std::string str;
	while(code != 0LL)
	{
		harlem::HashID next = code/40LL;
		harlem::HashID c = code - 40LL*next;
		assert(c < 40LL);
		if( c < 26LL )
		{
			c = c + 'a';
		}
		else if( c < 36LL )
		{
			c = c - 26LL + '0';
		}
		else
		{
			assert( c==36LL );
			c = '_';
		}

		str.append(1, (char)c);

		code = next;
	}

	return str;
}

namespace harlem {

HashMap::HashMap() 
{ 

}
		
HashMap::HashMap( const HashMap& ref ) 
{ 
	Copy(ref); 
}
	
HashMap::~HashMap() 
{

}

void HashMap::Copy( const HashMap& ref )
{
	this->i_params = ref.i_params;
	this->d_params = ref.d_params;
	this->s_params = ref.s_params;
	this->a_params = ref.a_params;
}
		
HashMap* HashMap::clone() const
{
	HashMap* popt = new HashMap();
	popt->Copy(*this);
	return popt;
}

void HashMap::set_i( const HashID& parmid, const int& value )
{
	i_params[parmid] = value;
}

void  HashMap::set_d( const HashID& parmid, const double& value )
{
	d_params[parmid] = value;
}

void  HashMap::set_s( const HashID& parmid, const std::string& value )
{
	s_params[parmid] = value;
}

void HashMap::set_a( const HashID& parmid, const boost::any& value )
{
	a_params[parmid] = value;
}

void  HashMap::set_i( const std::string& parmname, const int& value )
{
	 this->set_i( harlem::hash(parmname), value );
}

void  HashMap::set_d( const std::string& parmname, const double& value )
{
	this->set_d( harlem::hash(parmname), value );
}

void  HashMap::set_s( const std::string& parmname, const std::string& value )
{
	this->set_s( harlem::hash(parmname), value );
}

void HashMap::set_a( const std::string& parmname, const boost::any& value )
{
	this->set_a( harlem::hash(parmname), value );
}

int HashMap::get_i(const HashID& parmid) const
{
	std::map< HashID, int >::const_iterator itr = i_params.find(parmid);
	if( itr != i_params.end() ) return (*itr).second;
	return 0;
}
		
double HashMap::get_d(const HashID& parmid) const
{
	std::map< HashID, double >::const_iterator itr = d_params.find(parmid);
	if( itr != d_params.end() ) return (*itr).second;
	return 0.0;
}
		
std::string HashMap::get_s(const HashID& parmid) const
{
	std::map< HashID, std::string >::const_iterator itr = s_params.find(parmid);
	if( itr != s_params.end() ) return (*itr).second;
	return "";
}

boost::any HashMap::get_a(const HashID& parmid) const
{
	std::map< HashID, boost::any >::const_iterator itr = a_params.find(parmid);
	if( itr != a_params.end() ) return (*itr).second;
	return (int)0;
}

int HashMap::get_i(const std::string& parmname) const
{
	return this->get_i( hash(parmname) );
}
		
double HashMap::get_d(const std::string& parmname) const
{
	return this->get_d( hash(parmname) );
}

std::string HashMap::get_s(const std::string& parmname) const
{
	return this->get_s( hash(parmname) );
}

boost::any HashMap::get_a(const std::string& parmname) const
{
	return this->get_a( hash(parmname) );
}

bool HashMap::has_i(const HashID& parmid) const
{
	return ( i_params.count(parmid) > 0 );
}

bool HashMap::has_d(const HashID& parmid) const
{
	return ( d_params.count(parmid) > 0 );
}

bool HashMap::has_s(const HashID& parmid) const
{
	return ( s_params.count(parmid) > 0 );
}

bool HashMap::has_a(const HashID& parmid) const
{
	return ( a_params.count(parmid) > 0 );
}

bool HashMap::has_i(const std::string& parmname) const
{
	return this->has_i( hash(parmname));
}
		
bool HashMap::has_d(const std::string& parmname) const
{
	return this->has_d( hash(parmname));
}
		
bool HashMap::has_s(const std::string& parmname) const
{
	return this->has_s( hash(parmname));
}

bool HashMap::has_a(const std::string& parmname) const
{
	return this->has_a( hash(parmname));
}

}