/*! \file hastring.cpp 

  Classes and functions for string operations in HARLEM
 
   \author Igor Kurnikov
 
   \date 1998-2003 
*/
#define HASTRING_CPP

#include "haconst.h"
#include "hastring.h"
#include "haio.h"

#include <string>
#include <stdlib.h>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


//  already defined in DO_LIB
//char* my_sdup(const char* str)
//{
//	char* s = (char*) malloc(strlen(str)+1);
//	strcpy(s,str);
//	return s;
//}

int strcmp_trunc(const char* str1, const char* str2)
{
	int i,len;
	std::string str1_trunc, str2_trunc;
	len=strlen(str1);
	for(i=0; i< len; i++)
	{
		if( !isspace(str1[i])) str1_trunc+=str1[i]; 
	}
	len=strlen(str2);
	for(i=0; i< len; i++)
	{
		if( !isspace(str2[i])) str2_trunc+=str2[i]; 
	}
	return( strcmp(str1_trunc.c_str(), str2_trunc.c_str()) );
} 

int stricmp_trunc(const char* str1, const char* str2)
{
	int i,len;

	std::string str1_trunc, str2_trunc;
	len=strlen(str1);

	for(i=0; i< len; i++)
	{
		if( !isspace(str1[i])) str1_trunc+=str1[i]; 
	}

	len=strlen(str2);

	for(i=0; i< len; i++)
	{
		if( !isspace(str2[i])) str2_trunc+=str2[i]; 
	}
	return( stricmp_loc(str1_trunc.c_str(), str2_trunc.c_str()) );
}


//const char* HaString::trim()
//{
//	int i;
//	int ibeg= size(); int iend= -1;
//	for(i=0; i < size(); i++)
//	{
//		if(!isspace( (*this)[i])) 
//		{
//			ibeg=i;
//			break;
//		}
//	}
//	for(i= size()-1; i >=0 ; i--)
//	{
//		if(!isspace( (*this)[i])) 
//		{
//			iend=i;
//			break;
//		}
//	}
//
//	if(ibeg > iend) 
//	{
//		erase();
//	}
//	else
//	{
//		for(i=ibeg; i <= iend; i++)
//		{
//			(*this)[i-ibeg]= (*this)[i];
//		}
//		(*this).resize(iend-ibeg+1);
//	}
//	return (*this).c_str();
//}


int stricmp_loc(const char* str1, const char* str2)
{
#if defined(_MSC_VER)
	return _stricmp(str1,str2);
#else
	int len1= strlen(str1);
	int len2= strlen(str2);
	if(len1 != len2)
		return( (len1 > len2)? 1 : -1);
	for(int i=0; i < len1; i++)
	{
		int k1=  tolower(str1[i]);
		int k2=  tolower(str2[i]);
		if( k1 != k2)
			return ( (k1 > k2)? 1 : -1 );
	}
	return (0);
#endif
}

//HaString& HaString::to_upper()
//{
//	int len=size();
//	for(int i=0; i < len; i++)
//		(*this)[i]=toupper( (*this)[i]);
//	return *this;
//}
//
//HaString& HaString::to_lower()
//{
//	int len=size();
//	for(int i=0; i < len; i++)
//		(*this)[i]=tolower( (*this)[i]);
//	return *this;
//}

//int HaString::SplitTo(StrVec& str_vec, const char *strDelimit)
//{
//	str_vec.resize(0);
//	int len = this->size();
//	char* str = (char*)malloc(len+1);
//    char* token = strtok( str, strDelimit );
//    while( token != NULL )
//    {
//       std::string str_token = token;
//	   str_vec.push_back(token);
//       token = strtok( NULL, strDelimit);
//    }
//	free(str);
//	return TRUE;
//}

#if defined(linux) || defined (__DECCXX)

int _fstrnicmp(const char *s1, const char *s2, size_t n)
{
  return( strncasecmp(s1,s2,n));
} 

#endif

int strcpy_to_fort(char* str_fort, const char* c_str, int len_fort)
{
	int ires;
	int len_c = strlen(c_str);
	int i;
	if( len_c > len_fort) len_c = len_fort;
	for(i =0 ; i < len_c; i++)
	{
		str_fort[i] = c_str[i];
	}
	if(len_c < len_fort)
	{
		ires = len_c;
		for(int i= len_c; i < len_fort; i++)
		{
			str_fort[i] = ' ';
		}
	}
	else
	{
		ires = len_fort;
	}
	return ires;
}

StrDoubleMap_itr::StrDoubleMap_itr(StrDoubleMap& new_map): int_map(new_map)
{
	itr = int_map.begin();
}

StrDoubleMap_itr::StrDoubleMap_itr( StrDoubleMap_itr& ref):
int_map(ref.int_map),itr(ref.itr)
{
   	
}

std::string StrDoubleMap_itr::GetKey()
{
	return (*itr).first;
}

double StrDoubleMap_itr::GetVal()
{
	return (*itr).second;
}

int StrDoubleMap_itr::GetFirst()
{
   itr = int_map.begin();
   return( itr != int_map.end());
}

int StrDoubleMap_itr::GetNext()
{
   itr++;
   return( itr != int_map.end());
}

void StrDoubleMap::clear()
{
	std::map<std::string, double>::clear();
}

//int StrDoubleMap::count(const char* str)
//{
//	return map<std::string, double, less<std::string> >::count(str);
//}

int StrDoubleMap::size()
{
	return std::map<std::string, double>::size();
}

double StrDoubleMap::GetVal(const char* str)
{
	map<std::string, double>::iterator itr;
	itr = this->find(str);
	if(this->end() == itr) 
	{
		ierr = 1;
		return 0.0;
	}
	
	ierr = 0;
	return (*itr).second;	
}

void StrDoubleMap::SetVal(const char* str, double val)
{
	(*this)[str] = val;
}


StrStrMap_itr::StrStrMap_itr(StrStrMap& new_map): int_map(new_map)
{
	itr = int_map.begin();
}

StrStrMap_itr::StrStrMap_itr( StrStrMap_itr& ref):
int_map(ref.int_map),itr(ref.itr)
{
   	
}

const char* StrStrMap_itr::GetKey()
{
	return (*itr).first.c_str();
}

const char* StrStrMap_itr::GetVal()
{
	return (*itr).second.c_str();
}

int StrStrMap_itr::GetFirst()
{
   itr = int_map.begin();
   return( itr != int_map.end());
}

int StrStrMap_itr::GetNext()
{
   itr++;
   return( itr != int_map.end());
}


static std::string empty_str;

std::string StrStrMap::GetVal(const std::string& key)
{
	std::map<std::string, std::string>::iterator itr;
	itr = this->find(key);
	if(this->end() == itr) 
	{
		ierr = 1;
		return empty_str;
	}
	
	ierr = 0;
	return (*itr).second;
}

void StrStrMap::SetVal(const std::string& key, const std::string& value)
{
	(*this)[key] = value;
}

int StrIntMap::count(const char* str)
{
    return std::map<std::string, int>::count((std::string)str);
}

int StrIntMap::GetVal(const char* str)
{
	std::map<std::string, int>::iterator itr;
	itr = this->find(str);
	if(this->end() == itr) 
	{
		ierr = 1;
		return -9999999;
	}
	
	ierr = 0;
	return (*itr).second;
}

void StrIntMap::SetVal(const char* str, int val)
{
	(*this)[str] = val;
}


const char* StrVec::GetVal(size_t idx)
{
	return ((*this)[idx]).c_str();
}

void StrVec::SetVal(size_t idx, const char* val)
{
	(*this)[idx] = val;
}
 
std::string harlem::GetDirFromFileName(const std::string& fname)
{
	std::string dir_name = fname;
   
	size_t ip1 = fname.find_last_of('/');
	size_t ip2 = fname.find_last_of('\\');

	size_t ip = MaxFun(ip1,ip2);
	if( ip != -1) return fname.substr(0,ip);
	
	return "";
}

 
std::string harlem::GetPrefixFromFullName(const std::string& fname)
{
	std::string prefix = fname;

	size_t ip1 = prefix.find_last_of('/');
	size_t ip2 = prefix.find_last_of('\\');

	size_t ip;
	if( ip1 == std::string::npos ) ip = ip2;
	else if( ip2 == std::string::npos ) ip = ip1;
	else ip = MaxFun(ip1,ip2);

	if( ip != std::string::npos ) 
	{
		prefix = prefix.substr(ip+1);
	}

	ip = prefix.find_last_of('.');
	if( ip != std::string::npos )
	{
		prefix = prefix.substr(0,ip);
	}
	return ( prefix );
}
 
std::string harlem::GetExtFromFileName(const std::string& fname)
{	
	size_t ip = fname.find_last_of('.');

	if( ip != -1) return fname.substr(ip+1); 
	
	return "";
}

bool harlem::IsFloat(const std::string& str)
{
	std::istringstream is(str);
	double dval;
	is >> dval;
	if( is.fail() ) return false;
	if( IsInt(str) ) return false;
	return true;

//	std::string::const_iterator citr = str.begin();
//	for( ; citr != str.end(); citr++)
//	{
//		if( !isdigit(*citr) && (*citr != '-') && (*citr != ' ') && ) return false;
//	}
//	return true;


//	boost::trim(str);
//	try
//	{
//		boost::lexical_cast<float>(str);
//	}
//	catch (boost::bad_lexical_cast &)
//	{
//		return false;
//	}
//	if( harlem::IsInt(str) ) return false;
//
//	return true;	
}

bool harlem::IsInt(const std::string& str)
{
	std::string::const_iterator citr = str.begin();
	for( ; citr != str.end(); citr++)
	{
		if( !isdigit(*citr) && (*citr != '-') && (*citr != ' ') ) return false;
	}
	return true;
}

std::string harlem::ToString(int ival)
{
	std::stringstream os;
	os << ival;
	return os.str();
}

double harlem::ToDouble(std::string str)
{
	double dval = 0.0;
	std::istringstream is(str);
	is >> dval;
	return dval;
}

std::string harlem::StdXMLHeader()
{
	return "<?xml version=\"1.0\" encoding=\"ISO-8859-15\" ?>";
}

std::string harlem::HarlemDataHeader()
{
	return "<HARLEM_DATA version=\"01/07/2013\">";
}

std::string harlem::HarlemDataFooter()
{
	return "</HARLEM_DATA>";
}

#if defined(__INTEL_COMPILER) && defined(__x86_64)

extern const char* _intel_sse2_strchr ( const char * str, int character );
extern            char* _intel_sse2_strchr (            char * str, int character );

const char* _intel_sse2_strchr ( const char * str, int character )
{
	return strchr( str, character ); 
}

	  char* _intel_sse2_strchr (              char * str, int character )
{
	return strchr( str, character ); 
}
#endif