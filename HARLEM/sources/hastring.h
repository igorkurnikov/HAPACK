/*! \file hastring.h
 
    Character String class in HARLEM 
    use STL library if supported
 
    \author Igor Kurnikov 
    \date 1997-2009
*/
#ifndef HASTRING_H
#define HASTRING_H

#if defined(_MSC_VER)
#pragma warning (disable:4786)
#endif

// Use STL string
#include <string>

#if defined(_MSC_VER)

#if !defined(DllImport)
	#define DllImport   __declspec( dllimport )
#endif
#if !defined(DllExport)
	#define DllExport   __declspec( dllexport )
#endif
#endif 

#include "hastl.h"

class StrVec;

//class HaString: public std::string
////! Character string object in HARLEM based on STL string
//{
//public:
//	HaString() {}
//	HaString(const string& bas_str): std::string(bas_str) {}
//	HaString(const HaString& str): std::string(str) {}
//	HaString(const char* s): std::string(s) {}
//	HaString(const char* s, int n): std::string(s,n) {}
//	virtual ~HaString() {}
//
//	const char* trim();  //!< truncate spaces before and after the string 
//	HaString& to_upper();      //!< Converts all characters to upper case and returns the result. 
//	HaString& to_lower();      //!< Converts all characters to lower case and returns the result.
//
////	operator const char*() { return c_str(); }
//    
//	int SplitTo(StrVec& str_vec, const char *strDelimit= ""); 
//};


namespace harlem
{
	std::string GetDirFromFileName(const std::string& fname);     //!< Extract directory name from the full file name
	std::string GetPrefixFromFullName(const std::string& fname);  //!< Get prefix from the file name 
	std::string GetExtFromFileName(const std::string& fname);     //!< Get extention from the file name

	bool IsFloat(const std::string& str);    //!< Check if the string represent floating point number (return false if integer ) 
	bool IsInt  (const std::string& str);    //!< Check if the string integer number
	std::string ToString(int ival);   //!< Convert integer value to String
	double ToDouble(std::string str); //!< Convert string to double  

	std::string StdXMLHeader();      //!< first line of XML files created by HARLEM 
	std::string HarlemDataHeader();  //!< main element (HARLEM_DATA) header of XML files ceated by HARLEM
	std::string HarlemDataFooter();  //!< main element (HARLEM_DATA) footer of XML files ceated by HARLEM

};

class StrDoubleMap;

class StrDoubleMap_itr
//! Iterator class of a map of strings to real numbers
{
public:
	StrDoubleMap_itr( StrDoubleMap& new_map);
    StrDoubleMap_itr( StrDoubleMap_itr& ref);

	std::string GetKey();
	double   GetVal();
	int GetFirst();
	int GetNext();

protected:
	map<std::string, double, less<std::string> >::iterator itr;
	StrDoubleMap& int_map;
};

class StrDoubleMap: public map<std::string, double, less<std::string> >
//! Map of strings to real numbers
{
public:

	typedef StrDoubleMap_itr iterator;

	void clear();
//	int count(const char* str);
	int size();

	double GetVal(const char* str);
	void SetVal(const char* str, double val);
	int ierr; //!< indicate error
};

class StrStrMap;

class StrStrMap_itr
//! Iterator class of a map of strings to strings
{
public:
	StrStrMap_itr( StrStrMap& new_map);
	StrStrMap_itr( StrStrMap_itr& ref);
	
	const char* GetKey();
	const char* GetVal();
	int GetFirst();
	int GetNext();

protected:
	map<std::string, std::string, less<std::string> >::iterator itr;
	StrStrMap& int_map;
};

class StrStrMap: public map<std::string, std::string, less<std::string> >
//! Map of strings to strings
{
public:

	typedef StrStrMap_itr iterator;

	void clear();
	int count(const char* str);
	int size();

	const char* GetVal(const char* str);
	void SetVal(const char* str, const char* val);
	int ierr; //!< indicate error

};

class StrIntMap: public map<std::string, int, less<std::string> >
//! Map of strings to integer numbers
{
public:

#if defined(SWIG)
	void clear();
	int size();
#endif

	int count(const char* str);

	int GetVal(const char* str);
	void SetVal(const char* str, int val);
	int ierr; //!< indicate error
};

class StrVec: public std::vector<std::string>
//! class for a vector of strings 
{
public:
	StrVec() {}
	StrVec(size_t n): std::vector<std::string>(n) {}
	StrVec(const std::vector<std::string>& str_vec): std::vector<std::string>(str_vec) {}

#if defined(SWIG) 	
  void reserve(size_t n); 
  void resize(size_t n, const char* str);
  size_t size();
  void push_back(const char* str);
#endif  
	
    const char* GetVal(size_t idx);
	void SetVal(size_t idx, const char* val);
    
};

typedef map<std::string,void*,less<std::string> > StrPtrMap;
typedef map<int,std::string,less<int> > IntStrMap;

#if defined(linux) || defined(__DECCXX)
extern int _fstrnicmp(const char *s1, const char *s2, size_t n);
#endif

extern int stricmp_loc( const char* str1, const char* str2); //!< low case comparison of two strings 
                                                         
extern int strcmp_trunc(const char* str1, const char* str2); //!< compare two strings excluding all white spaces

extern int strcpy_to_fort(char* str_fort, const char* c_str, int len_fort); //!< Fill fortran string (with length len_fort) from c_str                                                              

extern int stricmp_trunc(const char* str1, const char* str2); //!< compare two strings excluding white spaces and converting to lower case

#endif /* !HASTRING_H */
