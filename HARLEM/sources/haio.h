/*! \file haio.h

   Include Files for standard C++ and C IO libraries

   \author  Igor Kurnikov 
   \date 1997-2002

*/
#ifndef HAIO_H
#define HAIO_H

// Use new C++ standard header for iostream  
#ifdef _MSC_VER
	#pragma warning (disable:4786)
	#define DllImport   __declspec( dllimport ) 
	#define DllExport   __declspec( dllexport )

#else
	#define DllExport 
	#define DllImport  extern
#endif

#include <stdio.h>
#include <iostream>
#include <fstream>
#if defined(_MSC_VER)
#include <sstream>
#endif
#include <strstream>
#include <iomanip>

using namespace std;

#if defined(_MSC_VER)
  
  typedef struct DIR DIR;

  struct dirent
  {
      char *d_name; //!< file name
  };

  DIR           *opendir(const char *); //!< opens the directory specified by name
  int           closedir(DIR *);        //!< closes the directory stream associated with dir
  struct dirent *readdir(DIR *);        //!< function is used to iterate through the directory stream dir
  void          rewinddir(DIR *);       //!< function can be used to reset the directory stream dir to the start. 
#else
#include <dirent.h>
#endif 

#if defined(_MSC_VER)
const std::string path_sep = "\\";
#else
const std::string path_sep = "/";
#endif

// typedef ostream ostream_withassign;
//Mikola temporal fix
#ifndef SWIG
extern bool find_line_in_file(FILE* fp, const char* str_comp, char* cur_line, const int len, const bool rew= false);
#endif

const std::string FLOAT_E16_8= "%16.8e";
const std::string FLOAT_F12_7= "%12.7f";
const std::string FLOAT_F15_7= "%15.7f";
const std::string FLOAT_F8_3 = "%15.7f";

//const std::string FLOAT_E16_8= "(E16.8)";
//const std::string FLOAT_F12_7= "(F12.7)";
//const std::string FLOAT_F15_7= "(F15.7)";
//const std::string FLOAT_F8_3 = "(F8.3)";

extern "C" 
{
extern int ErrorMessage(const char* str); //!< Print error message
extern int PrintMessage(const char* str); //!< Print Message to a standard output and control window on HARLEM panel
extern int ha_copy_file(const char* src, const char* tgt, const int mode = 0); //!< copy file named src to file named tgt
extern int ha_delete_file( const char* fname);  //!< delete file with the name fname
extern void RedirectIOToConsole();              //!< Set STDOUT (standart output channel) to Console
extern int RedirectIOToFile(const char* fname); //!< Redirect STDOUT (standart output channel) to file
extern int RedirectIOToMultipleFilesMPI(const char* fname); //!< Redirect STDOUT (standart output channel) to multiple files useful then run harlem in parallel
extern void write_log_(const char* str, int n); //!< Fortran function to write string to STDOUT 

extern int RestoreIOToConsole();                //!< Restore STDOUT to console that been redirected by RedirectIOToFile()


#if !defined(HAIO_CPP)
extern int PrintLog(const char* str, ... );  //!< HARLEM standard function to print info into STDOUT (syntax as in printf) 
extern int ErrorInMod(const char* module, const char* msg);  //!< HARLEM standard function to print info about error in function
#else
	int PrintLog(const char* str, ... );  //!< HARLEM standard function to print info into STDOUT (syntax as in printf) 
	int ErrorInMod(const char* module, const char* msg);  //!< HARLEM standard function to print info about error in function
#endif

}

#endif /* !HAIO_H */


