/*! \file gaufile.cpp

    Definitions of Classes to interact with Gaussian RWF Files

    \author Igor Kurnikov 

    \date 1997-2003
*/

#include <stdlib.h>

#include "haconst.h"
#include "haio.h"
#include "gaufile.h"

#include <string.h>

int
GauFile::set_open_mode(const char* mode)
{
  // Set mode to open Gaussian file 
  if(!strcmp(mode,"new"))
    open_mode=1;
  else if(!strcmp(mode,"old"))
    open_mode=2;
  else if(!strcmp(mode,"unknown"))
    open_mode=3;
  else if(!strcmp(mode,"ro_shared"))
    open_mode=4;
  else if(!strcmp(mode,"ro_shared_nowrite"))
    open_mode=5;
  else
    {
      std::cerr << "Gaufile::set_open_mode : unknown Gaussian file open mode: ";
      std::cerr << mode << std::endl;
      std::cerr << "set open mode to UNKNOWN " << std::endl;
    open_mode=3;
	return 1;
    }
	return 0;
}

GauFile::GauFile()
{
  fname= "";
  open_mode=3;
  iunit=1;
  init_alloc=0;
  flag_open=0;
}

GauFile::GauFile(const char* NewFname, const int NewIunit, const char* mode)
{
  open_mode=3;
  init_alloc=0;
  flag_open=0;  
  set_file_name(NewFname);
  set_gau_iunit(NewIunit);
  set_open_mode(mode); 
}

int 
GauFile::open()
{
#if defined(GAUSSVER)
	
  // Call Gaussian fopen subroutine for a file
  //    iunit = fileio/ntran unit to open.
  //   open_mode  1=new.
  //              2=old.
  //              3=unknown.
  //              4=old, readonly and shared.
  //              5=old, readonly, others can't write (nyi).
  //              negative values prohibit sharing (nyi).
  //  fname   ... name of the file.
  //  ialoci  ... requested initial allocation in wp words.
  //  ialoco  ... returned initial allocation in wp words.
  int fort_iunit=(int)iunit;
  int imode=(int)open_mode;
  int init_alloc_ask=(int)init_alloc;
  int init_alloc_get;
  char fname_for[1024]={1024*0};
  
  for(int i=0; i < fname.length(); i++) 
    fname_for[i]=fname[i];
  
  fname_for[fname.length()]=0;
  
  fopen_(&fort_iunit, &imode, fname_for, &init_alloc_ask, 
         &init_alloc_get,(int) 1024);
  if((int)init_alloc_get)
    {
      init_alloc=(int)init_alloc_get;
      flag_open=1;
      return 1;
    }
  else
    return 0;
#else
	PrintLog("Error in GauFile::open() \n");
	PrintLog("GAUSSIAN Library is not present \n");
	return FALSE;
#endif 

}

int 
GauFile::close(char* mode)
{
#if defined(GAUSSVER)
  int idisp;
  int fort_iunit=(int)iunit;

  if(!strcmp(mode,"keep"))
    idisp=0;
  else if(!strcmp(mode,"delete"))
    idisp=1;
  else
    {
    cerr << "error in Gaufile::close " << endl;
    cerr << "Unknown close mode: " << mode << endl;
    return 1;
   }
  fclose_(&fort_iunit,&idisp);

#else
	PrintLog("Error in GauFile::close() \n");
	PrintLog("GAUSSIAN Library is not present \n");
	return FALSE;
#endif 
  return 0;
}

void 
GauFile::dump_info()
{
#if defined(GAUSSVER)
  // Dump summary starage info on open Gaussian files
  fdump_();
#else
	PrintLog("Error in GauFile::dump_info() \n");
	PrintLog("GAUSSIAN Library is not present \n");
#endif 
}

int 
GauFile::fileio(const char* operation, const int sub_file, 
              int length, void* target, int position)
{
#if defined(GAUSSVER)
  int len   = (int)length;
  int isubf;
  int ipos  = (int)position;
  int ioper;
  int len_file;

  if(!strcmp(operation,"write"))
    ioper=1;
  else if(!strcmp(operation,"read"))
    ioper=2;
  else if(!strcmp(operation,"async_write"))
    ioper=-1;
  else if(!strcmp(operation,"async_read"))
    ioper=-2;
  else
    {
      cerr << "error in GauFile::fileio :" << endl;
      cerr << "unknown IO operation: " << operation << endl;
      return False;
    }
  if(sub_file == 0)
    {
      cerr << "error in GauFile::fileio :" << endl;
      cerr << "subfile number is equal to 0 " << endl;
      return False;
    } 

  int isign=sub_file/abs(sub_file);
 
  if(abs(sub_file) > 10000)
    isubf=(int)sub_file;
  else
    {
      isubf=(int)(isign*iunit*10000 + sub_file);
    }
  if(abs(ioper) == 2)
  {
	len_file=itqry_(&isubf);
	if(len_file < len) return False;
  }
	fileio_(&ioper,&isubf,&len,target,&ipos);
	return True;
#else
	PrintLog("Error in GauFile::fileio() \n");
	PrintLog("GAUSSIAN Library is not present \n");
	return FALSE;
#endif 
}






