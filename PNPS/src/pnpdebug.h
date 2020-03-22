//
// C++ Interface: debug
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PNPDEBUG
#define PNPDEBUG
#include <string>
#include <stdexcept>
#include <stdarg.h>
#include <stdio.h>
#include <time.h>
#if !defined(_MSC_VER) 
#include <sys/time.h>
#endif
#ifdef USE_WXLOG
#include <wx/wx.h>
#endif


//static char pnptmpdir[PNP_MAP_IO_STRING_LENGTH];
#define HA_EXIT_SUCCESS 1

#define PNPS_EXIT_FAILURE 0
#define PNPS_EXIT_SUCCESS 1

#define DeleteCArray(CArr) if(CArr!=NULL){delete [] CArr;CArr=NULL;}
#define DeleteCVecArray(CVArr,Nelem) if(CVArr!=NULL){int i;for(i=0;i<Nelem;i++){if(CVArr[i]!=NULL){delete [] CVArr[i];CVArr[i]=NULL;}}}
#define DeleteObjByPnt(CArr) if(CArr!=NULL){delete CArr;CArr=NULL;}

#define PNP_EXIT_FAIL(Message) pnpThrowErrorRaw(__LINE__,__FILE__,Message);
#define PNP_EXIT_FAIL_NULL(Pointer,Message) if(Pointer==NULL){pnpThrowErrorRaw(__LINE__,__FILE__,Message);}
#define PNP_EXIT_NULL_ON_NULL(Pointer,Message) if(Pointer==NULL){pnpThrowErrorRaw(__LINE__,__FILE__,Message);}
#define PNP_EXIT_ON_FAIL(Status) if(Status==EXIT_FAILURE){pnpThrowErrorRaw(__LINE__,__FILE__,"Something failed");}
#define PNP_EXIT_ON_FAIL_MES(Status,Message) if(Status==EXIT_FAILURE){pnpThrowErrorRaw(__LINE__,__FILE__,Message);}
#define PNP_EXIT_ON_FAIL(Status) if(Status==EXIT_FAILURE){pnpThrowErrorRaw(__LINE__,__FILE__,"Something failed");}
#define PNP_EXIT_FAIL_NULL1(Pointer,Message,Var1) if(Pointer==NULL){pnpThrowErrorRaw(__LINE__,__FILE__,Message,Var1);}

#if 0
#define PNP_EXIT_FAIL(Message) {pnpError(Message);return EXIT_FAILURE;}
#define PNP_EXIT_FAIL_NULL(Pointer,Message) if(Pointer==NULL){pnpError(Message);return EXIT_FAILURE;}
#define PNP_EXIT_NULL_ON_NULL(Pointer,Message) if(Pointer==NULL){pnpError(Message);return NULL;}
#define PNP_EXIT_ON_FAIL(Status) if(Status==EXIT_FAILURE){return EXIT_FAILURE;}
#define PNP_EXIT_ON_FAIL_MES(Status,Message) if(Status==EXIT_FAILURE){pnpError(Message);return EXIT_FAILURE;}
#define PNP_EXIT_ON_FAIL(Status) if(Status==EXIT_FAILURE){return EXIT_FAILURE;}
#define PNP_EXIT_FAIL_NULL1(Pointer,Message,Var1) if(Pointer==NULL){pnpError(Message,Var1);return EXIT_FAILURE;}
#endif



extern int TotalPNPWarningMassage;
extern int TotalPNPErrorMassage;

inline void InfoPrint(const char* str,...)
{
  va_list arg_list;
  va_start(arg_list,str);
  printf("Info: ");
  vprintf(str,arg_list);
  va_end(arg_list);
}
inline void pnpWarning(const char* str,...)
{
	va_list arg_list;
	va_start(arg_list,str);
	printf("WARNING: ");
	vprintf(str,arg_list);
	va_end(arg_list);
}
inline void pnpError(const char* str,...)
{
	va_list arg_list;
	va_start(arg_list,str);
	printf("ERROR: ");
	vprintf(str,arg_list);
	va_end(arg_list);
}

#define pnpThrowError(Message) pnpThrowErrorRaw(__LINE__,__FILE__,Message)

inline void pnpThrowErrorRaw(const int line, const char* filename, const char* str, ...)
{
	int i;
	char *s = new char[10240];

	i = sprintf(s, "(%s:%d): ", filename, line);

	va_list arg_list;
	va_start(arg_list, str);
	i += vsprintf(s + i, str, arg_list);
	va_end(arg_list);

	std::runtime_error e(s);
	delete s;

	throw e;
}

#if __GNUC__ < 3 && !defined(_MSC_VER)
inline int roundf(float f)
{
	return (int)(f+0.5);
}

#endif

#if defined(_MSC_VER)

#ifndef M_PI
#  define M_PI        3.14159265358979323846
#endif

#if _MSC_VER <= 1600
inline int roundf(float f)
{
	return (int)(f+0.5);
}
#endif

#ifdef USE_WXLOG
	#define DbgPrint0 wxLogMessage
	#define DbgPrint1 wxLogMessage
	#define DbgPrint2 wxLogMessage
#else

#endif


#endif

#if !defined(MPI_PARALLEL) 
		inline void pnpPrint(const char* str,...)
		{
			va_list arg_list;
			va_start(arg_list,str);
			vprintf(str,arg_list);
			va_end(arg_list);
		}
		inline void pnpPrint0(const char* str,...)
		{
			va_list arg_list;
			va_start(arg_list,str);
			vprintf(str,arg_list);
			va_end(arg_list);
		}
		inline void pnpPrintMain0(const char* str,...)
		{
			va_list arg_list;
			va_start(arg_list,str);
			vprintf(str,arg_list);
			va_end(arg_list);
		}
		inline void pnpPrintGroup0(const char* str,...)
		{
			va_list arg_list;
			va_start(arg_list,str);
			vprintf(str,arg_list);
			va_end(arg_list);
		}
#  else
		inline void pnpPrint(const char* str,...)
		{
			va_list arg_list;
			va_start(arg_list,str);
			vprintf(str,arg_list);
			va_end(arg_list);
		}
		inline void pnpPrint0(const char* str,...)
		{
			va_list arg_list;
			va_start(arg_list,str);
			if(pnpsapp!=NULL)
				printf("(%d/%d): ",pnpsapp->GetMyAbsRank(), pnpsapp->GetTotalProc());
			vprintf(str,arg_list);
			va_end(arg_list);
		}
		inline void pnpPrintMain0(const char* str,...)
		{
			if(pnpsapp->AmIBigBrother())
			{
				va_list arg_list;
				va_start(arg_list,str);
				printf("(%d/%d): ",pnpsapp->GetMyAbsRank(), pnpsapp->GetTotalProc());
				vprintf(str,arg_list);
				va_end(arg_list);
			}
		}
		inline void pnpPrintGroup0(const char* str,...)
		{
			if(pnpsapp->AmIGroupLeader())
			{
				va_list arg_list;
				va_start(arg_list,str);
				printf("(%d/%d): ",pnpsapp->GetMyAbsRank(), pnpsapp->GetTotalProc());
				vprintf(str,arg_list);
				va_end(arg_list);
			}
		}
#  endif

#ifdef DBG0
  #ifndef MPI_PARALLEL
    inline void DbgPrint0(const char* str,...)
    {
      va_list arg_list;
      va_start(arg_list,str);
      printf("DBG0: ");
      vprintf(str,arg_list);
      va_end(arg_list);
    }
  #else
    inline void DbgPrint0(const char* str,...)
    {
      va_list arg_list;
      va_start(arg_list,str);
			printf("DBG0 (%d/%d): ", pnpsapp->GetMyAbsRank(), pnpsapp->GetTotalProc());
      vprintf(str,arg_list);
      va_end(arg_list);
    }
  #endif
  #ifdef DBG1
    #ifndef MPI_PARALLEL
      inline void DbgPrint1(const char* str,...)
      {
        va_list arg_list;
        va_start(arg_list,str);
        printf("DBG1: ");
        vprintf(str,arg_list);
        va_end(arg_list);
      }
    #else
      inline void DbgPrint1(const char* str,...)
      {
        va_list arg_list;
        va_start(arg_list,str);
				//printf("DBG1: ");
				printf("DBG1 (%d/%d): ", pnpsapp->GetMyAbsRank(), pnpsapp->GetTotalProc());
        vprintf(str,arg_list);
        va_end(arg_list);
      }
    #endif
    #ifdef DBG2
      #ifndef MPI_PARALLEL
        inline void DbgPrint2(const char* str,...)
        {
          va_list arg_list;
          va_start(arg_list,str);
          printf("DBG2: ");
          vprintf(str,arg_list);
          va_end(arg_list);
        }
      #else
        inline void DbgPrint2(const char* str,...)
        {
          va_list arg_list;
          va_start(arg_list,str);
					//printf("DBG2: ");
					printf("DBG2 (%d/%d): ", pnpsapp->GetMyAbsRank(), pnpsapp->GetTotalProc());
          vprintf(str,arg_list);
          va_end(arg_list);
        }
      #endif
    #else
      #define DbgPrint2(...)
    #endif
  #else
    #define DbgPrint1(...)
    #define DbgPrint2(...)
  #endif
  
#else
  #define DbgPrint0(...)
  #define DbgPrint1(...)
  #define DbgPrint2(...)
  //#define StartClock0
  //#define StopClock0
  //#define StopClockWMes0(Massege)
#endif

#if !defined(_MSC_VER) // TEMPORAL FIX IGOR for MS VS compiler
#define DefClock0 clock_t time0;timeval tvtime0,tvtime1;
#define StartClock0 time0=clock ();gettimeofday(&tvtime0,NULL);
#define StopClock0 gettimeofday(&tvtime1,NULL);DbgPrint0("Time : %g s(CPU Time) %g  s(Wall Time)\n",((double)(clock ()-time0))/CLOCKS_PER_SEC,double(tvtime1.tv_sec)+(double(tvtime1.tv_usec)/1000000.0)-double(tvtime0.tv_sec)-(double(tvtime0.tv_usec)/1000000.0));
#define StopClockWMes0(Massege) gettimeofday(&tvtime1,NULL);pnpPrint0("Time for %s is %.5g s (CPU Time) %g  s(Wall Time)\n",(Massege),((double)(clock ()-time0))/CLOCKS_PER_SEC,double(tvtime1.tv_sec)+(double(tvtime1.tv_usec)/1000000.0)-double(tvtime0.tv_sec)-(double(tvtime0.tv_usec)/1000000.0));
#else
#include <windows.h>
#define DefClock0 LARGE_INTEGER tvcount0,tvcount1,tvfreq;FILETIME tvtime0,tvtime1;double elapsedtime;
#define StartClock0    QueryPerformanceCounter(&tvcount0);GetSystemTimeAsFileTime(&tvtime0);
#define StopClock0Only QueryPerformanceCounter(&tvcount1);GetSystemTimeAsFileTime(&tvtime1);QueryPerformanceFrequency(&tvfreq);elapsedtime=(double)(((ULONGLONG)(tvtime1.dwHighDateTime-tvtime0.dwHighDateTime)<<32)|(ULONGLONG)(tvtime1.dwLowDateTime-tvtime0.dwLowDateTime))/10000000.0;
#define StopClock0  StopClock0Only; DbgPrint0("Time : %.3f s(CPU Time) %.3f  s(Wall Time)\n",((double)(tvcount1.QuadPart - tvcount0.QuadPart))/tvfreq.QuadPart,elapsedtime);
#define StopClockWMes0(Massege) StopClock0Only; pnpPrint0("Time for %s is %.5g s (CPU Time) \n\t%.5g  s(Wall Time) %.6f MCounts\n",(Massege),((double)(tvcount1.QuadPart - tvcount0.QuadPart))/tvfreq.QuadPart,elapsedtime,((double)(tvcount1.QuadPart - tvcount0.QuadPart))/1000000.0);
#endif

#endif