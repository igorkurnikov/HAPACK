//
// C++ Interface: pnpsappess
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PNPPROCESS_H
#define PNPPROCESS_H

#ifdef MPI_PARALLEL
#  include "mpi.h"
#endif

#include <stdarg.h>
#include <stdio.h>
#include <time.h>
#if !defined(_MSC_VER) 
#include <sys/time.h>
#else
#include <winsock.h>
#endif

#include "vector"
#include <zlib.h>

#include <string>

#include "pnps.h"

class TiXmlElement;

/**
Class For Parallelization
	@author Nikolay Simakov <nsimakov@andrew.cmu.edu>
*/
class PNPSApp
{
	private:
		PNPSApp();
		int init(int np,int ng,int nppg,const char *t_TempDir, int SuperGroupNumber=0, int FirstSuperGroupProc=0);
		~PNPSApp();
	private:
		int MyAbsRank;///< Absolute rank of process
		int TotalProcs;///< Total number of proceses
		int MyRankInGroup;///< Rank of process in group
		int MyGroupNumber;///< Get group number
		int NProcsInMyGroup;///< Total number of process within group
		int NGroups;///< Total number of groups
		char TempDir[PNP_MAP_IO_STRING_LENGTH];///<TempDirectory to store temp results
#ifdef MPI_PARALLEL
	public:
		MPI::Group MyGroup;
		MPI::Intracomm MyComGroup;
	//	MPI::Group MySuperGroup;
		MPI::Intracomm MyComSuperGroup;
#endif
	public:
		inline int GetMyAbsRank(){return MyAbsRank;}
		inline int GetTotalProc(){return TotalProcs;}
		
		//!return number of the group to which current proccess belong
		inline int GetMyGroupNumber(){return MyGroupNumber;}
		//!return total number of groups
		inline int GetNumberOfGroups(){return NGroups;}
		
		//!return number of current proccess in the group
		inline int GetMyRankInGroup(){return MyRankInGroup;}
		//!return number of procceses in curent process group
		inline int GetNumProcsInMyGroup(){return NProcsInMyGroup;}
		inline int GetMyGroupLeader(){return 0;}
		inline int GetMaster(){return 0;}
		//!return number of procceses in GroupNumber group
		int GetNumProcsInGroup(int GroupNumber);
		
		inline bool AmIGroupLeader(){if(MyRankInGroup==0)return true;else return false;}
		inline bool AmIBigBrother(){if(MyAbsRank==0)return true;else return false;}
		
		inline int GetMasterProc(){return 0;}
		
		int SendCStr(int dest,const char *CStr);
		int RecvCStr(int dest,char *CStr);
    
		int BcastCStr(char *CStr,int root=0);
		TiXmlElement*  BcastTiXmlElement(TiXmlElement *Elt);
		TiXmlElement*  BcastTiXmlElementWithinGroup(TiXmlElement *Elt);
		
		int GetNumOfThreads() const;
		int SetNumOfThreads(const int m_NumOfThread);
	public:
		static int InitPNPSApp();
		static int InitPNPSApp(int np,int ng,int nppg,const char *t_TempDir, int SuperGroupNumber=0, int FirstSuperGroupProc=0);
		static int DeletePNPSApp();
		static PNPSApp* GetPNPSApp();
		/** save array of floats in compressed way function return id which can be used for retrival of array using PullFloatArray*/
		//int DropFloatArray(float *V,int N);
		//float* PullFloatArray(int id);
		int DropArray(Bytef *source,uLong sourceLen);
		int PullArray(Bytef* Out,int id);
		//! compressing pull out from zlib, level=6 default, diffrence is that comppress3 will allocate dest memory 
		static int compress3(Bytef **dest, uLongf *destLen, const Bytef *source, uLong sourceLen, int level);
		//!will add Process Number To filename, i.e. was file.ext -> filepref1.ext
		static void AddProcessNumberToFileName(char * out,const char *in,const char *pref,int pnum,int totnum);
		
		void AddMyAbsProcNumberToFileName(char * out,const char *in);
		void AddMyGroupNumberToFileName(char * out,const char *in);
		void AddMyGroupNumberToFileNameStdStr(std::string *filename);
		void AddMyGroupStringNameToFileNameStdStr(std::string *filename);
		void AddGroupStringNameToFileNameStdStr(std::string *filename,int GroupNumber);
		char *AddGroupStringNameToFileNameCStr(const char *filename,int GroupNumber);
	public:
		
	protected:
		std::vector<Bytef*> Vcom;
		std::vector<uLong> OrigLens;
		std::vector<uLong> CompLens;
		std::vector<std::string> GroupStringNames;
	public:
		static int RedirectOutputToFile(const char *FileName,bool AddGroupStringName);
		static int RedirectOutputToStd();
		
	private:
		
		int PrintPNPSAppInfo();
	public:
		bool HaveTimeToRun();//!<return trun if application still have time to run, return false then time close to end(only mpi)
		
		bool bTimerOn;
		int SetWallTimeToRun(const char *walltime);
		timeval StartTime;
		double SecondsToRun;
	public:
};
		const int MaxNumOfParts=24;
		const int MaxNPartSize=104857600;
		extern PNPSApp * pnpsapp;
#endif
