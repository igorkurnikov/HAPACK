//
// C++ Implementation: pnpsapp
//
// Description: 
//
//
// Author: Nikolay Simakov <nsimakov@andrew.cmu.edu>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifdef MPI_PARALLEL
#  include "mpi.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>


#include "pnpsapp.h"
#include "tinyxml.h"
#include "zlib.h"
#include "pnpdebug.h"
#include "mapio.h"

int TotalPNPWarningMassage=0;
int TotalPNPErrorMassage=0;
PNPSApp * pnpsapp=NULL;

FILE *FileOut=NULL;
FILE *PNPStdOut;
FILE *PNPStdErr;
PNPSApp::PNPSApp()
{
	SetNumOfThreads(-2);
}
int PNPSApp::init(int np,int ng,int nppg,const char *t_TempDir, int SuperGroupNumber, int FirstSuperGroupProc)
{
	int i;
	MyAbsRank=-1;
	TotalProcs=np;
	MyRankInGroup=-1;
	MyGroupNumber=-1;
	NProcsInMyGroup=nppg;
	NGroups=ng;
	bTimerOn=false;
	
#ifdef MPI_PARALLEL
	MPI::Status	status;
	//Create MySuperGroup
	MyAbsRank = MPI::COMM_WORLD.Get_rank();
	MyRankInGroup = MyAbsRank-FirstSuperGroupProc;
	MPI::COMM_WORLD.Barrier();
	MyComSuperGroup = MPI::COMM_WORLD.Split(SuperGroupNumber,MyRankInGroup);
	MPI::COMM_WORLD.Barrier();
	//InfoPrint("MyComSuperGroup.Get_size()=%d MyComSuperGroup.Get_rank()=%d MyGroupNumber=%d MyRankInGroup=%d\n", MyComSuperGroup.Get_size(), MyComSuperGroup.Get_rank(),SuperGroupNumber,MyRankInGroup);
	//Create MyGroup
	
	MyAbsRank = MyComSuperGroup.Get_rank();
	if(TotalProcs!=MyComSuperGroup.Get_size())
	{
		pnpError("Total number of processors is not coinside with given\n");
		pnpError("TotalProcs=%d -np=%d\n",MyComSuperGroup.Get_size(),TotalProcs);
		return EXIT_FAILURE;
	}
	if(NProcsInMyGroup*NGroups!=TotalProcs)
	{
		pnpError("Total number of processors, groups and proccesses per group is incorrect\n");
		return EXIT_FAILURE;
	}
	MyComSuperGroup.Barrier();
	
	MyRankInGroup=MyAbsRank%NProcsInMyGroup;
	MyGroupNumber=MyAbsRank/NProcsInMyGroup;
	
	MyComGroup= MyComSuperGroup.Split(MyGroupNumber,MyRankInGroup);
	/* Extract the original group handle */
	//MPI::Group orig_group=MPI::COMM_WORLD.Get_group();
	
	/* Divide tasks into distinct groups based upon rank */
	//int *ranks;
	//ranks=new int[NProcsInMyGroup];
	//for(i=0;i<NProcsInMyGroup;i++)
	//	ranks[i]=i+NProcsInMyGroup*MyGroupNumber;
	
	//MyGroup=orig_group.Incl(NProcsInMyGroup, ranks);
	
	/* Create new new communicator and then perform collective communications */
	//MyComGroup= MPI::COMM_WORLD.Create(MyGroup);
	MyRankInGroup=MyComGroup.Get_rank();
	NProcsInMyGroup=MyComGroup.Get_size();
	//MPI_Allreduce(&sendbuf, &recvbuf, 1, MPI_INT, MPI_SUM, new_comm);
	//MPI_Group_rank (new_group, &new_rank);
	//printf("rank= %d newrank= %d recvbuf= %d\n",rank,new_rank,recvbuf);
	
	//MyComGroup;
#else
	//Single CPU version
	MyAbsRank=0;
	TotalProcs=1;
	MyRankInGroup=0;
	MyGroupNumber=0;
	NProcsInMyGroup=1;
	NGroups=1;
#endif
	sprintf(TempDir,"%s\0",t_TempDir);
	return EXIT_SUCCESS;
}
int PNPSApp::PrintPNPSAppInfo()
{
	//Print information about this particular problem
	pnpPrint("<AboutPMFPNPSolver\n");
	pnpPrint("\tversion=\"%s\"\n",PNPSVERSION);
#ifdef MPI_PARALLEL
	pnpPrint("\tsubversion=\"MPI parallel\"\n");
#else
	pnpPrint("\tsubversion=\"Single CPU\"\n");
#endif
	pnpPrint("\tcompiler=\"\"\n");
	pnpPrint("\ttime_of_compilation=\"%s\"\n",__TIME__);
	pnpPrint("\tdate_of_compilation=\"%s\"\n",__DATE__);
	pnpPrint("\tMyAbsRank=\"%d\" \tTotalProcs = \"%d\"\n",
			GetMyAbsRank(),GetTotalProc());
	pnpPrint("\tMyRankInGroup=\"%d\" \tMyGroup=\"%d\" NumProcsInMyGroup=\"%d\" \tNumberOfGroups=\"%d\"\n",
					GetMyRankInGroup(), GetMyGroupNumber(),
					GetNumProcsInMyGroup(),GetNumberOfGroups());
#ifdef MPI_PARALLEL
	char *name;
	int namelen;
	//MPI::Get_processor_name(name, namelen);
	//pnpPrint("\tprocessor_name=\"%s\"\n",name);
#endif
	pnpPrint("\t/>\n");
	return EXIT_SUCCESS;
}
int PNPSApp::RedirectOutputToFile(const char *FileName,bool AddGroupStringName)
{
	//Std and Err Output
	InfoPrint("OutputFile=%s\n",FileName);
	char OutputFileLocal[PNP_MAP_IO_STRING_LENGTH];
	if(FileName!=NULL)
	{
#ifdef MPI_PARALLEL
// 		if(AddGroupStringName&&pnpsapp!=NULL)
// 		{
// 			char OutputFileLocal2[PNP_MAP_IO_STRING_LENGTH];
// 			if(pnpsapp->GetNumProcsInMyGroup()>1)
// 			{
// 				pnpsapp->AddMyGroupStringNameToFileNameStdStr(OutputFileLocal2,FileName);
// 				PNPSApp::AddProcessNumberToFileName(OutputFileLocal,OutputFileLocal2,"_p",pnpsapp->GetMyRankInGroup(),pnpsapp->GetNumProcsInMyGroup());
// 			}
// 			else
// 			{
// 				pnpsapp->AddMyGroupStringNameToFileNameStdStr(OutputFileLocal2,FileName);
// 			}
// 		}
// 		else
// 		{
		if(pnpsapp==NULL)
			PNPSApp::AddProcessNumberToFileName(OutputFileLocal,FileName,"_p",MPI::COMM_WORLD.Get_rank(),MPI::COMM_WORLD.Get_size());
		else
			PNPSApp::AddProcessNumberToFileName(OutputFileLocal,FileName,"_p",pnpsapp->MyComSuperGroup.Get_rank(),pnpsapp->MyComSuperGroup.Get_size());
// 		}
#else
		sprintf(OutputFileLocal,"%s\0",FileName);
#endif
		FileOut=fopen(OutputFileLocal,"w");
		PNPStdOut=stdout;
		PNPStdErr=stderr;
#if !defined(_MSC_VER) // TMP FIX IGOR
		stdout=FileOut;
		stderr=FileOut;
#endif
		setvbuf( stdout, NULL, _IONBF, 0 );
	}
	else
	{
		sprintf(OutputFileLocal,"NULL\0");
	}
	InfoPrint("OutputFile=%s\n",OutputFileLocal);
	return EXIT_SUCCESS;
}
int PNPSApp::RedirectOutputToStd()
{
	if(FileOut!=NULL)
	{
		fclose(FileOut);
		FileOut=NULL;
#if !defined(_MSC_VER) // TMP FIX IGOR
		stdout=PNPStdOut;
		stderr=PNPStdErr;
#endif
	}
	return EXIT_SUCCESS;
}
PNPSApp::~PNPSApp()
{
	//if(FileOut!=NULL)fclose(FileOut);
}
int PNPSApp::GetNumProcsInGroup(int GroupNumber)
{
	return 1;
}
PNPSApp* PNPSApp::GetPNPSApp()
{
	if (pnpsapp == NULL)
		InitPNPSApp();
	return pnpsapp;
}
int PNPSApp::InitPNPSApp()
{
	return InitPNPSApp(1,1,1,"/tmp");
}
int PNPSApp::InitPNPSApp(int np,int ng,int nppg,const char *t_TempDir, int SuperGroupNumber, int FirstSuperGroupProc)
{
	static bool Initiated=false;
	if((!Initiated)||(pnpsapp==NULL)||(pnpsapp->GetTotalProc()==1))
	{
		pnpsapp = new PNPSApp();
		int status = pnpsapp->init(np, ng, nppg,t_TempDir,SuperGroupNumber,FirstSuperGroupProc);
		if(status==EXIT_SUCCESS)
		{
			Initiated=true;
		}
		if(status==EXIT_FAILURE)
		{
			DeleteObjByPnt(pnpsapp);
			Initiated=false;
			return EXIT_FAILURE;
		}
	}
	return EXIT_SUCCESS;
}
int PNPSApp::DeletePNPSApp()
{
	DeleteObjByPnt(pnpsapp);
	return EXIT_SUCCESS;
}
int PNPSApp::SendCStr(int dest,const char *CStr)
{
#ifdef MPI_PARALLEL
	int count=0;
	while(CStr[count]!='\0')count++;
	count++;
	fprintf(stdout,"Str=%s Count=%d\n",CStr,count);
	pnpsapp->MyComSuperGroup.Send(&count, 1, MPI::INT, dest, 0);
	pnpsapp->MyComSuperGroup.Send(CStr, count, MPI::CHAR, dest, 0);
#endif
	return EXIT_SUCCESS;
}
int PNPSApp::RecvCStr(int dest,char *CStr)
{
#ifdef MPI_PARALLEL
	int count=0;
	MPI::Status  status;
	pnpsapp->MyComSuperGroup.Recv(&count, 1, MPI::INT, dest, 0, status);
	pnpsapp->MyComSuperGroup.Recv(CStr, count, MPI::CHAR, dest, 0, status);
#endif
	return EXIT_SUCCESS;
}
int PNPSApp::BcastCStr(char *CStr,int root)
{
#ifdef MPI_PARALLEL
	int count=0;
	while(CStr[count]!='\0')count++;
	count++;
	pnpsapp->MyComSuperGroup.Bcast(&count, 1, MPI::INT, root);
	pnpsapp->MyComSuperGroup.Bcast(CStr, count, MPI::CHAR, root);
#endif
	return EXIT_SUCCESS;
}

int PNPSApp::GetNumOfThreads() const
{
#ifdef _OPENMP
	return omp_get_max_threads();
#else
	return 1;
#endif
}
int PNPSApp::SetNumOfThreads(const int m_NumOfThread)
{
#ifdef _OPENMP
	//printf("Set threads number to %d\n\n", num_threads);
	if (m_NumOfThread <= 0) {
		omp_set_num_threads(omp_get_num_procs());
	}
	else {
		omp_set_num_threads(m_NumOfThread);
	}
#else
	if (m_NumOfThread>1)
		printf("Warning: This version was compiled without openmp support!\n");
#endif
	return GetNumOfThreads();
}

TiXmlElement* PNPSApp::BcastTiXmlElement(TiXmlElement *Elt)
{
#ifdef MPI_PARALLEL
	int count=0,i;
	char * CStr;
	std::string str;
	if(GetMyAbsRank()==GetMaster())
	{
		str<<*Elt;
		count=str.size();
		
		pnpsapp->MyComSuperGroup.Bcast(&count, 1, MPI::INT, GetMaster());
		CStr=new char[count+1];
		for(i=0;i<count;i++)CStr[i]=str[i];
		CStr[count]='\0';
		pnpsapp->MyComSuperGroup.Bcast(CStr, count+1, MPI::CHAR, GetMaster());
		delete [] CStr;
		return Elt;
	}
	else
	{
		pnpsapp->MyComSuperGroup.Bcast(&count, 1, MPI::INT, GetMaster());
		CStr=new char[count+1];
		pnpsapp->MyComSuperGroup.Bcast(CStr, count+1, MPI::CHAR, GetMaster());
		std::istringstream ins(CStr);
		if(Elt==NULL)Elt=new TiXmlElement("DontKnowName");
		ins>>*Elt;
		delete [] CStr;
		return Elt;
	}
#else
	return Elt;
#endif
}
TiXmlElement* PNPSApp::BcastTiXmlElementWithinGroup(TiXmlElement *Elt)
{
#ifdef MPI_PARALLEL
	int count=0,i;
	char * CStr;
	std::string str;
	if(GetMyRankInGroup()==GetMyGroupLeader())
	{
		str<<*Elt;
		count=str.size();
		
		MyComGroup.Bcast(&count, 1, MPI::INT, GetMyGroupLeader());
		CStr=new char[count+1];
		for(i=0;i<count;i++)CStr[i]=str[i];
		CStr[count]='\0';
		MyComGroup.Bcast(CStr, count+1, MPI::CHAR, GetMyGroupLeader());
		delete [] CStr;
		return Elt;
	}
	else
	{
		MyComGroup.Bcast(&count, 1, MPI::INT, GetMyGroupLeader());
		CStr=new char[count+1];
		MyComGroup.Bcast(CStr, count+1, MPI::CHAR, GetMyGroupLeader());
		std::istringstream ins(CStr);
		if(Elt==NULL)Elt=new TiXmlElement("DontKnowName");
		ins>>*Elt;
		delete [] CStr;
		return Elt;
	}
#else
	return Elt;
#endif
}
int PNPSApp::DropArray(Bytef *source,uLong sourceLen)
{
	int i;
	pnpPrint0("DropArray %p %d\n",source,sourceLen);

	DefClock0;
	StartClock0;
	uLongf destLen=0;
	Bytef *dest=NULL;
	this->compress3(&dest,&destLen, source, sourceLen,6);
	
	pnpPrint0("DropArray was %d become %d\n",sourceLen,destLen);
	StopClockWMes0("DropArray");
	int id=Vcom.size();
	
	Vcom.push_back(dest);
	OrigLens.push_back(sourceLen);
	CompLens.push_back(destLen);
	return id;
}
void PNPSApp::AddMyAbsProcNumberToFileName(char * out,const char *in)
{
	AddProcessNumberToFileName(out,in,"_p",MyAbsRank,TotalProcs);
}
void PNPSApp::AddProcessNumberToFileName(char * out,const char *in,const char *pref,int pnum, int totnum)
{
	char in_main[PNP_MAP_IO_STRING_LENGTH];
	char *in_ext,*in_det;
	strcpy(in_main,in);
	in_det=strrchr(in_main,'.');
	if(in_det!=NULL)
	{
		*in_det='\0';
		in_ext=in_det+1;
		if(totnum<=10)
			sprintf(out,"%s%s%.1d.%s\0",in_main,pref,pnum,in_ext);
		else if(totnum<=100)
			sprintf(out,"%s%s%.2d.%s\0",in_main,pref,pnum,in_ext);
		else
			sprintf(out,"%s%s%.3d.%s\0",in_main,pref,pnum,in_ext);
	}
	else
	{
		if(totnum<=10)
			sprintf(out,"%s%s%.1d\0",in_main,pref,pnum);
		else if(totnum<=100)
			sprintf(out,"%s%s%.2d\0",in_main,pref,pnum);
		else
			sprintf(out,"%s%s%.3d\0",in_main,pref,pnum);
	}
}
void PNPSApp::AddMyGroupNumberToFileName(char * out,const char *in)
{
	char in_main[PNP_MAP_IO_STRING_LENGTH];
	char *in_ext,*in_det;
	strcpy(in_main,in);
	in_det=strrchr(in_main,'.');
	if(in_det!=NULL)
	{
		*in_det='\0';
		in_ext=in_det+1;
		if(NGroups<=10)
			sprintf(out,"%s%s%.1d.%s\0",in_main,"_g",MyGroupNumber,in_ext);
		else if(NGroups<=100)
			sprintf(out,"%s%s%.2d.%s\0",in_main,"_g",MyGroupNumber,in_ext);
		else
			sprintf(out,"%s%s%.3d.%s\0",in_main,"_g",MyGroupNumber,in_ext);
	}
	else
	{
		if(NGroups<=10)
			sprintf(out,"%s%s%.1d\0",in_main,"_g",MyGroupNumber);
		else if(NGroups<=100)
			sprintf(out,"%s%s%.1d\0",in_main,"_g",MyGroupNumber);
		else
			sprintf(out,"%s%s%.1d\0",in_main,"_g",MyGroupNumber);
	}
}
void PNPSApp::AddMyGroupNumberToFileNameStdStr(std::string *filename)
{
	char c_str[PNP_MAP_IO_STRING_LENGTH];
	AddMyGroupNumberToFileName(c_str,filename->c_str());
	filename->assign(c_str);
}
void PNPSApp::AddGroupStringNameToFileNameStdStr(std::string *filename,int GroupNumber)
{
	if(GroupStringNames.size()>GroupNumber)
	{
		char in_main[PNP_MAP_IO_STRING_LENGTH];
		char out[PNP_MAP_IO_STRING_LENGTH];
		char *in_ext,*in_det;
		strcpy(in_main,filename->c_str());
		in_det=strrchr(in_main,'.');
		if(in_det!=NULL)
		{
			*in_det='\0';
			in_ext=in_det+1;
			sprintf(out,"%s_%s.%s\0",in_main,GroupStringNames[GroupNumber].c_str(),in_ext);
		}
		else
		{
			sprintf(out,"%s_%s%\0",in_main,GroupStringNames[GroupNumber].c_str());
		}
		filename->assign(out);
	}
	else
	{
		char c_str[PNP_MAP_IO_STRING_LENGTH];
		AddMyGroupNumberToFileName(c_str,filename->c_str());
		filename->assign(c_str);
	}
}
char * PNPSApp::AddGroupStringNameToFileNameCStr(const char *filename,int GroupNumber)
{
	char *out=new char[PNP_MAP_IO_STRING_LENGTH];
	DbgPrint0("PNPSApp::AddGroupStringNameToFileNameCStr filename=%s GroupNumber=%d\n",filename,GroupNumber);
	if(GroupStringNames.size()>GroupNumber)
	{
		char in_main[PNP_MAP_IO_STRING_LENGTH];
		
		char *in_ext,*in_det;
		strcpy(in_main,filename);
		in_det=strrchr(in_main,'.');
		if(in_det!=NULL)
		{
			*in_det='\0';
			in_ext=in_det+1;
			sprintf(out,"%s_%s.%s\0",in_main,GroupStringNames[GroupNumber].c_str(),in_ext);
		}
		else
		{
			sprintf(out,"%s_%s%\0",in_main,GroupStringNames[GroupNumber].c_str());
		}
		return out;
	}
	else
	{
		pnpWarning("GroupStringNames.size()>GroupNumber will add group number\n");
		AddMyGroupNumberToFileName(out,filename);
		return out;
	}
}
void PNPSApp::AddMyGroupStringNameToFileNameStdStr(std::string *filename)
{
	DbgPrint0("PNPSApp::PNPSApp::AddMyGroupStringNameToFileNameStdStr filename=%s\n",filename->c_str());
	if(MyGroupNumber<GroupStringNames.size())
	{
		char in_main[PNP_MAP_IO_STRING_LENGTH];
		char out[PNP_MAP_IO_STRING_LENGTH];
		char *in_ext,*in_det;
		strcpy(in_main,filename->c_str());
		in_det=strrchr(in_main,'.');
		if(in_det!=NULL)
		{
			*in_det='\0';
			in_ext=in_det+1;
			sprintf(out,"%s_%s.%s\0",in_main,GroupStringNames[GetMyGroupNumber()].c_str(),in_ext);
		}
		else
		{
			sprintf(out,"%s_%s%\0",in_main,GroupStringNames[GetMyGroupNumber()].c_str());
		}
		filename->assign(out);
	}
	else
	{
		char c_str[PNP_MAP_IO_STRING_LENGTH];
		AddMyGroupNumberToFileName(c_str,filename->c_str());
		filename->assign(c_str);
	}
}
int PNPSApp::PullArray(Bytef* Out,int id)
{
	uLongf destLen=OrigLens[id];
	Bytef *source=Vcom[id];
	uLong sourceLen=CompLens[id];
	DefClock0;
	StartClock0;
	uncompress (Out, &destLen, source, sourceLen);
	StopClockWMes0("PullArray");
	return EXIT_SUCCESS;
}
int PNPSApp::compress3(Bytef **dest, uLongf *destLen, const Bytef *source, uLong sourceLen, int level)
{
	z_stream stream;
	int err;
	uLong i;
	int countParts=0;
	Bytef* comp[MaxNumOfParts];
	for(i=0;i<MaxNumOfParts;i++)comp[i]=NULL;
	//*destLen=sourceLen*1.02+12;
	
	
	//for(i=0;i<destLen;i++)
	//	comp[i]=dest[i];
	
	stream.next_in = (Bytef*)source;
	stream.avail_in = (uInt)sourceLen;
	
	stream.next_out = comp[countParts];
	stream.avail_out = (uInt)*destLen;
	if ((uLong)stream.avail_out != *destLen) return Z_BUF_ERROR;
	
	stream.zalloc = (alloc_func)0;
	stream.zfree = (free_func)0;
	stream.opaque = (voidpf)0;
	
	*destLen=0;
	err = deflateInit(&stream, level);
	if (err != Z_OK) return err;
	
	do
	{
		comp[countParts]=new Bytef[MaxNPartSize];
		stream.next_out = comp[countParts];
		stream.avail_out = (uInt)MaxNPartSize;
		
		err = deflate(&stream, Z_FULL_FLUSH);
		*destLen+=stream.total_out;
		countParts++;
	}
	while (stream.avail_out == 0);
	
	DbgPrint0("Was allocated %d parts, destLen=%d\n",countParts,*destLen);
	
	*dest=new Bytef[*destLen];
	
	uLong j,k;
	for(i=0;i<countParts;i++)
	{
		if((i+1)*MaxNPartSize>*destLen)k=*destLen-i*MaxNPartSize;
		else k=MaxNPartSize;
		
		DbgPrint0("copy part %d part elms  %d\n",i,j);
		for(j=0;j<k;j++)
			(*dest)[i*MaxNPartSize+j]=comp[i][j-i*MaxNPartSize];
		delete [] (comp[i]);		
	}
	
	if (err != Z_STREAM_END) {
		deflateEnd(&stream);
		return err == Z_OK ? Z_BUF_ERROR : err;
	}
	
	err = deflateEnd(&stream);
	//DeleteCVecArray(comp,countParts);
	return err;
}
bool PNPSApp::HaveTimeToRun()
{
	if(!bTimerOn)//e.i. no time constrictions
		return true;
	timeval tvtime1;
#if !defined(_MSC_VER)
	gettimeofday(&tvtime1,NULL);
#endif
	double secondsrunning;
	secondsrunning=double(tvtime1.tv_sec)+(double(tvtime1.tv_usec)/1000000.0)-double(StartTime.tv_sec)-(double(StartTime.tv_usec)/1000000.0);
	
	if(secondsrunning>=SecondsToRun)
		return false;
	else
		return true;
}
int PNPSApp::SetWallTimeToRun(const char *walltime)
{
#if !defined(_MSC_VER)
	gettimeofday(&StartTime,NULL);
#endif
	int iH,iM,iS;
	bTimerOn=true;
	DbgPrint0("time to run %s\n",walltime);
	if(sscanf(walltime,"%d:%d:%d",&iH,&iM,&iS)!=3)
	{
		iH=0;
		if(sscanf(walltime,"%d:%d",&iM,&iS)!=2)
		{
			iM=0;
			if(sscanf(walltime,"%d",&iS)!=1)
			{
				pnpError("can not read wall time to run: %s\n",walltime);
				bTimerOn=false;
				return EXIT_FAILURE;
			}
		}
	}
	DbgPrint0("time to run %d:%d:%d\n",iH,iM,iS);
	SecondsToRun=double(iH)*60.0*60.0+double(iM)*60.0+double(iS);
	DbgPrint0("time to run %g s\n",SecondsToRun);
	return EXIT_SUCCESS;
}

