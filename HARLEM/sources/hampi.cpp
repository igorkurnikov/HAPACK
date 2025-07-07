#include <mpi.h>


#include "haconst.h"
#include "haio.h"
#include "wx/wxprec.h"

#include <sstream>

#include <boost/algorithm/string.hpp>

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "wx/dir.h"
#include "wx/filename.h"

#include <wx/filefn.h>
#include <wx/stream.h>
#include <wx/file.h>

#include "rapidxml.hpp"
#include "hastring.h"

#include "harlemapp.h"
#include "harlemapp_wx.h"
#include "hamolset.h"
#include "hampi.h"
#include "hawx_add.h"

HaMPI::HaMPI()
{
	myrank = 0;
	nprocs = 1;
	world_group = MPI_GROUP_NULL;

	msg_buffer.resize(10000);

	int ires;

	ires = MPI_Init(&(pApp->argc_loc), &(pApp->argv_loc)); 
	if(ires == 0)
	{
		ires = MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
		ires = MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
		ires = MPI_Comm_group(MPI_COMM_WORLD,&world_group);
	
//		std::string fname = "test_mpi" + harlem::ToString(myrank);
//		FILE* ftest = fopen(fname.c_str(),"w");
//	    fprintf(ftest,"MYRANK rank = %d : Num Processors= %d : after MPI_Init \n",myrank,nprocs);
//		fclose(ftest);
//		printf("MYRANK rank = %d : Num Processors= %d : after MPI_Init \n", myrank,nprocs);
	}
}

HaMPI::~HaMPI()
{
	if( world_group != MPI_GROUP_NULL) MPI_Group_free(&world_group);
	MPI_Finalize();
}

int HaMPI::Listen()
{
	using namespace rapidxml; 

	MPI_Status status;
	int i;
	int done = TRUE;
    int ires;
	int flag;

	HaMolSet* pmset = new HaMolSet();
	std::string mset_name = "MOLSET_MIRROR_" + harlem::ToString(myrank);
	pmset->SetName(mset_name.c_str());
	HaMolMechMod* p_mm_mod = pmset->GetMolMechMod(true);

	std::vector< wxCommandEvent* > wx_cmd_events;

//	PrintLog("HaMPI::Listen() start  myrank = %d\n",myrank );
	
	while(true)
	{
		try
		{
			ires = MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
//			ires = MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
//			PrintLog("HaMPI::Listen() pt 1  flag = %d\n", flag );
			if(!flag)
			{
				wxThread::Sleep(100);
				continue;
			}
//			PrintLog("HaMPI::Listen() pt 2 \n");
			//		ires = MPI_Bcast((void*)&basic_signal[0],BASIC_SIGNAL_DIM,MPI_INT,0,MPI_COMM_WORLD);
			ires = MPI_Recv((void*)&basic_signal[0],BASIC_SIGNAL_DIM,MPI_INT,
				MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);

			if(basic_signal[0] == XML_SIGNAL)
			{
//				PrintLog("HaMPI::Listen() pt 3  myrank = %d\n",myrank );
				int len = basic_signal[1];
				if( len <= 0) throw std::runtime_error(" invalid msg string len =" + harlem::ToString(len) );
				if( msg_buffer.size() < (len+1) ) msg_buffer.resize(len+1);
				ires = MPI_Bcast((void*)&msg_buffer[0],len,MPI_CHAR,0,MPI_COMM_WORLD);
				msg_buffer[len] = 0;

//				PrintLog("HaMPI::Listen() pt 4  msg recieved \n%s\n",msg_buffer.c_str() );

				xml_document<> doc;
				doc.parse<0>(&msg_buffer[0]);  
				xml_node<>* node1 = doc.first_node();
				if( node1 == NULL ) throw std::runtime_error(" No Head node in the XML message ");
				std::string tag = node1->name(); 
				boost::to_lower(tag);
				if( tag == "wxcommandevent" )
				{
					int type = 0;
					int id   = 0;
					xml_attribute<>* attr;
					for( attr = node1->first_attribute(); attr; attr = attr->next_attribute() )
					{
						std::string tag_attr = attr->name();
						boost::to_lower(tag_attr);
						if( tag_attr == "type" ) type = atoi( attr->value() );
						else if( tag_attr == "id" ) id = atoi( attr->value() );
					}
//					PrintLog("HaMPI::Listen() pt 5 \n type = %d  id = %d \n",type,id);  
					wx_cmd_events.push_back(new wxCommandEvent(type,id));
				}

				int nc = wx_cmd_events.size();
				int ic; 
				for(ic = 0; ic < nc; ic++)
				{
//					PrintLog("\n HaMPI::Listen() pt 6: wxCommandEvent Type = %d \n", wx_cmd_events[ic]->GetEventType());
//					PrintLog(" wxCommandEvent ID   = %d \n", wx_cmd_events[ic]->GetId());
					HarlemAppWX* p_wx_app = (HarlemAppWX*) pApp;
					p_wx_app->ProcessEvent(*wx_cmd_events[ic]);
					delete wx_cmd_events[ic];
				}
				wx_cmd_events.clear();
			}
			else if( basic_signal[0] == KILL_APP_SIGNAL )
			{
				break;
			}
		}
		catch( const std::exception& ex)
		{
			PrintLog("Error in HaMPI::Listen() \n");
			PrintLog("%s\n",ex.what());
		}
	}
	return TRUE;
}

std::string HaMPI::BuildXMLwxCmdEventBasic(int type, int id, bool add_header)
{
	std::stringstream os;

	if(add_header) 
	{
		os << harlem::StdXMLHeader() << std::endl;
	}
	os << "  <wxcommandevent type=\"" << type << "\" ";
	os << " id=\"" << id << "\" ";
	os << "/>" << std::endl;       

	return os.str();
}

int HaMPI::SendXmlMsgAllProc(const char* str)
{
//	PrintLog(" HaMPI::SendXmlMsgAllProc() \n %s\n",str);
	if( myrank != 0) 
	{
		PrintLog("Error in HaMPI::MPI_SendSignal(): Only Master can call it \n");
		return FALSE;
	}  
	int ierr;
	int len = strlen(str);

	basic_signal[0] = XML_SIGNAL;
	basic_signal[1] = len;
	basic_signal[2] = 0;
	basic_signal[3] = 0;

	vector<MPI_Request> req_vec;
	req_vec.resize(nprocs);

//	ierr = MPI_Bcast((void*)&basic_signal[0],BASIC_SIGNAL_DIM,MPI_INT,0,MPI_COMM_WORLD);
	int rank;
	for( rank = 1; rank < nprocs; rank++)
	{
		ierr = MPI_Isend(&basic_signal[0],BASIC_SIGNAL_DIM,MPI_INT,rank,0,MPI_COMM_WORLD,&req_vec[rank]);
	}

	if( len <= 0 ) 
	{
		PrintLog("Error in HaMPI::SendXmlMsgAllProc(): invalid msg string len = %d \n",len);
		return FALSE;
	}
	msg_buffer = str;
	
	ierr = MPI_Bcast((void*)&msg_buffer[0],len,MPI_CHAR,0,MPI_COMM_WORLD);
	
	return ierr;
}

int HaMPI::SendKillAppMsgAllProc()
{
	if( myrank != 0) 
	{
		PrintLog("Error in HaMPI::SendKillMsgAllProc(): Only Master can call it \n");
		return FALSE;
	}
	int ierr;

	vector<MPI_Request> req_vec;
	req_vec.resize(nprocs);

	basic_signal[0] = KILL_APP_SIGNAL;
	basic_signal[1] = 0;
	basic_signal[2] = 0;
	basic_signal[3] = 0;

//	ierr = MPI_Bcast((void*)&basic_signal[0],BASIC_SIGNAL_DIM,MPI_INT,0,MPI_COMM_WORLD);
	int rank;
	for( rank = 1; rank < nprocs; rank++)
	{
		ierr = MPI_Isend(&basic_signal[0],BASIC_SIGNAL_DIM,MPI_INT,rank,0,MPI_COMM_WORLD,&req_vec[rank]);
	}

	return ierr;
}


void HaMPI::ExecuteCommandProcArray(HaVec_int& proc_array, const char* cmd)
{
	
}

void HaMPI::ExecuteCommandAllProc(const char* cmd)
{
	
}
