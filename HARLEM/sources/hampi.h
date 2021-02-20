/*! \file hampi.h

	Class for parallel computing 

    \date 2005
    \author Igor Kurnikov 

*/

#if !defined(HARLEM_MPI_H)
#define HARLEM_MPI_H

#include <mpi.h>
#include <string>

class HaVec_int;

//!
//! Class to support parallel calculation in HARLEM
//!
class HaMPI
{
public:
	HaMPI();
	virtual ~HaMPI();	
	
	int Listen();                           //!< Wait for signal ( on slave)
	int SendXmlMsgAllProc(const char* msg); //!< Send XML Message To All Processors (from master)
	int SendKillAppMsgAllProc();            //!< Send Kill Application message to All Processors
	
	void ExecuteCommandProcArray(HaVec_int& proc_array, const char* cmd); //!< Send command to array of processors, return TRUE on success
	void ExecuteCommandAllProc(const char* cmd); //!< Execute Command on all processors

	int myrank; //!< Rank in  MPI_COMM_WORLD
	int nprocs; //!< Number of processes in MPI_COMM_WORLD

    MPI_Group world_group; //!< MPI_Group corresponding to MPI_COMM_WORLD

    static const int BASIC_SIGNAL_DIM = 4; 

	enum HARLEM_MPI_SIGNAL{ KILL_APP_SIGNAL = 0, XML_SIGNAL, WX_EVENT_SIGNAL };

    int basic_signal[BASIC_SIGNAL_DIM];
	std::string msg_buffer;

	static std::string BuildXMLwxCmdEventBasic(int type, int id, bool add_header = true); //!< Build XML string corresponding to wxCommandEvent no Object 

protected:
};

#endif // !defined(HARLEM_MPI) 


