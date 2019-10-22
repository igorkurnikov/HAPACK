/*! \file molmech_evt_handler.cpp

    wxWidgets event handler class for HaMolMechMod
 
    \author Igor Kurnikov
    \date 2009-

*/  

#define MOLMECH_EVT_HANDLER_CPP

#include <mpi.h>

#include <wx/event.h>
#include "hamolset.h"
#include "hamolmech.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "mm_driver_amber.h"
#include "hawx_add.h"
#include "molmech_evt_handler.h"

const wxEventType wxEVT_MOL_MECH = wxNewEventType();

BEGIN_EVENT_TABLE( MolMechEvtHandler,  wxEvtHandler)
   EVT_COMMAND( wxID_ANY, wxEVT_MOL_MECH, MolMechEvtHandler::OnTestMolMechEvent)
END_EVENT_TABLE()

MolMechEvtHandler::MolMechEvtHandler(HaMolMechMod* p_mm_mod_new)
{
	p_mm_mod = p_mm_mod_new;
	MolSet* pmset = p_mm_mod->GetMolSet();

	wxEvtHandler* evt_h = (wxEvtHandler*) pmset->p_evt_h;
	wxEvtHandler* evt_h_last = evt_h;
	wxEvtHandler* evt_h_next = evt_h->GetNextHandler();

	while( evt_h != NULL )
	{
		evt_h_last = evt_h;
		evt_h = evt_h->GetNextHandler();
	}
	if(evt_h_last != NULL) 
	{	
		evt_h_last->SetNextHandler(this);
		this->SetPreviousHandler(evt_h_last);
	}
}

MolMechEvtHandler::~MolMechEvtHandler(void)
{
	wxEvtHandler* evt_h_prev = this->GetPreviousHandler();
	wxEvtHandler* evt_h_next = this->GetNextHandler();
	if( evt_h_prev != NULL ) evt_h_prev->SetNextHandler(evt_h_next);
	if( evt_h_next != NULL ) evt_h_next->SetPreviousHandler(evt_h_prev);
}

bool MolMechEvtHandler::TryParent(wxEvent& event)
{
	return false;
}

void MolMechEvtHandler::OnTestCommand(wxCommandEvent& event)
{
	MolSet* pmset = p_mm_mod->GetMolSet();
	PrintLog(" In: MolMechEvtHandler::OnTestCommand() \n");
	PrintLog(" MolSet Name=%s \n",pmset->GetName() );
	event.Skip();
}

void MolMechEvtHandler::OnTestMolMechEvent(wxCommandEvent& event)
{
	int id = event.GetId();
	MolSet* pmset = p_mm_mod->GetMolSet();
//	PrintLog(" In: MolMechEvtHandler::OnTestMolMechEvent() \n");
//	PrintLog(" MolSet Name=%s \n",pmset->GetName() );
//	PrintLog(" wxCommandEvent ID =%d \n",id);

	try 
	{
		switch(id)
		{
		case MM_INIT_SIMULATIONS_STEP_2:
			PrintLog("About to Execute InitSimulationsStep2() \n");
			p_mm_mod->p_amber_driver->InitSimulationsStep2();
			break;
		case MM_DRIVER_AMBER_RUN_INTERNAL:
			PrintLog("About to Execute RunInternal_node() \n");
			p_mm_mod->RunInternal_node();
			break;
		case MM_SET_MPI_COMM_ALL_PROCS:
			PrintLog("About to Execute SetMPICommAllProcs() \n");
			p_mm_mod->p_amber_driver->SetMPICommAllProcs();
			p_mm_mod->p_amber_driver->InitParallelDatMod();
			break;
		case MM_MOD_SET_MPI_COMM_SPLIT_2:
			PrintLog("About to Execute SetMPICommSplit2() \n");
			p_mm_mod->SetMPICommSplit2();
			p_mm_mod->p_amber_driver->InitParallelDatMod();
			break;
		case MM_MOD_INIT_MIXED_HAMILTONIAN:
			PrintLog("About to Execute HaMolMechMod::InitMixedHamSimulations_node()  \n");
			p_mm_mod->InitMixedHamSimulations_node( p_mm_mod->p_mm_model );
			break;
		case MM_UPDATE_CONSTR_2:
			PrintLog("About to Execute MolMechModel::UpdateConstraints_2() \n");
			p_mm_mod->GetMolMechModel()->UpdateConstraints_2();
			break;
		default:
			PrintLog("Error in MolMechEvtHandler::OnTestMolMechEvent() Unknown ID= %d \n",id); 
		}

	}
	catch(std::exception& e) 
	{
		PrintLog(" Exception executing MolMechEvtHandler::OnTestMolMechEvent() \n %s \n",e.what() );
	}
	catch(...) 
	{ 
		PrintLog(" Exception executing MolMechEvtHandler::OnTestMolMechEvent() \n"); 
	} 

	event.Skip();
}
