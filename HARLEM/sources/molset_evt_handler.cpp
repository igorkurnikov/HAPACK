/*! \file molset_evt_handler.cpp

    wxWidgets event handler class for HaMolSet
 
    \author Igor Kurnikov
    \date 2009-

*/  

#define MOLSET_EVT_HANDLER_CPP

#include "wx/event.h"
#include "hamolset.h"
#include "hawx_add.h"
#include "molset_evt_handler.h"

const wxEventType wxEVT_MOLSET = wxNewEventType();

BEGIN_EVENT_TABLE( MolSetEvtHandler,  wxEvtHandler)
   EVT_COMMAND( wxID_ANY, wxEVT_MOLSET, MolSetEvtHandler::OnTestCommand)
END_EVENT_TABLE()

MolSetEvtHandler::MolSetEvtHandler(HaMolSet* pmset_new)
{	
	pmset = pmset_new;
}

MolSetEvtHandler::~MolSetEvtHandler()
{

}

bool MolSetEvtHandler::TryParent(wxEvent& event)
{
	return false;
}

void MolSetEvtHandler::OnTestCommand(wxCommandEvent& event)
{
	PrintLog(" In: MolSetEvtHandler::OnTestCommand() \n");
	PrintLog(" MolSet Name=%s \n",pmset->GetName() );
	event.Skip();
}
