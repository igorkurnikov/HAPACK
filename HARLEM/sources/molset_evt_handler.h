/*!  \file molset_evt_handler.h
  
     Class to define set of molecules in HARLEM and accompaniing iterator classes
 
    \author Igor Kurnikov  
    \date 1999-2002
*/
#if !defined(MOLSET_EVT_HANDLER_H)
#define MOLSET_EVT_HANDLER_H

#include <wx/event.h>

class MolSetEvtHandler: public wxEvtHandler
{
public:
	MolSetEvtHandler(MolSet* pmset);
	virtual ~MolSetEvtHandler();

	virtual bool TryParent(wxEvent& event); //!< Override wxEvtHandler to not calling wxApp::ProcessEvent()

	void OnTestCommand(wxCommandEvent& event); 

	MolSet* pmset;
private:
	DECLARE_EVENT_TABLE();
};

extern const wxEventType wxEVT_MOLSET;


#endif // !defined(MOLSET_EVT_HANDLER_H)


