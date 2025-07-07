/*!  \file molmech_evt_handler.h
  
     Class to define wxWdigets command processing for HaMolMechMod class
 
    \author Igor Kurnikov  
    \date 2009-
*/

#if !defined(MOLMECH_EVT_HANDLER)
#define MOLMECH_EVT_HANDLER

#include "wx/event.h"

class MolMechEvtHandler : public wxEvtHandler
{
public:
	MolMechEvtHandler(HaMolMechMod* p_mm_mod_new);
	virtual ~MolMechEvtHandler();

	virtual bool TryParent(wxEvent& event); //!< Override wxEvtHandler to not calling wxApp::ProcessEvent()

	void OnTestCommand(wxCommandEvent& event); //!< Test MolMech command 
	void OnTestMolMechEvent(wxCommandEvent& event); //!< Test wxEVT_MOL_MECH interception

	HaMolMechMod* p_mm_mod;

private:
	DECLARE_EVENT_TABLE();
};

extern const wxEventType wxEVT_MOL_MECH;

const int MOL_MECH_ID_TEST1 = 3000;
const int MOL_MECH_ID_TEST2 = 3001;
const int MM_DRIVER_AMBER_RUN_INTERNAL   = 3002;
const int MM_SET_MPI_COMM_ALL_PROCS      = 3003;
const int MM_MOD_SET_MPI_COMM_SPLIT_2    = 3004;
const int MM_INIT_SIMULATIONS_STEP_2     = 3005;
const int MM_MOD_INIT_MIXED_HAMILTONIAN  = 3006;
const int MM_UPDATE_CONSTR_2             = 3007;

#endif //! defined(MOLMECH_EVT_HANDLER)


