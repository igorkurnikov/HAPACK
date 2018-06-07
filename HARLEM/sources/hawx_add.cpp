/*! file hawx_add.cpp

    Additional modified wxWidget classes for HARLEM application       

    \author Igor Kurnikov 
    \date 2003-

*/

#define HAWX_ADD_CPP

#include "hawx_add.h"

#if wxUSE_GUI == 0 

IMPLEMENT_DYNAMIC_CLASS(wxCommandEvent, wxEvent)

wxCommandEvent::wxCommandEvent(wxEventType commandType, int theId)
              : wxEvent(theId, commandType)
{
    m_clientData = (char *) NULL;
    m_clientObject = (wxClientData *) NULL;
    m_extraLong = 0;
    m_commandInt = 0;
    m_isCommandEvent = true;

    // the command events are propagated upwards by default
    m_propagationLevel = wxEVENT_PROPAGATE_MAX;
}

wxString wxCommandEvent::GetString() const
{
        return m_cmdString;
}

#endif
