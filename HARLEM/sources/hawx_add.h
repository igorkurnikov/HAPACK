/*! \file hawx_add.h

    Additional modified wxWidget classes for HARLEM application
    \author Igor Kurnikov  
    \date 2009-
*/

#if !defined(HAWX_ADD_H)
#define HAWX_ADD_H

#include "wx/string.h"
#include "wx/event.h"

#if wxUSE_GUI == 0
// Define wxCommandEvent class that is absent in base wxWidgets

class WXDLLIMPEXP_CORE wxCommandEvent : public wxEvent
{
public:
    wxCommandEvent(wxEventType commandType = wxEVT_NULL, int winid = 0);

    wxCommandEvent(const wxCommandEvent& event)
        : wxEvent(event),
          m_cmdString(event.m_cmdString),
          m_commandInt(event.m_commandInt),
          m_extraLong(event.m_extraLong),
          m_clientData(event.m_clientData),
          m_clientObject(event.m_clientObject)
        { }

    // Set/Get client data from controls
    void SetClientData(void* clientData) { m_clientData = clientData; }
    void *GetClientData() const { return m_clientData; }

    // Set/Get client object from controls
    void SetClientObject(wxClientData* clientObject) { m_clientObject = clientObject; }
    wxClientData *GetClientObject() const { return m_clientObject; }

    // Get listbox selection if single-choice
    int GetSelection() const { return m_commandInt; }

    // Set/Get listbox/choice selection string
    void SetString(const wxString& s) { m_cmdString = s; }
    wxString GetString() const;

    // Get checkbox value
    bool IsChecked() const { return m_commandInt != 0; }

    // true if the listbox event was a selection.
    bool IsSelection() const { return (m_extraLong != 0); }

    void SetExtraLong(long extraLong) { m_extraLong = extraLong; }
    long GetExtraLong() const { return m_extraLong; }

    void SetInt(int i) { m_commandInt = i; }
    int GetInt() const { return m_commandInt; }

    virtual wxEvent *Clone() const { return new wxCommandEvent(*this); }

protected:
    wxString          m_cmdString;     // String event argument
    int               m_commandInt;
    long              m_extraLong;     // Additional information (e.g. select/deselect)
    void*             m_clientData;    // Arbitrary client data
    wxClientData*     m_clientObject;  // Arbitrary client object

private:
    DECLARE_DYNAMIC_CLASS_NO_ASSIGN(wxCommandEvent)
};

typedef void (wxEvtHandler::*wxCommandEventFunction)(wxCommandEvent&);

#define wxCommandEventHandler(func) \
    (wxObjectEventFunction)(wxEventFunction)wxStaticCastEvent(wxCommandEventFunction, &func)

#endif

#endif // !defined(HAWX_ADD_H)
