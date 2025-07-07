/*! \file ctrl_wx.h

    Custom wxWidgets controls in HARLEM 
 
    \author Igor Kurnikov  
    \date 2005-
*/

#include "hastl.h"

class CmdTextCtrl: public wxTextCtrl
{
public:
	CmdTextCtrl(wxWindow* parent, wxWindowID id, const wxString& value = "",
		const wxPoint& pos = wxDefaultPosition,const wxSize& size = wxDefaultSize,
		long style = 0, const wxValidator& validator = wxDefaultValidator,
		const wxString& name = wxTextCtrlNameStr);

	void OnChar(wxKeyEvent& event);

private:
    DECLARE_EVENT_TABLE()
};

class HaAtom;

class wxAtomEdit: public wxTextCtrl
{
public:
	wxAtomEdit(wxWindow* parent, wxWindowID id, const wxString& value = "",
		const wxPoint& pos = wxDefaultPosition,const wxSize& size = wxDefaultSize,
		long style = 0, const wxValidator& validator = wxDefaultValidator,
		const wxString& name = wxTextCtrlNameStr);

	virtual ~wxAtomEdit();

	static list<wxAtomEdit*> active_controls;
	static int BroadCastPickedAtom(HaAtom* PickAtom);

	int OnAtomPicked(HaAtom* PickAtom);
	
	void OnSetFocus(wxFocusEvent& event);


	enum PICK_MODE { PICK_ATOM, PICK_RESIDUE, PICK_MOLECULE, PICK_GROUP} pick_mode;

private:
    DECLARE_EVENT_TABLE()
};
