/*! \file dialogs_wx_3.h

    wxWidgets Dialogs 3 in HARLEM
    Dialogs created without wxDesigner

    \author Igor Kurnikov
    \date 2026-
*/

#if !defined(DIALOGS_WX_3_H)
#define DIALOGS_WX_3_H

#include <wx/listctrl.h>

class HaMolView;
class Object3D;
class MolSet;

class Object3DDlgWX : public wxDialog
{
public:
	Object3DDlgWX(HaMolView* new_pview, wxWindow *parent);
	virtual ~Object3DDlgWX();

	static int dlg_open;

	virtual bool TransferDataToWindow();

protected:
	HaMolView* pview;
	std::vector<Object3D*> obj_vec;

	wxListCtrl* obj_list_ctrl;
	wxTextCtrl* name_ctrl;
	wxTextCtrl* transp_ctrl;

	void CreateControls();
	void DDX_obj_list();

	void OnSetTransp(wxCommandEvent &event );
	void OnDelete(wxCommandEvent &event );
	void OnDisplay(wxCommandEvent &event );
	void OnUnDisplay(wxCommandEvent &event );
	void OnUpdate(wxCommandEvent &event );
	void OnClose(wxCloseEvent &event);
};

#endif // DIALOGS_WX_3_H
