/*! \file dialogs_wx_3.cpp

    wxWidgets Dialogs 3 in HARLEM implementation
    Dialogs created without wxDesigner

    \author Igor Kurnikov
    \date 2026-
*/

#include "hampi.h"

#include "wx/wx.h"
#include "wx/listctrl.h"

#include "hamolset.h"
#include "hamolview.h"
#include "hamolecule.h"
#include "object3d.h"
#include "ha_wx_res_wdr.h"

#include "dialogs_wx_3.h"


// ============================================================
//  Object3DDlgWX  –  "Display Molecules" dialog
// ============================================================

int Object3DDlgWX::dlg_open = FALSE;

Object3DDlgWX::Object3DDlgWX(HaMolView* new_pview, wxWindow *parent):
	wxDialog( parent, -1, "Display Molecules", wxDefaultPosition, wxDefaultSize,
			  wxDEFAULT_DIALOG_STYLE | wxRESIZE_BORDER )
{
	pview = new_pview;
	dlg_open = TRUE;
	obj_list_ctrl = NULL;
	CreateControls();
}

Object3DDlgWX::~Object3DDlgWX()
{
	dlg_open = FALSE;
}

void Object3DDlgWX::CreateControls()
{
	wxBoxSizer *topSizer = new wxBoxSizer( wxVERTICAL );

	wxBoxSizer *mainSizer = new wxBoxSizer( wxHORIZONTAL );

	// Left side: list + name edit
	wxBoxSizer *leftSizer = new wxBoxSizer( wxVERTICAL );

	leftSizer->Add( new wxStaticText(this, wxID_ANY, "Molecules and Objects:"),
		0, wxALIGN_LEFT | wxALL, 5 );

	obj_list_ctrl = new wxListCtrl( this, IDC_OBJ3D_OBJ_LIST,
		wxDefaultPosition, wxSize(280, 200),
		wxLC_REPORT | wxLC_NO_HEADER );
	obj_list_ctrl->InsertColumn(0, "Name", wxLIST_FORMAT_LEFT, 260);
	leftSizer->Add( obj_list_ctrl, 1, wxGROW | wxALL, 5 );

	leftSizer->Add( new wxTextCtrl(this, IDC_OBJ3D_NAME, "",
		wxDefaultPosition, wxSize(80, -1)),
		0, wxGROW | wxALL, 5 );

	mainSizer->Add( leftSizer, 1, wxGROW | wxALL, 0 );

	// Right side: buttons
	wxBoxSizer *btnSizer = new wxBoxSizer( wxVERTICAL );

	btnSizer->Add( new wxButton(this, IDC_OBJ3D_DELETE, "Delete"),
		0, wxALIGN_CENTER | wxALL, 5 );
	btnSizer->Add( new wxButton(this, IDC_OBJ3D_DISPLAY, "Display"),
		0, wxALIGN_CENTER | wxALL, 5 );
	btnSizer->Add( new wxButton(this, IDC_OBJ3D_UNDISPLAY, "Undisplay"),
		0, wxALIGN_CENTER | wxALL, 5 );
	btnSizer->Add( new wxButton(this, IDC_OBJ3D_SET_TRANSP, "Set Transparency"),
		0, wxALIGN_CENTER | wxALL, 5 );

	btnSizer->Add( new wxStaticText(this, wxID_ANY, "Transparency:"),
		0, wxALIGN_CENTER | wxALL, 5 );
	btnSizer->Add( new wxTextCtrl(this, IDC_OBJ3D_TRANSP, "",
		wxDefaultPosition, wxSize(80, -1)),
		0, wxALIGN_CENTER | wxALL, 5 );

	mainSizer->Add( btnSizer, 0, wxALIGN_CENTER | wxALL, 5 );

	topSizer->Add( mainSizer, 1, wxGROW | wxALL, 5 );

	topSizer->Add( new wxButton(this, IDC_OBJ3D_UPDATE, "Update Object List"),
		0, wxALIGN_CENTER_VERTICAL | wxALL, 5 );

	SetSizer( topSizer );
	topSizer->SetSizeHints( this );

	// Bind events (no event table needed)
	Bind(wxEVT_BUTTON, &Object3DDlgWX::OnSetTransp, this, IDC_OBJ3D_SET_TRANSP);
	Bind(wxEVT_BUTTON, &Object3DDlgWX::OnDelete,    this, IDC_OBJ3D_DELETE);
	Bind(wxEVT_BUTTON, &Object3DDlgWX::OnDisplay,   this, IDC_OBJ3D_DISPLAY);
	Bind(wxEVT_BUTTON, &Object3DDlgWX::OnUnDisplay,  this, IDC_OBJ3D_UNDISPLAY);
	Bind(wxEVT_BUTTON, &Object3DDlgWX::OnUpdate,    this, IDC_OBJ3D_UPDATE);
	Bind(wxEVT_CLOSE_WINDOW, &Object3DDlgWX::OnClose, this);
}

bool
Object3DDlgWX::TransferDataToWindow()
{
	if(pview != NULL)
	{
		DDX_obj_list();
	}
	wxDialog::TransferDataToWindow();
	return true;
}


void
Object3DDlgWX::DDX_obj_list()
{
	if(pview == NULL) return;
	MolSet* pmset = pview->GetMolSet();
	if(pmset == NULL) return;

	obj_list_ctrl->DeleteAllItems();
	obj_vec.clear();

	long idx = 0;

	// Add molecules
	for(int i = 0; i < pmset->GetNMol(); i++)
	{
		HaMolecule* pmol = pmset->GetMolByIdx(i);
		obj_list_ctrl->InsertItem(idx, pmol->GetObjName());
		if( !pmol->IsDisplayed() )
		{
			obj_list_ctrl->SetItemBackgroundColour(idx, wxColour(220, 220, 220));
		}
		obj_vec.push_back(pmol);
		idx++;
	}

	// Add other 3D objects
	std::list<Object3D*>::iterator oitr;
	for(oitr = (pmset->ViewObjects).begin(); oitr != (pmset->ViewObjects).end(); oitr++)
	{
		obj_list_ctrl->InsertItem(idx, (*oitr)->GetObjName());
		if( !(*oitr)->IsDisplayed() )
		{
			obj_list_ctrl->SetItemBackgroundColour(idx, wxColour(220, 220, 220));
		}
		obj_vec.push_back(*oitr);
		idx++;
	}
}

void Object3DDlgWX::OnClose(wxCloseEvent &event )
{
	event.Skip();
	dlg_open = FALSE;
}

void Object3DDlgWX::OnSetTransp(wxCommandEvent &event )
{
	double transp;
	wxTextCtrl* edit_transp = (wxTextCtrl*) FindWindow(IDC_OBJ3D_TRANSP);
	wxString str = edit_transp->GetValue();
	int ires = str.ToDouble(&transp);

	if( !ires) return;
	if(pview == NULL)
		return;

	MolSet* pmset = pview->GetMolSet();
	if(pmset == NULL)
		return;

	long item = -1;
	while( (item = obj_list_ctrl->GetNextItem(item, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != -1 )
	{
		obj_vec[item]->SetTransparency(transp);
	}

	TransferDataToWindow();
	pmset->RefreshAllViews(RFRefresh);
}

void Object3DDlgWX::OnDelete(wxCommandEvent &event )
{
	if(pview == NULL)
		return;
	MolSet* pmset = pview->GetMolSet();
	if(pmset == NULL)
		return;

	long item = -1;
	while( (item = obj_list_ctrl->GetNextItem(item, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != -1 )
	{
		wxString obj_name = obj_list_ctrl->GetItemText(item);
		pmset->DeleteObject3D( obj_name.ToStdString() );
	}
	TransferDataToWindow();
	pmset->RefreshAllViews(RFRefresh);
}

void Object3DDlgWX::OnDisplay(wxCommandEvent &event )
{
	if(pview == NULL)
		return;
	MolSet* pmset = pview->GetMolSet();
	if(pmset == NULL)
		return;

	long item = -1;
	while( (item = obj_list_ctrl->GetNextItem(item, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != -1 )
	{
		obj_vec[item]->SetDisplayed(true);
	}

	TransferDataToWindow();
	pmset->RefreshAllViews(RFRefresh);
}

void Object3DDlgWX::OnUnDisplay(wxCommandEvent &event )
{
	if(pview == NULL)
		return;
	MolSet* pmset = pview->GetMolSet();
	if(pmset == NULL)
		return;

	long item = -1;
	while( (item = obj_list_ctrl->GetNextItem(item, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != -1 )
	{
		obj_vec[item]->SetDisplayed(false);
	}

	TransferDataToWindow();
	pmset->RefreshAllViews(RFRefresh);
}

void Object3DDlgWX::OnUpdate(wxCommandEvent &event )
{
	TransferDataToWindow();
}
