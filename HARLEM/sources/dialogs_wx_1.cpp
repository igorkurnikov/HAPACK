/*! \file dialogs_wx_1.cpp

    wxWindows Dialogs 1 in HARLEM implementation
 
    \author Igor Kurnikov  
    \date 2003-
*/

#define DIALOGS_WX_1_CPP

#include <mpi.h>

#include "vec3d.h"
#include "g94_globals.h"
#include "g94_protos.h"

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "wx/valgen.h"
#include "wx/grid.h"
#include "wx/dir.h"
#include "wx/filename.h"
#include "wx/listctrl.h"
#include "wx/tokenzr.h"
#include "wx/notebook.h"
#include "wx/evtloop.h"
#include "wx/spinbutt.h"

#ifdef HAOGL
#include <wx/glcanvas.h>
#endif

#include "hampi.h"

#include "rapidxml.hpp"

#include "hamainframe_wx.h"
#include "ctrl_wx.h"
#include "ha_wx_aux_1.h"
#include "canvas3d.h"
#include "dialogs_wx_1.h"
#include "dialogs_wx_2.h"
#include "qc_dialogs_wx.h"
#include "ha_wx_res_wdr.h"
#include "abstree.h"
#include "halocorb.h"
#include "hacoord.h"
#include "haintcrd.h"
#include "haqchem.h"
#include "etcoupl.h"
#include "electrostmod.h"
#include "protonredox.h"
#include "hagaussian.h"
#include "hazindo.h"
#include "hadalton.h"
#include "harlemapp.h"
#include "hamolset.h"
#include "moleditor.h"
#include "hamolview.h"
#include "hamolecule.h"
#include "mm_params.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "hamolmech.h"
#include "mm_force_field.h"
#include "mm_traj_anal.h"
#include "mm_driver_amber.h"
#include "mm_driver_tinker.h"
#include "haintermol.h"
#include "haresdb.h"
#include "nuclacidmod.h"
#include "hascattermod.h"
#include "haempirical.h"
#include "edit_mut_map_dlg_wx.h"

#include "gaufile.h"
//>mikola to make some wxstuff to work under win
//#include "planview.h"
//#include "wxpnp_wdr.h"
#include <wx/dcbuffer.h>
#include <wx/wx.h>
#include <wx/spinctrl.h>
//#include <wxPlot/wxPlot.h>
#include <wx/scrolwin.h>
//<mikola
//#include "mysql.h"

#if defined(__WXMSW__)
// undefined MFC MACROS conflicting WX functions
#ifdef DrawText
#undef DrawText
#endif

#ifdef StartDoc
#undef StartDoc
#endif

#ifdef GetCharWidth
#undef GetCharWidth
#endif

#ifdef FindWindow
#undef FindWindow
#endif

#endif

#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/format.hpp>

CmdTextCtrl::CmdTextCtrl(wxWindow* parent, wxWindowID id, const wxString& value,
		const wxPoint& pos,const wxSize& size, long style, const wxValidator& validator,
		const wxString& name):
	wxTextCtrl(parent,id,value, pos, size,style,validator,name)
{

}

BEGIN_EVENT_TABLE(CmdTextCtrl, wxTextCtrl)
   EVT_KEY_DOWN(CmdTextCtrl::OnChar)
END_EVENT_TABLE()

void CmdTextCtrl::OnChar(wxKeyEvent& event)
{
	int code = event.GetKeyCode();
	std::string cmd_line;

	switch(code)
	{
	case(WXK_UP): // Arrow UP
		cmd_line= pApp->cmd_pr.RollHistory(1);
		if(!cmd_line.empty())
			this->SetValue(cmd_line.c_str());
		break;
	case( WXK_DOWN): // Arrow DOWN
		cmd_line= pApp->cmd_pr.RollHistory(-1);
		if(!cmd_line.empty())
			this->SetValue(cmd_line.c_str());
		break;
	default:
		event.Skip();
		break;
	}
}


//----------------------------------------------------------------------------
// HaFileDlg1
//----------------------------------------------------------------------------

// WDR: event table for HaFileDlg1
BEGIN_EVENT_TABLE(HaFileDlg1,wxDialog)
    EVT_BUTTON( IDC_CHOOSE_DIR, HaFileDlg1::ChooseDir )
    EVT_LIST_ITEM_SELECTED( IDC_FILE_LIST, HaFileDlg1::OnSelectFile )
	EVT_LIST_ITEM_ACTIVATED( IDC_FILE_LIST, HaFileDlg1::OnActivateFile )
    EVT_CHOICE( IDC_FILE_TYPE, HaFileDlg1::OnChangeFileType )
	EVT_CLOSE(HaFileDlg1::OnClose)
END_EVENT_TABLE()

HaFileDlg1::HaFileDlg1( wxWindow *parent, wxWindowID id, const wxString &title,
    const wxPoint &position, const wxSize& size, long style ) :
    wxDialog( parent, id, title, position, size, style )
{
    nsubdir = 0;
	file_types_ch = NULL; 
}

//HaFileDlg1::~HaFileDlg1()
//{
//	//wxDialog::~wxDialog();
//}

void HaFileDlg1::OnClose(wxCloseEvent& event)
{
	event.Skip();
//	this->EndModal(0);
//	this->Close();
}

void HaFileDlg1::OnChangeFileType( wxCommandEvent &event )
{
    OnChangeDir();
}

void HaFileDlg1::OnSelectFile( wxListEvent &event )
{
	wxListCtrl* file_list = (wxListCtrl*)FindWindow(IDC_FILE_LIST); 
	wxTextCtrl* file_name_edt = (wxTextCtrl*)FindWindow(IDC_FILE_NAME);
	file_name_edt->SetValue("");
	long idx= -1;
	while ((idx = file_list->GetNextItem(idx, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED)) != wxNOT_FOUND) 
	{
		if (idx < nsubdir) continue;
		wxString fname = file_list->GetItemText(idx);
		wxString fname_old = file_name_edt->GetValue();
		if (!fname_old.empty())
		{
			fname = fname_old + ";" + fname;
		}
		file_name_edt->SetValue(fname);
		PrintLog("%s\n",file_list->GetItemText(idx).ToStdString().c_str());
	}
}

void HaFileDlg1::OnActivateFile( wxListEvent &event )
{
	int idx = event.GetIndex();
	if( idx < nsubdir)
	{
		wxString sub_dir_name = event.GetText();
	    wxTextCtrl* dir_name = (wxTextCtrl*) FindWindow( IDC_DIR_NAME );
       
		wxFileName path( dir_name->GetValue(), "" );
		path.AppendDir( sub_dir_name );
        path.Normalize();

	    wxString path_str;
		path_str = path.GetPath();

	    dir_name->SetValue(path_str);
	
	    OnChangeDir();
	}
}

void HaFileDlg1::ChooseDir( wxCommandEvent &event )
{
    wxTextCtrl* dir_name = (wxTextCtrl*) FindWindow( IDC_DIR_NAME );
    wxString path = dir_name->GetValue();
    path = ::wxDirSelector("Choose a directory",path);
    dir_name->SetValue(path);
     
    OnChangeDir();
}

void HaFileDlg1::OnChangeDir()
{
    int isel = file_types_ch->GetSelection();
	if (isel == wxNOT_FOUND) return;
    wxString filter_str = (const char*) file_types_ch->GetClientData(isel);

    wxString path;
    wxTextCtrl* dir_name = (wxTextCtrl*) FindWindow( IDC_DIR_NAME );
    path = dir_name->GetValue();

    wxArrayString file_names;
	wxArrayString filters;
	wxStringTokenizer tkz(filter_str,";");
	while ( tkz.HasMoreTokens())
	{
		wxString token = tkz.GetNextToken();
		filters.Add(token);
	}

	wxArrayString tmp_arr;
	wxArrayString dir_arr;
	
	wxString dir_str;

	wxDir cur_dir(path);
//	wxImageList dir_img(16, 16, true);
//	dir_img.Add( wxIcon(_T("iconsmall"), wxBITMAP_TYPE_ICO_RESOURCE) );

//	wxIcon ico_dir;
//  ico_dir.CopyFromBitmap(MainToolBarBitmaps(0));
//	ico_dir.CopyFromBitmap(DirsBitmap(0));
//	dir_img.Add(ico_dir);

	bool bres = cur_dir.GetFirst(&dir_str, wxEmptyString, wxDIR_DIRS | wxDIR_DOTDOT );
	if( bres && dir_str != ".") dir_arr.Add(dir_str);
	while( cur_dir.GetNext(&dir_str) )
	{
		if(dir_str == ".") continue;
		dir_arr.Add(dir_str);
	}

	nsubdir = dir_arr.GetCount();

	int nfilters = filters.Count();
	int i;
	for(i =0; i < nfilters; i++)
	{
		tmp_arr.Clear();
		wxDir::GetAllFiles(path,&tmp_arr,filters[i], wxDIR_FILES  );

		int nfiles = tmp_arr.Count();
		int j;
		for(j = 0; j < nfiles; j++)
		{
		   file_names.Add(tmp_arr[j]);
		}
	}
	
    wxListCtrl* file_list = (wxListCtrl*) FindWindow( IDC_FILE_LIST);
    if(file_list->GetItemCount() > 0 ) file_list->ClearAll();
//	file_list->SetImageList(&dir_img,wxIMAGE_LIST_SMALL);

    int nd = dir_arr.GetCount();

    for(i=0; i < nsubdir; i++)
    {
		file_list->InsertItem(i, dir_arr[i],0);
		wxListItem item;
        item.m_itemId = i;
		wxFont item_font = item.GetFont();
#if defined _MSC_VER
		item_font.SetWeight(wxFONTWEIGHT_BOLD);
#endif
		item.SetFont(item_font);
        item.SetTextColour(*wxBLUE);
        file_list->SetItem( item );
	}

	int nf = file_names.GetCount();
    wxFileName fname;

    for(i=0; i < nf; i++)
    {
       fname.Assign(file_names[i]);
       file_list->InsertItem(i+nd, fname.GetFullName());
    }
	
}

void HaFileDlg1::FillFileTypes()
{

}

const int IDC_LOAD_FILE = 31001;

//----------------------------------------------------------------------------
// ChooseMolFileDlg
//----------------------------------------------------------------------------

// WDR: event table for ChooseMolFileDlg
BEGIN_EVENT_TABLE(ChooseMolFileDlg,HaFileDlg1)
    EVT_BUTTON( IDC_LOAD_FILE, ChooseMolFileDlg::OnLoadFile )
	EVT_CLOSE( ChooseMolFileDlg::OnClose )
END_EVENT_TABLE()

ChooseMolFileDlg::ChooseMolFileDlg( wxWindow *parent, wxWindowID id, const wxString &title,
    const wxPoint &position, const wxSize& size, long style ) :
    HaFileDlg1( parent, id, title, position, size, style )
{
	 sizer_main_v = new wxBoxSizer( wxVERTICAL );
     wxBoxSizer* sizer_h_dir = new wxBoxSizer( wxHORIZONTAL );

     wxStaticText *item2 = new wxStaticText( this, ID_TEXT, wxT("Directory:"), wxDefaultPosition, wxDefaultSize, 0 );
     sizer_h_dir->Add( item2, 0, wxALIGN_CENTER|wxALL, 5 );

     dir_name_txt = new wxTextCtrl( this, IDC_DIR_NAME, wxT(""), wxDefaultPosition, wxSize(250,-1), 0 );
	 wxString cur_dir = ::wxGetCwd();
     dir_name_txt->SetValue(cur_dir);

     sizer_h_dir->Add( dir_name_txt, 0, wxALIGN_CENTER|wxALL, 5 );

     wxButton *item4 = new wxButton( this, IDC_CHOOSE_DIR, wxT("Choose Dir"), wxDefaultPosition, wxDefaultSize, 0 );
     sizer_h_dir->Add( item4, 0, wxALIGN_CENTER|wxALL, 5 );

     sizer_main_v->Add( sizer_h_dir, 0, wxALL, 5 );

     wxListCtrl *item5 = new wxListCtrl( this, IDC_FILE_LIST, wxDefaultPosition, wxSize(400,120), wxLC_LIST| wxLC_NO_HEADER | wxSUNKEN_BORDER );
     sizer_main_v->Add( item5, 0, wxALIGN_CENTER|wxALL, 5 );

     wxBoxSizer *item6 = new wxBoxSizer( wxHORIZONTAL );

     wxStaticText *item7 = new wxStaticText( this, ID_TEXT, wxT("File Name"), wxDefaultPosition, wxDefaultSize, 0 );
     item6->Add( item7, 0, wxALIGN_CENTER|wxALL, 5 );

     file_name_txt = new wxTextCtrl( this, IDC_FILE_NAME, wxT(""), wxDefaultPosition, wxSize(250,-1), 0 );
     item6->Add( file_name_txt, 0, wxALIGN_CENTER|wxALL, 5 );

     choose_file_btn = new wxButton( this, IDC_LOAD_FILE, "Load File", wxDefaultPosition, wxDefaultSize, 0 );
     item6->Add( choose_file_btn, 0, wxALIGN_CENTER|wxALL, 5 );

     sizer_main_v->Add( item6, 0, wxALL, 5 );

     wxBoxSizer *item10 = new wxBoxSizer( wxHORIZONTAL );

     wxStaticText *item11 = new wxStaticText( this, ID_TEXT, wxT("File Type"), wxDefaultPosition, wxDefaultSize, 0 );
     item10->Add( item11, 0, wxALIGN_CENTER|wxALL, 5 );

     file_types_ch = new wxChoice( this, IDC_FILE_TYPE, wxDefaultPosition, wxDefaultSize, 0, NULL, 0 );
	 FillFileTypes();
     item10->Add( file_types_ch, 0, wxGROW|wxALL, 5 );

     sizer_main_v->Add( item10, 0, wxALL, 5 );

 //    wxBoxSizer *item16 = new wxBoxSizer( wxHORIZONTAL );

 //    wxStaticText *item17 = new wxStaticText( this, ID_TEXT, wxT("Molecule Name"), wxDefaultPosition, wxDefaultSize, 0 );
 //    item16->Add( item17, 0, wxALIGN_CENTER|wxALL, 5 );

 //    wxTextCtrl *item18 = new wxTextCtrl( this, IDC_OFDLG_MOLNAME, wxT(""), wxDefaultPosition, wxSize(80,-1), 0 );
 //    item16->Add( item18, 0, wxALIGN_CENTER|wxALL, 5 );

 //    sizer_main_v->Add( item16, 0, wxALIGN_CENTER_VERTICAL|wxALL, 5 );

     this->SetSizer( sizer_main_v );
     sizer_main_v->SetSizeHints( this );

     OnInitDialog();                          
}

//ChooseMolFileDlg::~ChooseMolFileDlg()
//{
//	HaFileDlg1::~HaFileDlg1();
//}

const char* pdb_filters = "*.pdb;*.ent";
const char* sybyl_filters = "*.syb;*.mol;*.mol2";
const char* mdl_filters = "*.mdl;*.mol";
const char* xyz_filters = "*.xyz";
const char* hin_filters = "*.hin";
const char* gauss_filters = "*.rwf;*.chk";
const char* harlem_filters = "*.hlm";
const char* amber_prep_filters = "*.in;*.prep";
const char* amber_top_filters = "*.top;*.crd;*.rst";
const char* amber_off_filters = "*.off;*.lib";
const char* nrg_filters = "*.nrg";

void ChooseMolFileDlg::FillFileTypes()
{
    file_types_ch->Clear();

    file_types_ch->Append("Protein Data Bank  (*.pdb;*.ent)", (void*)pdb_filters);
    file_types_ch->Append("TRIPOS MOL2 Format  (*.syb;*.mol;*.mol2)",(void*)sybyl_filters);
    file_types_ch->Append("MDL Mol File Format(*.mdl;*.mol)",(void*)mdl_filters);
    file_types_ch->Append("TINKER XYZ Format    (*.xyz)",(void*)xyz_filters);
    file_types_ch->Append("Gaussian RWF File Format (*.rwf;*.chk)",(void*)gauss_filters);
	file_types_ch->Append("HARLEM File Format (*.hlm)",(void*)harlem_filters);
    file_types_ch->Append("AMBER PREP File Format (*.in;*.prep)",(void*)amber_prep_filters);
	file_types_ch->Append("AMBER TOP File Format (*.top)",(void*)amber_top_filters);
	file_types_ch->Append("AMBER OFF File Format (*.off;*.lib)", (void*)amber_off_filters);
	file_types_ch->Append("ARBALEST File Format (*.hin)",(void*)hin_filters);
	file_types_ch->Append("MBX Program File Format (*.nrg)", (void*)nrg_filters);

    file_types_ch->SetSelection(5);
}

void ChooseMolFileDlg::OnInitDialog()
{ 	
	file_format = FormatHarlem;

	OnChangeDir();    
}

// WDR: handler implementations for ChooseMolFileDlg

void ChooseMolFileDlg::OnLoadFile( wxCommandEvent &event )
{
//	PrintLog(" ChooseMolFileDlg::OnLoadFile() pt 1 \n");
    
	file_name = file_name_txt->GetValue().ToStdString();
    dir_name = dir_name_txt->GetValue().ToStdString();
   
    int isel = file_types_ch->GetSelection();

	switch( isel )
    {
        case(0): file_format = FormatPDB;     break;
        case(1): file_format = FormatMol2;    break;
        case(2): file_format = FormatMDL;     break;
        case(3): file_format = FormatXYZ;     break;
        case(4): file_format = FormatRWF;     break;
		case(5): file_format = FormatHarlem;  break;
        case(6): file_format = FormatAmberPrep;  break;
		case(7): file_format = FormatAmberTop;   break;
		case(8): file_format = FormatAmberOff;   break;
		case(9): file_format = FormatHIN;      break;
		case(10): file_format = FormatNRG;      break;
    }
	this->Close();
}

void ChooseMolFileDlg::OnClose(wxCloseEvent& event)
{
	event.Skip();
}

//----------------------------------------------------------------------------
// SaveMolFileDlg
//----------------------------------------------------------------------------

// WDR: event table for SaveMolFileDlg
BEGIN_EVENT_TABLE(SaveMolFileDlg,HaFileDlg1)
    EVT_BUTTON( IDC_SAVE_FILE, SaveMolFileDlg::OnSaveFile )
END_EVENT_TABLE()

SaveMolFileDlg::SaveMolFileDlg( MolSet* pmset_new, wxWindow *parent, wxWindowID id, const wxString &title,
    const wxPoint &position, const wxSize& size, long style ) :
    HaFileDlg1( parent, id, title, position, size, style )
{
    pmset = pmset_new;

    save_mol_file_dlg( this, TRUE );

	file_types_ch = (wxChoice*) FindWindow(IDC_FILE_TYPE);
	wxCheckBox* chk_box;
	chk_box = (wxCheckBox*) FindWindow( IDC_SAVE_TRANSFORMED );
	if(pmset->save_opt_default.save_transform)
		chk_box->SetValue(true);
	chk_box  = (wxCheckBox*) FindWindow( IDC_SAVE_CONNECT );
	if(pmset->save_opt_default.save_connect)
		chk_box->SetValue(true);
	chk_box  = (wxCheckBox*) FindWindow( IDC_SAVE_AMBER_PDB );
	if(pmset->save_opt_default.save_amber_pdb)
		chk_box->SetValue(true);
	chk_box = (wxCheckBox*)FindWindow(IDC_SAVE_SEP_WAT_MOL);
	if (pmset->save_opt_default.save_sep_solv_mol)
		chk_box->SetValue(true);

    FillFileTypes();
    wxString cur_dir = ::wxGetCwd();
    wxTextCtrl* dir_name = (wxTextCtrl*) FindWindow( IDC_DIR_NAME );
    dir_name->SetValue(cur_dir);
     
    OnChangeDir(); 
}

const char* bmp_filters = "*.bmp";
const char* gif_filters = "*.gif";
const char* pcx_filters = "*.pcx";
const char* png_filters = "*.png";
const char* raw_port_filters   = "*.ppm";
const char* ascii_port_filters = "*.ppm";
const char* pict_filters = "*.pic;*.pict";
const char* rgb_filters  = "*.rgb";
const char* jpeg_filters = "*.jpg;*.jpeg";
const char* tiff_filters = "*.tiff;*.tif";
const char* ps_filters   = "*.ps;*.eps";
const char* pov_filters  = "*.pov";
const char* vrml_filters = "*.vrml;*.wrl";

bool SaveMolFileDlg::TransferDataToWindow()
{
	return wxDialog::TransferDataToWindow();
}

bool SaveMolFileDlg::TransferDataFromWindow()
{
	return wxDialog::TransferDataFromWindow();
}


void SaveMolFileDlg::FillFileTypes()
{ 
    file_types_ch->Clear();

    file_types_ch->Append("Protein Data Bank  (*.pdb;*.ent)", (void*)pdb_filters);
    file_types_ch->Append("HARLEM File Format (*.hlm)",(void*)harlem_filters);
	file_types_ch->Append("OLD HARLEM File Format (*.hlm)",(void*)harlem_filters);
    file_types_ch->Append("TINKER XYZ Format    (*.xyz)",(void*)xyz_filters);
    file_types_ch->Append("MSMS XYZR Format   (*.xyz)",(void*)xyz_filters);
	file_types_ch->Append("ARBALEST HIN Format   (*.hin)",(void*)hin_filters);
	file_types_ch->Append("MBX Program NRG Format   (*.nrg)", (void*)nrg_filters);

    file_types_ch->SetSelection(1);
}

// WDR: handler implementations for SaveMolFileDlg

void SaveMolFileDlg::OnSaveFile( wxCommandEvent &event )
{
	wxCheckBox* chk_box;

	chk_box = (wxCheckBox*) FindWindow( IDC_SAVE_TRANSFORMED );
	pmset->save_opt_default.save_transform = chk_box->IsChecked() ? TRUE : FALSE;
	
	chk_box   = (wxCheckBox*) FindWindow( IDC_SAVE_CONNECT );
	pmset->save_opt_default.save_connect = chk_box->IsChecked() ? TRUE : FALSE;

	chk_box   = (wxCheckBox*) FindWindow( IDC_SAVE_AMBER_PDB );
	pmset->save_opt_default.save_amber_pdb = chk_box->IsChecked() ? TRUE : FALSE;

	chk_box = (wxCheckBox*)FindWindow(IDC_SAVE_SEP_WAT_MOL);
	pmset->save_opt_default.save_sep_solv_mol = chk_box->IsChecked() ? TRUE : FALSE;
		
    wxTextCtrl* file_name_edt =  (wxTextCtrl*) FindWindow( IDC_FILE_NAME );
    wxString file_name = file_name_edt->GetValue();

    wxTextCtrl* dir_name_edt = (wxTextCtrl*) FindWindow( IDC_DIR_NAME );
    wxString dir_name = dir_name_edt->GetValue();

    wxFileName fname_full;
    fname_full.Assign(dir_name,file_name);
 	
	if(!file_name.empty() && pmset != NULL )
	{
       std::string file_name_full = fname_full.GetFullPath().ToStdString();
       ::wxSetWorkingDirectory(dir_name);

        int isel = file_types_ch->GetSelection();

		switch( isel )
		{
			case(0): pmset->SavePDBFile(file_name_full);     break;
			case(1): pmset->SaveHarlemFile(file_name_full);  break;
			case(2): pmset->SaveOldHarlemFile(file_name_full);  break;
			case(3): pmset->SaveXYZFile(file_name_full);     break;
			case(4): pmset->SaveXYZRadFile(file_name_full);  break;	
			case(5): pmset->SaveHINFile(file_name_full);  break;	
			case(6): pmset->SaveNRGFile(file_name_full);  break;
		}		
	}
	delete this;
}


//----------------------------------------------------------------------------
// SaveImageFileDlg
//----------------------------------------------------------------------------

// WDR: event table for SaveImageFileDlg
BEGIN_EVENT_TABLE(SaveImageFileDlg,HaFileDlg1)
    EVT_BUTTON( IDC_SAVE_FILE, SaveImageFileDlg::OnSaveFile )
END_EVENT_TABLE()

SaveImageFileDlg::SaveImageFileDlg( HaMolView* pview, wxWindow *parent, wxWindowID id, const wxString &title,
    const wxPoint &position, const wxSize& size, long style ) :
    HaFileDlg1( parent, id, title, position, size, style )
{
	mol_view = pview;
    save_image_file_dlg( this, TRUE );

	file_types_ch = (wxChoice*) FindWindow(IDC_FILE_TYPE);
    FillFileTypes();
    wxString cur_dir = ::wxGetCwd();
    wxTextCtrl* dir_name = (wxTextCtrl*) FindWindow( IDC_DIR_NAME );
    dir_name->SetValue(cur_dir);
     
    OnChangeDir();                              
}


void SaveImageFileDlg::FillFileTypes()
{    
    file_types_ch->Clear();

    file_types_ch->Append("Microsoft Bitmap  (*.bmp)", (void*)bmp_filters);
    file_types_ch->Append("CompuServe GIF  (*.gif)", (void*)gif_filters);
    file_types_ch->Append("PCX  ", (void*)pcx_filters);
    file_types_ch->Append("PNG  ", (void*)png_filters);
    file_types_ch->Append("Raw Portable Pixmap  (*.ppm)", (void*)raw_port_filters);
    file_types_ch->Append("ASCII Portable Pixmap  (*.ppm)", (void*)ascii_port_filters);
    file_types_ch->Append("Apple Macintosh PICT  (*.pic;*.pict)", (void*)pict_filters);
    file_types_ch->Append("Silicon Graphics RGB  (*.rgb)", (void*)rgb_filters);
    file_types_ch->Append("JPEG  (*.jpg;*.jpeg)", (void*)jpeg_filters);
    file_types_ch->Append("TIFF  (*.tiff;*.tif)", (void*)tiff_filters);
    file_types_ch->Append("Vector PostScript  (*.ps;*.eps)", (void*)ps_filters);
    file_types_ch->Append("POV RAY (*.pov)", (void*)pov_filters);
    file_types_ch->Append("VRML  (*.vrml;*.wrl)", (void*)vrml_filters);

    file_types_ch->SetSelection(8);
}


// WDR: handler implementations for SaveImageFileDlg

void SaveImageFileDlg::OnSaveFile( wxCommandEvent &event )
{
	TransferDataFromWindow();
    wxTextCtrl* file_name_edt =  (wxTextCtrl*) FindWindow( IDC_FILE_NAME );
    wxString file_name = file_name_edt->GetValue();

    wxTextCtrl* dir_name_edt = (wxTextCtrl*) FindWindow( IDC_DIR_NAME );
    wxString dir_name = dir_name_edt->GetValue();

    wxFileName fname_full;
    fname_full.Assign(dir_name,file_name);
 
	
	if(!file_name.empty() && mol_view != NULL )
	{
       wxString file_name_full = fname_full.GetFullPath();

        int isel = file_types_ch->GetSelection();

		switch( isel )
		{
		case(0): mol_view->WriteBMPFile(file_name_full.c_str());       break;
		case(1): mol_view->WriteGIFFile(file_name_full.c_str());       break;
		case(2): mol_view->WritePCXFile(file_name_full.c_str());       break;
		case(3): mol_view->WritePNGFile(file_name_full.c_str());       break;		
		case(4): mol_view->WritePPMFile(file_name_full.c_str(),True);  break;		
		case(5): mol_view->WritePPMFile(file_name_full.c_str(),False); break;		
		case(6): mol_view->WritePICTFile(file_name_full.c_str());  break;		
		case(7): mol_view->WriteIRISFile(file_name_full.c_str());  break;		
		case(8): mol_view->WriteJPEGFile(file_name_full.c_str());  break;		
		case(9): mol_view->WriteTIFFFile(file_name_full.c_str());  break;		
		case(10): mol_view->WriteVectPSFile(file_name_full.c_str());  break;		
		case(11): mol_view->WritePOVRayFile(file_name_full.c_str());  break;		
		case(12): mol_view->WriteVRMLFile(file_name_full.c_str());  break;			
		}		
	}
	delete this;
}

int PrintMessage(const char* str)
{
	HaMainFrameWX* frame_main = GetHaMainFrameWX();
	if(frame_main)
	{
		wxTextCtrl* log_win = (wxTextCtrl*) frame_main->FindWindow(IDC_LOG_WIN);
		log_win->SetInsertionPoint(0);
		log_win->WriteText(str);
		log_win->WriteText("\n");
	}
//	printf("%s\n", str);
	PrintLog("%s\n", str); 
	return TRUE;
}

//----------------------------------------------------------------------------
// LoadScriptDlgWX
//----------------------------------------------------------------------------

int LoadScriptDlgWX::dlg_open = False;

// WDR: event table for LoadScriptDlgWX
BEGIN_EVENT_TABLE(LoadScriptDlgWX,HaFileDlg1)
    EVT_BUTTON( IDC_EXEC_SCRIPT_FILE,    LoadScriptDlgWX::OnExecScriptFile )
	EVT_BUTTON( IDC_EXEC_SCRIPT_WIN,     LoadScriptDlgWX::OnExecScriptWindow )
	EVT_BUTTON( IDC_SAVE_SCRIPT_TO_FILE, LoadScriptDlgWX::OnSaveScriptToFile )
	EVT_LIST_ITEM_SELECTED ( IDC_FILE_LIST, LoadScriptDlgWX::OnSelectFile )
	EVT_CLOSE(LoadScriptDlgWX::OnClose)
END_EVENT_TABLE()

LoadScriptDlgWX::LoadScriptDlgWX( wxWindow *parent, wxWindowID id, const wxString &title,
    const wxPoint &position, const wxSize& size, long style ) :
    HaFileDlg1( parent, id, title, position, size, style )
{
	dlg_open = True;
    load_script_dlg( this, TRUE );

	file_types_ch = (wxChoice*) FindWindow(IDC_FILE_TYPE);
    FillFileTypes();
    script_dir = ::wxGetCwd();
    wxTextCtrl* dir_name = (wxTextCtrl*) FindWindow( IDC_DIR_NAME );
    dir_name->SetValue(script_dir);
     
    OnChangeDir();                              
}

const char* python_filters = "*.py";
const char* rasmol_filters = "*.*";

void LoadScriptDlgWX::FillFileTypes()
{  
    file_types_ch->Clear();

    file_types_ch->Append("Python Script (*.py)", (void*)python_filters);
    file_types_ch->Append("RASMOL Script (*.*)",  (void*)rasmol_filters);

    file_types_ch->SetSelection(0);
}

void LoadScriptDlgWX::OnSelectFile(wxListEvent &event )
{
	HaFileDlg1::OnSelectFile(event);

    wxTextCtrl* file_name_edt =  (wxTextCtrl*) FindWindow( IDC_FILE_NAME );
    wxString file_name = file_name_edt->GetValue();

    wxTextCtrl* dir_name_edt = (wxTextCtrl*) FindWindow( IDC_DIR_NAME );
    wxString dir_name = dir_name_edt->GetValue();

	wxTextCtrl* script_text_edt = (wxTextCtrl*) FindWindow( IDC_SCRIPT_TEXT );

	if( file_name.IsEmpty() ) return; 

    wxFileName fname_full;
    fname_full.Assign(dir_name,file_name);

    wxString file_name_full = fname_full.GetFullPath();
	script_text_edt->LoadFile(file_name_full);
}

// WDR: handler implementations for LoadScriptDlgWX

void LoadScriptDlgWX::OnExecScriptFile( wxCommandEvent &event )
{    
    wxTextCtrl* file_name_edt =  (wxTextCtrl*) FindWindow( IDC_FILE_NAME );
    wxString file_name = file_name_edt->GetValue();

    wxTextCtrl* dir_name_edt = (wxTextCtrl*) FindWindow( IDC_DIR_NAME );
    wxString dir_name = dir_name_edt->GetValue();

	wxTextCtrl* script_text_edt = (wxTextCtrl*) FindWindow( IDC_SCRIPT_TEXT );

    wxFileName fname_full;
    fname_full.Assign(dir_name,file_name);

    wxString fname_full_str = fname_full.GetFullPath();
    
    int isel = file_types_ch->GetSelection();

	wxString msg = "Execute Script from file ";
	msg += file_name;
	PrintMessage(msg.c_str());
    
	switch( isel )
    {
        case(0): pApp->ExecuteScriptFromFile(fname_full_str.c_str());   break;
        case(1): pApp->ExecRasMolScript(fname_full_str.c_str());        break;
    }
}

void LoadScriptDlgWX::OnExecScriptWindow( wxCommandEvent &event )
{    
	wxString cur_dir = ::wxGetCwd();
	wxString file_name = "ha_tmp_script.py";
    wxFileName fname_full;
    fname_full.Assign(cur_dir,file_name);
	wxString fname_full_str = fname_full.GetFullPath();
    
	wxTextCtrl* script_text_edt = (wxTextCtrl*) FindWindow( IDC_SCRIPT_TEXT );
    script_text_edt->SaveFile(fname_full_str);

    wxChoice* file_types_ch = (wxChoice*) FindWindow(IDC_FILE_TYPE);
    int isel = file_types_ch->GetSelection();

	wxString msg = "Execute Script from the text window ";
	PrintMessage(msg.c_str());
    
	switch( isel )
    {
        case(0): pApp->ExecuteScriptFromFile(fname_full_str.c_str());   break;
        case(1): pApp->ExecRasMolScript(fname_full_str.c_str());        break;
    }
}

void LoadScriptDlgWX::OnSaveScriptToFile( wxCommandEvent &event )
{
    wxTextCtrl* file_name_edt =  (wxTextCtrl*) FindWindow( IDC_FILE_NAME );
    wxString file_name = file_name_edt->GetValue();

    wxTextCtrl* dir_name_edt = (wxTextCtrl*) FindWindow( IDC_DIR_NAME );
    wxString dir_name = dir_name_edt->GetValue();

    wxFileName fname_full;
    fname_full.Assign(dir_name,file_name);

    wxString fname_full_str = fname_full.GetFullPath();

	wxTextCtrl* script_text_edt = (wxTextCtrl*) FindWindow( IDC_SCRIPT_TEXT );
    script_text_edt->SaveFile(fname_full_str);

	OnChangeDir();
}

void LoadScriptDlgWX::OnClose( wxCloseEvent &event )
{
	event.Skip();
	dlg_open = False;	
}

/////////////////////////////////////////////////////////////////////
//  Object3D Manipulation Dialog:

int Object3DDlgWX::dlg_open = FALSE;

Object3DDlgWX::Object3DDlgWX(HaMolView* new_pview,wxWindow *parent ):
wxDialog( parent, -1, "Manipulate 3D Objects", wxDefaultPosition, wxDefaultSize, 
		   wxDEFAULT_DIALOG_STYLE )
{
	pview= new_pview;
    dlg_open = TRUE;
	object3d_dlg( this, TRUE );
}

Object3DDlgWX::~Object3DDlgWX()
{
     dlg_open = FALSE;
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


BEGIN_EVENT_TABLE(Object3DDlgWX,wxDialog)
    EVT_BUTTON( IDC_OBJ3D_SET_TRANSP,    Object3DDlgWX::OnSetTransp )
    EVT_BUTTON( IDC_OBJ3D_DELETE,        Object3DDlgWX::OnDelete )
    EVT_BUTTON( IDC_OBJ3D_DISPLAY,       Object3DDlgWX::OnDisplay )
    EVT_BUTTON( IDC_OBJ3D_UNDISPLAY,     Object3DDlgWX::OnUnDisplay )
    EVT_BUTTON( IDC_OBJ3D_UPDATE,        Object3DDlgWX::OnUpdate )
    EVT_CLOSE ( Object3DDlgWX::OnClose)
END_EVENT_TABLE()

void
Object3DDlgWX::DDX_obj_list()
{
	if(pview == NULL) return;
 	wxListBox* obj_list = (wxListBox*) FindWindow( IDC_OBJ3D_OBJ_LIST );
	
	MolSet* pmset = pview->GetMolSet();
	if(pmset == NULL) return;

	obj_list->Clear();
	std::list<Object3D*>::iterator oitr;
	obj_vec.clear();
	for(oitr = (pmset->ViewObjects).begin(); oitr != (pmset->ViewObjects).end(); oitr++)
	{
		obj_list->Append((*oitr)->GetObjName());
		obj_vec.push_back(*oitr);
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

	wxListBox* obj_list = (wxListBox*) FindWindow( IDC_OBJ3D_OBJ_LIST );

	int nitem = obj_list->GetCount();
	if(nitem == 0)
		return;
	
	for(int i=0; i < nitem; i++)
	{
		if( obj_list->IsSelected(i))
		{
			obj_vec[i]->SetTransparency(transp);		
		}
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

	wxListBox* obj_list = (wxListBox*) FindWindow( IDC_OBJ3D_OBJ_LIST );
	int nitem = obj_list->GetCount();
	if(nitem <= 0)
		return;
	for(int i=0; i < nitem; i++)
	{
		if( obj_list->IsSelected(i))
		{
			wxString obj_name = obj_list->GetString(i);
			pmset->DeleteObject3D( obj_name.ToStdString() );		
		}
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

	wxListBox* obj_list = (wxListBox*) FindWindow( IDC_OBJ3D_OBJ_LIST );

	int nitem = obj_list->GetCount();
	if(nitem == 0)
		return;
	
	for(int i=0; i < nitem; i++)
	{
		if( obj_list->IsSelected(i))
		{
			obj_vec[i]->SetDisplayed(true);		
		}
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

	wxListBox* obj_list = (wxListBox*) FindWindow( IDC_OBJ3D_OBJ_LIST );

	int nitem = obj_list->GetCount();
	if(nitem == 0)
		return;

	for(int i=0; i < nitem; i++)
	{
		if( obj_list->IsSelected(i))
		{
			obj_vec[i]->SetDisplayed(false);		
		}
	}

	TransferDataToWindow();
	pmset->RefreshAllViews(RFRefresh);		
}



void Object3DDlgWX::OnUpdate(wxCommandEvent &event )
{
	TransferDataToWindow();
}


// SolvateDlgWX Dialog

int SolvateDlgWX::dlg_open = FALSE;

SolvateDlgWX::SolvateDlgWX(MolSet* new_pmset, wxWindow* parent):
wxDialog( parent, -1, "Solvate Molecular Set", wxDefaultPosition, wxDefaultSize, 
		   wxDEFAULT_DIALOG_STYLE )
{
	pmset= new_pmset;
	if( pmset == NULL) return;
	p_mol_editor = pmset->GetMolEditor();
	solvate_mset_dlg(this, TRUE);
    dlg_open = TRUE;
}

SolvateDlgWX::~SolvateDlgWX()
{
    dlg_open = FALSE;
}


void SolvateDlgWX::OnInitDialog(wxInitDialogEvent& event)
{
	dlg_open = TRUE;

	wxTextCtrl* edit_ctrl; 
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_SOLVENT);
	edit_ctrl->SetValidator( StdStringValidator(&p_mol_editor->solv_name) );

	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_SOLV_BUF);
	edit_ctrl->SetValidator( wxDoubleValidator(&p_mol_editor->solv_buffer_dist, "%6.1f") );

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_MIN_SOLUTE_SOLV_DIST);
	edit_ctrl->SetValidator(wxDoubleValidator(&p_mol_editor->min_solute_solv_dist, "%6.1f"));

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_NUM_NA);
	edit_ctrl->SetValidator(wxGenericValidator(&p_mol_editor->num_na_add));

	edit_ctrl = (wxTextCtrl*)FindWindow(IDC_NUM_CL);
	edit_ctrl->SetValidator(wxGenericValidator(&p_mol_editor->num_cl_add));

	event.Skip();
}

bool SolvateDlgWX::TransferDataToWindow()
{	
	return wxDialog::TransferDataToWindow();
}

bool SolvateDlgWX::TransferDataFromWindow()
{
	return wxDialog::TransferDataFromWindow();
}


BEGIN_EVENT_TABLE(SolvateDlgWX,wxDialog)
    EVT_INIT_DIALOG( SolvateDlgWX::OnInitDialog )
    EVT_BUTTON( IDC_SOLV_MOL,    SolvateDlgWX::OnSolvMol )
	EVT_BUTTON( IDC_ADD_IONS,    SolvateDlgWX::OnAddIons )
    EVT_CLOSE( SolvateDlgWX::OnClose)
END_EVENT_TABLE()


void  SolvateDlgWX::OnSolvMol(wxCommandEvent &event) 
{
	TransferDataFromWindow();
	p_mol_editor->Solvate(pmset);
}

void  SolvateDlgWX::OnAddIons(wxCommandEvent& event)
{
	TransferDataFromWindow();
	p_mol_editor->AddIons( pmset, p_mol_editor->num_na_add, p_mol_editor->num_cl_add);
	pmset->RefreshAllViews(RFApply);
}

void SolvateDlgWX::OnClose(wxCloseEvent& event)
{
	SolvateDlgWX::dlg_open = false;
	event.Skip();
}


/////////////////////////////////////////////////////////////////////
// RedirectIODlgWX Dialog:

int RedirectIODlgWX::dlg_open = FALSE;

RedirectIODlgWX::RedirectIODlgWX(wxWindow* parent):
wxDialog( parent, -1, "Redirect IO", wxDefaultPosition, wxDefaultSize, 
		   wxDEFAULT_DIALOG_STYLE )
{
	redirect_io_dlg(this,TRUE);
}

RedirectIODlgWX::~RedirectIODlgWX()
{
        dlg_open = FALSE;
}


void
RedirectIODlgWX::OnInitDialog(wxInitDialogEvent& event)
{
	dlg_open = TRUE;
	wxTextCtrl* edit_log_file = (wxTextCtrl*) FindWindow(IDC_IO_LOG_FILE_NAME);
	edit_log_file->SetValue("harlem.out");	
}


BEGIN_EVENT_TABLE(RedirectIODlgWX,wxDialog)
    EVT_INIT_DIALOG( RedirectIODlgWX::OnInitDialog )
    EVT_BUTTON( IDC_IO_STDOUT_FILE,    RedirectIODlgWX::OnStdoutFile )
	EVT_BUTTON( IDC_IO_STDOUT_CONSOLE,    RedirectIODlgWX::OnStdoutConsole )
    EVT_CLOSE( RedirectIODlgWX::OnClose)
END_EVENT_TABLE()



void RedirectIODlgWX::OnClose(wxCloseEvent& event) 
{
	dlg_open = FALSE;
	event.Skip();
}

void RedirectIODlgWX::OnStdoutFile(wxCommandEvent &event) 
{
	wxTextCtrl* edit_log_file = (wxTextCtrl*) FindWindow(IDC_IO_LOG_FILE_NAME);
	wxString fname = edit_log_file->GetValue();
	RedirectIOToFile(fname.c_str());
}

void RedirectIODlgWX::OnStdoutConsole(wxCommandEvent &event) 
{
	
}


/////////////////////////////////////////////////////////////////////
// PathwaysDlgWX Dialog:

int PathwaysDlgWX::dlg_open = FALSE;

PathwaysDlgWX::PathwaysDlgWX(ETCouplMod* new_etmod, wxWindow* parent):
wxDialog( parent, -1, "ET Pathways Calculations", wxDefaultPosition, wxDefaultSize, 
		   wxDEFAULT_DIALOG_STYLE )
{
	etmod= new_etmod;
	sel_thresh = 0.5;
	pathways_dlg(this,TRUE);
}

PathwaysDlgWX::~PathwaysDlgWX()
{
       dlg_open = FALSE;
}

void
PathwaysDlgWX::OnInitDialog(wxInitDialogEvent& event)
{
	dlg_open = TRUE;

	if(etmod == NULL)
	{
		ErrorInMod("PathwaysDlgWX::OnInitDialog() ", 
		" the pointer to the ET Coupling Module is not valid  ");
		return;
	}

	wxCheckBox* calc_s_hbond_chk = (wxCheckBox*) FindWindow(IDC_PATH_S_HBOND);
    wxCheckBox* hbond_term_chk = (wxCheckBox*) FindWindow(IDC_PATH_HBOND);

	MolSet* pmset = etmod->GetMolSet();
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);

	calc_s_hbond_chk->SetValidator(wxGenericValidator(&p_mol_editor->m_calc_s_hbonds_flag));
	hbond_term_chk->SetValidator(wxGenericValidator(&etmod->m_hbond_paths_flag));
	
	event.Skip();
}

void
PathwaysDlgWX::OnSHBond(wxCommandEvent& event)
{
	MolSet* pmset = etmod->GetMolSet();
		
	pmset->HBonds_found = false;
	
	TransferDataFromWindow();
}

void
PathwaysDlgWX::OnHBondTerm(wxCommandEvent& event)
{
	etmod->pathways_graph_init_flag = false;
	TransferDataFromWindow();
}

void
PathwaysDlgWX::SetModuleDataFromControls(wxCommandEvent &event)
{
	TransferDataFromWindow();
}


bool
PathwaysDlgWX::TransferDataFromWindow()
{
	double old_nb_decay = etmod->pw_nb_decay;
	double old_hb_decay = etmod->pw_hb_decay;
	double old_cov_decay = etmod->pw_ln_cov_decay;
    double old_sel_thresh = this->sel_thresh;

	wxString str;
	bool bres;

	wxTextCtrl* nb_dist_decay = (wxTextCtrl*) FindWindow( IDC_PATH_NB_DECAY );
	str = nb_dist_decay->GetValue();
	bres = str.ToDouble(&etmod->pw_nb_decay);
	if( !bres )
	{
		etmod->pw_nb_decay = old_nb_decay;
		TransferDataToWindow();
		return false;
	}

	wxTextCtrl* hb_dist_decay = (wxTextCtrl*) FindWindow( IDC_PATH_HB_DECAY );
	str = hb_dist_decay->GetValue();
	bres = str.ToDouble(&etmod->pw_hb_decay);
	if( !bres )
	{
		etmod->pw_hb_decay = old_hb_decay;
		TransferDataToWindow();
		return false;
	}

	wxTextCtrl* sel_thresh_edt = (wxTextCtrl*) FindWindow( IDC_PATH_SEL_THRESH );
	str = sel_thresh_edt->GetValue();
	bres = str.ToDouble(&this->sel_thresh);
	if( !bres)
	{
		this->sel_thresh = old_sel_thresh;
		TransferDataToWindow();
		return false;
	}
	
	wxTextCtrl* cov_decay_edit = (wxTextCtrl*) FindWindow(IDC_PATH_COV_DECAY);
	str = cov_decay_edit->GetValue();

	bres = str.ToDouble(&etmod->pw_ln_cov_decay);
	if(!bres)
	{
		etmod->pw_ln_cov_decay = old_cov_decay;
		TransferDataFromWindow();
	}
	etmod->pw_ln_cov_decay = log(fabs(etmod->pw_ln_cov_decay));
	
	if( fabs(old_nb_decay - etmod->pw_nb_decay) > 0.01 )
				etmod->pathways_graph_init_flag = false;

	if( fabs(old_hb_decay - etmod->pw_hb_decay) > 0.01 )
				etmod->pathways_graph_init_flag = false;

	if(fabs( old_cov_decay - etmod->pw_ln_cov_decay) > 0.01)
		etmod->pathways_graph_init_flag = false;

	return wxDialog::TransferDataFromWindow();
}

bool
PathwaysDlgWX::TransferDataToWindow()
{
	wxTextCtrl* cov_decay_edit = (wxTextCtrl*) FindWindow(IDC_PATH_COV_DECAY);
	wxString str;
	str.Printf("%f",exp(etmod->pw_ln_cov_decay));
	str.Trim(false);
	cov_decay_edit->SetValue(str);

	wxTextCtrl* nb_dist_decay = (wxTextCtrl*) FindWindow( IDC_PATH_NB_DECAY );
	str.Printf("%f",etmod->pw_nb_decay);
	str.Trim(false);
	nb_dist_decay->SetValue(str);

	wxTextCtrl* hb_dist_decay = (wxTextCtrl*) FindWindow( IDC_PATH_HB_DECAY );
	str.Printf("%f",etmod->pw_hb_decay);
	str.Trim(false);
    hb_dist_decay->SetValue(str);
	
	wxTextCtrl* sel_thresh_edt = (wxTextCtrl*) FindWindow( IDC_PATH_SEL_THRESH );
	str.Printf("%f",this->sel_thresh);
	str.Trim(false);
	sel_thresh_edt->SetValue(str);

	return wxDialog::TransferDataToWindow();
}


BEGIN_EVENT_TABLE(PathwaysDlgWX,wxDialog)
    EVT_INIT_DIALOG( PathwaysDlgWX::OnInitDialog )
    EVT_BUTTON( IDC_COUPL_MAP,    PathwaysDlgWX::OnCouplMap )
    EVT_BUTTON( IDC_BEST_PATH,    PathwaysDlgWX::OnBestPath )
    EVT_BUTTON( IDC_PATH_SEL_COUPLED,  PathwaysDlgWX::OnSelCoupled )
    EVT_BUTTON( IDC_SET_MODULE_DATA,   PathwaysDlgWX::SetModuleDataFromControls )
    EVT_BUTTON( IDC_DUTTON_CALC,       PathwaysDlgWX::OnDutton )
	EVT_BUTTON( IDC_COLOR_MOL_SURF_BY_COUPL, PathwaysDlgWX::OnColorMolSurfETCoupl )
	EVT_CHECKBOX( IDC_PATH_S_HBOND,    PathwaysDlgWX::OnSHBond )
	EVT_CHECKBOX( IDC_PATH_HBOND,      PathwaysDlgWX::OnHBondTerm )
    EVT_CLOSE( PathwaysDlgWX::OnClose)
END_EVENT_TABLE()

void PathwaysDlgWX::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}


void PathwaysDlgWX::OnBestPath(wxCommandEvent &event)
{
	TransferDataFromWindow();
	wxBusyCursor wait;

	etmod->pathways_calc_type = BEST_PATH;
	etmod->path_coupl_calc();
	
	MolSet* pmset = etmod->GetMolSet();
	HaMolView* pView = pmset->GetActiveMolView();
	if(pView == NULL) return;

	pView->DrawBestPath = True;
	pmset->RefreshAllViews();
}

void PathwaysDlgWX::OnCouplMap(wxCommandEvent &event)
{
	TransferDataFromWindow();
    wxBusyCursor wait;

	etmod->pathways_calc_type = COUPL_MAP;
	etmod->path_coupl_calc();

	MolSet* pmset = etmod->GetMolSet();
	HaMolView* pView = pmset->GetActiveMolView();
	if(pView == NULL) return;
		
	pmset->SelectAtomsAll();
	pView->ScaleColourAttrib(TempAttr);
	pmset->RefreshAllViews();
}

void PathwaysDlgWX::OnSelCoupled(wxCommandEvent &event)
{
	TransferDataFromWindow();
	wxBusyCursor wait;
	etmod->select_important(this->sel_thresh);
}

void PathwaysDlgWX::OnDutton(wxCommandEvent &event)
{
	TransferDataFromWindow();
	wxBusyCursor wait;
    etmod->DuttonModelCalc();	
}

void PathwaysDlgWX::OnColorMolSurfETCoupl( wxCommandEvent &event)
{
	TransferDataFromWindow();
	wxBusyCursor wait;
    etmod->ColorMolSurfETCoupl();	
}

//////////////////////////////////////////////////////////////////////
// ElectrostDlgWX dialog

int ElectrostDlgWX::dlg_open = FALSE;

ElectrostDlgWX::ElectrostDlgWX(ElectrostMod* new_electrost_mod, wxWindow* parent):
wxDialog( parent, -1, "Continuum Electrostatics Calculations", wxDefaultPosition, wxDefaultSize, 
		   wxDEFAULT_DIALOG_STYLE )
{
	electrost_mod= new_electrost_mod;	
    dlg_open = TRUE;
	cont_electr_dlg(this,TRUE);
}

ElectrostDlgWX::~ElectrostDlgWX()
{
    dlg_open = FALSE;
}

void ElectrostDlgWX::OnInitDialog(wxInitDialogEvent& event)
{
	wxChoice* ref_state_cb = (wxChoice*) FindWindow(IDC_DELPHI_REF_STATE);
	ref_state_cb->Append("Vaccuum");
	ref_state_cb->Append("Solvent");
	ref_state_cb->SetSelection(0);

	wxTextCtrl* mol_name_edit= (wxTextCtrl*) FindWindow(IDC_DELPHI_MSET_NAME);
	MolSet* pmset= electrost_mod->GetMolSet();
	if(pmset != NULL && mol_name_edit != NULL)
	{
		mol_name_edit->SetValue(pmset->GetName());
	}

	wxTextCtrl* edt_ctrl;

	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_PERFIL);
	edt_ctrl->SetValidator( wxDoubleValidator(&electrost_mod->perfil, "%5.1f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_OFFSET_X);
	edt_ctrl->SetValidator( wxDoubleValidator(&electrost_mod->offset[0], "%8.3f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_OFFSET_Y);
	edt_ctrl->SetValidator( wxDoubleValidator(&electrost_mod->offset[1], "%8.3f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_OFFSET_Z);
	edt_ctrl->SetValidator( wxDoubleValidator(&electrost_mod->offset[2], "%8.3f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_EPSI);
	edt_ctrl->SetValidator( wxDoubleValidator(&electrost_mod->epsi, "%6.2f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_EPSOUT);
	edt_ctrl->SetValidator( wxDoubleValidator(&electrost_mod->epsout, "%6.2f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_RIONST);
	edt_ctrl->SetValidator( wxDoubleValidator(&electrost_mod->rionst, "%6.2f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_ION_RAD);
	edt_ctrl->SetValidator( wxDoubleValidator(&electrost_mod->exrad, "%6.2f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_PROBE_RAD);
	edt_ctrl->SetValidator( wxDoubleValidator(&electrost_mod->radprb, "%6.2f") );

	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_POT_VAL);
	edt_ctrl->SetValidator( wxFloatValidator(&electrost_mod->pot_isolevel, "%6.2f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_DOTS_NUMB);
	edt_ctrl->SetValidator( wxGenericValidator(&electrost_mod->dots_number) );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_TOT_ENE);
	edt_ctrl->SetValidator( wxDoubleValidator(&electrost_mod->tot_ene,"%12.6f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_NX);
	edt_ctrl->SetValidator( wxGenericValidator(&electrost_mod->nx) );

	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_NLIT);
	edt_ctrl->SetValidator( wxGenericValidator(&electrost_mod->nlit) );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_NNIT);
	edt_ctrl->SetValidator( wxGenericValidator(&electrost_mod->nnit) );

	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_POT_LOW_VAL);
	edt_ctrl->SetValidator( wxDoubleValidator(&electrost_mod->elpot_low_val,"%8.3f") );

	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DELPHI_POT_HIGH_VAL);
	edt_ctrl->SetValidator( wxDoubleValidator(&electrost_mod->elpot_high_val,"%8.3f") );

	event.Skip();
}

bool
ElectrostDlgWX::TransferDataFromWindow()
{	
	bool bres = wxDialog::TransferDataFromWindow();
	if(!bres) return false;
	
	electrost_mod->ny = electrost_mod->nx;
	electrost_mod->nz = electrost_mod->nx;
	
	wxChoice* choice_ctrl = (wxChoice*) FindWindow(IDC_DELPHI_BOUNDARY);
	electrost_mod->boundary = choice_ctrl->GetSelection();

	choice_ctrl = (wxChoice*) FindWindow(IDC_DELPHI_REF_STATE );
	electrost_mod->rdx_shft_mode  = choice_ctrl->GetSelection() + 1;

	wxTextCtrl* text_ctrl = (wxTextCtrl*) FindWindow(IDC_ELFIELD_FILE);
	electrost_mod->elfield_fname = (text_ctrl->GetValue()).c_str();

	wxCheckBox* check_ctrl;
	check_ctrl = (wxCheckBox*) FindWindow(IDC_DELPHI_IPER_X);
	electrost_mod->iper[0] = check_ctrl->GetValue();
	check_ctrl = (wxCheckBox*) FindWindow(IDC_DELPHI_IPER_Y);
	electrost_mod->iper[1] = check_ctrl->GetValue();
	check_ctrl = (wxCheckBox*) FindWindow(IDC_DELPHI_IPER_Z);
	electrost_mod->iper[2] = check_ctrl->GetValue();

	check_ctrl = (wxCheckBox*) FindWindow(IDC_DELPHI_PHIWRT);
	electrost_mod->phiwrt = check_ctrl->GetValue();

	return true; 
}

bool ElectrostDlgWX::TransferDataToWindow()
{	
	bool bres = wxDialog::TransferDataToWindow();
	if(!bres) return false;
		
	wxChoice* choice_ctrl = (wxChoice*) FindWindow(IDC_DELPHI_BOUNDARY);
    choice_ctrl->SetSelection(electrost_mod->boundary);

	choice_ctrl = (wxChoice*) FindWindow(IDC_DELPHI_REF_STATE );
	choice_ctrl->SetSelection(electrost_mod->rdx_shft_mode -1 );

	wxTextCtrl* text_ctrl = (wxTextCtrl*) FindWindow(IDC_ELFIELD_FILE);
	text_ctrl->SetValue(electrost_mod->elfield_fname.c_str());

	wxCheckBox* check_ctrl;
	check_ctrl = (wxCheckBox*) FindWindow(IDC_DELPHI_IPER_X);
	if(electrost_mod->iper[0]) check_ctrl->SetValue(true);
	else check_ctrl->SetValue(false);
	check_ctrl = (wxCheckBox*) FindWindow(IDC_DELPHI_IPER_Y);
	if(electrost_mod->iper[1]) check_ctrl->SetValue(true);
	else check_ctrl->SetValue(false);
	check_ctrl = (wxCheckBox*) FindWindow(IDC_DELPHI_IPER_Z);
	if(electrost_mod->iper[2]) check_ctrl->SetValue(true);
	else check_ctrl->SetValue(false);

	check_ctrl = (wxCheckBox*) FindWindow(IDC_DELPHI_PHIWRT);
	if(electrost_mod->phiwrt) check_ctrl->SetValue(true);
	else check_ctrl->SetValue(false);

	return true;
}


BEGIN_EVENT_TABLE(ElectrostDlgWX, wxDialog)
	EVT_INIT_DIALOG( ElectrostDlgWX::OnInitDialog )
	EVT_BUTTON  (IDC_DELPHI_COLOR_SURF_POT, ElectrostDlgWX::OnColorSurfPot)
	EVT_BUTTON  (IDC_ELFIELD_LOAD, ElectrostDlgWX::OnElFieldLoad)
	EVT_BUTTON  (IDC_ELFIELD_SAVE, ElectrostDlgWX::OnElFieldSave)
	EVT_BUTTON  (IDC_AVG_POT_DON, ElectrostDlgWX::OnAvgPotDon)
	EVT_CLOSE   (ElectrostDlgWX::OnClose)
	EVT_BUTTON  (IDC_DELPHI_SETPAR, ElectrostDlgWX::OnSetParamMod)
	EVT_BUTTON  (IDC_EDIT_ATOM_PARAMS, ElectrostDlgWX::OnEditAtomParams)
	EVT_BUTTON  (IDC_DELPHI_SAVE_INPUT, ElectrostDlgWX::OnSaveInpFiles)
	EVT_BUTTON  (IDC_DELPHI_ET_REORG_ENE, ElectrostDlgWX::OnETReorgEne)
	EVT_BUTTON  (IDC_DELPHI_REDOX_POT, ElectrostDlgWX::OnCalcRedoxPotShft)
	EVT_BUTTON  (IDC_DELPHI_RUN, ElectrostDlgWX::OnRun)
	EVT_BUTTON  (IDC_DELPHI_POT_ISOLVL, ElectrostDlgWX::OnCalcPotIsoLvl)
	EVT_BUTTON  (IDC_IND_CHARGE, ElectrostDlgWX::OnCalcIndCharge)
	EVT_BUTTON  (IDC_PLOT_IND_CHARGE, ElectrostDlgWX::OnPlotIndCharge)
	EVT_BUTTON  (IDC_DELPHI_POT_VDWDOTS, ElectrostDlgWX::OnCalcPotVdwDots) 
	EVT_BUTTON  (IDC_POTPLANEVIEW, ElectrostDlgWX::OnPotPlaneView)
	EVT_BUTTON  (IDC_CONCPLANEVIEW, ElectrostDlgWX::OnConcPlaneView)
END_EVENT_TABLE()


void ElectrostDlgWX::OnClose(wxCloseEvent& event)
{
	ElectrostDlgWX::dlg_open = false;
	event.Skip();
}

void ElectrostDlgWX::OnSetParamMod(wxCommandEvent& event)
{
	TransferDataFromWindow();
}

void
ElectrostDlgWX::OnPotPlaneView(wxCommandEvent& event)
{
	if(electrost_mod == NULL) return;
	if( electrost_mod->el_pot_map.GetNx() <= 0 || electrost_mod->el_pot_map.GetNy() <= 0 || electrost_mod->el_pot_map.GetNz() <= 0)
		return;
	if(electrost_mod->GetMolSet() == NULL) return;
	TransferDataFromWindow();
	
	wxBusyCursor wait;
	
	HaField3D* field=&(electrost_mod->el_pot_map);
	wxString Name=electrost_mod->GetMolSet()->GetName();
	Name<<"El.Pot.";
	//field->SetName(Name.c_str());
			
	PlaneViewOfHaField3D* PlaneV=new PlaneViewOfHaField3D(field,Name.c_str());
	PlaneV->SetMinMax(electrost_mod->elpot_low_val,electrost_mod->elpot_high_val);
	PlaneV->SetHideZeroValues(true);
	
	MolSet* pmset = electrost_mod->GetMolSet();
	pmset->AddObject3D(PlaneV);
	pmset->RefreshAllViews(RFRefresh);
	
	wxFieldPlaneView *PV=new wxFieldPlaneView(PlaneV,pmset,NULL,-1,Name);
	PV->Show(true);
}
void
ElectrostDlgWX::OnConcPlaneView(wxCommandEvent& event)
{
	if(electrost_mod == NULL) return;
	
	//ConcMap0
	if( electrost_mod->ConcMap0.GetNx() <= 0 || electrost_mod->ConcMap0.GetNy() <= 0 || electrost_mod->ConcMap0.GetNz() <= 0)
		return;
	if(electrost_mod->GetMolSet() == NULL) return;
	TransferDataFromWindow();
	
	wxBusyCursor wait;
	
	HaField3D* field=&(electrost_mod->ConcMap0);
	wxString Name=electrost_mod->GetMolSet()->GetName();
	Name<<"Conc.0";
	//field->SetName(Name.c_str());
			
	PlaneViewOfHaField3D* PlaneV=new PlaneViewOfHaField3D(field,Name.c_str());
	PlaneV->SetMinMax(0.0,2.0);
	PlaneV->SetHideZeroValues(true);
	
	MolSet* pmset = electrost_mod->GetMolSet();
	pmset->AddObject3D(PlaneV);
	pmset->RefreshAllViews(RFRefresh);
	
	wxFieldPlaneView *PV=new wxFieldPlaneView(PlaneV,pmset,NULL,-1,Name);
	PV->Show(true);
	
	//ConcMap1
	if( electrost_mod->ConcMap1.GetNx() <= 0 || electrost_mod->ConcMap1.GetNy() <= 0 || electrost_mod->ConcMap1.GetNz() <= 0)
		return;
	if(electrost_mod->GetMolSet() == NULL) return;
	TransferDataFromWindow();
	
	
	HaField3D* field1=&(electrost_mod->ConcMap1);
	wxString Name1=electrost_mod->GetMolSet()->GetName();
	Name1<<"Conc.1";
	//field->SetName(Name.c_str());
			
	PlaneViewOfHaField3D* PlaneV1=new PlaneViewOfHaField3D(field1,Name1.c_str());
	PlaneV1->SetMinMax(0.0,2.0);
	PlaneV1->SetHideZeroValues(true);
	
	pmset->AddObject3D(PlaneV1);
	pmset->RefreshAllViews(RFRefresh);
	
	wxFieldPlaneView *PV1=new wxFieldPlaneView(PlaneV1,pmset,NULL,-1,Name1);
	PV1->Show(true);
}

void ElectrostDlgWX::OnETReorgEne(wxCommandEvent& event)
{
	if(electrost_mod == NULL) return;
	TransferDataFromWindow();
	double reorg_ene = electrost_mod->CalcETReorgEne();
}

void ElectrostDlgWX::OnCalcRedoxPotShft(wxCommandEvent& event)
{
	if(electrost_mod == NULL) return;
	TransferDataFromWindow();
	wxBusyCursor wait;
	electrost_mod->CalcRedoxPotShft();
}


void
ElectrostDlgWX::OnRun(wxCommandEvent& event)
{
	if(electrost_mod == NULL) return;
	wxBusyCursor wait;
	TransferDataFromWindow();
	electrost_mod->run(RUN_BACKGROUND);
	TransferDataToWindow();
}

void ElectrostDlgWX::OnSaveInpFiles(wxCommandEvent& event)
{
	if(electrost_mod == NULL) return;
	TransferDataFromWindow();
	electrost_mod->SaveParamFile();
	electrost_mod->SaveChargeFile();
	electrost_mod->SaveRadiusFile();
	electrost_mod->SaveCoordFile();
}

void ElectrostDlgWX::OnCalcPotIsoLvl(wxCommandEvent& event)
{
	if(electrost_mod == NULL) return;
	TransferDataFromWindow();
	wxBusyCursor wait;
	int result= electrost_mod->BuildPotIsoSurface();
	if(result)
	{
		MolSet* pmset= electrost_mod->GetMolSet();
		pmset->RefreshAllViews(RFRefresh);
	}
}


void ElectrostDlgWX::OnCalcPotVdwDots(wxCommandEvent& event)
{
	if(electrost_mod == NULL) return;
	TransferDataFromWindow();
	bool result;
	wxBusyCursor wait;
	result= electrost_mod->BuildPotVdwDots();
	if(result)
	{
		MolSet* pmset= electrost_mod->GetMolSet();
		pmset->RefreshAllViews(RFRefresh);
	}
}

void ElectrostDlgWX::OnCalcIndCharge(wxCommandEvent& event)
{
    if(electrost_mod == NULL) return;
	TransferDataFromWindow();
	wxBusyCursor wait;
	electrost_mod->CalcIndCharge();
}

void ElectrostDlgWX::OnElFieldLoad(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	wxBusyCursor wait;
	electrost_mod->LoadElPotFromFile();
	TransferDataToWindow();
}

void ElectrostDlgWX::OnElFieldSave(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	wxBusyCursor wait;
	electrost_mod->el_pot_map.SaveToFile(electrost_mod->elfield_fname.c_str());	
}


void ElectrostDlgWX::OnPlotIndCharge(wxCommandEvent& event)
{
	if(electrost_mod == NULL) return;
	TransferDataFromWindow();
	
	wxBusyCursor wait;
	int result= electrost_mod->PlotIndCharge();
	if(result)
	{
		MolSet* pmset= electrost_mod->GetMolSet();
		pmset->RefreshAllViews(RFRefresh);
	}
}
void ElectrostDlgWX::OnColorSurfPot(wxCommandEvent& event) 
{
	if(electrost_mod == NULL) return;
	TransferDataFromWindow();
	bool result;
	wxBusyCursor wait;
	result=electrost_mod->ColorMolSurfElPot();
	if(result)
	{
		MolSet* pmset= electrost_mod->GetMolSet();
		pmset->RefreshAllViews(RFRefresh);
	}
}

void ElectrostDlgWX::OnAvgPotDon(wxCommandEvent& event) 
{
	MolSet* pmset= electrost_mod->GetMolSet();
	AtomGroup* at_grp= pmset->GetAtomGroupByID("DONOR");
	double phi = electrost_mod->CalcAvgPotOn(at_grp);
	PrintLog("Average potential on donor atoms is %10.4f kT \n",phi);
}

void ElectrostDlgWX::OnEditAtomParams(wxCommandEvent& event)
{
	if(electrost_mod == NULL) return;
	MolSet* pmset = electrost_mod->GetMolSet();
	if(pmset == NULL) return;
	
	AtomParamsDlgWX::ResetEditFlags();
	AtomParamsDlgWX::charge_edit_flag = TRUE;
	AtomParamsDlgWX::radius_edit_flag = TRUE;
	AtomParamsDlgWX::temperature_edit_flag = TRUE;

	if(AtomParamsDlgWX::dlg_open) return;
	AtomParamsDlgWX* ptr_atom_params_dlg = new AtomParamsDlgWX( pmset, this );
	ptr_atom_params_dlg->Show(TRUE);	
}


/////////////////////////////////////////////////////////////////////
// QChemParDlgWX Dialog:

int QChemParDlgWX::dlg_open = FALSE;

QChemParDlgWX::QChemParDlgWX(HaQCMod* new_phost_qcmod, wxWindow* parent):
wxFrame( parent, -1, "Quantum Chemical Calculations")
{
	p_qc_mod= new_phost_qcmod;
	dlg_open = TRUE;

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);
 
	wxMenuBar* qchem_1_menu_bar = qchem_menu_1();
    this->SetMenuBar(qchem_1_menu_bar);    

	qchem_param_dlg(this,TRUE);

	OnInitDialog();
}

QChemParDlgWX::~QChemParDlgWX()
{
    dlg_open = FALSE;
}

BEGIN_EVENT_TABLE(QChemParDlgWX, wxFrame)
//    EVT_INIT_DIALOG( QChemParDlgWX::OnInitDialog )
	EVT_MENU(IDC_IPACK_TEST_1, QChemParDlgWX::OnIpackTest1)
	EVT_MENU(IDC_IPACK_TEST_2, QChemParDlgWX::OnIpackTest2)
	EVT_MENU(IDC_TEST_RANDOM, QChemParDlgWX::OnTestRandom)
	EVT_MENU(IDC_CONVERT_CHK_TO_FCHK, QChemParDlgWX::OnConvertChkToFchk)
//	EVT_BUTTON(IDC_RUN_CNDO, QChemParDlgWX::OnRunCndo)
	EVT_BUTTON(IDC_STOP_CALC, QChemParDlgWX::OnStopCalc)
	EVT_MENU(IDC_SET_INT_ENGINE_IPACK, QChemParDlgWX::OnSetIntEngineIPACK)
	EVT_MENU(IDC_SET_INT_ENGINE_GAUSS, QChemParDlgWX::OnSetIntEngGauss)
	EVT_COMBOBOX(IDC_BASIS_SET, QChemParDlgWX::OnSelChangeBasisSet)
	EVT_COMBOBOX(IDC_ACTIVE_ORB, QChemParDlgWX::OnSelChangeActiveOrb)
	EVT_CHOICE(IDC_COMBO_METHOD, QChemParDlgWX::OnChangeMethod)
	EVT_CLOSE( QChemParDlgWX::OnClose )
	EVT_BUTTON   (IDC_QCPAR_SET_PAR, QChemParDlgWX::OnSetPar)
	EVT_BUTTON   (IDC_SAVE_INPFILE,  QChemParDlgWX::OnSaveInpFile)
	EVT_BUTTON   (IDC_QC_RUN_CALC,   QChemParDlgWX::OnRunCalc)
END_EVENT_TABLE()


void QChemParDlgWX::OnInitDialog()
{
	wxComboBox* bas_list= (wxComboBox*) FindWindow(IDC_BASIS_SET);
	bas_list->Append("3-21G");
	bas_list->Append("3-21G_M1");
	bas_list->Append("3-21++G");
	bas_list->Append("6-31G");
	bas_list->Append("6-311G");
	bas_list->Append("HAY1_DZ");
	bas_list->Append("HAY1_DZPP");
	bas_list->Append("MINB6G");
	
	wxChoice* actorb_list = (wxChoice*) FindWindow(IDC_ACTIVE_ORB);
	actorb_list->Append("VALENCE_AO");
	actorb_list->Append("FULL_ATOMIC_BASIS");
	
	wxChoice* method_list= (wxChoice*) FindWindow(IDC_COMBO_METHOD);
	method_list->Append("Hartree-Fock");
	method_list->Append("CNDO/2");
	method_list->Append("INDO/2");
	method_list->Append("ZINDO/1");
	method_list->Append("ZINDO/S");
	method_list->Append("EXTENDED_HUCKEL");
	method_list->Append("MP2");
	method_list->Append("DFT/BLYP3");

	if(p_qc_mod == NULL)
	{
		ErrorInMod(" Error in QChemParDlgWX::OnInitDialog() ",
		           " the pointer to the QChem Module is not valid \n");
		return;
	}

	wxRadioBox* radio_foreground = (wxRadioBox*) FindWindow(IDC_FOREGROUND);
	radio_foreground->SetSelection(0);

	wxRadioBox* qchem_engine = (wxRadioBox*) FindWindow(IDC_QCHEM_ENGINE);
	qchem_engine->SetSelection(1);

	wxTextCtrl* edt_ctrl;

	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_QCPAR_EDIT_CHARGE);
	edt_ctrl->SetValidator( wxGenericValidator(&p_qc_mod->charge) );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_QCPAR_EDIT_MULT);
	edt_ctrl->SetValidator( wxGenericValidator(&p_qc_mod->mult) );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_QCPAR_SOLV_CHRG_EDT);
	edt_ctrl->SetValidator( wxDoubleValidator(&p_qc_mod->distr_ext_chrg, "%6.2f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_MAX_SCF_ITER);
	edt_ctrl->SetValidator( wxGenericValidator(&p_qc_mod->max_scf_iter) );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_NITER_AVG_2);
	edt_ctrl->SetValidator( wxGenericValidator(&p_qc_mod->max_it_avg) );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_FERMI_TEMP);
	edt_ctrl->SetValidator( wxDoubleValidator(&p_qc_mod->temp0_fermi, "%8.2f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_DMAT_CONV);
	edt_ctrl->SetValidator( wxDoubleValidator(&p_qc_mod->conv_dm, "%9.3e") );

	TransferDataToWindow();
}

void QChemParDlgWX::DDX_bas_name(bool from_window )
{
	wxComboBox* bas_list = (wxComboBox*) FindWindow( IDC_BASIS_SET );
	wxString str;
    wxString bas_name = p_qc_mod->GetBasName();
    bas_name.Trim(false);
    bas_name.Trim(true);
	if(from_window)
	{
        str = bas_list->GetValue();
        str.Trim(false);
        str.Trim(true);

		if( p_qc_mod->GetNBfunc() == 0 || !str.Matches(bas_name.c_str()) ) // if basis name changed
		{
			p_qc_mod->InitBasis( str.ToStdString() );
		}
	}
	else
	{
		if(bas_name.IsEmpty())
		{
			bas_list->SetValue("");
		}
		else
		{
			int ipos = bas_list->FindString(bas_name);
			if(ipos >= 0) 
			{
				bas_list->SetSelection(ipos);
			}
			else
			{
				bas_list->Append( bas_name );
				ipos = bas_list->FindString( bas_name );
				bas_list->SetSelection(ipos);
			}
		}
	}
}

void QChemParDlgWX::DDX_active_orb(bool from_window )
{
	wxComboBox* actorb_list = (wxComboBox*) FindWindow( IDC_ACTIVE_ORB );
	if(from_window)
	{
		wxString str = actorb_list->GetValue();
		wxString str_old = p_qc_mod->GetLocOrbSetID();
		if(!str.IsEmpty())
		{
			if( p_qc_mod->GetNActiveOrb() == 0 || str != str_old ) // if basis name changed
			{
				p_qc_mod->InitLocOrb(str.c_str());
			}
		}
	}
	else
	{
		wxString act_orb_name_id = p_qc_mod->GetLocOrbSetID();
		act_orb_name_id.Trim(false);
		act_orb_name_id.Trim(true);
		if(act_orb_name_id.IsEmpty())
		{
			actorb_list->SetValue("");
		}
		else
		{
			int ipos = actorb_list->FindString(act_orb_name_id);
			if(ipos >= 0) 
			{
				actorb_list->SetSelection(ipos);
			}
			else
			{
				actorb_list->Append( act_orb_name_id );
				ipos = actorb_list->FindString(act_orb_name_id);
				actorb_list->SetSelection(ipos);
			}
		}
	}	
}

void QChemParDlgWX::DDX_wave_fun_type(bool from_window )
{
	wxChoice* method_list = (wxChoice*) FindWindow( IDC_COMBO_METHOD );
	wxString str;
	if(from_window)
	{
		str = method_list->GetStringSelection();
		if(!str.IsEmpty())
		{
			p_qc_mod->SetWaveFunType(str.c_str());
		}
	}
	else
	{
		int index;
		if( p_qc_mod->wave_fun_type == harlem::qc::HARTREE_FOCK) index = 0;
		if( p_qc_mod->wave_fun_type == harlem::qc::NDO)
		{
			if(p_qc_mod->ndo_method == harlem::qc::CNDO_2) index = 1;
			if(p_qc_mod->ndo_method == harlem::qc::INDO_2) index = 2;
			if(p_qc_mod->ndo_method == harlem::qc::ZINDO_1) index = 3;
			if(p_qc_mod->ndo_method == harlem::qc::ZINDO_S) index = 4;
		}
		if( p_qc_mod->wave_fun_type == harlem::qc::EXTENDED_HUCKEL) index = 5;
		if( p_qc_mod->wave_fun_type == harlem::qc::MP2) index = 6;
		if( p_qc_mod->wave_fun_type == harlem::qc::DFT_B3LYP) index = 7;

		method_list->SetSelection(index);
	}
}

bool QChemParDlgWX::TransferDataFromWindow()
{	
	MolSet* pmset= p_qc_mod->GetMolSet();
	DDX_bas_name  (true);
	DDX_active_orb(true);
	DDX_wave_fun_type(true);

	wxCheckBox* check_ctrl;
	check_ctrl = (wxCheckBox*) FindWindow(IDC_QCPAR_SOLV_CHRG_FLAG);
	// p_qc_mod->add_ext_chrg_flag = check_ctrl->GetValue();
	check_ctrl = (wxCheckBox*) FindWindow(IDC_UHF);
	p_qc_mod->iuhf = check_ctrl->GetValue();
	check_ctrl = (wxCheckBox*) FindWindow(IDC_GUESS_ONLY);
	p_qc_mod->guess_only = check_ctrl->GetValue();
	check_ctrl = (wxCheckBox*) FindWindow(IDC_GUESS_FROM_MO);
	p_qc_mod->set_guess_from_mos = check_ctrl->GetValue();


	wxTextCtrl* txt_ctrl;

	long ival;
	wxString str;
	HaGaussMod* p_gauss_mod = pmset->GetGaussMod(true);

	check_ctrl = (wxCheckBox*) FindWindow(IDC_QC_SAVE_BASIS_GEN);
	p_gauss_mod->SetSaveBasisSetGen( check_ctrl->GetValue() );

	check_ctrl = (wxCheckBox*) FindWindow(IDC_QC_NO_STD_ORIENT);
	p_gauss_mod->SetNoStdOrient( check_ctrl->GetValue() );

	txt_ctrl = (wxTextCtrl*) FindWindow(IDC_QC_NUM_COMP_NODES);
	str = txt_ctrl->GetValue();
	if( str.ToLong(&ival) ) 
	{
		p_gauss_mod->SetNumProc( (int) ival);
	}
	else
	{
		txt_ctrl->SetValue( wxString::Format("%d",p_gauss_mod->GetNumProc()) );
	}

	txt_ctrl = (wxTextCtrl*) FindWindow(IDC_QC_NUM_SH_MEM_CORES);
	str = txt_ctrl->GetValue();
	if( str.ToLong(&ival) ) 
	{
		p_gauss_mod->SetNumSharedMemCores( (int) ival);
	}
	else
	{
		txt_ctrl->SetValue( wxString::Format("%d",p_gauss_mod->GetNumSharedMemCores()) );
	}

	return wxFrame::TransferDataFromWindow();
}

bool QChemParDlgWX::TransferDataToWindow()
{	
	MolSet* pmset= p_qc_mod->GetMolSet();
	if(pmset != NULL)
	{
		wxTextCtrl* edit_mset_name = (wxTextCtrl*) FindWindow( IDC_MOLSET_NAME );
		wxTextCtrl* edit_inpfile_prefix = (wxTextCtrl*) FindWindow( IDC_INP_FILE_PREFIX );
		edit_mset_name->SetValue(pmset->GetName());
		edit_inpfile_prefix->SetValue(pmset->GetName());
	}
	DDX_bas_name  (false);
	DDX_active_orb(false);
	DDX_wave_fun_type(false);

	wxCheckBox* check_ctrl;
	check_ctrl = (wxCheckBox*) FindWindow(IDC_QCPAR_SOLV_CHRG_FLAG);
	if(!p_qc_mod->ext_chrg.empty()) check_ctrl->SetValue(true);
	else check_ctrl->SetValue(false);
	check_ctrl = (wxCheckBox*) FindWindow(IDC_UHF);
	if(p_qc_mod->iuhf) check_ctrl->SetValue(true);
	else check_ctrl->SetValue(false);
	check_ctrl = (wxCheckBox*) FindWindow(IDC_GUESS_ONLY);
	if(p_qc_mod->guess_only) check_ctrl->SetValue(true);
	else check_ctrl->SetValue(false);
	check_ctrl = (wxCheckBox*) FindWindow(IDC_GUESS_FROM_MO);
	if(p_qc_mod->set_guess_from_mos) check_ctrl->SetValue(true);
	else check_ctrl->SetValue(false);
	
	wxTextCtrl* txt_ctrl;

	long ival;
	wxString str;
	HaGaussMod* p_gauss_mod = pmset->GetGaussMod(true);

	check_ctrl = (wxCheckBox*) FindWindow(IDC_QC_SAVE_BASIS_GEN);
	check_ctrl->SetValue( p_gauss_mod->save_basis_set_gen );

	check_ctrl = (wxCheckBox*) FindWindow(IDC_QC_NO_STD_ORIENT);
	check_ctrl->SetValue( p_gauss_mod->no_std_orient );
	
	txt_ctrl = (wxTextCtrl*) FindWindow(IDC_QC_NUM_COMP_NODES);
	txt_ctrl->SetValue( wxString::Format("%d",p_gauss_mod->GetNumProc()) );
	txt_ctrl = (wxTextCtrl*) FindWindow(IDC_QC_NUM_SH_MEM_CORES);
	txt_ctrl->SetValue( wxString::Format("%d",p_gauss_mod->GetNumSharedMemCores()) );

	return wxFrame::TransferDataToWindow();
}

void QChemParDlgWX::OnSelChangeBasisSet(wxCommandEvent& event) 
{
	PrintLog(" QChemParDlgWX::OnSelChangeBasisSet() \n");
	TransferDataFromWindow();
}

void QChemParDlgWX::OnSelChangeActiveOrb(wxCommandEvent& event) 
{
	PrintLog(" QChemParDlgWX::OnSelChangeActiveOrb()  \n");
	TransferDataFromWindow();
}


void QChemParDlgWX::OnClose(wxCloseEvent& event)
{
    dlg_open = FALSE;
	event.Skip();
}

void QChemParDlgWX::OnSetPar(wxCommandEvent& event)
{
	TransferDataFromWindow();
}

void QChemParDlgWX::OnSaveInpFile(wxCommandEvent& event)
{
    TransferDataFromWindow();

    wxRadioBox* check_meth = (wxRadioBox*) FindWindow(IDC_QCHEM_ENGINE);
	int isel = check_meth->GetSelection();
	
	wxTextCtrl* edit_inpfile_prefix = (wxTextCtrl*) FindWindow( IDC_INP_FILE_PREFIX );
	wxString str = edit_inpfile_prefix->GetValue();

    MolSet* pmset = p_qc_mod->GetMolSet();
	
	if( isel == 1 )
	{
		HaGaussMod* gauss_mod = pmset->GetGaussMod(true);
        gauss_mod->SetFilePrefix(str.c_str());
		gauss_mod->SaveInpFile();
	}
	else if( isel == 2  )
	{
		HaDaltonMod* dalton_mod = pmset->GetDaltonMod(true);
	    dalton_mod->SetInpFilePrefix(str.c_str());
		dalton_mod->Save_mol_file();
		dalton_mod->Save_param_file();	
	}
	else if( isel == 3  )
	{
		HaZindoMod* zindo_mod = pmset->GetZindoMod(true);
	    zindo_mod->SaveParamFile("zindo.inp");
	}
}

void QChemParDlgWX::OnRunCalc(wxCommandEvent& event)
{
	TransferDataFromWindow();

	OnSaveInpFile(event);
    wxRadioBox* check_meth = (wxRadioBox*) FindWindow(IDC_QCHEM_ENGINE);
    wxRadioBox* check_foregr = (wxRadioBox*) FindWindow(IDC_FOREGROUND);

	wxTextCtrl* edit_inpfile_prefix = (wxTextCtrl*) FindWindow( IDC_INP_FILE_PREFIX );
	wxString str = edit_inpfile_prefix->GetValue();

    MolSet* pmset = p_qc_mod->GetMolSet();

	wxChoice* method_list= (wxChoice*) FindWindow(IDC_COMBO_METHOD);

	if(p_qc_mod->wave_fun_type == harlem::qc::NDO)
	{
		TransferDataFromWindow();
		p_qc_mod->RunCNDO();
		TransferDataToWindow();
		return;
	}

    int imeth = check_meth->GetSelection();
    int ifg   = check_foregr->GetSelection();

	if(p_qc_mod->wave_fun_type == harlem::qc::EXTENDED_HUCKEL)
	{
		wxBusyCursor wait;
		p_qc_mod->RunExtHuckel();
		TransferDataToWindow();
		return;
	}

	if( imeth == 1)
	{
		HaGaussMod* gauss_mod = pmset->GetGaussMod(true);
		gauss_mod->SaveInpFile();
		harlem::RunOptions ropt;
		if( ifg == 0 ) ropt.SetRunSync( false );
		else ropt.SetRunSync( true );
		
        gauss_mod->Run( &ropt );
		
		
	}
	else if( imeth == 2 )
	{
        HaDaltonMod* dalton_mod = pmset->GetDaltonMod(true);
		dalton_mod->Save_mol_file();
		dalton_mod->Save_param_file();
		if( ifg == 0)
		{
			dalton_mod->run(RUN_FOREGROUND);
		}
		else if( ifg == 1)
		{
            dalton_mod->run(RUN_BACKGROUND);
		}
	}
	else if( imeth == 3 )
	{
        HaZindoMod* zindo_mod = pmset->GetZindoMod(true);
		zindo_mod->SaveParamFile("zindo.inp");
	}
	TransferDataToWindow();
}

void QChemParDlgWX::OnConvertChkToFchk(wxCommandEvent& event)
{
	MolSet* pmset = p_qc_mod->GetMolSet();
	HaGaussMod* gauss_mod = pmset->GetGaussMod(true);
	gauss_mod->RunFormChk(gauss_mod->GetCHKFileName().c_str(),gauss_mod->GetFCHKFileName().c_str());
}

void QChemParDlgWX::OnStopCalc(wxCommandEvent& event) 
{
	p_qc_mod->stop_calc_flag = 1;
}


void QChemParDlgWX::OnRunCndo(wxCommandEvent& event) 
{
	wxBusyCursor wait;
	p_qc_mod->RunCNDO();
	TransferDataToWindow();
}


void QChemParDlgWX::OnIpackTest1(wxCommandEvent& event) 
{
#ifdef USE_IPACK
	p_qc_mod->TestIPack1();
#endif
}

void QChemParDlgWX::OnIpackTest2(wxCommandEvent& event) 
{
#ifdef USE_IPACK
	p_qc_mod->TestIPack2();
#endif
}


void QChemParDlgWX::OnTestRandom(wxCommandEvent& event) 
{
	p_qc_mod->TestRandomGen();
}

void QChemParDlgWX::OnChangeMethod(wxCommandEvent& event) 
{
	wxChoice* method_list= (wxChoice*) FindWindow(IDC_COMBO_METHOD);
	wxString str = method_list->GetStringSelection();
	p_qc_mod->SetWaveFunType(str.c_str());

	wxComboBox* bas_list= (wxComboBox*) FindWindow(IDC_BASIS_SET);
	wxComboBox* actorb_list = (wxComboBox*) FindWindow(IDC_ACTIVE_ORB);

	if(p_qc_mod->wave_fun_type == harlem::qc::NDO || p_qc_mod->wave_fun_type == harlem::qc::EXTENDED_HUCKEL )
	{
		bas_list->SetSelection(7);
        actorb_list->SetSelection(1);
	}
	TransferDataFromWindow();
}

void QChemParDlgWX::OnSetIntEngineIPACK(wxCommandEvent& event) 
{
	HaQCMod::SetQCIntEngine(QCIntEngineType::INT_ENGINE_IPACK);
}

void QChemParDlgWX::OnSetIntEngGauss(wxCommandEvent& event) 
{
	HaQCMod::SetQCIntEngine(QCIntEngineType::INT_ENGINE_GAUSS);
}

/////////////////////////////////////////////////////////////////////////////
// ETEffHamDlgWX dialog

int ETEffHamDlgWX::dlg_open = FALSE;

ETEffHamDlgWX::ETEffHamDlgWX(ETCouplMod* new_ptr_et_mod, wxWindow* parent):
wxFrame( parent, -1, "Donor/Acceptor Electronic Coupling Calculations")
{
	ptr_et_mod= new_ptr_et_mod;
    dlg_open = TRUE;

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);
 
	wxMenuBar* et_1_menu_bar = et_menu_1();
    this->SetMenuBar(et_1_menu_bar);    

	da_coupl_calc_dlg(this,TRUE);

	OnInitDialog();
}

ETEffHamDlgWX::~ETEffHamDlgWX()
{
    dlg_open = FALSE;
}


void
ETEffHamDlgWX::OnInitDialog()
{
	  wxTextCtrl* db_file  = (wxTextCtrl*) FindWindow(IDC_ETEFFH_DB_FILE);
	  wxTextCtrl* da_field = (wxTextCtrl*) FindWindow(IDC_ETEFFH_DA_FIELD);
	  wxTextCtrl* step = (wxTextCtrl*) FindWindow(IDC_ETEFFH_STEP);
	  wxTextCtrl* fin_field = (wxTextCtrl*) FindWindow(IDC_ETEFFH_FIN_FIELD);
      wxTextCtrl* firsteig = (wxTextCtrl*) FindWindow(IDC_ETEFFH_FIRSTEIG);
	  wxTextCtrl* lasteig = (wxTextCtrl*) FindWindow(IDC_ETEFFH_LASTEIG);
	  db_file->SetValue(ptr_et_mod->db_file_name.c_str());
	  da_field->SetValue("0.0");
	  fin_field->SetValue("0.01");
	  step->SetValue("0.001");
      firsteig->SetValue("1");
      lasteig->SetValue("1");
	  wxTextCtrl* edit_eigvec_isolvl = (wxTextCtrl*) FindWindow( IDC_ETEFFH_ISOLVL);
	  edit_eigvec_isolvl->SetValue("0.05");

	  HaQCMod* ptr_qcmod = ptr_et_mod->GetQCMod();
	  wxTextCtrl* grid_size = (wxTextCtrl*) FindWindow(IDC_ETEFFH_GRID_SIZE);
	  if(ptr_qcmod) grid_size->SetValidator(wxGenericValidator(&ptr_qcmod->m_grid_size));

	  TransferDataToWindow();
}

void ETEffHamDlgWX::DDX_eig_val(int set_var)
{
	int i;
	wxString str;
	
	wxListBox* eigval_list= (wxListBox*) FindWindow(IDC_ETEFFH_EIGVAL_LIST);
	
	if(set_var)
	{
		
	}
	else
	{
		eigval_list->Clear();
		HaVec_double& enel= ptr_et_mod->enel;
		for(i=1; i <= enel.size(); i++)
		{
			str.Printf("%4d :%16.9f",i,enel(i));
			eigval_list->Append(str);
		}
	}
}


bool
ETEffHamDlgWX::TransferDataToWindow()
{
	wxString str;

	DDX_tun_ener(false);
	DDX_da_field(false);
	DDX_eig_val(false);
	DDX_list_loc_orb_gf(false);

	FillListBoxHeffEigene();
	wxListBox* da_val_list= (wxListBox*) FindWindow(IDC_DA_COUPL);
	da_val_list->Clear();
	int id,ia;
	int ndon = ptr_et_mod->donor_orbs.GetNOrbs();
	int nacc = ptr_et_mod->acc_orbs.GetNOrbs();
	if( ptr_et_mod->da_coupl_val.num_rows() >= ndon &&  ptr_et_mod->da_coupl_val.num_cols() >= nacc)
	{
		for(id = 0;  id < ndon; id++)
		{
			for(ia = 0; ia < nacc; ia++)
			{
				str.Printf("%3d %3d :  %12.6e ",id,ia,ptr_et_mod->da_coupl_val.GetVal_idx0(id,ia));
				da_val_list->Append(str);
			}
		}
	}
	HaQCMod* ptr_qcmod = ptr_et_mod->GetQCMod();
	if(ptr_qcmod != NULL)
	{
		int nae = ptr_qcmod->GetNumAlphaEl(1);
		wxTextCtrl* et_homo_edt = (wxTextCtrl*) FindWindow(IDC_ET_HOMO);
		wxString str;
		str.Printf("%d",nae);
		et_homo_edt->SetValue(str);
	}
	wxCheckBox* check_ctrl = (wxCheckBox*) FindWindow(IDC_ETEFFH_DAB_INT_HUCK);
	if(ptr_et_mod->set_dab_huck_inter) check_ctrl->SetValue(true);
	else check_ctrl->SetValue(false);

	return wxFrame::TransferDataToWindow();
}

bool ETEffHamDlgWX::TransferDataFromWindow()
{
	DDX_tun_ener(true);
	DDX_da_field(true);
	DDX_eig_val(true);
	DDX_list_loc_orb_gf(true);
	
	wxCheckBox* check_ctrl;
	check_ctrl = (wxCheckBox*) FindWindow(IDC_ETEFFH_DAB_INT_HUCK);
	ptr_et_mod->set_dab_huck_inter = check_ctrl->GetValue();
	return wxFrame::TransferDataFromWindow();
}


BEGIN_EVENT_TABLE(ETEffHamDlgWX, wxFrame)
	EVT_BUTTON(IDC_ETEFFH_SET_PAR, ETEffHamDlgWX::OnSetParam)
	EVT_BUTTON(IDC_ETEFFH_CLEAR_HEFF, ETEffHamDlgWX::OnClearHeff)
	EVT_BUTTON(IDC_ETEFFH_CALC_HEFF, ETEffHamDlgWX::OnCalcHeff)
    EVT_BUTTON(IDC_ETEFFH_SCAN_ENE, ETEffHamDlgWX::OnCalcScanEne)	
	EVT_BUTTON(IDC_ETEFFH_DIAG, ETEffHamDlgWX::OnDiagHeff)
	EVT_BUTTON(IDC_MIN_SPLIT, ETEffHamDlgWX::OnMinSplit)
	EVT_MENU(IDC_ETEFFH_SAVE_DB, ETEffHamDlgWX::OnSaveDB)
	EVT_MENU(IDC_ETEFFH_BUILD_DB, ETEffHamDlgWX::OnBuildDB)
	EVT_MENU(IDC_SAVE_HEFF_XML, ETEffHamDlgWX::OnSaveHeffXml)
	EVT_MENU(IDC_LOAD_HEFF_XML, ETEffHamDlgWX::OnLoadHeffXml)
	EVT_MENU(IDC_SAVE_REDOX_ORBS_XML, ETEffHamDlgWX::OnSaveRedoxOrbsXml)
	EVT_MENU(IDC_LOAD_DONOR_ORBS_XML, ETEffHamDlgWX::OnLoadDonorOrbsXml)
	EVT_MENU(IDC_LOAD_ACC_ORBS_XML, ETEffHamDlgWX::OnLoadAccOrbsXml)
	EVT_BUTTON(IDC_CHOOSE_DONOR_ORB, ETEffHamDlgWX::OnChooseRedoxOrbs)
	EVT_BUTTON(IDC_CHOOSE_ACC_ORB, ETEffHamDlgWX::OnChooseRedoxOrbs)
	EVT_BUTTON(IDC_ETEFFH_RESET_SRC_ORB, ETEffHamDlgWX::OnResetSrcLorb)
	EVT_BUTTON(IDC_ETEFFH_RESET_TGT_ORB, ETEffHamDlgWX::OnResetTgtLorb)
	EVT_BUTTON(IDC_ETEFFH_CALC_GF_MO,    ETEffHamDlgWX::OnCalcGFfromMO)
	EVT_BUTTON(IDC_ETEFFH_CALC_GF_HEFF,  ETEffHamDlgWX::OnCalcGFfromHeff)
	EVT_MENU(IDC_ETEFFH_PROTECT_MAT,  ETEffHamDlgWX::OnPrintProtectMat)
	EVT_BUTTON(IDC_ETEFFH_PLOT_EIGVEC,  ETEffHamDlgWX::OnPlotEigenVec)
	EVT_BUTTON(IDC_ETEFFH_HDA_FROM_GF,  ETEffHamDlgWX::OnCalcHDAfromGF)
	EVT_BUTTON(IDC_ETEFFH_HDA_PERT,     ETEffHamDlgWX::OnCalcHDAPert)
	EVT_BUTTON(IDC_ETEFFH_PRINT_EIGVEC, ETEffHamDlgWX::OnPrintEigVec)
	EVT_BUTTON(IDC_ETEFFH_ZERO_LONG, ETEffHamDlgWX::OnZeroLong)
	EVT_MENU(IDC_ETEFFH_PR_OVLP_EL, ETEffHamDlgWX::OnPrintOvlpElem)
	EVT_MENU(IDC_ETEFFH_PR_HEFF_EL, ETEffHamDlgWX::OnPrintHeffElem)
	EVT_BUTTON(IDC_HEFF_EDIT_DIAG, ETEffHamDlgWX::OnHeffEditDiag)
	EVT_BUTTON(IDC_ET_SET_DONOR_TRUNC_MO, ETEffHamDlgWX::OnSetDonorTruncMo)
	EVT_CLOSE(ETEffHamDlgWX::OnClose)
	EVT_BUTTON(IDC_ET_CALC_HEFF_HUCK, ETEffHamDlgWX::OnCalcHeffHuck)
//	EVT_RIGHT_DOWN( ETEffHamDlgWX::OnRightDown )
END_EVENT_TABLE()

//void
//ETEffHamDlgWX::OnRightDown(wxMouseEvent& event)
//{
//	PrintLog("Mouse Event ID: %d \n", event.GetId());
//	PrintLog("IDC_DONOR_ORB = %d \n", IDC_DONOR_ORB);
//	PrintLog("IDC_ACC_ORB = %d \n", IDC_ACC_ORB);
//}


void ETEffHamDlgWX::OnSetParam(wxCommandEvent& event)
{
	if(ptr_et_mod == NULL) return;
	wxBusyCursor wait;
	TransferDataFromWindow();
}

void ETEffHamDlgWX::OnClearHeff(wxCommandEvent& event)
{
	if(ptr_et_mod == NULL) return;
	ptr_et_mod->heff_mat.newsize(0,0);
	wxListBox* eigval_list= (wxListBox*) FindWindow(IDC_ETEFFH_EIGVAL_LIST);
	eigval_list->Clear();
}

void ETEffHamDlgWX::OnCalcHeff(wxCommandEvent& event)
{
	if(ptr_et_mod == NULL) return;
	wxBusyCursor wait;
	TransferDataFromWindow();
	ptr_et_mod->RecalcHeff();
	
}

void ETEffHamDlgWX::OnCalcScanEne(wxCommandEvent& event)
{
	wxTextCtrl* step;
	step = (wxTextCtrl*) FindWindow(IDC_ETEFFH_STEP);
    wxTextCtrl* da_field = (wxTextCtrl*) FindWindow(IDC_ETEFFH_DA_FIELD); 
    wxTextCtrl* fin_field = (wxTextCtrl*) FindWindow(IDC_ETEFFH_FIN_FIELD);
    wxTextCtrl* firsteig = (wxTextCtrl*) FindWindow(IDC_ETEFFH_FIRSTEIG);
	wxTextCtrl* lasteig = (wxTextCtrl*) FindWindow(IDC_ETEFFH_LASTEIG);

	double step_val;
	double ifield_val;
	double ffield_val;
	long first_eig_val;
	long last_eig_val;

	wxString str;

	str = step->GetValue();
	str.ToDouble(&step_val);
    str = da_field->GetValue();
    str.ToDouble(&ifield_val); 
	str = fin_field->GetValue();
	str.ToDouble(&ffield_val); 
	str = firsteig->GetValue();
    str.ToLong(&first_eig_val);
    str = lasteig->GetValue();
    str.ToLong(&last_eig_val);

	wxBusyCursor wait;
    ptr_et_mod->ScanEigEneField(first_eig_val, last_eig_val, ifield_val, ffield_val, step_val);
}

void ETEffHamDlgWX::OnDiagHeff(wxCommandEvent& event)
{
	if(ptr_et_mod == NULL) return;
	wxBusyCursor wait;
	TransferDataFromWindow();
	ptr_et_mod->DiagHeff();
	TransferDataToWindow();
	
}

void ETEffHamDlgWX::OnCalcHDAfromGF(wxCommandEvent& event)
{
	if(ptr_et_mod == NULL) return;
	wxBusyCursor wait;
	TransferDataFromWindow();
	ptr_et_mod->CalcHDAfromGF();
	TransferDataToWindow();
	
}

void ETEffHamDlgWX::OnCalcHDAPert(wxCommandEvent& event)
{
	if(ptr_et_mod == NULL) return;
	wxBusyCursor wait;
	TransferDataFromWindow();
	double et_coupl;
	ptr_et_mod->CalcHDAPert(et_coupl);
	TransferDataToWindow();
	
}

void ETEffHamDlgWX::FillListBoxHeffEigene()
{
	int i;
	wxString str;
	
	wxListBox* eigval_list= (wxListBox*) FindWindow(IDC_ETEFFH_EIGVAL_LIST);
	if(eigval_list== NULL) return;
	eigval_list->Clear();
	if(ptr_et_mod == NULL) return;
	HaVec_double& enel= ptr_et_mod->enel;
	
	for(i=1; i <= enel.size(); i++)
	{
		str.Printf("%4d : %16.9f",i,enel(i));
		eigval_list->Append(str);
	}
}

void ETEffHamDlgWX::OnMinSplit(wxCommandEvent& event)
{
//	double fld,splt;
//  wxString str;
//	wxTextCtrl* da_field= (wxTextCtrl*) FindWindow(IDC_ETEFFH_DA_FIELD);
//	str = da_field->GetValue();
//	bool bres = str.ToDouble(&fld);
	wxBusyCursor wait;
	ptr_et_mod->CalcHDAEneSplit();
	TransferDataToWindow();
}

void ETEffHamDlgWX::OnSaveDB(wxCommandEvent& event)
{
	wxString str;
	if(ptr_et_mod == NULL) return;
	wxTextCtrl* db_file = (wxTextCtrl*) FindWindow(IDC_ETEFFH_DB_FILE);
	str = db_file->GetValue();
	ptr_et_mod->db_file_name = str.c_str();
	wxBusyCursor wait;
	ptr_et_mod->PutSubMatToDB();
}

void ETEffHamDlgWX::OnBuildDB(wxCommandEvent& event)
{
	wxString str;
	if(ptr_et_mod == NULL) return;
	wxTextCtrl* db_file = (wxTextCtrl*) FindWindow(IDC_ETEFFH_DB_FILE);
	str = db_file->GetValue();	
	ptr_et_mod->db_file_name = str.c_str();
	wxBusyCursor wait;
	ptr_et_mod->GetSubMatFromDB();
}

void ETEffHamDlgWX::OnSaveHeffXml(wxCommandEvent& event)
{
	wxString str;
	if(ptr_et_mod == NULL) return;

	wxString path = wxFileSelector("Select Xml File to Save Heff");
	if ( !path )
        return;

	FILE* fout = fopen(path.c_str(),"w");
	if(fout != NULL)
	{
		ptr_et_mod->SaveHeffXml(fout);
		fclose(fout);
	}
}

void ETEffHamDlgWX::OnLoadHeffXml(wxCommandEvent& event)
{
	wxString str;
	if(ptr_et_mod == NULL) return;

	wxString path = wxFileSelector("Select Xml File to Load Submatrix of Heff");
	
	if ( !path )
        return;

	FILE* finp = fopen(path.c_str(),"r");
	if(finp != NULL)
	{
		ptr_et_mod->LoadFragmHeffXml(finp);
		fclose(finp);
	}
}

void ETEffHamDlgWX::OnSaveRedoxOrbsXml(wxCommandEvent& event)
{
	wxString str;
	if(ptr_et_mod == NULL) return;

	wxString path = wxFileSelector("Select Xml File to Save Redox Orbitals");
	if ( !path )
        return;

	int ndon = ptr_et_mod->donor_orbs.GetNOrbs();
	int nacc = ptr_et_mod->acc_orbs.GetNOrbs();

	if( ndon == 0 && nacc == 0)
	{
		PrintLog("\nError in ETEffHamDlgWX::OnSaveRedoxOrbsXml() \n");
		PrintLog("No Redox Orbitals are defined \n\n");
		return;
	}

	FILE* fout = fopen(path.c_str(),"w");
	if(fout != NULL)
	{
		if(ndon > 0)
		{
			ptr_et_mod->donor_orbs.SaveXML(fout);
		}
		else
		{
			ptr_et_mod->acc_orbs.SaveXML(fout);
		}
		fclose(fout);
	}
}

void ETEffHamDlgWX::OnLoadDonorOrbsXml(wxCommandEvent& event)
{
	wxString str;
	if(ptr_et_mod == NULL) return;

	wxString path = wxFileSelector("Select Xml File to Load Donor Orbs");
	
	if ( !path )
        return;

	FILE* finp = fopen(path.c_str(),"r");
	if(finp != NULL)
	{
		ptr_et_mod->donor_orbs.LoadXmlFile(finp);
		fclose(finp);
	}
}

void ETEffHamDlgWX::OnLoadAccOrbsXml(wxCommandEvent& event)
{
	wxString str;
	if(ptr_et_mod == NULL) return;

	wxString path = wxFileSelector("Select Xml File to Load Acceptor Orbs");
	
	if ( !path )
        return;

	FILE* finp = fopen(path.c_str(),"r");
	if(finp != NULL)
	{
		ptr_et_mod->acc_orbs.LoadXmlFile(finp);
		fclose(finp);
	}
}


void ETEffHamDlgWX::OnPlotEigenVec(wxCommandEvent& event)
{
	wxString str;
	TransferDataFromWindow();
	wxListBox* eigval_list= (wxListBox*) FindWindow(IDC_ETEFFH_EIGVAL_LIST);

	wxArrayInt sel_arr;

	int neig = eigval_list->GetCount();
	int idx;
	int nsel = eigval_list->GetSelections(sel_arr);
	if(nsel == 0 || neig == 0)
	{
		PrintLog(" Error in ETEffHamDlgWX::OnPlotEigenVec() \n"); 
		PrintLog(" Empty eigen vector list or No selection \n"); 
		return;
	}

	wxTextCtrl* edit_isolvl = (wxTextCtrl*) FindWindow( IDC_ETEFFH_ISOLVL);
	str = edit_isolvl->GetValue();
		
	double flvl;
	bool bres;

	bres = str.ToDouble(&flvl);

	if(!bres)
	{
		std::cerr << " Error in ETEffHamDlgWX::OnPlotEigenVec() " << std::endl;
		std::cerr << " Error Reading EigenVector isolevel value " << std::endl;
		return;
	}

	HaQCMod* ptr_qcmod = ptr_et_mod->GetQCMod();
	if( ptr_qcmod == NULL)
		return;

	wxBusyCursor wait;
	for( idx = 0; idx < nsel; idx++)
	{	int isel = sel_arr[idx];
		ptr_et_mod->CreateEigVecContour(isel+1, flvl, ptr_qcmod->m_grid_size);
	}
	(ptr_et_mod->GetMolSet())->RefreshAllViews(RFRefresh);
	
}

void ETEffHamDlgWX::OnPrintEigVec(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	
	wxListBox* eigval_list= (wxListBox*) FindWindow(IDC_ETEFFH_EIGVAL_LIST);
	
	wxArrayInt sel_arr;

	int neig = eigval_list->GetCount();
	int idx;
	int nsel = eigval_list->GetSelections(sel_arr);
	if(nsel == 0 || neig == 0)
	{
		PrintLog(" Error in ETEffHamDlgWX::OnPrintEigenVec() \n"); 
		PrintLog(" Empty eigen vector list or No selection \n"); 
		return;
	}

	HaQCMod* ptr_qcmod = ptr_et_mod->GetQCMod();
	if( ptr_qcmod == NULL)
		return;

	wxBusyCursor wait;
	for( idx = 0; idx < nsel; idx++)
	{	int isel = sel_arr[idx];
		bool result= ptr_et_mod->PrintEigVecCoef(isel);
	}
}

void ETEffHamDlgWX::OnChooseRedoxOrbs(wxCommandEvent& event)
{
	wxArrayString actions;
	
	int choose_donor = FALSE;

	if( event.GetId() == IDC_CHOOSE_DONOR_ORB )
	{
		actions.Add("Pick Donor Orbitals From The Basis");
		actions.Add("Retrieve Donor Orbitals From a Fragment");
		choose_donor = TRUE;
	}
	if( event.GetId() == IDC_CHOOSE_ACC_ORB )
	{
		actions.Add("Pick Acceptor Orbitals From The Basis");
		actions.Add("Retrieve Acceptor Orbitals From a Fragment");
		choose_donor = FALSE;
	}

	wxString sres = wxGetSingleChoice("Get Redox Orbitals ", "Available Methods: ", actions);

	if( sres == "Pick Donor Orbitals From The Basis")
	{
		OnPickDonorOrbBasis(event);
	}
	else if( sres == "Pick Acceptor Orbitals From The Basis" )
	{
		OnPickAccOrbBasis(event);
	}
	else if( sres == "Retrieve Donor Orbitals From a Fragment" )
	{
		MolSet* pmset = this->ptr_et_mod->GetMolSet();
		int nf = pmset->Fragments.size();
		int i;
		wxArrayString frag_names;
		frag_names.Alloc(nf);
		for(i = 0; i < nf; i++)
		{
			MolSet* pfrag = (MolSet*) pmset->Fragments[i];
			frag_names.Add(pfrag->GetName());
		}

		int idx_frag = wxGetSingleChoiceIndex("Choose Fragment to Retrieve Redox Orbitals", "List of Fragments: ", frag_names);
		
		if( idx_frag > -1)
		{
			MolSet* pfrag = (MolSet*) pmset->Fragments[idx_frag];
			ETCouplMod* etmod_frag = pfrag->GetETCouplMod(0);
			int ires = ptr_et_mod->GetRedoxOrbsFromFrag(pfrag);
		}
	}
}


void ETEffHamDlgWX::OnPickDonorOrbBasis(wxCommandEvent& event)
{
	if(ptr_et_mod == NULL) return;

	SelectLocOrbDlgWX* ptr_sel_lorb_dlg= new SelectLocOrbDlgWX(ptr_et_mod,0, this);
	ptr_sel_lorb_dlg->phost_dialog = this;
	ptr_sel_lorb_dlg->Show(TRUE);
}

void ETEffHamDlgWX::OnPickAccOrbBasis(wxCommandEvent& event)
{
	if(ptr_et_mod == NULL) return;

	SelectLocOrbDlgWX* ptr_sel_lorb_dlg= new SelectLocOrbDlgWX(ptr_et_mod,1, this);
	ptr_sel_lorb_dlg->phost_dialog = this;
	ptr_sel_lorb_dlg->Show(TRUE);

}

void ETEffHamDlgWX::OnResetSrcLorb(wxCommandEvent& event)
{
	if(ptr_et_mod == NULL) return;

	wxListBox* list_lorb_src = (wxListBox*) FindWindow(IDC_DONOR_ORB);
	list_lorb_src->Clear();
	ptr_et_mod->donor_orbs.Clear();
}

void ETEffHamDlgWX::OnResetTgtLorb(wxCommandEvent& event)
{
	if(ptr_et_mod == NULL) return;

	wxListBox* list_lorb_tgt = (wxListBox*) FindWindow(IDC_ACC_ORB);
	list_lorb_tgt->Clear();
	ptr_et_mod->acc_orbs.Clear();
}



void ETEffHamDlgWX::DDX_list_loc_orb_gf(int set_var)
{
	wxListBox* list_lorb_src = (wxListBox*) FindWindow(IDC_DONOR_ORB);
	wxListBox* list_lorb_tgt = (wxListBox*) FindWindow(IDC_ACC_ORB);

	if(ptr_et_mod == NULL)
		return;

	HaQCMod* ptr_qc_mod = ptr_et_mod->GetQCMod();
		
	if(ptr_qc_mod == NULL)
		return;

	if(set_var) // Set lists of local orbitals pointers in ET module
	{

	}
	else // Set lists of local orbitals references from list of source and target orbitals of the ET module
	{
		list_lorb_src->Clear();
		list_lorb_tgt->Clear();

		int nsrc = ptr_et_mod->donor_orbs.GetNOrbs();
		int ntgt = ptr_et_mod->acc_orbs.GetNOrbs();
		int i;
		for(i = 0; i < nsrc; i++)
		{
			list_lorb_src->Append((ptr_et_mod->donor_orbs.GetLabel(i)).c_str());
		}
		for(i = 0; i < ntgt; i++)
		{
			list_lorb_tgt->Append((ptr_et_mod->acc_orbs.GetLabel(i)).c_str());
		}
	}
}

void ETEffHamDlgWX::DDX_tun_ener(int set_var)
{
	wxString str;
	wxTextCtrl* edit_tun_ener = (wxTextCtrl*) FindWindow(IDC_ETEFFH_TUN_ENER);
	if(set_var)
	{
		str = edit_tun_ener->GetValue();
		double fval;
		bool bres = str.ToDouble(&fval);
		if(bres) ptr_et_mod->SetTunEne(fval);
	}
	else
	{
		str.Printf("%12.6f",ptr_et_mod->GetTunEne());
		str.Trim(false);
		edit_tun_ener->SetValue(str);
	}
}

void ETEffHamDlgWX::DDX_da_field(int set_var)
{
	double tmp;
	wxString str;

	wxBusyCursor wait;
	if(ptr_et_mod == NULL) return;

	wxTextCtrl* da_field = (wxTextCtrl*) FindWindow(IDC_ETEFFH_DA_FIELD);
	if(set_var)
	{
		str = da_field->GetValue();
		bool bres = str.ToDouble(&tmp);
		if(!bres)
		{
			PrintLog(" Error reading Donor-acceptor Electrical field field \n");
		}
		else
		{
			ptr_et_mod->SetDAfield(tmp);
		}
	}
	else
	{
		tmp= ptr_et_mod->GetDAfield();
		str.Printf("%16.9f",tmp);
		str.Trim(false);
		da_field->SetValue(str);
	}	
}

void ETEffHamDlgWX::OnCalcGFfromMO(wxCommandEvent& event)
{
	wxBusyCursor wait;
	TransferDataFromWindow();
	ptr_et_mod->CalcGFDonAccOrb();	
}

void ETEffHamDlgWX::OnCalcGFfromHeff(wxCommandEvent& event)
{
	wxBusyCursor wait;
	TransferDataFromWindow();
	ptr_et_mod->CalcGFDonAccOrbHeff();
}

void ETEffHamDlgWX::OnPrintOvlpElem(wxCommandEvent& event) 
{
	wxBusyCursor wait;
	TransferDataFromWindow();
	ptr_et_mod->PrintOvlpElem();
}

void ETEffHamDlgWX::OnPrintHeffElem(wxCommandEvent& event) 
{
	wxBusyCursor wait;
	TransferDataFromWindow();
    ptr_et_mod->PrintHeffElem();
}


void ETEffHamDlgWX::OnPrintProtectMat(wxCommandEvent& event)
{
	ptr_et_mod->PrintProtectMat();
}


void ETEffHamDlgWX::OnZeroLong(wxCommandEvent& event) 
{
	double cutoff;
	wxString str;
	wxTextCtrl* cutoff_edit = (wxTextCtrl*) FindWindow(IDC_ETEFFH_CUT_DIST);
	str = cutoff_edit->GetValue();
	
	bool bres = str.ToDouble(&cutoff);

    if(!bres)
	{
		PrintLog("Error in ETEffHamDlgWX::OnZeroLong() \n");
		PrintLog(" Invalid cutoff distance \n");
		return;
	}

	ptr_et_mod->ZeroLongInter(cutoff);
}

void ETEffHamDlgWX::OnHeffEditDiag(wxCommandEvent& event) 
{
	EditHuckHamDlg* heff_edit_dlg = new EditHuckHamDlg(ptr_et_mod,NULL,-1,"Edit 1e Hamiltonian");	
	heff_edit_dlg->Show();
}

void ETEffHamDlgWX::OnSetDonorTruncMo(wxCommandEvent& event) 
{
    wxListBox* eigval_list= (wxListBox*) FindWindow(IDC_ETEFFH_EIGVAL_LIST);
	
	wxArrayInt sel_arr;

	int neig = eigval_list->GetCount();
	int idx;
	int nsel = eigval_list->GetSelections(sel_arr);

	if(nsel == 0) 
	{
		PrintLog(" No MOs are selected \n");
		return;
	}

	HaVec_int idx_sel(nsel);

	int ic = 0;
	for(idx = 0; idx < nsel; idx++)
	{
		idx_sel[idx] = sel_arr[idx]+1;	
	}
	wxBusyCursor wait;
	ptr_et_mod->AddRedoxOrbFromEigVec(idx_sel);
	TransferDataToWindow();
}

void ETEffHamDlgWX::OnCalcHeffHuck(wxCommandEvent& event) 
{
	wxBusyCursor wait;
    HaQCMod* ptr_qc_mod = ptr_et_mod->GetQCMod();
    ptr_qc_mod->InitLocOrb("VALENCE_AO");
    HaMat_double& ss = ptr_et_mod->GetActBasOvlpMat();
	ptr_qc_mod->InitHuckHam(ptr_et_mod->heff_mat,ss,*(ptr_qc_mod->ActBas));
}

void ETEffHamDlgWX::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}

/////////////////////////////////////////////////////////////////////////////
// sel_loc_orb dialog

int SelectLocOrbDlgWX::dlg_open = FALSE;

SelectLocOrbDlgWX::SelectLocOrbDlgWX(ETCouplMod* new_ptr_et_mod, int new_orb_type, wxWindow* parent):
wxDialog( parent, -1, "Select Localized Orbitals",wxDefaultPosition, wxDefaultSize, 
		   wxDEFAULT_DIALOG_STYLE )
{
	ptr_et_mod = new_ptr_et_mod;
	ptr_qc_mod = ptr_et_mod->GetQCMod();
	phost_dialog=NULL;
	orb_type = new_orb_type;
    dlg_open = TRUE;
	choose_sel_orb_dlg(this,TRUE);
}

SelectLocOrbDlgWX::~SelectLocOrbDlgWX()
{
     dlg_open = FALSE;
}


void SelectLocOrbDlgWX::DDX_lorb_list(int id)
{
	wxString label;
	wxListBox* pBox= (wxListBox*) FindWindow(IDC_LORB_LIST);
	pBox->Clear();
	int nlorb= ptr_qc_mod->GetNActiveOrb();	
	if( nlorb <= 0) return;
	for(int i= 0; i < nlorb; i++)
	{
		label = (ptr_qc_mod->ActBas->GetLabel(i)).c_str();
		pBox->Append(label);
	}
}


bool SelectLocOrbDlgWX::TransferDataToWindow()
{
	DDX_lorb_list(IDC_LORB_LIST);
	return true;
}

void SelectLocOrbDlgWX::OnExtractSel(wxCommandEvent& event)
{
	wxListBox* pBox= (wxListBox*) FindWindow(IDC_LORB_LIST);

	int i,j;
	if(phost_dialog == NULL) return;

	wxArrayInt sel_arr;

	int nn = pBox->GetSelections(sel_arr);
	if(nn == 0) return;

	LinCombOrb3D* prdx_orb = NULL;
	if( orb_type == 0) // donor
	{
		prdx_orb = &(ptr_et_mod->donor_orbs);
	}
	else
	{
		prdx_orb = &(ptr_et_mod->acc_orbs);
	}

	prdx_orb->Clear();
	prdx_orb->CreateEmptyOrbs(nn,ptr_qc_mod->ActBas);
	prdx_orb->ids.resize(nn,"ORB");

	for(i=0; i < nn; i++)
	{	
		int idx = sel_arr[i];
		for(j=0; j < nn; j++)
		{	
			prdx_orb->coef.SetVal_idx0(idx,i,1.0);
			prdx_orb->ids[i] = ptr_qc_mod->ActBas->GetLabel(idx);
		}
	}
	phost_dialog->TransferDataToWindow();
}

void SelectLocOrbDlgWX::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}

BEGIN_EVENT_TABLE(SelectLocOrbDlgWX, wxDialog)
   EVT_CLOSE(SelectLocOrbDlgWX::OnClose)
   EVT_BUTTON (IDC_EXTRACT_SEL, SelectLocOrbDlgWX::OnExtractSel)
END_EVENT_TABLE()

//////////////////////////////////////////////////////////////////////
// EditFragmDlgWX dialog

int EditFragmDlgWX::dlg_open = FALSE;

EditFragmDlgWX::EditFragmDlgWX(MolSet* new_pmset, wxWindow* parent ):
wxFrame( parent, -1, "Edit Fragments")
{
	pmset= new_pmset;

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

    edit_fragm_dlg(this,TRUE);
	OnInitDialog();
}

bool
EditFragmDlgWX::TransferDataToWindow()
{
	wxTextCtrl* parent_name_ctrl = (wxTextCtrl*) FindWindow(IDC_PARENT_MOL_SET);
	parent_name_ctrl->SetValue(pmset->GetName());

	int nf = pmset->Fragments.size();

	wxListBox* frag_list= (wxListBox*) FindWindow(IDC_FRAGM_LIST);
	frag_list->Clear();	

	int i;
	for(i = 0; i < nf; i++)
	{
		MolSet* frag = (MolSet*) pmset->Fragments[i];
		frag_list->Append(frag->GetName(), (void*) frag);
	}

	return wxFrame::TransferDataToWindow();
}

bool EditFragmDlgWX::TransferDataFromWindow()
{
	return wxFrame::TransferDataFromWindow();;
}

EditFragmDlgWX::~EditFragmDlgWX()
{
    dlg_open = FALSE;
}

void EditFragmDlgWX::OnInitDialog()
{
	dlg_open = TRUE;
	TransferDataToWindow();
}



BEGIN_EVENT_TABLE(EditFragmDlgWX, wxFrame)
	EVT_BUTTON  (IDC_EDTFRG_CREATE_FRAGM, EditFragmDlgWX::OnCreateFragm)
	EVT_BUTTON  (IDC_ADD_FRAGMENT, EditFragmDlgWX::OnAddFragment)
	EVT_BUTTON  (IDC_SELECT_ATOMS_MATCH_FRAG, EditFragmDlgWX::OnSelectAtomsMatchFrag)
	EVT_CLOSE   (EditFragmDlgWX::OnClose)
END_EVENT_TABLE()

void EditFragmDlgWX::OnCreateFragm(wxCommandEvent& event)
{
	wxTextCtrl* edit_frag_name = (wxTextCtrl*) FindWindow(IDC_EDTFRG_FRAG_NAME);
	wxString frg_name = edit_frag_name->GetValue();
	if(pmset)
	{
		MolSet* pfrag = pmset->CreateFragmentFromSelection( frg_name.ToStdString() );
		HaMainFrameWX* frame_main = GetHaMainFrameWX();
		frame_main->CreateMolView(pfrag);
	}
	TransferDataToWindow();
}

void EditFragmDlgWX::OnAddFragment(wxCommandEvent& event)
{
	wxArrayString opened_msets_names;
	VecPtr        opened_msets;

	int nm = pApp->molset_vec.size();
	int i;

	MolSet* pmset_open;

	for(i = 0; i < nm; i++)
	{
		pmset_open = (MolSet*) pApp->molset_vec[i];
		if( pmset_open == pmset || pmset->IsFragment(pmset_open)) continue;

		wxString ms_name = pmset_open->GetName();

		if( ms_name == "RESDB") continue;
		
		opened_msets.push_back((void*)pmset_open);
		opened_msets_names.Add(pmset_open->GetName());
	}

	wxMultiChoiceDialog add_fragm_dlg(NULL,"Choose Molecular Set To Add As a Fragment",
		"Opened Molecular Sets",opened_msets_names);

	add_fragm_dlg.ShowModal();

	wxArrayInt sel_idx = add_fragm_dlg.GetSelections();
	
	int nsel = sel_idx.Count();
	for(i = 0; i < nsel; i++)
	{
		int idx = sel_idx[i];
		pmset_open = (MolSet*) opened_msets[idx];
		pmset->AssociateFragment(pmset_open);
	}

	TransferDataToWindow();
}

void EditFragmDlgWX::OnSelectAtomsMatchFrag(wxCommandEvent& event)
{
	wxListBox* frag_list= (wxListBox*) FindWindow(IDC_FRAGM_LIST);

	wxArrayInt sel_idx; 
	frag_list->GetSelections(sel_idx);

	pmset->UnSelectAtomsAll();

	int nf = sel_idx.Count();
	int i;
	for(i = 0; i < nf; i++)
	{
		void* frag_data = frag_list->GetClientData(sel_idx[i]);

		MolSet* frag = (MolSet*) frag_data;
		if(frag != NULL)
		{
			pmset->SelectAtomsMatchingFragment(frag);
		}
	}
}

void EditFragmDlgWX::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}


//////////////////////////////////////////////////////////////////////
// BuildFilmDlgWX dialog

int BuildFilmDlgWX::dlg_open = FALSE;

BuildFilmDlgWX::BuildFilmDlgWX(MolSet* new_pmset, wxWindow* parent):
wxFrame( parent, -1, "Build Film")
{
	pmset= new_pmset;
	dx= 4.87;
	dy= 4.87;
	nx = 3;
	ny = 3;
	alpha= (PI/3)*RAD_TO_DEG;
	tilt = 0.0;
	nunit = 12;
	num_surf_layers = 4;

	add_surf_below_flag = false;
	add_surf_top_flag = false;
	add_atom_top_flag = false;
	add_atom_below_flag = false;

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

	build_film_dlg(this,TRUE);
	OnInitDialog();
}

BuildFilmDlgWX::~BuildFilmDlgWX()
{
        dlg_open = FALSE;
}

void BuildFilmDlgWX::OnInitDialog()
{
	wxTextCtrl* edt_ctrl;
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_FILM_DX);
	edt_ctrl->SetValidator( wxDoubleValidator(&dx, "%8.3f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_FILM_DY);
	edt_ctrl->SetValidator( wxDoubleValidator(&dy, "%8.3f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_FILM_NX);
	edt_ctrl->SetValidator( wxGenericValidator(&nx) );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_FILM_NY);
	edt_ctrl->SetValidator( wxGenericValidator(&ny) );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_FILM_ALPHA);
	edt_ctrl->SetValidator( wxDoubleValidator(&alpha, "%8.3f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_FILM_TILT);
	edt_ctrl->SetValidator( wxDoubleValidator(&tilt, "%8.3f") );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_FILM_NUNIT);
	edt_ctrl->SetValidator( wxGenericValidator(&nunit) );
	edt_ctrl = (wxTextCtrl*) FindWindow(IDC_FILM_NUM_SURF_LAYERS);
	edt_ctrl->SetValidator( wxGenericValidator(&num_surf_layers) );
	wxCheckBox* check_ctrl;
	check_ctrl = (wxCheckBox*) FindWindow(IDC_FILM_BLD_SURF_1);
	check_ctrl->SetValidator( wxGenericValidator(&add_surf_below_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_FILM_BLD_SURF_2);
	check_ctrl->SetValidator( wxGenericValidator(&add_surf_top_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_FILM_PUT_ATOM_TOP);
	check_ctrl->SetValidator( wxGenericValidator(&add_atom_top_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_FILM_PUT_ATOM_BELOW);
	check_ctrl->SetValidator( wxGenericValidator(&add_atom_below_flag) );
	
	dlg_open = TRUE;
	TransferDataToWindow();
}

BEGIN_EVENT_TABLE(BuildFilmDlgWX, wxFrame)
	EVT_BUTTON (IDC_CREATE_FILM,   BuildFilmDlgWX::OnCreateFilm)
	EVT_BUTTON (IDC_FILM_CREATE_SURF, BuildFilmDlgWX::OnCreateSurf)
	EVT_BUTTON  (IDC_CLOSE, BuildFilmDlgWX::OnCloseBtn)
	EVT_CLOSE  ( BuildFilmDlgWX::OnClose )
END_EVENT_TABLE()


void BuildFilmDlgWX::OnCreateFilm(wxCommandEvent& event)
{
	TransferDataFromWindow();
	if(pmset != NULL)
	{
		MolEditor* p_mol_editor = pmset->GetMolEditor(true);
		HaMolecule* pMol= p_mol_editor->CreateTransAlk(pmset,nunit);
		p_mol_editor->Create2DMolArray(pmset,pMol, dx, dy, nx, ny,
			                    alpha*DEG_TO_RAD, tilt*DEG_TO_RAD);

		if( add_surf_below_flag || add_surf_top_flag ||
			add_atom_top_flag || add_atom_below_flag)
		{
		   p_mol_editor->AddElectrSurf( pmset, add_surf_below_flag,add_surf_top_flag,
			                     add_atom_top_flag, add_atom_below_flag );
		}

		HaMolView* pView= pmset->GetActiveMolView();
		if(!pView)
			return;
		pView->ReDrawFlag |= RFInitial;
		pView->InitialTransform();
		pView->DefaultRepresentation();	
		pmset->RefreshAllViews();
	}
}

void BuildFilmDlgWX::OnCreateSurf(wxCommandEvent& event)
{
	TransferDataFromWindow();
	if(pmset != NULL)
	{
		MolEditor* p_mol_editor = pmset->GetMolEditor(true);
		HaMolecule* pMol = p_mol_editor->CreateSurf(pmset,num_surf_layers);

		HaMolView* pView= pmset->GetActiveMolView();
		if(!pView)
			return;
		pView->ReDrawFlag |= RFInitial;
		pView->InitialTransform();
		pView->DefaultRepresentation();	
		pmset->RefreshAllViews();
	}
}


void BuildFilmDlgWX::OnCloseBtn(wxCommandEvent& event)
{
	Close();
}

void BuildFilmDlgWX::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}

//////////////////////////////////////////////////////////////////////
// CrdSnapshotDlg dialog

int CrdSnapshotDlg::dlg_open = FALSE;

CrdSnapshotDlg::CrdSnapshotDlg(MolSet* new_pmset, wxWindow* parent):
wxFrame( parent, -1, "Coordinate Snapshots")
{
	pmset= new_pmset;
	
	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

	crd_snapshot_dlg(this,TRUE);
	OnInitDialog();
}

CrdSnapshotDlg::~CrdSnapshotDlg()
{
        dlg_open = FALSE;
}

void CrdSnapshotDlg::OnInitDialog()
{
	dlg_open = TRUE;
	snap_list = (wxGrid*) FindWindow( IDC_LIST_SNAP );
	snap_list->SetSelectionMode(wxGrid::wxGridSelectRows);

	wxCheckBox* check_ctrl;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_EN_QM);
	prop_show_flags["EN_QM"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["EN_QM"]));
	prop_col_num["EN_QM"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_EN_MM);
	prop_show_flags["EN_MM"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["EN_MM"]));
	prop_col_num["EN_MM"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_ES_QM);
	prop_show_flags["ES_QM"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["ES_QM"]));
	prop_col_num["ES_QM"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_ES_MM);
	prop_show_flags["ES_MM"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["ES_MM"]));
	prop_col_num["ES_MM"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_DS_QM);
	prop_show_flags["DS_QM"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["DS_QM"]));
	prop_col_num["DS_QM"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_DS_MM);
	prop_show_flags["DS_MM"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["DS_MM"]));
	prop_col_num["DS_MM"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_IN_QM);
	prop_show_flags["IN_QM"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["IN_QM"]));
	prop_col_num["IN_QM"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_IN_MM);
	prop_show_flags["IN_MM"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["IN_MM"]));
	prop_col_num["IN_MM"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_EN_MM_QM_DIFF);
	prop_show_flags["EN_MM_QM_DIFF"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["EN_MM_QM_DIFF"]));
	prop_col_num["EN_MM_QM_DIFF"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_ES_MM_QM_DIFF);
	prop_show_flags["ES_MM_QM_DIFF"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["EN_MM_QM_DIFF"]));
	prop_col_num["ES_MM_QM_DIFF"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_DS_MM_QM_DIFF);
	prop_show_flags["DS_MM_QM_DIFF"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["DS_MM_QM_DIFF"]));
	prop_col_num["DS_MM_QM_DIFF"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_EX_MM_QM_DIFF);
	prop_show_flags["EX_MM_QM_DIFF"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["EX_MM_QM_DIFF"]));
	prop_col_num["EX_MM_QM_DIFF"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_IN_MM_QM_DIFF);
	prop_show_flags["IN_MM_QM_DIFF"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["IN_MM_QM_DIFF"]));
	prop_col_num["IN_MM_QM_DIFF"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_BOND_ENE);
	prop_show_flags["BOND_ENE"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["BOND_ENE"]));
	prop_col_num["BOND_ENE"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_ANG_ENE);
	prop_show_flags["ANG_ENE"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["ANG_ENE"]));
	prop_col_num["ANG_ENE"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_TORS_ENE);
	prop_show_flags["TORS_ENE"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["TORS_ENE"]));
	prop_col_num["TORS_ENE"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_CNT_DIST);
	prop_show_flags["CNT_DIST"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["CNT_DIST"]));
	prop_col_num["CNT_DIST"] = -1;

	check_ctrl = (wxCheckBox*)FindWindow(IDC_SNAP_TEMP);
	prop_show_flags["TEMP"] = false;
	check_ctrl->SetValidator(wxGenericValidator(&prop_show_flags["TEMP"]));
	prop_col_num["TEMP"] = -1;

	sel_snap_id = (wxTextCtrl*) FindWindow( IDC_SEL_SNAP_ID );

	wxTextCtrl* text_ctrl;
	text_ctrl = (wxTextCtrl*) FindWindow( IDC_SNAP_MSET_NAME );
	text_ctrl->SetValue( pmset->GetName());

//	move_snap_btn = (wxSpinButton*) FindWindow( IDC_SNAP_MOVE_BTN );

	SetColumns();
	TransferDataToWindow();
}

void CrdSnapshotDlg::OnChangeProp(wxCommandEvent& event)
{
	TransferDataFromWindow();
	SetColumns();
	// FillAtomGroup();
}


void CrdSnapshotDlg::SetColumns()
{
	for ( auto& pn : prop_col_num)
	{
		pn.second = -1;
	}
	int num_cols = 1;

	for (auto& pf : prop_show_flags)
	{
		if (pf.second) num_cols++;
	}
	int tot_width = 550;
	int height = 400;
	snap_list->GetClientSize(&tot_width, &height);
	int col_width = (tot_width * 2 / 3) / num_cols;
	
	snap_list->AutoSizeColumns();
	snap_list->SetColLabelSize(50);

	int ncol_act = snap_list->GetNumberCols();
	if (ncol_act != num_cols)
	{
		snap_list->DeleteCols(0, ncol_act);
		snap_list->InsertCols(0, num_cols);
	}
	snap_list->SetColLabelValue(0, "CH");
	snap_list->SetColFormatBool(0);

	int icol = -1;
	for (auto& pf : prop_show_flags)
	{
		if (pf.second)
		{
			icol++;
			std::string prop_lbl = pf.first;
			snap_list->SetColLabelValue(icol, prop_lbl);
			snap_list->SetColMinimalWidth(icol, col_width);
			prop_col_num[prop_lbl] = icol;
		}
	}
}

BEGIN_EVENT_TABLE(CrdSnapshotDlg, wxFrame)
	EVT_BUTTON (IDC_ADD_SNAPSHOT,   CrdSnapshotDlg::OnAddSnapshot)
	EVT_BUTTON (IDC_DEL_SNAPSHOT,   CrdSnapshotDlg::OnDelSnapshot)
	EVT_BUTTON (IDC_DEL_ALL_SNAPSHOTS,   CrdSnapshotDlg::OnDelAllSnapshots)
	EVT_BUTTON (IDC_SET_CRD_FROM_SNAPSHOT,  CrdSnapshotDlg::OnSetCrdFromSnapshot)
	EVT_BUTTON (IDC_SAVE_CRD_TO_SNAPSHOT,   CrdSnapshotDlg::OnSaveCrdToSnapshot)
	EVT_BUTTON (IDC_LOAD_SNAP_FROM_MOL_FILE,   CrdSnapshotDlg::OnLoadSnapshotFromMolFile )
	EVT_BUTTON (IDC_LOAD_SNAPS_FROM_XML_FILE, CrdSnapshotDlg::OnLoadSnapshotsFromXMLFile )
	EVT_BUTTON (IDC_SAVE_SNAPS_TO_XML_FILE,   CrdSnapshotDlg::OnSaveSnapshotsToXMLFile )
	EVT_CHECKBOX(IDC_SNAP_EN_QM, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_EN_MM, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_ES_QM, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_ES_MM, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_DS_QM, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_DS_MM, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_IN_QM, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_IN_MM, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_CNT_DIST, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_TEMP, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_EN_MM_QM_DIFF, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_ES_MM_QM_DIFF, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_DS_MM_QM_DIFF, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_EX_MM_QM_DIFF, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_IN_MM_QM_DIFF, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_BOND_ENE, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_ANG_ENE, CrdSnapshotDlg::OnChangeProp)
	EVT_CHECKBOX(IDC_SNAP_TORS_ENE, CrdSnapshotDlg::OnChangeProp)

//	EVT_SPIN_UP( IDC_SNAP_MOVE_BTN, CrdSnapshotDlg::OnSnapMoveBtnUp )
//	EVT_SPIN_DOWN( IDC_SNAP_MOVE_BTN, CrdSnapshotDlg::OnSnapMoveBtnDown )
//	EVT_TEXT( IDC_SEL_SNAP_DESC, CrdSnapshotDlg::OnEditDescription )
	EVT_TEXT( IDC_SEL_SNAP_ID,   CrdSnapshotDlg::OnEditSnapID )
//	EVT_LISTBOX(IDC_LIST_SNAP, CrdSnapshotDlg::OnChangeSelSnapshot )
	EVT_CLOSE  ( CrdSnapshotDlg::OnClose )
END_EVENT_TABLE()


void CrdSnapshotDlg::OnAddSnapshot(wxCommandEvent& event)
{
	wxString snap_id = sel_snap_id->GetValue();
	CrdSnapshot* psnap = pmset->AddCrdSnapshot( snap_id.ToStdString() );
//	int ns = snap_list->GetItemCount();
//	if( psnap != NULL )
	//{
	//	long idx = snap_list->InsertItem( ns + 1, (wxString) psnap->GetName() );
	//	snap_list->SetItemPtrData(idx, (wxUIntPtr) psnap);

	//	snap_list->SetItemState(idx, wxLIST_STATE_SELECTED, wxLIST_STATE_SELECTED );
	//	OnChangeSelSnapshot(event);
	//}
}
  
void CrdSnapshotDlg::OnDelSnapshot(wxCommandEvent& event)
{
	CrdSnapshot* psnap = GetSelSnapshot();
	if( psnap == NULL )
	{
		PrintLog(" No Snapshot Selected() \n");
		return;
	}
	int ires = pmset->DeleteCrdSnapshot(psnap);
	if( ires )
	{
		wxArrayInt sel_idx_arr = snap_list->GetSelectedRows();
		int idx = sel_idx_arr[0];
		snap_list->DeleteRows(idx);
		int n = snap_list->GetNumberRows(); 
		if( idx < n )
		{
			snap_list->SelectRow( idx );
		}
		else if( n > 0 )
		{
			snap_list->SelectRow( idx - 1 );
		}
		OnChangeSelSnapshot( event );
	}
	TransferDataToWindow();
}

void CrdSnapshotDlg::OnDelAllSnapshots(wxCommandEvent& event)
{
	pmset->DeleteCrdSnapshots();
	TransferDataToWindow();
}

void CrdSnapshotDlg::OnSetCrdFromSnapshot(wxCommandEvent& event)
{
	CrdSnapshot* psnap = GetSelSnapshot();
	if( psnap != NULL )
	{
		psnap->SetAtomCrd();
		pmset->RefreshAllViews(RFApply | RFRefresh); 
	}
}

void CrdSnapshotDlg::OnSaveCrdToSnapshot(wxCommandEvent& event)
{
	CrdSnapshot* psnap = GetSelSnapshot();
	if( psnap != NULL )
	{
		psnap->SaveCurrentAtomCrd(); 
	}
}

void  CrdSnapshotDlg::OnLoadSnapshotFromMolFile(wxCommandEvent& event)
{
	ChooseMolFileDlg load_dlg(NULL, -1,  "Choose Molecular File");
	load_dlg.ShowModal();

	int na = pmset->GetNAtoms();

	CrdSnapshot* psnap = NULL;
    if( !load_dlg.file_name.empty() )
    {
		MolSet* p_cur_mset_old = GetCurMolSet();
		std::vector<std::string> tokens;
		boost::split(tokens, load_dlg.file_name, boost::is_any_of(";"));
		for (std::string fname : tokens)
		{
			boost::trim(fname);
			std::filesystem::path dir_path(load_dlg.dir_name);
			std::filesystem::path file_local_path(fname);
			std::filesystem::path file_full_path = dir_path / file_local_path;

			std::shared_ptr<MolSet> pmset_load(new MolSet());
			
			pmset_load->FetchFile(load_dlg.file_format, file_full_path.string().c_str() );
			int na_load = pmset_load->GetNAtoms();
			if (na_load != na)
			{
				PrintLog(" Error in CrdSnapshotDlg::OnLoadSnapshotFromMolFile() \n");
				PrintLog(" Number of atoms in the loaded file = %d  not equal to number of atoms in the current molset = %d \n", na_load, na);
				return;
			}
			HaVec_double crd(3 * na);
			AtomIteratorMolSet aitr(pmset_load.get());
			int i = 0;
			HaAtom* aptr;
			for (aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
			{
				crd[3 * i] = aptr->GetX();
				crd[3 * i + 1] = aptr->GetY();
				crd[3 * i + 2] = aptr->GetZ();
				i++;
			}

			std::string snap_name = fname;
			size_t lastindex = fname.find_last_of(".");
			if (lastindex > 0) snap_name = fname.substr(0, lastindex);
			psnap = pmset->AddCrdSnapshot(snap_name);

			PrintLog(" Coordinate Snapshot %s is set from molecular file %s \n", psnap->GetName().c_str(), fname.c_str());

			psnap->SaveCrd(crd);
			if (pmset_load->per_bc->IsSet())
			{
				HaVec_double pbox(6);
				pbox[0] = pmset_load->per_bc->GetA();
				pbox[1] = pmset_load->per_bc->GetB();
				pbox[2] = pmset_load->per_bc->GetC();
				pbox[3] = pmset_load->per_bc->GetAlpha();
				pbox[4] = pmset_load->per_bc->GetBeta();
				pbox[5] = pmset_load->per_bc->GetGamma();
				psnap->SavePBox(pbox);
			}
		}
		SetCurMolSet(p_cur_mset_old);
		TransferDataToWindow();
    }
}

void CrdSnapshotDlg::OnLoadSnapshotsFromXMLFile(wxCommandEvent& event)
{
	PrintLog(" CrdSnapshotDlg::OnLoadSnapshotsFromXMLFile() \n");
	ChooseMolFileDlg load_dlg(NULL, -1,  "Choose File with Coordinate Snapshots");
	load_dlg.file_types_ch->Clear();
	load_dlg.file_types_ch->Append("HARLEM format(*.hlm)",(void*)"*.hlm");
	load_dlg.file_types_ch->SetSelection(0);
	load_dlg.choose_file_btn->SetLabel("Load File");
	load_dlg.ShowModal();
	if( !load_dlg.file_name.empty() )
    {
		PrintLog("Loading Coordinate Snapshots from file %s \n",load_dlg.file_name.c_str());
		pmset->LoadCrdSnapshots( load_dlg.file_name );
		TransferDataToWindow();
	}
}

void CrdSnapshotDlg::OnSaveSnapshotsToXMLFile(wxCommandEvent& event)
{
	PrintLog(" CrdSnapshotDlg::OnSaveSnapshotsToXMLFile() \n");
	ChooseMolFileDlg save_dlg(NULL, -1,  "Choose Molecular File");
	save_dlg.file_types_ch->Clear();
	save_dlg.file_types_ch->Append("HARLEM format(*.hlm)",(void*)"*.hlm");
	save_dlg.file_types_ch->SetSelection(0);
	save_dlg.choose_file_btn->SetLabel("Save to File");
	save_dlg.ShowModal();
	
	if( !save_dlg.file_name.empty() )
    {
		pmset->SaveCrdSnapshots( save_dlg.file_name );
	}
}

void  CrdSnapshotDlg::OnSnapMoveBtnUp( wxSpinEvent& event )
{
	/*int ig_cur = snap_list->GetSelection();
	if( ig_cur < 1 ) return;
	if( ig_cur >= pmset->crd_snapshots.size() ) return;
	CrdSnapshot* snap1 = pmset->crd_snapshots[ig_cur-1];
	CrdSnapshot* snap2 = pmset->crd_snapshots[ig_cur];

	pmset->crd_snapshots[ig_cur-1] = snap2;
	pmset->crd_snapshots[ig_cur]   = snap1;

	TransferDataToWindow();
	snap_list->SetSelection( ig_cur-1 );
	wxCommandEvent event2;
	OnChangeSelSnapshot(event2);*/
}

void  CrdSnapshotDlg::OnSnapMoveBtnDown( wxSpinEvent& event )
{
	/*int ig_cur = snap_list->GetSelection();
	if( ig_cur < 0 ) return;
	if( ig_cur >= (snap_list->GetCount() - 1)) return;
	if( ig_cur >= pmset->crd_snapshots.size() ) return;

	CrdSnapshot* snap1 = pmset->crd_snapshots[ig_cur];
	CrdSnapshot* snap2 = pmset->crd_snapshots[ig_cur+1];

	pmset->crd_snapshots[ig_cur] = snap2;
	pmset->crd_snapshots[ig_cur+1]   = snap1;

	TransferDataToWindow();
	snap_list->SetSelection( ig_cur+1 );
	wxCommandEvent event2;
	OnChangeSelSnapshot(event2);*/
}

void CrdSnapshotDlg::OnChangeNumSnapshots()
{
	int nsnap = pmset->crd_snapshots.size();
	if( nsnap > 0 )
	{
		move_snap_btn->SetRange(0,nsnap-1);
	}
	else
	{
		move_snap_btn->SetRange(0,0);
	}
}

void CrdSnapshotDlg::OnChangeSelSnapshot(wxCommandEvent& event)
{
	char buf[80];
	long item = -1;
	//long last_sel_idx = snap_list->GetNextItem(item, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);
	//for (;;)
	//{
	//	long idx = snap_list->GetNextItem(item, wxLIST_NEXT_ALL, wxLIST_STATE_SELECTED);
	//	if (idx == -1) break;
	//	last_sel_idx = idx;
	//}
	//move_snap_btn->SetValue(last_sel_idx);
	//CrdSnapshot* psnap = GetSelSnapshot();
	//if( psnap != NULL )
	//{
	//	sel_snap_id->ChangeValue( psnap->GetName() );
	//	sel_snap_desc->ChangeValue( psnap->GetDesc() );

	//	HaVec_double crd = psnap->GetCrd();
	//	HaVec_double avg_crd(3);
	//	int n = crd.size()/3;
	//	if( n > 0 )
	//	{
	//		int i;
	//		for( i = 0; i < n; i++ )
	//		{
	//			avg_crd[0] += crd[3*i];
	//			avg_crd[1] += crd[3*i+2];
	//			avg_crd[2] += crd[3*i+2];
	//		}
	//		avg_crd[0] = avg_crd[0]/n;
	//		avg_crd[1] = avg_crd[1]/n;
	//		avg_crd[2] = avg_crd[2]/n;

	//		sprintf(buf," %16.9f  %16.9f  %16.9f ", avg_crd[0], avg_crd[1], avg_crd[2] ); 
	//		avg_crd_txt->SetValue( buf );
	//	}
	//	else
	//	{
	//		avg_crd_txt->SetValue(" Invalid avg crd: nat = 0 ");
	//	}
	//}
	//else
	//{
	//	sel_snap_id->ChangeValue( "" );
	//	sel_snap_desc->ChangeValue( "" );
	//	avg_crd_txt->SetValue("");
	//}
}

void CrdSnapshotDlg::OnEditSnapID( wxCommandEvent& event )
{
	CrdSnapshot* psnap = GetSelSnapshot();
	if( psnap != NULL )
	{
		wxArrayInt sel_idx_arr = snap_list->GetSelectedRows();
		int ig_cur = sel_idx_arr[0];
		std::string id = sel_snap_id->GetValue().ToStdString();
		psnap->SetName( id );
		snap_list->SetRowLabelValue(ig_cur, id );
	}
}

CrdSnapshot*  CrdSnapshotDlg::GetSelSnapshot()
{
	if(pmset == NULL) return NULL;
	wxArrayInt sel_idx_arr = snap_list->GetSelectedRows();
	if (sel_idx_arr.GetCount() < 1) return NULL;

	int sel_idx = sel_idx_arr[0];
	CrdSnapshot* psnap = pmset->crd_snapshots[sel_idx];
	return psnap;
}

void CrdSnapshotDlg::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}

bool CrdSnapshotDlg::TransferDataFromWindow()
{
	return wxFrame::TransferDataFromWindow();
}

bool CrdSnapshotDlg::TransferDataToWindow()
{
	//snap_list->ClearGrid();
	int i;
	int ns = pmset->crd_snapshots.size();
	int ncol = snap_list->GetNumberCols();
	int nrow = snap_list->GetNumberRows();

	if (nrow != ns )
	{
		if (nrow > 0) snap_list->DeleteRows(0, nrow);
		snap_list->AppendRows(ns);
	}

	for( i = 0; i < ns; i++ )
	{
		CrdSnapshot* psnap = pmset->crd_snapshots[i];
		snap_list->SetRowLabelValue(i, (wxString)psnap->GetName());
	}
	snap_list->SetRowLabelSize(wxGRID_AUTOSIZE);
	snap_list->AutoSize();
	// snap_list->SetColLabelSize(50);
	return wxFrame::TransferDataToWindow();
}

//////////////////////////////////////////////////////////////////////
// AtomParamsDlgWX dialog

int   AtomParamsDlgWX::dlg_open = FALSE;
bool  AtomParamsDlgWX::coord_edit_flag = false;
bool  AtomParamsDlgWX::element_edit_flag = true;
bool  AtomParamsDlgWX::charge_edit_flag = true;
bool  AtomParamsDlgWX::ff_symbol_edit_flag = false;
bool  AtomParamsDlgWX::mass_edit_flag = false;
bool  AtomParamsDlgWX::atname_edit_flag = true;
bool  AtomParamsDlgWX::radius_edit_flag = false;
bool  AtomParamsDlgWX::vdwrad_edit_flag = false;
bool  AtomParamsDlgWX::vdwene_edit_flag = false;
bool  AtomParamsDlgWX::temperature_edit_flag = false;
bool  AtomParamsDlgWX::hybridization_edit_flag = false;
bool  AtomParamsDlgWX::hb_status_edit_flag = false;
bool  AtomParamsDlgWX::solv_access_area_flag = false;


AtomParamsDlgWX::AtomParamsDlgWX(MolSet* new_pmset, wxWindow* parent):
wxFrame( parent, -1, "Edit Atom Parameters")
{
	pmset= new_pmset;
	p_mol_editor   = pmset->GetMolEditor(true);
	p_prot_rdx_mod = pmset->GetProtonRedoxMod(true); 

	n_x_coord = -1;
	n_y_coord = -1;
	n_z_coord = -1;
	n_element = -1;
	n_charge = -1;
	n_ff_symbol = -1;
	n_mass = -1;
	n_atname = -1;
	n_radius = -1;
	n_vdwrad = -1;
	n_vdwene = -1;
	n_temp = -1;
	n_hybrid = -1;
	n_hb_status = -1;

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

	wxMenuBar* atom_params_menu_bar = atom_params_menu();
    this->SetMenuBar(atom_params_menu_bar);    

	atom_params_dlg(this,TRUE);
	m_atom_lctrl = (wxGrid*) FindWindow( IDC_ATOM_LIST );
	OnInitDialog();     
	dlg_open = TRUE;
	
}

AtomParamsDlgWX::~AtomParamsDlgWX()
{
    dlg_open = FALSE;
}


void AtomParamsDlgWX::OnInitDialog()
{	
	wxCheckBox* check_ctrl;
	
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_COORD);
	check_ctrl->SetValidator( wxGenericValidator(&coord_edit_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_ELEMENT);
	check_ctrl->SetValidator( wxGenericValidator(&element_edit_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_CHARGE);
	check_ctrl->SetValidator( wxGenericValidator(&charge_edit_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_FF_SYMB);
	check_ctrl->SetValidator( wxGenericValidator(&ff_symbol_edit_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_MASS);
	check_ctrl->SetValidator( wxGenericValidator(&mass_edit_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_NAME);
	check_ctrl->SetValidator( wxGenericValidator(&atname_edit_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_RADIUS);
	check_ctrl->SetValidator( wxGenericValidator(&radius_edit_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_COORD);
	check_ctrl->SetValidator( wxGenericValidator(&coord_edit_flag) );	
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_VDWRAD);
	check_ctrl->SetValidator( wxGenericValidator(&vdwrad_edit_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_VDWENE);
	check_ctrl->SetValidator( wxGenericValidator(&vdwene_edit_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_TEMP);
	check_ctrl->SetValidator( wxGenericValidator(&temperature_edit_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_HYBRID);
	check_ctrl->SetValidator( wxGenericValidator(&hybridization_edit_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_HB_STATUS);
	check_ctrl->SetValidator( wxGenericValidator(&hb_status_edit_flag) );
	check_ctrl = (wxCheckBox*) FindWindow(IDC_EDTAT_SOLV_ACCESS_AREA);
	check_ctrl->SetValidator( wxGenericValidator(&solv_access_area_flag) );
	wxTextCtrl* edit_ph = (wxTextCtrl*) FindWindow(IDC_EDTAT_PH);
	edit_ph->SetValidator( wxDoubleValidator(&p_prot_rdx_mod->ph,"%6.2f") );

	SetColumns();
	FillAtomGroup();
	TransferDataToWindow();
}

int AtomParamsDlgWX::ResetEditFlags()
{
	coord_edit_flag = FALSE;
	element_edit_flag = FALSE;
	charge_edit_flag = FALSE;
	ff_symbol_edit_flag = FALSE;
	mass_edit_flag = FALSE;
	atname_edit_flag = FALSE;
	radius_edit_flag = FALSE;
	vdwrad_edit_flag = FALSE;
	vdwene_edit_flag = FALSE;
	temperature_edit_flag = FALSE;
	hybridization_edit_flag = FALSE;
	hb_status_edit_flag = FALSE;
	solv_access_area_flag = FALSE;

	return TRUE;
}

void AtomParamsDlgWX::SetColumns()
{
	n_atname = -1;
	n_element = -1;
	n_ff_symbol = -1;
	n_mass = -1;
	n_charge = -1;
	n_radius = -1;
	n_vdwrad = -1;
	n_vdwene = -1;
	n_temp = -1;
	n_x_coord = -1;
	n_y_coord = -1;
	n_z_coord = -1;
	n_hybrid = -1;
	n_hb_status = -1;
	n_solv_access = -1;

	int num_cols = 0;

	if( coord_edit_flag ) num_cols += 3;
	if( element_edit_flag ) num_cols++;
	if( charge_edit_flag ) num_cols++;
	if( ff_symbol_edit_flag ) num_cols++;
	if( mass_edit_flag ) num_cols++;
	if( atname_edit_flag ) num_cols++;
	if( radius_edit_flag ) num_cols++;
	if( vdwrad_edit_flag ) num_cols++;
	if( vdwene_edit_flag ) num_cols++;
	if( temperature_edit_flag ) num_cols++;
	if( hybridization_edit_flag ) num_cols++;
	if( hb_status_edit_flag ) num_cols++;
	if( solv_access_area_flag ) num_cols++;

	int tot_width = 550;
	int height = 400;
	m_atom_lctrl->GetClientSize(&tot_width, &height);

	int col_width = (tot_width*2 /3) /num_cols;

	int icol = -1;
	
	m_atom_lctrl->AutoSizeColumns();
	m_atom_lctrl->SetColLabelSize(21);

	int ncol_act = m_atom_lctrl->GetNumberCols();
	if(ncol_act != num_cols)
	{
		m_atom_lctrl->DeleteCols(0,ncol_act);
		m_atom_lctrl->InsertCols(0,num_cols);
	}
 
	if( atname_edit_flag )
	{
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "Name");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_atname = icol;
	}

	if( element_edit_flag )
	{
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "Elem");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_element = icol;
	}

	if( ff_symbol_edit_flag )
	{
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "FF Smbl");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_ff_symbol = icol;
	}

	if( mass_edit_flag )
	{
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "Mass");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_mass = icol;
	}

	if( charge_edit_flag )
	{
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "Charge");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_charge = icol;
	}

	if( radius_edit_flag )
	{
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "Radius");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_radius = icol;
	}

	if( vdwrad_edit_flag )
	{
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "VdW rad");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_vdwrad = icol;
	}

	if( vdwene_edit_flag )
	{
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "VdW ene");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_vdwene = icol;
	}

	if( temperature_edit_flag )
	{
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "Temp");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_temp = icol;
	}

	if( hybridization_edit_flag )
	{
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "Hybrid");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_hybrid = icol;
	}

	if( hb_status_edit_flag )
	{
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "HB status");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_hb_status = icol;
	}

	if( solv_access_area_flag )
	{
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "Solv Access");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_solv_access = icol;
	}


	if( coord_edit_flag )
	{
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "X");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_x_coord = icol;
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "Y");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_y_coord = icol;
		icol++;
        m_atom_lctrl->SetColLabelValue(icol, "Z");
	    m_atom_lctrl->SetColMinimalWidth(icol, col_width);
		n_z_coord = icol;
	}
}

void AtomParamsDlgWX::FillAtomGroup()
{
	char buf[256];

	m_atom_lctrl->ClearGrid();

	if(pmset == NULL) return;
	HaAtom* aptr;

	int nat_sel = 0;
	AtomIteratorMolSet aitr(pmset);
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if(!aptr->Selected()) continue;
		nat_sel++;
	}

	int ncol = m_atom_lctrl->GetNumberCols();
	int nrow = m_atom_lctrl->GetNumberRows();

	int nv = at_ptrs.size();

	if( nrow != nat_sel || nv != nat_sel )
	{
		if( nrow > 0) m_atom_lctrl->DeleteRows(0, nrow);
        m_atom_lctrl->AppendRows(nat_sel);
		at_ptrs.resize(nat_sel);
	}

	int itm = -1;
	int max_lbl = 0;

	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		if(!aptr->Selected()) continue;
		itm++;
		aptr->FillRef(buf);
		at_ptrs[itm] = (void*) aptr;

		int len = strlen(buf);
		if(len > max_lbl) max_lbl = len;
	
		m_atom_lctrl->SetRowLabelValue(itm, buf); 

		if( n_charge >= 0 )
		{
			sprintf(buf,"%9.4f",aptr->GetCharge());
			m_atom_lctrl->SetCellValue(itm, n_charge, buf); 
		}

		if( n_radius >= 0 )
		{
			sprintf(buf,"%9.4f",aptr->radius );
			m_atom_lctrl->SetCellValue(itm, n_radius, buf); 
		}

		if( n_vdwrad >= 0 )
		{
			sprintf(buf,"%9.4f",aptr->GetVdWRad() );
			m_atom_lctrl->SetCellValue(itm, n_vdwrad, buf); 
		}

		if( n_vdwene >= 0 )
		{
			sprintf(buf,"%9.4f",aptr->GetVdWEne());
			m_atom_lctrl->SetCellValue(itm, n_vdwene, buf);
		}

		if( n_temp >= 0 )
		{
			sprintf(buf,"%9.4f",aptr->tempf);
			m_atom_lctrl->SetCellValue(itm, n_temp, buf);
		}

		if( n_element >= 0 )
		{
			sprintf(buf,"%d",aptr->GetElemNo() );
			m_atom_lctrl->SetCellValue(itm, n_element, buf);
		}

		if( n_mass >= 0 )
		{
			sprintf(buf,"%9.4f",aptr->GetMass() );
			m_atom_lctrl->SetCellValue(itm, n_mass, buf);
		}

		if( n_x_coord >= 0 )
		{
			sprintf(buf,"%9.4f",aptr->GetX_Ang());
			m_atom_lctrl->SetCellValue(itm, n_x_coord, buf);

		}

		if( n_y_coord >= 0 )
		{
			sprintf(buf,"%9.4f",aptr->GetY_Ang() );
			m_atom_lctrl->SetCellValue(itm, n_y_coord, buf);
		}

		if( n_z_coord >= 0 )
		{
			sprintf(buf,"%9.4f",aptr->GetZ_Ang() );
			m_atom_lctrl->SetCellValue(itm, n_z_coord, buf);
		}

		if( n_ff_symbol >= 0 )
		{
			sprintf(buf,"%s", aptr->GetFFSymbol().c_str() );
			m_atom_lctrl->SetCellValue(itm, n_ff_symbol, buf);
		}
		
		if( n_atname >= 0 )
		{
			strcpy(buf,aptr->GetName());
			m_atom_lctrl->SetCellValue(itm, n_atname, buf);
		}

		if( n_hybrid >= 0 )
		{
			std::string hybrid_str;
			aptr->GetHybridTextStr(hybrid_str);
			strcpy(buf,hybrid_str.c_str() );
			m_atom_lctrl->SetCellValue(itm, n_hybrid, buf);
		}

		if( n_hb_status >= 0 )
		{
			std::string hb_status_str;
			aptr->GetHBStatusTextStr(hb_status_str);
			strcpy(buf,hb_status_str.c_str());
			m_atom_lctrl->SetCellValue(itm, n_hb_status, buf);
		}

		if( n_solv_access >= 0 )
		{
			sprintf(buf,"%9.4f",aptr->solv_access_area );
			m_atom_lctrl->SetCellValue(itm, n_solv_access, buf);
		}
	}
	m_atom_lctrl->SetRowLabelSize((int)(max_lbl*8.5));
	m_atom_lctrl->AutoSizeColumns();
}

bool AtomParamsDlgWX::TransferDataFromWindow()
{
	wxCheckBox* check_pbox = (wxCheckBox*) FindWindow(IDC_CHECK_PER_BOX);
	wxTextCtrl* txt_pbox_a = (wxTextCtrl*) FindWindow(IDC_PER_BOX_A);
	wxTextCtrl* txt_pbox_b = (wxTextCtrl*) FindWindow(IDC_PER_BOX_B);
	wxTextCtrl* txt_pbox_c = (wxTextCtrl*) FindWindow(IDC_PER_BOX_C);
	wxTextCtrl* txt_pbox_alpha  = (wxTextCtrl*) FindWindow(IDC_PER_BOX_ALPHA);
	wxTextCtrl* txt_pbox_beta   = (wxTextCtrl*) FindWindow(IDC_PER_BOX_BETA);
	wxTextCtrl* txt_pbox_gamma  = (wxTextCtrl*) FindWindow(IDC_PER_BOX_GAMMA);

	wxString str;
	if( check_pbox->GetValue() )
	{
		double a,b,c,alpha,beta,gamma;

		int num_valid = TRUE;
		bool bres;

		str = txt_pbox_a->GetValue();
		bres = str.ToDouble(&a);
		if (!bres || a < 0.0001) num_valid = FALSE;

		str = txt_pbox_b->GetValue();
		bres = str.ToDouble(&b);
		if (!bres || b < 0.0001) num_valid = FALSE;

		str = txt_pbox_c->GetValue();
		bres = str.ToDouble(&c);
		if (!bres || c < 0.0001) num_valid = FALSE;

		str = txt_pbox_alpha->GetValue();
		bres = str.ToDouble(&alpha);
		if (!bres || alpha < 0.0001) num_valid = FALSE;

		str = txt_pbox_beta->GetValue();
		bres = str.ToDouble(&beta);
		if (!bres || beta < 0.0001) num_valid = FALSE;

		str = txt_pbox_gamma->GetValue();
		bres = str.ToDouble(&gamma);
		if (!bres || gamma < 0.0001) num_valid = FALSE;

		if( num_valid )
		{
			alpha *= DEG_TO_RAD;
			beta  *= DEG_TO_RAD;
			gamma *= DEG_TO_RAD;

			int set_new = FALSE;

			if( fabs( a - pmset->per_bc->GetA()) > 0.00000001) set_new = TRUE;
			if( fabs( b - pmset->per_bc->GetB()) > 0.00000001) set_new = TRUE;
			if( fabs( c - pmset->per_bc->GetC()) > 0.00000001) set_new = TRUE;
			if( fabs( alpha - pmset->per_bc->GetAlpha()) > 0.00000001) set_new = TRUE;
			if( fabs( beta  - pmset->per_bc->GetBeta())  > 0.00000001) set_new = TRUE;
			if( fabs( gamma - pmset->per_bc->GetGamma()) > 0.00000001) set_new = TRUE;

			if( set_new )
			{
				pmset->per_bc->SetBox(a,b,c,alpha,beta,gamma);
			}
		}
	}
	return wxFrame::TransferDataFromWindow();
}

bool AtomParamsDlgWX::TransferDataToWindow()
{
	wxComboBox* chmap_list= (wxComboBox*) FindWindow(IDC_EDTAT_CHARGE_MAPS);
	chmap_list->Clear();
	if(pmset != NULL)
	{
		int nm = pmset->ChargeMaps.size();
		int i;
		for(i=0; i < nm; i++)
		{
			chmap_list->Append(pmset->ChargeMaps[i].GetName());
		}
	}
	wxCheckBox* check_pbox = (wxCheckBox*) FindWindow(IDC_CHECK_PER_BOX);
	wxTextCtrl* txt_pbox_a = (wxTextCtrl*) FindWindow(IDC_PER_BOX_A);
	wxTextCtrl* txt_pbox_b = (wxTextCtrl*) FindWindow(IDC_PER_BOX_B);
	wxTextCtrl* txt_pbox_c = (wxTextCtrl*) FindWindow(IDC_PER_BOX_C);
	wxTextCtrl* txt_pbox_alpha  = (wxTextCtrl*) FindWindow(IDC_PER_BOX_ALPHA);
	wxTextCtrl* txt_pbox_beta   = (wxTextCtrl*) FindWindow(IDC_PER_BOX_BETA);
	wxTextCtrl* txt_pbox_gamma  = (wxTextCtrl*) FindWindow(IDC_PER_BOX_GAMMA);

	wxTextCtrl* txt_num_atoms  = (wxTextCtrl*) FindWindow(IDC_NUM_ATOMS);
	txt_num_atoms->SetValue( wxString::Format("%d",pmset->GetNAtoms()) );

	if( pmset->per_bc->IsSet() )
	{
		check_pbox->SetValue(true);
		txt_pbox_a->SetEditable(true);
		txt_pbox_a->SetBackgroundColour(*wxWHITE);
		txt_pbox_a->SetValue( wxString::Format("%16.9f",pmset->per_bc->GetA()) );
		txt_pbox_b->SetEditable(true);
		txt_pbox_b->SetBackgroundColour(*wxWHITE);
		txt_pbox_b->SetValue( wxString::Format("%16.9f",pmset->per_bc->GetB()) );
		txt_pbox_c->SetEditable(true);
		txt_pbox_c->SetBackgroundColour(*wxWHITE);
		txt_pbox_c->SetValue( wxString::Format("%16.9f",pmset->per_bc->GetC()) );
		txt_pbox_alpha->SetEditable(true);
		txt_pbox_alpha->SetBackgroundColour(*wxWHITE);
		txt_pbox_alpha->SetValue( wxString::Format("%16.9f",pmset->per_bc->GetAlpha() * RAD_TO_DEG) );
		txt_pbox_beta->SetEditable(true);
		txt_pbox_beta->SetBackgroundColour(*wxWHITE);
		txt_pbox_beta->SetValue( wxString::Format("%16.9f",pmset->per_bc->GetBeta() * RAD_TO_DEG) );
		txt_pbox_gamma->SetEditable(true);
		txt_pbox_gamma->SetBackgroundColour(*wxWHITE);
		txt_pbox_gamma->SetValue( wxString::Format("%16.9f",pmset->per_bc->GetGamma() * RAD_TO_DEG) );
	}
	else
	{
		check_pbox->SetValue(false);
		txt_pbox_a->SetEditable(false);
		txt_pbox_a->SetBackgroundColour(*wxLIGHT_GREY);
		txt_pbox_a->SetValue("");
		txt_pbox_b->SetEditable(false);
		txt_pbox_b->SetBackgroundColour(*wxLIGHT_GREY);
		txt_pbox_b->SetValue("");
		txt_pbox_c->SetEditable(false);
		txt_pbox_c->SetBackgroundColour(*wxLIGHT_GREY);
		txt_pbox_c->SetValue("");
		txt_pbox_alpha->SetEditable(false);
		txt_pbox_alpha->SetBackgroundColour(*wxLIGHT_GREY);
		txt_pbox_alpha->SetValue("");
		txt_pbox_beta->SetEditable(false);
		txt_pbox_beta->SetBackgroundColour(*wxLIGHT_GREY);
		txt_pbox_beta->SetValue("");
		txt_pbox_gamma->SetEditable(false);
		txt_pbox_gamma->SetBackgroundColour(*wxLIGHT_GREY);
		txt_pbox_gamma->SetValue("");
	}
	return wxFrame::TransferDataToWindow();
}


void AtomParamsDlgWX::OnCopyToChmap(wxCommandEvent& event) 
{
	if(pmset == NULL)
		return;

	wxComboBox* chmap_list= (wxComboBox*) FindWindow(IDC_EDTAT_CHARGE_MAPS);
	wxString str = chmap_list->GetValue();
	
	str.Trim(false);
	str.Trim(true);

	if( str.IsEmpty() )
	{
		ErrorInMod("AtomParamsDlgWX::OnCopyToChmap()",
			       "Empty name of Charge Map");
		return;
	}

	pmset->SetChargeMapByCurrentCharges(str.c_str());
	TransferDataToWindow();
}


void AtomParamsDlgWX::OnSetFromChmap(wxCommandEvent& event) 
{
	if(pmset == NULL)
	return;


	wxComboBox* chmap_list= (wxComboBox*) FindWindow(IDC_EDTAT_CHARGE_MAPS);
	wxString str = chmap_list->GetValue();
	
	str.Trim(false);
	str.Trim(true);
	
	if(str.IsEmpty())
	{
		ErrorInMod("AtomParamsDlgWX::OnSetFromChmap()",
			       "Empty name of Charge Map");
		return;
	}

	AtomDoubleMap* charge_map = pmset->GetChargeMapByName(str.c_str());
	if(charge_map == NULL)
	{
		PrintLog(" No Charge Map with the name %s \n",str.ToStdString().c_str());
		return;
	}
	pmset->SetChargesFromChargeMap(charge_map);
	TransferDataToWindow();	
	FillAtomGroup();
}


BEGIN_EVENT_TABLE(AtomParamsDlgWX, wxFrame)
    EVT_MENU   (IDC_SET_STD_CHARGES, AtomParamsDlgWX::OnEdtAtAmberCh)
	EVT_MENU   (IDC_EDTAT_SET_HBDA_STATUS, AtomParamsDlgWX::OnEdtAtHBondStatus)
	EVT_MENU   (IDC_CALC_SOLV_ACCESS_AREA, AtomParamsDlgWX::OnCalcSolvAccessArea)
	EVT_MENU   (IDC_CALC_SOLV_ACCESS_AREA_GEOBALL, AtomParamsDlgWX::OnCalcSolvAccessAreaGeoball)
	EVT_MENU   (IDC_SAVE_CRD_EXCL_VOL_ARB, AtomParamsDlgWX::OnSaveCrdExclVolArb)
	EVT_MENU   (IDC_CREATE_EXCL_VOL_MOL, AtomParamsDlgWX::OnCreateExclVolMol)
	EVT_MENU   (IDC_EDTAT_AMBER_FF_SYMB, AtomParamsDlgWX::OnEdtAtAmberFFSymb)
	EVT_MENU   (IDC_EDTAT_CLEAR_ATOM_FF_PARAMS, AtomParamsDlgWX::OnClearAtomFFParams)
	EVT_MENU   (IDC_EDTAT_ELEM_MASS, AtomParamsDlgWX::OnEdtAtElemMass)
	EVT_MENU   (IDC_EDTAT_SET_ATOM_ELEM_FROM_TEMPL, AtomParamsDlgWX::OnSetAtElemFromTempl)
	EVT_BUTTON (IDC_EDTAT_COPY_TO_CHMAP, AtomParamsDlgWX::OnCopyToChmap)
	EVT_BUTTON (IDC_EDTAT_SET_FROM_CHMAP, AtomParamsDlgWX::OnSetFromChmap)
	EVT_MENU   (IDC_EDTAT_PH_CH, AtomParamsDlgWX::OnEditAtPhCh)
	EVT_MENU   (IDC_EDTAT_FORM_CH_TEMPL, AtomParamsDlgWX::OnEditAtFormChTempl)
	EVT_MENU   (IDC_MIN_DIST, AtomParamsDlgWX::OnMinDist)
	EVT_MENU   (IDC_SET_SPEC_RAD, AtomParamsDlgWX::OnSetSpecRad)
	EVT_CHECKBOX (IDC_EDTAT_CHARGE, AtomParamsDlgWX::OnChangeProp)
	EVT_CHECKBOX (IDC_EDTAT_SOLV_ACCESS_AREA, AtomParamsDlgWX::OnChangeProp)
	EVT_CHECKBOX (IDC_EDTAT_COORD, AtomParamsDlgWX::OnChangeProp)
	EVT_CHECKBOX (IDC_EDTAT_ELEMENT, AtomParamsDlgWX::OnChangeProp) 
	EVT_CHECKBOX (IDC_EDTAT_FF_SYMB, AtomParamsDlgWX::OnChangeProp)
	EVT_CHECKBOX (IDC_EDTAT_MASS, AtomParamsDlgWX::OnChangeProp)
	EVT_CHECKBOX (IDC_EDTAT_NAME, AtomParamsDlgWX::OnChangeProp)
	EVT_CHECKBOX (IDC_EDTAT_RADIUS, AtomParamsDlgWX::OnChangeProp)
	EVT_CHECKBOX (IDC_EDTAT_VDWRAD, AtomParamsDlgWX::OnChangeProp)
	EVT_CHECKBOX (IDC_EDTAT_VDWENE, AtomParamsDlgWX::OnChangeProp)
	EVT_CHECKBOX (IDC_EDTAT_TEMP, AtomParamsDlgWX::OnChangeProp)
	EVT_CHECKBOX (IDC_EDTAT_HYBRID, AtomParamsDlgWX::OnChangeProp)
	EVT_CHECKBOX (IDC_EDTAT_HB_STATUS, AtomParamsDlgWX::OnChangeProp)
	EVT_CHECKBOX (IDC_CHECK_PER_BOX, AtomParamsDlgWX::OnCheckPerBox)
	EVT_MENU     (IDM_LOAD_COORDS, AtomParamsDlgWX::OnLoadCoords)
	EVT_MENU     (IDM_SET_STD_ZMAT,  AtomParamsDlgWX::OnSetStdZMat)
	EVT_MENU     (IDM_PRINT_ZMAT,  AtomParamsDlgWX::OnPrintZMat)
	EVT_MENU     (IDM_SET_ZMAT_INT_CRD_FROM_ATCRD, AtomParamsDlgWX::OnSetZMatIntCrdFromAtCrd)
	EVT_MENU     (IDM_SET_ATCRD_FROM_ZMAT, AtomParamsDlgWX::OnSetAtCrdFromZMat)
	EVT_MENU     (IDC_CALC_DIPOLE,   AtomParamsDlgWX::OnCalcDipole)
	EVT_MENU     (IDC_SET_VDW_RAD,   AtomParamsDlgWX::OnSetVdwRad)
	EVT_MENU     (IDC_SET_PARSE_RAD, AtomParamsDlgWX::OnSetParseRad)
	EVT_MENU     (IDC_SET_HPP_RAD,   AtomParamsDlgWX::OnSetHPPRad)
	EVT_BUTTON   (IDC_CLOSE, AtomParamsDlgWX::OnCloseBtn)
	EVT_CLOSE    (AtomParamsDlgWX::OnClose)
	EVT_MENU     (IDC_EDTAT_BACKBONE_CH, AtomParamsDlgWX::OnSetStdCharges) // changed from std amber charges
	EVT_MENU     (IDC_EDTAT_ZERO_CHARGES, AtomParamsDlgWX::OnSetZeroCharges)
	EVT_BUTTON   (IDC_EDTAT_UPDATE_ATLIST, AtomParamsDlgWX::OnUpdateAtomGroup)
	EVT_BUTTON   (IDC_EDTAT_TRANSFER_TO_WIN,   AtomParamsDlgWX::OnTransferToWin)
	EVT_BUTTON   (IDC_EDTAT_TRANSFER_FROM_WIN, AtomParamsDlgWX::OnTransferFromWin)
	EVT_GRID_CELL_CHANGED( AtomParamsDlgWX::OnEndLabelEdit )
END_EVENT_TABLE()


void AtomParamsDlgWX::OnEndLabelEdit(wxGridEvent& event)
{
	PrintLog(" Enter AtomParamsDlgWX::OnEndLabelEdit() \n");
	
	int icol = event.GetCol();
	int irow = event.GetRow();
	
	wxString str_val = m_atom_lctrl->GetCellValue(irow,icol);

	long ival;
	bool bres;
	double dval;
	char buf[256];
	char buf2[256];
	HaAtom* aptr;

	if(pmset == NULL) return;
    
	wxString at_id = m_atom_lctrl->GetRowLabelValue(irow);
//	aptr = pmset->GetAtomByRef(at_id.c_str());
	aptr = (HaAtom*) at_ptrs[irow];
	if(aptr == NULL)
	{
		sprintf(buf2," Can't find atom with Ref %s",at_id.ToStdString().c_str());
		PrintMessage(buf2);
		return;
	}

	if(icol == n_atname)
	{
		if(str_val.IsEmpty())
		{
			sprintf(buf," Invalid Input for atom name: %s ", str_val.ToStdString().c_str());
			PrintMessage(buf);
			sprintf(buf,"%s",aptr->GetName());
			m_atom_lctrl->SetCellValue(irow, icol, buf);
		}
		else
		{
			aptr->SetName( str_val.ToStdString() );
			aptr->FillRef(buf);
			m_atom_lctrl->SetRowLabelValue(irow, buf); 
		}
	}

	if(icol == n_hybrid)
	{
		aptr->SetHybrid( str_val.ToStdString() );
	}

	if(icol == n_hb_status)
	{
		aptr->SetHBStatus(str_val.c_str());
	}


	if(icol == n_ff_symbol)
	{
		if(str_val.IsEmpty())
		{
			sprintf(buf," Invalid Input for atom Force Field Symbol: %s ", str_val.ToStdString().c_str());
			PrintMessage(buf);
			sprintf(buf,"%s",aptr->GetFFSymbol().c_str());
			m_atom_lctrl->SetCellValue(irow, icol, buf);
		}
		else
		{
			aptr->SetFFSymbol( str_val.ToStdString() );
		}
	}

	if(icol == n_charge)
	{
		bres = str_val.ToDouble(&dval);

		if( !bres)
		{
			sprintf(buf," Invalid Input for atom charge: %s ", str_val.ToStdString().c_str());
			PrintMessage(buf);
			sprintf(buf,"%9.4f",aptr->GetCharge());
			m_atom_lctrl->SetCellValue(irow, icol, buf);
		}
		else
		{
			aptr->SetCharge(dval);
		}
	}
	if(icol == n_mass)
	{
		bres = str_val.ToDouble(&dval);
		
		if(!bres)
		{
			sprintf(buf," Invalid Input for atom mass: %s ", str_val.ToStdString().c_str());
			PrintMessage(buf);
			sprintf(buf,"%9.4f",aptr->GetMass() );
			m_atom_lctrl->SetCellValue(irow, icol, buf);
		}
		else
		{
			aptr->SetMass(dval);
		}
	}
	else if( icol == n_radius)
	{
		bres = str_val.ToDouble(&dval);
		
		if(!bres)
		{
			sprintf(buf," Invalid Input for atom radius: %s ", str_val.ToStdString().c_str());
			PrintMessage(buf);
			sprintf(buf,"%9.4f",aptr->radius);
			m_atom_lctrl->SetCellValue(irow, icol, buf);
		}
		else
		{
			aptr->radius = dval;
		}
	}
	else if( icol == n_temp)
	{
		bres = str_val.ToDouble(&dval);
		
		if(!bres)
		{
			sprintf(buf," Invalid Input for atom temperature field: %s ", str_val.ToStdString().c_str());
			PrintMessage(buf);
			sprintf(buf,"%9.4f",aptr->tempf);
			m_atom_lctrl->SetCellValue(irow, icol, buf);
		}
		else
		{
			aptr->tempf = dval;
		}
	}
	else if( icol == n_element)
	{
		bres = str_val.ToLong(&ival);
		
		if(!bres)
		{
			sprintf(buf," Invalid Input for atom element: %s ", str_val.ToStdString().c_str());
			PrintMessage(buf);
			sprintf(buf,"%d",aptr->GetElemNo() );
			m_atom_lctrl->SetCellValue(irow, icol, buf);
		}
		else
		{
			aptr->SetElemNo(ival);
		}
	}
	else if( icol == n_x_coord)
	{
		bres = str_val.ToDouble(&dval);
		
		if(!bres)
		{
			sprintf(buf," Invalid Input for X coordinate: %s ", str_val.ToStdString().c_str());
			PrintMessage(buf);
			sprintf(buf,"%9.4f",aptr->GetX_Ang() );
			m_atom_lctrl->SetCellValue(irow, icol, buf);
		}
		else
		{
			aptr->SetX_Ang( dval );
		}
	}
	else if( icol == n_y_coord)
	{
		bres = str_val.ToDouble(&dval);

		if(!bres)
		{
			sprintf(buf," Invalid Input for Y coordinate: %s ", str_val.ToStdString().c_str());
			PrintMessage(buf);
			sprintf(buf,"%9.4f", aptr->GetY_Ang());
			m_atom_lctrl->SetCellValue(irow, icol, buf);
		}
		else
		{
			aptr->SetY_Ang( dval );
		}
	}
	else if( icol == n_z_coord)
	{
		bres = str_val.ToDouble(&dval);

		if(!bres)
		{
			sprintf(buf," Invalid Input for Z coordinate: %s ", str_val.ToStdString().c_str());
			PrintMessage(buf);
			sprintf(buf,"%9.4f",aptr->GetZ_Ang() );
			m_atom_lctrl->SetCellValue(irow, icol, buf);
		}
		else
		{
			aptr->SetZ_Ang( dval );
		}
	}
}

void AtomParamsDlgWX::OnCalcDipole(wxCommandEvent& event)
{
	if(pmset != NULL)
	{
		pmset->CalcDipole();
	}
}

void AtomParamsDlgWX::OnSetVdwRad(wxCommandEvent& event)
{
	if(pmset != NULL)
	{
		pmset->SetVdwRadii();
		FillAtomGroup();
	}
}

void AtomParamsDlgWX::OnSetParseRad(wxCommandEvent& event)
{
	if(pmset != NULL)
	{
		pmset->SetParseRadii();
		FillAtomGroup();
	}
}

void AtomParamsDlgWX::OnSetHPPRad(wxCommandEvent& event)
{
	if(pmset != NULL)
	{
		pmset->SetHPPRadii();
		FillAtomGroup();
	}
}

void AtomParamsDlgWX::OnSetSpecRad(wxCommandEvent& event) 
{
	
	wxTextCtrl* rad_edit = (wxTextCtrl*) FindWindow(IDC_EDTAT_RAD);
	wxString str = rad_edit->GetValue();
	double rad;
	str.ToDouble(&rad);
	if(pmset != NULL)
	{
		AtomIteratorMolSet aitr(pmset);
		HaAtom* aptr = aitr.GetFirstAtom();
		for(; aptr; aptr = aitr.GetNextAtom())
		{
			if(aptr->Selected())
			{
				aptr->radius = rad;
			}
		}
		FillAtomGroup();
	}
}

void AtomParamsDlgWX::OnMinDist(wxCommandEvent& event) 
{
	HaAtom* aptr1 = NULL;
	HaAtom* aptr2 = NULL;
	double dist_min;
	if(pmset != NULL)
	{
		dist_min = pmset->FindClosestContact(aptr1,aptr2);
	}	
}


void AtomParamsDlgWX::OnSetStdCharges(wxCommandEvent& event)
{
	if(pmset == NULL) return;
	p_mol_editor->SetStdAtomicParams(pmset,BACKBONE_CHRG | PROT_CHARGED_GROUPS_CHRG);
    FillAtomGroup();
}

void AtomParamsDlgWX::OnEditAtPhCh(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	if(pmset == NULL) return;
	p_prot_rdx_mod->SetChargesForCurrentPH();
	FillAtomGroup();
}

void AtomParamsDlgWX::OnEditAtFormChTempl(wxCommandEvent& event)
{
	TransferDataFromWindow();
	if(pmset == NULL) return;
	p_mol_editor->SetFormalAtChrgFromTempl(pmset);
	FillAtomGroup();
}

void AtomParamsDlgWX::OnSetZeroCharges(wxCommandEvent& event)
{
	if(pmset == NULL) return;
	p_mol_editor->SetStdAtomicParams(pmset,ZERO_CHRG);
	FillAtomGroup();
}

void AtomParamsDlgWX::OnEdtAtAmberCh(wxCommandEvent& event) 
{
	if(pmset == NULL) return;
	p_mol_editor->SetStdAtomicParams(pmset,AMBER_ALL_ATOM_CHRGS);	
	FillAtomGroup();
}

void AtomParamsDlgWX::OnEdtAtHBondStatus(wxCommandEvent& event) 
{
	if(pmset == NULL) return;
	p_mol_editor->SetStdAtomicParams(pmset,ATOM_HBOND_DA_STATUS);	
	FillAtomGroup();
}

void AtomParamsDlgWX::OnCalcSolvAccessArea(wxCommandEvent& event) 
{
	if(pmset == NULL) return;
	pmset->CalcMolSurface();
	FillAtomGroup();
}

void AtomParamsDlgWX::OnCalcSolvAccessAreaGeoball(wxCommandEvent& event) 
{
	if(pmset == NULL) return;
	pmset->CalcSolventAccessArea();
	FillAtomGroup();
}

void AtomParamsDlgWX::OnSaveCrdExclVolArb(wxCommandEvent& event) 
{
	if(pmset == NULL) return;
	pmset->SaveCrdExclVolArb();
	FillAtomGroup();
}

void AtomParamsDlgWX::OnCreateExclVolMol(wxCommandEvent& event) 
{
	if(pmset == NULL) return;
	pmset->CreateExcludedVolumeMol();
	FillAtomGroup();
}



void AtomParamsDlgWX::OnEdtAtAmberFFSymb(wxCommandEvent& event) 
{
	p_mol_editor->SetStdAtomicParams(pmset,AMBER_ALL_ATOM_FF_SYMBOLS);	
	FillAtomGroup();	
}

void  AtomParamsDlgWX::OnClearAtomFFParams(wxCommandEvent& event)
{
	MolEditor::ClearAtomFFParams(pmset);
	FillAtomGroup();
}


void AtomParamsDlgWX::OnUpdateAtomGroup(wxCommandEvent& event)
{
	FillAtomGroup();
}

void AtomParamsDlgWX::OnTransferToWin(wxCommandEvent& event)
{
	TransferDataToWindow();
}

void AtomParamsDlgWX::OnTransferFromWin(wxCommandEvent& event)
{
	TransferDataFromWindow();
	pmset->RefreshAllViews(RFRefresh);
}

void
AtomParamsDlgWX::OnClose(wxCloseEvent& event)
{
	AtomParamsDlgWX::dlg_open = false;
	event.Skip();
}

void
AtomParamsDlgWX::OnCloseBtn(wxCommandEvent& event)
{
	Close();
}

void AtomParamsDlgWX::OnChangeProp(wxCommandEvent& event)
{
	TransferDataFromWindow();
	SetColumns();
	FillAtomGroup();
}

void AtomParamsDlgWX::OnCheckPerBox(wxCommandEvent& event)
{
	wxCheckBox* check_pbox = (wxCheckBox*) FindWindow(IDC_CHECK_PER_BOX);

	if(check_pbox->GetValue())
	{
		if( pmset->per_bc->IsValid() )
		{
			pmset->per_bc->Set(true);
		}
		else
		{
			pmset->per_bc->SetStdBox(pmset);
		}
	}
	else
	{
		pmset->per_bc->Set(false);
	}
	pmset->OnChangePeriodicity();
	TransferDataToWindow();
}

void AtomParamsDlgWX::OnEdtAtElemMass(wxCommandEvent& event) 
{
	if(pmset == NULL) return;
	p_mol_editor->SetStdAtomicParams(pmset,ATOM_MASSES_ELEMENT);
	FillAtomGroup();
}

void AtomParamsDlgWX::OnSetAtElemFromTempl(wxCommandEvent& event) 
{
	if(pmset == NULL) return;
	p_mol_editor->SetStdAtomicParams(pmset,ATOM_ELEM_FROM_TEMPL);
	FillAtomGroup();
	HaMolView* pview = pmset->GetActiveMolView();
	if(pview)pview->CPKColourAttrib();
	pmset->RefreshAllViews( RFRefresh | RFColour );
}

void AtomParamsDlgWX::OnLoadCoords(wxCommandEvent& event) 
{
    std::string file_ext;
    std::string fname = wxFileSelector("Select file with atom coordinates to load").ToStdString();

    if ( fname.empty() ) return;

	std::filesystem::path p(fname);
	std::string ext = p.extension().string();
	if( boost::starts_with(ext,".")) ext = ext.substr(1);  // remove leading "." from extension
	boost::to_upper(ext);

	int format = FormatXYZ;
	if( file_ext == "PDB" || file_ext == "ENT") format = FormatPDB;
	if( file_ext == "HLM" ) format = FormatHarlem;
	if( file_ext == "HIN" ) format = FormatHIN;
	pmset->SetCoordFromFile(fname, format);
}

void AtomParamsDlgWX::OnSetStdZMat(wxCommandEvent& event)
{
	ZMatCrd* pzm = pmset->GetZMat();
	pzm->InitStdZMat();
}

void AtomParamsDlgWX::OnPrintZMat(wxCommandEvent& event)
{
	ZMatCrd* pzm = pmset->GetZMat();
	if( pzm->IsEmpty() )
	{
		PrintLog(" Z-matrix is not set \n");
	}
	else
	{
		std::string str = pzm->SaveToString();
		PrintLog(" Z-matrix: \n%s\n",str.c_str());
	}
}

void AtomParamsDlgWX::OnSetZMatIntCrdFromAtCrd(wxCommandEvent& event)
{
	ZMatCrd* pzm = pmset->GetZMat();
	if( pzm->IsEmpty() )
	{
		PrintLog(" Error in AtomParamsDlgWX::OnSetZMatIntCrdFromAtCrd() \n");
		PrintLog(" Z-matrix is not set \n");
	}
	else
	{
		int ires = pzm->SetFromAtomCrd();
		if( ires )
		{
			PrintLog(" Z-Matrix internal coordinates were set from current cartesian atom coordinates"); 
		}
	}
}

void AtomParamsDlgWX::OnSetAtCrdFromZMat(wxCommandEvent& event)
{
	ZMatCrd* pzm = pmset->GetZMat();
	if( pzm->IsEmpty() )
	{
		PrintLog(" Error in AtomParamsDlgWX::OnSetAtCrdFromZMat() \n");
		PrintLog(" Z-matrix is not set \n");
	}
	else
	{
		int ires = pzm->SetAtomCrd();
		pmset->RefreshAllViews(RFApply | RFRefresh);
		if( ires )
		{
			PrintLog(" Atom Cartesian Coordinates were set from Z-Matrix internal coordinates ");  	 
		}
	}
}

//////////////////////////////////////////////////////////////////////
// ResidueParamsDlgWX dialog

int   ResidueParamsDlgWX::dlg_open = FALSE;

ResidueParamsDlgWX::ResidueParamsDlgWX(MolSet* new_pmset, wxWindow* parent):
wxFrame( parent, -1, "Edit Residue Parameters")
{
	pmset= new_pmset;
	p_mol_editor = pmset->GetMolEditor(true);

	n_res_name = -1;
	n_res_name_modifier = -1;
	n_res_num = -1;
	n_chain_name = -1;
	n_mutation = -1;
	
	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

	res_params_menu_bar = res_params_menu();
    this->SetMenuBar(res_params_menu_bar);    

	res_params_dlg(this,TRUE);
	residue_grid = (wxGrid*) FindWindow( IDC_RESIDUE_LIST );
	OnInitDialog();     
}

ResidueParamsDlgWX::~ResidueParamsDlgWX()
{
    dlg_open = FALSE;
}


void ResidueParamsDlgWX::OnInitDialog()
{
	wxCheckBox* check_ctrl;

	res_name_check = (wxCheckBox*) FindWindow(IDC_EDTRES_NAME);
	res_name_check->Bind(wxEVT_CHECKBOX, &ResidueParamsDlgWX::OnChangeProp,this);
	res_name_check->SetValue(true);

	res_name_modifier_check = (wxCheckBox*) FindWindow(IDC_EDTRES_MODIFIER);
	res_name_modifier_check->Bind(wxEVT_CHECKBOX, &ResidueParamsDlgWX::OnChangeProp, this);
	res_name_modifier_check->SetValue(true);

	res_num_check = (wxCheckBox*) FindWindow(IDC_EDTRES_RES_NUM);
	res_num_check->Bind(wxEVT_CHECKBOX, &ResidueParamsDlgWX::OnChangeProp, this);
	res_num_check->SetValue(true);

	chain_check = (wxCheckBox*) FindWindow(IDC_EDTRES_CHAIN);
	chain_check->Bind(wxEVT_CHECKBOX, &ResidueParamsDlgWX::OnChangeProp, this);
	chain_check->SetValue(true);

	mutation_check = (wxCheckBox*)FindWindow(IDC_EDTRES_MUT_RES);
	mutation_check->Bind(wxEVT_CHECKBOX, &ResidueParamsDlgWX::OnChangeProp, this);
	mutation_check->SetValue(true);

	residue_grid->Bind(wxEVT_GRID_SELECT_CELL, &ResidueParamsDlgWX::OnSelectResidueRow, this);

	p_renum_start_n = (wxTextCtrl*) FindWindow(IDC_RESPAR_RENUM_START_NUM);
	p_renum_start_n->SetValue("1");

	p_sel_res_text = (wxTextCtrl*)FindWindow(IDC_RESPAR_SELECTED_RES);
	p_sel_res_text->SetValue("");

	p_mut_res_type_text = (wxTextCtrl*)FindWindow(IDC_RESPAR_MUTATED_RES_TYPE);
	p_mut_res_type_text->SetValue("");

	p_mutate_res_btn = (wxButton*)FindWindow(IDC_RESPAR_MUTATE_RES);
	p_mutate_res_btn->Bind(wxEVT_BUTTON, &ResidueParamsDlgWX::OnMutateResidue, this);

	p_set_transform_res_btn = (wxButton*)FindWindow(IDC_RESPAR_SET_TRANSFORMATION);
	p_set_transform_res_btn->Bind(wxEVT_BUTTON, &ResidueParamsDlgWX::OnSetTransformation, this);
	
	SetColumns();	
	TransferDataToWindow();
	dlg_open = TRUE;
}

void ResidueParamsDlgWX::SetColumns()
{
	n_res_name = -1;
	n_res_name_modifier = -1;
	n_chain_name = -1;
	n_res_num = -1;
	n_mutation = -1;
	
	int num_cols = 0;

    if(res_name_check->IsChecked() ) num_cols += 1;
    if(res_name_modifier_check->IsChecked()) num_cols += 1;
	if(res_num_check->IsChecked()) num_cols += 1;
    if(chain_check->IsChecked() ) num_cols += 1; 
	if(mutation_check->IsChecked()) num_cols += 1;
    	
	int tot_width = 550;
	int height = 400;
	residue_grid->GetClientSize(&tot_width, &height);

	if(num_cols == 0) num_cols = 1;
	int col_width = (tot_width*2 /3)/num_cols;

	int icol = -1;
	residue_grid->SetColLabelSize(21);
    
	int ncol_act = residue_grid->GetNumberCols();
	if(ncol_act != num_cols)
	{
		residue_grid->DeleteCols(0,ncol_act);
		residue_grid->InsertCols(0,num_cols);
	}

	if(res_name_check->IsChecked())
	{
		icol++;
		residue_grid->SetColLabelValue(icol, "Name");
		n_res_name = icol;
	}

	if(res_name_modifier_check->IsChecked())
	{
		icol++;
		residue_grid->SetColLabelValue(icol, "Modifier");
		n_res_name_modifier = icol;
	}

	if(res_num_check->IsChecked())
	{
		icol++;
		residue_grid->SetColLabelValue(icol, "Num");
		n_res_num = icol;
	}

	if(chain_check->IsChecked() )
	{
		icol++;
		residue_grid->SetColLabelValue(icol, "Ch");
		n_chain_name = icol;
	}

	if (mutation_check->IsChecked())
	{
		icol++;
		residue_grid->SetColLabelValue(icol, "Mut");
		n_mutation = icol;
	}
}

void ResidueParamsDlgWX::FillResidueGrid()
{
	residue_grid->ClearGrid();

	if(pmset == NULL) return;
	
	int nres_sel = 0;

	MolSet::ResidueIterator ritr(pmset);
	for(HaResidue* pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes())
	{
		if(!pres->HasSelectedAtoms()) continue;
		nres_sel++;
	}

	const int ncol = residue_grid->GetNumberCols();
	const int nrow = residue_grid->GetNumberRows();

	const int nv = res_ptrs.size();

	if(nrow != nres_sel || nv != nres_sel)
	{
		if( nrow > 0) residue_grid->DeleteRows(0, nrow);
        residue_grid->AppendRows(nres_sel);
		res_ptrs.resize(nres_sel);
	}

	int itm = -1;
	int max_lbl = 0;

	for(HaResidue* pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes())
	{
		if(!pres->HasSelectedAtoms()) continue;
		itm++;

		std::string rref = pres->GetRef();
		res_ptrs[itm] = (void*) pres;

		if(rref.size() > max_lbl) max_lbl = rref.size();
	
		residue_grid->SetRowLabelValue(itm, rref); 

		if( n_res_name >= 0 )
		{
			std::string ss = pres->GetName();
			residue_grid->SetCellValue(itm, n_res_name, ss); 
		}

		if( n_res_name_modifier >= 0 )
		{
			std::string ss = pres->NameModifier;
			residue_grid->SetCellValue(itm, n_res_name_modifier, ss); 
		}

		if( n_res_num >= 0 )
		{
			std::string ss = (boost::format("%d") % pres->GetSerNo()).str();
			residue_grid->SetCellValue(itm, n_res_num, ss); 
		}

		if( n_chain_name >= 0 )
		{
			HaChain* chain = pres->GetHostChain();
			std::string ss = std::string(1, chain->ident);
			residue_grid->SetCellValue(itm, n_chain_name, ss); 
		}

		if (n_mutation >= 0)
		{
			std::string ss;
			if (pres->IsAlchemicalTransformationSet())
			{
				ss = pres->p_res_transform->res_name_b;
			}
			residue_grid->SetCellValue(itm, n_mutation, ss);
		}
	}
	residue_grid->SetRowLabelSize((int)(max_lbl*8.5));
	residue_grid->AutoSizeColumns();
}

bool ResidueParamsDlgWX::TransferDataToWindow()
{
	FillResidueGrid();

	pres_sel = this->GetSelectedResidue();

	if (pres_sel)
	{
		this->p_sel_res_text->SetValue(pres_sel->GetRef());
	}
	else
	{
		this->p_sel_res_text->SetValue("");
	}

	return wxFrame::TransferDataToWindow();
}

BEGIN_EVENT_TABLE(ResidueParamsDlgWX, wxFrame)
	EVT_BUTTON(IDC_EDTRES_UPDATE_RESLIST, ResidueParamsDlgWX::OnUpdateResidueList)
	EVT_BUTTON(IDC_RESPAR_RENUM, ResidueParamsDlgWX::OnResidueRenumber)
	EVT_BUTTON(IDC_CLOSE, ResidueParamsDlgWX::OnCloseBtn)
	EVT_CLOSE (ResidueParamsDlgWX::OnClose)
	EVT_MENU(IDC_RESPAR_CHECK_STRUCT, ResidueParamsDlgWX::OnCheckStruct)
	EVT_MENU(IDC_RESPAR_DEL_EXTRA_AT, ResidueParamsDlgWX::OnDelExtraAt)
	EVT_MENU(IDC_RESPAR_ADD_MISS_AT, ResidueParamsDlgWX::OnAddMissAt)
	EVT_MENU(IDC_FIX_BONDS_USING_TEMPL, ResidueParamsDlgWX::OnFixBondsUsingTempl)
	EVT_MENU(IDC_ORDER_ATOMS_RES, ResidueParamsDlgWX::OnOrderAtomsInRes)
	EVT_MENU(IDC_RENAME_ATOMS_TO_AMBER,   ResidueParamsDlgWX::OnRenameAtomsToAmber)
	EVT_MENU(IDC_RENAME_ATOMS_TO_GROMACS, ResidueParamsDlgWX::OnRenameAtomsToGromacs)
	EVT_MENU(IDC_CONVERT_WAT_ARROW_VB, ResidueParamsDlgWX::OnConvertWaterArrowVB)
	EVT_MENU(IDC_CONVERT_WAT_FAST_AMBER, ResidueParamsDlgWX::OnConvertWaterFastAmber)
	EVT_GRID_CELL_CHANGED( ResidueParamsDlgWX::OnEndLabelEdit)
END_EVENT_TABLE()

void ResidueParamsDlgWX::OnUpdateResidueList(wxCommandEvent& event)
{
	FillResidueGrid();
}

HaResidue* ResidueParamsDlgWX::GetSelectedResidue()
{
	int nv = res_ptrs.size();
	int irow_selected = -1;

	wxArrayInt selectedRows = residue_grid->GetSelectedRows();
	for (int i = 0; i < selectedRows.size(); ++i) {
		irow_selected = selectedRows[i];
		break;
	}

	if (irow_selected >= 0)
	{
		HaResidue* pres_sel = (HaResidue*)res_ptrs[irow_selected];
		return pres_sel;
	}
	return nullptr;
}

void ResidueParamsDlgWX::OnResidueRenumber(wxCommandEvent& event)
{
	wxTextCtrl* p_renum_start_n = (wxTextCtrl*) FindWindow(IDC_RESPAR_RENUM_START_NUM);
	wxString n_start_str = p_renum_start_n->GetValue();
	long n_start;
	if( n_start_str.ToLong(&n_start) )
	{
		pmset->RenumberSelectedRes(n_start);
	}
	else
	{
		PrintLog("Invalid Starting Number for Residue Renumbering\n");
		p_renum_start_n->SetValue("1");
	}
	FillResidueGrid();
}

void ResidueParamsDlgWX::OnMutateResidue(wxCommandEvent& event)
{
	HaResidue* p_res_sel = this->GetSelectedResidue();
	if (!pres_sel)
	{
		PrintLog("No Residue Selected \n");
	}
	std::string mut_res_type = p_mut_res_type_text->GetValue();
	if (!mut_res_type.empty())
	{
		PrintLog("Mutate Residue %s to type: %s \n", pres_sel->GetRef(), mut_res_type.c_str() );
		p_res_sel->MutateTo(mut_res_type);
	}
	TransferDataToWindow();
}

void ResidueParamsDlgWX::OnSetTransformation(wxCommandEvent& event)
{
	PrintLog("ResidueParamsDlgWX::OnSetTransformation() \n");
	pres_sel = this->GetSelectedResidue();
	if (!pres_sel)
	{
		PrintLog("No Residue Selected \n");
		return;
	}
	std::string mut_res_type = p_mut_res_type_text->GetValue();
	if (!mut_res_type.empty())
	{
		PrintLog("Set Alchemical Transformation for Residue %s to type: %s \n", pres_sel->GetRef(), mut_res_type);
		pres_sel->SetAlchemicalTransformation(mut_res_type);
	}
	TransferDataToWindow();
}


void ResidueParamsDlgWX::OnCloseBtn(wxCommandEvent& event)
{
	Close();
}

void ResidueParamsDlgWX::OnClose(wxCloseEvent& event)
{
	ResidueParamsDlgWX::dlg_open = FALSE;
	event.Skip();

}

void ResidueParamsDlgWX::OnChangeProp(wxCommandEvent& event)
{
	TransferDataFromWindow();
	SetColumns();
	FillResidueGrid();
}

void ResidueParamsDlgWX::OnEndLabelEdit(wxGridEvent& event)
{
	PrintLog( " Enter ResidueParamsDlgWX::OnEndLabelEdit() \n");
	
	int icol = event.GetCol();
	int irow = event.GetRow();
	
    wxString str_val = residue_grid->GetCellValue(irow,icol);

	str_val.Strip(wxString::both);
	
	char buf[256];
	char buf2[256];

	if(pmset == NULL) return;
    
	HaResidue* pres = (HaResidue*) res_ptrs[irow];

	if(icol == n_res_name)
	{
		if(str_val.IsEmpty())
		{
			sprintf(buf," Invalid Input for residue name: %s ", str_val.ToStdString().c_str());
			PrintMessage(buf);
			sprintf(buf,"%s",pres->GetName());
			residue_grid->SetCellValue(irow, icol, buf);
		}
		else
		{
			pres->SetName( str_val.ToStdString() );
			pres->FillRef(buf);
			residue_grid->SetRowLabelValue(irow, buf); 
		}
	}
	if(icol == n_res_name_modifier)
	{
		pres->SetNameModifier( str_val.ToStdString() );
	}

	if(icol == n_res_num)
	{
		if(str_val.IsEmpty())
		{
			sprintf(buf," Invalid Input for Residue Number: %s ", str_val.ToStdString().c_str());
			PrintMessage(buf);
			sprintf(buf,"%d", pres->GetSerNo());
			residue_grid->SetCellValue(irow, icol, buf);
		}
		else
		{
//			pres->Set(str.c_str());
		}
	}

	std::string stri;
	std::string str_ha = str_val.ToStdString();

	if(icol == n_chain_name)
	{
		
	}
}

void ResidueParamsDlgWX::OnSelectResidueRow(wxGridEvent& event)
{
	int irow_selected = event.GetRow();

	pres_sel = nullptr;
	if (irow_selected >= 0 && irow_selected < res_ptrs.size())
	{
		pres_sel = (HaResidue*)res_ptrs[irow_selected];
	}
	if (pres_sel)
	{
		p_sel_res_text->SetValue(pres_sel->GetRef());
	}
	else
	{
		p_sel_res_text->SetValue("");
	}
}

void ResidueParamsDlgWX::OnCheckStruct(wxCommandEvent& event) 
{
	MolSet::ResidueIterator ritr(pmset);
	HaResidue* pres;
	int all_valid = TRUE;
	for( pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes())
	{
		if(pres->HasSelectedAtoms())
		{
			int ires = pres->CheckStruct();
			if( !ires ) all_valid = FALSE;
		}
	}
	if( all_valid )
	{
		PrintLog("\n All residues are valid \n");
	}
	else
	{
		PrintLog("\n Not all residues have valid structures \n");
	}
}

void ResidueParamsDlgWX::OnDelExtraAt(wxCommandEvent& event) 
{
	p_mol_editor->DeleteExtraAtoms(pmset);	
}

void ResidueParamsDlgWX::OnAddMissAt(wxCommandEvent& event) 
{
	p_mol_editor->AddMissingAtoms(pmset);
}

void ResidueParamsDlgWX::OnFixBondsUsingTempl(wxCommandEvent& event)
{
	p_mol_editor->FixBondsUsingTempl(pmset);
}

void ResidueParamsDlgWX::OnOrderAtomsInRes(wxCommandEvent& event)
{
	p_mol_editor->OrderAtomsInRes(pmset);
}

void ResidueParamsDlgWX::OnRenameAtomsToAmber(wxCommandEvent& event)
{
	p_mol_editor->RenameAtomsToAmber(pmset);
}

void ResidueParamsDlgWX::OnRenameAtomsToGromacs(wxCommandEvent& event)
{
	p_mol_editor->RenameAtomsToGromacs(pmset);
}
void ResidueParamsDlgWX::OnConvertWaterArrowVB(wxCommandEvent& event)
{
	p_mol_editor->ConvertWaterArrowVB(pmset);
	pmset->RefreshAllViews(RFApply);
}
void ResidueParamsDlgWX::OnConvertWaterFastAmber(wxCommandEvent& event)
{
	p_mol_editor->ConvertWaterFastAmber(pmset);
	pmset->RefreshAllViews(RFApply);
}

/////////////////////////////////////////////////////////////////////
// EditGeomDlgWX Dialog:

int EditGeomDlgWX::dlg_open = FALSE;
EditGeomDlgWX* EditGeomDlgWX::active_dlg_ptr = NULL;

EditGeomDlgWX::EditGeomDlgWX(MolSet* new_pmset, wxWindow* parent):
wxFrame( parent, -1, "Edit Molecular Geometry")

{
	pmset= new_pmset;
	p_mol_editor = pmset->GetMolEditor(true);
    dlg_open = TRUE;
    EditGeomDlgWX::active_dlg_ptr = this;

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);
	
	edit_geom_dlg(this,TRUE);

	OnInitDialog();
}

EditGeomDlgWX::~EditGeomDlgWX()
{
    dlg_open = FALSE;
    EditGeomDlgWX::active_dlg_ptr = NULL;
}

void EditGeomDlgWX::OnInitDialog()
{
	aptr1 = NULL;
	aptr2 = NULL;
	aptr3 = NULL;
	aptr4 = NULL;

	edit_at1= (wxTextCtrl*) FindWindow(IDC_EDTORS_AT1);
	edit_at2= (wxTextCtrl*) FindWindow(IDC_EDTORS_AT2);
	edit_at3= (wxTextCtrl*) FindWindow(IDC_EDTORS_AT3);
	edit_at4= (wxTextCtrl*) FindWindow(IDC_EDTORS_AT4);

	edit_at1->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( EditGeomDlgWX::SetIntFocusAt ), NULL, this ); 
	edit_at2->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( EditGeomDlgWX::SetIntFocusAt ), NULL, this );
	edit_at3->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( EditGeomDlgWX::SetIntFocusAt ), NULL, this );
	edit_at4->Connect( wxEVT_SET_FOCUS, wxFocusEventHandler( EditGeomDlgWX::SetIntFocusAt ), NULL, this );

	active_edit_at = edit_at1;
}

BEGIN_EVENT_TABLE(EditGeomDlgWX,wxFrame)
    EVT_BUTTON(IDC_EDTG_CALC_GEOM, EditGeomDlgWX::OnCalcGeom)
	EVT_BUTTON(IDC_EDTG_SET_GEOM, EditGeomDlgWX::OnSetGeom)
	EVT_BUTTON(IDC_EDTG_INCREM_TORS, EditGeomDlgWX::OnIncremTors)
	EVT_BUTTON(IDC_EDTG_DEL_BOND, EditGeomDlgWX::OnDelBond)
	EVT_BUTTON(IDC_EDTG_CREATE_BOND, EditGeomDlgWX::OnCreateBond)
	EVT_BUTTON(IDC_EDTG_DEL_ATOM, EditGeomDlgWX::OnDelAtom)
	EVT_BUTTON(IDC_EDTG_MERGE_MOL, EditGeomDlgWX::OnMergeMol)
	EVT_BUTTON(IDC_EDTG_ATTACH_MOL, EditGeomDlgWX::OnAttachMol)
    EVT_BUTTON(IDU_ATOM_PICK, EditGeomDlgWX::OnAtomSelect)
	EVT_CLOSE(EditGeomDlgWX::OnClose)
	EVT_SET_FOCUS(EditGeomDlgWX::SetIntFocusAt)
END_EVENT_TABLE()

void EditGeomDlgWX::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
    EditGeomDlgWX::active_dlg_ptr = NULL;
	event.Skip();
}

void EditGeomDlgWX::SetIntFocusAt(wxFocusEvent& event)
{
//	PrintLog(" EditGeomDlgWX::SetIntFocusAt()  \n");
	wxTextCtrl* text_ctrl;

	if(event.GetId() == IDC_EDTORS_AT1)
		active_edit_at = edit_at1;
	else if( event.GetId() == IDC_EDTORS_AT2)
		active_edit_at = edit_at2;
	else if( event.GetId() == IDC_EDTORS_AT3)
		active_edit_at = edit_at3;
	else if( event.GetId() == IDC_EDTORS_AT4)
		active_edit_at = edit_at4;

	int i;
	wxColour col_inact(255,255,255);
	for(i=0; i < 4; i++)
	{
		if(i == 0) text_ctrl = edit_at1;
		if(i == 1) text_ctrl = edit_at2;
		if(i == 2) text_ctrl = edit_at3;
		if(i == 3) text_ctrl = edit_at4;

		wxColour col_sel(255,200,50);
		if( active_edit_at == text_ctrl )
		{
			text_ctrl->SetBackgroundColour(col_sel);
			text_ctrl->Refresh();
		}
		else
		{
			text_ctrl->SetBackgroundColour(col_inact);
			text_ctrl->Refresh();
		}
	}
	event.Skip();
}

void EditGeomDlgWX::OnAtomSelect(wxCommandEvent& event)
{
//	PrintLog("EditGeomDlgWX OnAtomSelect got IDU_ATOM_PICK message \n");
	char buf[256],buf2[256];
	if(PkAtom)
	{
		PkAtom->FillRef(buf2);
		MolSet* pmset_at = PkAtom->GetHostMolSet();
		if(pmset_at != this->pmset)
		{
			sprintf(buf,"!%s!%s",pmset_at->GetName(),buf2);
		}
		else
		{
			strcpy(buf,buf2);
		}

		if( active_edit_at == edit_at1 )
		{
			edit_at1->SetValue(buf);
			aptr1 = PkAtom;
			edit_at2->SetFocus();
		}
		else if( active_edit_at == edit_at2 )
		{
			edit_at2->SetValue(buf);
			aptr2 = PkAtom;
			edit_at3->SetFocus();
		}
		else if( active_edit_at == edit_at3 )
		{
			edit_at3->SetValue(buf);
			aptr3 = PkAtom;
			edit_at4->SetFocus();
		}
		else if( active_edit_at == edit_at4 )
		{
			edit_at4->SetValue(buf);
			aptr4 = PkAtom;
			edit_at1->SetFocus();
		}
	}
}

void EditGeomDlgWX::GetAtoms()
{
   char buf[256];
   std::string str;

   std::string mset_name;
   wxString str_wx;

   str_wx = edit_at1->GetValue();
   strcpy(buf,str_wx.c_str());
   if(buf[0] == '!') 
      aptr1 = pApp->GetAtomByRef(buf);
   else
	  aptr1 = pmset->GetAtomByRef(buf);
 
   str_wx = edit_at2->GetValue();
   strcpy(buf,str_wx.c_str());   
   if(buf[0] == '!') 
      aptr2 = pApp->GetAtomByRef(buf);
   else
	  aptr2 = pmset->GetAtomByRef(buf);
  
   str_wx = edit_at3->GetValue();
   strcpy(buf,str_wx.c_str());   
   str = buf;
   boost::trim(str);

   if(str.empty())
   {
	   aptr3 = NULL;
   }
   else
   {
	   if(buf[0] == '!')
	   {
		   aptr3 = pApp->GetAtomByRef(buf);
	   }
	   else
	   {
		   aptr3 = pmset->GetAtomByRef(buf);
	   }
   }
   
   str_wx = edit_at4->GetValue();
   strcpy(buf,str_wx.c_str());   

   str = buf;
   boost::trim(str);

   if(str.empty())
   {
	   aptr4 = NULL;
   }
   else
   {
	   if(buf[0] == '!')
	   {
		   aptr4 = pApp->GetAtomByRef(buf);
	   }
	   else
	   {
		   aptr4 = pmset->GetAtomByRef(buf);
	   }
   }
}

void EditGeomDlgWX::OnCalcGeom(wxCommandEvent& event) 
{
   char buf[256];

   wxTextCtrl* edit_dist12 = (wxTextCtrl*) FindWindow(IDC_EDTG_DIST12);
   wxTextCtrl* edit_dist23 = (wxTextCtrl*) FindWindow(IDC_EDTG_DIST23);
   wxTextCtrl* edit_dist34 = (wxTextCtrl*) FindWindow(IDC_EDTG_DIST34);

   wxTextCtrl* edit_ang123 = (wxTextCtrl*) FindWindow(IDC_EDTG_ANG123);
   wxTextCtrl* edit_ang234 = (wxTextCtrl*) FindWindow(IDC_EDTG_ANG234);

   wxTextCtrl* edit_tors_val  = (wxTextCtrl*) FindWindow(IDC_EDTG_TORS_VAL);
   wxTextCtrl* edit_tors_incr = (wxTextCtrl*) FindWindow(IDC_EDTG_TORS_INCR);

   edit_dist12->SetValue("");
   edit_dist23->SetValue("");
   edit_dist34->SetValue("");

   edit_ang123->SetValue("");
   edit_ang234->SetValue("");
   
   edit_tors_val->SetValue("");

   double dist,ang,tors;

   this->GetAtoms();

   if(aptr1 != NULL && aptr2 != NULL)
   {
	   dist = Vec3D::CalcDistance(aptr1,aptr2,ANGSTROM_U);
	   sprintf(buf,"%9.4f",dist);
	   edit_dist12->SetValue(buf);

	   if( aptr3 != NULL )
	   {
          ang = Vec3D::CalcAngle(aptr1,aptr2,aptr3)*RAD_TO_DEG;
		  sprintf(buf,"%9.4f",ang);
		  edit_ang123->SetValue(buf);

		  if( aptr4 != NULL )
		  {
              tors = Vec3D::CalcTorsion(aptr1,aptr2,aptr3,aptr4)*RAD_TO_DEG;
		      sprintf(buf,"%9.4f",tors);
		      edit_tors_val->SetValue(buf);
		  }
	   }
   }

   if(aptr2 != NULL && aptr3 != NULL)
   {
	   dist = Vec3D::CalcDistance(aptr2,aptr3,ANGSTROM_U);
	   sprintf(buf,"%9.4f",dist);
	   edit_dist23->SetValue(buf);

	   if( aptr4 != NULL )
	   {
          ang = Vec3D::CalcAngle(aptr2,aptr3,aptr4)*RAD_TO_DEG;
		  sprintf(buf,"%9.4f",ang);
		  edit_ang234->SetValue(buf);
	   }
   }

   if(aptr3 != NULL && aptr4 != NULL)
   {
	   dist = Vec3D::CalcDistance(aptr3,aptr4,ANGSTROM_U);
	   sprintf(buf,"%9.4f",dist);
	   edit_dist34->SetValue(buf);
   }
   
}

void EditGeomDlgWX::OnSetGeom(wxCommandEvent& event) 
{	
   char buf[256];

   wxTextCtrl* edit_dist12 = (wxTextCtrl*) FindWindow(IDC_EDTG_DIST12);
   wxTextCtrl* edit_dist23 = (wxTextCtrl*) FindWindow(IDC_EDTG_DIST23);
   wxTextCtrl* edit_dist34 = (wxTextCtrl*) FindWindow(IDC_EDTG_DIST34);

   wxTextCtrl* edit_ang123 = (wxTextCtrl*) FindWindow(IDC_EDTG_ANG123);
   wxTextCtrl* edit_ang234 = (wxTextCtrl*) FindWindow(IDC_EDTG_ANG234);

   wxTextCtrl* edit_tors_val  = (wxTextCtrl*) FindWindow(IDC_EDTG_TORS_VAL);
   wxTextCtrl* edit_tors_incr = (wxTextCtrl*) FindWindow(IDC_EDTG_TORS_INCR);

   int ires;
//   bool bres;
   double val;
   double vold;

   this->GetAtoms();

   wxString str;
   
   str = edit_dist12->GetValue();
   strcpy(buf,str.c_str());

   ires = sscanf(buf,"%lf",&val);
   if(ires != 0 && ires != EOF)
   {
	   if(aptr1 != NULL && aptr2 != NULL)
	   {
		   vold = Vec3D::CalcDistance(aptr1,aptr2,ANGSTROM_U);
		   if( fabs(val - vold) > 0.0001)
		   {
				p_mol_editor->SetBondDist(aptr1,aptr2,val);
		   }
	   }
   }

   str = edit_dist23->GetValue();
   strcpy(buf,str.c_str());

   ires = sscanf(buf,"%lf",&val);
   if(ires != 0 && ires != EOF)
   {
	   if(aptr2 != NULL && aptr3 != NULL)
	   {
		   vold = Vec3D::CalcDistance(aptr2,aptr3,ANGSTROM_U);
		   if( fabs(val - vold) > 0.0001)
		   {
				p_mol_editor->SetBondDist(aptr2,aptr3,val);
		   }
	   }
   }

   str = edit_dist34->GetValue();
   strcpy(buf,str.c_str());

   ires = sscanf(buf,"%lf",&val);
   if(ires != 0 && ires != EOF)
   {
	   if(aptr3 != NULL && aptr4 != NULL)
	   {
		   vold = Vec3D::CalcDistance(aptr3,aptr4,ANGSTROM_U);
		   if( fabs(val - vold) > 0.0001)
		   {
				p_mol_editor->SetBondDist(aptr3,aptr4,val);
		   }
	   }
   }

   str = edit_ang123->GetValue();
   strcpy(buf,str.c_str());

   ires = sscanf(buf,"%lf",&val);
   if(ires != 0 && ires != EOF)
   {
	   if(aptr1 != NULL && aptr2 != NULL && aptr3 != NULL)
	   {
		   vold = Vec3D::CalcAngle(aptr1,aptr2,aptr3)*RAD_TO_DEG;
		   if( fabs(val - vold) > 0.0001)
		   {
				p_mol_editor->SetAngle(aptr1,aptr2,aptr3,val*DEG_TO_RAD);
		   }
	   }
   }

   str = edit_ang234->GetValue();
   strcpy(buf,str.c_str());

   ires = sscanf(buf,"%lf",&val);
   if(ires != 0 && ires != EOF)
   {
	   if(aptr2 != NULL && aptr3 != NULL && aptr4 != NULL)
	   {
		   vold = Vec3D::CalcAngle(aptr2,aptr3,aptr4)*RAD_TO_DEG;
		   if( fabs(val - vold) > 0.0001)
		   {
				p_mol_editor->SetAngle(aptr2,aptr3,aptr4,val*DEG_TO_RAD);
		   }
	   }
   }

   str = edit_tors_val->GetValue();
   strcpy(buf,str.c_str());

   ires = sscanf(buf,"%lf",&val);
   if(ires != 0 && ires != EOF)
   {
	   if(aptr1 != NULL && aptr2 != NULL && aptr3 != NULL && aptr4 != NULL)
	   {
		   vold = Vec3D::CalcTorsion(aptr1,aptr2,aptr3,aptr4)*DEG_TO_RAD;
		   if( fabs(val - vold) > 0.0001)
		   {
				p_mol_editor->SetTorsion(aptr1,aptr2,aptr3,aptr4,val*DEG_TO_RAD);
		   }
	   }
   }

   OnCalcGeom(event);
   pmset->AnnounceGeomChange();
//   pmset->RefreshAllViews(RFRefresh);
}

void EditGeomDlgWX::OnIncremTors(wxCommandEvent& event) 
{
	
}

void EditGeomDlgWX::OnDelAtom(wxCommandEvent& event) 
{
	wxTextCtrl* edit_at1= (wxTextCtrl*) FindWindow(IDC_EDTORS_AT1);
	this->GetAtoms();
	if(aptr1 != NULL)
	{
		pmset->DeleteAtom(aptr1);
		edit_at1->SetValue("");
		aptr1 = NULL;
	}
	pmset->AnnounceGeomChange();
}

void EditGeomDlgWX::OnDelBond(wxCommandEvent& event) 
{
	this->GetAtoms();
	if(aptr1 != NULL && aptr2 != NULL )
	{ 
		HaAtom::DeleteBond(aptr1,aptr2);
		pmset->AnnounceGeomChange();
	}
}

void EditGeomDlgWX::OnCreateBond(wxCommandEvent& event) 
{
	this->GetAtoms();
	if(aptr1 != NULL && aptr2 != NULL )
	{
		HaAtom::CreateBond(aptr1,aptr2);
		pmset->AnnounceGeomChange();
	}
}

void EditGeomDlgWX::OnMergeMol(wxCommandEvent& event) 
{
	this->GetAtoms();
	if(aptr1 != NULL && aptr2 != NULL )
	{
		HaMolecule* pMol1 = aptr1->GetHostMol();
		HaMolecule* pMol2 = aptr2->GetHostMol();
		if(pMol1 != pMol2)
		{
			p_mol_editor->MergeMolecules(pMol1,pMol2);
			pmset->AnnounceGeomChange();
		}
		wxTextCtrl* edit_at1= (wxTextCtrl*) FindWindow(IDC_EDTORS_AT1);
		wxTextCtrl* edit_at2= (wxTextCtrl*) FindWindow(IDC_EDTORS_AT2);
		edit_at1->SetValue("");
		aptr1 = NULL;
        edit_at2->SetValue("");
		aptr2 = NULL;
	}
}

void EditGeomDlgWX::OnAttachMol(wxCommandEvent& event) 
{
	PrintLog(" EditGeomDlgWX::OnAttachMol() pt 1 \n");
	this->GetAtoms();
	if(aptr1 != NULL && aptr2 != NULL )
	{
		HaMolecule::AttachFragment(aptr1,aptr2);
		pmset->AnnounceGeomChange();
		
		edit_at1->SetValue("");
		aptr1 = NULL;
        edit_at2->SetValue("");
		aptr2 = NULL;
	}
}

/////////////////////////////////////////////////////////////////////
// EditGroupsDlg Dialog:

int EditGroupsDlg::dlg_open = FALSE;
EditGroupsDlg* EditGroupsDlg::active_dlg_ptr = NULL;

EditGroupsDlg::EditGroupsDlg(MolSet* new_pmset, int itype, wxWindow* parent):
wxFrame( parent, -1, "Edit Atom Groups")
{
	this->SetExtraStyle(wxWS_EX_VALIDATE_RECURSIVELY);
	pmset= new_pmset;
	igrp_type =  itype;
    dlg_open = TRUE;
    active_dlg_ptr = this;

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

    wxMenuBar* edit_groups_menu_bar = edit_groups_menu();
    SetMenuBar(edit_groups_menu_bar);    
	
	edit_groups_dlg(this,TRUE);

	modal_run_flag = FALSE;
	p_loc_event_loop = NULL;

	OnInitDialog();
}

EditGroupsDlg::~EditGroupsDlg()
{
    dlg_open = FALSE;
    active_dlg_ptr = NULL;
}

void EditGroupsDlg::OnInitDialog()
{
	wxRadioBox* pick_level_box = (wxRadioBox*) FindWindow(IDC_PICKING_LEVEL);

	pick_level_box->SetSelection(0);

	if(igrp_type > 2 || igrp_type < 0)
		igrp_type = 1;
	
	wxRadioBox* group_type_box = (wxRadioBox*) FindWindow(IDC_TYPE_GROUP);
	group_type_box->SetSelection(igrp_type);

//	for( i = 0; i < 3; i++)
//	{
//		wxRadioButton* btn = (wxRadioButton*) FindWindow( ids_type[i]);
//		if(igrp_type == i) 
//			btn->SetValue(true);
//		else
//			btn->SetValue(false);
//	}


	wxListBox* grp_list= (wxListBox*) FindWindow(IDC_EDTGRP_LIST_GRP);
	if( grp_list->GetCount() > 0)
	{
		grp_list->SetSelection(0);
		OnChangeSelGroup();
	}
	TransferDataToWindow();
}


bool EditGroupsDlg::TransferDataToWindow()
{
	FillGroupList();

	wxTextCtrl* txt_mset_name = (wxTextCtrl*) FindWindow(IDC_EDTGRP_MOLSET_NAME);
	txt_mset_name->SetValue(pmset->GetName());

	return wxFrame::TransferDataToWindow();
}

int EditGroupsDlg::ShowModal()
{
	modal_run_flag = TRUE;
	wxFrame* frame_p = (wxFrame*)this;
	frame_p->Show();

	p_loc_event_loop = new wxEventLoop();
	p_loc_event_loop->Run();
	return TRUE;
}

BEGIN_EVENT_TABLE(EditGroupsDlg, wxFrame)
    EVT_RADIOBOX(IDC_TYPE_GROUP, EditGroupsDlg::OnGroupType)
	EVT_MENU(IDC_SET_SEL_FROM_GRP, EditGroupsDlg::OnSetSelFromGrp)
	EVT_MENU(IDC_ADD_SEL_TO_GRP, EditGroupsDlg::OnAddSelToGrp)
	EVT_MENU(IDC_DEL_SEL_FROM_GRP, EditGroupsDlg::OnDelSelFromGrp)
	EVT_MENU(IDC_CALC_DON_ACC_DIST, EditGroupsDlg::OnCalcDonAccDist)
	EVT_BUTTON(IDC_EDTGRP_RENAME_GRP, EditGroupsDlg::OnRenameGrp)
	EVT_BUTTON(IDC_EDTGRP_COPY_GRP, EditGroupsDlg::OnCopyGrp)
	EVT_BUTTON(IDC_EDTGRP_SET_FROM_EXPR, EditGroupsDlg::OnSetFromExpr)
	EVT_BUTTON(IDC_EDTGRP_ADD_FROM_EXPR, EditGroupsDlg::OnAddFromExpr)
	EVT_BUTTON(IDC_EDTGRP_DEL_EXPR, EditGroupsDlg::OnDelFromExpr)
	EVT_BUTTON(IDC_EDTGRP_KEEP_ONLY_EXPR, EditGroupsDlg::OnKeepOnlyExpr)

//	EVT_RADIOBOX(IDC_PICKING_LEVEL, EditGroupsDlg::OnPickLevel)
	EVT_MENU(IDC_CREATE_EXT_CH_MOL, EditGroupsDlg::OnCreateExtChMol)
	EVT_BUTTON (IDC_CLOSE, EditGroupsDlg::OnCloseBtn)
	EVT_CLOSE( EditGroupsDlg::OnClose)
	EVT_BUTTON (IDC_EDTGRP_NEW_GRP, EditGroupsDlg::OnNewGroup)
	EVT_BUTTON (IDC_EDTGRP_DEL_GRP, EditGroupsDlg::OnDelGroup)
	EVT_BUTTON (IDU_ATOM_PICK, EditGroupsDlg::OnAtomSelect)
	EVT_BUTTON (IDC_EDTGRP_RESET_MEMB, EditGroupsDlg::OnResetMemb)
	EVT_BUTTON (IDC_EDTGRP_DELETE_MEMB, EditGroupsDlg::OnDelMemb)
	EVT_BUTTON (IDC_EDTGRP_SET_PROT, EditGroupsDlg::OnSetProt)
	EVT_BUTTON (IDC_EDTGRP_SAVE_XYZ_FILE, EditGroupsDlg::OnSaveXYZFile )
	EVT_BUTTON(IDC_EDTGRP_SAVE_NDX_FILE, EditGroupsDlg::OnSaveNDXFile)
	EVT_BUTTON(IDC_EDTGRP_SORT_IDX, EditGroupsDlg::OnSortGrpIdx)
	EVT_MENU (IDC_STD_GROUPS, EditGroupsDlg::OnStdGroups)
	EVT_MENU(IDC_STD_PROTEIN_GROUPS, EditGroupsDlg::OnStdProteinGroups)
	EVT_MENU(IDC_CHECK_CHEM_GROUPS, EditGroupsDlg::OnCheckChemGroups)
	EVT_MENU(IDC_RENUMBER_GRP, EditGroupsDlg::OnRenumberGrp)
	EVT_MENU(IDC_COLOR_RIGID_CLUSTERS, EditGroupsDlg::OnColorRigidClusters)
	EVT_LISTBOX(IDC_EDTGRP_LIST_GRP, EditGroupsDlg::OnChangeSelGroup2)
END_EVENT_TABLE()

void EditGroupsDlg::OnCloseBtn(wxCommandEvent& event)
{
	Close();
}

void EditGroupsDlg::OnChangeSelGroup2(wxCommandEvent& event)
{
	OnChangeSelGroup();
}


void EditGroupsDlg::OnGroupType(wxCommandEvent& event)
{
	PrintLog("EditGroupsDlg::OnGroupType() \n");
	wxRadioBox* group_type_box = (wxRadioBox*) FindWindow(IDC_TYPE_GROUP);
	igrp_type = group_type_box->GetSelection();
//  wxRadioButton* radio_chem_grp  = (wxRadioButton*) FindWindow(IDC_RADIO_CHEM_GRP);
//	wxRadioButton* radio_redox_grp = (wxRadioButton*) FindWindow(IDC_RADIO_REDOX_GRP);
//	wxRadioButton* radio_named_grp = (wxRadioButton*) FindWindow(IDC_RADIO_NAMED_ATOM_GRP);

//	bool bval[3] = { false, false, false };
//	if( event.GetId() == IDC_RADIO_CHEM_GRP) bval[0] = true;
//	if( event.GetId() == IDC_RADIO_REDOX_GRP) bval[1] = true;
//	if( event.GetId() == IDC_RADIO_NAMED_ATOM_GRP) bval[2] = true;

//	radio_chem_grp->SetValue( bval[0]);
//	radio_redox_grp->SetValue( bval[1]);
//	radio_named_grp->SetValue( bval[2]);
		
	FillGroupList();
}

void EditGroupsDlg::OnTestRadio(wxCommandEvent& event)
{
	PrintLog("EditGroupsDlg::OnTestRadio() \n");
}

void EditGroupsDlg::OnPickLevel(wxCommandEvent& event)
{
	UpdateNumAtomCounter();
	DisplaySelectedGroup();
	event.Skip();
}


void EditGroupsDlg::FillGroupList()
{
	wxListBox* grp_list= (wxListBox*) FindWindow(IDC_EDTGRP_LIST_GRP);
	grp_list->Clear();
	int idx = 0;
	if( igrp_type == CHEM_GRP_TYPE )
	{
		ChemGroup* gptr;
		ChemGroupIterator gitr(pmset);
		for(gptr= gitr.GetFirst(); gptr; gptr= gitr.GetNext())
		{
			grp_list->Append(gptr->GetID() );
			grp_list->SetClientData(idx, (void*)gptr);
			idx++;
		}
	}
	else if( igrp_type == NAMED_ATOM_GROUP_TYPE )
	{
		AtomGroup* at_arr;
		AtomGroupIteratorMolSet litr(pmset);
		for(at_arr= litr.GetFirst(); at_arr; at_arr= litr.GetNext())
		{
			grp_list->Append(at_arr->GetID() );
			grp_list->SetClientData(idx,(void*)at_arr);
			idx++;
		}
	}
	else if( igrp_type == REDOX_CNT_TYPE )
	{
        AtomGroup* pdon = pmset->GetAtomGroupByID("DONOR");
		if(pdon == NULL) pdon = pmset->AddAtomGroup("DONOR");
        AtomGroup* pacc = pmset->GetAtomGroupByID("ACCEPTOR");
		if(pacc == NULL) pacc = pmset->AddAtomGroup("ACCEPTOR");
		grp_list->Append("DONOR");
		grp_list->Append("ACCEPTOR");
        grp_list->SetClientData(0,pdon);
		grp_list->SetClientData(1,pacc);
	}

	wxListBox* memb_list= (wxListBox*) FindWindow( IDC_EDTGRP_LIST_MEMB );
	memb_list->Clear();

	int n = grp_list->GetCount();
	if( n > 0) 
	{
		grp_list->SetSelection(0);
		OnChangeSelGroup();
	}
}

void EditGroupsDlg::UpdateNumAtomCounter()
{
	wxTextCtrl* txt_natoms = (wxTextCtrl*) FindWindow(IDC_EDTGRP_NATOMS);
	AtomGroup* sel_at_arr = GetSelGroup(); 
	if( sel_at_arr )
	{
		txt_natoms->SetValue( wxString::Format("%d",sel_at_arr->GetNAtoms()) );
	}
	else
	{
		txt_natoms->SetValue("");
	}
}

void EditGroupsDlg::OnClose(wxCloseEvent& event)
{
	if( modal_run_flag )
	{
		p_loc_event_loop->Exit();
//		delete p_loc_event_loop;
		p_loc_event_loop = NULL;
//		this->MakeModal(false);
		modal_run_flag = FALSE;
	}
    dlg_open = FALSE;
    active_dlg_ptr = NULL;

//	HaMolView* pview= pmset->GetActiveMolView();
//	if(pview)
//	{
//		pmset->SelectAtomsAll();
//		pview->DisableSpacefill();
//		pmset->RefreshAllViews();
//	}
	event.Skip();
}

void EditGroupsDlg::OnNewGroup(wxCommandEvent& event)
{
	wxListBox* grp_list= (wxListBox*) FindWindow(IDC_EDTGRP_LIST_GRP);
	wxTextCtrl* edt_grp_id = (wxTextCtrl*) FindWindow(IDC_GRP_ID);

	wxString str_wx = edt_grp_id->GetValue();
	
	std::string grp_id_str = str_wx.ToStdString();
	boost::to_upper(grp_id_str);
	boost::trim(grp_id_str);
    edt_grp_id->SetValue(grp_id_str.c_str());
	if(grp_id_str.length() == 0) return;

	AtomGroup* set_at_arr;

	if(igrp_type == CHEM_GRP_TYPE)
	{
	    set_at_arr = pmset->AddBlankChemGroup();
	}
	else if( igrp_type == NAMED_ATOM_GROUP_TYPE )
	{
	    set_at_arr = pmset->AddAtomGroup(grp_id_str.c_str());
	}
	else if( igrp_type == REDOX_CNT_TYPE )
	{
        return;
	}

	if(set_at_arr)
	{
	   grp_list->Append(set_at_arr->GetID() );
       int n = grp_list->GetCount();
	   grp_list->SetClientData(n-1,set_at_arr);
	   grp_list->SetSelection(n-1);
	}
	OnChangeSelGroup();
}

void EditGroupsDlg::OnDelGroup(wxCommandEvent& event)
{
	if(pmset)
	{
		wxListBox* grp_list= (wxListBox*) FindWindow(IDC_EDTGRP_LIST_GRP);
		int ig_cur= grp_list->GetSelection();
		AtomGroup* at_grp = GetSelGroup();

		wxString str_wx;
		
		if(ig_cur >= 0 && at_grp != NULL)
		{
			str_wx = grp_list->GetStringSelection();
			std::string gid  = str_wx.ToStdString();
			if( igrp_type == CHEM_GRP_TYPE)
			{
			   pmset->DeleteChemGroupPtr((ChemGroup*)at_grp);
			}
			else if( igrp_type == NAMED_ATOM_GROUP_TYPE)
			{
               pmset->DeleteAtomGroupPtr(at_grp);
			}
			TransferDataToWindow();
			int n=grp_list->GetCount();
			if(ig_cur  >= n) ig_cur--;
			grp_list->SetSelection(ig_cur);
			OnChangeSelGroup();
		}	
	}
}

AtomGroup* EditGroupsDlg::GetSelGroup()
{
	if(pmset == NULL) return NULL;
	wxListBox* grp_list  = (wxListBox*) FindWindow(IDC_EDTGRP_LIST_GRP);
	int ig_cur          = grp_list->GetSelection();
	if(ig_cur < 0)
		return NULL;
	
	AtomGroup* at_grp = (AtomGroup*) grp_list->GetClientData(ig_cur);

//	grp_list->GetText(ig_cur, buf);
//	std::string gid     = buf;
//  if( igrp_type == CHEM_GRP_TYPE )
//	{
//	   ChemGroup* gcur_ptr = pmset->GetChemGroupByID(gid);
//	   return(gcur_ptr);
//	}
//	else if( igrp_type == NAMED_ATOM_GROUP_TYPE || igrp_type == REDOX_CNT_TYPE )
//	{
//		AtomGroup* alist = pmset->GetAtomGroupByID(gid);
//		return alist;
//	}
	return at_grp;
}

void EditGroupsDlg::OnAtomSelect(wxCommandEvent& event)
{
	PrintLog("EditGroupsDlg OnAtomSelect got IDU_ATOM_PICK message \n");
	char buf[256];
	if(PkAtom)
	{
	    wxRadioBox* pick_level_box = (wxRadioBox*) FindWindow(IDC_PICKING_LEVEL);
	    int isel = pick_level_box->GetSelection();

		wxTextCtrl* txt_expr = (wxTextCtrl*) FindWindow(IDC_EDTGRP_EXPR_TXT);

		MolSet* pmset = PkAtom->GetHostMolSet();

		if(isel == 1)
		{
			HaResidue* pres=PkAtom->GetHostRes();
			pres->FillRef(buf);
		}
		else if( isel == 2)
		{
			ChemGroup* chem_grp = pmset->GetChemGroupByAtom(PkAtom);
			if(chem_grp != NULL)
			{
				chem_grp->FillRef(buf);
			}
			else
			{
				strcpy(buf,"");
			}
		}
		else
		{
			PkAtom->FillRef(buf);
		}
		txt_expr->SetValue(buf);
	}
}

void EditGroupsDlg::OnResetMemb(wxCommandEvent& event)
{
	wxListBox* memb_list= (wxListBox*) FindWindow( IDC_EDTGRP_LIST_MEMB );
	memb_list->Clear();
	AtomGroup* p_sel_grp = this->GetSelGroup();
	if(p_sel_grp) p_sel_grp->clear();

	UpdateNumAtomCounter();
	DisplaySelectedGroup();
}

void EditGroupsDlg::OnDelMemb(wxCommandEvent& event)
{
	wxListBox* memb_list= (wxListBox*) FindWindow( IDC_EDTGRP_LIST_MEMB );
	int n = memb_list->GetCount();
	AtomGroup* p_sel_grp = this->GetSelGroup();

	for(int i=n-1; i >= 0; i--)
	{
		if( memb_list->IsSelected(i) ) 
		{
			HaAtom* aptr = (HaAtom*) memb_list->GetClientData(i);
			if( p_sel_grp )
			{
				p_sel_grp->DeleteAtom(aptr);
			}
			memb_list->Delete(i);
		}
	}

	UpdateNumAtomCounter();
	DisplaySelectedGroup();
}

void EditGroupsDlg::DisplaySelectedGroup()
{
	wxListBox* memb_list= (wxListBox*) FindWindow( IDC_EDTGRP_LIST_MEMB );
	AtomGroup* gcur_ptr= GetSelGroup();
	if(gcur_ptr == NULL) return;
	int n = gcur_ptr->size();

	int i;

	wxRadioBox* group_type_box = (wxRadioBox*) FindWindow(IDC_TYPE_GROUP);
	int isel_gr_type = group_type_box->GetSelection();

	wxRadioBox* pick_level_box = (wxRadioBox*) FindWindow(IDC_PICKING_LEVEL);
	int isel_pick = pick_level_box->GetSelection();

	HaMolView* pView=pmset->GetActiveMolView();
	pmset->SelectAtomsAll();
	pView->ReDrawFlag |= RFRefresh;
	pView->DisableSpacefill();
	if(isel_gr_type == CHEM_GRP_TYPE || isel_pick == 2)
	{
		pView->GroupsColourAttrib();
	}
	else
	{
		pView->CPKColourAttrib();
	}
	pView->ReDrawFlag |= RFColour;
	pmset->UnSelectAtomsAll();

	wxString str_wx;

	for(i=0; i< n; i++)
	{
		(*gcur_ptr)[i]->Select();
	}

	pView->SetAtomScreenRadVal(SelectRad);
	pmset->RefreshAllViews();
}

void EditGroupsDlg::OnSetProt(wxCommandEvent& event)
{
	bool bres;
	wxTextCtrl* edit_prot= (wxTextCtrl*) FindWindow(IDC_EDTGRP_EDT_PROT);
	edit_prot->SetValue("");
	if( igrp_type != CHEM_GRP_TYPE ) return;
	
	ChemGroup* gcur_ptr= (ChemGroup*) GetSelGroup();
	wxString str;
	if(gcur_ptr)
	{
		str = edit_prot->GetValue();
		double prot;
		bres = str.ToDouble(&prot);
		gcur_ptr->SetProtect(prot);
		str.Printf("%6.4f",gcur_ptr->GetProtect());
		edit_prot->SetValue(str);
	}
}

void EditGroupsDlg::OnSaveXYZFile(wxCommandEvent& event)
{
	AtomGroup* p_at_arr = GetSelGroup();
	if( p_at_arr == NULL) return;

	std::string fname_init = pmset->GetName();
	fname_init += "_";
	fname_init += p_at_arr->GetID();
	fname_init += ".xyz";

	wxString fname_out = ::wxFileSelector("Select XYZ file to save Atom Group coordinates",
		::wxGetCwd(), fname_init,"xyz","*.xyz");

	if( fname_out.empty() ) return;

	int ires = p_at_arr->SaveXYZFile( fname_out.ToStdString() );
}

void EditGroupsDlg::OnSaveNDXFile(wxCommandEvent& event)
{
	AtomGroup* p_atgrp = GetSelGroup();
	if (p_atgrp == NULL) return;

	wxString fname_init = pmset->GetName();
	fname_init += "_";
	fname_init += p_atgrp->GetID();
	fname_init += ".ndx";

	wxString fname_out = ::wxFileSelector("Select NDX file to save Atom Group coordinates",
		::wxGetCwd(), fname_init, "ndx", "*.ndx");

	if (fname_out.empty()) return;
	pmset->SaveAtomGroupToNDXFile(p_atgrp, fname_out.ToStdString());
}

void EditGroupsDlg::OnSortGrpIdx(wxCommandEvent& event)
{
	AtomGroup* p_atgrp = GetSelGroup();
	if (p_atgrp == NULL) return;
	pmset->SortAtomGroupByIdx(p_atgrp);

	OnChangeSelGroup();
}


void EditGroupsDlg::OnChangeSelGroup()
{
//	cout << " EditGroupsDlg::OnChangeSelGroup() " << endl;
	char buf[256];
	wxListBox* grp_list = (wxListBox*) FindWindow( IDC_EDTGRP_LIST_GRP );
    AtomGroup* sel_at_arr = GetSelGroup(); 
	
	int idx = 0;
	if( grp_list->GetCount() && sel_at_arr )
	{
		wxListBox* memb_list= (wxListBox*) FindWindow( IDC_EDTGRP_LIST_MEMB );
		memb_list->Clear();
		HaAtom* aptr;
		AtomIteratorAtomGroup aitr(sel_at_arr);
		for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
		{
			aptr->FillRef(buf);
			memb_list->Append(buf);
			memb_list->SetClientData(idx, (void*)aptr);
			idx++;
		}
		
		wxTextCtrl* edit_prot= (wxTextCtrl*) FindWindow(IDC_EDTGRP_EDT_PROT);
        edit_prot->SetValue("");
        if(igrp_type == CHEM_GRP_TYPE)
		{
		   ChemGroup* chem_grp = (ChemGroup*) sel_at_arr;
		   sprintf(buf,"%6.4f", chem_grp->GetProtect());
		   edit_prot->SetValue(buf);
		}
	}
	
	UpdateNumAtomCounter();
	DisplaySelectedGroup();
}


void EditGroupsDlg::OnStdGroups(wxCommandEvent& event)
{
	if(pmset)pmset->SetStdChemGroups();
	TransferDataToWindow();
}

void EditGroupsDlg::OnStdProteinGroups(wxCommandEvent& event)
{
	if (pmset)pmset->SetStdProteinGroups();
	TransferDataToWindow();
}

void EditGroupsDlg::OnCheckChemGroups(wxCommandEvent& event)
{
	if (pmset)pmset->CheckChemGroups();
	TransferDataToWindow();
}

void EditGroupsDlg::OnRenumberGrp(wxCommandEvent& event)
{
	PrintLog(" EditGroupsDlg::OnRenumberGrp() \n"); 
	if(pmset)pmset->RenumberGrp();
	TransferDataToWindow();
}
  
void EditGroupsDlg::OnColorRigidClusters(wxCommandEvent& event)
{
	PrintLog(" EditGroupsDlg::OnRenumberGrp() \n"); 
	if(!pmset) return;

	HaMolView* pview = pmset->GetActiveMolView();
	if(!pview) return;

	pview->RigidClusterColourAttrib();
	pview->UpdateThisView(RFRefresh);
}


void EditGroupsDlg::OnSetSelFromGrp(wxCommandEvent& event) 
{
	AtomGroup* at_grp = GetSelGroup();
	if(at_grp == NULL) return;
	AtomIteratorMolSet aitr(pmset);
	HaAtom* aptr;
	for(aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		aptr->UnSelect();
	}
    AtomIteratorAtomGroup alitr(at_grp);
    for( aptr = alitr.GetFirstAtom(); aptr; aptr = alitr.GetNextAtom())
	{
		aptr->Select();
	}
}

void EditGroupsDlg::OnAddSelToGrp(wxCommandEvent& event) 
{
	AtomGroup* at_grp = GetSelGroup();
	if(at_grp == NULL) return;
	AtomIteratorMolSet aitr(pmset);
	HaAtom* aptr;
	for(aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		if(aptr->Selected() && !at_grp->HasAtom(aptr))
			at_grp->InsertAtom(aptr);
	}
	OnChangeSelGroup();
}

void EditGroupsDlg::OnDelSelFromGrp(wxCommandEvent& event) 
{
	AtomGroup* at_grp = GetSelGroup();
	if(at_grp == NULL) return;
	at_grp->DelSelAtoms();
	OnChangeSelGroup();	
}


void EditGroupsDlg::OnCalcDonAccDist(wxCommandEvent& event) 
{
	ETCouplMod* etmod = pmset->GetETCouplMod(1);
	etmod->calc_edge_dist();
}

void EditGroupsDlg::OnRenameGrp(wxCommandEvent& event) 
{
	AtomGroup* at_grp = GetSelGroup();
	if(at_grp == NULL) return;
	wxTextCtrl* edit_grp_id = (wxTextCtrl*) FindWindow(IDC_GRP_ID);
	wxString str = edit_grp_id->GetValue();
	std::string grp_id = str.ToStdString();
	boost::to_upper(grp_id);
	boost::trim(grp_id);
    edit_grp_id->SetValue(grp_id.c_str());
	if(grp_id.length() == 0) return;
	at_grp->SetID(grp_id.c_str());

	wxListBox* grp_list= (wxListBox*) FindWindow(IDC_EDTGRP_LIST_GRP);
	int ig_cur= grp_list->GetSelection();
    FillGroupList();	
	grp_list->SetSelection(ig_cur);
	OnChangeSelGroup();
}

void EditGroupsDlg::OnCopyGrp(wxCommandEvent& event) 
{
	AtomGroup* old_atlist = GetSelGroup();
	if(old_atlist == NULL) return;
    OnNewGroup(event);
	AtomGroup* new_atlist = GetSelGroup();
	if(new_atlist == NULL) return;

	HaAtom* aptr;
	AtomIteratorAtomGroup aitr(old_atlist);
	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		new_atlist->InsertAtom(aptr);
	}
    OnChangeSelGroup();
}

void EditGroupsDlg::OnSetFromExpr(wxCommandEvent& event)
{
	AtomGroup* p_at_arr = GetSelGroup();
	if(p_at_arr == NULL) return;

	wxTextCtrl* txt_expr = (wxTextCtrl*) FindWindow(IDC_EDTGRP_EXPR_TXT);
	std::string expr_str = txt_expr->GetValue().ToStdString();

	p_at_arr->SetFromExprStr(expr_str.c_str(),pmset);
	OnChangeSelGroup();
}

void EditGroupsDlg::OnAddFromExpr(wxCommandEvent& event)
{
	AtomGroup* p_at_arr = GetSelGroup();
	if(p_at_arr == NULL) return;

	wxTextCtrl* txt_expr = (wxTextCtrl*) FindWindow(IDC_EDTGRP_EXPR_TXT);
	std::string expr_str = txt_expr->GetValue().ToStdString();

	p_at_arr->AddFromExprStr(expr_str.c_str(),pmset);
	OnChangeSelGroup();
}

void EditGroupsDlg::OnDelFromExpr(wxCommandEvent& event)
{
	AtomGroup* p_at_arr = GetSelGroup();
	if(p_at_arr == NULL) return;

	wxTextCtrl* txt_expr = (wxTextCtrl*) FindWindow(IDC_EDTGRP_EXPR_TXT);
	std::string expr_str = txt_expr->GetValue().ToStdString();

	p_at_arr->DeleteAtomsExprStr(expr_str.c_str(),pmset);
	OnChangeSelGroup();
}

void EditGroupsDlg::OnKeepOnlyExpr(wxCommandEvent& event)
{
	AtomGroup* p_at_arr = GetSelGroup();
	if(p_at_arr == NULL) return;

	wxTextCtrl* txt_expr = (wxTextCtrl*) FindWindow(IDC_EDTGRP_EXPR_TXT);
	std::string expr_str = txt_expr->GetValue().ToStdString();

	p_at_arr->KeepOnlyAtomsExprStr(expr_str.c_str(),pmset);
	OnChangeSelGroup();
}

void EditGroupsDlg::OnCreateExtChMol(wxCommandEvent& event) 
{
	AtomGroup* atgrp = GetSelGroup();
	if(atgrp == NULL)
	{
		PrintLog("Error in EditGroupsDlg::OnCreateExtChMol() \n");
		PrintLog(" No named group selected \n");
		return;
	}
	std::string grp_name = atgrp->GetID();
	pmset->CreateAxxMol("EXTERNAL_CHARGES",grp_name.c_str());
	PrintLog("Create Axxiliary Molecule with external charges");
}

////////////////////////////////////////////////////////////////////
// CompAccountsDlg Dialog:

int CompAccountsDlg::dlg_open = FALSE;

CompAccountsDlg::CompAccountsDlg(wxWindow* parent):
wxFrame( parent, -1, "Remote Computer Accounts")
{
    dlg_open = TRUE;
	
	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

	comp_acc_dlg(this,TRUE);

}

CompAccountsDlg::~CompAccountsDlg()
{
        dlg_open = FALSE;
}



bool
CompAccountsDlg::TransferDataToWindow()
{
	char buf[256];
	if(pApp != NULL)
	{
		wxListBox* acc_list=  (wxListBox*) FindWindow(IDC_ACCOUNTS_LIST);
		
		acc_list->Clear();
		int i;
		int na = pApp->comp_accounts.size();
		for(i= 0 ; i < na; i++)
		{
			sprintf(buf,"%-15s %-20s %-15s ", pApp->comp_accounts[i].acc_id.c_str(),
				pApp->comp_accounts[i].login_str.c_str(),
				pApp->comp_accounts[i].interm_acc_ID.c_str());
			acc_list->Append(buf);
		}
	}
	return true;
}

BEGIN_EVENT_TABLE(CompAccountsDlg, wxFrame)
    EVT_BUTTON (IDC_ACC_SHOW_LOAD, CompAccountsDlg::OnAccShowLoad)
	EVT_BUTTON (IDC_CLOSE,  CompAccountsDlg::OnCloseBtn)
END_EVENT_TABLE()



void CompAccountsDlg::OnClose(wxCloseEvent& event) 
{
	dlg_open = FALSE;
	event.Skip();
}

void 
CompAccountsDlg::OnCloseBtn(wxCommandEvent& event) 
{
	Close();
}


void CompAccountsDlg::OnAccShowLoad(wxCommandEvent& event) 
{
	pApp->ShowAccountsLoad();
}

//////////////////////////////////////////////////////////////////////
// ResDBDlg dialog

int ResDBDlg::dlg_open = FALSE;

ResDBDlg::ResDBDlg(wxWindow* parent):
wxFrame( parent, -1, "Residue Templates DB")
{
    dlg_open = TRUE;
	
	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

    wxMenuBar* res_db_menu_bar = res_db_menu();
    SetMenuBar(res_db_menu_bar);    

	res_db_dlg(this,TRUE);
	
	sel_templ = NULL;
	OnInitDialog();
}

ResDBDlg::~ResDBDlg()
{
    dlg_open = FALSE;
}

void
ResDBDlg::OnInitDialog()
{
	UpdateTemplList();
}

BEGIN_EVENT_TABLE(ResDBDlg, wxFrame)
	EVT_BUTTON(IDC_CLOSE, ResDBDlg::OnCloseBtn)
	EVT_CLOSE( ResDBDlg::OnClose )
	EVT_MENU(IDC_RESDB_UPDATE_LIST, ResDBDlg::OnUpdateTemplList)
	EVT_MENU(IDC_RESDB_LOAD, ResDBDlg::OnResDBLoad)
	EVT_LISTBOX(IDC_TEMPL_LIST, ResDBDlg::OnChangeSelTempl)
END_EVENT_TABLE()

void
ResDBDlg::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}

void 
ResDBDlg::OnCloseBtn(wxCommandEvent& event) 
{
	Close();
}

void 
ResDBDlg::OnUpdateTemplList(wxCommandEvent& event) 
{
	UpdateTemplList();
}

void 
ResDBDlg::UpdateTemplList() 
{
	wxListBox* templ_list = (wxListBox*) FindWindow( IDC_TEMPL_LIST );
	templ_list->Clear();
	HaResDB* p_res_db = HaResDB::GetDefaultResDB();
	int nt = p_res_db->GetNMol();
	int it;
	int idx = -1;
	for(it = 0 ; it < nt; it++)
	{
		HaMolecule* ptempl = p_res_db->GetMolByIdx(it);
		idx++;
		templ_list->Append(ptempl->GetObjName());
		templ_list->SetClientData(idx,(void*)ptempl);
	}
}


void ResDBDlg::OnChangeSelTempl(wxCommandEvent& event)
{
	wxListBox* templ_list = (wxListBox*) FindWindow( IDC_TEMPL_LIST );
	wxTextCtrl* sel_templ_edit = (wxTextCtrl*) FindWindow(IDC_SELECTED_TEMPL);  
	
	sel_templ_edit->SetValue("");
	sel_templ = NULL;

	int icur_sel = templ_list->GetSelection();
	if(icur_sel < 0)
		return;
	void* ptr = templ_list->GetClientData(icur_sel);
	if( ptr == NULL ) return;
	sel_templ = (HaMolecule*) ptr;

	sel_templ_edit->SetValue(sel_templ->GetObjName());

	HaResDB* p_res_db = HaResDB::GetDefaultResDB();

	AtomIteratorMolSet aitr_mset( p_res_db );
	AtomIteratorMolecule aitr_mol(sel_templ);
	HaAtom* aptr;
	for( aptr = aitr_mset.GetFirstAtom(); aptr; aptr = aitr_mset.GetNextAtom())
	{
        aptr->UnSelect();
	}

	for( aptr = aitr_mol.GetFirstAtom(); aptr; aptr = aitr_mol.GetNextAtom())
	{
        aptr->Select();
	}

	HaMolView* pview = p_res_db->GetActiveMolView();
	if(pview == NULL) return;
	pview->RestrictSelected();
	pview->CenterSelected();

	pview->ReDrawFlag |= (RFRefresh | RFColour | RFApply);
	p_res_db->RefreshAllViews();
}

void 
ResDBDlg::OnResDBLoad(wxCommandEvent& event) 
{
	HaResDB* p_res_db = HaResDB::GetDefaultResDB();	
	UpdateTemplList();
	HaMolView* pView = p_res_db->GetActiveMolView();
    pView->ReDrawFlag |= (RFInitial| RFApply | RFColour);
    pView->InitialTransform();
    pView->DefaultRepresentation();	
	p_res_db->RefreshAllViews();
}

/////////////////////////////////////////////////////////////////////
// MolViewParDlg Dialog:

int MolViewParDlg::dlg_open = FALSE;

MolViewParDlg::MolViewParDlg(HaMolView* new_pview, wxWindow* parent):
wxFrame( parent, -1, "Molecular View Parameters")
{
	pview= new_pview;

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

	mol_view_param_dlg(this,TRUE);

	OnInitDialog();
    dlg_open = TRUE;
}

MolViewParDlg::~MolViewParDlg()
{
        dlg_open = FALSE;
}

void
MolViewParDlg::OnInitDialog()
{
	if(pview != NULL)
	{
		DDX_Text_double(this, IDC_MVPAR_CENX, pview->CenX,"%9.4f");
		DDX_Text_double(this, IDC_MVPAR_CENY, pview->CenY,"%9.4f");
		DDX_Text_double(this, IDC_MVPAR_CENZ, pview->CenZ,"%9.4f");

		DDX_Text_double(this, IDC_MVPAR_ORIGX, pview->Orig(1),"%9.4f");
		DDX_Text_double(this, IDC_MVPAR_ORIGY, pview->Orig(2),"%9.4f");
		DDX_Text_double(this, IDC_MVPAR_ORIGZ, pview->Orig(3),"%9.4f");

		DDX_Text_double(this, IDC_MVPAR_ROT11, pview->Rot(1,1),"%9.4f");
		DDX_Text_double(this, IDC_MVPAR_ROT12, pview->Rot(1,2),"%9.4f");
		DDX_Text_double(this, IDC_MVPAR_ROT13, pview->Rot(1,3),"%9.4f");
		DDX_Text_double(this, IDC_MVPAR_ROT21, pview->Rot(2,1),"%9.4f");
		DDX_Text_double(this, IDC_MVPAR_ROT22, pview->Rot(2,2),"%9.4f");
		DDX_Text_double(this, IDC_MVPAR_ROT23, pview->Rot(2,3),"%9.4f");
		DDX_Text_double(this, IDC_MVPAR_ROT31, pview->Rot(3,1),"%9.4f");
		DDX_Text_double(this, IDC_MVPAR_ROT32, pview->Rot(3,2),"%9.4f");
		DDX_Text_double(this, IDC_MVPAR_ROT33, pview->Rot(3,3),"%9.4f");
		
		DDX_Text_int(this, IDC_MVPAR_SLAB_MODE,   pview->pCanv->m_SlabMode);
		DDX_Text_int(this, IDC_MVPAR_SLAB_VALUE,  pview->pCanv->m_SlabValue);
		DDX_Text_int(this, IDC_MVPAR_SLAB_INTEN,  pview->pCanv->m_SlabInten);
		DDX_Text_int(this, IDC_MVPAR_SLICE_VALUE, pview->pCanv->m_SliceValue);
		DDX_Text_int(this, IDC_MVPAR_IMAGE_SIZE,  pview->pCanv->m_ImageSize);
		DDX_Text_int(this, IDC_MVPAR_ZOFFSET, pview->pCanv->m_ZOffset);
		DDX_Text_double(this, IDC_MVPAR_DSCALE, pview->DScale ,"%9.4f");
		DDX_Text_double(this, IDC_MVPAR_ZOOM, pview->Zoom,"%9.4f");
	}
	TransferDataToWindow();
}

bool
MolViewParDlg::TransferDataFromWindow()
{
	wxCheckBox* check_box;

	check_box = (wxCheckBox*) FindWindow( IDC_MVPAR_USE_SLAB_PLANE );
	pview->pCanv->m_UseSlabPlane = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_MVPAR_DRAW_BEST_PATH );
	pview->DrawBestPath = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_MVPAR_DRAW_UNIT_CELL );
	pview->DrawUnitCell = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_MVPAR_DRAW_CONT_SURF );
	pview->DrawContourSurf = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_MVPAR_DRAW_SOLID_SURF );
	pview->DrawSolidSurfaces = check_box->GetValue();

	return wxFrame::TransferDataFromWindow();
}

bool
MolViewParDlg::TransferDataToWindow()
{
	wxCheckBox* check_box;
	
	check_box = (wxCheckBox*) FindWindow( IDC_MVPAR_USE_SLAB_PLANE );
	if( pview->pCanv->m_UseSlabPlane ) check_box->SetValue(true);
	else check_box->SetValue(false);

	check_box = (wxCheckBox*) FindWindow( IDC_MVPAR_DRAW_BEST_PATH );
	if( pview->DrawBestPath ) check_box->SetValue(true);
	else check_box->SetValue(false);

	check_box = (wxCheckBox*) FindWindow( IDC_MVPAR_DRAW_UNIT_CELL );
	if( pview->DrawUnitCell ) check_box->SetValue(true);
	else check_box->SetValue(false);

	check_box = (wxCheckBox*) FindWindow( IDC_MVPAR_DRAW_CONT_SURF );
	if( pview->DrawContourSurf ) check_box->SetValue(true);
	else check_box->SetValue(false);

	check_box = (wxCheckBox*) FindWindow( IDC_MVPAR_DRAW_SOLID_SURF );
	if( pview->DrawSolidSurfaces ) check_box->SetValue(true);
	else check_box->SetValue(false);

	return wxFrame::TransferDataToWindow();
}

void
MolViewParDlg::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}

void 
MolViewParDlg::OnCloseBtn(wxCommandEvent& event) 
{
	Close();
}

// WDR: event table for MolViewParDlg
BEGIN_EVENT_TABLE(MolViewParDlg, wxFrame)
	EVT_BUTTON (IDC_CLOSE,         MolViewParDlg::OnCloseBtn)
	EVT_CLOSE ( MolViewParDlg::OnClose)
	EVT_BUTTON (IDC_MVPAR_RESET_VIEW_ROT,  MolViewParDlg::OnResetViewRot)
	EVT_BUTTON (IDC_MVPAR_SETCNTR,         MolViewParDlg::OnSetControlPar)
	EVT_BUTTON (IDC_MVPAR_SETVIEW,         MolViewParDlg::OnSetViewPar)
END_EVENT_TABLE()


void
MolViewParDlg::OnSetControlPar(wxCommandEvent& event)
{
	TransferDataToWindow();
}

void
MolViewParDlg::OnSetViewPar(wxCommandEvent& event)
{
	TransferDataFromWindow();
	if(pview != NULL)
	{
		MolSet* pmset = pview->GetMolSet();
		pmset->RefreshAllViews(RFRefresh | RFColour);
	}
}

void MolViewParDlg::OnResetViewRot(wxCommandEvent& event) 
{
	if(pview != NULL)
	{
		pview->ReDrawFlag |= RFInitial;
		pview->InitialTransform();
		MolSet* pmset = pview->GetMolSet();
		pmset->RefreshAllViews(RFRefresh | RFMagnify);
	}
	
}

/////////////////////////////////////////////////////////////////////
// MolSetParDlg Dialog:

int MolSetParDlg::dlg_open = FALSE;

MolSetParDlg::MolSetParDlg(MolSet* new_pmset, wxWindow* parent):
wxFrame( parent, -1, "Molecular Set Parameters")
{
	pmset= new_pmset;

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

	mol_set_param_dlg(this,TRUE);

	OnInitDialog();
    dlg_open = TRUE;
}

MolSetParDlg::~MolSetParDlg()
{
     dlg_open = FALSE;
}

void MolSetParDlg::OnInitDialog()
{
	TransferDataToWindow();
}

bool MolSetParDlg::TransferDataFromWindow()
{
	if(pmset != NULL)
	{	
		wxTextCtrl*    name_edt=  (wxTextCtrl*)  FindWindow(IDC_MSPAR_MSET_NAME);
		
		wxString str = name_edt->GetValue();
		pmset->SetName(str.c_str());
	}
	return wxFrame::TransferDataFromWindow();
}

bool MolSetParDlg::TransferDataToWindow()
{
	wxListBox* pmset_list = (wxListBox*) FindWindow(IDC_MSPAR_MOLSET_LIST);
	pmset_list->Clear();
	int nm = pApp->molset_vec.size();
	int i;
	MolSet* pmset_c = GetCurMolSet();
	int idx_mset_c = 0;
	for(i = 0; i < nm; i++)
	{
		MolSet* pmset = (MolSet*) pApp->molset_vec[i];
		if( pmset == pmset_c ) idx_mset_c = i;
		pmset_list->Append( pmset->GetName(), (void*) pmset );
	}
	if( pmset_list->GetCount() > 0) pmset_list->SetSelection(idx_mset_c);

	if(pmset != NULL)
	{
		wxTextCtrl*    name_edt=  (wxTextCtrl*)  FindWindow(IDC_MSPAR_MSET_NAME);
		
		wxListBox* frag_list= (wxListBox*) FindWindow(IDC_MSPAR_FRAGM_NAMES);
		wxListBox* mol_list=  (wxListBox*) FindWindow(IDC_MSPAR_MOL_NAMES);

		name_edt->SetValue(pmset->GetName());
		
		frag_list->Clear();
		mol_list->Clear();
		int nf = pmset->Fragments.size();
		int i;

		for(i=0; i < nf; i++)
		{
			MolSet* frag = (MolSet*) pmset->Fragments[i];
			frag_list->Append( frag->GetName() );
		}
		std::vector<HaMolecule*>::iterator mol_itr;
		for(mol_itr= pmset->HostMolecules.begin(); mol_itr != pmset->HostMolecules.end(); mol_itr++ )
		{
			mol_list->Append( (*mol_itr)->GetObjName() );
		}	
	}
	return wxFrame::TransferDataToWindow();
}

void MolSetParDlg::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}

void MolSetParDlg::OnCloseBtn(wxCommandEvent& event) 
{
	Close();
}

// WDR: event table for MolViewParDlg
BEGIN_EVENT_TABLE(MolSetParDlg, wxFrame)
	EVT_BUTTON (IDC_CLOSE,         MolSetParDlg::OnCloseBtn)
	EVT_CLOSE ( MolSetParDlg::OnClose)
	EVT_BUTTON (IDC_MSPAR_MULT_MOL,  MolSetParDlg::OnMultMol)
	EVT_BUTTON (IDC_MSPAR_SETCNTR,   MolSetParDlg::OnSetControlPar)
	EVT_BUTTON (IDC_MSPAR_SETMOL,    MolSetParDlg::OnSetMolPar)
END_EVENT_TABLE()


void MolSetParDlg::OnSetControlPar(wxCommandEvent& event)
{
	TransferDataToWindow();
}

void MolSetParDlg::OnSetMolPar(wxCommandEvent& event)
{
	TransferDataFromWindow();
}

void MolSetParDlg::OnMultMol(wxCommandEvent& event) 
{
	if(pmset != NULL)
	{
		MolEditor* p_mol_editor = pmset->GetMolEditor();
		p_mol_editor->ReplicatePeriodBox(pmset,2,2,2);
	}	
}



/////////////////////////////////////////////////////////////////////
// AtomPropColorDlg Dialog:

int AtomPropColorDlg::dlg_open = FALSE;

AtomPropColorDlg::AtomPropColorDlg(HaMolView* pview_new, wxWindow* parent):
wxFrame( parent, -1, "Color Atoms By Property")
{
	pview = pview_new;
	pmset= pview->GetMolSet();

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

	atom_prop_colors_dlg(this,TRUE);

	OnInitDialog();
    dlg_open = TRUE;
}

AtomPropColorDlg::~AtomPropColorDlg()
{
     dlg_open = FALSE;
}

void AtomPropColorDlg::OnInitDialog()
{
	n_col_mval  = 0;
	n_col_cname = 1; 
	n_col_r = 2;
	n_col_g = 3;
	n_col_b = 4;
	
	p_atom_prop = (wxChoice*) FindWindow(IDC_CHOICE_ATOM_PROP);
	p_atom_prop->Clear();
	p_atom_prop->Append("Temperature");
	p_atom_prop->Append("Atom RMSF");
	p_atom_prop->Append("Atom RMSD");
	p_atom_prop->SetSelection(0);

	p_color_names = (wxChoice*) FindWindow( IDC_CHOICE_COLORS );
	p_color_names->Clear();

	std::map<std::string, ColorVal>::iterator citr;
	citr = HaColor::name_col_map.begin();
	for( ; citr != HaColor::name_col_map.end(); citr++ )
	{
		std::string col_name = (*citr).first;
		p_color_names->Append(col_name.c_str());
	}
	p_color_names->SetSelection(1);

	p_txt_minval = (wxTextCtrl*) FindWindow( IDC_MIN_PROP_VAL );
	p_txt_minval->SetValue("0.0");

	p_txt_r = (wxTextCtrl*) FindWindow( IDC_COLOR_R );
	p_txt_g = (wxTextCtrl*) FindWindow( IDC_COLOR_G );
	p_txt_b = (wxTextCtrl*) FindWindow( IDC_COLOR_B );
	r_val = 255;
    g_val = 255;
	b_val = 255;

	num_interpol_colors = 0;
	p_txt_num_interpol_colors = (wxTextCtrl*) FindWindow( IDC_NUM_INTERPOL_COLORS );

	p_color_map_grid = (wxGrid*) FindWindow(IDC_GRID_COLOR_MAP);

	int ncol_act = p_color_map_grid->GetNumberCols();
	if(ncol_act > 0 )
	{
		p_color_map_grid->DeleteCols(0,ncol_act);
	}

	p_color_map_grid->InsertCols(0,5);
	p_color_map_grid->SetColLabelValue(n_col_mval, "Min Val");
	p_color_map_grid->SetColMinimalWidth(n_col_mval,30);
	p_color_map_grid->SetColLabelValue(n_col_cname, "Color Name");
	p_color_map_grid->SetColMinimalWidth(n_col_cname,30);
	p_color_map_grid->SetColLabelValue(n_col_r, "R");
	p_color_map_grid->SetColMinimalWidth(n_col_r,10);
	p_color_map_grid->SetColLabelValue(n_col_g, "G");
	p_color_map_grid->SetColMinimalWidth(n_col_g,10);
	p_color_map_grid->SetColLabelValue(n_col_b, "B");
	p_color_map_grid->SetColMinimalWidth(n_col_b,10);

	color_map.AddColorAndMinVal(0,0,255,-100.0);
	color_map.AddColorAndMinVal(0,255,0,1.0);
	color_map.AddColorAndMinVal(255,255,0,1.5);
	color_map.AddColorAndMinVal(255,0,0,2.0);

	TransferDataToWindow();
}

bool AtomPropColorDlg::TransferDataFromWindow()
{
	wxString str_min_val = p_txt_minval->GetValue();
	int ires = str_min_val.ToDouble(&min_prop_val);

	wxString str_r_val = p_txt_r->GetValue();
	wxString str_g_val = p_txt_g->GetValue();
	wxString str_b_val = p_txt_b->GetValue();

	ires = str_r_val.ToLong(&r_val);
	if(!ires || r_val < 0 || r_val > 255) r_val = -1;
	ires = str_g_val.ToLong(&g_val);
	if(!ires || g_val < 0 || g_val > 255) g_val = -1;
	ires = str_b_val.ToLong(&b_val);
	if(!ires || b_val < 0 || b_val > 255) b_val = -1;

	wxString val_str = p_txt_num_interpol_colors->GetValue();
	ires = val_str.ToLong(&num_interpol_colors);

	if(!ires || num_interpol_colors < 0 ) num_interpol_colors = 0;

	return wxFrame::TransferDataFromWindow();
}

bool AtomPropColorDlg::TransferDataToWindow()
{
//	wxListBox* pmset_list = (wxListBox*) FindWindow(IDC_MSPAR_MOLSET_LIST);
//	pmset_list->Clear();
	int nrow = p_color_map_grid->GetNumberRows();
	if( nrow > 0 ) p_color_map_grid->DeleteRows(0,nrow);

	int nc = color_map.GetNColors();

	wxString lbl;
	lbl.Printf("%3d",r_val);
	p_txt_r->SetValue(lbl);
	lbl.Printf("%3d",g_val);
	p_txt_g->SetValue(lbl);
	lbl.Printf("%3d",b_val);
	p_txt_b->SetValue(lbl);

	if( nc > 0 )
	{
		p_color_map_grid->AppendRows(nc);
	}

	lbl.Printf("%d",num_interpol_colors);
	p_txt_num_interpol_colors->SetValue(lbl);

	int i;
	for( i = 0; i < nc; i++ )
	{
		HaColor& color = color_map.GetColorByIdx(i);
		
		lbl.Printf("%3d",color.r);
		p_color_map_grid->SetCellValue(i,2,lbl);
		lbl.Printf("%3d",color.g);
		p_color_map_grid->SetCellValue(i,3,lbl);
		lbl.Printf("%3d",color.b);

		std::string col_name = HaColor::GetColorName(color.r,color.g,color.b);
		if( !col_name.empty() )
		{
			p_color_map_grid->SetCellValue(i,n_col_cname,col_name.c_str());
		}
		wxColour color_wx(color.r,color.g,color.b);
		p_color_map_grid->SetCellBackgroundColour(i,n_col_cname,color_wx);

		p_color_map_grid->SetReadOnly(i,n_col_cname);

		p_color_map_grid->SetCellValue(i,4,lbl);
		if( i == 0 || (i >= color_map.min_values.size()) ) 
		{
			lbl = "";
		}
		else
		{
			lbl.Printf("%9.4f",color_map.min_values[i]);
		}

		p_color_map_grid->SetCellValue(i,0,lbl);
	}
	
	return wxFrame::TransferDataToWindow();
}

void AtomPropColorDlg::OnClose(wxCloseEvent& event)
{
	dlg_open = FALSE;
	event.Skip();
}

void AtomPropColorDlg::OnColorAtomsProp(wxCommandEvent& event)
{
	wxString str_prop = p_atom_prop->GetStringSelection();
	pview->ColorAtomsByProp( str_prop.ToStdString() , &color_map );
}

void AtomPropColorDlg::OnClearColorMap(wxCommandEvent& event)
{
	color_map.Clear();
	TransferDataToWindow();
}

void AtomPropColorDlg::OnAddColorByName(wxCommandEvent& event)
{
	TransferDataFromWindow();

	wxString col_name = p_color_names->GetStringSelection();
	if( !HaColor::ColorNameExist( col_name.ToStdString() ) ) return;
	ColorVal col_val = HaColor::GetColorVal(col_name.ToStdString() );
	HaColor color( col_val );

	color_map.AddColorAndMinVal( color.r, color.g, color.b, min_prop_val, num_interpol_colors );
	TransferDataToWindow();
}

void AtomPropColorDlg::OnAddColorByRGB(wxCommandEvent& event)
{
	TransferDataFromWindow();
	if( r_val < 0 || g_val < 0 || b_val < 0 )
	{
		PrintLog(" Error in AtomPropColorDlg::OnAddColorByRGB() \n");
		PrintLog(" Invalid R G B values \n");
		return;
	}
	color_map.AddColorAndMinVal( r_val, g_val, b_val, min_prop_val, num_interpol_colors );
	TransferDataToWindow();
}

void AtomPropColorDlg::OnSaveColorMapFile(wxCommandEvent& event)
{
	wxString fname_save = ::wxFileSelector("Select txt file to save Color Map",
		                     ::wxGetCwd(),"color_map.txt","txt","*.txt");

	int ires = color_map.SaveToTxtFile( fname_save.ToStdString() );
	TransferDataToWindow();
}

void AtomPropColorDlg::OnLoadColorMapFile(wxCommandEvent& event)
{
	wxString fname_load = ::wxFileSelector("Select txt file to load Color Map",
		                     ::wxGetCwd(),"color_map.txt","txt","*.txt");

	int ires = color_map.LoadFromTxtFile( fname_load.ToStdString() );
	TransferDataToWindow();
}

void AtomPropColorDlg::OnGridCellChange( wxGridEvent &event )
{
	// PrintLog(" Enter AtomPropColorDlg::OnGridEndEdit \n");
	
	int col = event.GetCol();
	int row = event.GetRow();	

	bool bres;
	char buf[256];
	wxString lbl;

    // PrintLog("End Edit on row %d, col %d\n", row, col );

	wxString str_val = p_color_map_grid->GetCellValue(row,col);

	if( col == n_col_r || col == n_col_g || col == n_col_b )
	{
		HaColor& color = color_map.GetColorByIdx(row);
		long ival;
		bres = str_val.ToLong(&ival);
		if( row >= color_map.GetNColors() ) return; 
		
		if(!bres || ival < 0 || ival > 255 ) 
		{
			if( col == n_col_r )
			{
				ival = color.r;
			}
			else if( col == n_col_g )
			{
				ival = color.g;
			}
			else if( col == n_col_b )
			{
				ival = color.b;
			} 
		}

		if( col == n_col_r )
		{
			color.r = (unsigned char) ival;
		}
		else if( col == n_col_g )
		{
			color.g = (unsigned char) ival;
		}
		else if( col == n_col_b )
		{
			color.b = (unsigned char) ival;
		}

		sprintf(buf,"%3d", ival);
		p_color_map_grid->SetCellValue(row, col, buf);

		std::string col_name = HaColor::GetColorName(color.r,color.g,color.b);
		if( !col_name.empty() )
		{
			p_color_map_grid->SetCellValue(row,n_col_cname,col_name.c_str());
		}
		else
		{
			p_color_map_grid->SetCellValue(row,n_col_cname,"");
		}
		wxColour color_wx(color.r,color.g,color.b);
		p_color_map_grid->SetCellBackgroundColour(row,n_col_cname,color_wx);
	}
	if( col == n_col_mval )
	{
		double dval;
		bres = str_val.ToDouble(&dval);

		PrintLog( " dval = %9.4f \n",dval );
		if(  row >= color_map.min_values.size()  ) return;

		if(!bres || row == 0 ||  
			((row > 0) && color_map.min_values[row-1] > dval) ||  
			((row + 1) < color_map.min_values.size() && color_map.min_values[row+1] < dval) ) 
		{
			dval = color_map.min_values[row];
		}
		color_map.min_values[row] = dval;
		
		lbl.Printf("%9.4f",dval);
		p_color_map_grid->SetCellValue(row,n_col_mval,lbl);
	}

	return;
}

BEGIN_EVENT_TABLE(AtomPropColorDlg, wxFrame)
	EVT_CLOSE ( AtomPropColorDlg::OnClose)
	EVT_BUTTON( IDC_COLOR_ATOMS_PROP, AtomPropColorDlg::OnColorAtomsProp )
	EVT_BUTTON( IDC_CLEAR_COLOR_MAP,  AtomPropColorDlg::OnClearColorMap )
	EVT_BUTTON( IDC_ADD_COLOR_BY_NAME, AtomPropColorDlg::OnAddColorByName )
	EVT_BUTTON( IDC_ADD_COLOR_BY_RGB,  AtomPropColorDlg::OnAddColorByRGB )
	EVT_BUTTON( IDC_SAVE_COLOR_MAP_FILE, AtomPropColorDlg::OnSaveColorMapFile )
	EVT_BUTTON( IDC_LOAD_COLOR_MAP_FILE, AtomPropColorDlg::OnLoadColorMapFile )
	EVT_GRID_CELL_CHANGED( AtomPropColorDlg::OnGridCellChange )
END_EVENT_TABLE()


int NuclAcidDlgWX::dlg_open = FALSE;

NuclAcidDlgWX::NuclAcidDlgWX(NuclAcidMod* new_nucl_acid_mod, wxWindow* parent):
wxFrame( parent, -1, "Nuclear Acids Manipulation")
{
	this->SetExtraStyle(wxWS_EX_VALIDATE_RECURSIVELY);
	nucl_acid_mod = new_nucl_acid_mod;
	pmset = nucl_acid_mod->GetMolSet();

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);

    wxMenuBar* nucl_acid_menu_bar = nucl_acid_menu();
    SetMenuBar(nucl_acid_menu_bar);    

	nucl_acid_dlg( this, TRUE );

	display_mode = JUMNA_HLX_CRD;
	m_Grid = (wxGrid*) FindWindow(IDC_GRID);

	OnInitDialog();
	
}

NuclAcidDlgWX::~NuclAcidDlgWX()
{
        dlg_open = FALSE;
}

void
NuclAcidDlgWX::OnInitDialog()
{
	wxTextCtrl* edit_ctrl;
	
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_RES_INDEX);
    edit_ctrl->SetValue("123456789012345678901234567890123456789012345678901234567890");

	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_CHAIN_1);
	edit_ctrl->SetValidator( StdStringValidator(&nucl_acid_mod->seq[0]) );
	edit_ctrl = (wxTextCtrl*) FindWindow(IDC_CHAIN_2);
	edit_ctrl->SetValidator( StdStringValidator(&nucl_acid_mod->seq[1]) );

	DDX_Text_double(this,IDC_SUP_HELIX_RAD,nucl_acid_mod->sup_helix_rad,"%9.4f");
	DDX_Text_double(this,IDC_SUP_HELIX_PIT,nucl_acid_mod->sup_helix_pit,"%9.4f");
    DDX_Text_double(this,IDC_FST_TWIST,    nucl_acid_mod->fst_twist,    "%9.4f");
	DDX_Text_int (this, IDC_SYMM_UNIT_1,   nucl_acid_mod->nsym_unit[0] );
	DDX_Text_int (this, IDC_SYMM_UNIT_2,   nucl_acid_mod->nsym_unit[1] );
	DDX_Text_int (this, IDC_SYMM_UNIT_3,   nucl_acid_mod->nsym_unit[2] );
	DDX_Text_int (this, IDC_BREAK_RES_1,   nucl_acid_mod->nbreak_1     );
	DDX_Text_int (this, IDC_BREAK_RES_2,   nucl_acid_mod->nbreak_2     );
	DDX_Text_int (this, IDC_INTERS_SYMM_OFF_1, nucl_acid_mod->homon_symm_offs[0]  );
	DDX_Text_int (this, IDC_INTERS_SYMM_OFF_2, nucl_acid_mod->homon_symm_offs[1]  );
	DDX_Text_int (this, IDC_INTERS_SYMM_OFF_3, nucl_acid_mod->homon_symm_offs[2]  );
	DDX_Text_int (this, IDC_MAX_ITER,          nucl_acid_mod->max_iter );

	TransferDataToWindow();
	dlg_open = TRUE;
}

bool NuclAcidDlgWX::TransferDataToWindow()
{
	int i;
	for(i= 0; i < 4; i++)
	{
		boost::trim(nucl_acid_mod->seq[i]);
		boost::to_upper(nucl_acid_mod->seq[i]);
	}

	wxCheckBox* check_box;
	
	check_box = (wxCheckBox*) FindWindow( IDC_ENE_PER_UNIT );
	if( nucl_acid_mod->ene_per_unit_flag ) check_box->SetValue(true);
	else check_box->SetValue(false);

	check_box = (wxCheckBox*) FindWindow( IDC_HOMO_SYMM_1 );
	if( nucl_acid_mod->homon_symm_flags[0] ) check_box->SetValue(true);
	else check_box->SetValue(false);

	check_box = (wxCheckBox*) FindWindow( IDC_HOMO_SYMM_2 );
	if( nucl_acid_mod->homon_symm_flags[1] ) check_box->SetValue(true);
	else check_box->SetValue(false);

	check_box = (wxCheckBox*) FindWindow( IDC_HOMO_SYMM_3 );
	if( nucl_acid_mod->homon_symm_flags[2] ) check_box->SetValue(true);
	else check_box->SetValue(false);

	wxChoice* choice_ctrl;

	choice_ctrl = (wxChoice*) FindWindow( IDC_FF_TYPE );
	choice_ctrl->SetSelection(nucl_acid_mod->force_field);

	wxString str_wx;

    wxTextCtrl*     edt_sigma     = (wxTextCtrl*) FindWindow(IDC_SIGM_SLOPE); 
	str_wx.Printf("%9.4f",nucl_acid_mod->diel_slope);
	edt_sigma->SetValue(str_wx);

    wxTextCtrl*     edt_phos_chrg     = (wxTextCtrl*) FindWindow(IDC_PHOS_CHRG); 
	str_wx.Printf("%9.4f",nucl_acid_mod->phos_chrg);
	edt_phos_chrg->SetValue(str_wx);

	OnUpdateContent();

	return wxFrame::TransferDataToWindow();
}

bool
NuclAcidDlgWX::TransferDataFromWindow()
{
	wxCheckBox* check_box;

	check_box = (wxCheckBox*) FindWindow( IDC_ENE_PER_UNIT );
	nucl_acid_mod->ene_per_unit_flag = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_HOMO_SYMM_1 );
	nucl_acid_mod->homon_symm_flags[0] = check_box->GetValue();

	check_box = (wxCheckBox*) FindWindow( IDC_HOMO_SYMM_2 );
	nucl_acid_mod->homon_symm_flags[1] = check_box->GetValue();
	
	check_box = (wxCheckBox*) FindWindow( IDC_HOMO_SYMM_3 );
	nucl_acid_mod->homon_symm_flags[2] = check_box->GetValue();

    int isel;

	wxChoice* choice_ctrl;

	choice_ctrl = (wxChoice*) FindWindow( IDC_FF_TYPE );
	isel = choice_ctrl->GetSelection();
	nucl_acid_mod->SetFFtypeIdx(isel);


    bool bres;
	double tmp;
	wxString str_wx;

    wxTextCtrl*     edt_sigma     = (wxTextCtrl*) FindWindow(IDC_SIGM_SLOPE); 
	str_wx = edt_sigma->GetValue();
	bres = str_wx.ToDouble(&tmp);
    if(bres) nucl_acid_mod->SetDielSlope(tmp);

    wxTextCtrl*     edt_phos_chrg     = (wxTextCtrl*) FindWindow(IDC_PHOS_CHRG); 
	str_wx = edt_phos_chrg->GetValue();
	bres = str_wx.ToDouble(&tmp);
    if(bres) nucl_acid_mod->SetPhosChrg(tmp);

	return wxFrame::TransferDataFromWindow();
}



BEGIN_EVENT_TABLE(NuclAcidDlgWX, wxFrame)
    EVT_BUTTON(IDC_UPDATE_CONTROLS, NuclAcidDlgWX::OnUpdateControls)
	EVT_BUTTON(IDC_SAVE_CHANGES, NuclAcidDlgWX::OnSaveChanges)
	EVT_BUTTON(IDC_BUILD_NUCL_ACID, NuclAcidDlgWX::OnBuildNuclAcid)
	EVT_BUTTON(IDC_MIN_ENE, NuclAcidDlgWX::OnMinEne)
	EVT_BUTTON(IDC_NA_CALC_ENE, NuclAcidDlgWX::OnCalcEne)
	EVT_BUTTON(IDC_GEN_COMPL_STR, NuclAcidDlgWX::OnGenComplStr)
	EVT_MENU(IDC_PDIST_TO_CNT, NuclAcidDlgWX::OnPdistToCnt)
	EVT_MENU(IDC_LOC_BASE_STEP_PARS, NuclAcidDlgWX::OnLocBaseStepPars)
	EVT_MENU(IDC_LOC_HLX_PARS_BP, NuclAcidDlgWX::OnLocHlxParsBp)
	EVT_MENU(IDC_CALC_AXIS, NuclAcidDlgWX::OnCalcAxis)
	EVT_MENU(IDC_GLOB_HLX_PAR, NuclAcidDlgWX::OnGlobHlxPar)
	EVT_MENU(IDC_CALC_BB_CRD, NuclAcidDlgWX::OnCalcBbCrd)
	EVT_MENU(IDC_UPDATE_XYZ, NuclAcidDlgWX::OnUpdateXYZ)
	EVT_MENU(IDC_LOCK_SEL_CRD,NuclAcidDlgWX::OnLockSelCrd)
	EVT_MENU(IDC_BASE_BASE_PARS, NuclAcidDlgWX::OnBaseBasePars)
	EVT_BUTTON (IDC_CLOSE, NuclAcidDlgWX::OnCloseBtn)
	EVT_BUTTON (IDC_SET_FROM_JUMNA, NuclAcidDlgWX::OnSetFromJumna)
	EVT_BUTTON (IDC_SET_TO_JUMNA,   NuclAcidDlgWX::OnSetToJumna)
	EVT_BUTTON (IDC_UPDATE_GRID,    NuclAcidDlgWX::OnUpdateControls)
    EVT_CLOSE (NuclAcidDlgWX::OnClose)
	EVT_GRID_CELL_CHANGED( NuclAcidDlgWX::OnGridEndEdit )
	EVT_NOTEBOOK_PAGE_CHANGED(IDC_NUCL_ACID_NOTEB, NuclAcidDlgWX::OnPageChange)
END_EVENT_TABLE()


void NuclAcidDlgWX::OnClose(wxCloseEvent& event)
{
	event.Skip();
}

void NuclAcidDlgWX::OnCloseBtn(wxCommandEvent& event)
{
	Close();
}

void NuclAcidDlgWX::OnPageChange(wxNotebookEvent& event)
{
    TransferDataToWindow();	
}


void NuclAcidDlgWX::OnUpdateControls(wxCommandEvent& event) 
{
    TransferDataToWindow();	
}

void NuclAcidDlgWX::OnSaveChanges(wxCommandEvent& event) 
{
    TransferDataFromWindow();	
}

void NuclAcidDlgWX::OnBuildNuclAcid(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	nucl_acid_mod->BuildNuclAcid();

	HaMolView* pView = pmset->GetActiveMolView();
	if (!pView) return;

	pView->ReDrawFlag |= RFInitial;
	pView->InitialTransform();
	pView->DefaultRepresentation();
	pmset->RefreshAllViews();
}

void NuclAcidDlgWX::OnPdistToCnt(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	nucl_acid_mod->CalcPdistToShlxCnt();
}


void NuclAcidDlgWX::OnMinEne(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	wxBusyCursor wait;
	nucl_acid_mod->MinEne();
}

void NuclAcidDlgWX::OnCalcEne(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	wxBusyCursor wait;
	nucl_acid_mod->CalcEne(); 	
}

void NuclAcidDlgWX::OnGenComplStr(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	nucl_acid_mod->GenComplStrand();
	TransferDataToWindow();
}


void NuclAcidDlgWX::OnLocBaseStepPars(wxCommandEvent& event) 
{
	nucl_acid_mod->CalcLocHlxCrd(FALSE);
}

void NuclAcidDlgWX::OnLocHlxParsBp(wxCommandEvent& event) 
{
	nucl_acid_mod->CalcLocHlxCrd(TRUE);

}

void NuclAcidDlgWX::OnGlobHlxPar(wxCommandEvent& event) 
{
	nucl_acid_mod->CalcGlobHlxCrd();
}

void NuclAcidDlgWX::OnCalcAxis(wxCommandEvent& event) 
{
	nucl_acid_mod->CalcAxis();
	nucl_acid_mod->CalcBend();
}

void NuclAcidDlgWX::OnCalcBbCrd(wxCommandEvent& event) 
{
	nucl_acid_mod->CalcBBCrd();	
}


void NuclAcidDlgWX::OnUpdateContent()
{
	HaMolecule* pmol = nucl_acid_mod->p_dna_mol;
	if(pmol == NULL) return;

	HaChain* ch;
	ChainIteratorMolecule chitr(pmol);
	HaChain* ch1 = chitr.GetFirstChain();
	HaChain* ch2 = chitr.GetNextChain();
	int nr = ch1->GetNRes();

	if(display_mode == BASE_BASE_CRD)
	{
		if(nucl_acid_mod->bs_bs_pars.num_cols() < nr) nucl_acid_mod->CalcGlobHlxCrd();
	}

	int m_nCols = 0;
	int m_nRows = 0;

	if(display_mode == JUMNA_HLX_CRD && nucl_acid_mod->hel_crd.size() > 6*nr)
	{
		m_nCols = 6;
		m_nRows = nucl_acid_mod->GetNRes();
	}
	else if(display_mode == BASE_BASE_CRD)
	{
		m_nCols = 6;
		m_nRows = nr;		
	}

	StrVec lbl_res(nr*4);

	HaResidue* pres;

	int m = -1;
	int k;
	for(k = 0; k < 2; k++)
	{
		if(k == 0) ch = ch1;
		if(k == 1) ch = ch2;
		ResidueIteratorChain ritr_ch(ch);
		pres = ritr_ch.GetFirstRes();
		for(; pres; pres = ritr_ch.GetNextRes())
		{
			m++;
			char buf[256];
			pres->FillRef(buf,1);
			lbl_res[m] = buf;
		}
	}

	int row,col;
	
	int nrows = m_Grid->GetNumberRows();
	int ncols = m_Grid->GetNumberCols();

	if( nrows != m_nRows  || ncols != m_nCols )
	{
        
        if( ncols > 0) m_Grid->DeleteCols(0,ncols);
		if( nrows > 0) m_Grid->DeleteRows(0,nrows);
		m_Grid->InsertRows(0,m_nRows);
		m_Grid->InsertCols(0,m_nCols);
	}
	// fill rows/cols with text
    //  Titles:

	for (col = 0; col < m_Grid->GetNumberCols(); col++)
	{	 
		wxString str;
				
		if(display_mode == JUMNA_HLX_CRD)
		{
			if(col == 0) str = "Xdis";
			if(col == 1) str = "Ydis";
			if(col == 2) str = "Rise";
			if(col == 3) str = "Inc";
			if(col == 4) str = "Tip";
			if(col == 5) str = "Twist";
		}
		if(display_mode == BASE_BASE_CRD)
		{
			if(col == 0) str = "Shear";
			if(col == 1) str = "Stretch";
			if(col == 2) str = "Stagger";
			if(col == 3) str = "Buckle";
			if(col == 4) str = "Propel";
			if(col == 5) str = "Opening";
		}
		m_Grid->SetColLabelValue(col,str);
	}

	int ir;
	row = 0;

	for(k=0; k < 2; k++)
	{
		if(k >0 && display_mode == BASE_BASE_CRD) continue;
		if( display_mode == JUMNA_HLX_CRD && nucl_acid_mod->hel_crd.size() < 6*nr)
			continue;

		for(ir = 0; ir < nr; ir++)
		{
			wxString str;
			if(display_mode == JUMNA_HLX_CRD)
			{
				str.Printf("%s", lbl_res[ir + k*nr].c_str());
			}
			else if(display_mode == BASE_BASE_CRD)
			{
				str.Printf("%s-%s", lbl_res[ir].c_str(),lbl_res[ir+nr].c_str());
			}
			m_Grid->SetRowLabelValue(row,str);
			
			for (col = 0; col < m_Grid->GetNumberCols(); col++)
			{	 
				if(display_mode == JUMNA_HLX_CRD)
				{
					str.Printf("%7.3f",nucl_acid_mod->hel_crd[6*ir + col]);
				}
				else if(display_mode == BASE_BASE_CRD)
				{
					str.Printf("%7.3f",nucl_acid_mod->bs_bs_pars(col+1,ir+1));
				}

				int locked;
				
				if(display_mode == JUMNA_HLX_CRD)
				{
					locked = nucl_acid_mod->IsHelCoordLocked(row, col);

					if(locked)
					{
						wxColour clr_b(255, 255, 0);
						wxColour clr_f(255,   0, 0);
						m_Grid->SetCellBackgroundColour(row,col,clr_b);
						m_Grid->SetCellTextColour(row,col,clr_f);
					}
				
				}
				m_Grid->SetCellValue(row,col,str);
			}
			row++;
		}
	}
	m_Grid->AutoSizeColumns();
}

void NuclAcidDlgWX::OnLockSelCrd(wxCommandEvent& event)
{
	int col = m_Grid->GetGridCursorCol();
	int row = m_Grid->GetGridCursorRow();

	if( col >= 0 && row >= 0)
	{
		int ich   = nucl_acid_mod->ChainIdxOfRes(row);
		int ir_ch = nucl_acid_mod->IdxResInChain(row) + 1;			
		nucl_acid_mod->LockHelCoord(ich,ir_ch,col+1,1);
	}

	wxGridCellCoordsArray sel_cells = m_Grid->GetSelectedCells();
	 
//	int min_row = sel_cells.GetMinRow();
//	int min_col = sel_cells.GetMinCol();
//	int max_row = sel_cells.GetMaxRow();
//	int max_col = sel_cells.GetMaxCol();
    
//
//	for(row = min_row; row <= max_row; row++)
//	{
//		for( col = min_col; col <= max_col; col++)
//		{
	int n = sel_cells.Count();
	int i;
	for( i=0; i < n; i++)
	{
		row = sel_cells[i].GetRow();
		col = sel_cells[i].GetCol();

		if( row < 0 && col < 0) continue;
		int ich   = nucl_acid_mod->ChainIdxOfRes(row);
		int ir_ch = nucl_acid_mod->IdxResInChain(row) + 1;
			
		nucl_acid_mod->LockHelCoord(ich,ir_ch,col+1,1);
	}
	nucl_acid_mod->UpdateVarCoord();
	TransferDataToWindow();
}

void NuclAcidDlgWX::OnBaseBasePars(wxCommandEvent& event) 
{
	// TODO: Add your command handler code here
	display_mode = BASE_BASE_CRD;
	TransferDataToWindow();
}


void NuclAcidDlgWX::OnSetFromJumna(wxCommandEvent& event)
{
	if(!nucl_acid_mod) return;
	nucl_acid_mod->SetCoordsFromJumna();
	TransferDataToWindow();
}

void NuclAcidDlgWX::OnSetToJumna(wxCommandEvent& event)
{
	if(!nucl_acid_mod) return;
	nucl_acid_mod->SetIntCoordsToJumna();
	TransferDataFromWindow();
	
}

void NuclAcidDlgWX::OnUpdateXYZ(wxCommandEvent& event) 
{
	if(!nucl_acid_mod) return;
	nucl_acid_mod->SetIntCoordsToJumna();
	nucl_acid_mod->UpdateXYZ();
	TransferDataFromWindow();	
}


void NuclAcidDlgWX::OnGridEndEdit(wxGridEvent& event)
{
	PrintLog(" Enter NuclAcidDlgWX::OnGridEndEdit \n");
	
	int col = event.GetCol();
	int row = event.GetRow();	

    PrintLog("End Edit on row %d, col %d\n", row, col );

	wxString str_val = m_Grid->GetCellValue(row,col);

	double dval;
	bool bres;
	bres = str_val.ToDouble(&dval);

	if(bres)
	{
		nucl_acid_mod->hel_crd[6* row + col] = dval;
	}
	else
	{
		dval = nucl_acid_mod->hel_crd[6* row + col];
		str_val.Printf("%9.4f",dval);
		m_Grid->SetCellValue(row,col,str_val);
	}
}

int ScatterDlg::dlg_open = FALSE;

ScatterDlg::ScatterDlg(HaScatterMod* new_ptr_sc_mod, wxWindow* parent):
wxFrame( parent, -1, "Electron Scattering Module")
{
	ptr_sc_mod = new_ptr_sc_mod;
	dlg_open = TRUE;
	x1= y1 = z1 = 0.0;
	x2= y2 = z2 = 0.0;

	wxColour back_colour = wxSystemSettings::GetColour(wxSYS_COLOUR_BTNFACE);
 	SetBackgroundColour(back_colour);
	
	scatter_dlg(this,TRUE);

	OnInitDialog();
}

ScatterDlg::~ScatterDlg()
{
	dlg_open = FALSE;	
}

void ScatterDlg::OnInitDialog()
{
	DDX_Text_double(this,IDC_PTGRID_X1,x1,"%8.3f");
	DDX_Text_double(this,IDC_PTGRID_Y1,y1,"%8.3f");
	DDX_Text_double(this,IDC_PTGRID_Z1,z1,"%8.3f");
	DDX_Text_double(this,IDC_PTGRID_X2,x2,"%8.3f");
	DDX_Text_double(this,IDC_PTGRID_Y2,y2,"%8.3f");
	DDX_Text_double(this,IDC_PTGRID_Z2,z2,"%8.3f");
	DDX_Text_double(this,IDC_PTGRID_GAUSS_EXP,ptr_sc_mod->pt_exp,"%8.3f");
	DDX_Text_double(this,IDC_PTGRID_PSP_THRESH,ptr_sc_mod->psp_thresh_val,"%8.3f");

	TransferDataToWindow();
}

bool ScatterDlg::TransferDataToWindow()
{
	wxString str;
	
	wxTextCtrl* xmin_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_XMIN);
	wxTextCtrl* ymin_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_YMIN);
	wxTextCtrl* zmin_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_ZMIN);
	wxTextCtrl* xmax_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_XMAX);
	wxTextCtrl* ymax_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_YMAX);
	wxTextCtrl* zmax_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_ZMAX);
	wxTextCtrl* nx_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_NX);
	wxTextCtrl* ny_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_NY);
	wxTextCtrl* nz_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_NZ);

	double xmin_ang, ymin_ang, zmin_ang;
	double xmax_ang, ymax_ang, zmax_ang;
	int nx, ny, nz;

	nx= ptr_sc_mod->pot_grid->GetNx();
	ny= ptr_sc_mod->pot_grid->GetNy();
	nz= ptr_sc_mod->pot_grid->GetNz();

	str.Printf("%d",nx);
	nx_edt->SetValue(str);
	str.Printf("%d",ny);
	ny_edt->SetValue(str);
	str.Printf("%d",nz);
	nz_edt->SetValue(str);

	xmin_ang = ptr_sc_mod->pot_grid->GetXmin();
	ymin_ang = ptr_sc_mod->pot_grid->GetYmin();
	zmin_ang = ptr_sc_mod->pot_grid->GetZmin();
	xmax_ang = ptr_sc_mod->pot_grid->GetXmax();
	ymax_ang = ptr_sc_mod->pot_grid->GetYmax();
	zmax_ang = ptr_sc_mod->pot_grid->GetZmax();

	str.Printf("%8.3f",xmin_ang);
	xmin_edt->SetValue(str);
	str.Printf("%8.3f",ymin_ang);
	ymin_edt->SetValue(str);
	str.Printf("%8.3f",zmin_ang);
	zmin_edt->SetValue(str);

	str.Printf("%8.3f",xmax_ang);
	xmax_edt->SetValue(str);
	str.Printf("%8.3f",ymax_ang);
	ymax_edt->SetValue(str);
	str.Printf("%8.3f",zmax_ang);
	zmax_edt->SetValue(str);

	wxTextCtrl* psp_fname_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_PSP_FILE);
	psp_fname_edt->SetValue(ptr_sc_mod->psp_fname.c_str());

	return wxFrame::TransferDataToWindow();
}

bool ScatterDlg::TransferDataFromWindow()
{
	wxString str;

	wxTextCtrl* xmin_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_XMIN);
	wxTextCtrl* ymin_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_YMIN);
	wxTextCtrl* zmin_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_ZMIN);
	wxTextCtrl* xmax_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_XMAX);
	wxTextCtrl* ymax_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_YMAX);
	wxTextCtrl* zmax_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_ZMAX);
	wxTextCtrl* nx_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_NX);
	wxTextCtrl* ny_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_NY);
	wxTextCtrl* nz_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_NZ);

	double xmin_ang, ymin_ang, zmin_ang;
	double xmax_ang, ymax_ang, zmax_ang;
	long nx, ny, nz;

	str = xmin_edt->GetValue();
	str.ToDouble(&xmin_ang);
	str = ymin_edt->GetValue();
	str.ToDouble(&ymin_ang);
	str = zmin_edt->GetValue();
	str.ToDouble(&zmin_ang);
	str = xmax_edt->GetValue();
	str.ToDouble(&xmax_ang);
	str = ymax_edt->GetValue();
	str.ToDouble(&ymax_ang);
	str = zmax_edt->GetValue();
	str.ToDouble(&zmax_ang);

	ptr_sc_mod->pot_grid->SetGridCornersCoord(xmin_ang,ymin_ang,zmin_ang,xmax_ang,ymax_ang,zmax_ang);

	str = nx_edt->GetValue();
	str.ToLong(&nx);
	str = ny_edt->GetValue();
	str.ToLong(&ny);
	str = nz_edt->GetValue();
	str.ToLong(&nz);

	ptr_sc_mod->pot_grid->SetDimensions(nx,ny,nz);

	wxTextCtrl* psp_fname_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_PSP_FILE);
	str = psp_fname_edt->GetValue();
	ptr_sc_mod->psp_fname = str.c_str();

	return wxFrame::TransferDataFromWindow();
}


BEGIN_EVENT_TABLE(ScatterDlg, wxFrame)
    EVT_CLOSE(ScatterDlg::OnClose)
    EVT_BUTTON(IDC_PTGRID_CALC_PSP_XYZ, ScatterDlg::OnCalcPspXYZ)
	EVT_BUTTON(IDC_PTGRID_CALC_PSP_AO, ScatterDlg::OnCalcPsPAO)
	EVT_BUTTON(IDC_CALC_GRID_PSEUDOPOT, ScatterDlg::OnCalcPsPGrid)
	EVT_BUTTON(IDC_PTGRIG_SAVE_PSP_FILE, ScatterDlg::OnSavePSPFile)
	EVT_BUTTON(IDC_PTGRID_GRID_EIG, ScatterDlg::OnGridEigVal)
	EVT_BUTTON(IDC_PTGRID_SET_CORE_HAM, ScatterDlg::OnSetCoreHam)
	EVT_BUTTON(IDC_PTGRID_PSP_AO_XYZ, ScatterDlg::OnPSPGaussXYZ)
	EVT_BUTTON(IDC_PTGRID_KIN_ENE_XYZ, ScatterDlg::OnKinEneGaussXYZ)
END_EVENT_TABLE()


/////////////////////////////////////////////////////////////////////////////
// ScatterDlg message handlers
/////////////////////////////////////////////////////////////////////////////

void
ScatterDlg::OnClose(wxCloseEvent& event)
{
	event.Skip();
}


void ScatterDlg::OnCalcPspXYZ(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	if( ptr_sc_mod == NULL)
		return;

	double pp_val = ptr_sc_mod->GetPsP_xyz(x1,y1,z1,x2,y2,z2);

	PrintLog(" points: %8.3f %8.3f %8.3f  :   %8.3f %8.3f %8.3f \n",x1, y1, z1, x2, y2, z2);
	double dist = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );

	PrintLog(" dist = %8.3f  ppot = %12.6e \n\n", dist, pp_val);
}

void ScatterDlg::OnPSPGaussXYZ(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	if( ptr_sc_mod == NULL)
		return;

	double pp_val = ptr_sc_mod->GetPsP_Gauss_xyz(x1,y1,z1,x2,y2,z2);
	
    PrintLog(" points: %8.3f %8.3f %8.3f  :   %8.3f %8.3f %8.3f \n",
		      x1, y1, z1, x2, y2, z2);
	double dist = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );

	PrintLog(" dist = %8.3f  Gauss ppot = %12.6e \n\n", dist, pp_val);
}

void ScatterDlg::OnKinEneGaussXYZ(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	if( ptr_sc_mod == NULL)
		return;

	double kin_ene_val = ptr_sc_mod->GetKinEne_Gauss_xyz(x1,y1,z1,x2,y2,z2);
	
    PrintLog(" points: %8.3f %8.3f %8.3f  :   %8.3f %8.3f %8.3f \n",
		      x1, y1, z1, x2, y2, z2);
	double dist = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );

	PrintLog(" dist = %8.3f  Kinetic Energy Matrix Element between Gaussians  = %12.6e \n\n", dist, kin_ene_val);
}



void ScatterDlg::OnCalcPsPAO(wxCommandEvent& event) 
{
	if( ptr_sc_mod == NULL)
		return;
	wxBusyCursor wait;
	ptr_sc_mod->CalcPseudoPotAO();
}

void ScatterDlg::OnCalcPsPGrid(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	if( ptr_sc_mod == NULL)
		return;

	wxBusyCursor wait;
	ptr_sc_mod->CalcPseudoPotGrid();
}
		

void ScatterDlg::OnSavePSPFile(wxCommandEvent& event) 
{
	if(ptr_sc_mod == NULL)
		return;
	if(ptr_sc_mod->pot_grid == NULL)
		return;

	wxString str;
	wxTextCtrl* psp_fname_edt = (wxTextCtrl*) FindWindow(IDC_PTGRID_PSP_FILE);
	str = psp_fname_edt->GetValue();
	ptr_sc_mod->pot_grid->SaveField( str.ToStdString() );
}

void ScatterDlg::OnGridEigVal(wxCommandEvent& event) 
{
	wxBusyCursor wait;
	ptr_sc_mod->FindGridHamEigVec(3);
}

void ScatterDlg::OnSetCoreHam(wxCommandEvent& event) 
{
	TransferDataFromWindow();
	ptr_sc_mod->SetCoreHamOnGrid();
}

/////////////////////////////////////////////////////////////////////////////
// wxAtomEdit


std::list<wxAtomEdit*>  wxAtomEdit::active_controls;

wxAtomEdit::wxAtomEdit(wxWindow* parent, wxWindowID id, const wxString& value,
		const wxPoint& pos,const wxSize& size, long style, const wxValidator& validator,
		const wxString& name):
wxTextCtrl(parent,id,value, pos, size,style,validator,name)
{
//	active_controls.push_back(this);
	pick_mode = PICK_ATOM;
}

wxAtomEdit::~wxAtomEdit()
{
//	list<wxAtomEdit*>::iterator citr;
//	for(citr = active_controls.begin(); citr != active_controls.end(); citr++ )
//	{
//		if( (*citr) == this )
//			active_controls.erase(citr);
//	}
}

BEGIN_EVENT_TABLE(wxAtomEdit, wxTextCtrl)
   EVT_SET_FOCUS(wxAtomEdit::OnSetFocus)
END_EVENT_TABLE()

void wxAtomEdit::OnSetFocus(wxFocusEvent& event)
{
//	PrintLog("wxAtomEdit::OnSetFocus() \n");
//	wxWindow* parent = this->GetParent();
//	parent->ProcessEvent(event);
	event.Skip();	
}


int wxAtomEdit::BroadCastPickedAtom(HaAtom* PickAtom)
{
	std::list<wxAtomEdit*>::iterator citr;
	for(citr = active_controls.begin(); citr != active_controls.end(); citr++ )
	{
		(*citr)->OnAtomPicked(PickAtom);
	}
	return TRUE;
	
}

int wxAtomEdit::OnAtomPicked(HaAtom* PickAtom)
{
	if(PickAtom == NULL)
		return FALSE;

	if(wxWindow::FindFocus() == (wxWindow*) this)
	{
		char buf[256];
		if(pick_mode == PICK_ATOM)
		{
			PickAtom->FillRef(buf);
		}
		else if(pick_mode == PICK_RESIDUE)
		{
			HaResidue* pres = PickAtom->GetHostRes();
			pres->FillRef(buf);
		}
		else if( pick_mode == PICK_MOLECULE)
		{
			HaMolecule* pMol = PickAtom->GetHostMol();
			pMol->FillRef(buf);
		}
		SetValue(buf);
	}
	return TRUE;
}

int HaMolView::BroadcastCurrAtom()
//! This function send message IDU_ATOM_SELECT to Dialog classes that are able 
//! to prosess this info and call 
//! AtomEdit::BroadCastPickedAtom(PkAtom) - for atom select edit boxes
{
	MolSet* mset= GetMolSet();

	if( EditGroupsDlg::dlg_open)
	{
	    wxCommandEvent at_event( wxEVT_COMMAND_BUTTON_CLICKED, IDU_ATOM_PICK );
		at_event.SetEventObject( EditGroupsDlg::active_dlg_ptr );
	//	EditGroupsDlg::active_dlg_ptr->ProcessEvent(at_event);
		EditGroupsDlg::active_dlg_ptr->ProcessWindowEvent( at_event );
	}

	if( EditGeomDlgWX::dlg_open)
	{
	    wxCommandEvent at_event(wxEVT_COMMAND_BUTTON_CLICKED,IDU_ATOM_PICK);
		at_event.SetEventObject( EditGeomDlgWX::active_dlg_ptr );
//		EditGeomDlgWX::active_dlg_ptr->ProcessEvent(at_event);:
		EditGeomDlgWX::active_dlg_ptr->ProcessWindowEvent( at_event );
	}

	if (EditMutMapDlg::dlg_open)
	{
		wxCommandEvent at_event(wxEVT_COMMAND_BUTTON_CLICKED, IDU_ATOM_PICK);
		at_event.SetEventObject(EditMutMapDlg::active_dlg_ptr);
		//		EditGeomDlgWX::active_dlg_ptr->ProcessEvent(at_event);:
		EditMutMapDlg::active_dlg_ptr->ProcessWindowEvent(at_event);
	}
	
//	wxAtomEdit::BroadCastPickedAtom(PkAtom);

	return True;
}

void
HaMainFrameWX::OnShowResDb(wxCommandEvent &event)
{
	HaResDB* p_res_db = HaResDB::GetDefaultResDB();
	MolViewWX* mol_vew_wx = CreateMolView(p_res_db);

	if(ResDBDlg::dlg_open) return; 
	ResDBDlg* ptr_res_db_dlg = new ResDBDlg(this);
	ptr_res_db_dlg->Show(TRUE);
}
