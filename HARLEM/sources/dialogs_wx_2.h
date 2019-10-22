/*! \file dialogs_wx_2.h

    wxWidgets Dialogs 2 in HARLEM 
 
    \author Igor Kurnikov   
    \date 2005-
*/

#if !defined(DIALOGS_WX_2_H)
#define DIALOGS_WX_2_H
#ifndef SWIG
#include <wx/spinctrl.h>

class HaInterMolMod;
class HaEmpiricalMod;
class InterMolMCSimulator;
class TrajIOAgent;

class InterMolDlgWX : public wxFrame
{
public:
	InterMolDlgWX(HaInterMolMod* ptr_immod_new, wxWindow *parent);   
    virtual ~InterMolDlgWX();

    static int dlg_open;

	virtual bool TransferDataToWindow();
	virtual bool TransferDataFromWindow();

	void OnInitDialog();

protected:

	HaInterMolMod* ptr_im_mod;
    HaEmpiricalMod* ptr_emp_mod;
	InterMolMCSimulator* p_mc_sim;
	TrajIOAgent* p_io_ag;

	void OnIntermolInit(wxCommandEvent& event);
	void OnIntermolElStat(wxCommandEvent& event);
	void OnIntermolTotEne(wxCommandEvent& event);
	void OnIntermolSetIntCoord(wxCommandEvent& event);
	void OnIntermolCompIntCoord(wxCommandEvent& event);
	void OnMcDockRun(wxCommandEvent& event);
	void OnMcDockPlaybackTraj(wxCommandEvent& event);
	void OnMcDockPause(wxCommandEvent& event);
	void OnMcDockResume(wxCommandEvent& event);
	void OnMcDockStop(wxCommandEvent& event);
	void OnClose(wxCloseEvent &event);
protected:
	DECLARE_EVENT_TABLE()
};

class PlaneViewOfHaField3D;
class MolSet;
//wxFieldPlaneView

class wxFieldPlaneView : public wxDialog
{
	public:
		wxFieldPlaneView(PlaneViewOfHaField3D* PlView, MolSet* _MolSet, wxWindow* parent,wxWindowID id = -1, const wxString& title = "wxOGLPlaneViewerDlg", const wxPoint& pos = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxDEFAULT_DIALOG_STYLE, const wxString& name = "wxFieldPlaneView");
		~wxFieldPlaneView();

		bool TransferDataToWindow();
		bool TransferDataFromWindow();
		int OwnerOfView;
	protected:
		void OnScroll( wxScrollEvent& event );
		void OnLevelChange( wxSpinEvent& event );
		void OnEnter(wxCommandEvent& event);
		void OnChoosePlane( wxCommandEvent& event );
		void OnCopy1D1stVar( wxCommandEvent& event );
		void OnCopy1D2ndVar( wxCommandEvent& event );
		void OnCopy2D( wxCommandEvent& event );
		void OnHideZeroValue( wxCommandEvent& event );
		void OnHide( wxCommandEvent& event );
		void OnClose(wxCloseEvent& event );
		
	private:
		MolSet* m_MolSet;
		PlaneViewOfHaField3D *m_PlaneView;
	protected:
		DECLARE_EVENT_TABLE()
};
#endif
void CreatewxFieldPlaneView(PlaneViewOfHaField3D* PlView, MolSet* _MolSet,const char *title,int OwnerOfView);


#endif  // end !defined(DIALOGS_WX_2_H)


