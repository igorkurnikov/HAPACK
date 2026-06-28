/*! \file hamainframe_wx.h

    wxWidgets MainFrame and MolView classes for HARLEM application
    \author Igor Kurnikov
    \date 2003-
*/

#if !defined(HAMAINFRAME_WX_H)
#define HAMAINFRAME_WX_H

#include "hamolview.h"
#include "main_menu_base.h"

class MolViewWX;
class HaMolView;
class MolSet;
class HaDialogBar;
class mVueFrame;//mikola 09/19/2006

class wxSashLayoutWindow;

class HaMainFrameWX : public HaMainMenuBase
{
  public:

    HaMainFrameWX();

	void OnExecuteCommand( wxCommandEvent &event );

    void OnClose( wxCloseEvent& event);
// File Menu
    void OnFileNew   (wxCommandEvent &event) override;
    void OnFileOpen  (wxCommandEvent &event) override;
    void OnMolSave   (wxCommandEvent &event) override;
    void OnMolSaveAs (wxCommandEvent &event) override;
    void DoLoadScriptDialog(wxCommandEvent &event) override;
    void DoRedirectIODialog (wxCommandEvent &event) override;
    void SaveOutputFile (wxCommandEvent &event) override;
    void OnSaveImageClipboard (wxCommandEvent &event) override;
    void DoCompAccountsDialog (wxCommandEvent &event) override;
    void OnInfo (wxCommandEvent &event) override;
    void OnFileClose (wxCommandEvent &event) override;
    void OnPrint (wxCommandEvent &event) override;
    void OnSetup (wxCommandEvent &event) override;
    void OnExit(wxCommandEvent &event) override;
    void OnPyMod(wxCommandEvent &event);//<mikola Jul 19, 2006
// Edit Menu
    void OnSelectAll ( wxCommandEvent &event ) override;
    void DoAtomParamsDialog ( wxCommandEvent &event ) override;
    void DoResidueParamsDialog ( wxCommandEvent &event ) override;
    void DoMolSetParamDialog ( wxCommandEvent &event ) override;
	void DoAtomPropColorDialog ( wxCommandEvent &event ) override;
    void DoEditGroupsDialog ( wxCommandEvent &event ) override;
	void DoCrdSnapshotDialog ( wxCommandEvent &event ) override;
    void DoEditGeomDialog ( wxCommandEvent &event ) override;
    void DoEditMutationMapDialog(wxCommandEvent& event) override;
    void OnFindHbonds ( wxCommandEvent &event ) override;
    void DoSolvateDialog ( wxCommandEvent &event ) override;
	void OnWrapIntoUnitCell( wxCommandEvent &event ) override;
	void OnCenterSolute( wxCommandEvent &event ) override;
	void OnCenterMolInPBox( wxCommandEvent &event ) override;
    void OnAddMissingAtoms ( wxCommandEvent &event ) override;
    void OnAddPolarHydrogens ( wxCommandEvent &event ) override;
    void OnAddHHybrid ( wxCommandEvent &event ) override;
    void OnAddHydrogens ( wxCommandEvent &event ) override;
    void OnDelSelAtoms ( wxCommandEvent &event ) override;
	void OnDelOvlpMols ( wxCommandEvent &event ) override;
    void DoEditFragmDialog ( wxCommandEvent &event ) override;
    void DoBuildFilmDialog ( wxCommandEvent &event ) override;
    void OnClearPicked ( wxCommandEvent &event ) override;
	void OnSelAtomsInBoundBox ( wxCommandEvent &event ) override;
	void OnRevertAtomSelection ( wxCommandEvent &event ) override;
    void OnExpandAtomSelectionBonded(wxCommandEvent& event) override;
    void DoNuclAcidDialog ( wxCommandEvent &event ) override;
	void OnDescribeSecStruct ( wxCommandEvent &event ) override;
	void OnPrintHBonds ( wxCommandEvent &event ) override;
	void OnSetAlphaHelix ( wxCommandEvent &event ) override;
    void OnShowResDb ( wxCommandEvent &event ) override;

// Display Menu
//  Display->Display Mode menu
    void OnWireFrame ( wxCommandEvent &event ) override;
    void OnBackBone  ( wxCommandEvent &event );
    void OnSticks    ( wxCommandEvent &event ) override;
    void OnSpheres   ( wxCommandEvent &event ) override;
    void OnBallStick ( wxCommandEvent &event ) override;
    void OnRibbons   ( wxCommandEvent &event ) override;
    void OnStrands   ( wxCommandEvent &event ) override;
    void OnCartoons  ( wxCommandEvent &event ) override;
//  Display->Colours menu
    void OnMono       (wxCommandEvent &event ) override;
    void OnCPK        ( wxCommandEvent &event ) override;
    void OnShapely    ( wxCommandEvent &event ) override;
    void OnColGroups  ( wxCommandEvent &event ) override;
    void OnColResidues( wxCommandEvent &event ) override;
    void OnColChain   ( wxCommandEvent &event ) override;
    void OnColTemper  ( wxCommandEvent &event ) override;
    void OnStruct     ( wxCommandEvent &event ) override;
//  Display->Options menu
    void OnSlab       ( wxCommandEvent &event ) override;
    void OnHydrogen   ( wxCommandEvent &event ) override;
    void OnHetero     ( wxCommandEvent &event ) override;
    void OnSpecular   ( wxCommandEvent &event ) override;
    void OnStereo     ( wxCommandEvent &event ) override;
    void OnLabels     ( wxCommandEvent &event ) override;
	void OnDisplayPBox ( wxCommandEvent &event ) override;
//
    void DoViewParamDlg   ( wxCommandEvent &event ) override;
    void DoObject3DDialog ( wxCommandEvent &event ) override;
// Display->Labels menu
    void OnLabelsAtomId     ( wxCommandEvent &event ) override;
    void OnLabelsAtomNames  ( wxCommandEvent &event ) override;
    void OnLabelsAtomSeqNum ( wxCommandEvent &event ) override;
	void OnLabelsAtomFFSymb ( wxCommandEvent &event ) override;
    void OnLabelsAtomCharges( wxCommandEvent &event ) override;
    void OnLabelsOff        ( wxCommandEvent &event ) override;
// Display->Windows menu
    void OnWindowsCascade   ( wxCommandEvent &event ) override;
    void OnWindowsTileVert  ( wxCommandEvent &event ) override;
    void OnCenterViewSel    ( wxCommandEvent &event ) override;
// Measure menu
    void OnShowAtomID    ( wxCommandEvent &event ) override;
    void OnMeasureDist   ( wxCommandEvent &event ) override;
    void OnMeasureAngle  ( wxCommandEvent &event ) override;
    void OnMeasureDihed  ( wxCommandEvent &event ) override;
    void OnAddDistMonitor( wxCommandEvent& event ) override;
    void OnShowAtomLabel ( wxCommandEvent& event ) override;
// ET menu
    void OnEditRedox     ( wxCommandEvent &event ) override;
    void DoPathwaysDialog( wxCommandEvent &event ) override;
    void DoETEffHamDialog( wxCommandEvent &event ) override;
    void OnResetETModule ( wxCommandEvent &event ) override;
    void DoScatterDialog ( wxCommandEvent &event ) override;
// QChem menu
    void DoQChemParamDialog ( wxCommandEvent &event ) override;
    void DoLoadQCDatDialog ( wxCommandEvent &event ) override;
    void DoWaveFunAnalDialog ( wxCommandEvent &event ) override;
    void OnCalcPolarGcontr ( wxCommandEvent &event ) override;
    void OnSaveGrpOperMat ( wxCommandEvent &event ) override;
    void OnCalcPolarContrF ( wxCommandEvent &event ) override;
    void OnReadPolarContr ( wxCommandEvent &event ) override;
    void OnCalcBetaContr2idx ( wxCommandEvent &event ) override;
    void OnReadBetaContr2idx ( wxCommandEvent &event ) override;
// MMech menu
    void DoMolMechDialog    ( wxCommandEvent &event ) override;
    void DoInterMolDialog   ( wxCommandEvent &event ) override;
	void DoProtRedoxDialog  ( wxCommandEvent &event ) override;
    void DoContElectrDialog ( wxCommandEvent &event ) override;
    void DoPNPDialog        ( wxCommandEvent &event ) override;
    void DoContElectrPNPDialog( wxCommandEvent &event );
    void DoSmlPNPDialog(wxCommandEvent &event) override;
    void DoRFDialog(wxCommandEvent &event) override;
	void DoFlexModDialog(wxCommandEvent& event) override;
	void DoEDDialog(wxCommandEvent& event) override;
    //Test Menu
    void OnTestOper1( wxCommandEvent &event ) override;
    void OnTestOper2( wxCommandEvent &event ) override;
    void OnTestQCMod1( wxCommandEvent &event ) override;
    void OnCalcTMatr1( wxCommandEvent &event ) override;
    void OnTestMap( wxCommandEvent &event ) override;
    void OnTestMin1( wxCommandEvent &event ) override;
    void OnTestGraph1( wxCommandEvent &event ) override;
	void OnTestPython1( wxCommandEvent &event ) override;
	void OnTestAvgPopFunc( wxCommandEvent &event ) override;
	void OnSaveAmoebaTopFile1( wxCommandEvent &event ) override;
	void OnSaveAmoebaTopFile2( wxCommandEvent &event ) override;
	void OnSaveResFFMort( wxCommandEvent &event );
	void OnTestFFT1( wxCommandEvent &event ) override;
    void OnDumpOverlap2( wxCommandEvent &event ) override;
    void OnDumpOverlap( wxCommandEvent &event ) override;
    void OnDumpGaussBCommon( wxCommandEvent &event ) override;
    void OnDumpMolInfo( wxCommandEvent &event ) override;
//
    void OnMove(wxMoveEvent& event);
    void OnSize(wxSizeEvent& event);
//
    void OnLoadManual(wxCommandEvent& event) override;
    void OnMolConnect(wxCommandEvent& event) override;
    void OnWorldConnect(wxCommandEvent& event) override;
		//Test Menu->PNPS
		void OnOpenVectorFieldInHaMolView(wxCommandEvent& event) override;
		void OnOpenDielConstFromNindexInHaMolView(wxCommandEvent& event) override;
		void OnOpenDiffConstFromNindexInHaMolView(wxCommandEvent& event) override;
		void OnOpenQstFromNindexInHaMolView(wxCommandEvent& event) override;
    //Help menu
    void OnHelp( wxCommandEvent &event ) override;
    void OnAbout( wxCommandEvent &event ) override;
    //AddOn options

    std::shared_ptr<BaseMolView> CreateMolView(MolSet* pmset);

    int RedirectIOLogWindow(); //!< Create Log Window and redirect to it STDOUT and STDERR

	wxSashLayoutWindow* sash_win;

private:
//    DECLARE_CLASS(HaMainFrameWX);
    // mVueFrame *m_TheFrame;//mikola 09/19/2006
    DECLARE_EVENT_TABLE();
};

HaMainFrameWX* GetHaMainFrameWX();

// put here wxEvent ID used to commpunicate to MPI threads of HARLEM - have to switch to non wx messages
extern const wxEventType wxEVT_HARLEM_APP;
const int HARLEM_APP_EXIT_ID = 20001;

//! HaMainFrameWX starter for starting GUI from python
MSET_API void StartHaMainFrameWX();

class MolViewWX;

class MolViewFrame: public wxMDIChildFrame
{
public:
    MolViewFrame(wxMDIParentFrame* parent, wxString& mset_name);

    void OnClose(wxCloseEvent& event);
	void OnActivate(wxActivateEvent& event);

    MolViewWX* mol_view_wx;
private:
  DECLARE_EVENT_TABLE();
};

class MolViewWX: public wxWindow, public BaseMolView
{
public:
    MolViewWX( MolSet* pmset_new, MolViewFrame *frame, const wxPoint& pos,
               const wxSize& size, long style);
	virtual ~MolViewWX();
    virtual void OnDraw(wxDC& dc);

    int  SetWXImage(wxImage& wx_image); //!< Construct wxImage with from the Molecular View
    void OnMouseEvent(wxMouseEvent& event);
    void OnEraseBackground(wxEraseEvent& event);
    void OnPaint(wxPaintEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnClose(wxCloseEvent& event);
    void MouseMove( wxMouseEvent& event, int dx, int dy );

    void UpdateThisView(int lHint = 0) override;   //!< Update the View
    void SetTitle(std::string title_str) override; //!< Set View Title
    void BroadcastCurrAtom() override;             //!< Broadcast change of the selected atom

    int WriteBMPFile(std::string name);
    int WriteJPEGFile(std::string name);
    int WriteTIFFFile(std::string name);
    int WritePNGFile(std::string name);
    int WritePCXFile(std::string name);

    MolViewFrame* mol_frame;
    HaMolView* mol_view;
    MolSet*  pmset;

    int PointX; //!< X screen coordinate of the point of mouse click or release
    int PointY; //!< Y screen coordinate of the point of mouse click or release
    int InitX;  //!< reference X screen coordinate of the point of mouse click in mouse motion tracking
    int InitY;  //!< reference Y screen coordinate of the point of mouse click in mouse motion tracking

    int HeldButton;
  private:
    DECLARE_EVENT_TABLE();
};


#endif // !defined(HAMAINFRAME_WX_H)
