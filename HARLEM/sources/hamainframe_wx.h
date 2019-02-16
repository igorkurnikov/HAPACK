/*! \file hamainframe_wx.h

    wxWidgets MainFrame and MolView classes for HARLEM application
    \author Igor Kurnikov  
    \date 2003-
*/

#if !defined(HAMAINFRAME_WX_H)
#define HAMAINFRAME_WX_H

class MolViewWX;
class HaMolView;
class HaMolSet;
class HaDialogBar;
class mVueFrame;//mikola 09/19/2006

class wxSashLayoutWindow;

class HaMainFrameWX : public wxMDIParentFrame
{
  public:

    HaMainFrameWX();
 
	void OnExecuteCommand( wxCommandEvent &event );
	
	void OnHelp( wxCommandEvent &event );
    void OnAbout( wxCommandEvent &event );
    void OnTestGraph1( wxCommandEvent &event );
	void OnTestPython1( wxCommandEvent &event );
	void OnTestAvgPopFunc( wxCommandEvent &event );
	void OnSaveAmoebaTopFile1( wxCommandEvent &event );
	void OnSaveAmoebaTopFile2( wxCommandEvent &event );
	void OnSaveResFFMort( wxCommandEvent &event );
	void OnTestFFT1( wxCommandEvent &event );
    void OnTestMin1( wxCommandEvent &event );
    void OnTestMap( wxCommandEvent &event );
    void OnCalcTMatr1( wxCommandEvent &event );
    void OnTestQCMod1( wxCommandEvent &event );
    void OnTestOper2( wxCommandEvent &event );
    void OnTestOper1( wxCommandEvent &event );
    void OnDumpOverlap2( wxCommandEvent &event );
    void OnDumpOverlap( wxCommandEvent &event );
    void OnDumpGaussBCommon( wxCommandEvent &event );
    void OnDumpMolInfo( wxCommandEvent &event );
    void OnClose( wxCloseEvent& event);
// File Menu
    void OnFileNew   (wxCommandEvent &event);
    void OnFileOpen  (wxCommandEvent &event);
    void OnMolSave   (wxCommandEvent &event);
    void OnMolSaveAs (wxCommandEvent &event);
    void DoLoadScriptDialog(wxCommandEvent &event);
    void DoRedirectIODialog (wxCommandEvent &event);
    void SaveOutputFile (wxCommandEvent &event);
    void OnSaveImageClipboard (wxCommandEvent &event);
    void DoCompAccountsDialog (wxCommandEvent &event);
    void OnInfo (wxCommandEvent &event);
    void OnFileClose (wxCommandEvent &event);
    void OnPrint (wxCommandEvent &event);
    void OnSetup (wxCommandEvent &event);
    void OnExit(wxCommandEvent &event);
    void OnPyMod(wxCommandEvent &event);//<mikola Jul 19, 2006
// Edit Menu
    void OnSelectAll ( wxCommandEvent &event ); 
    void DoAtomParamsDialog ( wxCommandEvent &event );  
    void DoResidueParamsDialog ( wxCommandEvent &event );
    void DoMolSetParamDialog ( wxCommandEvent &event );
	void DoAtomPropColorDialog ( wxCommandEvent &event );
    void DoEditGroupsDialog ( wxCommandEvent &event );
	void DoCrdSnapshotDialog ( wxCommandEvent &event );
    void DoEditGeomDialog ( wxCommandEvent &event );
    void OnFindHbonds ( wxCommandEvent &event );
    void DoSolvateDialog ( wxCommandEvent &event );
	void OnWrapIntoUnitCell( wxCommandEvent &event );
	void OnCenterSolute( wxCommandEvent &event );
	void OnCenterMolInPBox( wxCommandEvent &event );
    void OnAddMissingAtoms ( wxCommandEvent &event );
    void OnAddPolarHydrogens ( wxCommandEvent &event );
    void OnAddHHybrid ( wxCommandEvent &event );
    void OnAddHydrogens ( wxCommandEvent &event );
    void OnDelSelAtoms ( wxCommandEvent &event );
	void OnDelOvlpMols ( wxCommandEvent &event );
    void DoEditFragmDialog ( wxCommandEvent &event );
    void DoBuildFilmDialog ( wxCommandEvent &event );
    void OnClearPicked ( wxCommandEvent &event );
	void OnSelAtomsInBoundBox ( wxCommandEvent &event );
	void OnRevertAtomSelection ( wxCommandEvent &event );
    void DoNuclAcidDialog ( wxCommandEvent &event );
	void OnDescribeSecStruct ( wxCommandEvent &event );
	void OnPrintHBonds ( wxCommandEvent &event );
	void OnSetAlphaHelix ( wxCommandEvent &event );
    void OnShowResDb ( wxCommandEvent &event );

// Display Menu
//  Display->Display Mode menu
    void OnWireFrame ( wxCommandEvent &event );
    void OnBackBone  ( wxCommandEvent &event );
    void OnSticks    ( wxCommandEvent &event );
    void OnSpheres   ( wxCommandEvent &event );
    void OnBallStick ( wxCommandEvent &event );
    void OnRibbons   ( wxCommandEvent &event );
    void OnStrands   ( wxCommandEvent &event );
    void OnCartoons  ( wxCommandEvent &event );
//  Display->Colours menu
    void OnMono       (wxCommandEvent &event );
    void OnCPK        ( wxCommandEvent &event );
    void OnShapely    ( wxCommandEvent &event );
    void OnColGroups  ( wxCommandEvent &event );
    void OnColResidues( wxCommandEvent &event );
    void OnColChain   ( wxCommandEvent &event );
    void OnColTemper  ( wxCommandEvent &event );
    void OnStruct     ( wxCommandEvent &event );
//  Display->Options menu
    void OnSlab       ( wxCommandEvent &event );
    void OnHydrogen   ( wxCommandEvent &event );
    void OnHetero     ( wxCommandEvent &event );
    void OnSpecular   ( wxCommandEvent &event );
    void OnStereo     ( wxCommandEvent &event );
    void OnLabels     ( wxCommandEvent &event );
	void OnDisplayPBox ( wxCommandEvent &event );
//
    void DoViewParamDlg   ( wxCommandEvent &event );
    void DoObject3DDialog ( wxCommandEvent &event );
// Display->Labels menu
    void OnLabelsAtomId     ( wxCommandEvent &event );
    void OnLabelsAtomNames  ( wxCommandEvent &event );
    void OnLabelsAtomSeqNum ( wxCommandEvent &event );
	void OnLabelsAtomFFSymb ( wxCommandEvent &event );
    void OnLabelsOff        ( wxCommandEvent &event );
// Display->Windows menu
    void OnWindowsCascade   ( wxCommandEvent &event );
    void OnWindowsTileVert  ( wxCommandEvent &event );
    void OnCenterViewSel    ( wxCommandEvent &event );  
// Measure menu
    void OnShowAtomID    ( wxCommandEvent &event );
    void OnMeasureDist   ( wxCommandEvent &event );
    void OnMeasureAngle  ( wxCommandEvent &event );
    void OnMeasureDihed  ( wxCommandEvent &event );
// ET menu
    void OnEditRedox     ( wxCommandEvent &event );
    void DoPathwaysDialog( wxCommandEvent &event );
    void DoETEffHamDialog( wxCommandEvent &event );
    void OnResetETModule ( wxCommandEvent &event );
    void DoScatterDialog ( wxCommandEvent &event );
// QChem menu
    void DoQChemParamDialog ( wxCommandEvent &event );
    void DoLoadQCDatDialog ( wxCommandEvent &event );
    void DoWaveFunAnalDialog ( wxCommandEvent &event );
    void OnCalcPolarGcontr ( wxCommandEvent &event );
    void OnSaveGrpOperMat ( wxCommandEvent &event );
    void OnCalcPolarContrF ( wxCommandEvent &event );
    void OnReadPolarContr ( wxCommandEvent &event );
    void OnCalcBetaContr2idx ( wxCommandEvent &event );
    void OnReadBetaContr2idx ( wxCommandEvent &event );
// MMech menu
    void DoMolMechDialog    ( wxCommandEvent &event );
    void DoInterMolDialog   ( wxCommandEvent &event );
	void DoProtRedoxDialog  ( wxCommandEvent &event );
    void DoContElectrDialog ( wxCommandEvent &event );
    void DoPNPDialog        ( wxCommandEvent &event );
    void DoContElectrPNPDialog( wxCommandEvent &event );
    void DoSmlPNPDialog(wxCommandEvent &event);
    void DoRFDialog(wxCommandEvent &event);
	void DoFlexModDialog(wxCommandEvent& event);
	void DoEDDialog(wxCommandEvent& event);
    //Test Menu
//
    void OnMove(wxMoveEvent& event);
    void OnSize(wxSizeEvent& event);
//
    void OnLoadManual(wxCommandEvent& event);
    void OnMolConnect(wxCommandEvent& event);
    void OnWorldConnect(wxCommandEvent& event);
		//Test Menu->PNPS
		void OnOpenVectorFieldInHaMolView(wxCommandEvent& event);
		void OnOpenDielConstFromNindexInHaMolView(wxCommandEvent& event);
		void OnOpenDiffConstFromNindexInHaMolView(wxCommandEvent& event);
		void OnOpenQstFromNindexInHaMolView(wxCommandEvent& event);
    //AddOn options

    MolViewWX* CreateMolView(HaMolSet* pmset);

	wxSashLayoutWindow* sash_win;

private:
//    DECLARE_CLASS(HaMainFrameWX);
    mVueFrame *m_TheFrame;//mikola 09/19/2006
    DECLARE_EVENT_TABLE();
};

HaMainFrameWX* GetHaMainFrameWX();

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

class MolViewWX: public wxWindow
{
public:    
    MolViewWX( HaMolSet* pmset_new, MolViewFrame *frame, const wxPoint& pos, 
               const wxSize& size, long style);
	virtual ~MolViewWX();
    virtual void OnDraw(wxDC& dc);
    void OnMouseEvent(wxMouseEvent& event);
    void OnEraseBackground(wxEraseEvent& event);
    void OnPaint(wxPaintEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnClose(wxCloseEvent& event);
    void MouseMove( wxMouseEvent& event, int dx, int dy );

    MolViewFrame* mol_frame;
    HaMolView* mol_view;
    HaMolSet*  pmset;

    int PointX; //!< X screen coordinate of the point of mouse click or release
    int PointY; //!< Y screen coordinate of the point of mouse click or release 
    int InitX;  //!< reference X screen coordinate of the point of mouse click in mouse motion tracking 
    int InitY;  //!< reference Y screen coordinate of the point of mouse click in mouse motion tracking
   
    int HeldButton;
  private:
    DECLARE_EVENT_TABLE();
};


#endif // !defined(HAMAINFRAME_WX_H)


