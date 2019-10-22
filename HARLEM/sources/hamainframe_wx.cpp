/*! file hamainframe_wx.cpp

    wxWidgets MainFrame and MolView classes for HARLEM application              

    \author Igor Kurnikov 
    \date 2003-

*/

#define HAMAINFRAME_WX_CPP

#include <mpi.h>

#if !defined(HARLEM_PYTHON_NO)
#include <Python.h>
#endif

#include "vec3d.h"

#include "wx/wxprec.h" 
#include <wx/app.h>
#include <wx/docmdi.h>
#include <wx/image.h>
#include <wx/tooltip.h>
#include <wx/filedlg.h>
#include "wx/clipbrd.h"
#include "wx/printdlg.h"

#include "wx/wx.h"
#include "wx/laywin.h"
#include "wx/sizer.h"
#include "wx/filename.h"
#include <wx/splitter.h>
#include "ha_wx_res_wdr.h"
#include "haconst.h"

#include "hamolview.h"
#include "hampi.h"
#include "harlemapp.h"
#include "hamainframe_wx.h"
#include "hamolecule.h"
#include "moleditor.h"
#include "harlemapp.h"
#include "stmmod.h"
#include "etcoupl.h"
#include "protonredox.h"
#include "hatests.h"
#include "hasurface.h"

#include "wxMolFlex.h"
#include "wxMolED.h"

#include "dialogs_wx_1.h"
#include "dialogs_wx_2.h"
#include "mm_dialogs_wx.h"
#include "qc_dialogs_wx.h"
#include "wx_prot_redox_dlg.h"
#include "harlemapp_wx.h"
#include "hamolmech.h"
#include "mm_driver_amber.h"

#include "contworld.h"//<mikola July 30, 2007
#include "mapio.h"//<mikola Jan 8, 2008
#include "hasvnrev.h"
#if defined(_MSC_VER)
#include <windows.h>
#endif

#define IsClose(u,v) (((u)>=(v)-1) && ((u)<=(v)+1))


HaMainFrameWX *m_HaMainFrameWX = nullptr;
HaMainFrameWX* GetHaMainFrameWX() { return m_HaMainFrameWX; }

void StartHaMainFrameWX()
{
	if (wxTheApp == nullptr) {
		PrintLog("Can not start HaMainFrameWX wxApp is not started\n");
	}
	HaMainFrameWX *m_mainFrame = new HaMainFrameWX();

	m_mainFrame->CreateToolBar();
	MainToolBarFunc(m_mainFrame->GetToolBar());

	//   m_mainFrame->Centre(wxBOTH);
	m_mainFrame->Show(TRUE);

	wxTheApp->SetTopWindow(m_mainFrame);

	MolSet* pmset = GetCurMolSet();
	if (pmset != NULL)
	{
		pmset->canvas_wx = m_mainFrame->CreateMolView(pmset);
		pmset->canvas_wx->mol_view->InitialTransform();
		pmset->canvas_wx->mol_view->DefaultRepresentation();
		pmset->RefreshAllViews();
	}
}

// For drawing lines in a canvas
long xpos = -1;
long ypos = -1;

BEGIN_EVENT_TABLE(MolViewWX, wxWindow)
    EVT_MOUSE_EVENTS(MolViewWX::OnMouseEvent)
    EVT_ERASE_BACKGROUND(MolViewWX::OnEraseBackground)
    EVT_PAINT(MolViewWX::OnPaint)
    EVT_SIZE(MolViewWX::OnSize)
    EVT_CLOSE(MolViewWX::OnClose)
END_EVENT_TABLE()

MolViewWX::MolViewWX(MolSet* pmset_new, MolViewFrame *frame, const wxPoint& pos, const wxSize& size, long style):
 wxWindow(frame, -1, pos, size, style)
{
  mol_frame = frame;

  PointX = 0;
  PointY = 0;
  InitX = 0;
  InitY = 0;
  pmset = pmset_new;

  mol_view = new HaMolView;
  mol_view->pCanv->resize(size.GetWidth(), size.GetHeight());
  mol_view->host_mol_set = pmset;
  mol_view->InitialTransform();
  mol_view->DefaultRepresentation();

//  mol_view = NULL;
//  if(pmset) mol_view = pmset->GetActiveMolView();

  HeldButton = 0;
#if defined(_MSC_VER)
  DragAcceptFiles(true);
#endif
}
 MolViewWX::~MolViewWX()
 {
	 pmset->canvas_wx  = NULL;
	 pmset->mset_pview = NULL;
 }
void MolViewWX::OnPaint(wxPaintEvent& event)
{
    wxPaintDC dc(this);
    PrepareDC(dc);    // This actually don't do anything for wxWindow
//  dc.SetAxisOrientation(TRUE,FALSE);

    OnDraw(dc);
}


// Define the repainting behaviour
void MolViewWX::OnDraw(wxDC& dc)
{
//	std::cerr << std::endl << "MolVieWX::OnDraw() pt 1 " << mol_view->pCanv->FBuffer << std::endl;

    if(mol_view == NULL || mol_view->pCanv->FBuffer == NULL) 
    {
        PrintLog("MolViewWX::OnDraw() \n");
        PrintLog("Current Molecular View not set \n");
        return;
    }

    mol_view->ReDrawFlag &= ~(RFTransZ|RFPoint);
    
    if( mol_view->ReDrawFlag )
    {       
        if( mol_view->ReDrawFlag & RFReSize )
        {
            mol_view->ReSizeScreen();
        }
        
        if( mol_view->ReDrawFlag & RFColour )
        {
            mol_view->ClearImage();
            mol_view->RefreshColors();
        }
        
        if( !(pmset->HostMolecules.empty()&& pmset->ViewObjects.size() == 0 ) )
        {
            if( mol_view->ReDrawFlag & RFApply ) mol_view->ApplyTransform();
            
            mol_view->DrawFrame();   // Fill the Fbuffer(and axxiliary DBuffer) that determine the image
        }
        else
        {
            mol_view->ClearBuffers();
        }
        mol_view->ReDrawFlag = 0;
    }

    unsigned int x_src = mol_view->pCanv->XRange();
    unsigned int y_src = mol_view->pCanv->YRange();

#if defined(_MSC_VER)
//#if 0    
    int size = sizeof(BITMAPINFOHEADER) + 256*sizeof(RGBQUAD);
    BITMAPINFO* BitInfo = (BITMAPINFO *)_fmalloc( size );

    BitInfo->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    BitInfo->bmiHeader.biCompression = BI_RGB;
    BitInfo->bmiHeader.biXPelsPerMeter = 0;
    BitInfo->bmiHeader.biYPelsPerMeter = 0;
    BitInfo->bmiHeader.biClrImportant = 0;
    BitInfo->bmiHeader.biSizeImage = 0;

    BitInfo->bmiHeader.biBitCount = 32;
    BitInfo->bmiHeader.biPlanes = 1;
 
    BitInfo->bmiHeader.biWidth  = x_src;
    BitInfo->bmiHeader.biHeight = y_src;

    HDC hDC = (HDC)dc.GetHDC();

    HBITMAP PixMap = CreateDIBitmap( hDC, (BITMAPINFOHEADER  *)BitInfo, 
                             CBM_INIT, mol_view->pCanv->FBuffer, BitInfo, DIB_RGB_COLORS);

    HDC memdc = ::CreateCompatibleDC( hDC );

    HGDIOBJ hOldBitmap = ::SelectObject( memdc, PixMap );
    ::BitBlt( hDC, 0, 0, x_src, y_src, memdc, 0, 0, SRCCOPY);
        
    ::SelectObject( memdc, hOldBitmap );
    ::DeleteDC( memdc );
    
    DeleteObject(PixMap);
    delete BitInfo;
#else

    wxImage img(x_src,y_src);
    mol_view->SetWXImage(img);

    wxBitmap btm_img(img);

    wxMemoryDC mem_dc;
    mem_dc.SelectObject(btm_img);

    dc.Blit(0, 0, x_src, y_src, &mem_dc, 0, 0);
    
    mem_dc.SelectObject(wxNullBitmap);

#endif
}

void MolViewWX::OnSize(wxSizeEvent& event)
{
    wxSize size_new = event.GetSize();

    int cx = size_new.GetWidth(); 
    int cy = size_new.GetHeight();

//  if(nType != SIZE_MINIMIZED)
        
    mol_view->pCanv->resize(cx,cy);
    mol_view->ReDrawFlag |= RFReSize;
    mol_view->ClearImage();
    
    pmset->RefreshAllViews();
}

void MolViewWX::OnMouseEvent(wxMouseEvent& event)
{
//  PrintLog( "MolViewWX::OnMouseEvent() \n");
    wxClientDC dc(this);
    PrepareDC(dc);

//  dc.SetPen(*wxBLACK_PEN);

    wxPoint point(event.GetLogicalPosition(dc));
      
    int dx,dy;
    
    if( event.LeftDown() || event.RightDown() || event.MiddleDown())
    {
       InitX = PointX = point.x;
       InitY = PointY = point.y;
       HeldButton = True;
    
       CaptureMouse();     // Capture the mouse until button up.
    }
    else if( event.LeftUp()) 
    {
       if (!HasCapture() )
          return; // If this window (view) didn't capture the mouse,
                   // then the user isn't drawing in this window.
    
       PointX = dx = point.x;
       PointY = dy = point.y;
       if( IsClose(dx,InitX) && IsClose(dy,InitY) )
       {   
          if( event.m_shiftDown || event.m_controlDown )
          {      
             mol_view->PickAtom(True,dx,mol_view->pCanv->YRange()-dy);
          } 
          else 
             mol_view->PickAtom(False,dx,mol_view->pCanv->YRange()-dy);
       }
       ReleaseMouse();   // Release the mouse capture established at
       // the beginning of the mouse drag.      
    }
    else if( event.RightUp() || event.MiddleUp()) 
    {
       if (!HasCapture() )
          return; 
	   ReleaseMouse(); 
	}
    else if(event.Dragging() && 
        (event.m_leftDown || event.m_middleDown || event.m_rightDown || event.m_shiftDown ))
    {
        if( IsClose((int)point.x, InitX) && IsClose((int)point.y, InitY) )
           return;  
        int dx,dy;
        dx = (int)point.x - PointX;
        dy = (int)point.y - PointY;
        MouseMove( event, dx, dy );
        PointX = point.x;
        PointY = point.y;        
        pmset->RefreshAllViews();
    }
}

void MolViewWX::MouseMove( wxMouseEvent& event, int dx, int dy )
{
//	PrintLog(" MolViewWX::MouseMove() dx=%d dy=%d \n",dx,dy); 
    if( mol_view->MouseMode == MMRasMol )
    {   
        if( event.m_shiftDown )
        {   
            if( event.m_leftDown )
            {   
                if( dy ) /* Zoom Vertical */
                {   
                    mol_view->Zoom += (double)dy/mol_view->pCanv->HRange();
                    mol_view->ReDrawFlag |= RFZoom;
//					std::cerr << " MouseMove() Zoom = " << mol_view->Zoom << std::endl;
                }
            } 
            else if( event.m_middleDown | event.m_rightDown )
                if( dx ) /* Z Rotation Horizontal */
                {   
                    mol_view->WrapShiftVal( 2, (double)-dx/mol_view->pCanv->WRange() );
                    mol_view->ReDrawFlag |= RFRotateZ;
                }
        } 
        else if( event.m_controlDown )
        {   
            if( event.m_leftDown )
            {   
                if( dy ) /* Slab Vertical */
                {   
                    mol_view->ClampShiftVal( 7, (double)dy/mol_view->pCanv->YRange() );
                    mol_view->ReDrawFlag |= RFSlab;
                }
            }
            
        } 
        else /* Unmodified! */
            if( event.m_leftDown )
            {   
                if( dx ) /* Rotate Y Horizontal */
                {   
                    mol_view->WrapShiftVal( 1, (double)dx/mol_view->pCanv->WRange() );
                    mol_view->ReDrawFlag |= RFRotateY;
                }
                
                if( dy ) /* Rotate X Vertical */
                {   
                    mol_view->WrapShiftVal( 0, (double)-dy/mol_view->pCanv->HRange() );
                    mol_view->ReDrawFlag |= RFRotateX;
                }
            } 
            else if( event.m_rightDown | event.m_middleDown )
            {   
                if( dx ) /* Translate X Horizontal */
                {   
                    mol_view->ClampShiftVal( 4, (double)dx/mol_view->pCanv->XRange() );
                    mol_view->ReDrawFlag |= RFTransX;
                }
                
                if( dy ) /* Translate Y Vertical */
                {   
                    mol_view->ClampShiftVal( 5, (double)-dy/mol_view->pCanv->YRange() );
                    mol_view->ReDrawFlag |= RFTransY;
                }
            }
    } 
    else if( mol_view->MouseMode == MMQuanta )
    {   
        if( event.m_shiftDown )
        {   
            if( event.m_leftDown )
            {   
                if( dy ) /* Slab Vertical */
                {   
                    mol_view->ClampShiftVal( 7, (double)dy/mol_view->pCanv->YRange() );
                    mol_view->ReDrawFlag |= RFSlab;
                }
            } 
            else if( event.m_middleDown | event.m_rightDown )
            {   
                if( dx ) /* Translate X Horizontal */
                {   
                    mol_view->ClampShiftVal( 4, (double)dx/mol_view->pCanv->XRange() );
                    mol_view->ReDrawFlag |= RFTransX;
                }
                
                if( dy ) /* Translate Y Vertical */
                {   
                    mol_view->ClampShiftVal( 5, (double)-dy/mol_view->pCanv->YRange() );
                    mol_view->ReDrawFlag |= RFTransY;
                }
            } 
            else /* No Mouse Buttons */
                if( dy ) /* Zoom Vertical */
                {   
                    mol_view->ClampShiftVal( 3, (double)dy/mol_view->pCanv->HRange() );
                    mol_view->ReDrawFlag |= RFZoom;
                }
        } 
        else if( event.m_middleDown | event.m_leftDown )
        {   
            if( dx ) /* Rotate Y Horizontal */
            {   
                mol_view->WrapShiftVal( 1, (double)dx/mol_view->pCanv->WRange() );
                mol_view->ReDrawFlag |= RFRotateY;
            }
            
            if( dy ) /* Rotate X Vertical */
            {   
                mol_view->WrapShiftVal( 0, (double)-dy/mol_view->pCanv->HRange() );
                mol_view->ReDrawFlag |= RFRotateX;
            }
        } 
        else if( event.m_rightDown )
        if( dx ) /* Z Rotation Horizontal */
        {   
            mol_view->WrapShiftVal( 2, (double)-dx/mol_view->pCanv->WRange() );
            mol_view->ReDrawFlag |= RFRotateZ;
        }
    
    } 
    else /* MMInsight */
    {   
        if( event.m_leftDown && !(event.m_middleDown | event.m_rightDown) )
        {
            if( dx ) /* Rotate Y Horizontal */
            {   
                mol_view->WrapShiftVal( 1, (double)dx/mol_view->pCanv->WRange() );
                mol_view->ReDrawFlag |= RFRotateY;
            }
            
            if( dy ) /* Rotate X Vertical */
            {   
                mol_view->WrapShiftVal( 0, (double)dy/mol_view->pCanv->HRange() );
                mol_view->ReDrawFlag |= RFRotateX;
            }
        }
        else if( event.m_middleDown && !(event.m_leftDown | event.m_rightDown))
        {
            if( dx ) /* Translate X Horizontal */
            {   
                mol_view->ClampShiftVal( 4, (double)dx/mol_view->pCanv->XRange() );
                mol_view->ReDrawFlag |= RFTransX;
            }
            
            if( dy ) /* Translate Y Vertical */
            {   
                mol_view->ClampShiftVal( 5, (double)dy/mol_view->pCanv->YRange() );
                mol_view->ReDrawFlag |= RFTransY;
            }
        }   
        else if( event.m_leftDown & event.m_middleDown && !(event.m_rightDown))
        {
            mol_view->ClampShiftVal( 3, (double)dx/mol_view->pCanv->WRange() - (double)dy/mol_view->pCanv->HRange() );
            mol_view->ReDrawFlag |= RFZoom;
        }   
        else if( (event.m_leftDown & event.m_rightDown) && !(event.m_middleDown) )
        {
            mol_view->WrapShiftVal( 2, (double)dx/mol_view->pCanv->WRange() - (double)dy/mol_view->pCanv->HRange() );
            mol_view->ReDrawFlag |= RFRotateZ;
        }   
        else if( event.m_leftDown & event.m_middleDown & event.m_rightDown )
        {
            mol_view->ClampShiftVal( 7, (double)dx/mol_view->pCanv->XRange() - (double)dy/mol_view->pCanv->YRange() );
            mol_view->ReDrawFlag |= RFSlab;
        }
    }   
}

void MolViewWX::OnEraseBackground(wxEraseEvent& event)
{
//  PrintLog("MolViewWX::OnEraseBackground() \n");
//    RECT rect;
//    ::GetClientRect(GetHwnd(), &rect);

//    COLORREF ref = PALETTERGB(m_backgroundColour.Red(),
//                              m_backgroundColour.Green(),
//                              m_backgroundColour.Blue());
//    HBRUSH hBrush = ::CreateSolidBrush(ref);
//    if ( !hBrush )
//        wxLogLastError(wxT("CreateSolidBrush"));

//    HDC hdc = (HDC)event.GetDC()->GetHDC();

//    int mode = ::SetMapMode(hdc, MM_TEXT);

//    ::FillRect(hdc, &rect, hBrush);
//    ::DeleteObject(hBrush);
//    ::SetMapMode(hdc, mode);
}


void MolViewWX::OnClose(wxCloseEvent& event)
{
    MolSet* pmset = this->mol_view->GetMolSet();
    wxString name = pmset->GetName();
    wxString msg = "Delete Molecular Set ";
    msg += pmset->GetName();
    msg += " associated with a view ?";

    int answer = ::wxMessageBox(msg,"Confirm", wxYES_NO | wxCANCEL );

	if( answer == wxCANCEL && event.CanVeto()) 
	{
		PrintLog(" answer = wxCANCEL \n");
		event.Veto();
		return;
	}

	delete this->mol_view;

    if( answer == wxYES)
    {
        delete pmset;
    }
	else
	{
		pmset->canvas_wx  = NULL;
		pmset->mset_pview = NULL;
	}

	SetCurMolSet(NULL);

//    PrintLog(" MolViewWX::OnClose() \n");
}

int HaMolView::SetWXImage(wxImage& wx_image) //!< Set wxImage with the current image data
{
    unsigned int x_src = pCanv->XRange();
    unsigned int y_src = pCanv->YRange();

    wx_image.Rescale(x_src,y_src);
    if( pCanv->FBuffer == NULL || x_src == 0 || y_src == 0) 
    {
        PrintLog("HaMolView::SetWXImage() \n");
        PrintLog(" Image data not set");
        return FALSE;
    }

    unsigned char* img_data = (unsigned char*) malloc( x_src*y_src*3);
    
    int n= x_src*y_src;
    unsigned int i;
    int j = 0;

    uint_4* ptr_src = pCanv->FBuffer + x_src*y_src;
    unsigned char* ptr_dst = img_data;
    unsigned char* ptr_cur;

    ptr_src = ptr_src + x_src;
    for(i = 0; i < n; i++)
    {
        if( (i % x_src) == 0) ptr_src = ptr_src - (2*x_src);
        ptr_cur = (unsigned char*)ptr_src;
        *ptr_dst = *(ptr_cur+2);  // setting red
        ptr_dst++; 
        *ptr_dst = *(ptr_cur+1);  // setting green
        ptr_dst++; 
        *ptr_dst = *ptr_cur;  // setting blue
        ptr_dst++;
        ptr_src++; 
    }
    wx_image.SetData(img_data); 
    return TRUE;
}

MolViewFrame::MolViewFrame(wxMDIParentFrame* parent,wxString& mset_name):
wxMDIChildFrame(parent, -1, mset_name, wxPoint(10, 10),
                wxDefaultSize, wxDEFAULT_FRAME_STYLE)
{

}

BEGIN_EVENT_TABLE(MolViewFrame, wxMDIChildFrame)
    EVT_CLOSE(MolViewFrame::OnClose)
	EVT_ACTIVATE(MolViewFrame::OnActivate)
END_EVENT_TABLE()

void MolViewFrame::OnClose(wxCloseEvent& event )
{
//  PrintLog(" MolViewFrame::OnClose() \n");
    bool bres = mol_view_wx->Close();
    if( bres)
	{
		DestroyChildren();
		Destroy();
	}
}


void MolViewFrame::OnActivate(wxActivateEvent& event)
{
	CurMolView = mol_view_wx->mol_view;

	MolSet* pmset = mol_view_wx->mol_view->GetMolSet();
//	PrintLog(" MolViewFrame::OnActivate()  set CurMolSet to: %s \n", pmset->GetName() ); 
	SetCurMolSet(pmset);
	
	event.Skip();
}

MolViewWX* HaMainFrameWX::CreateMolView(MolSet* pmset)
{
    MolViewFrame* subframe;
    wxString mset_name= "EMPTY MOLSET";
    if(pmset) mset_name = pmset->GetName();
    subframe = new MolViewFrame(this, mset_name);

#ifdef __WXMSW__
  subframe->SetIcon(wxString("notepad"));
#endif
#ifdef __X__
  subframe->SetIcon(wxIcon("doc.xbm"));
#endif

  int width, height;
  subframe->GetClientSize(&width, &height);

  MolViewWX* pview_wx = new MolViewWX(pmset, subframe, wxPoint(0,0), wxSize(width, height),0); 

  subframe->mol_view_wx = pview_wx;
  wxCursor curs(wxCURSOR_ARROW);  
  pview_wx->SetCursor(curs);
//  pview_wx->SetCursor(wxCursor(wxCURSOR_PENCIL));
  
  pmset->mset_pview = pview_wx->mol_view;
  pmset->canvas_wx = pview_wx;
  
#ifdef __X__
    // X seems to require a forced resize
    int x, y;
    GetSize(&x, &y);
    SetSize(-1, -1, x, y);
#endif

  subframe->Show(TRUE);
  pview_wx->mol_view->InitialTransform();
  pview_wx->mol_view->DefaultRepresentation();
  pmset->RefreshAllViews(RFRefresh | RFColour | RFApply);
  return subframe->mol_view_wx;
}

//IMPLEMENT_CLASS(HaMainFrameWX, wxMDIParentFrame)

// WDR: event table for HaMainFrameWX

BEGIN_EVENT_TABLE(HaMainFrameWX, wxMDIParentFrame)
    EVT_MOVE(HaMainFrameWX::OnMove)
	EVT_SIZE(HaMainFrameWX::OnSize)
	EVT_CLOSE(HaMainFrameWX::OnClose)
	EVT_BUTTON( IDC_CMD_EXEC, HaMainFrameWX::OnExecuteCommand )
    EVT_TEXT_ENTER( IDC_CMD_TXT, HaMainFrameWX::OnExecuteCommand )
//File Menu
    EVT_MENU( ID_FILE_NEW_WX,         HaMainFrameWX::OnFileNew )
    EVT_MENU( ID_FILE_OPEN_WX,        HaMainFrameWX::OnFileOpen )
    EVT_MENU( IDM_MOL_SAVE_WX,        HaMainFrameWX::OnMolSave )
    EVT_MENU( IDM_MOL_SAVE_AS_WX,     HaMainFrameWX::OnMolSaveAs )
    EVT_MENU( IDM_LOAD_SCRIPT_WX,     HaMainFrameWX::DoLoadScriptDialog )
    EVT_MENU( IDM_REDIRECT_IO_WX,     HaMainFrameWX::DoRedirectIODialog )
    EVT_MENU( IDM_SAVE_IMAGE_WX,      HaMainFrameWX::SaveOutputFile )
    EVT_MENU( IDM_SAVE_IMAGE_CLIPBOARD_WX, HaMainFrameWX::OnSaveImageClipboard )
    EVT_MENU( IDM_COMP_ACCOUNTS_WX,   HaMainFrameWX::DoCompAccountsDialog )
    EVT_MENU( IDM_INFO_WX,            HaMainFrameWX::OnInfo )
    EVT_MENU( ID_FILE_CLOSE_WX,       HaMainFrameWX::OnFileClose )
    EVT_MENU( IDM_PRINT_WX,           HaMainFrameWX::OnPrint )
    EVT_MENU( IDM_SETUP_WX,           HaMainFrameWX::OnSetup )
    EVT_MENU( IDM_EXIT_WX,            HaMainFrameWX::OnExit )
    EVT_MENU( IDM_PYMOD_WX,            HaMainFrameWX::OnPyMod )//<mikola Jul 19, 2006
// Edit Menu
    EVT_MENU( IDM_SELECT_WX,              HaMainFrameWX::OnSelectAll )
    EVT_MENU( IDM_EDIT_ATOM_PARAM_WX,     HaMainFrameWX::DoAtomParamsDialog )
    EVT_MENU( IDM_EDIT_RES_PARAM_WX,      HaMainFrameWX::DoResidueParamsDialog )
    EVT_MENU( IDM_EDIT_MOLSETS,           HaMainFrameWX::DoMolSetParamDialog )
	EVT_MENU( IDM_ATOM_PROP_COLORS,       HaMainFrameWX::DoAtomPropColorDialog )
    EVT_MENU( IDM_EDIT_GROUPS_WX,         HaMainFrameWX::DoEditGroupsDialog )
	EVT_MENU( IDM_CRD_SNAPSHOT,           HaMainFrameWX::DoCrdSnapshotDialog )
    EVT_MENU( IDM_EDIT_GEOM_WX,           HaMainFrameWX::DoEditGeomDialog )
    EVT_MENU( IDM_FIND_HBONDS_WX,         HaMainFrameWX::OnFindHbonds )
    EVT_MENU( IDM_SOLVATE_WX,             HaMainFrameWX::DoSolvateDialog )
	EVT_MENU( IDM_WRAP_UNIT_CELL,         HaMainFrameWX::OnWrapIntoUnitCell )
	EVT_MENU( IDM_CENTER_SOLUTE,          HaMainFrameWX::OnCenterSolute )
	EVT_MENU( IDM_CENTER_MOL_PBOX,        HaMainFrameWX::OnCenterMolInPBox )
    EVT_MENU( IDM_ADD_MISSING_ATOMS_WX,   HaMainFrameWX::OnAddMissingAtoms )
    EVT_MENU( IDM_ADD_POLAR_HYDROGENS_WX, HaMainFrameWX::OnAddPolarHydrogens )
    EVT_MENU( IDM_ADD_H_HYBRID_WX,        HaMainFrameWX::OnAddHHybrid )
    EVT_MENU( IDM_ADD_HYDROGENS_WX,       HaMainFrameWX::OnAddHydrogens )
    EVT_MENU( IDM_DEL_SEL_ATOMS_WX,       HaMainFrameWX::OnDelSelAtoms )
	EVT_MENU( IDM_DEL_OVLP_MOLS,          HaMainFrameWX::OnDelOvlpMols )
    EVT_MENU( IDM_EDIT_FRAGM_WX,          HaMainFrameWX::DoEditFragmDialog )
    EVT_MENU( IDM_BUILD_FILM_WX,          HaMainFrameWX::DoBuildFilmDialog )
    EVT_MENU( IDM_CLEAR_PICKED_WX,        HaMainFrameWX::OnClearPicked )
	EVT_MENU( IDM_SEL_ATOMS_IN_BOUND_BOX, HaMainFrameWX::OnSelAtomsInBoundBox )
	EVT_MENU( IDM_REVERT_SELECTION,       HaMainFrameWX::OnRevertAtomSelection )
    EVT_MENU( IDM_NUCL_ACID_WX,           HaMainFrameWX::DoNuclAcidDialog )
	EVT_MENU( IDM_DESCRIBE_SEC_STRUCT,    HaMainFrameWX::OnDescribeSecStruct )
	EVT_MENU( IDM_PRINT_HBONDS,           HaMainFrameWX::OnPrintHBonds )
	EVT_MENU( IDM_SET_ALPHA_HELIX,        HaMainFrameWX::OnSetAlphaHelix )
    EVT_MENU( IDM_SHOW_RES_DB_WX,         HaMainFrameWX::OnShowResDb )
// Display Menu
//   Display->Display Mode menu
    EVT_MENU( IDM_WIREFRAME_WX,    HaMainFrameWX::OnWireFrame )
    EVT_MENU( IDM_BACKBONE_WX,     HaMainFrameWX::OnBackBone )
    EVT_MENU( IDM_STICKS_WX,       HaMainFrameWX::OnSticks )
    EVT_MENU( IDM_SPHERES_WX,      HaMainFrameWX::OnSpheres )
    EVT_MENU( IDM_BALLSTICK_WX,    HaMainFrameWX::OnBallStick )
    EVT_MENU( IDM_RIBBONS_WX,      HaMainFrameWX::OnRibbons )
    EVT_MENU( IDM_STRANDS_WX,      HaMainFrameWX::OnStrands )
    EVT_MENU( IDM_CARTOONS_WX,     HaMainFrameWX::OnCartoons )
//   Display->Colours menu
    EVT_MENU( IDM_MONO_WX,         HaMainFrameWX::OnMono )
    EVT_MENU( IDM_CPK_WX,          HaMainFrameWX::OnCPK )
    EVT_MENU( IDM_SHAPELY_WX,      HaMainFrameWX::OnShapely )
    EVT_MENU( IDM_COL_GROUPS_WX,   HaMainFrameWX::OnColGroups )
    EVT_MENU( IDM_COL_RESIDUES_WX, HaMainFrameWX::OnColResidues )
    EVT_MENU( IDM_CHAIN_WX,        HaMainFrameWX::OnColChain )
    EVT_MENU( IDM_TEMPER_WX,       HaMainFrameWX::OnColTemper )
    EVT_MENU( IDM_STRUCT_WX,       HaMainFrameWX::OnStruct )
    EVT_MENU( IDM_COL_DON_COUPL_WX,HaMainFrameWX::OnColTemper )
//   Display->Options menu
    EVT_MENU( IDM_SLAB_WX,         HaMainFrameWX::OnSlab )
    EVT_MENU( IDM_HYDROGEN_WX,     HaMainFrameWX::OnHydrogen )
    EVT_MENU( IDM_HETERO_WX,       HaMainFrameWX::OnHetero )
    EVT_MENU( IDM_SPECULAR_WX,     HaMainFrameWX::OnSpecular )
    EVT_MENU( IDM_STEREO_WX,       HaMainFrameWX::OnStereo )
    EVT_MENU( IDM_LABELS_WX,       HaMainFrameWX::OnLabels )
	EVT_MENU( IDM_DISPLAY_PBOX,    HaMainFrameWX::OnDisplayPBox )
//
    EVT_MENU( IDM_VIEW_PARAM_WX,   HaMainFrameWX::DoViewParamDlg )
    EVT_MENU( IDM_OBJECT3D_WX,     HaMainFrameWX::DoObject3DDialog )
//   Display->Labels menu
    EVT_MENU( IDM_LABELS_ATOM_ID_WX,     HaMainFrameWX::OnLabelsAtomId )
    EVT_MENU( IDM_LABELS_ATOM_NAMES_WX,  HaMainFrameWX::OnLabelsAtomNames )
    EVT_MENU( IDM_LABELS_ATOM_SEQNUM_WX, HaMainFrameWX::OnLabelsAtomSeqNum )
	EVT_MENU( IDM_LABELS_ATOM_FF_SYMB,   HaMainFrameWX::OnLabelsAtomFFSymb )
    EVT_MENU( IDM_LABELS_OFF_WX ,        HaMainFrameWX::OnLabelsOff )
// Display->Windows menu
    EVT_MENU( ID_WINDOW_CASCADE_WX,      HaMainFrameWX::OnWindowsCascade )
    EVT_MENU( ID_WINDOW_TILE_VERT_WX,    HaMainFrameWX::OnWindowsTileVert )
    EVT_MENU( IDM_CENTER_VIEW_SEL_WX,    HaMainFrameWX::OnCenterViewSel )   
// Measure Menu
    EVT_MENU( IDM_SHOW_ATOMID_WX,      HaMainFrameWX::OnShowAtomID )
    EVT_MENU( IDM_MEASURE_DISTANCE_WX, HaMainFrameWX::OnMeasureDist )
    EVT_MENU( IDM_MEASURE_ANGLE_WX,    HaMainFrameWX::OnMeasureAngle )
    EVT_MENU( IDM_MEASURE_DIHEDRAL_WX, HaMainFrameWX::OnMeasureDihed )
// ET Menu
    EVT_MENU( IDM_EDIT_REDOX_WX,       HaMainFrameWX::OnEditRedox )
    EVT_MENU( IDM_RUN_PATHWAYS_WX,     HaMainFrameWX::DoPathwaysDialog )
    EVT_MENU( IDM_ET_EFF_HAM_WX,       HaMainFrameWX::DoETEffHamDialog )
    EVT_MENU( IDM_RESET_ET_WX,         HaMainFrameWX::OnResetETModule )
    EVT_MENU( IDM_EL_SCATTER_WX,       HaMainFrameWX::DoScatterDialog )
// QChem Menu
    EVT_MENU( IDM_QCHEM_PARAM_WX,       HaMainFrameWX::DoQChemParamDialog )
    EVT_MENU( IDM_QCHEM_LOAD_DATA_WX,   HaMainFrameWX::DoLoadQCDatDialog )
    EVT_MENU( IDM_WAVEFUN_ANAL_WX,      HaMainFrameWX::DoWaveFunAnalDialog )
    EVT_MENU( IDM_CALC_POLAR_GCONTR_WX, HaMainFrameWX::OnCalcPolarGcontr )
    EVT_MENU( IDM_SAVE_GROUP_OPER_WX,   HaMainFrameWX::OnSaveGrpOperMat )
    EVT_MENU( IDM_POL_GRP_CONTR_F_WX,   HaMainFrameWX::OnCalcPolarContrF )
    EVT_MENU( IDM_READ_POLAR_CONTR_WX,  HaMainFrameWX::OnReadPolarContr )
    EVT_MENU( IDM_CALC_ROTANG_GRP_CONTR2_WX, HaMainFrameWX::OnCalcBetaContr2idx )
    EVT_MENU( IDM_READ_ROTANG_GRP_CONTR2_WX, HaMainFrameWX::OnReadBetaContr2idx )
// MMech Menu
    EVT_MENU( IDM_MM_PARAM_WX,         HaMainFrameWX::DoMolMechDialog )
    EVT_MENU( IDM_INTER_MOL_WX,        HaMainFrameWX::DoInterMolDialog )
    EVT_MENU( IDM_CONT_ELECTR_WX,      HaMainFrameWX::DoContElectrDialog )
    EVT_MENU( IDM_PNP,                 HaMainFrameWX::DoPNPDialog )
    //EVT_MENU( IDM_CONT_ELECTR_WX_PNP,  HaMainFrameWX::DoContElectrPNPDialog )
    EVT_MENU( IDM_PNP_SMLPNP,          HaMainFrameWX::DoSmlPNPDialog )
    EVT_MENU( IDM_PNP_RF,          HaMainFrameWX::DoRFDialog )
	EVT_MENU( ID_FLEX_MOD,          HaMainFrameWX::DoFlexModDialog )
	EVT_MENU( ID_MENU_ED,              HaMainFrameWX::DoEDDialog )
	EVT_MENU( IDM_PROT_REDOX,          HaMainFrameWX::DoProtRedoxDialog )
//Test Menu
		EVT_MENU( IDM_OPENVECFIELD_HAVIEW,HaMainFrameWX::OnOpenVectorFieldInHaMolView )
		EVT_MENU( IDM_OPEN_EPS_NIND,HaMainFrameWX::OnOpenDielConstFromNindexInHaMolView )
		EVT_MENU( IDM_OPEN_DIFF_NIND,HaMainFrameWX::OnOpenDiffConstFromNindexInHaMolView )
		EVT_MENU( IDM_OPEN_QST_NIND,HaMainFrameWX::OnOpenQstFromNindexInHaMolView )
// Main Toolbar
    EVT_MENU( IDM_MANUAL_WX,           HaMainFrameWX::OnLoadManual )
    EVT_MENU( IDM_MOL_CONNECT_WX,      HaMainFrameWX::OnMolConnect )
    EVT_MENU( IDM_WORLD_CONNECT_WX,    HaMainFrameWX::OnWorldConnect )
    EVT_MENU( IDM_MOL_INFO_WX, HaMainFrameWX::OnDumpMolInfo )
    EVT_MENU( IDM_GAUSS_BCOMMON_WX, HaMainFrameWX::OnDumpGaussBCommon )
    EVT_MENU( IDM_DUMP_OVERLAP, HaMainFrameWX::OnDumpOverlap )
    EVT_MENU( IDM_DUMP_OVERLAP2_WX, HaMainFrameWX::OnDumpOverlap2 )
    EVT_MENU( IDM_TEST_OPER_1_WX, HaMainFrameWX::OnTestOper1 )
    EVT_MENU( IDM_TEST_OPER_2_WX, HaMainFrameWX::OnTestOper2 )
    EVT_MENU( IDM_TEST_QCMOD_1_WX, HaMainFrameWX::OnTestQCMod1 )
    EVT_MENU( IDM_CALC_MATR_1_WX, HaMainFrameWX::OnCalcTMatr1 )
    EVT_MENU( IDM_TEST_MAP_WX, HaMainFrameWX::OnTestMap )
    EVT_MENU( IDM_TEST_MIN_1_WX, HaMainFrameWX::OnTestMin1 )
    EVT_MENU( IDM_TEST_GRAPH_1_WX, HaMainFrameWX::OnTestGraph1 )
	EVT_MENU( IDM_TEST_PYTHON_1, HaMainFrameWX::OnTestPython1 )
	EVT_MENU( IDM_TEST_AVG_POP_FUNC, HaMainFrameWX::OnTestAvgPopFunc )
	EVT_MENU( IDM_SAVE_AMOEBA_TOP_FILE_1, HaMainFrameWX::OnSaveAmoebaTopFile1 )
	EVT_MENU( IDM_SAVE_AMOEBA_TOP_FILE_2, HaMainFrameWX::OnSaveAmoebaTopFile2 )
	EVT_MENU( IDM_TEST_FFT_1, HaMainFrameWX::OnTestFFT1 )
    EVT_MENU( IDM_ABOUT_WX, HaMainFrameWX::OnAbout )
    EVT_MENU( IDM_HELP_WX, HaMainFrameWX::OnHelp )
//
    
END_EVENT_TABLE()

HaMainFrameWX::HaMainFrameWX() :
    wxMDIParentFrame(NULL, -1, "HARLEM", wxPoint(0, 0), wxDefaultSize, wxDEFAULT_FRAME_STYLE | wxHSCROLL | wxVSCROLL , "HaMainFrameWX")
{
		m_HaMainFrameWX = this;
	int xs = wxSystemSettings::GetMetric(wxSYS_SCREEN_X);
    int ys = wxSystemSettings::GetMetric(wxSYS_SCREEN_Y);

    int xp = MinFun(xs*3/4,1200);
    int yp = MinFun(ys*3/4,900);

	this->SetInitialSize( wxSize(xp, yp) );

	 wxMenuBar* main_menu_bar = MainMenu();
	 SetMenuBar(main_menu_bar);    

//	 CreateToolBar();
//   MainToolBarFunc(this->GetToolBar());
    
	 wxSize size_w(1000,50);
	
	 sash_win = new wxSashLayoutWindow( this, -1, wxDefaultPosition, size_w, wxNO_BORDER | wxSW_3D );

	 sash_win->SetDefaultSize( wxSize(1000,100) );
	 sash_win->SetAlignment(wxLAYOUT_BOTTOM);
	 sash_win->SetOrientation( wxLAYOUT_HORIZONTAL );
	 sash_win->SetSashVisible( wxSASH_BOTTOM, true );
	 
	 bot_bar_dlg( sash_win, true, true );
	 sash_win->SetAutoLayout( true );
	 sash_win->Layout();

	 sash_win->Show();
	 sash_win->Refresh();

//     wxGraphs::Parent=this;//<mikola, July 20, 2006

}

void HaMainFrameWX::OnExecuteCommand( wxCommandEvent &event )
{
    wxTextCtrl* cmd_txt_edt = (wxTextCtrl*) FindWindow(IDC_CMD_TXT);
    wxString cmd = cmd_txt_edt->GetValue();
	pApp->cmd_pr.SetCmdLine( cmd.ToStdString() );
    pApp->ExecuteCommand();
}

void HaMainFrameWX::OnClose(wxCloseEvent& event)
{
	PrintLog("HaMainFrameWX::OnClose() \n");
	if( pApp->mpi_driver->nprocs > 1 )
	{
		std::string msg = HaMPI::BuildXMLwxCmdEventBasic(wxEVT_HARLEM_APP, HARLEM_APP_EXIT_ID);
	    pApp->mpi_driver->SendXmlMsgAllProc(msg.c_str());
	}

	this->DestroyChildren();
	this->Destroy();
	

//	if(pApp->python_thread != NULL)
//	{
//		((wxThread*)pApp->python_thread)->Kill();
//		pApp->python_thread = NULL;
//	}
//	PrintLog("HaMainFrameWX::OnClose() pt 2 \n");
}

// WDR: handler implementations for HaMainFrameWX

void HaMainFrameWX::OnHelp( wxCommandEvent &event )
{
#if defined(_MSC_VER)
	std::string help_file;
    help_file = pApp->harlem_home_dir;
    int len = help_file.length();
    if( (char)help_file[len-1] != '\\') help_file += '\\';
    help_file += "harlemwin.hlp";
    ::WinHelpA(NULL,help_file.c_str(),HELP_INDEX,0L);   
#endif
}

void HaMainFrameWX::OnAbout(wxCommandEvent& WXUNUSED(event))
{
	wxString Mes="HARLEM: HAmiltonians to Research LargE Molecules \n\n";      
	Mes+="Igor Kurnikov, Nikolay Simakov, Kirill Speransky, \n"; 
	Mes+="Maria Kurnikova 1997 - 2018 \n\n"; 
    Mes+="Graphical Interface Based on RASMOL 2.6 of Roger Sayle \n";
	Mes+="and wxWidgets library \n";
	Mes+="Quantum Chemical functionality based on IPACK library of Wolfgang Wenzel \n";
	Mes+="Command line/scripting interface is based on PYTHON 3.7 \n";
	Mes+="Also linked to VFLIB, TINYXML, LAPACK, BLAS and other libraries \n\n";
	Mes+="Build date -  ";Mes+=__DATE__;Mes+="\n";        
//	Mes+=HaSVNRevision();Mes+="\n";             
//	Mes+=HaSVNDate();              
	wxMessageBox(Mes, wxString("About HARLEM"));               
}    
    
void HaMainFrameWX::OnMove(wxMoveEvent& event)
{
	wxLayoutAlgorithm().LayoutMDIFrame(this);
	event.Skip();
}

void HaMainFrameWX::OnSize(wxSizeEvent& event)
{
	wxLayoutAlgorithm().LayoutMDIFrame(this);
}

void MolSet::SetName(const char* new_name)
{
    name_mset = new_name;
    if(canvas_wx != NULL) canvas_wx->mol_frame->SetTitle(new_name);
}

void HaMolView::UpdateThisView( int lHint)
{
    lHint |= RFRefresh;
    ReDrawFlag |= lHint;
    MolSet* pmset = GetMolSet();
    if(pmset->canvas_wx) pmset->canvas_wx->Refresh();
}

// File Menu
void HaMainFrameWX::OnFileNew(wxCommandEvent &event)
{
	MolSet* pmset = new MolSet();
	if(pmset)
	{
		pmset->canvas_wx = CreateMolView(pmset);
	}
}

void HaMainFrameWX::OnFileOpen(wxCommandEvent &event)
{
	ChooseMolFileDlg load_dlg(NULL, -1, wxString("Open Molecular File"));
	
	wxCheckBox* new_view_chk = new wxCheckBox( &load_dlg, -1, wxT("Load into a new View"), wxDefaultPosition, wxDefaultSize, 0 );
    load_dlg.sizer_main_v->Add( new_view_chk, 0, wxALL, 5 );

    wxCheckBox* center_chk = new wxCheckBox( &load_dlg, -1, wxT("Translate to the center of the screen"), wxDefaultPosition, wxDefaultSize, 0 );
    load_dlg.sizer_main_v->Add( center_chk, 0, wxALL, 5 );

    wxCheckBox* calc_bonds_chk = new wxCheckBox( &load_dlg, -1, wxT("Build bonds from atom-atom distance criteria"), wxDefaultPosition, wxDefaultSize, 0 );
//	calc_bonds_chk->SetValidator(wxGenericValidator(&(MolSet::p_load_opt_default->calc_bonds)));
    load_dlg.sizer_main_v->Add( calc_bonds_chk, 0, wxALL, 5 );
	load_dlg.sizer_main_v->SetSizeHints( &load_dlg );
	load_dlg.ShowModal();

	::wxSetWorkingDirectory(load_dlg.dir_name);
     
    MolSet* pmset = GetCurMolSet();
    if( pmset == NULL || pmset->canvas_wx == NULL)
    {
        pmset = new MolSet();
        pmset->canvas_wx = CreateMolView(pmset);        
    }

	pmset->FetchFile(load_dlg.file_format,load_dlg.file_name_full.c_str());
	pmset->canvas_wx->mol_view->InitialTransform();
    pmset->canvas_wx->mol_view->DefaultRepresentation();
    pmset->RefreshAllViews(RFRefresh | RFColour | RFApply);
}

void HaMainFrameWX::OnMolSave(wxCommandEvent &event)
{
	OnMolSaveAs(event);
}

void HaMainFrameWX::OnMolSaveAs(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet();
	if(pmset == NULL) return;
	SaveMolFileDlg* save_dlg = new SaveMolFileDlg(pmset, NULL, -1, wxString("Save Molecular File"));
    save_dlg->Show(TRUE);	
}

void HaMainFrameWX::DoLoadScriptDialog(wxCommandEvent &event)
{
	if(LoadScriptDlgWX::dlg_open) return;
	LoadScriptDlgWX* load_dlg = new LoadScriptDlgWX(NULL, -1, wxString("Load Script"));
    load_dlg->Show(TRUE);
}


void HaMainFrameWX::DoRedirectIODialog(wxCommandEvent &event)
{
	if(RedirectIODlgWX::dlg_open) return; 
	RedirectIODlgWX* ptr_redirect_io_dlg = new RedirectIODlgWX(this);
	ptr_redirect_io_dlg->Show(TRUE);
}

void HaMainFrameWX::SaveOutputFile(wxCommandEvent &event)
{
	HaMolView* pview = CurMolView;
	if( pview != NULL )
	{
	   SaveImageFileDlg* save_dlg = new SaveImageFileDlg(pview,NULL, -1, wxString("Save Image File"));
       save_dlg->Show(TRUE);
	}
}

void HaMainFrameWX::OnSaveImageClipboard(wxCommandEvent &event)
{
	HaMolView* pview = CurMolView;
	
	if(pview != NULL)
	{
		unsigned int x_src = pview->pCanv->XRange();
		unsigned int y_src = pview->pCanv->YRange();
		wxImage img(x_src,y_src);
		int ires = pview->SetWXImage(img);
		if(!ires) return;
		
		wxBitmap btm_img(img);

#if defined(_MSC_VER)
		if (wxTheClipboard->Open())
		{
			wxTheClipboard->SetData(new wxBitmapDataObject(btm_img));
			wxTheClipboard->Close();
		}
#endif
	}
}

void
HaMainFrameWX::DoCompAccountsDialog(wxCommandEvent &event)
{
	if(CompAccountsDlg::dlg_open) return; 
	CompAccountsDlg* ptr_comp_accounts_dlg = new CompAccountsDlg(this);
	ptr_comp_accounts_dlg->Show(TRUE);
}

void
HaMainFrameWX::OnInfo(wxCommandEvent &event)
{

}

void
HaMainFrameWX::OnFileClose(wxCommandEvent &event)
{
	
}

void
HaMainFrameWX::OnPrint(wxCommandEvent &event)
{
	HaMolView* pview = CurMolView;
	if(pview == NULL) return;

#if defined(_MSC_VER)
   
    char *device, *driver, *output;
    int xsize, xres, yres;
    int dx, dy, caps;
    char printer[80];

    DOCINFOA info;
    RECT rect;
    HDC hDC;

    GetProfileStringA("windows","device", "", printer, 80 );
    if( !(device = strtok(printer,",")) ) return;
    if( !(driver = strtok((char*)NULL,", ")) ) return;
    if( !(output = strtok((char*)NULL,", ")) ) return;

    hDC = CreateDCA(driver,device,output,NULL);
    if( !hDC ) return;

    caps = GetDeviceCaps( hDC, RASTERCAPS );
    if( !(caps & RC_STRETCHDIB) ) return;
    
    xres = GetDeviceCaps( hDC, LOGPIXELSX );
    yres = GetDeviceCaps( hDC, LOGPIXELSY );
    xsize = GetDeviceCaps( hDC, HORZRES );

    dx = xsize - xres;
    dy = (int)(((long)dx*pview->pCanv->YRange())/pview->pCanv->XRange());

    /* Should set printer abort procedure */
    /* Position Image on Printed Page */
    rect.top = yres;        rect.bottom = rect.top + dy;
    rect.left = xres>>1;    rect.right = rect.left + dx;
    Escape( hDC, SET_BOUNDS, sizeof(RECT), (char *)&rect, NULL );

    /* Start HARLEM Document */
    info.cbSize = sizeof(DOCINFO);
    info.lpszDocName = "HARLEM";
    info.lpszOutput = NULL;
    StartDocA( hDC, &info );
    StartPage( hDC );

    int size = sizeof(BITMAPINFOHEADER) + 256*sizeof(RGBQUAD);
    BITMAPINFO* BitInfo = (BITMAPINFO *)_fmalloc( size );

    BitInfo->bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    BitInfo->bmiHeader.biCompression = BI_RGB;
    BitInfo->bmiHeader.biXPelsPerMeter = 0;
    BitInfo->bmiHeader.biYPelsPerMeter = 0;
    BitInfo->bmiHeader.biClrImportant = 0;
    BitInfo->bmiHeader.biSizeImage = 0;

    BitInfo->bmiHeader.biBitCount = 32;
    BitInfo->bmiHeader.biPlanes = 1;

    BitInfo->bmiHeader.biWidth = pview->pCanv->XRange();
    BitInfo->bmiHeader.biHeight = pview->pCanv->YRange();

    StretchDIBits( hDC, xres>>1, yres, dx, dy, 
                        0, 0, pview->pCanv->XRange(), pview->pCanv->YRange(), 
                        pview->pCanv->FBuffer, BitInfo, DIB_RGB_COLORS, SRCCOPY );

    EndPage( hDC );
    EndDoc( hDC );

    DeleteDC( hDC );
	delete BitInfo;
#endif
}

void HaMainFrameWX::OnSetup(wxCommandEvent &event)
{
	wxPrintDialog print_dlg(this);
    print_dlg.ShowModal();
}

void HaMainFrameWX::OnExit(wxCommandEvent &event)
{
	this->Close();
}
//>mikola Jul 19, 2006
void HaMainFrameWX::OnPyMod(wxCommandEvent &event)
{
//   fprintf(stdout,"HaMainFrameWX::OnPyMod %d\n",(int)wxPyMod::dlg_open);
//   if(wxPyMod::dlg_open) return;
// 
//   wxPyMod* py_mod_dlg = new wxPyMod( this );
//   py_mod_dlg->Show(TRUE);
	 pApp->CreateCommandWindow();
}
//<mikola Jul 19, 2006
// Edit Menu
void HaMainFrameWX::OnSelectAll(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	if(pmset == NULL) return;
	HaMolView* pview = pmset->GetActiveMolView();
	if( pview == NULL) return;
	int mask = NormAtomFlag;
	if( pview->HetaGroups ) mask |= HeteroFlag;
	if( pview->Hydrogens )  mask |= HydrogenFlag;
	pmset->SelectAtomsMask(mask);	
}


void
HaMainFrameWX::DoAtomParamsDialog(wxCommandEvent &event)
{
	if(AtomParamsDlgWX::dlg_open) return;
	MolSet* pmset = GetCurMolSet(); 
	if( pmset == NULL) return;
	AtomParamsDlgWX* ptr_atom_params_dlg = new AtomParamsDlgWX( pmset, this );
	ptr_atom_params_dlg->Show(TRUE);	
}

void HaMainFrameWX::DoResidueParamsDialog(wxCommandEvent &event)
{
	if(ResidueParamsDlgWX::dlg_open) return;
	MolSet* pmset = GetCurMolSet(); 
	if( pmset == NULL) return;
	ResidueParamsDlgWX* ptr_res_params_dlg = new ResidueParamsDlgWX( pmset, this );
	ptr_res_params_dlg->Show(TRUE);
}

void HaMainFrameWX::DoMolSetParamDialog(wxCommandEvent &event)
{
	if(MolSetParDlg::dlg_open) return;
	MolSet* pmset = GetCurMolSet(); 
	if( pmset == NULL) return;
	MolSetParDlg* ptr_mset_par_dlg = new MolSetParDlg( pmset, this );
	ptr_mset_par_dlg->Show(TRUE);
}

void HaMainFrameWX::DoAtomPropColorDialog(wxCommandEvent &event)
{
	if(AtomPropColorDlg::dlg_open) return;
	MolSet* pmset =  GetCurMolSet(); 
	if( pmset == NULL) return;
	HaMolView* pview = pmset->GetActiveMolView();

	AtomPropColorDlg* ptr_atom_prop_color_dlg = new AtomPropColorDlg( pview, this );
	ptr_atom_prop_color_dlg->Show(TRUE);
}

void HaMainFrameWX::DoEditGroupsDialog(wxCommandEvent &event)
{
    if(EditGroupsDlg::dlg_open) return;
	MolSet* pmset = GetCurMolSet(); 
	if( pmset == NULL) return;
	EditGroupsDlg* edit_grp_dlg = new EditGroupsDlg( pmset, 2, this );
	edit_grp_dlg->Show(TRUE);
}

void HaMainFrameWX::DoCrdSnapshotDialog(wxCommandEvent &event)
{
	if(CrdSnapshotDlg::dlg_open) return;
	MolSet* pmset = GetCurMolSet(); 
	if( pmset == NULL) return;
	CrdSnapshotDlg* crd_snap_dlg = new CrdSnapshotDlg( pmset, this );
	crd_snap_dlg->Show(TRUE);
}

void HaMainFrameWX::DoEditGeomDialog(wxCommandEvent &event)
{
	if(EditGeomDlgWX::dlg_open) return;
	MolSet* pmset = GetCurMolSet(); 
	if( pmset == NULL) return;
	EditGeomDlgWX* edit_geom_dlg = new EditGeomDlgWX( pmset, this );
	edit_geom_dlg->Show(TRUE);
}

void HaMainFrameWX::OnFindHbonds(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	if(pmset != NULL) 
	{
		MolEditor* p_mol_editor = pmset->GetMolEditor(true);
		p_mol_editor->CalcHBonds(pmset,true);
	}
}

void HaMainFrameWX::DoSolvateDialog(wxCommandEvent &event)
{
	if(SolvateDlgWX::dlg_open) return;
	MolSet* pmset = GetCurMolSet();
	if( pmset == NULL) return;
	SolvateDlgWX* ptr_solvate_dlg = new SolvateDlgWX( pmset, this );
    ptr_solvate_dlg->Show(TRUE);
}

void HaMainFrameWX::OnWrapIntoUnitCell(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet();
	if(pmset == NULL) return;
	MolEditor* p_mol_editor = pmset->GetMolEditor();
	int ires = p_mol_editor->WrapToUnitCell(pmset,pmset->per_bc);
	if( ires )
	{
		pmset->RefreshAllViews(RFRefresh | RFColour | RFApply);
	}
}

void HaMainFrameWX::OnCenterSolute(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet();
	if(pmset == NULL) return;
	MolEditor* p_mol_editor = pmset->GetMolEditor();
	int ires = p_mol_editor->CenterSoluteInSolvent(pmset);
}

void HaMainFrameWX::OnCenterMolInPBox(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet();
	if(pmset == NULL) return;
	MolEditor* p_mol_editor = pmset->GetMolEditor();
	int ires = p_mol_editor->CenterMolInPBox(pmset);
}

void HaMainFrameWX::OnAddMissingAtoms(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet();
	MolEditor mol_editor;
	if(pmset != NULL) mol_editor.AddMissingAtoms(pmset);
	
}

void HaMainFrameWX::OnAddPolarHydrogens(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	MolEditor mol_editor;
	if(pmset != NULL) mol_editor.AddPolarHydrogens(pmset);
	
}

void
HaMainFrameWX::OnAddHHybrid(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet();
	MolEditor mol_editor;
	if(pmset != NULL) mol_editor.AddHydrogensHybrid(pmset);
	
}

void HaMainFrameWX::OnAddHydrogens(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	MolEditor mol_editor;
	if(pmset != NULL) mol_editor.AddHydrogens(pmset);
}

void HaMainFrameWX::OnDelSelAtoms(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	if(pmset == NULL) return;
	AtomIteratorMolSet aitr(pmset);
	HaAtom* aptr;
	AtomGroup atlist;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if(aptr->Selected())atlist.InsertAtom(aptr);
	}
	pmset->DeleteAtoms(atlist);
	pmset->RefreshAllViews(RFRefresh | RFApply);	
}

void HaMainFrameWX::OnDelOvlpMols(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	if(pmset == NULL) return;
	AtomIteratorMolSet aitr(pmset);
	HaAtom* aptr;
	AtomGroup at_arr;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if(aptr->Selected())at_arr.InsertAtom(aptr);
	}
	MolEditor* p_mol_editor = pmset->GetMolEditor();
	p_mol_editor->DeleteOverlapMols(pmset,at_arr);
	pmset->RefreshAllViews(RFRefresh | RFApply);
}

void
HaMainFrameWX::DoBuildFilmDialog(wxCommandEvent &event)
{
	if(BuildFilmDlgWX::dlg_open) return;
	MolSet* pmset = GetCurMolSet(); 
	if(pmset == NULL) return;
	BuildFilmDlgWX* ptr_build_film_dlg = new BuildFilmDlgWX( pmset, this );
	ptr_build_film_dlg->Show(TRUE);
}

void
HaMainFrameWX::DoEditFragmDialog(wxCommandEvent &event)
{
	if(EditFragmDlgWX::dlg_open) return;

	MolSet* pmset = GetCurMolSet(); 
	if(pmset == NULL) return;

	EditFragmDlgWX* ptr_edit_fragm_dlg = new EditFragmDlgWX( pmset, this );
	ptr_edit_fragm_dlg->Show(TRUE);
}

void HaMainFrameWX::OnClearPicked(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	if(pmset == NULL) return;
	pmset->ClearPickedAtoms();
	pmset->RefreshAllViews(RFRefresh);
}

void HaMainFrameWX::OnSelAtomsInBoundBox ( wxCommandEvent &event )
{
	MolSet* pmset = GetCurMolSet(); 
	if(pmset == NULL) return;
	pmset->SelectAtomsInBoundaryBox();
}

void HaMainFrameWX::OnRevertAtomSelection ( wxCommandEvent &event )
{
	MolSet* pmset = GetCurMolSet(); 
	if(pmset == NULL) return;
	pmset->RevertAtomSelection();
}

void
HaMainFrameWX::OnDescribeSecStruct(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	if(pmset == NULL) return;
	pmset->DescribeSecStruct();
}

void
HaMainFrameWX::OnPrintHBonds(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	if(pmset == NULL) return;
	pmset->PrintHBonds();
}

void
HaMainFrameWX::OnSetAlphaHelix(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	if(pmset == NULL) return;
	MolEditor* p_mol_editor = pmset->GetMolEditor(true);
	p_mol_editor->SetAlphaHelix(pmset);
	pmset->RefreshAllViews(RFApply | RFRefresh);
}

void
HaMainFrameWX::DoNuclAcidDialog(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	if( NuclAcidDlgWX::dlg_open) return;
	if(pmset == NULL) return;
	NuclAcidMod* ptr_nucl_acid_mod = pmset->GetNuclAcidMod(true);
	if(ptr_nucl_acid_mod == NULL)
		return;
		
	NuclAcidDlgWX* nucl_acid_dlg = new NuclAcidDlgWX(ptr_nucl_acid_mod, this);
	nucl_acid_dlg->Show(TRUE);
}




// Display Menu
//  Display->Display Mode menu
void HaMainFrameWX::OnWireFrame(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->DisableSpacefill();
		pview->EnableWireframe(WireFlag,0);
		pview->SetRibbonStatus(False,0,0);
		pview->DisableBackbone();
		pview->UpdateThisView(RFRefresh);
	}
}

void HaMainFrameWX::OnBackBone(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->DisableSpacefill();
		pview->DisableWireframe();
		pview->SetRibbonStatus(False,0,0);
		pview->EnableBackbone(CylinderFlag,0.32);
		pview->UpdateThisView(RFRefresh);
	}
}
void HaMainFrameWX::OnSticks(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->DisableSpacefill();
		if( pmset->GetNAtoms() < 256 )
		{  
			pview->EnableWireframe(CylinderFlag,0.16);
		} 
		else 
			pview->EnableWireframe(CylinderFlag,0.32);
		pview->SetRibbonStatus(False,0,0);
		pview->DisableBackbone();
		pview->UpdateThisView(RFRefresh);
	}
}

void
HaMainFrameWX::OnSpheres(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->SetAtomScreenRadVdW();
		pview->DisableWireframe();
		pview->SetRibbonStatus(False,0,0);
		pview->DisableBackbone();
		pview->UpdateThisView(RFRefresh);
	}
}

void HaMainFrameWX::OnBallStick(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->SetAtomScreenRadVal(0.48);
		pview->EnableWireframe(CylinderFlag,0.16);
		pview->SetRibbonStatus(False,0,0);
		pview->DisableBackbone();
		pview->UpdateThisView(RFRefresh);
	}
}

void HaMainFrameWX::OnRibbons(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->DisableSpacefill();
		pview->DisableWireframe();
		pview->SetRibbonStatus(True,RibbonFlag,0);
		pview->DisableBackbone();
		pview->UpdateThisView(RFRefresh);		
	}
}

void HaMainFrameWX::OnStrands(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->DisableSpacefill();
		pview->DisableWireframe();
		pview->SetRibbonStatus(True,StrandFlag,0);
		pview->DisableBackbone();
		pview->UpdateThisView(RFRefresh);
	}
}
void
HaMainFrameWX::OnCartoons(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->DisableSpacefill();
		pview->DisableWireframe();
		pview->SetRibbonCartoons();
		pview->DisableBackbone();
		pview->UpdateThisView(RFRefresh);
	}
}
//  Display->Colours menu
void
HaMainFrameWX::OnMono(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->MonoColourAttrib(255,255,255);
		pview->UpdateThisView(RFRefresh);
	}
}
void
HaMainFrameWX::OnCPK(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->CPKColourAttrib();
		pview->UpdateThisView(RFRefresh);
	}
}
void
HaMainFrameWX::OnShapely(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->ShapelyColourAttrib();
		pview->UpdateThisView(RFRefresh);
	}
}

void
HaMainFrameWX::OnColGroups(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->GroupsColourAttrib();
		pview->UpdateThisView(RFRefresh);
	}
}

void
HaMainFrameWX::OnColResidues(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->ScaleColourAttrib( ResidueAttr );
		pview->UpdateThisView(RFRefresh);
	}
}

void
HaMainFrameWX::OnColChain(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
        pview->ScaleColourAttrib( ChainAttr );
		pview->UpdateThisView(RFRefresh);
	}
}
void

HaMainFrameWX::OnColTemper(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
        pview->ScaleColourAttrib( TempAttr );
		pview->UpdateThisView(RFRefresh);
	}
}

void
HaMainFrameWX::OnStruct(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->StructColourAttrib();
		pview->UpdateThisView(RFRefresh);
	}
}

//  Display->Options menu
void
HaMainFrameWX::OnSlab(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->SetUseSlabPlane(!pview->UseSlabPlane());
		pview->UpdateThisView(RFRefresh);
	}
}

void
HaMainFrameWX::OnHydrogen(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		int mask = NormAtomFlag;
		if( pview->HetaGroups )
		{
			mask |= HeteroFlag;
		}
		pview->Hydrogens = !pview->Hydrogens;
		pview->ReDrawFlag |= RFRefresh;

		if( pview->Hydrogens )
		{      
			pmset->SelectAtomsMask(mask|HydrogenFlag);
		} 
		else 
		{
			pmset->SelectAtomsMask(mask);
			pview->RestrictSelected();
		}			
		pview->UpdateThisView(RFRefresh);
	}
}
void
HaMainFrameWX::OnHetero(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		int mask = NormAtomFlag;
		if( pview->Hydrogens )
			mask |= HydrogenFlag;
		pview->HetaGroups = !pview->HetaGroups;
		pview->ReDrawFlag |= RFRefresh;

		if( pview->HetaGroups )
		{      
			pmset->SelectAtomsMask(mask|HeteroFlag);
		} 
		else 
		{
			pmset->SelectAtomsMask(mask);
			pview->RestrictSelected();
		}		
		pview->UpdateThisView(RFRefresh);
	}
}
void
HaMainFrameWX::OnSpecular(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->FakeSpecular = !pview->FakeSpecular;
		pview->UpdateThisView(RFColour);
	}
}

void
HaMainFrameWX::OnStereo(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		if( pview->UseStereo )
		{   
			pview->SetStereoMode(False);
		} 
		else 
		{
			pview->SetStereoMode(True);
		}
		pview->UpdateThisView(RFRefresh);
	}
}

void
HaMainFrameWX::OnLabels(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->LabelOptFlag = !pview->LabelOptFlag;
		pview->DefaultLabels(pview->LabelOptFlag);
		pview->UpdateThisView(RFRefresh);
	}
}

void HaMainFrameWX::OnDisplayPBox(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->DrawUnitCell = (!pview->DrawUnitCell); 
		pview->UpdateThisView(RFRefresh);
	}
}

void
HaMainFrameWX::DoViewParamDlg(wxCommandEvent &event)
{
	if(MolViewParDlg::dlg_open) return;
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		MolViewParDlg* ptr_view_par_dlg = new MolViewParDlg( pview , this);                          
		ptr_view_par_dlg->Show(TRUE);
	}
}

void
HaMainFrameWX::DoObject3DDialog(wxCommandEvent &event)
{
	if(Object3DDlgWX::dlg_open) return;
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		Object3DDlgWX* ptr_object3d_dlg = new Object3DDlgWX( pview, this );  
		ptr_object3d_dlg->Show(TRUE);
	}
}

// Display->Labels menu
void
HaMainFrameWX::OnLabelsAtomId(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->LabelOptFlag = !pview->LabelOptFlag;
		pview->DefaultLabels(pview->LabelOptFlag);
		if( pmset->GetNChains() >1 )
		{   
			pview->DefineLabels("%n%r:%c.%a");
		} 
		else if( pmset->GetNRes() > 1 )
		{   
			pview->DefineLabels("%n%r.%a");
		} 
		else 
			pview->DefineLabels("%a");	
		pview->UpdateThisView(RFRefresh);
	}
}

void HaMainFrameWX::OnLabelsAtomNames(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->DefineLabels("%a");	
		pview->UpdateThisView(RFRefresh);
	}
}

void HaMainFrameWX::OnLabelsAtomSeqNum(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->DefineLabels("%e%i");
		pview->UpdateThisView(RFRefresh);
	}
}

void HaMainFrameWX::OnLabelsAtomFFSymb(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->DefineLabels("%f");	
		pview->UpdateThisView(RFRefresh);
	}
}

void HaMainFrameWX::OnLabelsOff(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->DeleteLabels();
		pview->UpdateThisView(RFRefresh);
	}
}
// Display->Windows menu
void
HaMainFrameWX::OnWindowsCascade(wxCommandEvent &event)
{

}

void
HaMainFrameWX::OnWindowsTileVert(wxCommandEvent &event)
{

}

void
HaMainFrameWX::OnCenterViewSel(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet(); 
	HaMolView* pview = NULL;
	if(pmset != NULL) pview = pmset->GetActiveMolView();
	if(pview)
	{
		pview->CenterSelected();
		pview->UpdateThisView(RFRefresh);
	}
}

// Measure Menu

void HaMainFrameWX::OnShowAtomID(wxCommandEvent &event)
{ 
	HaMolView* pview = CurMolView;
	if(pview) pview->SetPickMode(PickIdent);
}

void HaMainFrameWX::OnMeasureDist(wxCommandEvent &event)
{
	HaMolView* pview = CurMolView;
	if(pview) pview->SetPickMode(PickDist);
}

void HaMainFrameWX::OnMeasureAngle(wxCommandEvent &event)
{
	HaMolView* pview = CurMolView;
	if(pview) pview->SetPickMode(PickAngle);
}

void HaMainFrameWX::OnMeasureDihed(wxCommandEvent &event)
{ 
	HaMolView* pview = CurMolView;
	if(pview) pview->SetPickMode(PickTorsn);
}



// ET Menu 

void HaMainFrameWX::OnEditRedox( wxCommandEvent &event )
{
    if(EditGroupsDlg::dlg_open) return;
	MolSet* pmset = GetCurMolSet();
	if(pmset != NULL)
	{
	   EditGroupsDlg* edit_grp_dlg = new EditGroupsDlg( pmset, 1, this );
	   edit_grp_dlg->Show(TRUE);
	}
}


void HaMainFrameWX::DoPathwaysDialog( wxCommandEvent &event)
{
	if(PathwaysDlgWX::dlg_open) return;

    MolSet* pmset = GetCurMolSet();
	if(pmset != NULL)
	{
	  ETCouplMod* et_coupl_mod = pmset->GetETCouplMod(true);
	  if( et_coupl_mod == NULL)
	   	 return;

	  PathwaysDlgWX* ptr_pathways_dlg = new PathwaysDlgWX( et_coupl_mod, this );
	  ptr_pathways_dlg->Show(TRUE);
	}
}

void HaMainFrameWX::DoETEffHamDialog( wxCommandEvent &event)
{
	if(ETEffHamDlgWX::dlg_open) return;

    MolSet* pmset = GetCurMolSet();
	if(pmset != NULL)
	{
	   ETCouplMod* ptr_et_coupl_mod = pmset->GetETCouplMod(true);
	   if(ptr_et_coupl_mod == NULL)
		  return;

	   ETEffHamDlgWX* ptr_et_ham_dlg = new ETEffHamDlgWX( ptr_et_coupl_mod, this );
	   ptr_et_ham_dlg->Show(TRUE);
	}
}


void HaMainFrameWX::OnResetETModule(wxCommandEvent &event)
{
    MolSet* pmset = GetCurMolSet();
	if(pmset != NULL)
	{
		ETCouplMod* ptr_et_coupl_mod= pmset->GetETCouplMod(false);
		if(ptr_et_coupl_mod != NULL)
			ptr_et_coupl_mod->Clear();	
	}
}

void HaMainFrameWX::DoScatterDialog(wxCommandEvent &event)
{
    if( ScatterDlg::dlg_open) return;
	
    MolSet* pmset = GetCurMolSet();
	if(pmset != NULL)
	{
	   HaScatterMod* ptr_sc_mod = pmset->GetScatterMod(true);
	   if(ptr_sc_mod == NULL)
	 	  return;
		
	   ScatterDlg* scatter_dlg = new ScatterDlg(ptr_sc_mod, this);
	   scatter_dlg->Show(TRUE);
	}
}

// QChem Menu

void HaMainFrameWX::DoQChemParamDialog(wxCommandEvent &event)
{
	if(QChemParDlgWX::dlg_open) return;

    MolSet* pmset = GetCurMolSet();
	if(pmset != NULL)
	{
	   HaQCMod* ptr_qc_mod = pmset->GetQCMod(true);
	   if(ptr_qc_mod == NULL)
		  return;
	   
	   QChemParDlgWX* ptr_qchem_par_dlg = new QChemParDlgWX( ptr_qc_mod, this );
	   ptr_qchem_par_dlg->Show(TRUE);
	}
}


void HaMainFrameWX::DoLoadQCDatDialog(wxCommandEvent &event)
{
	if(LoadQCDatDlgWX::dlg_open) return;
    MolSet* pmset = GetCurMolSet();
	if(pmset != NULL)
	{
	  HaQCMod* qcmod = pmset->GetQCMod(true);
	  LoadQCDatDlgWX* ptr_load_qcdat_dlg = new LoadQCDatDlgWX( qcmod , this );
	  ptr_load_qcdat_dlg->Show(TRUE);
	}
}

void HaMainFrameWX::DoWaveFunAnalDialog(wxCommandEvent &event)
{
	if(WaveFunAnalDlgWX::dlg_open) return;

    MolSet* pmset = GetCurMolSet();
	if(pmset != NULL)
	{
	   HaQCMod* ptr_qc_mod = pmset->GetQCMod(true);
	   if(ptr_qc_mod == NULL)
		  return;

	   WaveFunAnalDlgWX* ptr_wfun_anal_dlg = new WaveFunAnalDlgWX( ptr_qc_mod, this );
	   ptr_wfun_anal_dlg->Show(TRUE);
	}
}

void
HaMainFrameWX::OnCalcPolarGcontr(wxCommandEvent &event)
{
	HaTests::calc_polar_gcontr();
}

void
HaMainFrameWX::OnSaveGrpOperMat(wxCommandEvent &event)
{
	HaTests::save_grp_oper_mat(); 
}

void
HaMainFrameWX::OnCalcPolarContrF(wxCommandEvent &event)
{
	HaTests::calc_polar_contr_f();
}

void
HaMainFrameWX::OnReadPolarContr(wxCommandEvent &event)
{
	HaTests::read_polar_contr();
}

void
HaMainFrameWX::OnCalcBetaContr2idx(wxCommandEvent &event)
{
	HaTests::calc_beta_contr_2idx();
}

void HaMainFrameWX::OnReadBetaContr2idx(wxCommandEvent &event)
{
	HaTests::read_beta_contr_2idx(); 
}	

// MMech Menu

void HaMainFrameWX::DoMolMechDialog(wxCommandEvent &event)
{
    if( MolMechDlgWX::dlg_open) return;
	
    MolSet* pmset = GetCurMolSet();
	if(pmset != NULL)
	{
	   HaMolMechMod* ptr_mm_mod = pmset->GetMolMechMod(true);
	   if(ptr_mm_mod == NULL)
		  return;
		
	   MolMechDlgWX* mol_mech_dlg = new MolMechDlgWX(ptr_mm_mod,this);
	   mol_mech_dlg->Show(TRUE);
	}
}

void HaMainFrameWX::DoInterMolDialog(wxCommandEvent &event)
{
    if( InterMolDlgWX::dlg_open) return;
	
    MolSet* pmset = GetCurMolSet();
	if(pmset != NULL)
	{
	  HaInterMolMod* ptr_im_mod = pmset->GetInterMolMod(true);
	  if(ptr_im_mod == NULL)
		  return;
		
	  InterMolDlgWX* inter_mol_dlg = new InterMolDlgWX(ptr_im_mod,this);
	  inter_mol_dlg->Show(TRUE);
	}
}

void HaMainFrameWX::DoProtRedoxDialog(wxCommandEvent &event)
{
	if( ProtonRedoxDlg::dlg_open) return;

	MolSet* pmset = GetCurMolSet();
	if(pmset != NULL)
	{
		ProtonRedoxMod* p_prot_rdx_mod = pmset->GetProtonRedoxMod(true);
		if(p_prot_rdx_mod == NULL) return;

		ProtonRedoxDlg* prot_rdx_dlg = new ProtonRedoxDlg(p_prot_rdx_mod,this);
		prot_rdx_dlg->Show(TRUE);
	}
}


void HaMainFrameWX::DoContElectrDialog(wxCommandEvent &event)
{
	if(ElectrostDlgWX::dlg_open) return;

    MolSet* pmset = GetCurMolSet();
	if(pmset != NULL)
	{
	   ElectrostMod* electrost_mod = pmset->GetElectrostMod(true);
	   ElectrostDlgWX* ptr_electrost_dlg = new ElectrostDlgWX( electrost_mod, this );
	   ptr_electrost_dlg->Show(TRUE);
	}
}

void HaMainFrameWX::DoFlexModDialog(wxCommandEvent &event)
{
//#if defined(MOL_FLEX)
	if(wxMolFlex::dlg_open) return;

	wxMolFlex* mol_mech_dlg = new wxMolFlex(this);
	mol_mech_dlg->Show(TRUE);
//#endif
}

void HaMainFrameWX::DoEDDialog(wxCommandEvent &event)
{
	if(wxMolED::dlg_open) return;
	
    MolSet* pmset = GetCurMolSet();
	if(pmset != NULL)
	{
	   CollectCrdAnalMod* p_ccrd_mod = pmset->GetCollectCrdAnalMod(true);
	   if(p_ccrd_mod == NULL) return;
		
	   wxMolED* cluster_anal_dlg = new wxMolED(p_ccrd_mod,this);
	   cluster_anal_dlg->Show(TRUE);
	}
}

void HaMainFrameWX::DoPNPDialog(wxCommandEvent &event)
{
  //if(wxPNPFrame::dlg_open) return;
  MolSet* pmset = GetCurMolSet();
  
	if(pmset != NULL)
	{
//      PNPMod* pnp_mod = pmset->GetPNPMod(true);
//      wxPNPFrame* frame = new wxPNPFrame(this, -1, wxT( "Poisson Nernst Planck GUI" ),wxDefaultPosition, wxSize(800,480));
//      pnp_mod->SetPNPFrame(frame);
//      frame->SetPNPMod(pnp_mod);
//      frame->Show(true);
	}
}
void HaMainFrameWX::DoRFDialog(wxCommandEvent &event)
{
//   if(wxRFSmlPnl::dlg_open) return;
//   MolSet* pmset = GetCurMolSet();
//   
//   if(pmset != NULL)
//   {
//     PNPMod* pnp_mod = pmset->GetPNPMod(true);
//     wxString Title("PNP. Mol.Set: ");
//     Title+=pmset->GetName();
//     wxRFSmlPnl* frame = new wxRFSmlPnl(pnp_mod,this, -1, Title);
//     frame->Show(true);
//   }
}
//Test Menu
void
HaMainFrameWX::OnOpenVectorFieldInHaMolView(wxCommandEvent& event)
{
#ifdef PNP_DEPRECATED
	wxFileDialog dialog(this,
			wxT("Load Vector Field3D to HaMolView"),
			wxEmptyString,
			wxT("smth.gz"),
			wxT("Vector Field3D (*.gz)|*.gz|Vector Field3D (*.bin)|*.bin|DX (*.dx)|*.dx"),
			wxFD_OPEN|wxFD_FILE_MUST_EXIST);
	if (dialog.ShowModal() == wxID_OK)
	{
		wxString FileNameTree;
		FileNameTree=dialog.GetPath();
		VectorField3D* VField=new VectorField3D(FileNameTree.c_str());

		if(GetCurMolSet()==NULL)
			OnFileNew(event);

		int i=0;
		for(i=0;i<VField->Nelem;i++)
		{
			HaField3D* field=VField->GetHaField3D(i);
			wxString Name="";
			Name<<"VF["<<i<<"]";
			field->SetName( Name.ToStdString() );
			
			PlaneViewOfHaField3D* PlaneV=new PlaneViewOfHaField3D(field,Name.c_str());
			PlaneV->SetHideZeroValues(true);
			
			MolSet* pmset = GetCurMolSet();
			pmset->AddObject3D(PlaneV);
			HaMolView* pView= pmset->GetActiveMolView();
			if(pView)
			{
				pView->ReDrawFlag |= RFInitial;
				pView->InitialTransform();
			}
			pmset->RefreshAllViews(RFRefresh);
			
			wxFieldPlaneView *PV=new wxFieldPlaneView(PlaneV,pmset,this,-1,Name);
			PV->Show(true);
		}
	}
#endif
}
void 
HaMainFrameWX::OnOpenDielConstFromNindexInHaMolView(wxCommandEvent& event)
{
#ifdef PNP_DEPRECATED
	wxFileDialog dialog(this,
			wxT("Open Diel.Conts. Maps From NodeIndex File"),
			wxEmptyString,
			wxT("smth.gz"),
			wxT("Node Indexing (*.gz)|*.gz"),
			wxFD_OPEN|wxFD_FILE_MUST_EXIST);
	if (dialog.ShowModal() == wxID_OK)
	{
		wxString FileNameTree;
		FileNameTree=dialog.GetPath();
		NodeIndexing *NIndex=new NodeIndexing();
		NIndex->ReadFromFile(FileNameTree.c_str());

		if(GetCurMolSet()==NULL)
			OnFileNew(event);
		
		int i=0;
		for(i=0;i<3;i++)
		{
			HaField3D* field=NIndex->GetHaField3D(NodeIndexing::DielConst,NodeIndexing::Epsilon0);
			wxString Name="";
			Name<<"Epsilon["<<i<<"]";
			field->SetName( Name.ToStdString() );
			
			PlaneViewOfHaField3D* PlaneV=new PlaneViewOfHaField3D(field,Name.c_str());
			PlaneV->SetHideZeroValues(true);
			
			MolSet* pmset = GetCurMolSet();
			pmset->AddObject3D(PlaneV);
			HaMolView* pView= pmset->GetActiveMolView();
			if(pView)
			{
				pView->ReDrawFlag |= RFInitial;
				pView->InitialTransform();
			}
			pmset->RefreshAllViews(RFRefresh);
			
			wxFieldPlaneView *PV=new wxFieldPlaneView(PlaneV,pmset,this,-1,Name);
			PV->Show(true);
		}
		delete NIndex;
	}
#endif
}
void 
HaMainFrameWX::OnOpenDiffConstFromNindexInHaMolView(wxCommandEvent& event)
{
#ifdef PNP_DEPRECATED
	wxFileDialog dialog(this,
			wxT("Open Diff.Conts. Maps From NodeIndex File"),
			wxEmptyString,
			wxT("smth.gz"),
			wxT("Node Indexing (*.gz)|*.gz"),
			wxFD_OPEN|wxFD_FILE_MUST_EXIST);
	if (dialog.ShowModal() == wxID_OK)
	{
		wxString FileNameTree;
		FileNameTree=dialog.GetPath();
		NodeIndexing *NIndex=new NodeIndexing();
		NIndex->ReadFromFile(FileNameTree.c_str());

		if(GetCurMolSet()==NULL)
			OnFileNew(event);
		
		int i=0;
		for(i=0;i<NIndex->NIonsTypes;i++)
		{
			HaField3D* field=NIndex->GetHaField3D(NodeIndexing::DiffConst,(NodeIndexing::NodeIndexDescriptor)NIndex->IonField[i]);
			wxString Name="";
			Name<<"Diffusion["<<i<<"]";
			field->SetName( Name.ToStdString() );
			
			PlaneViewOfHaField3D* PlaneV=new PlaneViewOfHaField3D(field,Name.c_str());
			PlaneV->SetHideZeroValues(true);
			
			MolSet* pmset = GetCurMolSet();
			pmset->AddObject3D(PlaneV);
			HaMolView* pView= pmset->GetActiveMolView();
			if(pView)
			{
				pView->ReDrawFlag |= RFInitial;
				pView->InitialTransform();
			}
			pmset->RefreshAllViews(RFRefresh);
			
			wxFieldPlaneView *PV=new wxFieldPlaneView(PlaneV,pmset,this,-1,Name);
			PV->Show(true);
		}
		delete NIndex;
	}
#endif
}
void 
HaMainFrameWX::OnOpenQstFromNindexInHaMolView(wxCommandEvent& event)
{
#ifdef PNP_DEPRECATED
	wxFileDialog dialog(this,
			wxT("Open Qstatic Maps From NodeIndex File"),
			wxEmptyString,
			wxT("smth.gz"),
			wxT("Node Indexing (*.gz)|*.gz"),
			wxFD_OPEN|wxFD_FILE_MUST_EXIST);
	if (dialog.ShowModal() == wxID_OK)
	{
		wxString FileNameTree;
		FileNameTree=dialog.GetPath();
		NodeIndexing *NIndex=new NodeIndexing();
		NIndex->ReadFromFile(FileNameTree.c_str());
		
		if(GetCurMolSet()==NULL)
			OnFileNew(event);

		HaField3D* field=NIndex->GetHaField3D(NodeIndexing::Charge,NodeIndexing::ChargeMask);
		wxString Name="";
		Name<<"Qst";
		field->SetName( Name.ToStdString() );
		
		PlaneViewOfHaField3D* PlaneV=new PlaneViewOfHaField3D(field,Name.c_str());
		PlaneV->SetHideZeroValues(true);
		
		MolSet* pmset = GetCurMolSet();
		pmset->AddObject3D(PlaneV);
		HaMolView* pView= pmset->GetActiveMolView();
		if(pView)
		{
			pView->ReDrawFlag |= RFInitial;
			pView->InitialTransform();
		}
		pmset->RefreshAllViews(RFRefresh);
		
		wxFieldPlaneView *PV=new wxFieldPlaneView(PlaneV,pmset,this,-1,Name);
		PV->Show(true);
		delete NIndex;
	}
#endif
}
void
HaMainFrameWX::DoSmlPNPDialog(wxCommandEvent &event)
{
//   if(wxPNPSmlPnl::dlg_open) return;
//   MolSet* pmset = GetCurMolSet();
//   
//   if(pmset != NULL)
//   {
//     PNPMod* pnp_mod = pmset->GetPNPMod(true);
//     wxString Title("PNP. Mol.Set: ");
//     Title+=pmset->GetName();
//     wxPNPSmlPnl* frame = new wxPNPSmlPnl(pnp_mod,this, -1, Title);
//     frame->Show(true);
//   }
}

void HaMainFrameWX::DoContElectrPNPDialog( wxCommandEvent &event )
{
//   if(wxContElectSmlPnl::dlg_open) return;
//   MolSet* pmset = GetCurMolSet();
//   
//   if(pmset != NULL)
//   {
//     PNPMod* pnp_mod = pmset->GetPNPMod(true);
//     wxString Title("Continuum Electrostatics. Mol.Set: ");
//     Title+=pmset->GetName();
//     wxContElectSmlPnl* frame = new wxContElectSmlPnl(pnp_mod,this, -1, Title);
//     frame->Show(true);
//   }
}
void HaMainFrameWX::OnLoadManual(wxCommandEvent &event)
{
	StrVec args;
	StrVec prog_output;

	args.push_back(pApp->manual_main_page);

    HarlemApp::RunExternalProgram(RUN_BACKGROUND, pApp->html_browser,args, prog_output, FALSE);	
}

void HaMainFrameWX::OnMolConnect(wxCommandEvent &event)
{
	if(CurMolView)CurMolView->SetPickMode(PickMolConnect);	
}

void HaMainFrameWX::OnWorldConnect(wxCommandEvent &event)
{
	MolSet* pmset = GetCurMolSet();
	if(!pmset) return;
 	MoleculesType::iterator mol_itr;
	for(mol_itr= pmset->HostMolecules.begin(); mol_itr != pmset->HostMolecules.end(); mol_itr++ )
	{
		(*mol_itr)->SetConnected(true);
	}
	CurMolView->m_screen_transform = true;
	CurMolView->CalcRotCenter();
	PrintMessage(" Transform Dials Connected to the World ");
}

void MolSet::RefreshAllViews( long lHint  )
{
	lHint |= RFRefresh;
    if( canvas_wx != NULL)
    {
        canvas_wx->mol_view->ReDrawFlag |= lHint;
        canvas_wx->Refresh();
    }
}


void HaMainFrameWX::OnTestGraph1( wxCommandEvent &event )
{
	HaTests::test_graph_1();
}

void HaMainFrameWX::OnTestPython1( wxCommandEvent &event )
{
    HaTests::test_python_1();
}

void HaMainFrameWX::OnTestAvgPopFunc( wxCommandEvent &event )
{
	ProtonRedoxMod::TestCalcPopFun();
}

void HaMainFrameWX::OnSaveAmoebaTopFile1( wxCommandEvent &event )
{
	HaMolMechMod::TestSaveAmoebaTopFile1();
}

void HaMainFrameWX::OnSaveAmoebaTopFile2( wxCommandEvent &event )
{
	HaMolMechMod::TestSaveAmoebaTopFile2();
}

void HaMainFrameWX::OnTestFFT1( wxCommandEvent &event )
{
	MMDriverAmber::TestFFT1();
}

void HaMainFrameWX::OnTestMin1( wxCommandEvent &event )
{
    HaTests::test_min_1();
}

void HaMainFrameWX::OnTestMap( wxCommandEvent &event )
{
    multimap< int, double, less<int> > test_map;
    multimap< int, double, less<int> >::iterator itr1;
    multimap< int, double, less<int> >::iterator itr2;
    
    multimap<int, double, less<int> >::value_type p1(1,2.0);    
    test_map.insert(p1);
    multimap<int, double, less<int> >::value_type p2(1,3.0);    
    test_map.insert(p2);
    multimap<int, double, less<int> >::value_type p3(2,5.0);    
    test_map.insert(p3);
    multimap<int, double, less<int> >::value_type p4(3,7.0);    
    test_map.insert(p4);
    multimap<int, double, less<int> >::value_type p5(3,8.0);    
    test_map.insert(p5);

    for( int i = 1 ; i <= 4; i++)
    {
        PrintLog("For i = %d \n",i);
        itr1 = test_map.lower_bound(i);
        itr2 = test_map.upper_bound(i);

        PrintLog( "For i = %d:\n",i);

        if( itr1 == test_map.end())
        {
            PrintLog(" itr1 == test_map.end() \n");
        }
        else
        {
            PrintLog(" itr1 -> %6.3f \n", (*itr1).second);
        }
        if( itr2 == test_map.end())
        {
            PrintLog(" itr2 == test_map.end() \n");
        }
        else
        {
            PrintLog(" itr2 -> %6.3f \n", (*itr2).second);
        }
        PrintLog("\n");
    }    
}

void HaMainFrameWX::OnCalcTMatr1( wxCommandEvent &event )
{
    MolSet* pmset = GetCurMolSet();
    if(pmset == NULL) return;
    StmMod* stm_mod = pmset->GetSTMMod(true);
    pApp->StartWait();
    stm_mod->CalcTMatr1();
    pApp->EndWait();
}   

void HaMainFrameWX::OnTestQCMod1( wxCommandEvent &event )
{
    HaTests::test_qcmod_1(); 
}

void HaMainFrameWX::OnTestOper2( wxCommandEvent &event )
{
    HaTests::test_oper_2();
}

void HaMainFrameWX::OnTestOper1( wxCommandEvent &event )
{
    HaTests::test_oper_1(); 
}

void HaMainFrameWX::OnDumpOverlap2( wxCommandEvent &event )
{
    HaTests::dump_overlap2();
}

void HaMainFrameWX::OnDumpOverlap( wxCommandEvent &event )
{
    HaTests::dump_overlap();
}

void HaMainFrameWX::OnDumpGaussBCommon( wxCommandEvent &event )
{
    HaTests::dump_gauss_bcommon(); 
}

void HaMainFrameWX::OnDumpMolInfo( wxCommandEvent &event )
{
    HaTests::dump_mol_info();
}
