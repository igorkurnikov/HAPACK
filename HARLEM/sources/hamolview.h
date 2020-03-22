/*!  \file hamolview.h

    Classes to define 3D view of molecular set in HARLEM

    based on
    RasMol2 Molecular Graphics
    Roger Sayle, August 1995
    Version 2.6

    \author Igor Kurnikov
    \date 1998-2002

*/
#ifndef HAMOLVIEW_H
#define HAMOLVIEW_H


#include "hastl.h"
#include "habond.h"
#include "haatgroup.h"
#include "canvas3d.h"
#include "command.h"
#include "halinalg.h"

#include "hamolset.h"

const int LOAD_AMBER_RESTART_TIMER_ID = 100301;

class HaChain;

class HaMolecule;
class AtomExpr;
class Object3D;
class wxThread;

class Monitor
{
public:
    Monitor();
    virtual ~Monitor();
    HaAtom  *src;
    HaAtom  *dst;

    bool operator ==( const Monitor& ref) const;
    bool operator < ( const Monitor& ref) const;

    double dist;
    short col;
};

const int ResidueAttr = 0x00;
const int ChainAttr   = 0x01;
const int TempAttr    =	0x02;
const int ChargeAttr  = 0x03;

const int RibColInside   = 0x01;
const int RibColOutside  = 0x02;
const int RibColBoth     = 0x03;

const int DrawKnotFlag =  0x7e;
const int WideKnotFlag =  0x0e;

typedef void* PSItemPtr;

class wxImage;


class HaMolView 
//!  Class for 3D representation of a set of molecules
{
public:

    HaMolView();
	virtual ~HaMolView();

public:
	MolSet* GetMolSet(); //!< Get Molecular Set that the view displays
	
	int debug_level;

	void ResetView();

	int UseDepthCue;
	int UseStereo;
	int UseClipping;
	
	int SSBondMode,HBondMode;
	static int ZoneBoth;       //!< Flag to select bond only if both atoms ofthe bond are selected

	int UseTransparent;
	int UseOutLine;

	double StereoAngle;

	int DrawBoundBox,DrawAxes;
	int DrawDoubleBonds;
	int DrawUnitCell;

	double CartoonHeight; //!< Cartoon representation height in Angstroms
	int SolventDots;
	double ProbeRadius;    //!< Probe Radius in Angstroms

	int DrawDots,DrawLabels;
	int DrawMonitDistance;
	int DrawBetaArrows;
	
  int DrawObj3D;//draw *Obj3D (BoxObj3D,etc) //mikola 30july06
public:
	MolSet* host_mol_set;

public:
	BoxPartition HashTable;
 	
	int XOffset, YOffset;
	int ideltx, idelty;  //!< Current translation values along X and Y axis

  double CenX, CenY, CenZ;          //!< Coordinates of the Center of the Zone

	double CurRX; //!< Current value of rotation of the view around X axis  ( in [-1.1] interval)
	double CurRY; //!< Current value of rotation of the view around Y axis  ( in [-1.1] interval)
	double CurRZ; //!< Current value of rotation of the view around Z axis  ( in [-1.1] interval)

	double CurTX; //!< Current translation value of the view along X axis ( in [-1.1] interval)
	double CurTY; //!< Current translation value of the view along Y axis ( in [-1.1] interval)

	double CurZoom;      //!< Current Zoom value ( in [-1.1] interval)
	double CurSlabValue; //!< Current slab value ( in [-1.1] interval)

	double LastRX,LastRY,LastRZ;  //!< Previous values of CurRX,CurRY,CurRZ
	                              //!< Rotations around X,Y and Z axis 

	double LastTX, LastTY;        //!< Previous values of CurTX,CurTY - 
	                              //!< Translations along X,Y axis 

	int PointX; //!< X screen coordinate of the point of mouse click or release
	int PointY; //!< Y screen coordinate of the point of mouse click or release 
	int InitX;  //!< reference X screen coordinate of the point of mouse click in mouse motion tracking 
	int InitY;  //!< reference Y screen coordinate of the point of mouse click in mouse motion tracking

	static int FakeSpecular,SpecPower;
	int UseLabelCol;
	static int UseBackFade;
	static double Ambient;

 // Transformation matricies connecting absolute and screen coordinate system:
	HaMat_double Rot;  //!< Coordinates of x, y, z directional vectors in the screen(transformed) coordinate system 
	HaVec_double Orig; //!< Coordinates of the origin of the absolute coordinate system in the screen(transformed) c.s. 

	double Zoom;           //!< Zoom value

	double Scale;          //!< coef of transformation from the World Coordinates in Ang to Screen integer coordinates 
	double DScale;         //!< 1.0/(Max size of the displayed portion of the molecular set)
	                       
	int UseScreenClip;
	int m_screen_transform; //!< Flag to set rotations and translations only on the screen, real atomic coordinates are not changed

	int ReDrawFlag;
	int Hydrogens,HetaGroups;
	int DrawAtoms;
	int DrawBonds;
	int DrawRibbon;
	int DrawBestPath;
	int DrawContourSurf;
	int DrawSolidSurfaces;

protected:
	int StereoView;

	int FBufX,FBufY; //!< Save previous values of XRange and YRange

public:

	virtual int ExecuteCommand(CmdParser& cmd_pr);
	int ExecuteSetCommand(CmdParser& cmd_pr);
	int ExecuteColourCommand(CmdParser& cmd_pr);
	void ConnectObject(Object3D* pObj);

	int BroadcastCurrAtom();
	int FillCurrAtomRef(char* buf);

	int GetImageSize();
	int GetImageRadius();

	void SetUseSlabPlane(int new_UseSlabPlane) { pCanv->m_UseSlabPlane= new_UseSlabPlane; }
	int UseSlabPlane() { return pCanv->m_UseSlabPlane; }
	void SetSlabMode(int new_SlabMode ) { pCanv->m_SlabMode= new_SlabMode; }
	int SlabMode() { return pCanv->m_SlabMode; }
	void SetSlabValue(int new_SlabValue ) { pCanv->m_SlabValue= new_SlabValue; }
	int SlabValue() { return pCanv->m_SlabValue; }

	void SetSlabInten( int new_SlabInten)   { pCanv->m_SlabInten= new_SlabInten; }
	void SetSliceValue( int new_SliceValue) { pCanv->m_SliceValue= new_SliceValue; }
	int SliceValue() { return pCanv->m_SliceValue; }
	
	void SetZOffset(int new_ZOffset) { pCanv->m_ZOffset = new_ZOffset; }
	int ZOffset() { return pCanv->m_ZOffset; }

#if !defined(SWIG)
#if !defined(_WIN32) & !defined(TWIN32)
	static int OpenDisplay( int, int );
#endif
#endif

	void ClampShiftVal(int  ival, double  value ); //!< Shift value of the current translation, zoom or slab with restriction [-1,1]
	void WrapShiftVal(int iaxis, double value ); //!< Shift current rotation value around X,Y or Z axis, wraping them to [-1,1] interval

public:

	virtual void UpdateThisView(int lHint=0);
	virtual void RefreshScreen() {}
	int CreateImage();
	virtual void ClearImage(){}
	void ClearBuffers();
	virtual void TransferImage(){}
	virtual int PrintImage(){ return 1;}
	virtual int ClipboardImage(){ return 1;}
	void ReSizeScreen();
	void ReAllocBuffers();

	void BuildHashTable();  //!< Build HashTable - partion of selected atoms to grid bins
 
	void DrawFrame();        //!< Fill FBuffer and DBuffer which determine the image call RenderFrame
	void ResetRenderer();
	void SetStereoMode( int );
	void SetPickMode( int );
	void PickAtom( int, int, int );

// Transform functions
public:
	void InitialTransform(); //!< Arrange the molecular set within the window in standard orientation, and set rotation matrix to unity
	void PrepareTransform();
	void ApplyTransform();	
	void ResetTransform();	
	void CalcRotCenter(int sel_atoms = FALSE); //!< Compute center of rotation of connected molecules ( or selected atoms)

	void CenterSelected(); //!< Center Selected portion of the molecular set

	void GetTransfCoord(double x_abs, double y_abs, double z_abs,
		                double& x_tr, double& y_tr, double& z_tr) const; //!< Get tranformed coordinates (to get screen coordinates multiply by Scale and + ixadd(iyadd)
	
public:
	void SetAtomScreenRadVal( double rad ); //!< Set radii of atoms as they seen on the screen
	void SetRadiusTemperature();
	void SetAtomScreenRadVdW(); //!< Set image atom radii to van der Waals radii 
	void DefaultRepresentation();
	void DisableSpacefill();
	void EnableWireframe( int mask, double rad );
	void DisableWireframe();
	void EnableBackbone( int mask, double rad );
	void DisableBackbone();
	void SetHBondStatus( int enable, double rad );  //!< Set H-Bond display status 
	void SetSSBondStatus( int enable, double rad ); //!< Set SS-Bond display status 
	void SetRibbonCartoons();
    void SetRibbonStatus( int enable, int flag, double width );
    void SetTraceTemperature();
	
	static void SetMouseMode( int );
	static int MouseMode;
	static int UseHourGlass;

// Color Manipulation functions:

	vector<int> min_color_map; //!< Array of minimal set of colors to describe the image
	IntIntMap rev_color_map; //!< Map from colors to indexes in min_color_map

	static void RefreshColors(); //!< Reset colors for current display parameters for changes in fade and others
	
	int ComputeRevColourMap(); //!< Find minimal set of colors to define the image (min_color_map) 
	                            //!< and compute reverse color map rev_color_map
	                            //!< return number of colors in min_color_map

    static HaColor BackColor;  //!< Color of the background
    static HaColor BoxColor;   //!< Color of the Box
    static HaColor LabelColor; //!< Color of labels

	int ColorAtomsByProp( const std::string& str_prop, DValColorMap* p_col_map = NULL ); //!< Color selected atoms by specified property and color map
	void ColourBackNone();
	void ColourBondNone();
	void ColourHBondType();
	void ColourHBondNone();
	void ColourSSBondNone();
	void ColourMonitNone();
	void ColourBackAttrib( int r, int g, int b);
	void ColourBondAttrib( int r, int g, int b);
	void ColourHBondAttrib( int r, int g, int b);
	void ColourSSBondAttrib( int r, int g, int b);
	void ColourMonitAttrib( int r, int g, int b);
	void ColourDotsAttrib( int r, int g, int b);
	void ColourDotsPotential();
    void ColourRibbonNone( int flag );
    void ColourRibbonAttrib( int flag, int r, int g, int b );

	void ScaleColourAttrib( int attr );
	void MonoColourAttrib( int r, int g, int b );
	void CPKColourAttrib();
    static int ColorAtomCPK( HaAtom* aptr); //! Color Atom According to chemical element
	void GroupsColourAttrib();
	void ShapelyColourAttrib();	
	void AminoColourAttrib();
	void StructColourAttrib();
	void RigidClusterColourAttrib();
	
	static int IsCPKColour( HaAtom* aptr );
	static HaColorMap cpk_col_map; //!< map of atoms CPK colors

//! \name Labels manipulation functions

	void FormatLabel(HaChain* chain, HaResidue* group, HaAtom* aptr, 
		             const char* label, char* ptr );
	int DeleteLabels();
	void DefineLabels( const char* label );
	void DefaultLabels( int enable );
	void DisplayLabels();

	int LabelOptFlag;
	
//! \name Monitors functions:

	list<Monitor> MonitList;

	void DeleteMonitors();
	void AddMonitors( HaAtom* src, HaAtom* dst );
	void CreateMonitor( int src, int dst );
	void DisplayMonitors();

//! \name Display ET related quantities:

	void DisplayETBestPath(); //!< Display PATHWAYS Best Path 

//! \name Display Electrostatic Related quantities:

	void DisplayContourSurf();

//! \name Eigen Vector Display: 

	int AnimateEigenVector( HaVec_double& evec, AtomContainer* at_coll ); //!< Run an animation of an atom eigen vector in a separate thread
	int AnimateEigenVectorInternal( HaVec_double& evec, AtomContainer* at_coll ); //!< Internal function to run animation of an atom eigen vector
	void StopAnimation(); //!< Stop an animation
	wxThread* p_anim_thread; //!< Thread to run an animation
	int to_stop_animation;   //!< Flag to stop an animation


//! \name Atom Selection Functions:
	void RestrictSelected(); //!< Restrict view to only selected portion of the molecular set

//! \name Dot Surface functions:

	int TestDot( double x, double y, double z, bool solvent_access);
	void InitElemDots();
	void AddElemDots( int elem, int density );
	void FreeElemDots();

	void DeleteDotSurfaces();
	void CalculateDotSurface( int density);
	void DisplayDotSurfaces();
  void DisplayObj3D();//<mikola 30July06
//! \name Ribbons and Cartoons functions:

	void DisplayRibbon( HaChain * ); //!< Function called in RenderFrame to display secondary structure elements as ribbons

//! \name Image Output related functions:

protected:
	void WriteByte( int );
	void WriteLSBShort( int );
	void WriteMSBShort( int );
	void WriteMSBLong( unsigned int );	
	void WriteGIFCode(int code );
	void WritePPMWord(int i );	

    int FindDepth( PSItemPtr item, int type);
    void DepthSort( PSItemPtr* data, char* type, int count);	
	int ClipVectSphere( HaAtom * );
	int ClipVectBond( HaAtom *, HaAtom * );
	void WriteVectColour( int col );
	void WriteVectSphere( PSItemPtr *data, char *type,int index);
	void WriteVectWire( HaAtom* src, HaAtom* dst, int col, int dash );
	void WriteVectStick(HaAtom* src,HaAtom* dst, int col, int rad );

	int CountPSItems();
	void FetchPSItems( PSItemPtr *, char * );
	void WritePSItems( PSItemPtr *, char *, int );
	void WriteVectDots();
	void WriteVectLabels();
	void WriteVectMonitors();
	void FlushPICTBuffer();
	void FlushPICTPacket();
	void WritePICTCode(int val );
	void WritePICTData();
	void FlushIRISBuffer();
	void FlushIRISPacket();
	void WriteIRISCode(int val );
	void DetermineIRISSizes(int* rowstart,short* rowsize, int* min, int* max );
	void WriteIRISHeader(int* rowstart, short* rowsize, int min, int max );
	void WriteIRISData();

public:
	void WriteImageFile(const char* name, int type );
	int WritePPMFile( const char* name,int raw );   //!< Write the Molecular image into PPM file
	int WriteGIFFile( const char* name );           //!< Write the Molecular image into GIF file
	int WriteBMPFile( const char* name );           //!< Write the Molecular image into BMP file
	int WritePNGFile( const char* name );           //!< Write the Molecular image into PNG file
	int WriteJPEGFile( const char* name );          //!< Write the Molecular image into JPEG file
	int WriteTIFFFile( const char* name );          //!< Write the Molecular image into TIFF file
	int WritePCXFile( const char* name );           //!< Write the Molecular image into PCX file
	int WritePICTFile(const char* name );           //!< Write the Molecular image into APPLE PICT file
	int WriteIRISFile(const char* name );           //!< Write the Molecular image into IRIS RGB file

	int WriteVectPSFile(const char* name );         //!< Write the Molecular image into Vector Postscript File
	int WritePOVRayFile(const char* name);          //!< Write the POV-ray file with molecular structure description
	int WriteVRMLFile(const char* name);            //!< Write the VRML file with molecular structure description

	int WriteScriptFile(const char *name);          //!< Write RASMOL-type script describing structure
	int WriteMolScriptFile(const char* name);       //!< Write MolScript file with molecular structure description

public:

	void WriteVRMLDots();
	void WriteVRMLColour( int indent, int shade ); //!< Write VRML specification for colour
	void WriteVRMLTriple( double x, double y, double z ); //!< Write coordinate in VRML format
	void WriteVRMLAtoms(); //!< Write atoms in VRML format
	void WriteVRMLLine( int src, int dst, int shade, int* flag );
	void WriteVRMLWireframe();

protected:
	void WriteScriptAtoms();
	void WriteScriptBonds();
	void WriteScriptBackbone();
	void WriteScriptRibbons();
	void WriteScriptLabels();
	void WriteScriptMonitors();
	void WriteScriptHBonds( char* obj );

public:

#if !defined(HA_NOGUI)
    int SetWXImage(wxImage& wx_image); //!< Set wxImage with the current image data
#endif

	Canvas3D* pCanv;  //!< Pointer to the canvas object that contains the PixelMap buffer to plot different objects 
	
	void DisplaySpaceFill();   //!< Plot atoms as spheres in the buffer
	void DisplayWireframe();   //!< function to plot wireframe molecular representation
	void DisplayCylinder(int x1,int y1,int z1,
						 int x2,int y2,int z2,
						 int c1,int c2,int rad);
	void DisplayDoubleBonds();
    void DisplayBackbone();     //!< Plot the molecule backbone 
	void DisplayHBonds();  //!< Display H-bonds
    void DisplaySSBonds(); //!< Display SS-bonds
	void DisplayBoxes();   //!< Plot boundary box in the buffer
	void DisplayOnScreenInfo(); //!< Print info on the display screen
	void DisplayPickedAtoms();  //!< Display circles around picked atoms 
	void RenderFrame();   //!< The Main function to plot object in the buffer pixelmap before output to the display

	void TestAtomProximity( HaAtom* ptr, int xpos, int ypos); //!< Set Atom as the Picked atom if it is closest so far to the (xpos,ypos) point
	void IdentifyAtom( int xpos, int ypos );               //!< Find the atom the mouse points to 
	void InitializeTables();
	void InitializeRenderer();

};

//
// Macros below require definition
// MoleculesType::iterator mol_itr;
//

#define ForEachMol_VIEW     for(mol_itr=(GetMolSet()->HostMolecules).begin(); mol_itr != (GetMolSet()->HostMolecules).end(); mol_itr++)

const int PickNone    =     0x00;
const int PickIdent   =     0x01;
const int PickDist    =     0x02;
const int PickAngle   =     0x03;
const int PickTorsn   =     0x04;
const int PickLabel   =     0x05;
const int PickMonit   =     0x06;
const int PickCentr   =     0x07;
const int PickMolConnect =  0x08;

const int DefaultWide = 480;
const int DefaultHigh = 480;

const int RFRotateX  = 0x0001;
const int RFRotateY  = 0x0002;
const int RFRotateZ  = 0x0004;
const int RFZoom     = 0x0008;
const int RFTransX   = 0x0010;
const int RFTransY   = 0x0020;
const int RFTransZ   = 0x0040;
const int RFSlab     = 0x0080;
const int RFReSize   = 0x0100;
const int RFColour   = 0x0200;
const int RFRefresh  = 0x0400;
const int RFPoint1   = 0x1000;
const int RFPoint2   = 0x2000;

const int RFTrans    = 0x0070;
const int RFRotate   = 0x0007;
const int RFApply    = 0x017F;
const int RFDials    = 0x00FF;
const int RFMagnify  = 0x0108;
const int RFInitial  = 0x01FF;
const int RFPoint    = 0x3000;


const int MMRasMol   = 0x00;
const int MMInsight  = 0x01;
const int MMQuanta   = 0x02;


const int ViewLeft   =   0;
const int ViewRight  =   1;
 
unsigned int isqrt( unsigned int );

#define ZValid_v(z)     ((!UseSlabPlane() ) || ((z) < SlabValue()))
#define XValid_v(x)     (((x)>=0)&&((x)< pCanv->View.xmax))
#define YValid_v(y)     (((y)>=0)&&((y)< pCanv->View.ymax))

#ifdef HAMOLVIEW_CPP

HaMolView*  CurMolView;

#else

extern HaMolView* CurMolView;

#endif


#endif // !defined HAMOLVIEW_H 


