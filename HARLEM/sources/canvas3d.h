/*!  \file canvas3d.h

    Classes to plot 3D images in a window in HARLEM
 
    \author Igor Kurnikov 
    \date 1998-2008

*/

#ifndef CANVAS3D_H
#define CANVAS3D_H

#include "haconst.h"
#include "hastring.h"

// These values set the sizes of the sphere rendering
// tables. The first value, maxrad, is the maximum
// sphere radius and the second value is the table
// size = (maxrad*(maxrad+1))/2 + 1
//


// #define MAXRAD    120   256   
// #define MAXTABLE  7261  32897 

typedef uint_4 ColorVal;

//#if defined(SIXTEENBIT)
//   typedef unsigned short ColorVal;

#define MAXRAD    255
#define MAXTABLE  32641


   class Knot 
   //! Element of a ribbon or cartoon
   {
   public:
	   int px, py, pz;      //!< Spline Control Co-ordinate - coordinates of corners will be px-wx-dx,py-wy-dy,etc
	   int tx, ty, tz;      //!< Spline Direction Vector    
	   int hnx, hny, hnz;   //!< Horizontal Normal Vector   
	   int vnx, vny, vnz;   //!< Vertical Normal Vector   
	   int dx, dy, dz;      //!< Ribbon Height Vector       
	   int wx, wy, wz;      //!< Ribbon Width Vector        
	   char hinten;         //!< Horizontal Intensity - for plotting color of the side surfaces of the cartoon elements
	   char vinten;         //!< Vertical Intensity   - for plotting upper and lower surface of the cartoon (and ribbon) elements)     
	   short hsize;         //!< Horizontal Vector Length   
	   short vsize;         //!< Vertical Vector Length     
	   short wide;          //!< Ribbon Width               
   };
   
#define MAXVERT 10
   
   class Edge 
   {
   public:
	   int dx,dz;
	   int dr,dg,db;
	   int x,z;
	   int r,g,b;
   };

   class Vert
   //! Vertice of a polygon for plotting
   {
   public:
       int x, y, z; //!< screen cordinates of the vertice
       int inten;   //!< Colour of the given vertice (including light intensity changes due to the turn of the surface element 
   };

   class Poly
   //! Polygon class for plotting
   {
   public:
	   Vert v[MAXVERT]; //!< Verticies of the polygon
	   int count; //!< Number of verticies in the polygon 
   };


   class ArcEntry
   {
   public:
                short dx,dy,dz;
                short inten;
                int offset;
   };



typedef struct {
        ColorVal  *fbuf;     /* Pixels of the image       */
        short  *dbuf;        /* z-coordinates of the points of the image used so  
		                        points with smaller z will not be plotted */
        int xmax, ymax;
        int yskip;
    } ViewStruct; /*  Structure with picture to be filled by pixel utils  */

#define ZValid(z)     ((!m_UseSlabPlane) || ((z) < m_SlabValue))   
#define XValid(x)     (((x)>=0)&&((x)<View.xmax))   
#define YValid(y)     (((y)>=0)&&((y)<View.ymax))    

const int LutSize = 20000;

const double DefaultAmbient = 0.4;  //!< Intensity of light coming from everywhere relative to intensity of front light
const int    ColourDepth    = 32;   //!< the number of intermediate colors to model fading intensity of reflected light
const int    ColourMask     = 31;   //!< ColourMask = ColourDepth-1

const int SlabReject  =  0x00;
const int SlabHalf    =  0x01;
const int SlabHollow  =  0x02;
const int SlabFinal   =  0x03;
const int SlabClose   =  0x04;
const int SlabSection =  0x05;

const int ColBits    =  24;

#define RComp(x)   (((x)>>16)&0xff)  //!< Red Component of color x
#define GComp(x)   (((x)>>8)&0xff)   //!< Green Component of color x
#define BComp(x)   ((x)&0xff)        //!< Blue Component of color x

class HaColor
//! \brief public specifying Color in HARLEM
//!
//! Derivatives of this color (lower intensity are in Lut[] array)
//!   
{
public:
	HaColor(int r, int g, int b);    //!< Constructor from RGB 
	HaColor(ColorVal cval);          //!< Set from unsigned int that is currently define pixel in FBuffer  
	HaColor(const HaColor& color);
	~HaColor();

	bool operator==( const HaColor& ref) const;	

	int SetColor(int r, int g, int b); 
	int SetColor( const std::string& col_str ); //!< Set Color By Name
	static bool ColorNameExist(const std::string& col_str ); //!< Check if color name exist 
	static ColorVal GetColorVal( const std::string& col_str ); //!< Get Color Value by Name
	static std::string GetColorName( int r, int g, int b ); //!< Get Color Name by RGB values
	static ColorVal GetColorVal(int r, int g, int b);  //!< Get Color Value corresponding to RGB values
    static int RegisterColor(int r, int g, int b);     //!< Set entry of the color and its shades in Lut[], cval_idx_map and used_colors 
    static ColorVal GetPureColorByCIdx(int cidx);       //!< Get Color from color index in lut[]
    static void GetPureRGBbyCIdx(int cidx, int& r, int& g, int& b ); //!< Get RGB values for a color with index cidx
	static void GetRGBFromColVal(ColorVal cval, int& r, int& g, int& b ); //!< Compute RGB components of color given by ColorVal
	static int GetNumColors(){ return used_colors.size(); }  //!< Get Number of registerd colors
	static int RegisterColorName(const std::string& col_name, int r, int g, int b); //!< Register Color Name
	

    unsigned char r;  //!< red   intensity of the color
    unsigned char g;  //!< green intensity of the color
    unsigned char b;  //!< blue  intensity of the color
	
	ColorVal cval;    //!< Pixel color value. RGB values of this color to be extracted by RComp(x),GComp(x),BComp(x)
	int  cidx;        //!< position of the this colour in Lut[] array 

	static IntIntMap cval_idx_map;        //!< Map of main color values as ColorVal numbers to indecies in Lut[] array

	static std::map<ColorVal,std::string, less<ColorVal> >    col_name_map; //!< map of color names  from color values
	static std::map<std::string, ColorVal,less<std::string> > name_col_map; //!< map of color values from color names
	
	static vector<ColorVal> used_colors;  //!< Array of used(initiated) colors
	static int color_names_init_flag;  //!< flag to indicate that color names have been initiated

 static void GetRGBbyCIdxFloat(int cidx, float *fcol)
 {
   int ir, ig, ib;
   HaColor::GetPureRGBbyCIdx(cidx,ir, ig, ib);
   fcol[0]=(float)ir/255.0;
   fcol[1]=(float)ig/255.0;
   fcol[2]=(float)ib/255.0;
 }
};

class HaColorMap
//! Class to represent an array of colors to represent some property 
{
public:
	HaColorMap();
	virtual ~HaColorMap();

	std::vector<HaColor> colors;

	HaColor& GetColorByIdx(int idx) { return colors[idx];  } //!< zero based index acces to color array
	int GetNColors() { return colors.size(); } 

	int AddColor(int r, int g, int b); //!< Add Color entry into the map
	
	int AddUniformRange(int num_cols, HaColor& col1, HaColor& col2); //!< Add colors interpolating between two colors  
	
	virtual int GetIdxByProp(void* prop) { return 0; } //!< function to map a property to index 
};

class StrColorMap: public HaColorMap
//! Class for mapping string to colors
{
public:
	StrColorMap() {}
    virtual ~StrColorMap() {}

	HaColor* GetColorForStr(const char* str);
	int AddStrColorPair( const char* str, int r, int g, int b);

	StrIntMap str_idx_map;
	IntIntMap cval_idx_map;
};

class IValColorMap: public HaColorMap
//! Class for mapping string to colors
{
public:
	IValColorMap();
    virtual ~IValColorMap();

	HaColor* GetColorForInt(int ival);
	int AddIntColorPair( int ival, int r, int g, int b);

	IntIntMap int_idx_map;
	IntIntMap cval_idx_map;
};

class DRangeColorMap: public HaColorMap
//! Class for mapping continjuos range of real values to colors
{
public:
	DRangeColorMap();
    virtual ~DRangeColorMap();

	virtual HaColor* GetColorForVal(double val);

	double min_val;
	double max_val;
};

class DValColorMap: public HaColorMap
//! Class for specified ranges of real values to colors
{
public:
	DValColorMap();
    virtual ~DValColorMap();

	virtual HaColor* GetColorForVal(double val);

	int SaveToTxtFile(const std::string& fname);     //!< Save Color Map to   text file 
	int LoadFromTxtFile(const std::string& fname);   //!< Load Color Map from text file

	void Clear(); //!< Clear Color Map
	bool IsValid() const; //!< Check if the color map is valid
	int AddColorAndMinVal(int r, int g, int b, double val_min, int n_interp_col = 0); //!< Add Color entry into the map and minimal value for the color (ignored for the first added color) 

	std::vector<double> min_values; // array of minimal property values (N Colors() - 1) of number of colors 
};



class Canvas3D
//!  Class to define Window to Plot 3D Images 
{
public:

	Canvas3D();
	~Canvas3D();

	void resize(int new_XRange, int new_YRange);
	ColorVal* AllocImage();
	void DestroyImage();
	short* AllocDBuffer();
	void DestroyDBuffer();

	ColorVal  *FBuffer;  //!< Arrays of Pixels With the image
	short  *DBuffer;     //!< Arrays of maximal Z-coordinates of the objects

	int FBClear; //!< Flag to indicate if FBuffer values set to zero 
	int DBClear; //!< Flag to indicate if DBuffer values set to zero 

	ViewStruct View; //!< Structure to describe active portion of the image to be filled
	                 //!< by plotting functions

	static ColorVal Lut[LutSize]; //!< palette pixel entries - set by a call to DefineColourMap()
	static int empty_lut_idx;            //!< current empty entry in Lut[] array

	int m_FontSize;
	int m_SplineCount;
	int m_UseSlabPlane;
	int m_SlabValue;
	int m_SlabMode;
	int m_SlabInten;
	int m_SliceValue;
	int m_ImageSize;    //!< Slab thickness of the Image
	int m_ImageRadius;  //!< 1/2 m_ImageSize
	int m_ZOffset;      //!< Slab position

	static unsigned char *LookUp[MAXRAD]; //!< LookUp shows starting point of 
	static unsigned char Array[MAXTABLE]; //!< Sqrt(rad**2 - i*i) array 
	                                      //!< in Array 
	                                      //!< MAXTABLE= MAXRAD*(MAXRAD-1)/2+1

	static unsigned int ColConstTable[MAXRAD];
	static unsigned int *ColConst;

	int XRange() const { return m_XRange; }
	int YRange() const { return m_YRange; }
	int HRange() const { return m_XRange/2; }
	int WRange() const { return m_YRange/2; }
	int Range () const { return MinFun(m_XRange,m_YRange); }

	void PlotDeepPoint( int x, int y, int z, int col);
	void ClipDeepPoint( int x, int y, int z, int col );                              
	//! Draw line between points (screen coords) with color indexes col1 and col2 - no depth cueing and no clipping
	void DrawTwinLine( int x1,int y1,int z1,int x2,int y2,int z2, int col1,int col2 );
	//! Draw line between points (screen coords) with color indexes col1 and col2 - clipping but no depth cueing
	void ClipTwinLine( int x1,int y1,int z1,int x2,int y2,int z2, int col1,int col2 );
	//! Draw line between points (screen coords) with color indexes col1 and col2 - depth cueing but not clipping
	void DrawTwinVector( int x1,int y1,int z1,int x2,int y2,int z2, int col1,int col2 );
	//! Draw line between points (screen coords) with color indexes col1 and col2 - clipping and depth cueing 
	void ClipTwinVector( int x1,int y1,int z1,int x2,int y2,int z2, int col1,int col2 );
	//! Draw dashed line between points (screen coords) with color indexes col1 and col2 - clipping and depth cueing 
	void ClipDashVector( int x1,int y1,int z1,int x2,int y2,int z2, int col1,int col2 );

	void DrawCylinder( int x1,int y1,int z1,int x2,int y2,int z2,int c1,int c2,int rad);
	void ClipCylinder( int x1,int y1,int z1,int x2,int y2,int z2,int c1,int c2,int rad );
	void DashRibbon( Knot* src, Knot* dst, int col1, int col2 );
	void StrandRibbon( Knot* src, Knot* dst, int col1, int col2 );
	void SolidRibbon2( Knot* src, Knot* dst, int col1, int col2 );
	void SolidRibbon( Knot* src, Knot* dst, int col );
	void RectRibbon ( Knot* src, Knot* dst, int col );
	void DrawSphere( int x,int y,int z,int rad,int col );
	void ClipSphere( int x,int y,int z,int rad,int col );

	void SetFontSize( int size );
	void DisplayTextString( int x, int y, int z, const char* label, int col  ); //!< Display string in a given position and with given color

protected:

	int m_XRange, m_YRange; //!< Image Dimensions
	char FontDimen[23];   

	int OutCode(int x, int y, int z);
	void PlotPoint(int x,int y,int z,int col);
	void ClipPoint(int x, int y,  int z, int col);
	void ClipVector(int x1,int y1,int z1,
					int x2,int y2,int z2,int col);
	void ClipLine(int x1,int y1,int z1,int x2,int y2,int z2,int col);
public:
//	void OutLinePolygon( Poly* p );
//	void DrawPolygon( Poly* p );
	void ClipPolygon( Poly* p, double transp = 0.0); //!< Function to plot a piice of surface bourders by a polygon
protected:
	int TestSphere( int x, int y, int z, int rad );
	void DrawArcAc(short* dbase, ColorVal* fbase, int z, int c);
	void DrawArcDn( short* dbase, ColorVal* fbase, int z, int c);
	void DrawCylinderCaps( int x1,int y1,int z1, 
							  int x2, int y2, int z2, 
							  int c1, int c2, int rad );
	int TestCylinder( int x1, int y1, int z1,
					  int x2, int y2, int z2, int rad );
	void ClipArcAc( short* dbase, ColorVal* fbase,
					int x, int y, int z, int c);
	void ClipArcDn( short* dbase, ColorVal* fbase,
					int x, int y, int z, int c);
	void ClipCharacter( int x, int y, int z, int glyph, int col );


};



#endif /* !CANVAS3D_H */
