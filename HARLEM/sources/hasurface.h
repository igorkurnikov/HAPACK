/*! \file hasurface.h

   Classes to define surfaces in HARLEM 
 
   \author  Igor Kurnikov

   \date 1999-2002
*/
#if !defined(HASURFACE_H)
#define HASURFACE_H


#include "hastl.h"
#include "haconst.h"
#include "halinalg.h"
#include "vec3d.h"
#include "object3d.h"
#include "haobject.h"

#define MAXVARS 10
#define MAXTIMES 10

const int MAX_ISO_VERTS = 1200000;
const double VERTEX_SCALE  = 10000.0;

class AtomContainer;

//! Class for a scalar field on a rectangular grid
class HaField3D:public HaObject
{
public:

    HaField3D(float *fvec=NULL,int new_Nx=0, int new_Ny=0, int new_Nz=0, bool deligate_control=false);
//	HaField3D();
 
	virtual ~HaField3D();

	virtual void clear();  //!< erase the field and deallocate the memory it takes
	void FillZeros();      //!< Fill the field with zeros
	void FillValues(float val);//!< Fill the field with some value
	void MultiplyByValues(float val);//!< MultiplyBuValues

    int SaveGZ(const char* fname,int Columns=2);//!<Save the field to gzip file
    int LoadGZ(const char* fname);//!<Load the field from gzip file
	int SaveToFile(const char* fname,      int binary=TRUE);       //!< Save the field to a file
	int RestoreFromFile(const char* fname, int binary=TRUE);  //!< Restor the field from a file

	virtual void copy_from(const HaField3D& ref_field);

	virtual bool SetDimensions(int new_Nx, int new_Ny, int new_Nz); //!< Set the grid dimensions

	int GetNx() const { return m_Nx ; } //!< return number of grid points along X axes 
	int GetNy() const { return m_Ny ; } //!< return number of grid points along Y axes 
	int GetNz() const { return m_Nz ; } //!< return number of grid points along Z axes 
	
	double GetXmin() const { return xmin; } //!< return Min X coordinate of the grid
	double GetYmin() const { return ymin; } //!< return Min Y coordinate of the grid
	double GetZmin() const { return zmin; } //!< return Min Z coordinate of the grid
	double GetXmax() const { return xmax; } //!< return Max X coordinate of the grid
	double GetYmax() const { return ymax; } //!< return Max Y coordinate of the grid
	double GetZmax() const { return zmax; } //!< return Max Z coordinate of the grid

	double GetXstep() const { return xstep; } //!< return the length of the step along X axes 
	double GetYstep() const { return ystep; } //!< return the length of the step along Y axes 
	double GetZstep() const { return zstep; } //!< return the length of the step along Z axes 
	
	void GetMinMaxValue(float *ValueMin,float *ValueMax) const;

	bool SetGridCornersCoord(double xmin_new, double ymin_new, double zmin_new,  
		                    double xmax_new, double ymax_new, double zmax_new); //!< set min and max expansions of the grid
  bool ShiftGridCorners(double xsht, double ysht, double zsht);//!< shift min and max expansions of the grid by Rsht, mikola 11-19-06
  bool SetCenterAsZero(float scale);//!< Set center of field to (0,0,0), scale grids/Born, mikola 12-06-06

	int GetLinIdx(int ix, int iy, int iz); //!< Get linear index of the point of the grid (0 - based)

	float* GetFieldPtr() { return m_field_data.begin(); } //!< Get the pointer to the beginning of the grid data

	//! get cartesian coordinates of the point of the grid with grid indexes (ix,iy,iz) 
	bool GetXYZ(float& x, float& y , float& z, const int ix, const int iy, const int iz);
	Vec3D GetGridPtCrd(int ix, int iy, int iz); //!< Get Cartesian Coordinates of the grid point as Vec3D

	float* GetValPtr( int ix, int iy, int iz); //!< get pointer to the value of the field at the grid point 
	float GetValue( int ix, int iy, int iz);   //!< get Value of the field at the grid point (0 - based indexing)
	                                           //!< first point of the grid is (0,0,0)

	void SetValue(int ix, int iy, int iz, float val); //!< Set Value of the grid point

	int GetClosestGridPoint( double x, double y, double z, int& ix, int& iy, int& iz); //!< get a grid vertex closest to 
	                                                                                   //!< the given point (x,y,z)

    double CalcLinInter(double V0,double V1,double x);
	double GetInterpolValAtPoint( double x, double y, double z); //!< Interpolate the value of the field in the point

	bool grid_to_xyz_float(const int numverts, const float* vr, const float* vc, const float* vl,
	                       float* xyz_coord);  //!< convert vertices from grid to (x,y,z) coordinates

  float CompareHaField3D(HaField3D *CompareWith,float Prec);
protected:
	float xmin, xmax;    //!< grid maximum and minimum expansions along X  axes
	float ymin, ymax;    //!< grid maximum and minimum expansions along Y axes
	float zmin, zmax;    //!< grid maximum and minimum expansions along Z axes
	
	float xstep;   //!< grid step along X axes
	float ystep;   //!< grid step along Y axes
	float zstep;   //!< grid step along Z axes

	int m_Nx, m_Ny, m_Nz;        //!< Grid dimensions along X,Y and Z axes 
	HaVec_float m_field_data;    //!< Linear array of the field data

};


class HaNonLocField3D_2 : public HaField3D
{
public:
	HaNonLocField3D_2();
	virtual ~HaNonLocField3D_2();

	bool SetDepth(int new_depth);
	int GetDepth() const;

	virtual void clear();
	virtual bool SetDimensions(int new_Nx, int new_Ny, int new_Nz);

	float GetValue_nloc(int ix, int iy, int iz,
		                int ir_shift, int ic_shift, int il_shift);

protected:
    int depth;
	int nl_cube_size;
	HaMat_float m_nonl_field_data;

};

class ValAtPoint
{
public:
	ValAtPoint() {}
	ValAtPoint(short ir_new,short ic_new, short il_new, double new_val);

	virtual ~ValAtPoint() {}
	
	short ix;
	short iy;
	short iz;
	double val;
};

class HaNonLocField3D : public HaField3D
{
public:
	HaNonLocField3D();
	virtual ~HaNonLocField3D();

	std::vector< std::list<ValAtPoint> > fvals;

	int SaveField(const std::string& fname);

protected:

};

class HaSurface
//! class to define surface build from triangles
{	
public:

	HaSurface();
	virtual ~HaSurface();
	virtual void clear();
    double GetSurfVolume() const { return surf_volume ;} ;

	int     valid;       //!< valid/initialized surface flag 
	float   isolevel;    //!< the isolevel of the surface 
	
	bool calc_isosurf( HaField3D* field, float iso_level); //!< Calculate isosurface  

	static void main_march( float *ptGRID, int NC, int NR, int NL,
                 int LOWLEV,
                 float GLEV, float ARX, float ARY, float ARZ,
                 int NVERTS, float *VX, float *VY, float *VZ,
                 float *NX,float *NY, float *NZ,
                 int NPTS, int *VPTS, int*IVERT, int *IPTS, int *IPOLY,
                 int *ITRI);

	//! Molecular Surface Types enum:
	enum MolSurfaceType { VDW_SURF=0,     //!< Van der Waals Surface 
		                  SACCESS_SURF=1, //!< Solvent Access Surface
						  SEXCL_SURF=2    //!< Solvent Exclusion surface
	                    };  

    static int CalcMolSurf( HaSurface* sptr, int surf_type, float solv_rad, AtomContainer& at_coll ); //!< Compute Molecular surface for atom collection

	static class GEPOLParams
	//! parameters to build molecular surface using GEPOL93 algorithm
	{
	public:
		GEPOLParams();
		~GEPOLParams();

		float rmin_axx_sph;    //!< minimal radius (Ang) of axxiliary spheres to build excluded solvent surface (smaller the better) (default = 0.5)
		float overlap_axx_sph; //!< overlapping factor (from 0 to 1) larger the better  (default = 0.8)
		int   ndiv_sph;        //!< level of subdivision of the original triangles (from 1 to 5) forming sphere surface (default = 3)
	} gepol_prm;
	
//	static Vec3DValArray vrt_std;  //!< verticies of main triangles for a standard triangulation of a sphere
//    static Vec3DValArray tcnt_std;  //!< positions of centers of spherical triangles for a triangulation of a sphere (for ndiv_sph)

	HaMat_float verts;      //!< array of vertices coordinates
	HaMat_float norms;      //!< array of normals
	HaMat_int   tr_indx;    //!< indexes of verticies of triangles (O-based)

	int GetNumTr()    { return tr_indx.num_cols(); }     //!< Get number of triangles
	int GetNumVerts() const { return verts.num_cols(); } //!< Get the number of vertices in the surface

	virtual int  SetNumVerts( const int new_num_verts);  //!< Set the number of verticies
	
// tesserae surface representation 
// surface presented as patches (parts of main triangles) of the trianglulated spheres
// with a given surface area and spheres they belong to (that determine normal vectors)

	HaVec_ptr sph;  //!< Array of pointers to centers of spheres used to identify tesserae
	HaVec_double   srad; //!< Radii of spheres 
	HaVec_short  flag_sph; //!< indicator to show if spheres are allocated externally or other prop
    Vec3DValArray   tess; //!< array of tesserae
	HaVec_double tsurf;    //!< surface area of tesserae
	HaVec_short  itsph;   //!< index of a sphere tesserae belongs to
	HaVec_short  itsph_r; //!< index of a real sphere tesserae assigned to
    
	bool Print_info(std::ostream& sout, const int level) const; //!< dump information about the surface
	HaVec_double surface_alpha;   //!< Array of atomic surfaces
	HaVec_double volume_alpha;    //!< Array of atomic volumes
	HaMat_double d_volume_alpha;  //!< Matrix of atomic volumes derivatives 
	HaMat_double d_surface_alpha; //!< Matrix of atomic surfaces derivatives
	double surface_alpha_total;   //!< Molecular surface area calculated using Alpha Shape Theory
    int CalcMolSurfAlpha(int calc_d,double solv_rad, HaMat_double& cnt_crd_alpha, HaVec_double& cnt_rad_alpha); //!< Calculate Molecular surface using Alpha Shape theory, calc_d (=0 or =1) - to calculate derivatives of surface area vs atom coordinates

protected:
	double surf_volume;

};

const int MC_SIZE = 100000;
const int MV_SIZE = 100000;

extern "C"{

void gsurf_(int* ksurf,  //!< type of the surface ( =0 - VdW surf, =1 - Solvent Access surf, =2 Solvent Excluded Surf)
	   freal* rmin,     //!< for solv excl surf - min axx sphere radii (smaller the better) def = 0.5
	   freal* ofac,     //!< for solv excl surf - overlapping factor (from 0 to 1) larger the better def = 0.8
	   freal* rd,       //!< Solvent radius
       int* ndiv,   //!< level of subdivision of the original triangles (from 1 to 5) def = 3
	   logical* ass1,   //!< for solv excl surf - assign surfaces of axx spheres to real ones (def = FALSE)
	   int* natom,  //!< total number of centers 
	   logical* ghost,  //!< flag to indicate the presence of ghost atoms
       int*   np,    //!< the number of triangle forming the surface (OUTPUT)
       double* surf_volume);      //!< the surface volume

typedef struct
{
  freal xe[MC_SIZE];   //!< X coordinates of the center of the sphere
  freal ye[MC_SIZE];   //!< Y coordinates of the center of the sphere
  freal ze[MC_SIZE];   //!< Z coordinates of the center of the sphere
  freal re[MC_SIZE];   //!< radii of the spheres
  int_2 iuse[MC_SIZE];  //!< if =1 - not use, if =6 - use, if =3 - ghost sphere??
} csfe_type; 

typedef struct
{
	int_4  jvt1[3*60];
	int_4  jvt2[3*4];
} penta_type;

typedef struct      // surface triangle tessarae description
{
   int_4 ito[MV_SIZE]; //!< triangle index on a standard sphere
   int_4 iso[MV_SIZE]; //!< indexes of spheres tringles belong to(in CSFE) (1-based) i
   int_4 isa[MV_SIZE]; //!< indexes of atoms (initials spheres) triangles are assigned to
   freal xp[MV_SIZE];  //!< X coordinates of the centers of tesserae (Averaged position of secondary triangles) 
   freal yp[MV_SIZE];  //!< Y coordinates of the centers of tesserae (Averaged position of secondary triangles)
   freal zp[MV_SIZE];  //!< Z coordinates of the centers of tesserae (Averaged position of secondary triangles)
   freal ap[MV_SIZE];  //!< surface area of tesserae
} pun_type;

typedef struct
{
   double cv[32*3];
   freal xc1[15360];
   freal yc1[15360];
   freal zc1[15360];
} poli_type;

extern csfe_type csfe_;    //!< Coordinates and radii of atoms to build the molecular surface
extern pun_type pun_;
extern penta_type penta_;
extern poli_type poli_;

}


class HaDisplayedSurface: public HaSurface, public Object3D
{
public:
	HaDisplayedSurface();
	virtual ~HaDisplayedSurface();
	
	virtual void clear();

// Object3D virtuals:
	
	virtual int RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt);

	virtual int Translate( const Vec3D& tr_vec ); //!< Translate the object by tr_vec

	virtual int SetTransparency(double transp_new); //!< Set transparency of the surface
	bool ColourUniform(int r, int g, int b );
	

	virtual int SetNumVerts(const int new_num_verts);  //!< Set the number of verticies
	HaVec_short  colors;      /* array [numverts] of color table indexes */
	

protected:

};


class HaDot
{
public:
    HaDot() {};
	HaDot(const double new_x, const double new_y,
		  const double new_z, const int col);

	virtual ~HaDot();
 //short col;
 int col;//mikola's modifying

	double GetX() const { return xpos; }
	double GetY() const { return ypos; }
	double GetZ() const { return zpos; }

	double GetX_Ang() const { return xpos; }
	double GetY_Ang() const { return ypos; }
	double GetZ_Ang() const { return zpos; }

	void SetX( const double new_x ) { xpos = new_x; }
	void SetY( const double new_y ) { ypos = new_y; }
	void SetZ( const double new_z ) { zpos = new_z; }

	void SetX_Ang( const double new_x ) { xpos = new_x; }
	void SetY_Ang( const double new_y ) { ypos = new_y; }
	void SetZ_Ang( const double new_z ) { zpos = new_z; }

protected:

	double xpos;
	double ypos;
	double zpos;
};



class DotStruct : public Object3D
{
public:

    DotStruct();
	virtual ~DotStruct();

	std::vector<HaDot> dots;
	int GetCount() { return dots.size(); }

	void AddDot(double x, double y, double z, int col );

// Object3D virtual functions:

	virtual const char* GetObjName() const { return "Dot Surface"; }

	virtual int RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt );

	virtual int Translate( const Vec3D& tr_vec); // Translate the object by tr_vec

};



#if defined(HASURFACE_CPP)
const int  MAX_ATOM_GEOBALL = 20000;
extern "C" {
 int GeoBallSurface(int* switchval, int* natom, double* solv_rad, double* cnt_crd, double* cnt_rad, double* surface , double* volume, double* deriv_surf , double* deriv_vol);

// typedef struct
//	{
//		double coord[MAX_ATOM_GEOBALL];
//		double radius[MAX_ATOM_GEOBALL];
//	} geoball_data;
// extern geoball_data coord_;
}
#endif




#endif // End of !defined(HASURFACE_H)
