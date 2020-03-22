/*!  \file object3d.h

    Classes to define general 3D objects that may be displayed in HARLEM
 
    \author Igor Kurnikov , University of Pittsburgh 
    \date 1999-2002

*/
#if !defined(OBJECT3D_H)
#define OBJECT3D_H

#include "hastring.h"

const int OBJ3D_MOLECULE = 1;     //!< Type of the object - molecule (HaMolecule class)
const int OBJ3D_SURFACE  = 2;     //!< Type of the object - surface  (HaSurface class)
const int OBJ3D_DOT_SURFACE = 3;  //!< Type of the object - dot surface (DotStruct class)
const int OBJ3D_BOX = 4;//<mikola 30July06
const int OBJ3D_MEMBRANEZ = 5;//<mikola 30July06
const int OBJ3D_MEMBRANE_TUBE = 6;//<mikola 30July06
const int OBJ3D_MATRIX = 7;//<mikola 30July06
const int OBJ3D_PlaneViewerOfHaField3D = 8;//<mikola 26March08
class Vec3D;
class HaMat_double;
class HaMolView;//<mikola 30July06
class Canvas3D;//<mikola 30July06
class HaField3D;//<mikola 26March08

class Object3D 
//! Class to define 3D objects that can be displayed and manipulated
{
public:

	Object3D(const int new_obj_type, const char* new_name= "GEN_Object3D" );
	Object3D(const Object3D& obj_ref);
	virtual ~Object3D();

	bool IsConnected() { return connect_flag; }
	void SetConnected(const bool status) { connect_flag = status; } 

	bool IsDisplayed() { return (displayed_flag) ; }
	void SetDisplayed(const bool status) { displayed_flag = status; }

	virtual const char* GetObjName() const { return name.c_str(); } 
	virtual bool SetObjName(const char* new_name) {  name = new_name; return true;}

	virtual int RotateX( double theta, const Vec3D& cnt );  //!< Rotate around X axis that go through cnt
	virtual int RotateY( double theta, const Vec3D& cnt );  //!< Rotate around Y axis that go through cnt
	virtual int RotateZ( double theta, const Vec3D& cnt );  //!< Rotate around Z axis that go through cnt

	virtual int RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt); //!< Rotate object around the center  

	virtual int Translate( const Vec3D& tr_vec); //!< Translate the object by (dx,dy,dz)

	int GetObjType() const { return obj_type; }
 
    virtual int DrawObj(HaMolView* molview){return TRUE;}//<mikola 30July06

	virtual int SetTransparency(double trasp_new); //!< SetObjectTransparency
	double transparency;      //!< transparency (from 0 to 1) works only for surfaces

	//! Get Min and Max coordinates of this object
	virtual bool GetObjectMinMaxCrd(double& MinX_v, double& MinY_v, double& MinZ_v,
		double& MaxX_v, double& MaxY_v, double& MaxZ_v) const{return false;}
protected:

	int obj_type; //!< Type of the object 
	bool connect_flag;
	bool displayed_flag; 	
	
	std::string name;
};

const int SCREEN_ORIENTATION = 0;
const int REFERENCE_ORIENTATION   = 1;
//>mikola 30July06
class BoxObj3D : public Object3D
{
  public:

    BoxObj3D(const char* new_name= "BoxObj3D",int r=128,int g=128,int b=255,float the_x0=-1.0,float the_y0=-1.0,float the_z0=-1.0,float the_x1=1.0,float the_y1=1.0,float the_z1=1.0);
    virtual ~BoxObj3D();

    float x0,y0,z0,x1,y1,z1;
    void SetBox(float the_x0,float the_y0,float the_z0,float the_x1,float the_y1,float the_z1);
    
    int Color;
// Object3D virtual functions:
    virtual int RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt );

    virtual int Translate( const Vec3D& tr_vec); // Translate the object by tr_vec

    virtual int DrawObj(HaMolView* molview);

	virtual bool GetObjectMinMaxCrd(double& MinX_v, double& MinY_v, double& MinZ_v,
		double& MaxX_v, double& MaxY_v, double& MaxZ_v) const;
};
class MatrixObj3D : public Object3D
{
  public:

    MatrixObj3D(const char* new_name= "MatrixObj3D",int r=255,int g=255,int b=255,float the_x0=-1.0,float the_y0=-1.0,float the_z0=-1.0,float the_dx1=1.0,float the_dy1=1.0,float the_dz1=0.0,float the_dx2=1.0,float the_dy2=1.0,float the_dz2=0.0);
    virtual ~MatrixObj3D();

    float x0,y0,z0,dx1,dy1,dz1,dx2,dy2,dz2;
    int Ni,Nj;
    void SetMatrixGuiders(float the_x0,float the_y0,float the_z0,float the_dx1,float the_dy1,float the_dz1,float the_dx2,float the_dy2,float the_dz2);
    void SetColors(int Nx,int Ny,float * fmap,double elpot_high_val,double elpot_low_val);
    int Color;
    HaMat_int ColMat;
// Object3D virtual functions:
    virtual int RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt );

    virtual int Translate( const Vec3D& tr_vec); // Translate the object by tr_vec

    virtual int DrawObj(HaMolView* molview);
};
class MembraneZObj3D : public Object3D
{
  public:

    MembraneZObj3D(const char* new_name= "MembraneZObj3D",int r=128,int g=255,int b=128,float the_x0=-1.0,float the_y0=-1.0,float the_z0=-1.0,float the_x1=1.0,float the_y1=1.0,float the_z1=1.0,float the_mem_z1=-1.0,float the_mem_z2=1.0);
    virtual ~MembraneZObj3D();

    float x0,y0,z0,x1,y1,z1;
    float mem_z1,mem_z2;
    float hstep;
    void SetBox(float the_x0,float the_y0,float the_z0,float the_x1,float the_y1,float the_z1);
    void SetMembraneZ(float the_mem_z1,float the_mem_z2);
    
    int Color;
// Object3D virtual functions:
    virtual int RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt );

    virtual int Translate( const Vec3D& tr_vec); // Translate the object by tr_vec

    virtual int DrawObj(HaMolView* molview);
};
class TubeObj3D : public Object3D
{
  public:

    TubeObj3D(const char* new_name= "TubeObj3D",int r=128,int g=255,int b=128,float the_x0=-1.0,float the_y0=-1.0,float the_z0=-1.0,float the_x1=1.0,float the_y1=1.0,float the_z1=1.0);
    virtual ~TubeObj3D();

    float x0,y0,z0,x1,y1,z1;
    
    float mem_z1,mem_z2,mem_x,mem_y;
    float hstep,phistep;
    float R1,R2;
    void SetBox(float the_x0,float the_y0,float the_z0,float the_x1,float the_y1,float the_z1);
    void SetTube3d(float the_mem_z1,float the_mem_z2,float the_mem_x,float the_mem_y,float the_R1,float the_R2);
    
    int Color;
    int Style;
// Object3D virtual functions:
    virtual int RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt );

    virtual int Translate( const Vec3D& tr_vec); // Translate the object by tr_vec

    virtual int DrawObj(HaMolView* molview);
	virtual bool GetObjectMinMaxCrd(double& MinX_v, double& MinY_v, double& MinZ_v,
		double& MaxX_v, double& MaxY_v, double& MaxZ_v) const;
};
//<mikola 30July06
//<mikola 26March08
class PlaneViewOfHaField3D : public Object3D
{
	public:

		PlaneViewOfHaField3D(HaField3D *_field,const char* new_name= "PlaneViewOfHaField3D", int OwnerOfData=0);
		virtual ~PlaneViewOfHaField3D();
		enum {PlaneXY=0,PlaneYZ=1,PlaneZX=2};
		int OwnerOfData;
		
		
		
		void SetMinMax(double m_Min,double m_Max);
		void SetLevel(int NewLevel);
		void SetPlaneXY();
		void SetPlaneYZ();
		void SetPlaneZX();
		void SetPlane(int newPlane);
		void SetHideZeroValues(bool newHideZeroValues);
		int GetPlane(){return Plane;}
		int GetLevel(){return Level;}
		double GetValueMin(){return ValueMin;}
		double GetValueMax(){return ValueMax;}
		HaField3D * GetHaField3D(){return field;}
		bool GetHideZeroValues(){return HideZeroValues;}
		
		void SetColorsOfPlane();
// Object3D virtual functions:
		virtual int RotateObj( const HaMat_double& rot_mat, const Vec3D& cnt );
		//! Translate the object by tr_vec
		virtual int Translate( const Vec3D& tr_vec); 
		virtual int DrawObj(HaMolView* molview);
		virtual bool GetObjectMinMaxCrd(double& MinX_v, double& MinY_v, double& MinZ_v,
			double& MaxX_v, double& MaxY_v, double& MaxZ_v) const;
	private:
		bool HideZeroValues;
		int Plane;
		int Level;
		double ValueMin,ValueMax;
		HaField3D *field;
		
		HaMat_int ColMat;
		float x0,y0,z0,dx1,dy1,dz1,dx2,dy2,dz2;
		int Ni,Nj;
};
//>mikola 26March08

#endif // End of !defined(OBJECT3D_H)
