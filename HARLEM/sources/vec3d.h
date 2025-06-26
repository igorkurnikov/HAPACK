/*! \file vec3d.h 

    Classes to define a vector or point in 3D space in HARLEM

    \author Igor Kurnikov  
    \date 2000-2002
*/
#if !defined(VEC3D_H)
#define VEC3D_H

enum CoordUnits{ BOHR_U,   ANGSTROM_U };
enum AngleUnits{ DEGREE_U, RADIAN_U };

#include "hastl.h"
#include "haconst.h"

class HaMat_double;
class HaVec_double;
class PtrPtrMap;
class AtomGroup;


//! class to define 3D vector or a point in 3D space
class Vec3D 
{
protected:
  double pos[3];
public:
  Vec3D();
  Vec3D(double x,double y,double z);
  virtual ~Vec3D();

  inline double GetX() const { return pos[0]; } //!< X Coordinate in Ang
  inline double GetY() const { return pos[1]; } //!< Y Coordinate in Ang
  inline double GetZ() const { return pos[2]; } //!< Z Coordinate in Ang

  inline double GetX_Bohr() const { return pos[0]/BOHR_TO_ANG; } //!<  X Coordinate in Bohr
  inline double GetY_Bohr() const { return pos[1]/BOHR_TO_ANG; } //!<  Y Coordinate in Bohr
  inline double GetZ_Bohr() const { return pos[2]/BOHR_TO_ANG; } //!<  Z Coordinate in Bohr

  inline double GetX_Ang() const { return pos[0]; } //!<  X Coordinate in Ang
  inline double GetY_Ang() const { return pos[1]; } //!<  Y Coordinate in Ang
  inline double GetZ_Ang() const { return pos[2]; } //!<  Z Coordinate in Ang

  Vec3D& operator=(const Vec3D& ref_vec);
  Vec3D& operator=(double val) { pos[0] = val; pos[1] = val; pos[2] = val; return (*this);}
  void SetCoordFrom(const Vec3D& ref_vec);

  double& operator[](int i) { return pos[i];}              //!< 0-index based coordinate access 
  const double& operator[](int i) const { return pos[i];}  //!< 0-index based coordinate access, const method

  Vec3D& operator+=(const Vec3D& v); //!< add vector to a current vector

  friend Vec3D operator+(const Vec3D& v1,const Vec3D& v2);  //!< sum of two vectors
  friend Vec3D operator-(const Vec3D& v1,const Vec3D& v2);  //!< difference of two vectors
  friend Vec3D operator*(const double factor,const Vec3D& v);    //!< multiply the vector by a scaler
  
  double length() const; //!< Calculate vector's length (Ang)
  double length2() const; //!< Calculate vector's length squared (Ang*Ang)

  int IsClose(const Vec3D& pt, double crit=0.001) const; //< Check if the point pt is close to the current (all coordinates within crit from the corordinates of the current point)

  void SetX(const double x_new); //!< Set X Coordinate of the vector, argument is in Ang
  void SetY(const double y_new); //!< Set Y Coordinate of the point, argument is in  Ang
  void SetZ(const double z_new); //!< Set Z Coordinate of the point, argument is in  Ang

  void SetX_Bohr( double x_new ); //!<  Set X Coordinate using Bohr units
  void SetY_Bohr( double y_new ); //!<  Set Y Coordinate using Bohr units
  void SetZ_Bohr( double z_new ); //!<  Set Z Coordinate using Bohr units

  void SetX_Ang( double x_new ); //!< Set X Coordinate using Ang units
  void SetY_Ang( double y_new ); //!< Set Y Coordinate using Ang units
  void SetZ_Ang( double z_new ); //!< Set Z Coordinate using Ang units

  int SetFromStr( const char* str); //!< Set from string, assume 3 float number in the string

  void Scale(double scale ); //!< Multiply components of the vector by a constant  
  int normalize();           //!< Normalize to length = 1.0
  void SetZeros();           //!< Set zero values to the vector components 

  int RotatePt(Vec3D& still_pt, Vec3D& n, double cosa, double sina); //!< rotate point around axis with unit vec n passing through point still_pt 
  int Rotate(Vec3D& n, double cosa, double sina); //!< rotate vector around unit vector n by angle with cosa and sina
	
// Geometry Calculations:

static double CalcTorsion( const Vec3D* atm1, const Vec3D* atm2,
						   const Vec3D* atm3, const Vec3D* atm4);  //!< Calculate Torsion Angle (radians) this is standard def as used in JUMNA
static double CalcDihedral( const Vec3D* atm1, const Vec3D* atm2,
						    const Vec3D* atm3, const Vec3D* atm4); //!< Calculate Dihedral Angle (pi - torsional angle) (radians)
static double CalcAngle( const Vec3D* atm1, const Vec3D* atm2, const Vec3D* atm3); //!< Calculate Valence Angle between three atoms (radians)
static double CalcDistance( const Vec3D* atm1, const Vec3D* atm2,  //!< Calculate Distance between two atom
						    const CoordUnits units = ANGSTROM_U);
static double CalcDistanceSq( const Vec3D* atm1, const Vec3D* atm2,  //!< Calculate Squared Distance between two atom
						    const CoordUnits units = ANGSTROM_U);

static int VecProduct( Vec3D& vprod, const Vec3D& vec1, const Vec3D& vec2); //!< Calculate vector product of two vectors
static double DotProduct( const Vec3D& vec1, const Vec3D& vec2); //!< Calculate dot product between two vectors
static int diff(Vec3D& c, const Vec3D& vec1, const Vec3D& vec2); //!< Calculate difference between vec1 and vec2 
static int sum (Vec3D& c, const Vec3D& vec1, const Vec3D& vec2); //!< Calculate sum of vec1 and vec2 

static int CalcHlxParams(const Vec3D& c0_1, const Vec3D& v1_1, const Vec3D& v2_1, const Vec3D& v3_1,
                         const Vec3D& c0_2, const Vec3D& v1_2, const Vec3D& v2_2, const Vec3D& v3_2,
						 double& shift, double& slide, double& rise,
						 double& tilt, double& roll,double& twist,
						 int idir=1); //!< Compute Helical transformation parameters 
                             //!< shift,slide,rise,tip,roll,twist -  to characterize
                             //!< relative position of two bodies given by position of the centers and directions of axes
                             //!< idir - direction of the helix from 1 -> 2 or reverse

static int SetAtomPos( Vec3D* pptr, const Vec3D* ptr1, const Vec3D* pptr2, const Vec3D* pptr3,
					   double dist, double val_angle, double dih_angle); //!< Set position of the point using distance, valence angle and dihedral angle 
                                                                         //!< dih_angle - should be computed by CalcTorsion() not CalcDihedral()
  double* GetPnt(){return pos;}
};   


class PointIterator
//! Abstract class for an iterator on collections of points in 3D space (3D vectors)
{
public:

	PointIterator() {}
	virtual ~PointIterator() {}

//	virtual PointIterator* clone() const = 0;  //!< Create a copy of the iterator with the same state 
	virtual Vec3D* GetFirstPt() = 0; //!< Get First Point in the collection 
	virtual Vec3D* GetNextPt() = 0;  //!< Get Next Point in the collection

//	virtual Vec3D* GetPtByIdx(int idx) = 0; //!< Get Point by index ( zero based ) 
//	virtual    int GetNumPt() = 0;          //!< Get total of number of points in the collection
};

class PointIterator_const
//! Abstract class for an iterator on collections of points in 3D space (3D vectors)
{
public:

	PointIterator_const() {}
	virtual ~PointIterator_const() {}

//	virtual PointIterator_const* clone() const = 0;  //!< Create a copy of the iterator with the same state 
	virtual const Vec3D* GetFirstPt() = 0; //!< Get First Point in the collection 
	virtual const Vec3D* GetNextPt() = 0;  //!< Get Next Point in the collection

//	virtual Vec3D* GetPtByIdx(int idx) = 0; //!< Get Point by index ( zero based ) 
//	virtual    int GetNumPt() = 0;          //!< Get total of number of points in the collection
};


class PointContainer
//! Abstract class of a container (group) of points in 3D space
{
public:
    virtual PointIterator*       GetPointIteratorPtr() = 0; //!< create Point Interator and return pointer (should be deleted when done)
    virtual PointIterator_const* GetPointIteratorPtr() const = 0; //!< create const Point Interator and return pointer (should be deleted when done)
	virtual int GetNumPt() const = 0;                        //!< Get total of number of points in the collection
    
//  virtual std::string GetRef() const = 0;  //!< Get the text reference for a point collection
//	virtual bool FillRef(char* buf,int mode = 0) const = 0; //!< Write text reference for a point collection to buffer 

	int IsWithinRadius(Vec3D* pptr, double dlimit2 ) const; //!< Check if the point is within radius limit from the points of the collection 

	bool GetMinMaxCrd(double& MinX_v, double& MinY_v, double& MinZ_v,
				      double& MaxX_v, double& MaxY_v, double& MaxZ_v) const; //!< Get Min and Max coordinates of the points of the container
	 
	bool GetAverageCoord(double& avx, double& avy, double& avz) const; //!<  Get average coordinates of the points of the container
	Vec3D GetAverageCoord() const; //!< Get average coordinates of the points of the container

	static int FindCoordMatch( PointContainer& g1, PointContainer& g2, PtrPtrMap& pt_pt_map); 

	
	static double GetSuperimposeMat( HaVec_double& ref_crd, PointContainer& g2,
                            HaMat_double& rot_mat, HaVec_double& transl_vec); //!< find the rotation matrix and translation vector that  
                                                                                           //!< superimpose collections of points g2 over g1
	                                                                                       //!< return eps - root mean square deviation (Ang) ( negative if error )
	static int GetSuperimposeMat( PointContainer& g1 , PointContainer& g2,
                            HaMat_double& rot_mat, HaVec_double& transl_vec, double& eps); //!< find the rotation matrix and translation vector that  
                                                                                           //!< superimpose collections of points g2 over coordinates in ref_crd
	                                                                                       //!< eps - root mean square deviation (Ang)
	static double CalcRMSD( PointContainer& g1, PointContainer& g2, int transform = TRUE ); //!< Calculate RMSD( root mean square deviation between two Point collections) finding best superimpose of g2 over g1 

	int Transform(const HaMat_double& rot_mat, const HaVec_double& transl_vec); //!< Transform coordinates of points using translation vector and rotation matrix 
	int SaveCrdToArray(HaVec_double& crd_arr);  //!< Save coordinates of the Point collection to a double array 
	HaVec_double GetCrdArray();           //!< Get Coordinate Array
	int SetCrdFromArray( const HaVec_double& crd_arr); //!< Set coordinates of the Point collection from a double array
};

namespace harlem
{
	namespace geom
	{
		static double GetSuperimposeMat( HaVec_double& ref_crd, PointContainer& g2, HaMat_double& rot_mat, HaVec_double& transl_vec)
		{
			return PointContainer::GetSuperimposeMat( ref_crd, g2, rot_mat, transl_vec);
		}
	}
}

class PointIteratorGen
{
public: 
	PointIteratorGen( PointContainer& pt_col) { ppt_col = &pt_col; ppt_itr = pt_col.GetPointIteratorPtr(); }
	PointIteratorGen( PointIteratorGen& pitr) { ppt_col = pitr.ppt_col; ppt_itr = ppt_col->GetPointIteratorPtr(); }
	virtual ~PointIteratorGen() { delete ppt_itr; }

	Vec3D* GetFirstPt() { return ppt_itr->GetFirstPt(); }
	Vec3D* GetNextPt()  { return ppt_itr->GetNextPt();  }
	int    GetNumPt()      { return ppt_col->GetNumPt(); }

private:
	PointIterator*   ppt_itr;
	PointContainer* ppt_col;

};

class PointIteratorGen_const
{
public: 
	PointIteratorGen_const( const PointContainer& pt_col) { ppt_col = &pt_col; ppt_itr = pt_col.GetPointIteratorPtr(); }
	PointIteratorGen_const( PointIteratorGen_const& pitr) { ppt_col = pitr.ppt_col; ppt_itr = ppt_col->GetPointIteratorPtr(); }
	virtual ~PointIteratorGen_const() { delete ppt_itr; }

	const Vec3D* GetFirstPt() { return ppt_itr->GetFirstPt(); }
	const Vec3D* GetNextPt()  { return ppt_itr->GetNextPt();  }
	int   GetNumPt()          { return ppt_col->GetNumPt(); }

private:
    PointIterator_const*   ppt_itr;
	const PointContainer* ppt_col;

};


class Vec3DValArray;

class Vec3DValArrayIterator: public PointIterator
{
public:
    Vec3DValArrayIterator(Vec3DValArray* pt_array);
	virtual ~Vec3DValArrayIterator();

	virtual Vec3D* GetFirstPt();
	virtual Vec3D* GetNextPt();

protected:
	int cur_idx;
	Vec3DValArray* cur_arr;
};

class Vec3DValArrayIterator_const: public PointIterator_const
{
public:
    Vec3DValArrayIterator_const(const Vec3DValArray* pt_array);
	virtual ~Vec3DValArrayIterator_const();

	virtual const Vec3D* GetFirstPt();
	virtual const Vec3D* GetNextPt();

protected:
	int cur_idx;
	const Vec3DValArray* cur_arr;
};

class Vec3DValArray: public std::vector<Vec3D>, public PointContainer
{
public:
    Vec3DValArray();
    Vec3DValArray(int n);
	virtual ~Vec3DValArray();

	virtual PointIterator*       GetPointIteratorPtr();
	virtual PointIterator_const* GetPointIteratorPtr() const;
    virtual int GetNumPt() const;
};

class Quaternion
//! Based on code by Laurent Schmalen and KoLoR
//! E-Mail:  Laurent.Schmalen@Ci.Educ.Lu
//!
{
  private:
    double W, X, Y, Z;      // components of a quaternion

  public:
    // default constructor
    Quaternion();
    Quaternion(const double, const double = 0.0, const double = 0.0, const double = 0.0);
	// default destructor
	virtual ~Quaternion();

    // quaternion multiplication
    friend Quaternion operator * (const Quaternion&, const Quaternion&);
    const Quaternion& operator *= (const Quaternion&);

    // conjugates Quaternion
    const Quaternion& operator ~ ();

    // this function inverts the quaternion
    const Quaternion& operator - ();

    // this normalizes a quaternion 
    const Quaternion& Normalize();
	int GetQuaternion(double& w, double& x, double& y, double& z);
	int SetQuaternion(double& w, double& x, double& y, double& z);

    const Quaternion& QuaternionFromAxis(const double, double, double, double);
    static int QuaternionToRotMat(const Quaternion& q, HaMat_double& rmat);
	static int RotMatToQuaternion(const HaMat_double& rmat, Quaternion& q);

    void Slerp(const Quaternion&, const Quaternion&, const double);

    // some additional Quaternion functions
    // getting the exponent of a quaternion and logarithm
    const Quaternion& exp();
    const Quaternion& log();
	void PrintOn();
};


class Rot3D
//! Class to describe rotations in 3D space
{
public:
// Assumptions for Euler angles:
// -PI <  phi, psi < PI,   -1 < cos_theta < 1 

	static int RotMatToEuler(const HaMat_double& rmat, double& phi, double& cos_theta, double& psi); //!< Compute rotation matrix for given Euler angles
	
	static int EulerToRotMat(const double phi, const double cos_theta, const double psi, 
		                        HaMat_double& rmat); //!< Compute rotation matrix corresponding to Euler angles

    static int IncrEulerAng( double& phi, double& cos_theta, double& psi, 
					         double delt_phi, double delt_cos_theta, double delt_psi); //!< increment Euler anlges taking into account limits and continuity

    static int NormalizeEulerAng( double& phi, double& cos_theta, double& psi); //!< assign Euler anlges values within
	static int RotMatToQuat(const HaMat_double& rmat, Quaternion& q);
	static int QuatToRotMat(const Quaternion& q, HaMat_double& rmat);
};



class BoxPartition;

class BoxRegionPointIterator: public PointIterator
{
public:

	BoxRegionPointIterator();
	BoxRegionPointIterator(const BoxRegionPointIterator& iter);
	virtual ~BoxRegionPointIterator();

// Virtual function of the PointIterator 

	virtual Vec3D* GetFirstPt();
	virtual Vec3D* GetNextPt();
	virtual int GetNumPt(); 

	int iax_min;
	int iay_min;
	int iaz_min;
	int iax_max;
	int iay_max;
	int iaz_max;

	BoxPartition* partition; 
};


class BoxPartition : public std::vector< VecPtr >
//!< Class for the BoxPartition of a collection of points into cubic sections to speed up some algorithm requiring cycles on nearby atoms
{
public:
	BoxPartition();

	int SetDimensions(int nx_new, int ny_new, int nz_new); //!< Set Dimensions of the partition of box
	
	int SetBoundaries(double xmin_new, double ymin_new, double zmin_new, 
		              double xmax_new, double ymax_new, double zmax_new); //!< Set geometrical boundaries of the box


	int DistributePointsToCells( PointContainer& pt_coll); //!< partition points of a point collection to cells of the box partition
	int AddPoint( Vec3D* pt); //!< Add Point to the partition

	double xmin;   //!< X coordinate of front-left-lower corner of the box
	double ymin;   //!< Y coordinate of front-left-lower corner of the box
	double zmin;   //!< Z coordinate of front-left-lower corner of the box
 
	double dx;     //!< Size of the interval of the partion in X direction 
	double dy;     //!< Size of the interval of the partion in Y direction 
	double dz;     //!< Size of the interval of the partion in Z direction 

	int nx;       //!< number of intervals in X direction 
	int ny;       //!< number of intervals in Y direction 
	int nz;       //!< number of intervals in Z direction 

	int GetIdx_x( double x_coord ); //!< Get index of the interval in X direction corresonding to X-coord (-1 if out of interval) 
	int GetIdx_y( double y_coord ); //!< Get index of the interval in Y direction corresonding to Y-coord (-1 if out of interval) 
	int GetIdx_z( double z_coord ); //!< Get index of the interval in Z direction corresonding to Z-coord (-1 if out of interval) 

	int GetPointCellIdx(const Vec3D* ppt);             //!< Get the index of the cell containing the point (-1) if out of the box
	int GetPointCellIdx(double x, double y, double z); //!< Get the index of the cell containing the point (-1) if out of the box
    int GetLinCellIdx( int ix, int iy, int iz); //!< Get Linear Cell index from coordinates indexes 

	int ixrad;  //!< extention indexes in X-direction ( +- ixrad) of the neibouring region of the cells    
	int iyrad;  //!< extention indexes in Y-direction ( +- iyrad) of the neibouring region of the cells 
	int izrad;  //!< extention indexes in Z-direction ( +- izrad) of the neibouring region of the cells 

	int SetRegionRad(double dist); //!< Set radius of the local region (in Bohr) to search for neighbors of the point(or cell)   
	
	int GetNeighbors(Vec3D& pt, AtomGroup& neighbors);

	BoxRegionPointIterator GetRegionPointIterator();


};


#endif /* !VEC3D_H */
