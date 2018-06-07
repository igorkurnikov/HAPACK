/*! \file rigidbodycoord.h

    Classes to define Rigid Body coordinates manipulations in HARLEM

    \author Igor Kurnikov 
 
    \date 2009

*/
#ifndef RIGIDBODYCOORD_H
#define RIGIDBODYCOORD_H

#include "atom_mapping.h"

namespace harlem
{
	//! \brief Class to define Rigid Body Coordinates for one or more 3D Objects
	//! \nosubgrouping
	class RigidBodyCoord: public harlem::Coord
	{
	public:
		RigidBodyCoord();
		RigidBodyCoord( const RigidBodyCoord& ref );
		virtual ~RigidBodyCoord();

		virtual std::string GetClassName() const { return "RigidBodyCoord"; } 

		virtual harlem::Coord* clone(); 
		virtual int SetFrom(const harlem::Coord* pcrd);
		virtual HaVec_double AsVecDouble() const;                  //!< Transform Coordinates to a Vector of double values 
		virtual int SetFromVecDouble(const HaVec_double& dbl_vec); //!< Set Coordinates from a vector of double values

		virtual int GetNumCrd() const;        //!< Get Total Number of Internal Coordinates
		virtual void FreezeCrd(int idx);         //!< Freeze Internal Coordinate with the index idx
		virtual int IsCrdFrozen(int idx) const; //!< Check if the coordinate is frozen
		virtual int LoadFromStream(std::istream& is, const harlem::HashMap* popt = NULL ); //!< Read Coordinates from stream
		virtual int SaveToStream(std::ostream&  os, const harlem::HashMap* popt = NULL ) const; //!< Write Coordinates to stream
		int SetFromCurrAtomCrd(AtomContainer* at_cont, int iobj);     //!< Set internal coordinates of iobj object from Current Cartesian Coordinates of Atom Collection 
		int SetFromCurrAtomCrd(vector<AtomContainer*> vec_at_cont);  //!< Set from Current Cartesian Coordinates of an Array of Atom Collections

		virtual void SetNumObj(int n_obj_new);    //!< Set Number of objects which positions are described by these coordinates
		int GetNumObj() const;            //!< Get the number of rigid body objects described by coordinates

		double GetPhi(int iobj) const;      //!< Get Euler Phi angle  of iobj-th object  (0-based index)
		double GetCosTheta(int iobj) const; //!< Get Euler Cos(theta) of iobj-th object  (0-based index)
		double GetPsi(int iobj) const;      //!< Get Euler Psi angle  of iobj-th object  (0-based index)
		double GetTransX(int iobj) const;   //!< Get Translation Along X axis of iobj-th object  (0-based index)
		double GetTransY(int iobj) const;   //!< Get Translation Along Y axis of iobj-th object  (0-based index)
		double GetTransZ(int iobj) const;   //!< Get Translation Along Z axis of iobj-th object  (0-based index)

		void SetPhi(int iobj, double phi_new);    //!< Set Euler Phi angle  of iobj-th object  (0-based index)
		void SetCosTheta(int iobj, double cos_theta_new); //!< Set Euler Cos(theta) of iobj-th object  (0-based index)
		void SetPsi(int iobj, double psi_new);    //!< Set Euler Psi angle  of iobj-th object  (0-based index)
		void SetTransX(int iobj, double x_new);   //!< Set Translation Along X axis of iobj-th object  (0-based index)
		void SetTransY(int iobj, double y_new);   //!< Set Translation Along Y axis of iobj-th object  (0-based index)
		void SetTransZ(int iobj, double z_new);   //!< Set Translation Along Z axis of iobj-th object  (0-based index)

		void FreezeObject(int iobj);         //!< Freeze motion of object iobj 
		int  IsObjectFrozen(int iobj) const; //!< Check if object iobj is frozen 

	protected:
		int n_obj;
		HaVec_double crd_v;  
		HaVec_int    frozen_idx;
	};

	class RigidBodyCoordDiscretized: public RigidBodyCoord
	{
	public: 
		RigidBodyCoordDiscretized();
		RigidBodyCoordDiscretized( const RigidBodyCoordDiscretized& ref );
		virtual ~RigidBodyCoordDiscretized();

		virtual harlem::Coord* clone(); 
		virtual int SetFrom(const harlem::Coord* pcrd); 
		virtual std::string GetClassName() const { return "RigidBodyCoordDiscretized"; } 

		virtual void SetNumObj(int n_obj_new);    //!< Set Number of objects which positions are described by these coordinates

		virtual void ConvertDiscrCrdToFloat();  //!< Compute values of indexes corresponding to float values of internal coordinates
		virtual void ConvertFloatCrdToDiscr();  //!< Compute float values of coordinate corresponding to current integer index values 

		void SetStandardLimits(); //!< Set Standard Limits on Internal coordinates 

		int SetDiscrNumForCrd(int idx, int npt);            //!< Set number of discretization points(npt) for internal coordinate with index idx 
		int SetLimits(int idx, double amin, double amax);   //!< Set limits for internal coordinate with index idx 

		HaVec_int crd_v_int; //!< Integers that determine coordinates of the system (from 0 - to dim_crd[i] )  
		HaVec_int dim_crd;   //!< Number of discretization points for coordinates ( number of intervals - 1 )  

	protected:
		HaVec_double  crd_min;  //!< Minimal values for the coordinates
		HaVec_double  crd_max;  //!< Maximal values for the coordinates
	};

};

#endif // !defined(RIGIDBODYCOORD_H) 


