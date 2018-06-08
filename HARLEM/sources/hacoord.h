/*! \file hacoord.h

    Abstract Classes to define general coordinates in HARLEM

    \author Igor Kurnikov 
 
    \date 2009

*/

#ifndef HACOORD_H
#define HACOORD_H

#include "hatypes.h"

class HaVec_double;

namespace harlem
{
	//! \brief Abstract Class to define Coordinates of the system
	//! \nosubgrouping
	class Coord
	{
	public:
		virtual std::string GetClassName() const = 0;   //!< Get Class Name of the Coordinate

		virtual Coord* clone() = 0;                   //!< Get a copy of coordinates
		virtual HaVec_double AsVecDouble() const = 0;   //!< Transform Coordinates to a Vector of double values 
		virtual int SetFrom(const Coord* pcrd) = 0;   //!< Set Coordinates from the other coordinate object   
		virtual int SetFromVecDouble(const HaVec_double& dbl_vec) = 0;   //!< Set Coordinates from a vector of double values
		virtual int LoadFromStream(std::istream& is, const harlem::HashMap* popt = NULL ) = 0; //!< Read Coordinates from stream
		virtual int SaveToStream(std::ostream&  os,  const harlem::HashMap* popt = NULL ) const = 0;  //!< Write Coordinates to stream 
	};

	class LoadCrdOptions: public harlem::HashMap
	{
	public:
		LoadCrdOptions() { SetStdOptions(); }
		LoadCrdOptions( const LoadCrdOptions& ref ) { SetStdOptions(); Copy(ref); }
		virtual ~LoadCrdOptions() {}

		void SetStdOptions() { SetLoadNotFrozenCrd(true); }

		virtual void Copy( const harlem::HashMap& ref ) 
		{  
			const LoadCrdOptions* pref = dynamic_cast<const LoadCrdOptions*>(&ref);
			if( pref != NULL )
			{
				SetLoadNotFrozenCrd( pref->ToLoadNotFrozenCrd() ); 
			}
		}

		virtual harlem::HashMap* clone() const { return new LoadCrdOptions(*this); }  

		void SetLoadNotFrozenCrd( bool set_par ) { load_not_frozen = set_par; } 
		bool ToLoadNotFrozenCrd() const { return load_not_frozen; }

		void SetLoadAllCrd( bool set_par ) { load_not_frozen = !set_par; }
		bool ToLoadAllCrd() const { return !load_not_frozen; }  

	protected:
		bool load_not_frozen;
	};

	class SaveCrdOptions: public harlem::HashMap
	{
	public:
		SaveCrdOptions() { save_not_frozen = true; }
		SaveCrdOptions( const SaveCrdOptions& ref ) { SetStdOptions(); Copy(ref); }
		virtual ~SaveCrdOptions() {}

		void SetStdOptions() { SetSaveNotFrozenCrd(true); }

		virtual void Copy( const harlem::HashMap& ref ) 
		{  
			const SaveCrdOptions* pref = dynamic_cast<const SaveCrdOptions*>(&ref);
			if( pref != NULL )
			{
				SetSaveNotFrozenCrd( pref->ToSaveNotFrozenCrd() ); 
			}
		}

		virtual harlem::HashMap* clone() const { return new SaveCrdOptions(*this); }  

		void SetSaveNotFrozenCrd( bool set_par ) { save_not_frozen = set_par; }
		bool ToSaveNotFrozenCrd() const { return save_not_frozen; }

		void SetSaveAllCrd( bool set_par ) { save_not_frozen = !set_par; }
		bool ToSaveAllCrd() const { return !save_not_frozen; }  

	protected:
		bool save_not_frozen;
	};

};

#endif // !defined(HACOORD_H) 


