/*! \file hatypes.h

    Some basic types to use in HARLEM 
 
    \author Igor Kurnikov  
    \date 2010-
*/

#if !defined(HATYPES_H)
#define HATYPES_H

#include <mpi.h>

#include "hastl.h"
#include "hastring.h"

#include <boost/any.hpp>

namespace harlem
{
	typedef long long HashID;

	std::string unhash(const HashID& id);
    HashID hash(const std::string& name);

	class HashMap
	{
	public:
		HashMap();
		HashMap( const HashMap& ref );
		virtual ~HashMap();

		virtual void Copy( const HashMap& ref );
		virtual HashMap* clone() const;

		static HashMap empty_map;

		// Set functions
		void set_i( const HashID& parmid, const int& value ); 
		void set_d( const HashID& parmid, const double& value );
		void set_s( const HashID& parmid, const std::string& value );
		void set_a( const HashID& parmid, const boost::any& value );

		void set_i( const std::string& parmname, const int& value );
		void set_d( const std::string& parmname, const double& value ); 
		void set_s( const std::string& parmname, const std::string& value ); 
		void set_a( const std::string& parmname, const boost::any& value ); 

		// Get functions

		int get_i(const HashID& parmid) const;
		double get_d(const HashID& parmid) const;
		std::string get_s(const HashID& parmid) const;
		boost::any get_a(const HashID& parmid) const;

		int get_i(const std::string& parmname) const;
		double get_d(const std::string& parmname) const;
		std::string get_s(const std::string& parmname) const;
		boost::any get_a(const std::string& parmname) const;

		// Test functions

		bool has_i(const HashID& parmid) const;
		bool has_d(const HashID& parmid) const;
		bool has_s(const HashID& parmid) const;
		bool has_a(const HashID& parmid) const;


		bool has_i(const std::string& parmname) const; 
		bool has_d(const std::string& parmname) const; 
		bool has_s(const std::string& parmname) const; 
		bool has_a(const std::string& parmname) const;

	protected:

		std::map< HashID, int >         i_params;
		std::map< HashID, double >      d_params;
		std::map< HashID, std::string > s_params;
		std::map< HashID, boost::any >  a_params;
	};

	class SaveOptions : public harlem::HashMap
	{
	public:

		SaveOptions() { SetSaveHeader(false); SetSaveFooter(false); }
		SaveOptions( const SaveOptions& ref ) { Copy(ref); }
		virtual ~SaveOptions() {}

		virtual void Copy( const harlem::HashMap& ref ) 
		{  
			harlem::HashMap::Copy(ref);
			const SaveOptions* pref = dynamic_cast<const SaveOptions*>(&ref);
			if( pref != NULL )
			{
				SetSaveHeader( pref->ToSaveHeader() ); 
				SetSaveFooter( pref->ToSaveFooter() ); 
			}
		}

		virtual harlem::HashMap* clone() const { return new SaveOptions(*this); }  
	
		bool ToSaveHeader() const { return save_header; }
		void SetSaveHeader( bool set_par = true ) { save_header = set_par; }

		bool ToSaveFooter() const { return save_footer; }
		void SetSaveFooter( bool set_par = true ) { save_footer = set_par; }

	protected:
		bool save_header;
		bool save_footer;
	};
}

class HaEnum
//!< Enum-like Class in HARLEM to be used as enumerated parameters
{
public:
	virtual int& value() = 0;
	virtual const char* label() const = 0;
	virtual int SetWithValue(int value) = 0; //!< Set using integer value corresponding to one of the enum values
	virtual int SetWithLabel(const char* label) = 0; //!< Set using label string

	HaEnum& operator=(int value) { SetWithValue(value); return (*this); }  

	int Bcast(MPI_Comm& comm, int root = 0); //!< Broadcast HaEnum over MPI communicator

	virtual std::vector<int> all_values() const = 0;
	virtual std::vector<std::string> GetAllLabels() const = 0;
	virtual std::vector<std::string> GetActiveLabels() { return GetAllLabels(); }
};

class HaEnum1: public HaEnum
//! Enum-like class with int to string map
{	
public:
	virtual std::vector<int> all_values() const;
	virtual std::vector<std::string> GetAllLabels() const;
protected:
	virtual IntStrMap& GetLabelsMap() const = 0 ;
};


#endif // end !defined(HATYPES_H)