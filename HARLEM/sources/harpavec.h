/*! \file harpavec.h

    Classes to define RPA vector

    \author Igor Kurnikov  
    \date 1998-2002
*/
#ifndef HARPAVEC_H
#define HARPAVEC_H

enum OPER_TYPE {REAL_OPER,IMAG_OPER};

#include "halocorb.h"
#include "halocexcit.h"
#include "haatgroup.h"
#include "hamatdb.h"


class HaRPAvec
//!< Class to define RPA vector 
{
public:
	HaRPAvec();
	HaRPAvec(HaQCMod & qc_mod);
	HaRPAvec(const HaLocExcit & excit,HaQCMod & qc_mod);
	virtual ~HaRPAvec();

	int GetNBfunc() const;
	
	bool Print_info(ostream &sout, const int level) const;

    int GetNumOccMO() const; // the number occupied MOs
	int GetNumVacMO() const; // the number of vacant MOs
	                         // used to define the excitation array 

    int SetFromLocExcit(const HaLocExcit & excit);
    // set RPA vector expanding local excitation

	int SetFromAOMat(const HaMat_double & aomat,
		             const OPER_TYPE optyp=REAL_OPER );
	// set RPA vector from AO density assumed premultiplied on bith sides by S^-1

	int SetFromLOGrpMat(const std::string& gid1,const std::string& gid2, 
		                    const HaMat_double & fmloc,
							const OPER_TYPE optyp=REAL_OPER);
	//!< set RPA vector from the submatrix of operator matrix fmat 
	//!< (premultiplied by S^-1) on active local orbitals of groups ig1 and ig2

	friend double
		   SProd(const HaRPAvec & left, const HaRPAvec & right); 
	// Calculate scalar product between two RPA vector
	
	friend HaMat_double
		   SProd(const vector<HaRPAvec> & left, const vector<HaRPAvec> & right); 

	friend HaVec_double
		   SProd(const HaRPAvec & RPAv, const vector<HaRPAvec> & RPAv_arr); 

	// Calculate scalar product between arrays of RPA vectors 
	
	friend double dot2(const HaRPAvec & left, const HaRPAvec & right);

	friend HaVec_double
		   dot2(const vector<HaRPAvec> & left, const vector<HaRPAvec> & right);

	friend double norm2(const HaRPAvec & RPAv);
	friend HaVec_double norm2(const vector<HaRPAvec> & RPAv_arr);

	friend class HaRPAHam;
	friend class HaRPAResolv;

	// multiply by a const:
	friend HaRPAvec 
	operator*(const double factor, const HaRPAvec & RPAv); 

	friend vector<HaRPAvec> 
	operator*(const HaVec_double & vfactor, const vector<HaRPAvec> & RPAv_arr); 


	friend HaRPAvec
	operator+(HaRPAvec & left, HaRPAvec & right); 

	friend vector<HaRPAvec>
	operator+(vector<HaRPAvec> & left, vector<HaRPAvec> & right); 

	friend vector<HaRPAvec> 
	operator-(vector<HaRPAvec> & left, vector<HaRPAvec> & right); 

	HaRPAvec &
	operator+=(HaRPAvec & rpav);

	friend vector<HaRPAvec> &
	operator+=(vector<HaRPAvec> & left,  vector<HaRPAvec> & right);

	HaRPAvec &
	operator-=(HaRPAvec & rpav);

	friend vector<HaRPAvec> &
	operator-=(vector<HaRPAvec> & left, vector<HaRPAvec> & right);

	HaQCMod* GetpHost();

	const bool Get_AO_dens(HaMat_double & dens, const int imat);
	
	HaMat_double Z_mat; //!< matricies of RPA Z and Y densities in MO basis 
	HaMat_double Y_mat; // 
	                     
	HaGrpOperID id;

protected:
	
	HaQCMod* phost; 

};




#endif /* !HALOCEXCIT_H */
