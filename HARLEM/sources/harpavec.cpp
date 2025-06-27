/*!  harpavec.cpp

    Classes to define RPA vector

    \author Igor Kurnikov 
	\date 1998-2002

*/

#define HARPAVEC_CPP

#include <mpi.h>

#define exception math_exception
#include <math.h>
#undef exception

#include "hamultipole.h"
#include "halocexcit.h"
#include "harpavec.h"
#include "hamolecule.h"
#include "haqchem.h"


HaRPAvec::HaRPAvec()
{
	phost=NULL;
}

HaRPAvec::HaRPAvec(HaQCMod & host)
{
	phost=&host;

}
	
HaRPAvec::HaRPAvec(const HaLocExcit & excit,HaQCMod & host)
{
	phost=&host;

}
	

HaRPAvec::~HaRPAvec()
{

}

HaQCMod* 
HaRPAvec::GetpHost() 
{
	return(phost);
}

int 
HaRPAvec::GetNBfunc() const
{
	assert(phost != NULL);
	return(phost->GetNBfunc());
}


int 
HaRPAvec::GetNumOccMO() const
{
	assert(phost != NULL);
	return(phost->GetNumOccMO());
}

int 
HaRPAvec::GetNumVacMO() const
{
	assert(phost != NULL);
	return(phost->GetNumVacMO() );
}

int
HaRPAvec::SetFromLocExcit(const HaLocExcit & excit)
{	
	assert(phost != NULL);

	int nb=this->GetNBfunc();
	int no=this->GetNumOccMO();
	int nv=this->GetNumVacMO();

    HaMat_double AODens(nb,nb);
	excit.Expand(AODens.begin(),(*phost));

	this->SetFromAOMat(AODens);

	return 1;
}


int
HaRPAvec::SetFromAOMat(const HaMat_double & aomat,
					   const OPER_TYPE optyp)
{	
	if(phost == NULL) return FALSE;

	int nb=this->GetNBfunc();
	int no=this->GetNumOccMO();
	int nv=this->GetNumVacMO();

	if(aomat.num_rows() != nb) return FALSE;
	if(aomat.num_cols() != nb) return FALSE;
	if(phost->MO_coef.num_cols() == 0) return FALSE;

	HaMat_double& ss = phost->GetOvlpMat();
	HaMat_double SC_occ;
	HaMat_double scr_occ;
	scr_occ.set_ext_alloc(phost->MO_coef.begin(),nb,no);
	matmult(SC_occ, ss,scr_occ);

	HaMat_double SC_vac;
	HaMat_double scr_vac;
	scr_vac.set_ext_alloc(phost->MO_coef.begin(),nb,nv);
	matmult(SC_vac, ss ,scr_vac);

    HaMat_double scr;
	 
	matmult_T1(scr,SC_vac,aomat);
	matmult   (Z_mat,scr,SC_occ);
	Y_mat.newsize(nv,no);
	if(optyp == REAL_OPER)
		mat_scale(Y_mat,Z_mat,-1.0);
	else if( optyp == IMAG_OPER)
		Y_mat=Z_mat;
	return 1;
}

int HaRPAvec::SetFromLOGrpMat(const std::string& gid1,const std::string& gid2, 
					    const HaMat_double & fmloc, const OPER_TYPE optyp) 
{
	int nb=phost->GetNBfunc();		

	int i;

	HaVec_int loindx1;
	phost->GetLocOrbIdxOfGrp(gid1,loindx1);
	int nlorb1=loindx1.size();
	
	HaVec_int loindx2;
	phost->GetLocOrbIdxOfGrp(gid2,loindx2);
	int nlorb2=loindx2.size();
	
	assert(nlorb1 == fmloc.num_rows() && nlorb2 == fmloc.num_cols());

	// building AO basis vectors of local orbitals:

	HaMat_double td; // transitional density matrix in atomic basis set functions
	
	std::string bas_type = phost->ActBas->GetClassName();

	if(bas_type == "LinCombOrb3D" )
	{
		LinCombOrb3D* lcmb = (LinCombOrb3D*) phost->ActBas;

		HaMat_double vl2(nb,nlorb2);
		int j;
		for(i=0; i < nlorb2; i++)
		{
			for(j =0; j < nb; j++)
			{
				vl2.SetVal_idx0(j,i, lcmb->coef.GetVal_idx0(j,loindx2[i]));
			}
		}
		HaMat_double vl1(nb,nlorb1);

		for(i=0; i < nlorb1; i++)
		{
			for(j =0; j < nb; j++)
			{
				vl1.SetVal_idx0(j,i, lcmb->coef.GetVal_idx0(j,loindx1[i]));
			}
		}

		HaMat_double scr;
		matmult(scr,vl1,fmloc);
		matmult_T2(td,scr,vl2);
	}
	else
	{
		td = fmloc;
	}
	// Add very small mumber to avoid zero density

	double thresh=1.0e-15;
	if( fabs(td(1,1)) < thresh) 
		td(1,1)=thresh;


	SetFromAOMat(td,optyp);

	return True;
}

double
dot2(const HaRPAvec & left, const HaRPAvec & right)
{
	int no= left.GetNumOccMO();
	int nv= left.GetNumVacMO();
	assert( no == right.GetNumOccMO() );
	assert( nv == right.GetNumVacMO() );
	double ss=0.0;
	int i,a;
    for(a=1; a <= nv; a++)
	{
		for(i=1; i <= no; i++)
		{
			ss+= left.Z_mat(a,i) * right.Z_mat(a,i);
			ss+= left.Y_mat(a,i) * right.Y_mat(a,i);
		}
	}
	return ss;
}

HaVec_double
dot2(const std::vector<HaRPAvec> & left, const std::vector<HaRPAvec> & right)
{
	int nvec=left.size();
	assert(nvec == right.size());
	HaVec_double sv(nvec);
	for(int iv=0; iv < nvec; iv++)
	{	
		sv[iv]=dot2(left[iv],right[iv]);
	}
	return sv;
}

double norm2(const HaRPAvec & RPAv)
{
	return(dot2(RPAv,RPAv));
}

HaVec_double
norm2(const std::vector<HaRPAvec> & RPAv_arr)
{
	return(dot2(RPAv_arr,RPAv_arr));
}

HaMat_double
SProd(const std::vector<HaRPAvec> & left, const std::vector<HaRPAvec> & right)
{
	int nvec1=left.size();
	int nvec2=right.size();
	assert( (nvec1 > 0) && (nvec2 > 0) );
	HaMat_double  ss(nvec1,nvec2);
	ss=0.0;
	for (int iv1=0; iv1 < nvec1; iv1++)
	{
		for (int iv2=0; iv2 < nvec2; iv2++)
		{
			ss(iv1+1,iv2+1)=  SProd(left[iv1].Z_mat,right[iv2].Z_mat); 
			ss(iv1+1,iv2+1)+= SProd(left[iv1].Y_mat,right[iv2].Y_mat); 
		}
	}
	return ss;
}

HaVec_double
SProd(const HaRPAvec & RPAv, const std::vector<HaRPAvec> & RPAv_arr2)
{
	std::vector<HaRPAvec> RPAv_arr;
	RPAv_arr.push_back(RPAv);
	HaMat_double ss_mat=SProd(RPAv_arr,RPAv_arr2);
	HaVec_double ss_vec(ss_mat.num_rows());
	for(int i=1; i <= ss_vec.size(); i++)
		ss_vec(i)=ss_mat(1,i);
	return(ss_vec);
}


HaRPAvec operator*(const double factor, const HaRPAvec & RPAv)
{
	HaRPAvec RPAv_new(RPAv);
	mat_scale(RPAv_new.Z_mat,RPAv_new.Z_mat,factor);
	mat_scale(RPAv_new.Y_mat,RPAv_new.Y_mat,factor);
	return RPAv_new;
}

std::vector<HaRPAvec> operator*(const HaVec_double & vfactor, const std::vector<HaRPAvec> & RPAv_arr)
{
	std::vector<HaRPAvec> RPAv_arr_new(RPAv_arr);
	int nvec=RPAv_arr.size();
	for(int iv=0; iv < nvec; iv++)
	{
		mat_scale(RPAv_arr_new[iv].Z_mat,RPAv_arr_new[iv].Z_mat,vfactor[iv]);
		mat_scale(RPAv_arr_new[iv].Y_mat,RPAv_arr_new[iv].Y_mat,vfactor[iv]);		
	}
	return(RPAv_arr_new);
}
 


HaRPAvec operator+(HaRPAvec & left, HaRPAvec & right)
{
	assert(left.GetpHost() != NULL);
	assert(left.GetpHost() == right.GetpHost() );
	HaRPAvec RPAv(left);
	mat_add( RPAv.Z_mat,RPAv.Z_mat,right.Z_mat);
	mat_add( RPAv.Y_mat,RPAv.Y_mat,right.Y_mat);
	return RPAv;
}

std::vector<HaRPAvec> operator+(std::vector<HaRPAvec> & left, std::vector<HaRPAvec> & right)
{
	int nvec=left.size();
	assert(nvec > 0);
	assert(nvec == right.size());
	assert(left[0].GetpHost() == right[0].GetpHost());
	std::vector<HaRPAvec> RPAv_arr(left);
	for(int iv=0; iv < nvec; iv++)
	{
		mat_add(RPAv_arr[iv].Z_mat,RPAv_arr[iv].Z_mat, right[iv].Z_mat);
		mat_add(RPAv_arr[iv].Y_mat,RPAv_arr[iv].Y_mat, right[iv].Y_mat);
	}
	return RPAv_arr;
}


HaRPAvec operator-(HaRPAvec & left, HaRPAvec & right)
{
	assert(left.GetpHost() != NULL);
	assert(left.GetpHost() == right.GetpHost() );
	HaRPAvec RPAv(left);
	mat_diff( RPAv.Z_mat,RPAv.Z_mat,right.Z_mat);
	mat_diff( RPAv.Y_mat,RPAv.Y_mat,right.Y_mat);
	return RPAv;
}

std::vector<HaRPAvec> operator-(std::vector<HaRPAvec> & left, std::vector<HaRPAvec> & right)
{
	int nvec=left.size();
	assert(nvec > 0);
	assert(nvec == right.size());
	assert(left[0].GetpHost() == right[0].GetpHost());
	std::vector<HaRPAvec> RPAv_arr(left);
	for(int iv=0; iv < nvec; iv++)
	{
		mat_diff(RPAv_arr[iv].Z_mat,RPAv_arr[iv].Z_mat, right[iv].Z_mat);
		mat_diff(RPAv_arr[iv].Y_mat,RPAv_arr[iv].Y_mat, right[iv].Y_mat);
	}
	return RPAv_arr;
}

HaRPAvec & HaRPAvec::operator+=(HaRPAvec & right)
{
	assert(this->GetpHost() != NULL);
	assert(this->GetpHost() == right.GetpHost() );
	mat_add( this->Z_mat,this->Z_mat,right.Z_mat);
	mat_add( this->Y_mat,this->Y_mat,right.Y_mat);
	return (*this);
}


HaRPAvec & HaRPAvec::operator-=(HaRPAvec & right)
{
	assert(this->GetpHost() != NULL);
	assert(this->GetpHost() == right.GetpHost() );
	mat_diff( this->Z_mat,this->Z_mat,right.Z_mat);
	mat_diff( this->Y_mat,this->Y_mat,right.Y_mat);
	return (*this);
}


std::vector<HaRPAvec> & operator+=(std::vector<HaRPAvec> & left, std::vector<HaRPAvec> & right)
{
	int nvec=left.size();
	assert(nvec > 0);
	assert(nvec == right.size());
	assert(left[0].GetpHost() == right[0].GetpHost());

	for(int iv=0; iv < nvec; iv++)
	{
		mat_add(left[iv].Z_mat,left[iv].Z_mat, right[iv].Z_mat);
		mat_add(left[iv].Y_mat,left[iv].Y_mat, right[iv].Y_mat);
	}
	return left;
}

std::vector<HaRPAvec> & operator-=(std::vector<HaRPAvec> & left, std::vector<HaRPAvec> & right)
{
	int nvec=left.size();
	assert(nvec > 0);
	assert(nvec == right.size());
	assert(left[0].GetpHost() == right[0].GetpHost());

	for(int iv=0; iv < nvec; iv++)
	{
		mat_diff(left[iv].Z_mat,left[iv].Z_mat, right[iv].Z_mat);
		mat_diff(left[iv].Y_mat,left[iv].Y_mat, right[iv].Y_mat);
	}
	return left;
}





