/*! haintengine.cpp
 
    Classes to calculate elemenatal integrals on basis function in HARLEM
    and make summations over densities 
 
    Implementation
 
    \author Igor Kurnikov 
    \date 1998-2004
*/

#include <mpi.h>

#include "haintengine.h"
#include "haqchem.h"
#include "haatombasis.h"
#include "gaufile.h"

HaIntEngine::HaIntEngine()
{
	ptr_qc_mod = NULL;
}

HaIntEngine::HaIntEngine(const HaQCMod & qc_mod)
{
	ptr_qc_mod = &qc_mod;
}

HaIntEngine::~HaIntEngine()
{

}

int HaIntEngine::Gen2eInt(HaMat_double & R1, HaMat_double & R2, HaMat_double & R3,
					  const TwoElComb &Type2e) const
{
#if defined(GAUSSVER)

	HaQCMod::LoadGauCom(ptr_qc_mod->AtBasis);
	int nb = ptr_qc_mod->GetNBfunc();
	
	int nbasis= nb;
	int ntt= nbasis*(nbasis+1)/2;
	int nbas6d; // size of the basis set with cartesian based d and f 
	                // functions
	getnb6_(&nbas6d);
	int ntt6d=nbas6d*(nbas6d+1)/2;
	
// Axxiliary constants for fofdir

	int c__0=0, c__1=1;	
	logical c_true=1, c_false=0;
    double zero= 0.0, one = 1.0;
	double xx=0;
	int iaxx;

// make copy of PRISM parameters (ipflag and allowp[50]) to be used in 
// in non-const Fortran subroutine

	int ipflag=0;
	logical allowp[50];
	int i;
	for(i=0; i< 50; i++)
		allowp[i]= 0;

	int iprtf=0;
	int nmat=1;
	int nsing=1;
	int ntripl=0;

	int iopcl=0;
		
	int Accuracy=30; // Accurate 2-e integrals 
	int icntrl=600+Accuracy;
	
	// allocate memory for fofdir
	int lwork=4000000; 
	double* pwork;
	pwork=(double*)malloc(4000000*sizeof(double));
    assert(pwork != NULL); 

	R3.newsize(ntt,ntt);
	R3=0.0;
    
	int iout = 6;
	int ihmeth = 1;

	fofdir_(&iout, &iprtf, &ihmeth,
		    &iopcl, &icntrl, 
		    &c__0, &ipflag, &allowp[0], 
		    &c__0, &c__0, 
		    &c_false, &c_true, &c_true, 
		    &zero, &one, 
		    &nmat, &nsing, &ntripl, 
		    &nbasis, 
		    &c__0, &c__1, &c__1,       // symmetry 
	    	&xx, &xx, &xx, &xx,        // is not used 
		    &xx, &xx,                  // ha is not added
		    &xx, &xx, &xx, &xx, 
		    &xx, &xx, &xx, &xx, 
		    &c__0, &c__0, &xx, &c_false, &c__0,     // coordinates are not used
			&c_false, &xx, &xx,
			&c__1, &iaxx,
		R3.begin(), &xx, &xx, 
		&lwork, pwork);

	int itype;

	if(Type2e == REG)
	{
		
	}
	else if(Type2e == RAFF1)
	{
		// load Raffenetti 1 integrals
		R1.newsize(ntt,ntt);
		R1=0.0;
		itype=0;
        regraf_(&itype,&nbasis,&R3(1,1),&R1(1,1),&xx,&xx,pwork,&lwork);

	}


	else if(Type2e == RAFF2)
	{
		// load Raffenetti 1 integrals
		R1.newsize(ntt,ntt);
		R2.newsize(ntt,ntt);
		R1=0.0; R2=0.0;
		itype=1;
        regraf_(&itype,&nbasis,&R3(1,1),&R1(1,1),&R2(1,1),&xx,pwork,&lwork);		
	}
	
	else if(Type2e == RAFF3)
	{
		// Load Raffenetti 1,3 integrals to R1, R2 and R3
		R1.newsize(ntt,ntt);
		R1=0.0; 
		itype=2;
        regraf_(&itype,&nbasis,&R3(1,1),&R1(1,1),&xx,&R3(1,1),pwork,&lwork);		
	}

	free(pwork);

#endif

	return TRUE;
}


int HaQCMod::LoadGauCom(const GauBasisSet& gbas)
//! Load necessary Gaussian commons for Gaussian basis set gbas
{
#if defined(GAUSSVER) 
	gbas.LoadToGaussianBCommon();
	io_.in=5;
	io_.iout=6;
	io_.ipunch=9;	
	kjob_.numprc=1;
#endif
	return true;
}

int HaIntEngine::IntCanonToSq(double *pcan, HaMat_double & SqInt, int N)
//! convert 2-e integrals from Gaussian canonical form to square form 
//! (ij|kl) as  ( (N*(N+1)/1) X (N*(N+1)/1) ) matrix
{
	if(N <= 0) return 0;
	if(pcan == NULL) return 0;
	size_t ntt=N*(N+1)/2;
    SqInt.newsize(ntt, ntt);

	size_t nn=0;
	size_t i,j,k,l;
	for( i=1; i <= N; i++)
	{
		for( j=1; j < i; j++)
		{
			for(k=1; k < i; k++)           // i != j != k 
			{
				for(l=1; l <= k; l++)
				{
					SqInt(ij_indx1(i,j), ij_indx1(k,l))=pcan[nn];
					SqInt(ij_indx1(k,l), ij_indx1(i,j))=pcan[nn];
					nn++;
				}
			}
			k=i;
			for(l=1; l <= j; l++)          // i == k != j 
			{
				SqInt(ij_indx1(i,j), ij_indx1(k,l))=pcan[nn];
				SqInt(ij_indx1(k,l), ij_indx1(i,j))=pcan[nn];
				nn++;
			}
		}
		j=i;
		for(k=1; k <= i; k++)               // i == j  
		{
			for(l=1; l <= k; l++)
			{
				SqInt(ij_indx1(i,j), ij_indx1(k,l))=pcan[nn];
				SqInt(ij_indx1(k,l), ij_indx1(i,j))=pcan[nn];
				nn++;
			}
		}
	}
	
	return 1;
}

int HaBasisSet::CalcOvlpMat(const HaBasisSet* pbas1, const HaBasisSet* pbas2, HaMat_double& ovlp_mat)
{
	if( pbas1 == NULL || pbas2 == NULL) return FALSE;
	std::string bas_type_1 = pbas1->GetClassName();
	std::string bas_type_2 = pbas2->GetClassName();

	if(bas_type_1 == "GauBasisSet" || bas_type_2 == "GauBasisSet" )
	{
		return GauBasisSet::CalcOvlpMat((GauBasisSet*)pbas1,(GauBasisSet*)pbas2,ovlp_mat);
	}
	if(bas_type_1 == "LinCombOrb3D" || bas_type_2 == "LinCombOrb3D" )
	{
        return LinCombOrb3D::CalcOvlpMat((LinCombOrb3D*)pbas1, (LinCombOrb3D*)pbas2, ovlp_mat);
	}
	return FALSE;
}


//int
//HaIntEngine::Eval1eOp(HaBasisSet* pbas1, HaBasisSet* pbas2, HaMat_double& oper_mat,int oper_type)
//{
//	if( pbas1 == NULL || pbas2 == NULL) return FALSE;
//	std::string bas_type_1  = pbas1->GetClassName();
//	std::string bas_type_2 = pbas2->GetClassName();
//
//	if(bas_type_1 == "GauBasisSet" || bas_type_2 == "GauBasisSet" )
//	{
//		return GauBasisSet::Eval1eOp((GauBasisSet*)pbas1,(GauBasisSet*)pbas2,oper_mat);
//	}
//	if(bas_type_1 == "LinCombOrb3D" || bas_type_2 == "LinCombOrb3D" )
//	{
//      return LinCombOrb3D::Eval1eOp((LinCombOrb3D*)pbas1, (LinCombOrb3D*)pbas2, oper_mat);
//	}
//	return FALSE;
//}

