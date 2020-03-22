/*! \file halocexcit.cpp

    Classes to define local excitation in HARLEM

    \author  Igor Kurnikov  
    \date 1997-2002
    
*/

#define HALOCEXCIT_CPP

#include <mpi.h>

#include "haqchem.h"
#include "hamultipole.h"
#include "halocexcit.h"
#include "hamolecule.h"
#include "hamolset.h"

#include "haatgroup.h"

#define exception math_exception
#include <math.h>
#undef exception


HaLocExcit::HaLocExcit()
{
	PartOrb = NULL;
	HoleOrb = NULL;
}

HaLocExcit::HaLocExcit(ArrayOrb3D & p1, ArrayOrb3D & p2)
{
	PartOrb = &p1;
	HoleOrb = &p2;
}

HaLocExcit::~HaLocExcit()
{

}

int
HaLocExcit::Expand(double* dmatr, HaQCMod & host) const
// Function expand local excitation to Nbasis x Nbasis Transition density matrix
// corresponding to the basis of the molecule host
//
// T_0= C_1 X C_2^T
// dmatr= (I - d_AO*S) * T_0 * S * d_AO
//
{
	int nb=host.GetNBfunc();
	if(nb <= 0) return false;

	int i,j;

	for(i=nb*nb; i !=0 ; i--)
	{
		dmatr[i-1]=0.0;
	}
	
	HaMat_double cf_p;
	HaMat_double cf_h;

	PartOrb->ExpandInBas(cf_p,host.AtBasis);
	HoleOrb->ExpandInBas(cf_h,host.AtBasis);

	for(i = 0 ; i < nb; i++)
	{
		int jbeg = nb*i;
		for(j = 0; j < nb; j++)
		{
			dmatr[jbeg + j] = cf_p.GetVal_idx0(i,0)*cf_h.GetVal_idx0(j,0);
		}
	}

//	HaMat_double dm(nb,nb,dmatr,1);

//	HaMat_double Dens1e=(HaMat_double)(*(host.pDens1e));  // d_AO
//	mat_scale(Dens1e,Dens1e,0.5); // Scale Density matrix because of 2 spin
//	HaMat_double ss =   host.ovlp_map;

//	HaMat_double scr1,scr2;
//	matmult(scr1,dm,ss);      
//	matmult(scr2,scr1,Dens1e);  // scr2 = T_0 * S * d_AO 
//	matmult(scr1,Dens1e,ss);       
//	mat_scale(scr1,scr1,-1.0); 
//	mat_add_unit(scr1, scr1, 1.0); // scr1 = I - d_AO*S
//	matmult(dm,scr1,scr2);         // d_tr^AO= (I-d_AO S) T_0 S d_AO 
    
	return true;	

}

