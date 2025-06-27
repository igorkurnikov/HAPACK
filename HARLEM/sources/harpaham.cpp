/*! \file harpaham.cpp

    Classes to define RPA Hamiltonian in HARLEM.

    \author Igor Kurnikov 
    \date 1998-2002

*/

#define HARPAHAM_CPP

#include <mpi.h>

#include "halinalg.h"
#include "halocexcit.h"
#include "haqchem.h"
#include "harpaham.h"
#include "g94_globals.h"
#include "g94_protos.h"


HaRPAHam::HaRPAHam()
{
	Ene=0.0;
	opmode=FULL;
}



HaRPAHam::~HaRPAHam()
{

}

bool
HaRPAHam::SetOpMode(int imode)
{
	assert(imode >= 0 && imode < 4);
	opmode=(OperMode)imode;
	return true;
}

bool
HaRPAHam::SetEnergy(const double NewEne)
{
	Ene=NewEne;
	return true;
}


bool
HaRPAHam::Apply_init(std::vector<HaRPAvec> & RPAvec) const
{

	int nvec= RPAvec.size(); 
	assert(nvec > 0);
	int no=RPAvec[0].GetNumOccMO();
    int nv=RPAvec[0].GetNumVacMO();
	
	const HaQCMod* phost = RPAvec[0].GetpHost();
	if(phost == NULL)
	{
		ErrorInMod("HaRPAHam::Apply_init()","RPA vector doesn't has a host module");
		return false;
	}
   	if(phost->MOene.size() == 0)
	{
		ErrorInMod("HaRPAHam::Apply_init()","MOs are not set");
		return false;	    
	}
	HaVec_double Evac(nv);
	HaVec_double Eocc(no); 

	int i,a;
	for( a =1; a <= nv; a++)
		Evac(a)= phost->MOene(no+a) ;
	for( i =1; i <= no; i++)
		Eocc(i)= phost->MOene(i) ;

	for (int iv=0; iv < nvec; iv++)
	{
		for(int a=1; a <= nv; a++)
		{
			for(int i=1; i <= no; i++)
			{
				RPAvec[iv].Z_mat(a,i) = RPAvec[iv].Z_mat(a,i)*(Evac(a) - Eocc(i));
				RPAvec[iv].Y_mat(a,i) = RPAvec[iv].Y_mat(a,i)*(Evac(a) - Eocc(i));
			}
		}
	}
	return true;
}

bool HaRPAHam::Apply(std::vector<HaRPAvec> & RPAvec)
// Calculate  Z'= A Z + B Y  
//            Y'= B Z + A Y 
{
#if defined(GAUSSVER) 

	int c__0=0, c__1=1;
	
	logical c_true=1, c_false=0;
    double zero= 0.0, one = 1.0;

	HaQCMod* phost = RPAvec[0].GetpHost();
	
	HaQCMod::LoadGauCom(phost->AtBasis);
	
	int nvec=RPAvec.size();
	assert(nvec > 0);
	int no=RPAvec[0].GetNumOccMO();
    int nv=RPAvec[0].GetNumVacMO();

	
	int nmat= 4*nvec;
	int nsing = 2*nvec;
	int ntripl = 0;
	int iraf= -1;
	int nbasis = phost->GetNBfunc();
	int nttfd = nbasis*nbasis;
	int nb = nbasis;
	int ntt= nbasis*(nbasis+1)/2;
	int nbas6d; // size of the basis set with cartesian based d and f 
	                // functions

	getnb6_(&nbas6d);
	int ntt6d=nbas6d*(nbas6d+1)/2;

	double xx=0;
// make copy of PRISM parameters (ipflag and allowp[50]) to be used in 
// in non-const Fortran subroutine

	int ipflag=0;
	logical allowp[50];
	int i;
	for(i=0; i< 50; i++)
		allowp[i]= 0;

	int ipflag_copy=ipflag;
	logical allowp_copy[50];
	for(i=0; i< 50; i++)
		allowp_copy[i]=allowp[i];
		
    HaVec_double v_td(4*nvec*ntt6d); // allocation 4*nv*ntt6d is necessary 
                                       // for fofdir
	         
	HaMat_double tdsymm;
	tdsymm.set_ext_alloc(v_td.begin(),ntt,2*nvec);
 	// HaMat_double tdsymm(ntt,2*nvec,v_td.begin(),1); 
	// array of symmetrical parts of Z and Y

	HaMat_double tdasymm;
	tdsymm.set_ext_alloc(&v_td(2*nvec*ntt+1),ntt,2*nvec);
    //HaMat_double tdasymm(ntt,2*nvec,&v_td(2*nvec*ntt+1),1); 
	// array of asymmetrical parts of Z and Y


    HaVec_double v_f(4*nvec*ntt6d);
	HaMat_double fa;
	fa.set_ext_alloc(v_f.begin(),ntt,2*nvec);
	HaMat_double fb;
	fa.set_ext_alloc(&v_f(2*nvec*ntt+1),ntt,2*nvec);

	// allocate memory for fofdir
	int lwork= HaQCMod::max_gauss_mem; 
	double* pwork;
	pwork=(double*)malloc(HaQCMod::max_gauss_mem*sizeof(double));
        assert(pwork != NULL); 

	HaMat_double C_occ;
	C_occ.set_ext_alloc(&(phost->MO_coef(1,1)),nb,no);
	HaMat_double C_vac;
	C_vac.set_ext_alloc(&(phost->MO_coef(1,no+1)),nb,nv);

	HaMat_double scr;
	HaMat_double X_AO(nb,nb);
	HaMat_double X_MO(nv,no);

	int iv;
	for(iv=0; iv < nvec; iv++)
	{
		matmult(scr,C_vac,RPAvec[iv].Z_mat);
		matmult_T2(X_AO,scr,C_occ);
//		cout << " X_AO " << endl;
//		X_AO.Print_format(cout,"%7.4f ");
		X_AO.GetSymmPart(&tdsymm(1,2*iv+1));
		X_AO.GetASymmPart(&tdasymm(1,2*iv+1));

		matmult(scr,C_vac,RPAvec[iv].Y_mat);
		matmult_T2(X_AO,scr,C_occ);

		X_AO.GetSymmPart(&tdsymm(1,2*iv+2));
		X_AO.GetASymmPart(&tdasymm(1,2*iv+2));
	}
//		cout << " tdsymm " << endl;
//		tdsymm.Print_format(cout,"%7.4f ");
//		cout << " tdasymm " << endl;
//		tdasymm.Print_format(cout,"%7.4f ");
	
	int iprtf=0, iopcl=0, icntrl=230;
	int iout = 6;
	int ihmeth = 1;
	int iaxx;


	fofdir_(&iout, &iprtf, &ihmeth, 
		    &iopcl, &icntrl, 
		&iraf, &ipflag_copy, &allowp_copy[0], 
		&c__0, &c__0, 
		&c_false, &c_true, &c_true, 
		&zero, &one, 
		&nmat, &nsing, &ntripl, 
		&nbasis, 
		&c__0, &c__1, &c__1,       // symmetry 
		&xx, &xx, &xx, &xx,        // is not used 
		&xx, &xx,                  // ha is not added
		&tdsymm(1,1), &tdsymm(1,1), &tdsymm(1,1), &tdsymm(1,1), 
		&fa(1,1), &fa(1,1), &xx, &xx, 
		&c__0, &c__0, &xx, &c_false, &c__0,       // coordinates are not used 
		&c_false, &xx, &xx,
		&c__1, &iaxx,
		&xx, &xx, &xx, 
		&lwork, pwork);
	
	free(pwork);
	v_td.newsize(0);

//		cout << " fa " << endl;
//		fa.Print_format(cout,"%7.4f ");
//		cout << " fb " << endl;
//		fb.Print_format(cout,"%7.4f ");
  

	this->Apply_init(RPAvec);     // MO diagonal approximation for RPA hamiltonian

	for(iv=0; iv < nvec; iv++)
	{
		HaMat_double AmB_Z(nb,nb); // (A-B)*Z
		AmB_Z=0.0;  
		AmB_Z.AddSymmPart(&fa(1,2*iv+1)); //FoFdir calculates 1/2 of (A(1) - B) X 
//		cout << " AmB_Z " << endl;
//		AmB_Z.Print_format(cout,"%6.3f ");

		HaMat_double ApB_Z(nb,nb);  // (A+B)*Z 
		ApB_Z=0.0;
		ApB_Z.AddASymmPart(&fb(1,2*iv+1));

		HaMat_double AmB_Y(nb,nb);  // (A-B)*Y
		AmB_Y=0.0;
		AmB_Y.AddSymmPart(&fa(1,2*iv+2));

		HaMat_double ApB_Y(nb,nb);  // (A+B)*Y  
		ApB_Y=0.0;
		ApB_Y.AddASymmPart(&fb(1,2*iv+2));

		X_AO=0.0;
		mat_add(X_AO,AmB_Z,ApB_Z);
		mat_add(X_AO,X_AO,ApB_Y);
		mat_diff(X_AO,X_AO,AmB_Y);

   	    matmult_T1(scr,C_vac,X_AO);
	    matmult   (X_MO,scr,C_occ);

		mat_add(RPAvec[iv].Z_mat,RPAvec[iv].Z_mat,X_MO);

		X_AO=0.0;
		mat_add(X_AO,AmB_Y,ApB_Y);
		mat_add(X_AO,X_AO,ApB_Z);
		mat_diff(X_AO,X_AO,AmB_Z);

   	    matmult_T1(scr,C_vac,X_AO);
	    matmult   (X_MO,scr,C_occ);

		mat_add(RPAvec[iv].Y_mat,RPAvec[iv].Y_mat,X_MO);

	}
#else
	PrintLog("Error in HaRPAHam::Apply() \n");
	PrintLog("GAUSSIAN Library is not available \n");
#endif

	return 1;
}

std::vector<HaRPAvec>
HaRPAHam::operator *(const std::vector<HaRPAvec> & RPAv)
{
	std::vector<HaRPAvec> RPAv_new(RPAv);
	HaMat_double tmpz,tmpy;
	int iv,nv;
	
	switch(opmode)
	{
	case FULL:
		(*this).Apply(RPAv_new);
		break;
	case MO_DIAG:
		(*this).Apply_init(RPAv_new);
		break;
	case E_MIN_H:
		
		(*this).Apply(RPAv_new);
		nv=RPAv.size();
		for( iv=0; iv < nv; iv++)
		{
			mat_copy_scale(tmpz, RPAv[iv].Z_mat,Ene);
			mat_copy_scale(tmpy, RPAv[iv].Y_mat,-Ene);
			
			mat_diff(RPAv_new[iv].Z_mat,tmpz,RPAv_new[iv].Z_mat);
			mat_diff(RPAv_new[iv].Y_mat,tmpy,RPAv_new[iv].Y_mat);
			
		}
		break;
	case E_MIN_H0:
		
		(*this).Apply_init(RPAv_new);
		nv=RPAv.size();
		for( iv=0; iv < nv; iv++)
		{
			mat_copy_scale(tmpz, RPAv[iv].Z_mat,Ene);
			mat_copy_scale(tmpy, RPAv[iv].Y_mat,-Ene);
			
			mat_diff(RPAv_new[iv].Z_mat,tmpz,RPAv_new[iv].Z_mat);
			mat_diff(RPAv_new[iv].Y_mat,tmpy,RPAv_new[iv].Y_mat);
		}
		break;
	default:
		;
	}
	return RPAv_new;
}


HaRPAResolv::HaRPAResolv()
{
	Ene=0.0;
	opmode=FULL;
	imag=false;
}


HaRPAResolv::~HaRPAResolv()
{

}

bool
HaRPAResolv::SetOpMode(const int imode)
{
	assert(imode >= 0 && imode < 2);
	if(imode == 0)
	{
		opmode=FULL;
		return true;
	}
	else if(imode == 1)
	{
		opmode=MO_DIAG;
		return true;
	}
	return false;
}

bool
HaRPAResolv::SetEnergy(const double NewEne)
{
	Ene=NewEne;
	return true;
}

bool
HaRPAResolv::SetImag(const bool new_imag)
{
	imag=new_imag;
	return true;
}


bool
HaRPAResolv::Apply(std::vector<HaRPAvec> & a1) const
{
	std::vector<HaRPAvec> x(a1);
	Apply_Init(x);
	
	HaRPAHam  h1;
	h1.SetEnergy(this->Ene);
	h1.SetOpMode(2);
	
	HaRPAResolv g0;
	g0.SetOpMode(1);
	g0.SetEnergy(this->Ene);

	std::cout << " HaRPAResolv::Apply " << std::endl;
	
	int maxiter=100;
	HaVec_double tol;
	tol.newsize(a1.size());
	tol=0.00000001;
	int result=CG_mult(h1,x,a1,g0,maxiter,tol);

//	cout << " Second call to CG_mult " << endl;
//	result=CG_mult(h1,x,a1,g0,maxiter,tol);


	std::cout << " Final accuracy is tol= " << tol << std::endl;
	std::cout << " The required number of iterations:  " << maxiter << std::endl;
	
	a1=x;
	
	return true;
}

bool
HaRPAResolv::Apply_Init(std::vector<HaRPAvec> & RPAvec) const
{
	int nvec= RPAvec.size(); 
	assert(nvec > 0);
	int no=RPAvec[0].GetNumOccMO();
    int nv=RPAvec[0].GetNumVacMO();
	
	const HaQCMod* phost = RPAvec[0].GetpHost();
	if(phost == NULL)
	{
		ErrorInMod("HaRPAResolv::Apply_init()","RPA vector doesn't has a host module");
		return false;
	}
   	if(phost->MOene.size() == 0)
	{
		ErrorInMod("HaRPAResolv::Apply_init()","MOs are not set");
		return false;	    
	}

	HaVec_double Evac(nv);
	HaVec_double Eocc(no); 

	int i,a;
	for( a =1; a <= nv; a++)
		Evac(a)= phost->MOene(no+a) ;
	for( i =1; i <= no; i++)
		Eocc(i)= phost->MOene(i) ;

	for (int iv=0; iv < nvec; iv++)
	{
		for(int a=1; a <= nv; a++)
		{
			for(int i=1; i <= no; i++)
			{		
				RPAvec[iv].Z_mat(a,i) = RPAvec[iv].Z_mat(a,i)/(Ene-(Evac(a) - Eocc(i)));
				RPAvec[iv].Y_mat(a,i) = RPAvec[iv].Y_mat(a,i)/(-Ene-(Evac(a) - Eocc(i)));
			}
		}
	}
	return true;
}

std::vector<HaRPAvec>
HaRPAResolv::solve(const std::vector<HaRPAvec> & RPAv) const
{
	std::vector<HaRPAvec> RPAv_new(RPAv);
	if(opmode == FULL)
	{
		(*this).Apply(RPAv_new);
	}
	else if(opmode == MO_DIAG)
	{
		(*this).Apply_Init(RPAv_new);
	}
	return RPAv_new;
}
