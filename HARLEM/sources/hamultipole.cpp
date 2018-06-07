/*! \file hamultipole.cpp
    
    implementation of classes defining 1e operators
    Electrical Multipole operators and 
    magnetic dipole operators
 
    \file Igor Kurnikov , University of Pittsburgh 
    \date 1998-2002
*/
#include <mpi.h>

#include "g94_globals.h"
#include "g94_protos.h"

#include "haqchem.h"
#include "hamultipole.h"
#include "harpaham.h"


#ifdef USE_IPACK  // IPACK headers
#include "basis.h"
#include "operators.h"
#endif

HaOperR::HaOperR()
{

}

HaOperR::~HaOperR(void)
{

}


HaOperRDelt::HaOperRDelt(HaQCMod* ptr_qc_mod_new)
{
	data.resize(3);
	i_lond=0;
	ptr_qc_mod = ptr_qc_mod_new;
}

HaOperRDelt::~HaOperRDelt(void)
{

}

int HaOperRDelt::EvalGauBasisSet(GauBasisSet* pbset, HaMat_doubleArr& rmats)
{
	int nbasis = pbset->GetNBfunc();
	if(HaQCMod::int_engine == HaQCMod::int_engine.INT_ENGINE_GAUSS)
	{
		PrintLog("Compute Angular Momentum Operator Matricies using GAUSSIAN Library \n");

#if defined(GAUSSVER) 

		int c0=0, c1=1;
		int iprint=0;
		double xx=0;
		HaMat_double coord;
		HaQCMod::LoadGauCom(*pbset);
		int nbasis = pbset->GetNBfunc();
		int nbas6d;
		getnb6_(&nbas6d);
		int ntt6d=nbas6d*(nbas6d+1)/2;

		int ntt= nbasis*(nbasis+1)/2;
		int natoms = pbset->GetNumCnt();
		pbset->GetCntCoord(coord);
		int nvec=6;
		HaVec_double scr_val( ntt6d*nvec);
		HaMat_double valmat;
	    valmat.set_ext_alloc(scr_val.begin(),ntt,nvec);

		double* pwork;
		int mdv=10*nbasis*nbasis;
		if(mdv < HaQCMod::max_gauss_mem) mdv=HaQCMod::max_gauss_mem;
	
		pwork=(double*)malloc(mdv*sizeof(double));
		assert(pwork != NULL);
		int* piwork=(int*)pwork;

		int icalc =12;        // Calc Del and R X Del
		if(i_lond == 1) icalc= 9; // Overlap, kinetic and angular momentum integrals for GIAO.
		int idgst=10; // return Operator Matricies

		int ipflag=0;
		logical allowp[50];
		int i;
		for(i=0; i< 50; i++)
			allowp[i]= 0;

		denbas_(&io_.iout,&icalc,&idgst,
				&c0,&c0,&c0,&c0,
				&c0,&c0,&c0,&xx,&xx,
				&c0,&nbasis,&c1, &c0,
				coord.begin(),&c0, &c0, &xx, &natoms,
				&xx,&c0,
				&xx,valmat.begin(),&xx,
				&c0,&c1,&c0,&xx,
	           &iprint,&ipflag,&allowp[0],&c0,
                   piwork,pwork,&mdv);
		free(pwork);

		int ishift=4;
		if(i_lond == 1) ishift=1;

		if(rmats.size() != 3) rmats.resize(3);

		for (i=0; i < 3; i++)
		{
			size_t nb= nbasis;
			data[i].newsize(nbasis,nbasis);
			HaVec_double fmat_lin(ntt,&valmat(1,i+ishift)); 
			mat_asymm_from_lin(data[i],fmat_lin,nb);
			rmats[i] = data[i];
		}
#else
		PrintLog(" GAUSSIAN Library is not available \n");
#endif
	}
	else if(HaQCMod::int_engine == QCIntEngineType::INT_ENGINE_IPACK)
	{
		PrintLog("Compute Angular Momentum Operator Matricies using IPACK Library \n");
    	HaQCMod::InitIPack();
		InternalBasis* ip_bas = pbset->CreateIPackBas();

		int nb = nbasis;
	    ARRAY<Mat> r_mat_arr(3);
		int i;
		for(i=0; i <3; i++)
		{
			r_mat_arr[0].reset(nb,nb);
			r_mat_arr[0].set(0);
		}
		int same_bas = TRUE;
		Location cnt(0.0,0.0,0.0);
		angular(*ip_bas, cnt, r_mat_arr);

		normalize(*ip_bas, *ip_bas,r_mat_arr,1);
		if(rmats.size() != 3) rmats.resize(3);
		for(i = 1; i < 4; i++)
		{
            rmats[i-1].set_from(r_mat_arr[i]);
		}
		free(ip_bas);
	}

	return true;
}


bool
HaOperRDelt::RecalcLondon(GauBasisSet* pbas)
{
#if defined(GAUSSVER) 

	int c0=0, c1=1;
	int iprint=0;
	double xx=0;
	HaMat_double coord;
	HaQCMod::LoadGauCom(*pbas);
	int nbasis = pbas->GetNBfunc();
	int nbas6d;
	getnb6_(&nbas6d);
	int ntt6d=nbas6d*(nbas6d+1)/2;

	int ntt= nbasis*(nbasis+1)/2;
	int natoms = pbas->GetNumCnt();
	pbas->GetCntCoord(coord);
	int nvec=6;
	HaVec_double scr_val( ntt6d*nvec);
	HaMat_double valmat;
	valmat.set_ext_alloc(scr_val.begin(),ntt,nvec);
	HaMat_double fangm(ntt,3);

	double* pwork;
	int mdv=10*nbasis*nbasis;
	if(mdv < 4000000) mdv=4000000;
	
	pwork=(double*)malloc(mdv*sizeof(double));
    assert(pwork != NULL);
	int* piwork=(int*)pwork;

	int icalc= 9; // Overlap, kinetic and angular momentum integrals for GIAO.
	int idgst=10; // return Operator Matricies

	int ipflag=0;
	logical allowp[50];
	int i;
	for(i=0; i< 50; i++)
		allowp[i]= 0;

    denbas_(&io_.iout,&icalc,&idgst,
		   &c0,&c0,&c0,&c0,
		   &c0,&c0,&c0,&xx,&xx,
		   &c0,&nbasis,&c1, &c0,
		   coord.begin(),&c0, &c0,&xx, &natoms,
		   &xx,&c1,                    // cgrid,ngrid
		   &xx,valmat.begin(),&xx,
		   &c0,&c1,&c0,&xx,
	       &iprint,&ipflag,&allowp[0],&c0,
           piwork,pwork,&mdv);

	HaMat_double TdB;
	TdB.set_ext_alloc(&valmat(1,4),ntt,3);
	
	mat_scale(TdB, TdB, 0.5);

#if 0
	cout << "Print Tdb matricies: " << endl;
	for(i=1; i <=3 ; i++)
	{
		cout << "TdB component " << i << endl;
		HaVec_double fmat_lin(ntt,&TdB(1,i));
		HaMat_double scr;
		mat_asymm_from_lin(scr,fmat_lin,nbasis);
		scr.Print_format(cout," %10.5f");
	}
#endif

	fangm=TdB;                // fill london angular momentum matrix 
	                         // by derivatives of kinetic energy matricies versus B  

	icalc= 10; // Potential energy integrals for GIAO. 
	iprint =3;

	natoms = pbas->GetNumCnt();
	int nb = pbas->GetNBfunc();
	Vec3D* pptr_old = NULL;
    Vec3D* pptr = NULL;

	HaVec_double tmp(nb);
	int nc = 0;
	
	for(i = 0; i < nb; i++)
	{
		pptr = 	pbas->GetHostPt(i);
		if(pptr != pptr_old)
		{	
			HaAtom* aptr = (HaAtom*) pptr;
			tmp[nc] = aptr->GetCharge();
			nc++;
			pptr_old = pptr;
		}
	}

	HaVec_double AtmChg(nc);
	for(i = 0; i < nc; i++)
		AtmChg[i] = tmp[i];

    denbas_(&io_.iout,&icalc,&idgst,
		   &c0,&c0,&c0,&c0,
		   &c0,&c0,&c0,&xx,&xx,
		   &c0,&nbasis,&c1,&c0, 
		   coord.begin(),&c0, &c0,AtmChg.begin(), &natoms,
		   &xx,&c1,
		   &xx,valmat.begin(),&xx,
		   &c0,&c1,&c0,&xx,
	           &iprint,&ipflag,&allowp[0],&c0,
                   piwork,pwork,&mdv);


	iprint = 3;

	HaMat_double UdB;
	UdB.set_ext_alloc(valmat.begin(),ntt,3);

#if 0	
	cout << "Print Udb matricies: " << endl;
	for(i=1; i <=3 ; i++)
	{
		cout << "UdB component " << i << endl;
		HaVec_double fmat_lin(ntt,&UdB(1,i));
		HaMat_double scr;
		mat_asymm_from_lin(scr,fmat_lin,nbasis);
		scr.Print_format(cout," %10.5f");;
	}
#endif
	mat_scale(UdB, UdB, 0.5);


	mat_diff(fangm,fangm,UdB); // substract derivatives of Potential energy versus B
                               // as in l10002 module

	int icntrl=6127; int iopcl=0; int iraf =0;
	logical c_true=1, c_false=0;
    double zero= 0.0, one = 1.0;	
	int nmat=1, nsing=0, ntripl = 0;
	logical DoPurF=c_false;
	double ScaHFX=one;
	double* pdena = NULL;
	double* pdenb=  NULL;

	int iout = 6;
	int ihmeth = 1;
	int c__0=0, c__1=1;	
	int iaxx;

	fofdir_( &iout, &iprint,  &ihmeth, 
		&iopcl, &icntrl, 
		&iraf, &ipflag, &allowp[0], 
		&c0, &c0, 
		&c_false, &c_true, &DoPurF, 
		&zero, &ScaHFX, 
		&nmat, &nsing, &ntripl, 
		&nbasis, 
		&c0, &c1, &c1,       // symmetry 
		&xx, &xx, &xx, &xx,        // is not used 
		&xx, &xx,                  // ha is not added
		pdena, pdenb, pdena, pdenb, 
		&xx, &xx, valmat.begin(), valmat.begin() , 
		&c__0, &c__0, &xx, &c_false, &c__0,       // coordinates are not used 
		&c_false, &xx, &xx,
		&c__1, &iaxx,
		&xx, &xx, &xx, 
		&mdv, pwork);

	free(pwork);

	HaMat_double f2eB;
	f2eB.set_ext_alloc(valmat.begin(),ntt,3); // 2e contributions to rX Grad operator fokian on London orbitals
	
	mat_scale(f2eB, f2eB, 0.25); // scale by 1/4 not by 1/2 as in 1002 because
	                             // dena is probably 2 times as it should be
	                             // thus I got agreement with l1002 

	mat_add(fangm,fangm,f2eB);
	mat_scale(fangm,fangm,-1.0);

	for ( i=0; i < 3; i++)
	{
	        size_t nb= nbasis;
		data[i].newsize(nbasis,nbasis);
		HaVec_double fmat_lin(ntt,&fangm(1,i+1)); 
		mat_asymm_from_lin(data[i],fmat_lin,nb);
	}

#endif

	return true;
}


bool 
HaOperRDelt::LondonDaltonCalc()
{
	int nb=ptr_qc_mod->GetNBfunc();
	int no=ptr_qc_mod->GetNumOccMO();
	int nv=ptr_qc_mod->GetNumVacMO();
	int nb2;

	ifstream mofile("MO.DAT");
	mofile >> nb2;
	assert( nb2 == nb*nb);

	int i;
	HaVec_double buf(nb2);
	for(i=1; i <= nb2; i++)
		mofile >> buf(i);

	mofile.close();

	HaMat_double cmo_dalt;
	cmo_dalt.set_ext_alloc(&buf(1),nb,nb);
	HaVec_double dalt_mo_sign(nb);

	HaMat_double& ss = ptr_qc_mod->GetOvlpMat();
	
	HaMat_double cmo= ptr_qc_mod->MO_coef;

	for(i=1; i <= nb; i++)
// fill array of sign correction for MOs from Dalton
	{
		dalt_mo_sign(i)= (cmo_dalt(1,i)*cmo(1,i) < 0) ? -1.0 : 1.0; 
	}

	int nsize;

	ifstream idfile("MAGLON.DAT");
	idfile >> nsize;
	assert( no*nv == nsize);

	HaMat_double SC_occ,SC_vac,scr ;
	HaMat_double scr_occ;
	scr_occ.set_ext_alloc(&cmo(1,1),nb,no);
	matmult(SC_occ, ss,scr_occ);

	HaMat_double scr_vac;
	scr_vac.set_ext_alloc(&cmo(1,no+1),nb,nv);
	matmult(SC_vac, ss,scr_vac);

	buf.newsize(nsize);
	HaMat_double Zm;
	Zm.set_ext_alloc(buf.begin(),nv,no);

	for(i=1; i <= nsize; i++)
		idfile >> buf(i);

	for(i=1; i <= no; i++)
	{
		for(int j=1; j <= nv; j++)
		{
			Zm(j,i)*=(-0.5)*dalt_mo_sign(i)*dalt_mo_sign(no+j);
		}
	}

	matmult(scr,SC_vac,Zm);
	matmult_T2(data[0],scr,SC_occ);

	for(i=1; i <= nsize; i++)
		idfile >> buf(i);


	for(i=1; i <= no; i++)
	{
		for(int j=1; j <= nv; j++)
		{
			Zm(j,i)*=(-0.5)*dalt_mo_sign(i)*dalt_mo_sign(no+j);
		}
	}

	matmult(scr,SC_vac,Zm);
	matmult_T2(data[1],scr,SC_occ);

	for(i=1; i <= nsize; i++)
		idfile >> buf(i);


	for(i=1; i <= no; i++)
	{
		for(int j=1; j <= nv; j++)
		{
			Zm(j,i)*=(-0.5)*dalt_mo_sign(i)*dalt_mo_sign(no+j);
		}
	}

	matmult(scr,SC_vac,Zm);
	matmult_T2(data[2],scr,SC_occ);
	
	return true;
}



bool
HaOperRDelt::RecalcFromHr()
{
	assert(ptr_qc_mod != NULL);
	assert(data.size() == 3);

	HaOperR r1;
	HaOperKinEner rT;
    HaMat_double sm1, tm;
	HaMat_doubleArr rm;
	HaMat_doubleArr gm;
	
	rT.EvalGauBasisSet(&(ptr_qc_mod->AtBasis),tm);
	sm1 = ptr_qc_mod->GetOvlpMat();
	HaMat_double::mat_inverse(sm1);

	r1.EvalGauBasisSet(&(ptr_qc_mod->AtBasis),rm);

	gm[0]= rm[0] * sm1 * tm - tm * sm1 * rm[0];
	gm[1]= rm[1] * sm1 * tm - tm * sm1 * rm[1];
	gm[2]= rm[2] * sm1 * tm - tm * sm1 * rm[2];

	data[0] = rm[1]*sm1*gm[2] - rm[2]*sm1*gm[1];
	data[1] = rm[2]*sm1*gm[0] - rm[0]*sm1*gm[2];
	data[2] = rm[0]*sm1*gm[1] - rm[1]*sm1*gm[0];

	data[0].ASymmetrize();
	data[1].ASymmetrize();
	data[2].ASymmetrize();


	return true;
}

bool
HaOperRDelt::RecalcFromHr2()
{
	assert(ptr_qc_mod != NULL);
	assert(data.size() == 3);

	HaMat_double& ss = ptr_qc_mod->GetOvlpMat();

	HaOperR r1;
	HaMat_double tmp;
	HaMat_doubleArr rmats;
    HaMat_double sm1 = ss;
	HaMat_double::mat_inverse(sm1);
	
	r1.FillMat(&(ptr_qc_mod->AtBasis),rmats);
	int i;
	for( i=0; i< 3 ; i++)
	{
		matmult(tmp,rmats[i],sm1);
		matmult(rmats[i],sm1,tmp);
	}
	return true;
	
	HaRPAHam h1;
	h1.SetEnergy(0.0);

	vector<HaRPAvec> x;
	
	HaMat_double OrbMat;
	int ic;
	for( ic=0; ic <3 ; ic++) 
	{
		OrbMat = rmats[ic];
		HaRPAvec RPAv(*ptr_qc_mod);
		RPAv.SetFromAOMat(OrbMat,REAL_OPER);
		x.push_back(RPAv);
	}
	h1.Apply(x);
	
    HaOperGrad rg1;

	HaMat_doubleArr gm;

	rg1.EvalGauBasisSet(&ptr_qc_mod->AtBasis,gm); 

	HaMat_double Hrmx, Hrmy, Hrmz;
	HaMat_double scr;
	
	int nb=ptr_qc_mod->GetNBfunc();
	int no=ptr_qc_mod->GetNumOccMO();
	int nv=ptr_qc_mod->GetNumVacMO();
	
	HaMat_double SC_occ;
	HaMat_double scr_occ;
	scr_occ.set_ext_alloc(ptr_qc_mod->MO_coef.begin(),nb,no);
	matmult(SC_occ, ss,scr_occ);

	HaMat_double SC_vac;
	HaMat_double scr_vac;
	scr_vac.set_ext_alloc(&(ptr_qc_mod->MO_coef(1,no+1)),nb,nv);
	matmult(SC_vac, ss,scr_vac);

	matmult(scr,SC_vac,x[0].Z_mat);
	matmult_T2(Hrmx,scr,SC_occ);

	matmult(scr,SC_vac,x[1].Z_mat);
	matmult_T2(Hrmy,scr,SC_occ);

	matmult(scr,SC_vac,x[2].Z_mat);
	matmult_T2(Hrmz,scr,SC_occ);

	int j;
	for(i=1; i <= nb; i++)
	{
		for(j=1; j <= nb ; j++)
		{
			const HaAtom* aprt1= ptr_qc_mod->GetAtomOfAO(i);
			const HaAtom* aprt2= ptr_qc_mod->GetAtomOfAO(j);
			double X= 0.5*(aprt1->GetX() + aprt2->GetX());
			double Y= 0.5*(aprt1->GetY() + aprt2->GetY());
			double Z= 0.5*(aprt1->GetZ() + aprt2->GetZ());

			data[0](i,j)-= Y*gm[2](i,j)- Z*gm[1](i,j);
			data[1](i,j)-= Z*gm[0](i,j)- X*gm[2](i,j);
			data[2](i,j)-= X*gm[1](i,j)- Y*gm[0](i,j);


			data[0](i,j)=0.0; 
			data[1](i,j)=0.0; 
			data[2](i,j)=0.0;


			data[0](i,j)-= Y*Hrmz(i,j)- Z*Hrmy(i,j);
			data[1](i,j)-= Z*Hrmx(i,j)- X*Hrmz(i,j);
			data[2](i,j)-= X*Hrmy(i,j)- Y*Hrmx(i,j);

		}
	}

	return true;
}


HaOperGrad::HaOperGrad()
{
	
}

HaOperGrad::~HaOperGrad(void)
{

}

int HaOperGrad::EvalGauBasisSet(GauBasisSet* pbset,HaMat_doubleArr& fmats)
{
	int nbasis = pbset->GetNBfunc();
	if(HaQCMod::int_engine == HaQCMod::int_engine.INT_ENGINE_GAUSS)
	{
		PrintLog("Compute Gradient Operator Matrix using GAUSSIAN Library \n");

#if defined(GAUSSVER) 

		int c0=0, c1=1;
		int iprint=0;
		double xx=0;
		HaMat_double coord;
		HaQCMod::LoadGauCom(*pbset);
		int nbasis = pbset->GetNBfunc();
		int nbas6d;
		getnb6_(&nbas6d);
		int ntt6d=nbas6d*(nbas6d+1)/2;

		int ntt= nbasis*(nbasis+1)/2;
		int natoms = pbset->GetNumCnt();
		pbset->GetCntCoord(coord);
		HaVec_double scr_val( ntt6d*6);
		HaMat_double valmat;
		valmat.set_ext_alloc(scr_val.begin(),ntt,6); 
		
		int icalc=12; // Calc Del anf R X Del
		int idgst=10; // return Operator Matricies

		double* pwork;
		int mdv=10*nbasis*nbasis;
		if(mdv < 4000000) mdv=4000000;
	
		pwork=(double*)malloc(mdv*sizeof(double));
		assert(pwork != NULL);
		int* piwork=(int*)pwork;

		int ipflag=0;
		logical allowp[50];
		int i;
		for(i=0; i< 50; i++)
			allowp[i]= 0;

		denbas_(&io_.iout,&icalc,&idgst,
				&c0,&c0,&c0,&c0,
				&c0,&c0,&c0,&xx,&xx,
				&c0,&nbasis,&c1,&c0,
				coord.begin(),&c0,&c0,&xx, &natoms,
				&xx,&c0,
				&xx,valmat.begin(),&xx,
				&c0,&c1,&c0,&xx,
				&iprint,&ipflag,&allowp[0],&c0,
                piwork,pwork,&mdv);

		free(pwork);

		if(fmats.size() != 3) fmats.resize(3);

		for (i=0; i < 3; i++)
		{
			size_t nb = nbasis;
			fmats[i].newsize(nbasis,nbasis);
			HaVec_double fmat_lin(ntt,&valmat(1,i+1)); 
			mat_asymm_from_lin(fmats[i],fmat_lin,nb);
		}
#else
		PrintLog("GAUSSIAN Library is not available \n");
#endif
	}
	else if(HaQCMod::int_engine == QCIntEngineType::INT_ENGINE_IPACK)
	{
		PrintLog("Compute Gradient Operator Matrix using IPACK Library \n");
		PrintLog("Incorrectly computed in IPACK!!! - substituted by R operator \n");
 
		HaQCMod::InitIPack();
		InternalBasis* ip_bas = pbset->CreateIPackBas();

		int nb = nbasis;
	    ARRAY<Mat> r_mat_arr(3);
		int i;
		for(i=0; i <3; i++)
		{
			r_mat_arr[0].reset(nb,nb);
			r_mat_arr[0].set(0);
		}
		int same_bas = TRUE;
		Location cnt(0.0,0.0,0.0);
		MomentData mom_data;
		mom_data.max_mu = 1;
		mom_data.loc.assign(cnt);
        
		moment(*ip_bas, mom_data, r_mat_arr);
		normalize(*ip_bas, *ip_bas,r_mat_arr,1);
		if(fmats.size() != 3) fmats.resize(3);
		for(i = 1; i < 4; i++)
		{
            fmats[i-1].set_from(r_mat_arr[i]);
		}
		free(ip_bas);
	}
	return true;
}



HaOperKinEner::HaOperKinEner()
{
	
}

HaOperKinEner::~HaOperKinEner(void)
{

}

int HaOperKinEner::EvalGauBasisSet(GauBasisSet* pbset,HaMat_double& fmat)
{
	if(HaQCMod::int_engine == HaQCMod::int_engine.INT_ENGINE_GAUSS)
	{
		PrintLog("Compute Kinetic Energy Operator T Matrix using GAUSSIAN Library \n");

#if defined(GAUSSVER) 

		int c0=0, c1=1;
		int iprint=0;
		double xx=0;
		HaMat_double coord;
		HaQCMod::LoadGauCom(*pbset);
		int nbasis = pbset->GetNBfunc();
		int nbas6d;
		getnb6_(&nbas6d);
		int ntt6d=nbas6d*(nbas6d+1)/2;

		int ntt= nbasis*(nbasis+1)/2;
		int natoms = pbset->GetNumCnt();
		pbset->GetCntCoord(coord);
		HaVec_double scr_val( ntt6d*2);
		HaMat_double valmat;
		valmat.set_ext_alloc(scr_val.begin(),ntt,2); 
		
		int icalc=1;  // Overlap an Kinetic Energy integrals
		int idgst=10; // return Operator Matricies

		double* pwork;
		int mdv=10*nbasis*nbasis;
		if(mdv < 4000000) mdv=4000000;
	
		pwork=(double*)malloc(mdv*sizeof(double));
		assert(pwork != NULL);
		int* piwork=(int*)pwork;

		int ipflag=0;
		logical allowp[50];
		int i;
		for(i=0; i< 50; i++)
			allowp[i]= 0;

		denbas_(&io_.iout,&icalc,&idgst,
		        &c0,&c0,&c0,&c0,
		        &c0,&c0,&c0,&xx,&xx,
	            &c0,&nbasis,&c1,&c0,
		        coord.begin(),&c0, &c0,&xx, &natoms,
		        &xx,&c0,
		        &xx,valmat.begin(),&xx,
		        &c0,&c1,&c0,&xx,
	            &iprint,&ipflag,&allowp[0],&c0,
                   piwork,pwork,&mdv);

	    free(pwork);
		
        size_t nb = nbasis;
		fmat.newsize(nbasis,nbasis);
		HaVec_double fmat_lin(ntt,&valmat(1,2)); 
		mat_symm_from_lin(fmat,fmat_lin,nb);
#else
		PrintLog("GAUSSIAN Library is not available \n");
#endif
	}
	else if( HaQCMod::int_engine == QCIntEngineType::INT_ENGINE_IPACK )
	{
		PrintLog("Compute Kinetic Energy Matrix using IPACK Library \n");
    	HaQCMod::InitIPack();
		InternalBasis* ip_bas = pbset->CreateIPackBas();
		
		int nb = pbset->GetNBfunc();

	    ARRAY<Mat> t_mat_arr(1);
        t_mat_arr[0].reset(nb,nb);
		t_mat_arr[0].set(0);

		kinetic(*ip_bas,t_mat_arr);
		normalize(*ip_bas,*ip_bas,t_mat_arr,0);
		fmat.set_from(t_mat_arr[0]);
		
		free(ip_bas);		
	}

	return true;
}

int HaOperR::FillMat(HaBasisSet* pbas1, HaMat_doubleArr& rmats)
{
	if(pbas1 == NULL) return FALSE;
    
	std::string bas_type = pbas1->GetClassName();
	if(bas_type == "GauBasisSet")
	{
		return EvalGauBasisSet((GauBasisSet*)pbas1,rmats);
	}
	if(bas_type == "LinCombOrb3D")
	{
        HaMat_doubleArr rmat_b;
		LinCombOrb3D* lcmb1 = (LinCombOrb3D*) pbas1;
        int ires = FillMat(lcmb1->bas,rmat_b);
		if(!ires || rmat_b.size() != 3) return FALSE;
		if(rmats.size() != 3) rmats.resize(3);
		int i;
		for(i = 0; i < 3; i++)
		{
			LinCombOrb3D::Eval1eOp(lcmb1, lcmb1, rmat_b[i], rmats[i]);
		}
		return TRUE;
	}
	return FALSE;
}

int
HaOperRDelt::FillMat(HaBasisSet* pbas1, HaMat_doubleArr& rmats)
{
	if(pbas1 == NULL) return FALSE;
	std::string bas_type_1 = pbas1->GetClassName();
	if(bas_type_1 == "GauBasisSet")
	{
		return EvalGauBasisSet((GauBasisSet*)pbas1,rmats);
	}
	if(bas_type_1 == "LinCombOrb3D")
	{
        HaMat_doubleArr rmat_b;
		LinCombOrb3D* lcmb1 = (LinCombOrb3D*)pbas1;
        int ires = FillMat(lcmb1->bas,rmat_b);
		if(!ires || rmat_b.size() != 3) return FALSE;
		if(rmats.size() != 3) rmats.resize(3);
		int i;
		for(i = 0; i < 3; i++)
		{
			LinCombOrb3D::Eval1eOp(lcmb1, lcmb1, rmat_b[i], rmats[i]);
		}
		return TRUE;
	}
	return FALSE;
}

int
HaOperGrad::FillMat(HaBasisSet* pbas1, HaMat_doubleArr& rmats)
{
	if(pbas1 == NULL) return FALSE;

	std::string bas_type_1 = pbas1->GetClassName();
	if(bas_type_1 == "GauBasisSet" )
	{
		return EvalGauBasisSet((GauBasisSet*)pbas1,rmats);
	}
	if(bas_type_1 == "LinCombOrb3D")
	{
        HaMat_doubleArr rmat_b;
		LinCombOrb3D* lcmb1 = (LinCombOrb3D*)pbas1;
        int ires = FillMat(lcmb1->bas,rmat_b);
		if(!ires || rmat_b.size() != 3) return FALSE;
		if(rmats.size() != 3) rmats.resize(3);
		int i;
		for(i = 0; i < 3; i++)
		{
			LinCombOrb3D::Eval1eOp(lcmb1, lcmb1, rmat_b[i], rmats[i]);
		}
		return TRUE;
	}
	return FALSE;
}


int HaOperR::EvalGauBasisSet(GauBasisSet* pbas1, HaMat_doubleArr& rmats)
{
	int nbasis = pbas1->GetNBfunc();
	if(HaQCMod::int_engine == HaQCMod::int_engine.INT_ENGINE_GAUSS)
	{
		PrintLog("Compute electric dipole Matrix using GAUSSIAN Library \n");

#if defined(GAUSSVER) 

		int c0=0, c1=1, c2=2;
		int iprint=0;
		double xx=0;
		HaMat_double coord;
		HaQCMod::LoadGauCom(*pbas1);

		int nbas6d;
		getnb6_(&nbas6d);
		int ntt6d=nbas6d*(nbas6d+1)/2;
		int ntt= nbasis*(nbasis+1)/2;
		int natoms = pbas1->GetNumCnt();
		pbas1->GetCntCoord(coord);
		HaVec_double scr_val( ntt6d*3);
		HaMat_double valmat;
		valmat.set_ext_alloc(scr_val.begin(),ntt,3); 
		
//	HaVec_double AtmChg(natoms,1.0);

		int icalc=7; // Multipole integrals
		int idgst=10; // return Operator Matricies
		int maxl=1; 
        int minmat=1; // these parameter can be used
		int maxmat=3; // to calculate only some of multipole integrals

		double* pwork;
		int mdv=10*nbasis*nbasis;
		if(mdv < HaQCMod::max_gauss_mem) mdv=HaQCMod::max_gauss_mem;

		pwork=(double*)malloc(mdv*sizeof(double));
		assert(pwork != NULL);
		int* piwork=(int*)pwork;

		int ipflag=0;
		logical allowp[50];
		int i;
		for(i=0; i< 50; i++)
			allowp[i]= 0;

		denbas_(&io_.iout,&icalc,&idgst,
			    &c0,&c0,&c1,&maxl,
				&minmat,&maxmat,&c0,&xx,&xx, 
			    &c0,&nbasis,&c1, &c0, 
			    coord.begin(), &c0,&c0,&xx, &natoms,
			    &xx,&c0,
			    &xx,valmat.begin(),&xx,
			    &c0,&c1,&c2,&xx,
		        &iprint,&ipflag,&allowp[0], &c0,
                piwork,pwork,&mdv);

	     free(pwork);

	     if(rmats.size() != 3) rmats.resize(3);

	     for (i=0; i < 3; i++)
		 {
	        size_t nb = nbasis;
		    rmats[i].newsize(nbasis,nbasis);
		    HaVec_double fmat_lin(ntt,&valmat(1,i+1)); 
		    mat_symm_from_lin(rmats[i],fmat_lin,nb);
		 }
#else
		PrintLog("GAUSSIAN Library is not available \n");
#endif
	}
	else if(HaQCMod::int_engine == QCIntEngineType::INT_ENGINE_IPACK)
	{
		PrintLog("Compute electric dipole Matrix using IPACK Library \n");
    	HaQCMod::InitIPack();
		InternalBasis* ip_bas = pbas1->CreateIPackBas();

		int nb = nbasis;
	    ARRAY<Mat> r_mat_arr(3);
		int i;
		for(i=0; i <3; i++)
		{
			r_mat_arr[0].reset(nb,nb);
			r_mat_arr[0].set(0);
		}
		int same_bas = TRUE;
		Location cnt(0.0,0.0,0.0);
		MomentData mom_data;
		mom_data.max_mu = 1;
		mom_data.loc.assign(cnt);
        
		moment(*ip_bas, mom_data, r_mat_arr);
		normalize(*ip_bas, *ip_bas,r_mat_arr,1);
		if(rmats.size() != 3) rmats.resize(3);
		for(i = 1; i < 4; i++)
		{
            rmats[i-1].set_from(r_mat_arr[i]);
		}
		free(ip_bas);
	}
	return true;
}



