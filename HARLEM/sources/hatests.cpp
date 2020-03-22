/*! file hatests.cpp

    Axxiliary classes to test functionality of HARLEM 

    \author Igor Kurnikov 
    \date 2006

*/

#define HATESTS_CPP
#define GAU_MAIN

#include "haconst.h"

#include <mpi.h>
#if !defined(HARLEM_PYTHON_NO)
#include "Python.h"
#endif
#include "wx/wfstream.h"

#include "haio.h"
#include "harlemapp.h"
#include "haatbasdb.h"
#include "hamolecule.h"
#include "hamolset.h"
#include "hamultipole.h"
#include "harpavec.h"
#include "harpaham.h"
#include "hamatdb.h"
#include "hadalton.h"
#include "hagaussian.h"
#include "haconst.h"
#include "command.h"
#include "abstree.h"
#include "haqchem.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "hamolmech.h"
#include "hapseudopot.h"
#include "math.h"
#include "argraph.h"
#include "randomip.h"
#include "argedit.h"
#include "vf2_sub_state.h"
#include "match.h"
#include "mm_force_field.h"
#include "hatests.h"




void HaTests::calc_polar_gcontr()
// Calculate Group-Group contributions to polarizability
{
	MolSet* pmset = GetCurMolSet();
	if(pmset == NULL) return;
	HaQCMod* ptr_qc_mod= pmset->GetQCMod(true);

	if(ptr_qc_mod == NULL) return;

  int nb=ptr_qc_mod->GetNBfunc();
  int nab=ptr_qc_mod->GetNActiveOrb();
  int ngrp=pmset->GetNChemGroups();

	assert(ngrp > 0 && ngrp < nb);

  HaOperR r1;

  HaMat_doubleArr rmats;
  HaMat_double sm1;
  HaMat_double tmp;
  HaMat_double xml,yml,zml;
  
  r1.FillMat(ptr_qc_mod->ActBas,rmats);
  HaBasisSet::CalcOvlpMat(ptr_qc_mod->ActBas,ptr_qc_mod->ActBas,sm1);
  HaMat_double::mat_inverse(sm1);

  int i;
  for(i=0 ; i < 3; i++)
  {
	matmult(tmp,sm1,rmats[i]);
	matmult(rmats[i],tmp,sm1);
  }

  xml = rmats[0];
  yml = rmats[1];
  zml = rmats[2];

  HaMat_double gsmat;

  vector<HaRPAvec> x;
  int ic,ig1,ig2;
  ic = 2;
  for( ic=1; ic <=3 ; ic++)
  {
	  for( ig1=1; ig1 <= ngrp; ig1++)
	  {
		  for( ig2=1; ig2 <= ngrp; ig2++)
		  {
			  HaRPAvec RPAv(*ptr_qc_mod);

			  std::string gid1=(pmset->GetChemGroupByIdx(ig1)).GetID();
			  std::string gid2=(pmset->GetChemGroupByIdx(ig2)).GetID();

 			  if(ic==1)ptr_qc_mod->ExtractLocOrbSubMat(gid1,gid2,xml,gsmat);
			  if(ic==2)ptr_qc_mod->ExtractLocOrbSubMat(gid1,gid2,yml,gsmat);
			  if(ic==3)ptr_qc_mod->ExtractLocOrbSubMat(gid1,gid2,zml,gsmat);
			  RPAv.SetFromLOGrpMat(gid1,gid2,gsmat);
			  x.push_back(RPAv);
		  }
	  }
  }


  HaRPAResolv g1;
  g1.SetEnergy(0.0773);

    g1.Apply(x);
//	g1.Apply_Init(x);

  HaMat_double rrt(3*ngrp*ngrp,3*ngrp*ngrp);

  int nn=0;
  for( ic=1; ic <=3 ; ic++)
  {
	  for( ig1=1; ig1 <= ngrp; ig1++)
	  {
		  for( ig2=1; ig2 <= ngrp; ig2++)
		  {
			  HaRPAvec RPAv(*ptr_qc_mod);
			  vector<HaRPAvec> RPAv_arr;
			  RPAv_arr.push_back(RPAv);

			  std::string gid1=(pmset->GetChemGroupByIdx(ig1)).GetID();
			  std::string gid2=(pmset->GetChemGroupByIdx(ig2)).GetID();
			
 			  if(ic==1)ptr_qc_mod->ExtractLocOrbSubMat(gid1,gid2,xml,gsmat);
			  if(ic==2)ptr_qc_mod->ExtractLocOrbSubMat(gid1,gid2,yml,gsmat);
			  if(ic==3)ptr_qc_mod->ExtractLocOrbSubMat(gid1,gid2,zml,gsmat);

			  RPAv_arr[0].SetFromLOGrpMat(gid1,gid2,gsmat);
			  nn++;
			  HaMat_double gg;
			  gg.set_ext_alloc(&rrt(1,nn),1,3*ngrp*ngrp);
			  gg=SProd(RPAv_arr,x);
		  }
	  }
  }
  mat_scale(rrt,rrt,-2.0);

  HaMat_double rr_tenz(3,3);
  rr_tenz=0.0;

  int ngsq=ngrp*ngrp;
  Index1D I1(1,ngsq);
  Index1D J1(1,ngsq);
	
  HaMat_double rrg(ngsq,ngsq);
  rrg(I1,I1)=rrt(I1,I1);
  cout << " Group contributions to xx component of Polariz Tensor " << endl;
  rrg.Print_format(cout," %12.6f ");

  int j;
  for( i=1; i <= ngsq; i++)
	  for( j=1; j <= ngsq; j++)
		rr_tenz(1,1)+=rrg(i,j);

  rrg(I1,I1)=rrt(I1+ngsq,I1+ngsq);
  cout << " Group contributions to yy component of Polariz Tensor " << endl;
  rrg.Print_format(cout," %12.6f ");

  for(i=1; i <= ngsq; i++)
	  for(j=1; j <= ngsq; j++)
		rr_tenz(2,2)+=rrg(i,j);


  rrg(I1,I1)=rrt(I1+2*ngsq,I1+2*ngsq);
  cout << " Group contributions to zz component of Polariz Tensor " << endl;
  rrg.Print_format(cout," %12.6f ");

  for(i=1; i <= ngsq; i++)
	  for(j=1; j <= ngsq; j++)
		rr_tenz(3,3)+=rrg(i,j);


  cout << " Total Polarization tenzor is: " << endl;
  cout << rr_tenz << endl;

}


void HaTests::save_grp_oper_mat()
{
//  Save to the database file Group-Group matricies of one-electron
//  operators for the active molecule
	MolSet* pmset = GetCurMolSet();
	if( pmset == NULL) return;
	HaQCMod* ptr_qc_mod= pmset->GetQCMod(true);
	if( ptr_qc_mod == NULL) return;

	int ngrp;
	HaMatDB f_goper_mat("group_oper_mat.hdb","w");

    HaOperR r1;
	HaMat_doubleArr rmats;
	HaMat_double sm1;
	HaMat_double tmp;
	r1.FillMat(ptr_qc_mod->ActBas,rmats);
	HaBasisSet::CalcOvlpMat(ptr_qc_mod->ActBas,ptr_qc_mod->ActBas,sm1);
    HaMat_double::mat_inverse(sm1);

	ngrp= pmset->GetNChemGroups();
	assert(ngrp > 0);

	HaMat_double ActOrbMat;
	HaMat_double ActOrbSubMat;
		
	int ic,ig1,ig2;
    for( ic=0 ; ic <3 ; ic++)
	{
	  matmult(tmp,sm1,rmats[ic]);
	  matmult(ActOrbMat,tmp,sm1);

          char buf[10];
	  for( ig1=1; ig1 <= ngrp; ig1++)
	  {
		  for( ig2=1; ig2 <= ngrp; ig2++)
		  {
			ChemGroup& g1= pmset->GetChemGroupByIdx(ig1);
			ChemGroup& g2= pmset->GetChemGroupByIdx(ig2);
	        ptr_qc_mod->ExtractLocOrbSubMat(g1.GetID(), g2.GetID(), ActOrbMat, ActOrbSubMat);
			std::string oper_str = "OPER_R_";
            sprintf(buf,"%1d",ic+1);
            oper_str += buf; 
			HaGrpOperID g_op_id(oper_str.c_str(),g1,g2);
			f_goper_mat.put(g_op_id,ActOrbSubMat);
		  }
	  }
  }

  HaOperRDelt rd1(ptr_qc_mod);
  rd1.FillMat(ptr_qc_mod->ActBas,rmats);
  HaBasisSet::CalcOvlpMat(ptr_qc_mod->ActBas,ptr_qc_mod->ActBas,sm1);
  HaMat_double::mat_inverse(sm1);
  
  int i;
  for(i=0 ; i < 3; i++)
  {
	matmult(tmp,sm1,rmats[i]);
	matmult(rmats[i],tmp,sm1);
  }
//	rd1.LondonDaltonCalc();
          char buf[10];

    for( ic=0; ic < 3 ; ic++)
	{
      ActOrbMat = rmats[i];
	  for( ig1=1; ig1 <= ngrp; ig1++)
	  {
		  for( ig2=1; ig2 <= ngrp; ig2++)
		  {
			ChemGroup& g1= pmset->GetChemGroupByIdx(ig1);
			ChemGroup& g2= pmset->GetChemGroupByIdx(ig2);
	        ptr_qc_mod->ExtractLocOrbSubMat(g1.GetID(), g2.GetID(), ActOrbMat, ActOrbSubMat);
			std::string oper_str = "OPER_RxDel_";
            sprintf(buf,"%1d",ic+1);
            oper_str += buf; 
			HaGrpOperID g_op_id(oper_str.c_str(),g1,g2);
			f_goper_mat.put(g_op_id,ActOrbSubMat);
		  }
	  }
  }

  f_goper_mat.close();
}

void HaTests::calc_polar_contr_f()
{
// Calculate Group-Group contributions to polarizability
// using Operator active orbital matricies saved in the file

	MolSet* pmset = GetCurMolSet();
	if( pmset == NULL) return;
	HaQCMod* ptr_qc_mod= pmset->GetQCMod(true);

	if(ptr_qc_mod == NULL) return;

  int nb=ptr_qc_mod->GetNBfunc();
  int nab=ptr_qc_mod->GetNActiveOrb();
  int ngrp= pmset->GetNChemGroups();

  assert(ngrp > 0 && ngrp < nb);

  HaMat_double ActOrbMat(nab,nab,0.0);
  HaMat_double ActOrbSubMat;

  assert(ngrp > 0 && ngrp < nb);

  HaMatDB f_goper_mat("group_oper_mat.hdb","r");

  vector<HaRPAvec> x;
  int ic,ig1,ig2;

  int nn1[3]={0,0,0}; // number of RPA vectors of the certain projection

  char buf[10];
  for( ic=1; ic <=3 ; ic++)
  {
	  for( ig1=1; ig1 <= ngrp; ig1++)
	  {
		  for( ig2=1; ig2 <= ngrp; ig2++)
		  {
			  HaRPAvec RPAv(*ptr_qc_mod);
			
			  ChemGroup& g1= pmset->GetChemGroupByIdx(ig1);
			  ChemGroup& g2= pmset->GetChemGroupByIdx(ig2);
			  std::string oper_str = "OPER_R_";
              sprintf(buf,"%1d",ic+1);
              oper_str += buf; 
			  HaGrpOperID g_op_id(oper_str.c_str(),g1,g2);

			  if(!f_goper_mat.get(g_op_id,ActOrbSubMat))
			  {
					cout << "not found: " << g_op_id.GetStr() << endl;
					continue;
			  }

			  cout << "found: " << g_op_id.GetStr() << endl;
			  RPAv.SetFromLOGrpMat(g1.GetID(),g2.GetID(),ActOrbSubMat);
			  RPAv.id=g_op_id;
			  x.push_back(RPAv);
			  nn1[ic-1]++;
		  }
	  }
  }
  int nt=nn1[0]+nn1[1]+nn1[2];

  HaRPAResolv g1;
  g1.SetEnergy(0.0773);

  g1.Apply(x);

  HaMat_double rrt(nt,nt);
  int nn2[3]={0,0,0};

  HaMatDB polar_contr_db("polar_contr_db.hdb","a");

  int m=0;
  for( ic=1; ic <=3 ; ic++)
  {
	  for( ig1=1; ig1 <= ngrp; ig1++)
	  {
		  for( ig2=1; ig2 <= ngrp; ig2++)
		  {
			  HaRPAvec RPAv(*ptr_qc_mod);
			  vector<HaRPAvec> RPAv_arr;
			  RPAv_arr.push_back(RPAv);
			
			  ChemGroup& g1 = pmset->GetChemGroupByIdx(ig1);
			  ChemGroup& g2 = pmset->GetChemGroupByIdx(ig2);
			  std::string oper_str = "OPER_R_";
              sprintf(buf,"%1d",ic);
              oper_str += buf;
			  HaGrpOperID g_op_id(oper_str.c_str(),g1,g2);

			  if(!f_goper_mat.get(g_op_id,ActOrbSubMat))
			  {
				  cout << "not found: " << g_op_id.GetStr() << endl;
				  continue;
			  }
			
			  cout << "found: " << g_op_id.GetStr() << endl;
			  RPAv_arr[0].SetFromLOGrpMat(g1.GetID(),g2.GetID(),ActOrbSubMat);
			  m++;
			  HaMat_double gg;
			  gg.set_ext_alloc(&rrt(1,m),1,nt);
			  gg=SProd(RPAv_arr,x);

			  for(int i=1; i <=  x.size(); i++)
			  {
				  HaGrp4MatID g4_mat_id(g_op_id,x[i-1].id);
				  HaMat_double dm(1,1,gg(1,i));
				  polar_contr_db.put(g4_mat_id,dm);
			  }
		  }
	  }
  }
  polar_contr_db.close();

  mat_scale(rrt,rrt,-2.0);

  HaMat_double rr_tenz(3,3);
  rr_tenz=0.0;

  Index1D I1(1,nn1[0]);
  Index1D J1(1,nn1[0]);
  Index1D I2(1,nn1[1]);
  Index1D J2(1,nn1[1]);
  Index1D I3(1,nn1[2]);
  Index1D J3(1,nn1[2]);
	
  HaMat_double rrg(nn1[0],nn1[0]);
  rrg(I1,I1)=rrt(I1,I1);
  cout << " Group contributions to xx component of Polariz Tensor " << endl;
  rrg.Print_format(cout," %12.6f ");

  int i,j;
  for( i=1; i <= nn1[0]; i++)
	  for( j=1; j <= nn1[0]; j++)
		rr_tenz(1,1)+=rrg(i,j);

  rrg.newsize(nn1[1],nn1[1]);
  rrg(I2,I2)=rrt(I2+nn1[0],I2+nn1[0]);
  cout << " Group contributions to yy component of Polariz Tensor " << endl;
  rrg.Print_format(cout," %12.6f ");

  for(i=1; i <= nn1[1]; i++)
	  for(j=1; j <= nn1[1]; j++)
		rr_tenz(2,2)+=rrg(i,j);


  rrg.newsize(nn1[2],nn1[2]);
  rrg(I3,I3)=rrt(I3+nn1[0]+nn1[1],I3+nn1[0]+nn1[1]);
  cout << " Group contributions to zz component of Polariz Tensor " << endl;
  rrg.Print_format(cout," %12.6f ");

  for(i=1; i <= nn1[2]; i++)
	  for(j=1; j <= nn1[2]; j++)
		rr_tenz(3,3)+=rrg(i,j);


  cout << " Total Polarization tenzor is: " << endl;
  cout << rr_tenz << endl;

}


void
HaTests::calc_polar_contr_2idx()
{
// Calculate Group-Group contributions to polarizability
// using Operator active orbital matricies saved in the file
	MolSet* pmset = GetCurMolSet();
	if( pmset == NULL) return;
	HaQCMod* ptr_qc_mod= pmset->GetQCMod(true);
	if( ptr_qc_mod == NULL) return;

  int nb=ptr_qc_mod->GetNBfunc();
  int nab=ptr_qc_mod->GetNActiveOrb();
  int ngrp= pmset->GetNChemGroups();

  assert(ngrp > 0 && ngrp < nb);

  int iprint=3;

  HaMat_double ActOrbMat(nab,nab,0.0);
  HaMat_double ActOrbSubMat;

  assert(ngrp > 0 && ngrp < nb);

  HaMatDB f_goper_mat("group_oper_mat.hdb","r");

  vector<HaRPAvec> x;
  int ic,ig1,ig2;

  HaMat_doubleArr rmats;
  HaMat_double sm1;
  HaMat_double tmp;

  HaOperR r1;

  r1.FillMat(ptr_qc_mod->ActBas,rmats);
  HaBasisSet::CalcOvlpMat(ptr_qc_mod->ActBas,ptr_qc_mod->ActBas,sm1);
  HaMat_double::mat_inverse(sm1);

  int i;
  for(i=0 ; i < 3; i++)
  {
	matmult(tmp,sm1,rmats[i]);
	matmult(rmats[i],tmp,sm1);
  }
  
  char buf[10];

  for( ic=0; ic <=2 ; ic++)
  {
	  std::string oper_str = "OPER_R_";
      sprintf(buf,"%1d",ic+1);  
	  HaGrpOperID g_op_id(oper_str.c_str() );
	  ActOrbMat = rmats[ic];
	  HaRPAvec RPAv(*ptr_qc_mod);
	  RPAv.SetFromAOMat(ActOrbMat,REAL_OPER);
	  if(iprint > 2)
	  {
		 cout << "Z_mat of RPA vector of R component " << ic << endl;
		 RPAv.Z_mat.Print_format(cout," %10.5f ");
		 cout << endl;
	  }
	  RPAv.id=g_op_id;
	  x.push_back(RPAv);
  }

  HaRPAResolv g1;
  g1.SetEnergy(0.0773);

  g1.Apply(x);

  HaMat_double rr_tenz(3,3,0.0);
  HaMat_double gg;

  HaMatDB polar_contr_db("polar_contr_2idx_db.hdb","a");

  for( ic=1; ic <=3 ; ic++)
  {
	  for( ig1=1; ig1 <= ngrp; ig1++)
	  {
		  for( ig2=1; ig2 <= ngrp; ig2++)
		  {
			  HaRPAvec RPAv(*ptr_qc_mod);
			  vector<HaRPAvec> RPAv_arr;
			  RPAv_arr.push_back(RPAv);
			
			  ChemGroup& g1= pmset->GetChemGroupByIdx(ig1);
			  ChemGroup& g2= pmset->GetChemGroupByIdx(ig2);
			  std::string oper_str = "OPER_R_";
              sprintf(buf,"%1d",ic);
              oper_str += buf;
			  HaGrpOperID g_op_id(oper_str.c_str(),g1,g2);

			  if(!f_goper_mat.get(g_op_id,ActOrbSubMat))
			  {
				  cout << "not found: " << g_op_id.GetStr() << endl;
				  continue;
			  }
			
			  cout << "found: " << g_op_id.GetStr() << endl;
			  RPAv_arr[0].SetFromLOGrpMat(g1.GetID(),g2.GetID(),ActOrbSubMat);
			  gg=SProd(RPAv_arr,x);

			  for(int i=1; i <=  x.size(); i++)
			  {
				  HaGrp4MatID g4_mat_id(g_op_id,x[i-1].id);
				  rr_tenz(ic,i)+=gg(1,i);
				  HaMat_double dm(1,1,gg(1,i));
				  polar_contr_db.put(g4_mat_id,dm);
			  }
		  }
	  }
  }
  polar_contr_db.close();

  mat_scale(rr_tenz,rr_tenz,-2.0);

  cout << " Total Polarization tenzor is: " << endl;
  cout << rr_tenz << endl;

}

void HaTests::calc_beta_contr_2idx()
{
// Calculate Group-Group contributions to polarizability
// using Operator active orbital matricies saved in the file

	MolSet* pmset = GetCurMolSet();
	if( pmset == NULL) return;
	HaQCMod* ptr_qc_mod= pmset->GetQCMod(true);
	if(ptr_qc_mod == NULL) return; 

  int nb=ptr_qc_mod->GetNBfunc();
  int nab=ptr_qc_mod->GetNActiveOrb();
  int ngrp= pmset->GetNChemGroups();

  assert(ngrp > 0 && ngrp < nb);

  int iprint=3;

  HaMat_double ActOrbMat(nab,nab,0.0);
  HaMat_double ActOrbSubMat;

  assert(ngrp > 0 && ngrp < nb);

  HaMatDB f_goper_mat("group_oper_mat.hdb","r");

  vector<HaRPAvec> x;
  int ic,ig1,ig2;

  HaOperRDelt r1(ptr_qc_mod);

//  rd1.RecalcLondon(&ptr_qc_mod->AtBasis);
//  rd1.RecalcFromHr2();
  r1.LondonDaltonCalc();

  HaMat_doubleArr rmats(3);
  HaMat_doubleArr rmats_l(3);

  HaMat_double sm1;
  HaMat_double tmp;

  int i;
  rmats = r1.data;

  for(i = 0; i < 3; i++)
  {
	  ptr_qc_mod->ProjMatToActBas(rmats[i], rmats_l[i]);
  }
 
  HaBasisSet::CalcOvlpMat(ptr_qc_mod->ActBas,ptr_qc_mod->ActBas,sm1);
  HaMat_double::mat_inverse(sm1);

  for(i=0 ; i < 3; i++)
  {
	 matmult(tmp,sm1,rmats_l[i]);
	 matmult(rmats_l[i],tmp,sm1);
  }

  char buf[10];
  for( ic=0; ic <3 ; ic++)
  {
	  std::string oper_str = "OPER_RxDel_";
          sprintf(buf,"%1d",ic);
          oper_str += buf;
	  HaGrpOperID g_op_id(oper_str.c_str() );
	  ActOrbMat = rmats_l[i];
	  HaRPAvec RPAv(*ptr_qc_mod);
	  RPAv.SetFromAOMat(ActOrbMat,IMAG_OPER);
	  if(iprint > 2)
	  {
		 cout << "Z_mat of RPA vector of RxGrad component " << ic << endl;
		 RPAv.Z_mat.Print_format(cout," %10.5f ");
		 cout << endl;
	  }
	  RPAv.id=g_op_id;
	  x.push_back(RPAv);
  }

  HaRPAResolv g1;
  g1.SetEnergy(0.0773);

  g1.Apply(x);

  if(iprint > 2)
  {
	  for(ic =0; ic< x.size(); ic++)
	  {
	     cout << "Z_mat of RPA vector of Solution (E-H)^-1|RxGrad) component " << ic << endl;
		 x[ic].Z_mat.Print_format(cout," %10.5f ");
		 cout << endl;
	  }
  }


  HaMat_double rr_tenz(3,3,0.0);
  HaMat_double gg;

  HaMatDB beta_contr_db("beta_contr_2idx_db.hdb","a");

  for( ic=1; ic <=3 ; ic++)
  {
	  for( ig1=1; ig1 <= ngrp; ig1++)
	  {
		  for( ig2=1; ig2 <= ngrp; ig2++)
		  {
			  HaRPAvec RPAv(*ptr_qc_mod);
			  vector<HaRPAvec> RPAv_arr;
			  RPAv_arr.push_back(RPAv);
			
			  ChemGroup& g1= pmset->GetChemGroupByIdx(ig1);
			  ChemGroup& g2= pmset->GetChemGroupByIdx(ig2);
			  std::string oper_str = "OPER_R_";
                          sprintf(buf,"%1d",ic);
                          oper_str += buf;
			  HaGrpOperID g_op_id(oper_str.c_str(),g1,g2);

			  if(!f_goper_mat.get(g_op_id,ActOrbSubMat))
			  {
				  cout << "not found: " << g_op_id.GetStr() << endl;
				  continue;
			  }
			
			  cout << "found: " << g_op_id.GetStr() << endl;
			  RPAv_arr[0].SetFromLOGrpMat(g1.GetID(),g2.GetID(),ActOrbSubMat,REAL_OPER);
			
			  if(iprint > 2)
			  {
				  cout << " Z_mat for R group pertubation comp g1 g2 " <<
					  ic << " " <<  g1.GetID() << " " << g2.GetID() << endl;
					RPAv_arr[0].Z_mat.Print_format(cout," %10.5f ");

			  }
			
			
			  gg=SProd(RPAv_arr,x);

			  for(int i=1; i <=  x.size(); i++)
			  {
				  HaGrp4MatID g4_mat_id(g_op_id,x[i-1].id);
				  rr_tenz(ic,i)+=gg(1,i);
				  HaMat_double dm(1,1,gg(1,i));
				  beta_contr_db.put(g4_mat_id,dm);
			  }
		  }
	  }
  }
  beta_contr_db.close();

  mat_scale(rr_tenz,rr_tenz,-1.0);

  cout << " Total G' tenzor is: " << endl;
  cout << rr_tenz << endl;

  double beta;
  beta=-(rr_tenz(1,1)+rr_tenz(2,2)+rr_tenz(3,3))/(3*0.0773);
  cout << " Beta is: " << endl;
  cout << beta << endl;

}


void HaTests::read_polar_contr()
{
  ofstream flog("harlem.log");

  MolSet* pmset = GetCurMolSet();
  if( pmset == NULL) return;
  HaQCMod* ptr_qc_mod= pmset->GetQCMod(true);
  if( ptr_qc_mod == NULL) return;

  int nb=ptr_qc_mod->GetNBfunc();
  int nab=ptr_qc_mod->GetNActiveOrb();
  int ngrp= pmset->GetNChemGroups();

  assert(ngrp > 0 && ngrp < nb);

  double beta=0.0;

  int ngsq=ngrp*ngrp;
 	
  HaMat_double rmt(ngsq,ngsq);

  HaMatDB polar_contr_db("polar_contr_db.hdb","r");

  HaMat_double gg;

  vector<HaMat_double> rmg_vec;

  char buf[10];
  for( int ic=1; ic <=3 ; ic++)
  {
	  std::string oper_str = "OPER_R_";
          sprintf(buf,"%1d",ic);
          oper_str += buf;
	  HaMat_double rmg(ngsq,ngsq);
	  for( int ig1=1; ig1 <= ngrp; ig1++)
	  {
		  for( int ig2=1; ig2 <= ngrp; ig2++)
		  {
			  ChemGroup& g1= pmset->GetChemGroupByIdx(ig1);
			  ChemGroup& g2= pmset->GetChemGroupByIdx(ig2);
			  HaGrpOperID g_op_id1(oper_str.c_str(),g1,g2);
			
			  for( int ig3=1; ig3 <= ngrp; ig3++)
			  {
				  for( int ig4=1; ig4 <= ngrp; ig4++)
				  {
					
					  ChemGroup& g3= pmset->GetChemGroupByIdx(ig3);
					  ChemGroup& g4= pmset->GetChemGroupByIdx(ig4);
					  HaGrpOperID g_op_id2(oper_str.c_str(),g3,g4);
					
					  HaGrp4MatID g4_mat_id(g_op_id1,g_op_id2);
					  if(!polar_contr_db.get(g4_mat_id,gg))
					  {
						  cout << "not found: " << g4_mat_id.GetStr() << endl;
						  flog << "not found: " << g4_mat_id.GetStr() << endl;
						  continue;
					  }
					  cout << "found: " << g4_mat_id.GetStr() << endl;
					  flog << "found: " << g4_mat_id.GetStr() << endl;
					  rmg(ig1+(ig2-1)*ngrp, ig3+(ig4-1)*ngrp)=gg(1,1);
				  }
			  }
		  }
	  }
	  mat_scale(rmg,rmg,-2.0);
	  rmg_vec.push_back(rmg);
  }

  cout << " Group contributions to xx component of Polariz Tensor " << endl;
  rmg_vec[0].Print_format(cout," %12.6f ");

  HaMat_double rr_tenz(3,3,0.0);

  int i,j;
  for( i=1; i <= ngsq; i++)
	  for( j=1; j <= ngsq; j++)
		rr_tenz(1,1)+=rmg_vec[0](i,j);

  cout << " Group contributions to yy component of Polariz Tensor " << endl;
  rmg_vec[1].Print_format(cout," %12.6f ");

  for( i=1; i <= ngsq; i++)
	  for( j=1; j <= ngsq; j++)
		rr_tenz(2,2)+=rmg_vec[1](i,j);

  cout << " Group contributions to xx component of Polariz Tensor " << endl;
  rmg_vec[2].Print_format(cout," %12.6f ");

  for( i=1; i <= ngsq; i++)
	  for( j=1; j <= ngsq; j++)
		rr_tenz(3,3)+=rmg_vec[2](i,j);


  cout << " Total Polarization tenzor is: " << endl;
  cout << rr_tenz << endl;

  flog.close();

}

void HaTests::read_polar_contr_2idx()
{
  ofstream flog("harlem.log");

  MolSet* pmset = GetCurMolSet();
  if( pmset == NULL) return;
  HaQCMod* ptr_qc_mod= pmset->GetQCMod(true);
  if( ptr_qc_mod == NULL) return;

  int nb=ptr_qc_mod->GetNBfunc();
  int nab=ptr_qc_mod->GetNActiveOrb();
  int ngrp=pmset->GetNChemGroups();

  assert(ngrp > 0 && ngrp < nb);

  HaMatDB polar_contr_db("polar_contr_2idx_db.hdb","r");

  HaMat_double rr_tenz(3,3,0.0);
  HaMat_double gg;

  char buf[10];

  for( int ic=1; ic <=3 ; ic++)
  {
	  for( int ig1=1; ig1 <= ngrp; ig1++)
	  {
		  for( int ig2=1; ig2 <= ngrp; ig2++)
		  {
			
			  ChemGroup& g1= pmset->GetChemGroupByIdx(ig1);
			  ChemGroup& g2= pmset->GetChemGroupByIdx(ig2);
			  std::string oper_str = "OPER_R_";
                          sprintf(buf,"%1d",ic);
                          oper_str += buf;
			  HaGrpOperID g_op_id1(oper_str.c_str(),g1,g2);
			
			  for( int ic2=1; ic2 <=3 ; ic2++)
			  {
			          oper_str = "OPER_R_";
                                  sprintf(buf,"%1d",ic2);	
                                  oper_str += buf;
				  HaGrpOperID g_op_id2(oper_str.c_str());
				
				  HaGrp4MatID g4_mat_id(g_op_id1,g_op_id2);
				  if(!polar_contr_db.get(g4_mat_id,gg))
				  {
					  cout << "not found: " << g4_mat_id.GetStr() << endl;
					  flog << "not found: " << g4_mat_id.GetStr() << endl;
					  continue;
				  }
				  cout << "found: " << g4_mat_id.GetStr() << endl;
				  flog << "found: " << g4_mat_id.GetStr() << endl;
				  rr_tenz(ic,ic2)+=gg(1,1);
			  }
		  }
	  }
  }

  mat_scale(rr_tenz,rr_tenz,-2.0);

  cout << " Total Polarization tenzor is: " << endl;
  cout << rr_tenz << endl;

  flog.close();

}

void HaTests::read_beta_contr_2idx()
{
  ofstream flog("harlem.log");

  MolSet* pmset = GetCurMolSet();
  if( pmset == NULL) return;
  HaQCMod* ptr_qc_mod= pmset->GetQCMod(true);
  if( ptr_qc_mod == NULL) return;

  int nb=ptr_qc_mod->GetNBfunc();
  int nab=ptr_qc_mod->GetNActiveOrb();
  int ngrp= pmset->GetNChemGroups();

  assert(ngrp > 0 && ngrp < nb);

  HaMatDB beta_contr_db("beta_contr_2idx_db.hdb","r");

  HaMat_double rr_tenz(3,3,0.0);
  HaMat_double gg;

  char buf[10];
  for( int ic=1; ic <=3 ; ic++)
  {
	  for( int ig1=1; ig1 <= ngrp; ig1++)
	  {
		  for( int ig2=1; ig2 <= ngrp; ig2++)
		  {
			
			  ChemGroup& g1= pmset->GetChemGroupByIdx(ig1);
			  ChemGroup& g2= pmset->GetChemGroupByIdx(ig2);
			  std::string oper_str = "OPER_R_";
              sprintf(buf,"%1d",ic);
              oper_str += buf;
			  HaGrpOperID g_op_id1(oper_str.c_str(),g1,g2);
			
			  for( int ic2=1; ic2 <=3 ; ic2++)
			  {
			          oper_str = "OPER_RxDel_";
                                  sprintf(buf,"%1d",ic2);	
                                  oper_str += buf;
				  HaGrpOperID g_op_id2(oper_str.c_str());
				
				  HaGrp4MatID g4_mat_id(g_op_id1,g_op_id2);
				  if(!beta_contr_db.get(g4_mat_id,gg))
				  {
					  cout << "not found: " << g4_mat_id.GetStr() << endl;
					  flog << "not found: " << g4_mat_id.GetStr() << endl;
					  continue;
				  }
				  cout << "found: " << g4_mat_id.GetStr() << endl;
				  flog << "found: " << g4_mat_id.GetStr() << endl;
				  rr_tenz(ic,ic2)+=gg(1,1);
			  }
		  }
	  }
  }

  mat_scale(rr_tenz,rr_tenz,-1.0);

  cout << " Total G' tenzor is: " << endl;
  cout << rr_tenz << endl;

  double beta;
  beta=-(rr_tenz(1,1)+rr_tenz(2,2)+rr_tenz(3,3))/(3*0.0773);
  cout << " Beta is: " << endl;
  cout << beta << endl;

  flog.close();

}

void HaTests::test_oper_1()
// function to check Analytical calculation of one-electron operators:
{
	MolSet* pmset = GetCurMolSet();
	if(!pmset) return ;

	HaQCMod* ptr_qc_mod= pmset->GetQCMod(true);

	if(ptr_qc_mod == NULL) return;

	HaMat_double& ss = ptr_qc_mod->GetOvlpMat();

	ss.Print_format(cout," %6.3f ");

	HaMat_double sm1 = ss;
	HaMat_double::mat_inverse(sm1);

	HaOperR r1;
	HaOperKinEner t1;
    HaOperGrad rg1;
	HaOperRDelt rd1(ptr_qc_mod);
	
//	rd1.RecalcLondon(&ptr_qc_mod->AtBasis);
//	cout << " r X Grad for london orbitals";
//	rd1.Print_info(cout,1);	
//	cout << endl;

	HaMat_doubleArr rm;
	HaMat_doubleArr gm;
	HaMat_doubleArr rdm;
	HaMat_doubleArr rdm_n;
	HaMat_double tmat;
	
	r1.FillMat(&ptr_qc_mod->AtBasis,rm);
    
	int ns = rm.size();
	int i;
	for(i = 0; i < ns; i++)
	{
		PrintLog(" \n");
		PrintLog(" The matrix of component %d of R-operator \n ",i+1);
		rm[i].Print_format(cout," %6.3f ");
	}

	t1.EvalGauBasisSet(&ptr_qc_mod->AtBasis,tmat);
	PrintLog(" \n");
	PrintLog(" The matrix of Kinietic Energy operator T \n");
	tmat.Print_format(cout," %6.3f ");
	
	rd1.EvalGauBasisSet(&ptr_qc_mod->AtBasis,rdm);
	ns = rdm.size();
	for(i = 0; i < ns; i++)
	{
		PrintLog(" \n");
		PrintLog(" The matrix of component %d of R X Delt operator \n ",i+1);
		rdm[i].Print_format(cout," %6.3f ");
	}

	return;

    rg1.EvalGauBasisSet(&ptr_qc_mod->AtBasis,gm);
	rd1.RecalcFromHr2();
	rdm = rd1.data;
		
	int nb=ptr_qc_mod->AtBasis.GetNBfunc();
	HaMat_double CMO;
	CMO.set_ext_alloc(ptr_qc_mod->MO_coef.begin(),nb,nb);

	cout << " MOs : " << endl;
	CMO.Print_format(cout," %10.5f ");


    HaMat_double scr, rdmx_mo;	
	matmult_T1(scr,CMO,rdm[0]);
	matmult   (rdmx_mo,scr,CMO);

	cout << "RxGrad_X in MO basis" << endl;
	rdmx_mo.Print_format(cout,"%10.5f ");
	cout << endl;

	cout << " RxGrad operator matricies with substracted ( 0.5*(R_1+R_2)*Grad_ij) " << endl;

	cout << "RxGrad_X corrected: " << endl;
	rdm[0].Print_format(cout,"%10.5f ");
	cout << "RxGrad_Y corrected: " << endl;
	rdm[1].Print_format(cout,"%10.5f ");
	cout << "RxGrad_Z corrected: " << endl;
	rdm[2].Print_format(cout,"%10.5f ");



	rdm_n[0] = rm[1]*sm1*gm[2] - rm[2]*sm1*gm[1];
	rdm_n[1] = rm[2]*sm1*gm[0] - rm[0]*sm1*gm[2];
	rdm_n[2] = rm[0]*sm1*gm[1] - rm[1]*sm1*gm[0];

	cout << " (R X Grad)_X mat calculated from R and Delt : " << endl;
	cout << rdm_n[0].Print_format(cout,"%10.5f ");
	cout << endl;

	cout << " (R X Grad)_Y mat calculated from R and Delt : " << endl;
	cout << rdm_n[1].Print_format(cout,"%10.5f ");
	cout << endl;

	cout << " (R X Grad)_Z mat calculated from R and Delt : " << endl;
	cout << rdm_n[2].Print_format(cout,"%10.5f ");
	cout << endl;

	HaOperKinEner rT;

	HaMat_double tm, tm_n;
	rT.EvalGauBasisSet(&ptr_qc_mod->AtBasis,tm);
	tm_n= gm[0]*sm1*gm[0]+gm[1]*sm1*gm[1]+gm[2]*sm1*gm[2]; 
	mat_scale(tm_n, tm_n, -0.5);

	cout << " Kinetic Energy mat calculated from Grad mat : " << endl;
	cout << tm_n.Print_format(cout,"%10.5f ");
	cout << endl;

    HaMat_doubleArr gm_n;

	gm_n[0]= rm[0] * sm1 * tm - tm * sm1 * rm[0];
	gm_n[1]= rm[1] * sm1 * tm - tm * sm1 * rm[1];
	gm_n[2]= rm[2] * sm1 * tm - tm * sm1 * rm[2];
	
	cout << " Grad_X mat calculated as [T,r] : " << endl;
	cout << gm_n[0].Print_format(cout,"%10.5f ");
	cout << endl;

	cout << " Grad_Y mat calculated as [T,r] : " << endl;
	cout << gm_n[1].Print_format(cout,"%10.5f ");
	cout << endl;

	cout << " Grad_Z mat calculated as [T,r] : " << endl;
	cout << gm_n[2].Print_format(cout,"%10.5f ");
	cout << endl;

}

void HaTests::test_oper_2()
{
	MolSet* pmset = GetCurMolSet();
	if(pmset == NULL) return;
	HaQCMod* ptr_qc_mod= pmset->GetQCMod(true);

	assert(ptr_qc_mod != NULL);

	HaMat_double sm1;
	HaMat_double tmp;

	HaBasisSet::CalcOvlpMat(&ptr_qc_mod->AtBasis,&ptr_qc_mod->AtBasis,sm1);
	HaMat_double::mat_inverse(sm1);

	HaOperR r1;
	HaMat_doubleArr rmats;
	r1.FillMat(&ptr_qc_mod->AtBasis,rmats);
	int i;
	for(i = 0; i <3; i++)
	{
		matmult(tmp,rmats[i],sm1);
		matmult(rmats[i],sm1,tmp);		
	}

	HaRPAHam h1;
	h1.SetEnergy(0.0);

	vector<HaRPAvec> x;
	
	HaMat_double OrbMat;
	for( int ic=0; ic < 3 ; ic++)
	{
		OrbMat = rmats[i];
		HaRPAvec RPAv(*ptr_qc_mod);
		RPAv.SetFromAOMat(OrbMat,REAL_OPER);
		x.push_back(RPAv);
	}
	h1.Apply(x);
	
    HaOperGrad rg1;

	HaMat_doubleArr gm;

	rg1.FillMat(&ptr_qc_mod->AtBasis,gm);

	HaMat_double Hrmx, Hrmy, Hrmz;
	HaMat_double scr;
	
	int nb= ptr_qc_mod->GetNBfunc();
	int no= ptr_qc_mod->GetNumOccMO();
	int nv= ptr_qc_mod->GetNumVacMO();

    HaMat_double& ss = ptr_qc_mod->GetOvlpMat();	
	
	HaMat_double SC_occ;
	HaMat_double scr_occ;
	scr_occ.set_ext_alloc(ptr_qc_mod->MO_coef.begin(),nb,no);
	matmult(SC_occ, ss,scr_occ);

	HaMat_double SC_vac;
	HaMat_double scr_vac;
	scr_vac.set_ext_alloc(ptr_qc_mod->MO_coef.begin()+no*nb,nb,nv);
	matmult(SC_vac, ss,scr_vac);

	matmult(scr,SC_vac,x[0].Z_mat);
	matmult_T2(Hrmx,scr,SC_occ);

	matmult(scr,SC_vac,x[1].Z_mat);
	matmult_T2(Hrmy,scr,SC_occ);

	matmult(scr,SC_vac,x[2].Z_mat);
	matmult_T2(Hrmz,scr,SC_occ);

	cout << "x[0].Z_mat" << endl;
	x[0].Z_mat.Print_format(cout,"%10.5f ");
	cout << "x[1].Z_mat" << endl;
	x[1].Z_mat.Print_format(cout,"%10.5f ");
	cout << "x[2].Z_mat" << endl;
	x[2].Z_mat.Print_format(cout,"%10.5f ");

	HaMat_double CMO;
	CMO.set_ext_alloc(ptr_qc_mod->MO_coef.begin(),nb,nb);

    HaMat_double hrmx_mo, hrmy_mo,hrmz_mo ;	
	matmult_T1(scr,CMO,Hrmx);
	matmult   (hrmx_mo,scr,CMO);
	matmult_T1(scr,CMO,Hrmy);
	matmult   (hrmy_mo,scr,CMO);
	matmult_T1(scr,CMO,Hrmz);
	matmult   (hrmz_mo,scr,CMO);

	cout << "[H,r]_X in MO basis" << endl;
	hrmx_mo.Print_format(cout,"%10.5f ");
	cout << endl;

	cout << "[H,r]_Y in MO basis" << endl;
	hrmy_mo.Print_format(cout,"%10.5f ");
	cout << endl;

	cout << "[H,r]_Z in MO basis" << endl;
	hrmz_mo.Print_format(cout,"%10.5f ");
	cout << endl;

}

void HaTests::test_qcmod_1()
{
	MolSet* pmset = GetCurMolSet();
	if(pmset == NULL) return;
	HaQCMod* ptr_qc_mod = pmset->GetQCMod(true);
	const HaPseudoPot* pot_ptr;

	GauAtomBasis* bptr;

	bptr=bas_db.Extract("HAY1_DZ","Fe");	
	if(bptr != NULL)bptr->SaveGaussianInp(std::cout);

	bptr=bas_db.Extract("HAY1_DZ","Ru");	
	if(bptr != NULL)bptr->SaveGaussianInp(std::cout);

	bptr=bas_db.Extract("HAY1_DZ","C");	
	if(bptr != NULL)bptr->SaveGaussianInp(std::cout);
			
	std::cout << " Pseudo potential data " << std::endl << std::endl;

	pot_ptr = pseudo_db.Extract("HAY_1","Fe");

	if(pot_ptr != NULL)
	{
		pot_ptr->SaveGaussInp(std::cout);
		std::cout << "********" << std::endl;
	}
	else
	{
		cout << "test_qcmod_1(): didn't find Fe in PseudoPot DB " << endl;
	}

	pot_ptr = pseudo_db.Extract("HAY_1","Ru");
	if(pot_ptr != NULL)
	{
		pot_ptr->SaveGaussInp(std::cout);
		std::cout << "********" << std::endl;
	}
	else
	{
		cout << "test_qcmod_1(): didn't find Ru in PseudoPot DB " << endl;
	}
}


void HaTests::dump_mol_info()
{
	MolSet* pmset = GetCurMolSet();
	if(pmset == NULL) return;
	pmset->Print_info(cout,1);
}

void HaTests::dump_gauss_bcommon()
{
	HaGaussMod::PrintCurBcommon();
}

void HaTests::dump_overlap()
{
	MolSet* pmset = GetCurMolSet();
	if( pmset == NULL) return;
	HaQCMod* ptr_qc_mod= pmset->GetQCMod(true);

	if(ptr_qc_mod == NULL)
	{
		cerr << " Error in dump_overlap() " << endl;
		cerr << " QChem Module is not set " << endl;
		return ;
	}

	HaMat_double ss;

	GauBasisSet::CalcOvlpMat( &ptr_qc_mod->AtBasis, &ptr_qc_mod->AtBasis,ss);

	PrintLog(" Overlap matrix \n"); 
	ss.Print_format(cout,"%10.5f ");
}

void HaTests::dump_overlap2()
{
	HaMat_double ss;
	GauBasisSet bas1,bas2;

	MolSet* pmset = GetCurMolSet();
	if( pmset == NULL) return;

//	HaQCMod::int_engine = QCIntEngineType::INT_ENGINE_GAUSS;
    HaQCMod::int_engine = QCIntEngineType::INT_ENGINE_IPACK;

	bas1.InitForMolSet("3-21G", pmset);
//	bas2.InitForMolSet("MINB6G",pmset);
//	bas1.InitForMolSet("sto-3g", pmset);
//	bas1.InitForMolSet("sto-6g", pmset);
	bas2.InitForMolSet("sto-6g", pmset);

//	GauBasisSet::CalcOvlpMat( &bas1, &bas2,ss);
	GauBasisSet::CalcOvlpMat( &bas1, &bas2,ss);

	PrintLog(" Overlap matrix between 3-21G and STO-6G basis \n"); 
//  PrintLog(" Overlap matrix between STO-3G and STO-6G basis \n"); 
	ss.Print_format(cout,"%10.5f ");
}

class HaMinQ1 : public HaMinimizer
{
public:
	HaMinQ1() { };
	~HaMinQ1() { };

	virtual int CalcValGrad(HaVec_double& x, double& val, HaVec_double& grad)
	{
		int ns = x.size();
		int i;
		val = 0.0;
		for(i=0; i < ns; i++)
		{
			val += i*(x[i] - i)*(x[i] - i)*(x[i] - i)*(x[i] - i);
			grad[i] = 4*i*(x[i] - i)*(x[i] - i)*(x[i] - i);
		}
		nitr++;

		if(nitr == 1)
		{
			glast.newsize(ns);
		}
	
		double df = 0.0;
		if( nitr > 1)
		{
			df = fabs( val - flast);
		}

		double gmax = 0.0;
		for(i=0; i < ns; i++)
		{
			gmax = MaxFun(fabs(grad[i]),gmax);
		}
	
		PrintLog(" iter= %5d  f= %12.6e  df = %12.6e  gmax = %12.6e \n",
			       nitr, val,df,gmax);
		
		flast = val;
		for(i=0; i < ns; i++)
			glast[i] = grad[i];

		return TRUE;
	}
};

void
HaTests::test_min_1()
//! function to test minimization class
{
	HaMinQ1 cmin;
	cmin.SetNVar(3);

	HaVec_double pt_ini(3,0.0);
	cmin.SetInitPoint(pt_ini);

	cmin.Minimize();
}

bool my_visitor(int n, node_id ni1[], node_id ni2[], void *usr_data)
{
	int* pnc = (int*) usr_data;
	(*pnc)++;

	int i;
	int nc = *pnc;

	PrintLog("\nMatch number %d \n",nc);
    for(i=0; i<n; i++)
	{
    PrintLog("\tNode %hd of graph 1 is paired with node %hd of graph 2\n",
               ni1[i], ni2[i]);
	}
	PrintLog("\n");
	return false;
}

void HaTests::test_python_1()
//! function to test HARLEM/PYTHON interactions
{
#if !defined(HARLEM_PYTHON_NO)
	PyRun_SimpleString("str1 = \"this is str one\"");
	PyObject* main_module = PyImport_AddModule("__main__");
	PyObject* global_dict = PyModule_GetDict(main_module);
	
    PyObject* str_old = PyDict_GetItemString(global_dict,"str1");

#if PY_VERSION_HEX >= 0x03000000
	PrintLog(" Here is the extracted str1 content: %s \n", PyBytes_AS_STRING(str_old));
#else
	PrintLog(" Here is the extracted str1 content: %s \n", PyString_AsString(str_old));
#endif
#endif
}


void HaTests::test_graph_1()
//! function to test graph matching library
{
    ARGEdit ed1,ed2;  // The object used to create the graph
    int i;

        // Insert the four nodes
        for(i=0; i<4; i++)
		{
          ed1.InsertNode(NULL); // The inserted node will have index i.
          ed2.InsertNode(NULL); // NULL stands for no semantic attribute.
		}

        // Insert the edges
		ed1.InsertEdge(0, 1, NULL);
		ed1.InsertEdge(0, 2, NULL);
		ed1.InsertEdge(0, 3, NULL);

		ed2.InsertEdge(1, 0, NULL);
		ed2.InsertEdge(1, 2, NULL);
		ed2.InsertEdge(1, 3, NULL);


        // Now the Graph can be constructed...
        Graph g1(&ed1);
		Graph g2(&ed2);

        VF2SubState s0(&g2, &g1);
		VF2SubState s1(&g2, &g1);

  int n;
  node_id ni1[100], ni2[100];

  if (!match(&s0, &n, ni1, ni2))
    { 
	  PrintLog("No matching found!\n");
      return;
    }

  PrintLog("Found a matching with %d nodes:\n", n);
  
  for(i=0; i<n; i++)
    PrintLog("\tNode %hd of graph 1 is paired with node %hd of graph 2\n",
               ni1[i], ni2[i]);

  int nc = 0;
  match(&s1,my_visitor, &nc);

}

static double fene1_old(double x)
{
	double amp = 1.0;
	double sig = 3.0;
	double x0  = 0.0;

//	double e1 = - amp*exp( - (x - x0)*(x - x0)/(2*sig*sig));

    double e1 = amp*(x - x0)*(x - x0);

	return e1;
}

static double fene1(double x)
{
	double amp = 0.3;
	double sig = 3.0;
	double x0  = -5.0;

//	double e1 = - amp*exp( - (x - x0)*(x - x0)/(2*sig*sig));

    double e1 = amp*(x - x0)*(x - x0);

	return e1;
}


void HaTests::model_mc_calc()
{
	Random rand_num_gen(240);

	int i;

    double ene_new;
	double step;
	double step_ap = 1.0;

	double x = 0.5;
	double ene = fene1(x);

	double kt = 1.0;

	FILE* fout = fopen("traj.dat","w");

	for(i = 0; i < 2000; i++)
	{
		step = step_ap*(rand_num_gen() - 0.5);
		
		double xnew = x + step;
        
		ene_new = fene1(xnew);
		
		double ksi = rand_num_gen();
		if( ene_new < ene || exp( - (ene_new-ene)/kt) > ksi ) 
		{
			x = xnew;
			ene = ene_new;
			fprintf(fout, " %9.4f %9.4f \n", x, ene); 
			PrintLog(" Point accepted x = %8.3f ene = %8.3f \n",x, ene); 
		}
	}

	fclose(fout);
}
