/*! \file ndo.cpp

    CNDO/INDO computation related functions

    \author Igor Kurnikov , University of Pittsburgh 
    \date 2004- 

*/

#include <mpi.h>

#include <string>
#include "harlemapp.h"
#include "hamatdb.h"
#include "haqchem.h"
#include "hamolset.h"
#include "g94_protos.h"
#include "math.h"

HaMatDB* HaQCMod::p_ndo_pars_db = NULL;

extern "C"{
extern void trmatd_(double* t, int* maxl, double* e);
extern void extrdm_(int* iout, int* iprint, int* icount,
					int* nbzdo, int* iuhf1, int* nttzdo,
					int* iflag, int* iextp,
                    double* rmsdp, double* pa, double* scr1, double* deltp);
}

//!< index of elements of symmetrical lower triangle stored matrix (1-based indexes i,j) return index in array 0-based)
int LInd(int i,int j){ return MaxFun(i,j)*(MaxFun(i,j)-1)/2 + MinFun(i,j) -1; }
//!< index of elements of symmetrical lower triangle stored matrix (0-based indexes i,j) return index in array 0-based)
int LInd0(int i,int j){ return MaxFun(i,j)*(MaxFun(i,j)+1)/2 + MinFun(i,j); }
int idx4s(int i,int j, int k, int l) 
//! index of coulomb integrals in atom coulom integral array (AO indexes 1-9) 0-based index in integral array
//! swap dx2_y2 <=> dxz   dxy <=> dyz
{
	if( i == 9) { i = 7; } else if( i == 7) { i = 9; } 
	if( i == 8) { i = 6; } else if( i == 6) { i = 8; } 
	if( k == 9) { k = 7; } else if( k == 7) { k = 9; } 
	if( k == 8) { k = 6; } else if( k == 6) { k = 8; } 
	if( j == 9) { j = 7; } else if( j == 7) { j = 9; } 
	if( j == 8) { j = 6; } else if( j == 6) { j = 8; } 
	if( l == 9) { l = 7; } else if( l == 7) { l = 9; } 
	if( l == 8) { l = 6; } else if( l == 6) { l = 8; } 

	if( k > l) { int ix = k; k = l ; l = ix; }	
   int kl = l*(l-1)/2 + k;
     if( i > j) { int ix = i; i = j ; j = ix; }	
   int ij = j*(j-1)/2 + i;
   
   return(LInd(kl,ij));
}

int idx4(int i,int j, int k, int l) 
//! index of coulomb integrals in atom coulom integral array (AO indexes 1-9) 0-based index in integral array
{
   int kl = l*(l-1)/2 + k;
   int ij = j*(j-1)/2 + i;
   
   return(LInd(kl,ij));
}


double HaQCMod::ZIndoGInt(double r, double gaa, double gbb)
//! function to compute electronic repulsion in CNDO approximating between coulomb and one center integral
//! Modified(by parameter 1.2) Motaga-Nishimoto Formula
{
      return 1.2/(r + 2.4/(gaa+gbb));
}


int HaQCMod::NDOExp(int ia , int l, int ovlp, int& n, HaVec_double& cf, HaVec_double& exp)
{
  double MuCI[48] = {1.2, 1.7,
		  0.65 , 0.975, 1.3,   1.625,  1.95, 2.275,  2.6,   2.925,
		  2.2/3, 0.95,  3.5/3, 4.15/3, 1.60, 5.45/3, 6.1/3, 2.25 ,
		  0.0  , 1.21,
		  1.23,  1.3,  1.3,   1.32, 1.36, 1.37, 1.423, 1.473, 1.482, 1.509,
		                0.0,   0.0,    0.0,  0.0,    0.0,   0.0,
		  0.0, 0.0, 
          1.275, 1.34, 1.397, 1.447,1.49, 1.47, 1.482, 1.497, 1.507, 1.699};

  double Coeff1[48] = {0.0, 0.0,
		  0.0 , 0.0, 0.0,   0.0 ,  0.0,  0.0,  0.0,   0.0,
		  0.0 , 0.0, 0.0,   0.0 ,  0.0,  0.0,  0.0,   0.0,
		  0.0 , 0.0,
          0.35922, 0.36461, 0.37378, 0.40714, 0.38984, 0.40379, 0.41333, 0.4212, 0.44729, 0.0,
		                0.0,   0.0,    0.0,  0.0,    0.0,   0.0,
		  0.0, 0.0, 
		  0.28655, 0.39924, 0.47775, 0.50051, 0.51241, 5.2248,  0.53881, 0.55017, 0.55763, 0.56165};
  
  double Coeff2[48] = {0.0, 0.0,
		  0.0 , 0.0, 0.0,   0.0 ,  0.0,  0.0,  0.0,   0.0,
		  0.0 , 0.0, 0.0,   0.0 ,  0.0,  0.0,  0.0,   0.0,
		  0.0 , 0.0,
		  0.76601, 0.75561, 0.74564, 0.73242, 0.72965, 0.71984, 0.71262, 0.70658, 0.69683, 0.0,
		                0.0,   0.0,    0.0,  0.0,    0.0,   0.0,
		  0.0, 0.0, 
		  0.82468, 0.71381, 0.65432, 0.62684, 0.59542, 0.58711, 0.57193, 0.56104, 0.55359, 0.54887};

  double MuDC[48] = {0.0, 0.0,
		  0.0 , 0.0, 0.0,   0.0 ,  0.0,  0.0,  0.0,   0.0,
		  0.0 , 0.0, 0.0,   0.0 ,  0.0,  0.0,  0.0,   0.0,
		  0.0 , 0.0,
          2.02,    2.497,   2.738,   2.966,   3.195,   3.383,  3.5803,  3.7765,  3.9721,  0.0, 
		                0.0,   0.0,    0.0,  0.0,    0.0,   0.0,
		  0.0, 0.0, 
		  1.879,   2.161,   2.389,   2.594,   2.795,   2.797,   2.978,   2.998,   3.321,   3.617};

  double MuDO1[48] = {0.0, 0.0,
		  0.0 , 0.0, 0.0,   0.0 ,  0.0,  0.0,  0.0,   0.0,
		  0.0 , 0.0, 0.0,   0.0 ,  0.0,  0.0,  0.0,   0.0,
		  0.0 , 0.0,
		  4.22244, 4.67,    5.05186, 5.13843, 5.76739, 6.06828, 6.38612, 6.70551, 6.79466, 0.0,
		                0.0,   0.0,    0.0,  0.0,    0.0,   0.0,
		  0.0, 0.0,
		  3.83698, 3.63923, 3.5908,  3.81446, 4.12429, 4.35709, 4.56078, 4.77173, 4.98896, 5.21258};

  double MuDO2[48] = {0.0, 0.0,
		  0.0 , 0.0, 0.0,   0.0 ,  0.0,  0.0,  0.0,   0.0,
		  0.0 , 0.0, 0.0,   0.0 ,  0.0,  0.0,  0.0,   0.0,
		  0.0 , 0.0,
		  1.74647, 1.98614, 2.17279, 2.07723, 2.50969, 2.61836, 2.74495, 2.87381, 2.76527, 0.0, 
		                0.0,   0.0,    0.0,  0.0,    0.0,   0.0,
		  0.0, 0.0, 
		  1.739,   1.80383, 1.7175,  1.86369, 2.15492, 2.26533, 2.36478, 2.47121, 2.58374, 2.70557};

      if( ia < 0 || ia > 48)
	  {
         ErrorInMod(" HaQCMod::NDOExp() ",
			         " Iat is out of range ");
		 return FALSE;
	  }
      else if(ia == 0) 
	  {
        n=0;
        cf[0]  = 0.0;
        cf[1]  = 0.0;
        exp[0] = 0.0;
        exp[1] = 0.0;
	  }
      else if(l <= 1) 
	  {
        n = 1;
        cf[0] = 1.0;
        cf[1] = 0.0;
        exp[0] = MuCI[ia-1];
        exp[1] = 0.0;
	  }
      else if(ndo_method <= 2 || ia <= 20)
	  {
        n=0;
        cf[0]  = 0.0;
        cf[1]  = 0.0;
        exp[0] = 0.0;
        exp[1] = 0.0;
	  }
      else if(ovlp)
	  {
        n = 2;
        cf[0]  = Coeff1[ia-1];
        cf[1]  = Coeff2[ia-1];
        exp[0] = MuDO1[ia-1];
        exp[1] = MuDO2[ia-1];
	  }
      else
	  {
		n = 1;
        cf[0] = 1.0;
        cf[1] = 0.0;
        exp[0] = MuDC[ia-1];
        exp[1] = 0.0;
      }
	  return TRUE;
}

static HaVec_double fact;       // fact[i] = i! - factorial of i
static HaVec_double bin_coef;   // binominal coefficients

static int psign[2] = { 1, -1 };

// scaling factor sigma/pi/delta for ZINDO/S, CNDO/S
      HaMat_double spi_scale; // scale factor for sigma/pi interations

	  double s_sigma = 1.0;
 	  double p_sigma = 1.267;
	  double d_sigma = 1.000;
	  double p_pi    = 0.585;
	  double d_pi    = 1.0;
	  double d_delta = 1.0;
 
	  const int nc[48] = { 1, 1,
	                 2, 2, 2, 2, 2, 2, 2, 2,
					 3, 3, 3, 3, 3, 3, 3, 3,
                     4, 4, 
					 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
					       4, 4, 4, 4, 4, 4,
					 5, 5,
					 5, 5, 5, 5, 5, 5, 5, 5, 5, 5}; // row number



int HaQCMod::CNDOInteg(int iprint,int iatom,int natoms, const HaVec_int& ian, const HaMat_double& c,
 				       double* ss, double* ss_scaled, double* gamma1, double* p_gamma2, double* gamma3,
				       const HaVec_int& ilst_bf_at, const HaVec_int& ifst_bf_at)
//
// iatom - force for which is computed
//
{	  
	  double cutoff = 1.0E-6;
	  HaMat_double gamma2;
	  gamma2.set_ext_alloc(p_gamma2,natoms,natoms);
    
      HaVec_int lc(9);
	  lc(1) = 0; lc(2) = 1; lc(3) = 1; lc(4) = 1; 
	  lc(5) = 2; lc(6)= 2; lc(7) = 2; lc(8) = 2; lc(9) = 2;
      
	  HaVec_int mc(9);
	  mc(1) = 0; mc(2) = 1; mc(3) = -1; mc(4) = 0; 
	  mc(5) = 0; mc(6) = 1; mc(7) = -1; mc(8) = 2; mc(9) = -2;

      if(fact.size() < 30)
	  {
		 fact.newsize(30);
         int i;
		 fact[0] = 1.0;
		 for( i = 1; i < 30; i++)
		 {
			fact[i] = i*fact[i-1];
		 }
	  }

	  if(bin_coef.size() < 465)
	  {
		  bin_coef.newsize(465);
		  cbincf_(bin_coef.begin(),&fact[0]);
	  }

	  HaVec_double GamSSS(48),GamSSD(48),GamSDD(48);

      GamSSS(1) = 12.85;  /* H */  GamSSD(1) = 0.0;  GamSDD(1) = 0.0;  // CNDO/S JAFFE_JCP(68)_48_4050
      GamSSS(2) = 0.0;    /* He */ GamSSD(2) = 0.0;  GamSDD(2) = 0.0; 
// 2-nd row	  
	  GamSSS(3) = 0.0;    /* Li */ GamSSD(3) = 0.0;  GamSDD(3) = 0.0; 
	  GamSSS(4) = 0.0;    /* Be */ GamSSD(4) = 0.0;  GamSDD(4) = 0.0; 
	  GamSSS(5) = 0.0;    /* B */  GamSSD(5) = 0.0;  GamSDD(5) = 0.0; 
	  GamSSS(6) = 11.11;  /* C */  GamSSD(6) = 0.0;  GamSDD(6) = 0.0;   // CNDO/S JAFFE_JCP(68)_48_4050
	  GamSSS(7) = 12.01;  /* N */  GamSSD(7) = 0.0;  GamSDD(7) = 0.0;   // CNDO/S JAFFE_JCP(68)_48_4050
	  GamSSS(8) = 13.0;   /* O */  GamSSD(8) = 0.0;  GamSDD(8) = 0.0;   // CNDO/S JAFFE_JCP(68)_48_4050
	  GamSSS(9) = 14.0;   /* F */  GamSSD(9) = 0.0;  GamSDD(9) = 0.0; 
      GamSSS(10) = 0.0;   /* Ne */ GamSSD(10) = 0.0;  GamSDD(10) = 0.0; 
// 3-rd row 
      GamSSS(11) = 0.0;           GamSSD(11) = 0.0;  GamSDD(11) = 0.0; 
	  GamSSS(12) = 4.79;          GamSSD(12) = 0.0;  GamSDD(12) = 0.0; 
	  GamSSS(13) = 0.0;           GamSSD(13) = 0.0;  GamSDD(13) = 0.0; 
	  GamSSS(14) = 0.0;           GamSSD(14) = 0.0;  GamSDD(14) = 0.0; 
	  GamSSS(15) = 0.0;           GamSSD(15) = 0.0;  GamSDD(15) = 0.0; 
	  GamSSS(16) = 10.09;         GamSSD(16) = 0.0;  GamSDD(16) = 0.0; 
	  GamSSS(17) = 0.0;           GamSSD(17) = 0.0;  GamSDD(17) = 0.0; 
	  GamSSS(18) = 0.0;           GamSSD(18) = 0.0;  GamSDD(18) = 0.0; 
// 4-th row
      GamSSS(19) = 0.0;  /* K */  GamSSD(19) = 0.0;  GamSDD(19) = 0.0; 
	  GamSSS(20) = 3.25; /* Ca */ GamSSD(20) = 4.00; GamSDD(20) = 6.03;    // ZINDO - ZERNER_JACS(80)_102_589
      GamSSS(21) = 3.89; /* Sc */ GamSSD(21) = 4.71; GamSDD(21) = 7.02;    // ZINDO - ZERNER_JACS(80)_102_589
	  GamSSS(22) = 4.5;  /* Ti */ GamSSD(22) = 5.38; GamSDD(22) = 7.98;    // ZINDO - ZERNER_JACS(80)_102_589
	  GamSSS(23) = 5.07; /* V  */ GamSSD(23) = 6.01; GamSDD(23) = 8.91;    // ZINDO - ZERNER_JACS(80)_102_589
	  GamSSS(24) = 5.6;  /* Cr */ GamSSD(24) = 6.6;  GamSDD(24) = 9.81;    // ZINDO - ZERNER_JACS(80)_102_589
	  GamSSS(25) = 6.09; /* Mn */ GamSSD(25) = 7.16; GamSDD(25) = 10.68;   // ZINDO - ZERNER_JACS(80)_102_589
	  GamSSS(26) = 6.54; /* Fe */ GamSSD(26) = 7.68;  GamSDD(26) = 11.52;  // ZINDO - ZERNER_JACS(80)_102_589
	  GamSSS(27) = 6.96; /* Co */ GamSSD(27) = 8.16;  GamSDD(27) = 12.32;  // ZINDO - ZERNER_JACS(80)_102_589
	  GamSSS(28) = 7.34; /* Ni */ GamSSD(28) = 8.61;  GamSDD(28) = 13.10;  // ZINDO - ZERNER_JACS(80)_102_589
	  GamSSS(29) = 7.68; /* Cu */ GamSSD(29) = 9.01;  GamSDD(29) = 13.84;  // ZINDO - ZERNER_JACS(80)_102_589
	  GamSSS(30) = 7.98; /* Zn */ GamSSD(30) = 9.39;  GamSDD(30) = 14.55;  // ZINDO - ZERNER_JACS(80)_102_589
// 
      GamSSS(31) = 0.0; /*   */     GamSSD(31) = 0.0;  GamSDD(31) = 0.0; 
      GamSSS(32) = 0.0; /*   */     GamSSD(32) = 0.0;  GamSDD(32) = 0.0;
	  GamSSS(33) = 0.0; /*   */     GamSSD(33) = 0.0;  GamSDD(33) = 0.0;
	  GamSSS(34) = 0.0; /*   */     GamSSD(34) = 0.0;  GamSDD(34) = 0.0;
	  GamSSS(35) = 0.0; /*   */     GamSSD(35) = 0.0;  GamSDD(35) = 0.0;
	  GamSSS(36) = 0.0; /*   */     GamSSD(36) = 0.0;  GamSDD(36) = 0.0;
// 4-th row
	  GamSSS(37) = 0.0; /* Y  */     GamSSD(37) = 0.0;  GamSDD(37) = 0.0;
	  GamSSS(38) = 0.0; /*    */     GamSSD(38) = 0.0;  GamSDD(38) = 0.0;
	  GamSSS(39) = 0.0; /*    */     GamSSD(39) = 0.0;  GamSDD(39) = 0.0;
	  GamSSS(40) = 0.0; /*    */     GamSSD(40) = 0.0;  GamSDD(40) = 0.0;
	  GamSSS(41) = 0.0; /*    */     GamSSD(41) = 0.0;  GamSDD(41) = 0.0;
	  GamSSS(42) = 0.0; /*    */     GamSSD(42) = 0.0;  GamSDD(42) = 0.0;
	  GamSSS(43) = 0.0; /*    */     GamSSD(43) = 0.0;  GamSDD(43) = 0.0;
	  GamSSS(44) = 0.0; /*    */     GamSSD(44) = 0.0;  GamSDD(44) = 0.0;
	  GamSSS(45) = 0.0; /*    */     GamSSD(45) = 0.0;  GamSDD(45) = 0.0;
	  GamSSS(46) = 0.0; /*    */     GamSSD(46) = 0.0;  GamSDD(46) = 0.0;
	  GamSSS(47) = 0.0; /*    */     GamSSD(47) = 0.0;  GamSDD(47) = 0.0;
	  GamSSS(48) = 0.0; /* Cd */     GamSSD(48) = 0.0;  GamSDD(48) = 0.0;

      double ToEv = HARTREE_TO_EV;

      int nattt = (natoms*(natoms+1))/2;
      int i,j;
	  int k,l;
      if(iatom == 0)
	  {
		  for(k = 0; k < nattt; k++)
			  gamma1[k] = 0.0;
		  if(ndo_method == harlem::qc::ZINDO_1 || ndo_method == harlem::qc::ZINDO_S ) 
		  {
			  for(k = 0; k < nattt; k++)
				  gamma3[k] = 0.0;
			  for(k = 1; k <= natoms; k++)
			  {
				  for(l = 1; l <= natoms; l++)
				  {
					  gamma2(k,l) = 0.0;
				  }
			  }
		  }
      }

//     Loop over pairs of atoms.

      HaMat_double pairs(9,9),pairs_sc(9,9),t(9,9),temp(9,9),temp_sc(9,9);

      for (k = 1 ; k <= natoms; k++)  
	  {
        if(ian(k) == 0) continue;
        for( l = k; l <= natoms; l++)
		{
		  if(iatom != 0 && k != iatom && l != iatom ) continue;
 
		  int ank = ian(k);
          int anl = ian(l);
          int nzk, nzl;

		  if(ian(l) == 0) anl = 1; // if dummy atom set as hydrogen

		  Vec3D c1,c2,e;
          for(i = 1; i <=3; i++)
		  {
            c1[i-1] = c(i,k);  // coordinates of the first atom
            c2[i-1] = c(i,l);  // coordinates of the second atom
		  }

//        Calculate unit vector along interatom axis.
		  e = c2-c1;
		  double r = e.length();
          e.normalize();

	      if(spi_scale.num_rows() != 9) // scale factor for sigma/pi interations
		  {
		     spi_scale.newsize(9,9);
		     spi_scale = 1.0;
	         spi_scale(1,1) = s_sigma;
//	         spi_scale(1,4) = (s_sigma + p_sigma)*0.5;  spi_scale(4,1) = spi_scale(1,4);	    
             spi_scale(2,2) = p_pi;
	         spi_scale(3,3) = p_pi;
             spi_scale(4,4) = p_sigma;
	         spi_scale(5,5) = d_sigma;
	         spi_scale(6,6) = d_pi;
		     spi_scale(7,7) = d_pi;
		     spi_scale(8,8) = d_delta;
	       	 spi_scale(9,9) = d_delta;
// It seems in ZINDO only scale s-s p-p and d-d not s-p p-d etc
//		     spi_scale(1,5) = (s_sigma + d_sigma)*0.5;  spi_scale(5,1) = spi_scale(1,5);
//		     spi_scale(4,5) = (p_sigma + d_sigma)*0.5;  spi_scale(5,4) = spi_scale(4,5);
//		     spi_scale(2,6) = (p_pi + d_pi)*0.5;        spi_scale(6,2) = spi_scale(2,6);
//		     spi_scale(3,7) = (p_pi + d_pi)*0.5;        spi_scale(7,3) = spi_scale(3,7);		
		  } 

		  HaVec_double czspk(2),ezspk(2),czspl(2),ezspl(2);
          HaVec_double czdk(2), ezdk(2), czdl(2), ezdl(2);

          if(ndo_method == harlem::qc::ZINDO_S) 
		  {
             double GKSS = GamSSS(ank)/ToEv;
             double GLSS = GamSSS(anl)/ToEv;
             double GKDD = GamSDD(ank)/ToEv;
             double GLDD = GamSDD(anl)/ToEv;
             gamma1[LInd(k,l)] = ZIndoGInt(r,GKSS,GLSS);
             if(ian(l) >= 21) gamma2(k,l) = ZIndoGInt(r,GKSS,GLDD);
             if(ian(k) >= 21) 
			 {
                if(k == l)
				{
                   gamma2(k,l) = GamSSD(ank)/ToEv;
				}
                else
				{
                   gamma2(l,k) = ZIndoGInt(r,GLSS,GKDD);
				}
                if(ian(l) >= 21) gamma3[LInd(k,l)]= ZIndoGInt(r,GKDD,GLDD);
             }
		  }
          else
		  {
//           Computation of 1-center Coulomb integrals over Slater s functions
             int n1sp = nc[ank-1];
             int n2sp = nc[anl-1];
			 int n1d  = n1sp-1;
			 int n2d  = n2sp-1;
	 
			 NDOExp(ank,0, FALSE, nzk, czspk,ezspk);
			 NDOExp(ank,2, FALSE, nzk, czdk,ezdk);
			 NDOExp(anl,0, FALSE, nzl, czspl,ezspl);
			 NDOExp(anl,2, FALSE, nzl, czdl,ezdl);

             gamma1[LInd(k,l)] = gint1_(&ezspk(1),&ezspl(1),&n1sp,&n2sp,&r,&fact[0],bin_coef.begin());
			 {
                if(ian(l) >= 21) gamma2(k,l) = gint1_(&ezspk(1),&ezdl(1),&n1sp,&n2d,&r,&fact[0],bin_coef.begin());
                if(ian(k) >= 21 && k != l) gamma2(l,k) = gint1_(&ezspl(1),&ezdk(1),&n2sp,&n1d,&r,&fact[0],bin_coef.begin());
                if(ian(k) >= 21 && ian(l) >= 21) gamma3[LInd(k,l)] = gint1_(&ezdk(1),&ezdl(1),&n1d,&n2d,&r,&fact[0],bin_coef.begin());
             }
		  }

          if(ian(l) == 0 ) continue;
               
          int llk = ifst_bf_at(k);
          int lll = ifst_bf_at(l);
          int ulk = ilst_bf_at(k);
          int ull = ilst_bf_at(l);
          int norbk = ulk - llk + 1;
          int norbl = ull - lll + 1;

//         Loop thru pairs of basis functions, one on each atom

          if(norbk == 0 || norbl == 0) continue;

          NDOExp(ank,0,  TRUE, nzk,czspk,ezspk); // Get SP exponents and coef of the first atom
		  NDOExp(ank,2,  TRUE, nzk,czdk,ezdk);   // Get D exponents and coef of the first atom
		  NDOExp(anl,0,  TRUE, nzl,czspl,ezspl); // Get SP exponents and coef of the second atom
		  NDOExp(anl,2,  TRUE, nzl,czdl,ezdl);   // Get D exponents and coef of the second atom

		  for(i=1; i <= norbk; i++)
		  {
             int nck = nc[ank-1]; //  ANK - first atom type, NC(ANK) - row number (?)
             if( i > 4 ) nck--;
//             nck = MinFun(nck,3);
             for(j = 1; j <= norbl; j++)
			 {
                int ncl = nc[anl-1];

                if( j > 4) ncl--;
//                ncl = MinFun(ncl,3);
                if(k == l)
				{
                  if(i == j) 
                    pairs(i,j) = 1.0;
                  else
                    pairs(i,j) = 0.0;
				}
                else if(mc(i) != mc(j))
                   pairs(i,j) = 0.0;
                else if(mc(i) < 0) 
                   pairs(i,j) = pairs(i-1,j-1);
                else
				{
		           pairs(i,j) = 0.0;
				   int im,jm;
				   double ezk_r, ezl_r, ck, cl;
				   for( im = 1; im <= 2; im++)
				   {
					   for(jm=1; jm <= 2; jm++)
					   {
						  if( i <= 4 && j <= 4)
						  {
                             ezk_r  = ezspk(im)*r; ck = czspk(im);
                             ezl_r  = ezspl(jm)*r; cl = czspl(jm);
						  }
						  else if( i > 4 && j <= 4)
						  {
                             ezk_r  = ezdk(im)*r; ck = czdk(im);
                             ezl_r  = ezspl(jm)*r; cl = czspl(jm);
						  }
						  else if( i <= 4 && j > 4)
						  {
                             ezk_r  = ezspk(im)*r; ck = czspk(im);
                             ezl_r  = ezdl(jm)*r; cl = czdl(jm);
						  }
						  else
						  {
                             ezk_r  = ezdk(im)*r; ck = czdk(im);
                             ezl_r  = ezdl(jm)*r; cl = czdl(jm);
						  }
						  if( fabs(ck) < 1.0E-6 || fabs(cl) < 1.0E-6 ) continue;

						   int lg = 0;
						   double ss1 = psign[ (lc(j)+mc(j)) % 2 ]*ssz_(&nck,&ncl,&lc(i),&lc(j),&mc(i),&ezk_r,&ezl_r,
							            fact.begin(),bin_coef.begin(),&lg); 
                           pairs(i,j) += ck*cl*ss1;  
					   }
				   }
				}
			 } // end of cycle on j
          } // end of cycle on i
          
          pairs_sc = pairs;
 		  int n1 = spi_scale.num_rows();
 		  int n2 = pairs_sc.num_rows();
 		  int nn = MinFun(n1,n2);
 
 		  if( k != l)
 		  {
 			for(i=1; i <= nn; i++)
 			{
 				for(j=1; j <= nn; j++)
 				{
 					pairs_sc(i,j) *= spi_scale(i,j);
 				}
 			}
 		  }

          int maxl = MaxFun(lc(norbk),lc(norbl));
          if(r > cutoff)
		  {
             trmatd_(t.v(),&maxl,&e[0]);  // Get transformation matrix T for direction E(3)

			 int kk;
             for(i = 1; i <= norbk; i++)
			 {
                for(j = 1; j <= norbl; j++)
				{
                   temp(i,j) = 0.0;
				   temp_sc(i,j) = 0.0;
                   for(kk = 1; kk <= norbl; kk++)
				   {
                      temp(i,j) += t(j,kk)*pairs(i,kk);
 					  temp_sc(i,j) += t(j,kk)*pairs_sc(i,kk);
				   }
				}
			 }
             for(i = 1; i <= norbk; i++)
			 {
				 for(j = 1; j <= norbl; j++)
				 {
                    pairs(i,j) = 0.0;
					pairs_sc(i,j) = 0.0;
                    for(kk = 1; kk <= norbk; kk++)
					{
                       pairs(i,j)    += t(i,kk)*temp(kk,j);
 					   pairs_sc(i,j) += t(i,kk)*temp_sc(kk,j);
					}
				 }
			 }
          }
//        Fill S matrix
          for(i = 1; i <= norbk; i++)
		  {
            int llkp = llk +i-1;
            for(j = 1; j <= norbl; j++)
			{
              int lllp = lll +j-1;
              ss[LInd(llkp,lllp)] = pairs(i,j);
			  ss_scaled[LInd(llkp,lllp)] = pairs_sc(i,j);
			}
          }
        }  // end of cycle on l
     } // end of cycle on k

    return TRUE;
}


// 1/2(ionization potential + electron affinity) 
// parameters from D(N-1) S(1) 
// S-orbitals:

const int MAX_AT_TYPE = 105;

//    CNDO/2 AND INDO/2 : 1/2(IONIZATION POTENTIAL + ELECTRON AFFINITY)

static double ie_ss_1[MAX_AT_TYPE] =
{
   -7.1761,                                            -100.0,     // 1. H-He
   -3.1055,   -5.94557,                                            // 2  Li-Be  
   -9.59407, -14.051,   -19.31637 ,-25.39017, -32.2724,-100.0,     //    B-Ne
   -2.805,    -5.222,                                              // 3. Na-Mg    
   -8.288,   -11.157,   -13.551,   -16.328,   -19.841, -100.0,     //    Al-Ar   
   -2.385,    -3.400,                                              // 4. K -Ca    
              -4.345,    -4.44,     -5.24,     -5.225,   -5.335,   //    Sc-Mn   
//              -5.725,    -5.12,     -3.84,     -3.985,   -6.220,   //    Fe-Zn  in Zindo
               100.0,    -5.12,     -3.84,     -3.985,   -6.220,   //    Fe-Zn  Igor playing    
 -100.0,     -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    Ga-Kr    
 -100.0,     -100.0,                                               // 5. Rb-Sr   
             -100.0,   -100.0,    -100.0,       -3.93, -100.0,     //    Y -Tc
             -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    Ru-Cd 
 -100.0,     -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    In-Xe  			 
 -100.0,     -100.0,                                               // 6. Cs-Ba   
             -100.0,                                               //    La
             -100.0, -100.0, -100.0, -100.0, -100.0, -100.0,-100.0,//    Ce-Gd
             -100.0, -100.0, -100.0, -100.0, -100.0, -100.0,-100.0,//    Tb-Lu
                       -100.0,    -100.0,     -100.0,  -100.0,     //    Hf-Re
             -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    Os-Hg                                     
 -100.0,     -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    Tl-Rn			 
 -100.0,     -100.0,                                               // 7. Fr-Ra   
             -100.0,                                               //    Ac
             -100.0, -100.0, -100.0, -100.0, -100.0, -100.0,-100.0,//    Th-Cm
             -100.0, -100.0, -100.0, -100.0, -100.0, -100.0,-100.0,//    Bk-Lr
			 -100.0,    -100.0};                                   //    Rf-Db

// P-orbitals:
static double ie_pp_1[MAX_AT_TYPE] =
{
  -100.0,                                              -100.0,     // 1. H-He
   -1.258,    -2.563,                                              // 2  Li-Be  
   -4.001,    -5.572,    -7.275,    -9.111,   -11.080, -100.0,     //    B-Ne
   -1.565,    -2.100,                                              // 3. Na-Mg
   -2.950,    -4.270,    -6.080,    -7.385,    -9.380, -100.0,     //    Al-Ar 
   -1.265,    -1.989,                                              // 4. K -Ca 
              -2.399,    -1.978,    -2.542,    -2.971,   -2.612,   //    Sc-Mn
//			  -3.655,    -2.550,    -1.317,    -1.165,   -1.640,   //    Fe-Zn
			 100.0,      -2.550,    -1.317,    -1.165,   -1.640,   //    Fe-Zn  igor playing
 -100.0,     -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    Ga-Kr    
 -100.0,     -100.0,                                               // 5. Rb-Sr   
             -100.0,   -100.0,    -100.0,       -0.71, -100.0,     //    Y -Tc
             -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    Ru-Cd 
 -100.0,     -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    In-Xe  			 
 -100.0,     -100.0,                                               // 6. Cs-Ba   
             -100.0,                                               //    La
             -100.0, -100.0, -100.0, -100.0, -100.0, -100.0,-100.0,//    Ce-Gd
             -100.0, -100.0, -100.0, -100.0, -100.0, -100.0,-100.0,//    Tb-Lu
                       -100.0,    -100.0,     -100.0,  -100.0,     //    Hf-Re
             -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    Os-Hg                                     
 -100.0,     -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    Tl-Rn			 
 -100.0,     -100.0,                                               // 7. Fr-Ra   
             -100.0,                                               //    Ac
             -100.0, -100.0, -100.0, -100.0, -100.0, -100.0,-100.0,//    Th-Cm
             -100.0, -100.0, -100.0, -100.0, -100.0, -100.0,-100.0,//    Bk-Lr
			 -100.0,    -100.0};                                   //    Rf-Db

// D-orbitals:
static double ie_dd_1[MAX_AT_TYPE] =
{
 -100.0,                                               -100.0,     // 1. H-He
 -100.0,    -100.0,                                                // 2  Li-Be  
 -100.0,    -100.0,    -100.0,    -100.0,    -100.0,   -100.0,     //    B-Ne
   -0.445,   -0.65,                                                // 3. Na-Mg
   -0.890,   +0.625,     -0.025,    -0.405,    +1.900, -100.000,   //    Al-Ar 
   -0.635,   -0.770,                                               // 4. K -Ca
             -2.15,      -2.70,     -3.245,    -3.545,   -3.410,   //    Sc-Mn
             -3.735,     -4.210,    -3.975,    -6.240,  -13.850,   //    Fe-Zn
 -100.0,     -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    Ga-Kr    
 -100.0,     -100.0,                                               // 5. Rb-Sr   
             -100.0,   -100.0,    -100.0,       -4.53, -100.0,     //    Y -Tc
             -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    Ru-Cd 
 -100.0,     -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    In-Xe  			 
 -100.0,     -100.0,                                               // 6. Cs-Ba   
             -100.0,                                               //    La
             -100.0, -100.0, -100.0, -100.0, -100.0, -100.0,-100.0,//    Ce-Gd
             -100.0, -100.0, -100.0, -100.0, -100.0, -100.0,-100.0,//    Tb-Lu
                       -100.0,    -100.0,     -100.0,  -100.0,     //    Hf-Re
             -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    Os-Hg                                     
 -100.0,     -100.0,   -100.0,    -100.0,     -100.0,  -100.0,     //    Tl-Rn			 
 -100.0,     -100.0,                                               // 7. Fr-Ra   
             -100.0,                                               //    Ac
             -100.0, -100.0, -100.0, -100.0, -100.0, -100.0,-100.0,//    Th-Cm
             -100.0, -100.0, -100.0, -100.0, -100.0, -100.0,-100.0,//    Bk-Lr
			 -100.0,    -100.0};                                   //    Rf-Db
		  
		  
static double beta0_sp_1[MAX_AT_TYPE] =  // off-diagonal interaction prefactor for S and P orbitals 
{
//  -11.8,                                                  -1.0,    // 1. H-He  ZINDO/1 parameters
    -12.0,                                                  -1.0,    // ZINDO/S Del Bene-Jaffe for spectra - default in ZINDO
	-5.0,      -13.0,                                               // 2  Li-Be  
//   -17.0,     -21.0,    -25.0,     -35.0,      -36.0,     -1.0,    //    B-Ne first appear in ZINDO PARAMB
   -8.0,     -17.0,    -26.0,     -34.0,      -44.0,     -1.0,    //    B-Ne Del Bene-Jaffe for spectra - default in ZINDO
   -5.0,       -6.0,                                               // 3. Na-Mg
//   -7.0,       -9.0,    -10.0,     -11.5,      -11.0,     -1.0,    //    Al-Ar  first appear in ZINDO PARAMB
     -7.0,       -9.0,    -15.0,     -15.0,      -11.0,     -1.0,    //    Al-Ar  Del Bene-Jaffe for spectra - default in ZINDO
   -1.0,        2.0,                                               // 4. K -Ca
//               -2.0,     -7.0,     -12.0,      -17.0,    -22.0,    //    Sc-Mn  From Zindo param - from Clark??
//              -26.0,    -29.0,     -32.0,      -35.0,     -1.0,    //    Fe-Zn  From Zindo param - from Clark??
               -1.0,     -1.0,      -1.0,       -1.0,     -1.0,    //    Sc-Mn  KIRCHNER-BACON-ZERNER VALUES OF 1977  ??
               -1.0,     -1.0,      -1.0,       -1.0,    -10.0,    //    Fe-Zn  KIRCHNER-BACON-ZERNER VALUES OF 1977  ??
   -1.0,       -1.0,     -1.0,      -1.0,       -1.0,     -1.0,    //    Ga-Kr    
   -1.0,       -1.0,                                               // 5. Rb(37)-Sr   
               -1.0,     -1.0,      -1.0,       -1.0,     -1.0,    //    Y -Tc
               -2.0,     -1.0,      -1.0,       -1.0,     -1.5,    //    Ru(44)-Cd 
   -1.0,       -1.0,     -1.0,      -1.0,       -1.0,     -1.0,    //    In-Xe  			 
   -1.0,       -1.0,                                               // 6. Cs-Ba   
               -1.0,                                               //    La
               -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,  -1.0,//    Ce-Gd
               -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,  -1.0,//    Tb-Lu
                         -1.0,      -1.0,       -1.0,     -1.0,    //    Hf-Re
               -1.0,     -1.0,      -1.0,       -1.0,     -1.0,    //    Os-Hg                                     
   -1.0,       -1.0,     -1.0,      -1.0,       -1.0,     -1.0,    //    Tl-Rn			 
   -1.0,       -1.0,                                               // 7. Fr-Ra   
               -1.0,                                               //    Ac
               -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,  -1.0,//    Th-Cm
               -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,  -1.0,//    Bk-Lr
			   -1.0,      -1.0};                                   //    Rf-Db

static double beta0_dd_1[MAX_AT_TYPE] =  // off-diagonal interaction prefactor for S and P orbitals
{
  -11.8,                                                  -1.0,    // 1. H-He
   -1.0,       -1.0,                                               // 2  Li-Be  
   -1.0,       -1.0,     -1.0,      -1.0,       -1.0,     -1.0,    //    B-Ne
   -7.4,      -11.4,                                               // 3. Na-Mg - Zindo pars -arbitrary
   -1.0,       -1.0,     -1.0,      -1.0,       -1.0,     -1.0,    //    Al-Ar 
   -1.0,       -1.0,                                               // 4. K -Ca
              -18.0,    -19.0,     -20.0,      -21.0,    -22.0,    //    Sc-Mn
			  -23.0,    -31.0,     -32.0,      -33.0,    -34.0,    //    Fe-Zn
   -1.0,       -1.0,     -1.0,      -1.0,       -1.0,     -1.0,    //    Ga-Kr    
   -1.0,       -1.0,                                               // 5. Rb-Sr   
//              -14.14,   -17.03,    -19.81,     -21.83,   -23.55,   //    Y(39) -Tc
//              -26.29,   -27.17,    -27.59,     -27.94,   -27.94,   //    Ru-Cd 
               -7.0,    -10.0,     -13.00,     -15.00,   -17.00,   //    Y(39) -Tc  Zindo default
              -20.0,    -21.0,     -21.50,     -22.00,   -22.00,   //    Ru-Cd      Zindo default      
   -1.0,       -129.0,     -1.0,      -1.0,       -1.0,     -1.0,  //    In-Xe  			 
   -1.0,       -1.0,                                               // 6. Cs-Ba   
               -1.0,                                               //    La
               -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,  -1.0,//    Ce-Gd
               -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,  -1.0,//    Tb-Lu
                         -1.0,      -1.0,       -1.0,     -1.0,    //    Hf-Re
               -1.0,     -1.0,      -1.0,       -1.0,     -1.0,    //    Os-Hg                                     
   -1.0,       -1.0,     -1.0,      -1.0,       -1.0,     -1.0,    //    Tl-Rn			 
   -1.0,       -1.0,                                               // 7. Fr-Ra   
               -1.0,                                               //    Ac
               -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,  -1.0,//    Th-Cm
               -1.0,   -1.0,   -1.0,   -1.0,   -1.0,   -1.0,  -1.0,//    Bk-Lr
			   -1.0,     -1.0};                                    //    Rf-Db

                           //    Rf-Db

static double g1sp_val[18] = 
{ 0.0,      0.0,
  0.092012, 0.1407,   0.199265, 0.267708, 0.346029, 0.43423,  0.532305, 0.0, 
  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0}; 

static double f2pp_val[18] = 
{ 0.0,      0.0,
  0.049865, 0.089125, 0.13041,  0.17372,  0.219055, 0.266415, 0.31580,   0.0, 
  0.0,      0.0,      0.0,      0.0,      0.0,      0.0,      0.0,       0.0}; 

static double beta0_val[18] =  // this looks as CNDO/2 values  CNDO/S ZINDO values are different beta_H=-12, beta_C=-17, beta_N=-26 see JAFFE_68 
{                              // beta_O = -45 eV
	-9.0,      0.0,            
   -9.0,    -13.0,    -17.0,    -21.0,    -25.0,    -31.0,    -39.0,       0.0, 
   -7.7203,  -9.4471, -11.3011, -13.065,  -15.070,  -18.150,  -22.330,     0.0};

int ll_bf[9] = {0,1,1,1,2,2,2,2,2};

HaMat_double ip_dnm1_s1;
HaMat_double ip_dnm2_s2;
HaMat_double ip_dn_s0;
HaMat_double ipea_dnm1_s1;
HaMat_double ipea_dnm2_s2;
HaMat_double ip_coef;
HaMat_double gf_integ;  // G1(S,P) F2(P,P) G2(S,D) G1(P,D) F2(P,D) G3(P,D) F2(D,D) F4(D,D)                   
HaMat_double rr_integ; // r1sppd(i), r2sddd(i), r2sdpp(i)
HaMat_double beta_zindo_1;
HaMat_double beta_zindo_s;

double
HaQCMod::GetNDOValEl(int elem, double& ns_val, double& np_val, double& nd_val,double& nf_val)
{
	double nval = 0.0;
	ns_val = 0.0;
	np_val = 0.0;
	nd_val = 0.0;
	nf_val = 0.0;
	
	double ncore;
	if( elem <= 2) ncore = 0;
	else if( elem <= 10) ncore = 2;
	else if( elem <= 18) ncore = 10;
	else if( elem <= 36) ncore = 18;
	else if( elem <= 54) ncore = 36;
	else if( elem <= 86) ncore = 54;
	else ncore = 86;
	
	nval = elem - ncore;
	if( elem == 1) 
	{
		ns_val = 1;
	}
	else if( elem < 20)
	{
		ns_val = MinFun(2.0,nval);
		np_val = MaxFun(nval-2.0,0.0);
	}
	else if( elem <= 30)
	{
		if( 1 ) // default config mixing (ZINDO param ISW2)
		{
			double cf1 = ip_coef.r0(elem-1,1);
			double cf2 = ip_coef.r0(elem-1,2);
		    np_val = 0.0;
		    ns_val = cf1 + 2.0* cf2;
            nd_val = nval - ns_val;
		}
		else
		{
			np_val = 0.0;
			ns_val = 1.0;  // S(1)D(N-1) config
//            ns_val = 2.0; // S(2)D(N-2) config
//            ns_val = 0.0; // S(0)D(N) config
            nd_val = nval - ns_val;
		}
	    
		if( elem == 30)
		{
// may have special treatment see ZINDO
		}
	}
	else if( elem <= 38)
	{
        ns_val = MinFun(2.0,nval);
        np_val = MaxFun(nval-2.0,0.0);
	}
	else if( elem <= 48 )
	{
		double cf1 = ip_coef.r0(elem-1,1);
		double cf2 = ip_coef.r0(elem-1,2);
		double cf3 = ip_coef.r0(elem-1,3);
		double cf4 = ip_coef.r0(elem-1,4);
		np_val = 0.0;
		ns_val = cf1 + 2.0* cf2 + cf4;
		np_val = cf4;
        nd_val = (nval - 1.0)*cf1 + (nval - 2.0)*cf2 + 
				  nval* cf3 + (nval - 1.0)*cf4;
		if( elem == 48)
		{
			np_val = 0.0;
            ns_val = cf1 + 2.0* cf2;
            nd_val = nval - ns_val;
		}
	}
	else if( elem <= 56 )
	{
        ns_val = MinFun(2.0,nval);
        np_val = MaxFun(nval-2.0,0.0);        
	}
	else if( elem <= 71)
	{
//		for the LANTHANIDES, F(N-3)D(1) S(2) dominates
        ns_val = 2.0;
		np_val = 0.0;
		nd_val = 1.0;
		nf_val = nval - 3.0;
	}
	else if( elem <= 88)
	{
       ns_val = MinFun(2.0,nval);
       np_val = MaxFun(nval-2.0,0.0);        
	}
	else
	{
		double cf1 = ip_coef.r0(elem-1,1);
		double cf2 = ip_coef.r0(elem-1,2);
		double cf3 = ip_coef.r0(elem-1,3);
		double cf4 = ip_coef.r0(elem-1,4);
        ns_val = 2.0*(cf1 + cf2 + cf3 + cf4);
		np_val = cf2;
		nd_val = cf1 + 2.0* cf3 + 3.0*cf4;
		nf_val = nval - ns_val - np_val - nf_val;
	}
	return nval;
} 


double
HaQCMod::RRIntSTO(int it, int n1, double a1, int n2, double a2, int n3, 
				  double a3, int n4, double a4)    
//!  general slater condon radial integral                          
//!  for STOs characterized by n,l,m  and exp a1
//!  it - type of the integral
//!
{
      double fn =sqrt( pow(2.0*a1, 2*n1+1) );                                     
      double fn1=sqrt( pow(2.0*a2, 2*n2+1) );                                    
      double fn2=sqrt( pow(2.0*a3, 2*n3+1) );                                    
      double fn3=sqrt( pow(2.0*a4, 2*n4+1) );
	  
      fn=fn*fn1*fn2*fn3/sqrt(fact(2*n1+1)*fact(2*n2+1)*fact(2*n3+1)*fact(2*n4+1));
	  
      int n=n1+n3;                                                           
      int m=n2+n4;                                                          
      double a=a1+a3;                                                          
      double b=a2+a4;                                                          
      int nl=n+ it +1;                                                          
      int ml=n-it;                                                           
      int mn=m+n+1;                                                         
      double fa=fact(nl)*fact(m-it)/( pow(a,nl)*pow(b,m-it) );                          
      double aa=1.0;                                                            
      double bb=pow(a+b,mn);                                                     
      int i;
	  double f;
	  for( i=1; i <= nl; i++)
	  {
         aa = aa*a;                                                           
         bb = bb/(a+b);                                                      
         f=fact(mn-i)/(aa*bb);                                              
         fa = fa - fact(nl)*f/fact(nl-i+1);                                     
         if(i <= ml) fa=fa+fact(ml)*f/fact(ml-i+1);                         
      }                                                          
                                                             
      f=fn*fa*2.0/3.0;                                                    
      return f;
}

int HaQCMod::SetNDOCore(int natoms, int nbasis, 
	               double* fmat,double* gss, double* gsd, double* gdd, int* ian,
                   int* ifst_bf_at, int* ilst_bf_at,
				   double* core_ch, int* iat_bf,double* fa)
{
    HaMat_double ene_d(3,MAX_AT_TYPE);    // (I+A)/2 for CNDO/INDO/ZINDO parameter setup
    HaMat_double beta0(3,MAX_AT_TYPE); // off-diagonal interaction coef for NDO methods 
	int i;

	if(p_ndo_pars_db == NULL) p_ndo_pars_db = new HaMatDB();
	if(!p_ndo_pars_db->is_open())
	{
		std::string db_name = pApp->harlem_home_dir + "basis/" + "ndo_pars_1.db";
		FILE* fdb = fopen(db_name.c_str(),"r");
		int fexist = FALSE;	
		if(fdb != NULL) 
		{
			fexist = TRUE;
		    fclose(fdb);
		}
		int ires = 0;
		if(fexist)
		{
		    ires = p_ndo_pars_db->open(db_name.c_str(),"r");
//			p_ndo_pars_db->ListKeys();
		}
		if(!ires)
		{
			PrintLog("Unable to open DB file %s \n",db_name.c_str());
		}
	}

	if( gf_integ.num_rows() == 0)
	{
		gf_integ.newsize(MAX_AT_TYPE,9);
		gf_integ = 0.0;
	}

	if(rr_integ.num_rows() == 0)
	{
		rr_integ.newsize(MAX_AT_TYPE,3);
		rr_integ = 0.0;
		int elem;
		for( elem = 21; elem <= 48; elem++)
		{
		   int nn,n1, n2, n3;
           HaVec_double cf_sp(2),exp_sp(2),cf_d(2),exp_d(2);

           NDOExp(elem,0, FALSE, nn, cf_sp, exp_sp);
		   n1 = nc[elem-1];
		   n2 = n1;
		   n3 = n1 - 1;
		   double es = exp_sp[0];
		   double ep = exp_sp[0];
           NDOExp(elem,2, FALSE, nn, cf_d,  exp_d);
		   double ed = exp_d[0];
		   
		   if(ed < 1.0e-6) continue;

           double r1sppd = RRIntSTO(1,n1,es,n2,ep,n2,ep,n3,ed);
		   double r2sddd = RRIntSTO(2,n1,es,n3,ed,n3,ed,n3,ed);
		   double r2sdpp = RRIntSTO(2,n1,es,n2,ep,n3,ed,n2,ep);

		   rr_integ.r0(elem-1,0) = r1sppd;
           rr_integ.r0(elem-1,1) = r2sddd;
		   rr_integ.r0(elem-1,2) = r2sdpp;
		}
	}

	if(p_ndo_pars_db == NULL) p_ndo_pars_db = new HaMatDB();
	if(p_ndo_pars_db->is_open())
	{
		p_ndo_pars_db->GetMat("IP_DNM1_S1",   ip_dnm1_s1);
		p_ndo_pars_db->GetMat("IP_DNM2_S2",   ip_dnm2_s2);
		p_ndo_pars_db->GetMat("IP_DN_S0",     ip_dn_s0);
		p_ndo_pars_db->GetMat("IPEA_DNM1_S1", ipea_dnm1_s1);
		p_ndo_pars_db->GetMat("IPEA_DNM2_S2", ipea_dnm2_s2);
        p_ndo_pars_db->GetMat("IP_COEF",      ip_coef);
		p_ndo_pars_db->GetMat("G_F_INTEG",    gf_integ);
		p_ndo_pars_db->GetMat("BETA_ZINDO_1", beta_zindo_1);
		p_ndo_pars_db->GetMat("BETA_ZINDO_S", beta_zindo_s);
        p_ndo_pars_db->close();
	}

	if( gf_integ.num_rows() == 0)
	{
		gf_integ.newsize(MAX_AT_TYPE,12);
		gf_integ = 0.0;
	}

	for(i = 0; i < MAX_AT_TYPE; i++)
	{
		ene_d.r0(0,i) = ie_ss_1[i]/HARTREE_TO_EV;
		ene_d.r0(1,i) = ie_pp_1[i]/HARTREE_TO_EV;
		ene_d.r0(2,i) = ie_dd_1[i]/HARTREE_TO_EV;
        beta0.r0(0,i) = beta0_sp_1[i]/HARTREE_TO_EV;
		beta0.r0(1,i) = beta0_sp_1[i]/HARTREE_TO_EV;
		beta0.r0(2,i) = beta0_dd_1[i]/HARTREE_TO_EV;
	}


	if(ndo_method == harlem::qc::CNDO_2)
	{
        gf_integ = 0.0;
	}

 	for(i = 0; i < 18; i++) // to check CNDO/2 still working
 	{
        beta0.SetVal_idx0(0,i, beta0_val[i]/HARTREE_TO_EV);
		beta0.SetVal_idx0(1,i, beta0_val[i]/HARTREE_TO_EV);
	}

	if(ndo_method == harlem::qc::ZINDO_1 || ndo_method == harlem::qc::ZINDO_S)
	{
		int n1 = ip_dnm1_s1.num_rows();
		int n2 = ip_dnm2_s2.num_rows();
		int n3 = ip_coef.num_rows();

		int nn = MinFun(n1,n2);
		nn = MinFun(nn,n3);
	
		for(i = 0; i < nn; i++)
		{
			double cf1 = ip_coef.r0(i,1);
			double cf2 = ip_coef.r0(i,2);
            double cf = cf1 + cf2;
            double e1 = cf1 * ip_dnm1_s1.r0(i,1);
			double e2 = cf2 * ip_dnm2_s2.r0(i,1);
			ene_d.r0(0,i) = ((e1+e2)/HARTREE_TO_EV)/cf;

			e1 = cf1 * ip_dnm1_s1.r0(i,2);
			e2 = cf2 * ip_dnm2_s2.r0(i,2);
			ene_d.r0(1,i) = ((e1+e2)/HARTREE_TO_EV)/cf;

			e1 = cf1 * ip_dnm1_s1.r0(i,3);
			e2 = cf2 * ip_dnm2_s2.r0(i,3);
			ene_d.r0(2,i) = ((e1+e2)/HARTREE_TO_EV)/cf;

            beta0.r0(0,i) = beta0_sp_1[i]/HARTREE_TO_EV;
		    beta0.r0(1,i) = beta0_sp_1[i]/HARTREE_TO_EV;
 	        beta0.r0(2,i) = beta0_dd_1[i]/HARTREE_TO_EV;
		}
	}
	
	int j;
	if( ndo_method == harlem::qc::ZINDO_1 && beta_zindo_1.num_rows() > 0)
    {
		for(i = 0; i < beta_zindo_1.num_rows(); i++)
		{
			for(j = 1; j < 3; j++)
			{
               beta0.r0(j-1,i) = beta_zindo_1.r0(i,j)/HARTREE_TO_EV;
			}
		}
	}
	if( ndo_method == harlem::qc::ZINDO_S && beta_zindo_s.num_rows() > 0 )
	{
		for(i = 0; i < beta_zindo_s.num_rows(); i++)
		{
			for(j = 1; j < 3; j++)
			{
               beta0.r0(j-1,i) = beta_zindo_s.r0(i,j)/HARTREE_TO_EV;
			}
		}
	}
	
	if(ndo_method == harlem::qc::CNDO_2 || ndo_method == harlem::qc::INDO_2)
	{
		int n1 = ipea_dnm1_s1.num_rows();
		int n2 = ipea_dnm2_s2.num_rows();
		int n3 = ip_coef.num_rows();

		int nn = MinFun(n1,n2);
		nn = MinFun(nn,n3);
	
		for(i = 0; i < nn; i++)
		{
 			double cf1 = ip_coef.r0(i,1);
			double cf2 = ip_coef.r0(i,2);
            double cf = cf1 + cf2;
             
            double e1 = cf1 * ipea_dnm1_s1.r0(i,1);
 			double e2 = cf2 * ipea_dnm2_s2.r0(i,1);
			ene_d.r0(0,i) =  ((e1+e2)/HARTREE_TO_EV)/cf;

			e1 = cf1 * ipea_dnm1_s1.r0(i,2);
			e2 = cf2 * ipea_dnm2_s2.r0(i,2);
			ene_d.r0(1,i) = ((e1+e2)/HARTREE_TO_EV)/cf;

			e1 = cf1 * ipea_dnm1_s1.r0(i,3);
			e2 = cf2 * ipea_dnm2_s2.r0(i,3);
			ene_d.r0(2,i) = ((e1+e2)/HARTREE_TO_EV)/cf;

//          beta0.SetVal_idx0(0,i, beta0_sp_1[i]/HARTREE_TO_EV);
//		    beta0.SetVal_idx0(1,i, beta0_sp_1[i]/HARTREE_TO_EV);
//	        beta0.SetVal_idx0(2,i, beta0_dd_1[i]/HARTREE_TO_EV);
		}
	}

	if(ndo_method == harlem::qc::INDO_2 || ndo_method == harlem::qc::ZINDO_1 ||  ndo_method == harlem::qc::ZINDO_S )
	{
 		if(gf_integ.num_rows() == 0)
		{
			gf_integ.newsize(MAX_AT_TYPE,9);
			gf_integ = 0.0;
            memcpy(&gf_integ.r0(0,1),g1sp_val,18*sizeof(double));
            memcpy(&gf_integ.r0(0,2),f2pp_val,18*sizeof(double));
		}
	}

	int ia,ja,iadr;
    for(i= 0; i < nbasis; i++)
	{
        ia = iat_bf[i] - 1;
        iadr = LInd0(i,i);
		int lli = 0;
        if(ian[ia] > 0)
		{
		   lli = ll_bf[i- ifst_bf_at[ia]+1]; // Angular momentum of orbital i
           double e0 = ene_d.r0(lli, ian[ia]-1);
		   double gamma_s = gss[LInd0(ia,ia)];
		   double gamma_p = gss[LInd0(ia,ia)];
           double gamma_sp = gss[LInd0(ia,ia)];
           double gamma_pd = gsd[ia + natoms*ia];
		   double gamma_sd = gsd[ia + natoms*ia];
		   double gamma_d = gdd[LInd0(ia,ia)];
 
           fa[iadr]   = e0;
 		   if(ndo_method == harlem::qc::CNDO_2)
		   {
			  if( lli < 2)  fmat[iadr] = e0 + (0.5 - core_ch[ia]) * gamma_s;
		      else fmat[iadr] = e0 + (0.5 - core_ch[ia]) * gamma_d;
//			  if( lli < 2)  fmat[iadr] = e0 + 0.5 * gamma_s;
//		      else fmat[iadr] = e0 + 0.5 * gamma_d;
		   }
           if(ndo_method == harlem::qc::INDO_2 || ndo_method == harlem::qc::ZINDO_1 ||  ndo_method == harlem::qc::ZINDO_S ) 
           {
               int elem = ian[ia];
			   double g1sp = gf_integ.r0(elem-1,1);
               double f2pp = gf_integ.r0(elem-1,2);
               double g2sd = gf_integ.r0(elem-1,3);
			   double g1pd = gf_integ.r0(elem-1,4);
               double f2pd = gf_integ.r0(elem-1,5);
               double g3pd = gf_integ.r0(elem-1,6);
               double f2dd = gf_integ.r0(elem-1,7);
               double f4dd = gf_integ.r0(elem-1,8);

			   double ns,np,nd,nf;
			   double nval = GetNDOValEl(elem,ns,np,nd,nf);
			   if( lli == 0)
			   {
                   fmat[iadr] = e0 - ((ns - 1.0)*gamma_s + np*gamma_sp + nd*gamma_sd);
                   fmat[iadr] += np*g1sp/6.0 + nd*g2sd/10.0;
			   }
			   else if(lli == 1)
			   {
				   if( np < 1.0E-6 && elem != 20)
				   {
                       np = 1.0;
					   ns -= 1.0;
				   }
                   fmat[iadr] = e0 - (ns*gamma_sp + (np-1.0)*gamma_p + nd*gamma_pd);
                   fmat[iadr] += (np-1.0)*2.0*f2pp/25.0 + ns*g1sp/6.0 + 
					              nd*( g1pd/15.0 + 3.0*g3pd/70.0);
			   }
			   else if( lli == 2)
			   {
                   fmat[iadr] = e0 - (ns*gamma_sd + np*gamma_pd + (nd-1.0)*gamma_d);
                   fmat[iadr] += (nd-1.0)*2.0*(f2dd+f4dd)/63.0 + ns*g2sd/10.0 + 
					             np*(g1pd/15.0 + 3.0*g3pd/70.0);
			   }         
		   }
		}

        for( ja = 0 ; ja <  natoms; ja++)
		{
		   if( ja == ia) continue;
		   int elem = ian[ja]; 
		   
		   double nval,ns,np,nd,nf;
		   if( elem > 0)
		   {
			   nval = GetNDOValEl(elem,ns,np,nd,nf);
           }
		   else
		   {
				ns = core_ch[ja];
				np = 0.0; nd = 0.0; nf = 0.0;
		   }

		   double gamma_sp = gss[LInd0(ia,ja)];
          

		   if( lli < 2) 
		   {
		       double gamma_sd = gsd[ia + natoms*ja];
               fmat[iadr] -= ((ns+np)*gamma_sp + nd*gamma_sd)  ;  
//			   fmat[iadr] -= core_ch[ja]*gss[LInd0(ja,ia)]; // f_ii = f_ii - sum_jat( ch_core(jat)*G(iat,jat) )
		   }
		   else
		   {
			   double gamma_sd = gsd[ja + natoms*ia];
			   double gamma_dd = gdd[LInd0(ia,ja)];
			   fmat[iadr] -= ((ns+np)*gamma_sd + nd*gamma_dd)  ;  
//               fmat[iadr] -= core_ch[ja]*gsd[ja+ ia*natoms]; // f_ii = f_ii - sum_jat( ch_core(jat)*G(iat,jat) )
		   }
		}
	}
//    Off diagonal terms -- beta0(ia,ja)*s(i,j).
//
    for( i= 1; i < nbasis; i++)
	{
		ia = iat_bf[i] - 1;
		int lli = ll_bf[i- ifst_bf_at[ia]+1];
        int iatyp = ian[ia] - 1;
		for( j = 0; j < i; j++)
		{
		  ja = iat_bf[j] - 1;
		  int llj = ll_bf[j- ifst_bf_at[ja]+1];
          int jatyp = ian[ja] - 1;
          iadr = LInd0(i,j);
		  fmat[iadr] *= 0.5 * ( beta0.GetVal_idx0(lli,iatyp) + beta0.GetVal_idx0(llj,jatyp));
//		  if( iatyp > 9 || jatyp > 9) fmat[iadr] *= 0.75;
		  fa[iadr] = fmat[iadr];
		}
	}
	return TRUE;
}


int HaQCMod::SetINDOAtCoulMat(int natoms, HaVec_double& cm, int ia, int elem, double* gss, double* gsd, double* gdd)
{
   int i;
   int j;

   if(elem < 1) return TRUE;

   int nt = 1;             // only S-orbitals
   if( elem > 2) nt = 2;    // S and P orbitals
   if( elem > 20) nt = 3;    // S,P and D orbitals

   if( nt == 1)
   {
	  cm.newsize(idx4(1,1,1,1) + 1);
   }
   if( nt == 2)
   {
	  cm.newsize(idx4(4,4,4,4) + 1);
   }
   if( nt == 3)
   {
	  cm.newsize(idx4(9,9,9,9) + 1);
   }
   cm = 0.0;

   double fp258 = 0.2581989;
   double fp218 = 0.2182179;      
   double sqrt3 = sqrt(3.0);

   double gamma_ss = gss[LInd0(ia,ia)];
   double gamma_sp = gss[LInd0(ia,ia)];
   double gamma_pp = gss[LInd0(ia,ia)];
   double gamma_sd = gsd[ia + natoms*ia];
   double gamma_pd = gsd[ia + natoms*ia];
   double gamma_dd = gdd[LInd0(ia,ia)];
   
   			   
   double g1sp = gf_integ.r0(elem-1,1)/3.0;
   double f2pp = gf_integ.r0(elem-1,2)/25.0;
   double g2sd = gf_integ.r0(elem-1,3)/5.0;
   double g1pd = gf_integ.r0(elem-1,4)/15.0;
   double f2pd = gf_integ.r0(elem-1,5)/35.0;
   double g3pd = gf_integ.r0(elem-1,6)/245.0;
   double f2dd = gf_integ.r0(elem-1,7)/49.0;
   double f4dd = gf_integ.r0(elem-1,8)/441.0;
   
   double r1sppd = rr_integ.r0(elem-1,0);  
   double r2sddd = rr_integ.r0(elem-1,1);
   double r2sdpp = rr_integ.r0(elem-1,2);
   
   i = 1; j = 1;
   cm[ idx4(i,j,i,j) ] = gamma_ss;    //  (ss|ss)
                                          
   if( nt > 1)
   {
      cm[ idx4(i,j,i+1,j+1) ] = gamma_sp;   //   (ss|pxpx)                                                 
      cm[ idx4(i,j,i+2,j+2) ] = gamma_sp;   //   (ss|pypy)
	  cm[ idx4(i,j,i+3,j+3) ] = gamma_sp;   //   (ss|pzpz)
   }

   if( nt > 2 )
   {
      cm[ idx4(i,j,i+4,j+4) ] = gamma_sd;   //   (ss|dz2 dz2)                                                 
      cm[ idx4(i,j,i+5,j+5) ] = gamma_sd;   //   (ss|dx2-y2 dx2-y2)                                                 
      cm[ idx4(i,j,i+6,j+6) ] = gamma_sd;   //   (ss|dxy dxy)                                                 
      cm[ idx4(i,j,i+7,j+7) ] = gamma_sd;   //   (ss|dxz dxz)                                                 
      cm[ idx4(i,j,i+8,j+8) ] = gamma_sd;   //   (ss|dyz dyz)                                                 
   }

   if( nt == 1) return TRUE;

   i = 1; j = 2;
   cm[ idx4(i,j,i,j) ] = g1sp;    //  (s px|s px)

   if( nt > 2)
   {
      cm[ idx4s(i,j,i+1,j+3) ] = -0.1490712*r1sppd;      //  (s px|px dz2)
      cm[ idx4s(i,j,i+1,j+4) ] =  fp258*r1sppd;          //  (s px|px dx2-y2)
      cm[ idx4s(i,j,i+2,j+5) ] =  fp258*r1sppd;          //  (s px|py dxy)
      cm[ idx4s(i,j,i+3,j+6) ] =  fp258*r1sppd;          //  (s px|pz dxz)
   }

   i = 2; j = 2;

   cm[ idx4s(i,j,i,j) ]     = gamma_pp + 4.0*f2pp;  //  (px px|px px)
   cm[ idx4s(i,j,i+1,j+1) ] = gamma_pp - 2.0*f2pp;  //  (px px|py py)
   cm[ idx4s(i,j,i+2,j+2) ] = gamma_pp - 2.0*f2pp;  //  (px px|pz pz)
   
   if( nt > 2)
   {
      cm[ idx4s(i,j,i+3,j+3) ]     = gamma_pd - 2.0*f2pd;  //  (px px|dz2 dz2)
	  cm[ idx4s(i,j,i+3,j+4) ]     = -3.4641016*f2pd;      //  (px px|dx2-y2 dz2)
	  cm[ idx4s(i,j,i+4,j+4) ]     = gamma_pd + 2.0*f2pd;  //  (px px|dx2-y2 dx2-y2)
	  cm[ idx4s(i,j,i+5,j+5) ]     = gamma_pd + 2.0*f2pd;  //  (px px|dxy dxy)
	  cm[ idx4s(i,j,i+6,j+6) ]     = gamma_pd + 2.0*f2pd;  //  (px px|dxz dxz)
	  cm[ idx4s(i,j,i+7,j+7) ]     = gamma_pd - 4.0*f2pd;  //  (px px|dyz dyz)
   }

   i = 1; j = 3;
   
   cm[ idx4(i,j,i,j) ]     = g1sp;  //  (s py|s py)
   
   if( nt > 2)
   {
       cm[ idx4s(i,j,i+1,j+4) ]     = fp258*r1sppd;       //  (s py|px dxy)
	   cm[ idx4s(i,j,i+2,j+2) ]     = -0.1490712*r1sppd;  //  (s py|py dz2)
       cm[ idx4s(i,j,i+2,j+3) ]     = -fp258*r1sppd;      //  (s py|py dx2-y2)
	   cm[ idx4s(i,j,i+3,j+6) ]     = fp258*r1sppd;      //  (s py|pz dyz)

   }
                                       
   i = 2; j = 3;
   
   cm[ idx4(i,j,i,j) ]     = 3.0*f2pp;  //  (px py|px py)

   if(nt > 2)
   {
      cm[ idx4s(i,j,i+3,j+4) ]     = -3.4641016*f2pd;  //  (px py|dz2 dxy)
	  cm[ idx4s(i,j,i+6,j+6) ]     = 3.0*f2pd;         //  (px py|dxz dyz)
   }

   i = 3; j = 3;

   cm[ idx4(i,j,i,j) ]     = gamma_pp + 4.0*f2pp;  //  (py py|py py)
   cm[ idx4(i,j,i+1,j+1) ] = gamma_pp - 2.0*f2pp;  //  (py py|pz pz)
                             
   if( nt > 2)
   {
      cm[ idx4s(i,j,i+2,j+2) ]     = gamma_pd - 2.0*f2pd;  //  (py py|dz2 dz2)
	  cm[ idx4s(i,j,i+2,j+3) ]     = 3.4641016*f2pd;       //  (py py|dz2 dx2-y2)
	  cm[ idx4s(i,j,i+3,j+3) ]     = gamma_pd + 2.0*f2pd;  //  (py py|dx2-y2 dx2-y2)
	  cm[ idx4s(i,j,i+4,j+4) ]     = gamma_pd + 2.0*f2pd;  //  (py py|dxy dxy)
	  cm[ idx4s(i,j,i+5,j+5) ]     = gamma_pd - 4.0*f2pd;  //  (py py|dxz dxz)
	  cm[ idx4s(i,j,i+6,j+6) ]     = gamma_pd + 2.0*f2pd;  //  (py py|dyz dyz)
   }
                                                                                           
   i = 1; j = 4;

   cm[ idx4(i,j,i,j) ]     = g1sp;  //  (s pz|s pz)

   if(nt > 2)
   {
	  cm[ idx4s(i,j,i+1,j+4) ]   = fp258*r1sppd;  //  (s pz|px dxz)
	  cm[ idx4s(i,j,i+2,j+5) ]   = fp258*r1sppd;  //  (s pz|py dyz)
	  cm[ idx4s(i,j,i+3,j+1) ]   = 0.2981424*r1sppd;  //  (s pz|pz dz2)
   }
                                          
   i = 2; j = 4;

   cm[ idx4(i,j,i,j) ]     = 3.0*f2pp;  //  (px pz|px pz)                                      
                             
   if(nt > 2)
   {
      cm[ idx4s(i,j,i+3,j+4) ]     = sqrt3*f2pd;  //  (px pz|dz2 dxz) 
	  cm[ idx4s(i,j,i+4,j+4) ]     = 3.0*f2pd;    //  (px pz|dx2-y2 dxz) 
	  cm[ idx4s(i,j,i+5,j+5) ]     = 3.0*f2pd;    //  (px pz|dxy dxz) 
   }
 
   i = 3; j = 4;
  
   cm[ idx4(i,j,i,j) ]     = 3.0*f2pp;  //  (py pz|py pz)          

   if(nt > 2)
   {
	  cm[ idx4s(i,j,i+2,j+5) ]     = sqrt3*f2pd;  //  (py pz|dz2 dxz) 
	  cm[ idx4s(i,j,i+3,j+5) ]     = -3.0*f2pd;   //  (py pz|dx2-y2 dxz) 
	  cm[ idx4s(i,j,i+4,j+4) ]     =  3.0*f2pd;   //  (py pz|dxy dxz) 
   }

   i = 4; j = 4;

   cm[ idx4(i,j,i,j) ]     = gamma_pp + 4.0*f2pp;  //  (pz pz|pz pz)         

   if(nt > 2)
   {
      cm[ idx4s(i,j,i+1,j+1) ]     = gamma_pd + 4.0*f2pd;  //  (pz pz|dz2 dz2)  
	  cm[ idx4s(i,j,i+2,j+2) ]     = gamma_pd - 4.0*f2pd;  //  (pz pz|dx2-y2 dx2-y2)    
	  cm[ idx4s(i,j,i+3,j+3) ]     = gamma_pd - 4.0*f2pd;  //  (pz pz|dxy dxy)    
	  cm[ idx4s(i,j,i+4,j+4) ]     = gamma_pd + 2.0*f2pd;  //  (pz pz|dxz dxz)    
	  cm[ idx4s(i,j,i+5,j+5) ]     = gamma_pd + 2.0*f2pd;  //  (pz pz|dyz dyz)    
   }
                                            
   if(nt == 2) return TRUE;

   i = 1; j = 5;
 
   cm[ idx4s(i,j,i+1,j-3) ] = -0.0894427*r2sdpp;  //  (s dz2|px px) 
   cm[ idx4s(i,j,i+2,j-2) ] = -0.0894427*r2sdpp;  //  (s dz2|py py)  
   cm[ idx4s(i,j,i+3,j-1) ] =  0.1788854*r2sdpp;  //  (s dz2|pz pz)   
   cm[ idx4s(i,j,i,j) ]     =  g2sd;              //  (s dz2|s dz2)   
   cm[ idx4s(i,j,i+4,j) ]   =  0.1277753*r2sddd;  //  (s dz2|dz2 dz2)   
   cm[ idx4s(i,j,i+5,j+1) ] =  -0.1277753*r2sddd; //  (s dz2|dx2-y2 dx2-y2) 
   cm[ idx4s(i,j,i+6,j+2) ] =  -0.1277753*r2sddd; //  (s dz2|dxy dxy)
   cm[ idx4s(i,j,i+7,j+3) ] =  0.0638877*r2sddd;  //  (s dz2|dxz dxz)
   cm[ idx4s(i,j,i+8,j+4) ] =  0.0638877*r2sddd;  //  (s dz2|dyz dyz)
   
   i = 2; j = 5;
       
   cm[ idx4s(i,j,i,j) ]     =  g1pd + 18.0*g3pd;               //  (px dz2|px dz2) 
   cm[ idx4s(i,j,i,j+1) ]   =  -sqrt3*g1pd - 5.1961524*g3pd;   //  (px dz2|px dx2-y2) 
   cm[ idx4s(i,j,i+1,j+2) ] =  -sqrt3*g1pd - 5.1961524*g3pd;   //  (px dz2|py dxy) 
   cm[ idx4s(i,j,i+2,j+3) ] =  -sqrt3*g1pd + 20.7846097*g3pd;  //  (px dz2|pz dxz) 
                                     
   i = 3; j = 5;
   
   cm[ idx4s(i,j,i,j) ]     =  g1pd + 18.0*g3pd;     //  (py dz2|py dz2) 
   cm[ idx4s(i,j,i,j+1) ]   =  sqrt3*g1pd + 5.1961524*g3pd;   //  (py dz2|py dx2-y2) 
   cm[ idx4s(i,j,i+1,j+4) ] =  -sqrt3*g1pd + 20.7846097*g3pd;   //  (py dz2|pz dyz)                                                       
     
   i = 4; j = 5;
                
   cm[ idx4s(i,j,i,j) ]     =  4.0*g1pd + 27.0*g3pd;     //  (pz dz2|pz dz2) 

   i = 5; j = 5;

   cm[ idx4s(i,j,i,j) ]     =  gamma_dd + 4.0*f2dd + 36.0*f4dd;     //  (dz2 dz2|dz2 dz2) 
   cm[ idx4s(i,j,i+1,j+1) ] =  gamma_dd - 4.0*f2dd +  6.0*f4dd;     //  (dz2 dz2|dx2-y2 dx2-y2) 
   cm[ idx4s(i,j,i+2,j+2) ] =  gamma_dd - 4.0*f2dd +  6.0*f4dd;     //  (dz2 dz2|dxy dxy) 
   cm[ idx4s(i,j,i+3,j+3) ] =  gamma_dd + 2.0*f2dd - 24.0*f4dd;     //  (dz2 dz2|dxz dxz)
   cm[ idx4s(i,j,i+4,j+4) ] =  gamma_dd + 2.0*f2dd - 24.0*f4dd;     //  (dz2 dz2|dyz dyz) 

   i = 1; j = 6;

   cm[ idx4s(i,j,i+1,j-4) ]     = 0.1549193*r2sdpp ;      //  (s dx2-y2|px px) 
   cm[ idx4s(i,j,i+2,j-3) ]     = -0.1549193*r2sdpp ;     //  (s dx2-y2|py py) 
   cm[ idx4s(i,j,i,j) ]         = g2sd;                   //  (s dx2-y2|s dx2-y2) 
   cm[ idx4s(i,j,i+4,j) ]       = -0.1277753*r2sddd ;     //  (s dx2-y2|dz2 dx2-y2)
   cm[ idx4s(i,j,i+7,j+2) ]     = 0.1106567*r2sddd ;      //  (s dx2-y2|dxz dxz)
   cm[ idx4s(i,j,i+8,j+3) ]     = -0.1106567*r2sddd ;     //  (s dx2-y2|dyz dyz)
  
   i = 2; j = 6;

   cm[ idx4s(i,j,i,j) ]         = 3.0*g1pd + 24.0*g3pd;  // (px dx2-y2 | px dx2-y2)
   cm[ idx4s(i,j,i+1,j+1) ]     = 3.0*g1pd - 21.0*g3pd;  // (px dx2-y2 | py dxy)
   cm[ idx4s(i,j,i+2,j+2) ]     = 3.0*g1pd - 6.0*g3pd;   // (px dx2-y2 | pz dxz)

   i = 3; j = 6;

   cm[ idx4s(i,j,i,j) ]         = 3.0*g1pd + 24.0*g3pd;  // (py dx2-y2 | py dx2-y2)
   cm[ idx4s(i,j,i+1,j+3) ]     = -3.0*g1pd - 6.0*g3pd;   // (px dx2-y2 | pz dyz)

   i = 4; j = 6;

   cm[ idx4s(i,j,i,j) ]         = 15.0*g3pd;  // (pz dx2-y2 | pz dx2-y2)

   i = 5; j = 6;

   cm[ idx4s(i,j,i,j)     ] = 4.0*f2dd + 15.0*f4dd;               // (dz2 dx2-y2 | dz2 dx2-y2)
   cm[ idx4s(i,j,i+3,j+2) ] = -3.4641016*f2dd + 17.3205081*f4dd;  // (dz2 dx2-y2 | dxz dxz)
   cm[ idx4s(i,j,i+4,j+3) ] =  3.4641016*f2dd - 17.3205081*f4dd;  // (dz2 dx2-y2 | dyz dyz)

   i = 6; j = 6;

   cm[ idx4s(i,j,i,j)     ] = gamma_dd + 4.0*f2dd + 36.0*f4dd;    // (dx2-y2 dx2-y2 | dx2-y2 dx2-y2)
   cm[ idx4s(i,j,i+1,j+1) ] = gamma_dd + 4.0*f2dd - 34.0*f4dd;    // (dx2-y2 dx2-y2 | dx2-y2 dx2-y2)
   cm[ idx4s(i,j,i+2,j+2) ] = gamma_dd - 2.0*f2dd -  4.0*f4dd;    // (dx2-y2 dx2-y2 | dx2-y2 dx2-y2)
   cm[ idx4s(i,j,i+3,j+3) ] = gamma_dd - 2.0*f2dd -  4.0*f4dd;    // (dx2-y2 dx2-y2 | dx2-y2 dx2-y2)

   i = 1; j = 7;

   cm[ idx4s(i,j,i+1,j-4) ]     = 0.1549193*r2sdpp ;      //  (s dxy|px py) 
   cm[ idx4s(i,j,i,j) ]         = g2sd;                   //  (s dxy|s dxy) 
   cm[ idx4s(i,j,i+4,j) ]       = -0.1277753*r2sddd ;     //  (s dxy|dz2 dxy)
   cm[ idx4s(i,j,i+7,j+2) ]     = 0.1106567*r2sddd ;      //  (s dxy|dxz dyz)

   i = 2; j = 7;

   cm[ idx4s(i,j,i,j) ]      = 3.0*g1pd + 24.0*g3pd;             // (px dxy | px dxy)
   cm[ idx4s(i,j,i+1,j-2) ]  = -1.7320508*g1pd -5.1961524*g3pd;  // (px dxy | py dz2)
   cm[ idx4s(i,j,i+1,j-1) ]  = 3.0*g1pd +21.0*g3pd;              // (px dxy | py dx2-y2)
   cm[ idx4s(i,j,i+2,j+2) ]  = 3.0*g1pd -6.0*g3pd;               // (px dxy | pz dyz)

   i = 3; j = 7;

   cm[ idx4s(i,j,i,j) ]      = 3.0*g1pd + 24.0*g3pd;             // (py dxy | py dxy)
   cm[ idx4s(i,j,i+1,j+1) ]  = 3.0*g1pd -  6.0*g3pd;             // (py dxy | pz dxz)

   i = 4; j = 7;

   cm[ idx4s(i,j,i,j) ]      = 15.0*g3pd;             // (pz dxy | pz dxy)

   i = 5; j = 7;

   cm[ idx4s(i,j,i,j) ]      = 4.0*f2dd + 15.0*f4dd;               // (dz2 dxy | dz2 dxy)
   cm[ idx4s(i,j,i+3,j+2) ]  = -3.4641016*f2dd + 17.3205081*f4dd;  // (dz2 dxy | dxz dyz)

   i = 6; j = 7;

   cm[ idx4s(i,j,i,j) ]      = 35.0*f4dd;                         // (dx2-y2 dxy | dx2-y2 dxy)

   i = 7; j = 7;

   cm[ idx4s(i,j,i,j) ]      = gamma_dd + 4.0*f2dd + 36.0*f4dd;  // (dxy dxy | dxy dxy)
   cm[ idx4s(i,j,i+1,j+1) ]  = gamma_dd - 2.0*f2dd - 4.0*f4dd;   // (dxy dxy | dxz dxz)
   cm[ idx4s(i,j,i+2,j+2) ]  = gamma_dd - 2.0*f2dd - 4.0*f4dd;   // (dxy dxy | dyz dyz)

   i = 1; j = 8;

   cm[ idx4s(i,j,i+1,j-4) ]     = 0.1549193*r2sdpp ;      //  (s dxz|px pz) 
   cm[ idx4s(i,j,i,j) ]         = g2sd;                   //  (s dxz|s dxz) 
   cm[ idx4s(i,j,i+4,j) ]       = 0.0638877*r2sddd ;      //  (s dxz|dz2 dxz) 
   cm[ idx4s(i,j,i+5,j) ]       = 0.1106567*r2sddd ;      //  (s dxz|dx2-y2 dxz) 
   cm[ idx4s(i,j,i+6,j+1) ]     = 0.1106567*r2sddd ;      //  (s dxz|dxy dyz) 

   i = 2; j = 8;

   cm[ idx4s(i,j,i,j) ]      = 3.0*g1pd + 24.0*g3pd;             // (px dxz | px dxz)
   cm[ idx4s(i,j,i+1,j+1) ]  = 3.0*g1pd -  6.0*g3pd;             // (px dxz | py dyz)
   cm[ idx4s(i,j,i+2,j-3) ]  = 3.4641016*g1pd - 15.5884573*g3pd; // (px dxz | pz dz2)
   cm[ idx4s(i,j,i+2,j-2) ]  = 15.0*g3pd;                        // (px dxz | pz dx2-y2)

   i = 3; j = 8;

   cm[ idx4s(i,j,i,j) ]      = 15.0*g3pd;             // (py dxz | py dxz)
   cm[ idx4s(i,j,i+1,j+1) ]  = 15.0*g3pd;             // (py dxz | pz dyz)

   i = 4; j = 8;
   
   cm[ idx4s(i,j,i,j) ]      = 3.0*g1pd + 24.0*g3pd;   // (pz dxz | pz dxz)

   i = 5; j = 8;

   cm[ idx4s(i,j,i,j) ]    = f2dd + 30.0*f4dd;                 // (dz2 dxz | dz2 dxz)
   cm[ idx4s(i,j,i+1,j) ]  = 1.7320508*f2dd - 8.6602540*f4dd;  // (dz2 dxz | dx2-y2 dxz)
   cm[ idx4s(i,j,i+2,j+1) ]= 1.7320508*f2dd - 8.6602540*f4dd;  // (dz2 dxz | dxy dyz)

   i = 6; j = 8;

   cm[ idx4s(i,j,i,j) ]      = 3.0*f2dd + 20.0*f4dd;   // ( dx2-y2 dxz | dx2-y2 dxz)
   cm[ idx4s(i,j,i+1,j+1) ]  = 3.0*f2dd - 15.0*f4dd;   // ( dx2-y2 dxz | dxy dyz)

   i = 7; j = 8;

   cm[ idx4s(i,j,i,j) ]      = 3.0*f2dd + 20.0*f4dd;   // ( dxy dxz | dxy dxz)

   i = 8;  j = 8;

   cm[ idx4s(i,j,i,j) ]      = gamma_dd + 4.0*f2dd + 36.0*f4dd;   // ( dxz dxz | dxz dxz)
   cm[ idx4s(i,j,i+1,j+1) ]  = gamma_dd - 2.0*f2dd -  4.0*f4dd;   // ( dxz dxz | dyz dyz)

   i = 1; j = 9;

   cm[ idx4s(i,j,i+2,j-5) ]     = 0.1549193*r2sdpp ;      //  (s dyz|py pz) 
   cm[ idx4s(i,j,i,j) ]         = g2sd;                   //  (s dyz|s dyz) 
   cm[ idx4s(i,j,i+4,j) ]       = 0.0638877*r2sddd ;      //  (s dyz|dz2 dyz) 
   cm[ idx4s(i,j,i+5,j) ]       = -0.1106567*r2sddd ;     //  (s dyz|dx2-y2 dyz) 
   cm[ idx4s(i,j,i+6,j-1) ]     = 0.1106567*r2sddd ;      //  (s dyz|dxy dxz) 

   i = 2; j = 9;
   
   cm[ idx4s(i,j,i,j) ]      = 15.0*g3pd;                 // (px dyz | px dyz)
   cm[ idx4s(i,j,i+1,j-1) ]  = 15.0*g3pd;                 // (px dxz | py dxz)
   cm[ idx4s(i,j,i+2,j-2) ]  = 15.0*g3pd;                 // (px dxz | pz dxy)

   i = 3; j = 9;

   cm[ idx4s(i,j,i,j) ]      = 3.0*g1pd + 24.0*g3pd;            // (py dyz | py dyz)
   cm[ idx4s(i,j,i+1,j-4) ]  = 3.4641016*g1pd -15.5884573*g3pd; // (py dyz | pz dz2)
   cm[ idx4s(i,j,i+1,j-3) ]  = -15.0*g3pd;                      // (py dyz | pz dx2-y2)

   i = 4; j = 9;

   cm[ idx4s(i,j,i,j) ]      = 3.0*g1pd + 24.0*g3pd;            // (pz dyz | pz dyz)

   i = 5; j = 9;

   cm[ idx4s(i,j,i,j) ]      = f2dd + 30.0*f4dd;                 // (dz2 dyz | dz2 dyz)
   cm[ idx4s(i,j,i+1,j) ]    = -1.7320508*f2dd + 8.6602540*f4dd; // (dz2 dyz | dx2-y2 dyz)
   cm[ idx4s(i,j,i+2,j-1) ]  =  1.7320508*f2dd - 8.6602540*f4dd; // (dz2 dyz | dxy dxz)

   i = 6; j = 9;

   cm[ idx4s(i,j,i,j) ]      = 3.0*f2dd + 20.0*f4dd;             // (dx2-y2 dyz | dx2-y2 dyz)
   cm[ idx4s(i,j,i+1,j-1) ]  = -3.0*f2dd + 15.0*f4dd;            // (dx2-y2 dyz | dxy dxz)

   i = 7; j = 9;

   cm[ idx4s(i,j,i,j) ]      = 3.0*f2dd + 20.0*f4dd;             // (dxy dyz | dxy dyz)

   i = 8; j = 9;

   cm[ idx4s(i,j,i,j) ]      = 3.0*f2dd + 20.0*f4dd;             // (dxz dyz | dxz dyz)

   i = 9; j = 9;

   cm[ idx4s(i,j,i,j) ]      = gamma_dd + 4.0*f2dd + 36.0*f4dd;  // (dyz dyz | dyz dyz)

// Rearrange d-orbitals integrals:

   

   return TRUE;
}

int HaQCMod::RunCNDOThread()
{
	stop_calc_flag = FALSE;

	MolSet* pmset = GetMolSet();
	InitBasis("MINB6G");
	InitLocOrb("FULL_ATOMIC_BASIS");
	wave_fun_type = harlem::qc::NDO;

	PrintLog(" \n\n");
	if( ndo_method == harlem::qc::CNDO_2)  PrintLog(" Run CNDO/2 calculations \n\n");
	if( ndo_method == harlem::qc::INDO_2)  PrintLog(" Run INDO/2 calculations \n\n");
	if( ndo_method == harlem::qc::ZINDO_1) PrintLog(" Run ZINDO/1 calculations \n\n");
	if( ndo_method == harlem::qc::ZINDO_S) PrintLog(" Run ZINDO/S calculations \n\n");

	int icntrl = 0; // do not compute gradient

    int natoms = pmset->GetNAtoms();
	int icharg = GetCharge();
	int multip = GetMult();

	int iout = 6;
	int iprint = 1;

	HaMat_double coord(3,natoms);
	HaVec_int ian(natoms);
	HaVec_int iattyp(natoms);  // if =-1 dummy atom
    HaVec_double core_ch(natoms); // charges of atoms cores and external charges
	AtomIteratorMolSet aitr(pmset);
	HaAtom* aptr;
	int i = 0;
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		if(aptr->IsDummy()) 
		{
			ian[i] = 0;
            iattyp[i] = -1;
			core_ch[i] = aptr->GetCharge();
		}
		else
		{
			iattyp[i] = 1;
            ian[i] = aptr->GetElemNo();
		}
		coord.SetVal_idx0(0,i, aptr->GetX_Bohr());
		coord.SetVal_idx0(1,i, aptr->GetY_Bohr());
		coord.SetVal_idx0(2,i, aptr->GetZ_Bohr());
		i++;
	}

	int iter  = 0; // out: number of iterations?
	double energy;   // out: the computed energy

	int nb = GetNBfunc();

// parameters
    int igrad = FALSE; 
    int usepsd = TRUE;

    int iuhf1 = iuhf + 1;
    int iextp = 0;
    if(max_it_noavg != 0) iextp = 2;
	int icount = 0;
       
//  Index of fist and last basis functions of the atom
    HaVec_int   ilst_bf_at(natoms); // index of the first basis function of the atom
    HaVec_int   ifst_bf_at(natoms);  // index of the last basis function of the atom

// Atom index of the basis function
    HaVec_int iat_bf(12000); // IIU  

    int convgd = FALSE;

	int nbzdo;  // Number of basis functions
    int naezdo; // Number of alpha electrons
	int nbezdo; // Number of beta electrons
    
	int iat;
	for(iat = 0; iat < natoms; iat++)
	{
		ifst_bf_at[iat]= 0;
		ilst_bf_at[iat]= -1;
	}

	int natb = AtBasis.GetNumAtBas();
	int ibs;
	int idx = 1;

	AtomIntMap atoms_idx = pmset->GetAtomSeqNumMap();

	for(ibs = 0; ibs < natb; ibs++)
	{
		GauAtomBasis& atb = AtBasis.GetAtBasByIdx(ibs);
		iat = atoms_idx[atb.GetAtHost()];
 		core_ch[iat] = atb.GetNumElectr();
		ifst_bf_at[iat] = idx;
		idx += atb.GetNBfunc();
        ilst_bf_at[iat] = idx - 1;
		int k;
		for(k= ifst_bf_at[iat]; k <= ilst_bf_at[iat]; k++)
		{
			iat_bf[k-1] = iat+1;
		}
	}

	nbzdo = AtBasis.GetNBfunc();
    int nae_act = GetNumAlphaEl();
	int nbe_act = GetNumBetaEl();

	naezdo = nae_act;
	nbezdo = nbe_act;

	PrintLog( " N alpha electrons = %d N beta electrons = %d \n", nae_act, nbe_act);
//	PrintLog( " naezdo  = %d nbezdo = %d \n", naezdo, nbezdo);

    int navzdo = nbzdo - nae_act; // Number of vacant alpha basis functions
    int nbvzdo = nbzdo - nbe_act; // Number of vacant beta basis functions
      
	HaVec_int itmp(nbzdo);
	for(int ii = 0; ii < nbzdo; ii++)
	{
		itmp[ii] = iat_bf[ii];
	}
	iat_bf.newsize(nbzdo);
	iat_bf = itmp;

	int nttzdo = (nbzdo*(nbzdo+1))/2;  
      
//    Allocate MO, gradient, and overlap arrays first.

	HaMat_double cmo_a(nbzdo,nbzdo); // IMOA
    HaMat_double cmo_b;              // IMOB
	if(iuhf) cmo_b.newsize(nbzdo,nbzdo);

    HaVec_double ene_mo_a(nbzdo);      // IEigA
	HaVec_double ene_mo_b;
    if(iuhf) ene_mo_b.newsize(nbzdo);  // IEigB

    int dopsda = usepsd && (naezdo > 0) && (naezdo < nbzdo);
    int dopsdb = usepsd && (nbezdo > 0) && (nbezdo < nbzdo);
      
	HaVec_double gss(natoms*(natoms+1));  // Gamma integrals between S and P orbitals
	HaMat_double gsd(natoms,natoms);  // Gamma integrals between S or P and D orbitals
	HaVec_double gdd(natoms*(natoms+1));  // Gamma integrals between D orbitals

    HaVec_double ss(nbzdo*(nbzdo+1));  //  Overlap Matrix and One-e(core) hamiltonian
	HaVec_double ss_scaled(nbzdo*(nbzdo+1));  // Overlap Matrix scaled by sigma/pi scaling factor (for CNNDO/S, ZINDO/S)

//    Integral scratch arrays.

   vector<HaVec_double> at_coul_int(natoms);

   HaVec_double f2(18);
   HaVec_double g1(18);

   HaVec_double ff(nbzdo*(nbzdo-1));  // IF - Fock matrix

//     Density matrices.

   HaMat_double pa(nbzdo,nbzdo); // IPA - alpha density matrix 
   HaMat_double pb;              // IPB - beta denstiy matrix
   double* pb_ptr = NULL;
   if(iuhf) 
   {
	   pb.newsize(nbzdo,nbzdo);
	   pb_ptr = pb.begin();
   }
   else
   {
       pb_ptr = pa.begin();
   }
   
   HaMat_double pas(nbzdo,nbzdo); // IPAS - save alpha density matrix 
   HaMat_double pbs;              // IPBS - save beta density matrix 
   if(iuhf) pbs.newsize(nbzdo,nbzdo);

//     Scratch space.
   size_t len_scr = nbzdo;
   if(usepsd) len_scr = MaxFun(nbzdo,nbzdo*naezdo);
	  
//      HaVec_double scr1(len_scr);  // ISCR1

      HaVec_double scr1(nbzdo*nbzdo*iuhf1);  // ISCR1

	  int lscr2 = 2;
	  HaVec_double scr2(lscr2*nbzdo);  // ISCR2

	  HaVec_double deltp; // IDeltp - previous density matricies for interpolation (3 matricies saved)
      deltp.newsize(3*iuhf1*nttzdo);

	  double zero = 0.0;
	  int izero = 0;

      CNDOInteg(iprint, izero,natoms,ian,coord,
				 ss.v(), ss_scaled.v(), gss.v(), gsd.v(), gdd.v(),
				 ilst_bf_at, ifst_bf_at);
		
	  HaMat_double smat(nbzdo,nbzdo);

	  int j;

//	  PrintLog("\n Overlap Matrix \n");
//	  HaMat_double::PrintSymmMat(cout,ss.begin(),nbzdo);

//	  PrintLog("\n Scaled Overlap Matrix \n");
//	  HaMat_double::PrintSymmMat(cout,ss_scaled.begin(),nbzdo);
      
	  HaVec_double* ph1 = &ss;
	  if( ndo_method == harlem::qc::ZINDO_S)
	  {
		 ph1 = &ss_scaled;
	  }


	  HaVec_double h1save, gss_save, gdd_save;
	  HaMat_double gsd_save;
	  if(igrad)
	  {
	        h1save = (*ph1);
	        gss_save = gss;
	        gsd_save = gsd;
	        gdd_save = gdd;
	  }

//     Convert the overlap matrix into the core Hamiltonian and
//    generate an initial guess.

      logical tr = 1;
	  logical fl = 0;
      
	  int calc_init_guess = TRUE; 

	  MO_coef.newsize(0,0);

	  if( set_guess_from_mos )
	  {
		  if(MO_coef.num_rows() != nbzdo)
		  {
			  PrintLog(" Dimensions of MO_coef do not equal to the number of Basis Functions \n");
			  PrintLog(" Initial guess will be computed \n");
		  }
		  else
		  {
			  PrintLog(" Setting initial guess from MO coef \n");
              calc_init_guess = FALSE;
			  cmo_a = MO_coef;
		  }
	  }

	  SetNDOCore(natoms,nbzdo,ph1->v(),gss.v(),gsd.v(),gdd.v(),ian.v(),
                 ifst_bf_at.v(),ilst_bf_at.v(),core_ch.v(),iat_bf.v(),ff.v());

//	  PrintLog("\n G(SP-SP) Coulomb Matrix \n");
//	  HaMat_double::PrintSymmMat(cout,gss.begin(),natoms);

//	  PrintLog("\n G(SP-D) Coulomb Matrix \n");
//	  gsd.Print_format(cout,"%10.5f");
	  
//	  PrintLog("\n G(D-D) Coulomb Matrix \n");
//	  HaMat_double::PrintSymmMat(cout,gdd.begin(),natoms);
	      
//	  PrintLog("\n Core Hamiltonian: \n");
//	  HaMat_double::PrintSymmMat(cout,ph1->begin(),nbzdo);

      int ia;
	  for( ia = 0; ia < natoms; ia++) // Set INDO electron repulsion matrix for atom orbitals
	  {
	     SetINDOAtCoulMat(natoms, at_coul_int[ia], ia, ian[ia], gss.begin(), gsd.begin(), gdd.begin());
	  }

      double nuc_rep_ene = CalcNucRepEne(natoms,ian,coord,core_ch,gss.begin());

      if(calc_init_guess )
	  {
		  for(i = 0; i < nbzdo; i++) // scale by a half diagonal elements for a better initial guess
		  {
              ff[i*(i+1)/2 + i] *= 0.5;
		  }
//	      PrintLog("\n 1e-Hamiltonian matrix for Initial Guess: \n");
//		  HaMat_double::PrintSymmMat(cout,ff.begin(),nbzdo);
			
		  HaMat_double smat(nbzdo,nbzdo);
		  mat_symm_from_lin(smat,ff,nbzdo);
		  HaMat_double::mat_sdiag(smat,cmo_a,ene_mo_a);

//		  PrintLog(" Eigen values of initial guess: \n");
//		  for(i = 0; i < nbzdo; i++)
//		  {
//			  PrintLog("%d %12.6f \n",i+1,ene_mo_a[i]);
//		  }
	  }
	  
      if(iuhf) cmo_b = cmo_a;
 
	  energy = zero;

      deltp = 0.0;

//     SCF loop.  Given MOs, form a density matrix and Fock matrix.
//    Evaluate the energy and test convergence.  Form the next set
//    of MOs.

      int niter = max_it_noavg;
      if(niter == 0) niter = max_scf_iter;
      niter = niter + max_it_avg;

      if(guess_only) 
	  {
		  PrintLog("Calculate Guess Only \n");
		  niter = 0;
	  }

	  int iflag = 0;
	  if(iter_temp == 0) iter_temp = 20;
	  double temp_fermi = 0.0;

	  for(iter = 1; iter < niter; iter++)
	  {
         energy = zero;
		 if(temp0_fermi > 1.0e-5 && iter < iter_temp) 
		 {
			temp_fermi = temp0_fermi*(iter_temp - iter)/iter_temp;
		 }
		 else
		 {
            temp_fermi = 0.0;
		 }

		 FormDenMat(cmo_a,pa.begin(),naezdo,ene_mo_a.begin(),temp_fermi);
		 if(iuhf)FormDenMat(cmo_b,pb.begin(),nbezdo,ene_mo_b.begin(),temp_fermi);

		 double rmsdp = 0.0;
		 
//		 iflag = 0;
         extrdm_(&iout,&iprint,&icount,&nbzdo,&iuhf1,&nttzdo,&iflag,&iextp,
                 &rmsdp,pa.begin(),scr1.begin(),deltp.begin());
		
		if((iter > 1) && (iter <= max_it_avg))
		{
			HaMat_double scr3;
			mat_scale(pa,pa,0.5);
			mat_scale(scr3,pas,0.5);
			mat_add(pa,pa,scr3);
			if(iuhf)
			{
				mat_scale(pb,pb,0.5);
				mat_scale(scr3,pbs,0.5);
				mat_add(pb,pb,scr3);
			}
		}
        if(iter < max_it_avg)
		{
		  pas = pa;
          if(iuhf) pbs = pb;
		}

//		mat_scale(pa,pa,0.5);
//		if(iuhf) mat_scale(pb,pb,0.5);

		double* p_pb = pa.begin();
        if(iuhf) p_pb = pa.begin();

		ff = *ph1;
        FormFockNDO(natoms, ian.begin(), ifst_bf_at.begin(), ilst_bf_at.begin(),
					pa.begin(), p_pb, ff.begin(), 
					gss.begin(),gsd.begin(),gdd.begin(), at_coul_int);
		
 
		int ij = 0;
		energy = 0.0;
		for( i= 0; i < nbzdo; i++)
		{
			for(j = 0; j < i; j++)
			{
		       energy += pa[ij]*( (*ph1)[ij] + ff[ij]);   
			   ij++;
			}
            energy += 0.5*pa[ij]*((*ph1)[ij] + ff[ij]);	
			ij++;
		}

//	    PrintLog("\n Fock Matrix \n");
//		HaMat_double::PrintSymmMat(cout,ff.begin(),nbzdo);

// Don't do pseudodiagonalization:

        dopsda = 0;
		dopsdb = 0;

        if(dopsda && (iter >= 3) &&  iter % 20 != 1)
		{
//           psudag_(&tr,&fl,&nbzdo,&naezdo,&navzdo,&zero,ff.begin(),
//                   cmo_a.begin(),ene_mo_a.begin(),scr1.begin(),&anorm);
		}
        else
		{
			HaMat_double smat(nbzdo,nbzdo);
			mat_symm_from_lin(smat,ff,nbzdo);
			HaMat_double::mat_sdiag(smat,cmo_a,ene_mo_a);
		}
        if(iuhf == 0)
		{
            energy = energy + energy;
        }
		else
		{ 
           FormFockNDO(natoms, ian.begin(), ifst_bf_at.begin(), ilst_bf_at.begin(),
					pb.begin(), pa.begin(), ff.begin(), 
					gss.begin(),gsd.begin(),gdd.begin(), at_coul_int);
          
	       if(dopsdb && (iter >= 3) && (iter % 20) != 1) 
		   {
//              psudag_(&tr,&fl,&nbzdo,&nbezdo,&nbvzdo,&zero,ff.begin(),
//                      cmo_b.begin(),ene_mo_b.begin(),scr1.begin(),&anorm);
		   }
		   else
		   { 
			  HaMat_double smat(nbzdo,nbzdo);
			  mat_symm_from_lin(smat,ff,nbzdo);
			  HaMat_double::mat_sdiag(smat,cmo_b,ene_mo_b);
		   }
        }
        if(iprint > 0)
		{
             PrintLog(" iter= %4d el_ene = %14.8f au, tot_ene= %14.8f rmsdp = %14.8f \n",
			       	 iter, energy, (energy + nuc_rep_ene), rmsdp);
		}
        convgd = rmsdp < conv_dm && (iter > max_it_avg) && (iter > 1) && iflag != 1 && iflag != 2;
        
		if(guess_only || stop_calc_flag || convgd) break;
      } // end cycle on iter

	  if(!convgd && !guess_only) 
	  {
		  PrintLog(" Failure to converge in RunCNDO() \n");
	  }
     
	  int imethod = ndo_method;

	  if(iprint > 0) PrintLog(" NDO ene = %12.6f au(Hartree) \n", (energy + nuc_rep_ene));

//     Optionally, compute the gradient by numerically differentiating the
//     integrals.

      if(igrad)
	  {
		  formp_ha_(&tr,&nbzdo,&nbzdo,&naezdo,cmo_a.begin(),pa.begin());
          if(iuhf) formp_ha_(&tr,&nbzdo,&nbzdo,&nbezdo,cmo_b.begin(),pb.begin());

		  double delta[2]= {0.0001,-0.0002};
		  double epm[2];

          HaMat_double cnew = coord;

	      HaMat_double force(3,natoms);  // IFX
          if(igrad != 0) force.newsize(3,natoms);
	
	      for( i = 0; i < natoms*3; i++)
		  {
             int iatom = i/3 + 1;
			 int istep;
             for( istep= 0; istep < 2; istep++)
			 {
				cnew[i] += delta[istep];
                epm[istep] = 0.0;

				(*ph1) = h1save;
				gss = gss_save;
				gsd = gsd_save;
				gdd = gdd_save;
                
                SetNDOCore(natoms,nbzdo,ph1->v(), gss.v(), gsd.v(), gdd.v(), ian.v(),
                           ifst_bf_at.v(), ilst_bf_at.v(),core_ch.v(), iat_bf.v(), ff.v());
								                
                FormFockNDO(natoms, ian.begin(), ifst_bf_at.begin(), ilst_bf_at.begin(),
         					pa.begin(), pb_ptr, ff.begin(), 
		         			gss.begin(),gsd.begin(),gdd.begin(), at_coul_int);
                
				if(iuhf == 0)
				{
                   epm[istep] = epm[istep] + epm[istep];
				}
                else
				{
                    FormFockNDO(natoms, ian.begin(), ifst_bf_at.begin(), ilst_bf_at.begin(),
					pb.begin(), pa.begin(), ff.begin(), 
					gss.begin(),gsd.begin(),gdd.begin(), at_coul_int);
				} 
			 }
             cnew[i] = cnew[i] + delta[0]; // restore original coordinates
             force[i] = (epm[0]-epm[1])/(2.0*delta[0]);
		 }
         zdofnm_(&natoms,core_ch.begin(),coord.begin(),force.begin());
	  }
	
	MO_coef = cmo_a;
	MOene = ene_mo_a;

	return TRUE;
}

int 
HaQCMod::FormFockNDO(int natoms, int* ian,const int* ifst_bf_at, const int* ilst_bf_at,
					 double* da, double* db, double* fm, 
					 const double* gss, const double* gsd,const double* gdd,
					 vector<HaVec_double>& at_coul_int)
//! ian[i] atom types if < 1 - point charge 
{
	int ia,ja;
	double ftr = 0.0;
	int itr = 32;
    int ipr = 0;
	for( ia= 0; ia < natoms; ia++)
	{
        int elem = ian[ia];
		if( elem < 1) continue;
		HaVec_double &ci = at_coul_int[ia];
		int norb = ilst_bf_at[ia] - ifst_bf_at[ia] + 1;
		int ioff = ifst_bf_at[ia] - 1;
		if( norb < 1) continue;
		int i,j,k,l;

		for(i= 1; i <= norb; i++)
		{
			for( j = i; j <= norb; j++)
			{
				for( k = i; k <= norb; k++)
				{
					for(l = k; l <= norb ; l++)
					{		
						double v = ci[idx4(i,j,k,l)]; // Coulomb integral (ij|kl)   
						if( fabs(v) < 1.0e-6) continue;
						
						int ig = ioff + i; // global indexes
						int jg = ioff + j;
						int kg = ioff + k;
						int lg = ioff + l;

						if(ipr) PrintLog(" V(%2d,%2d|%2d,%2d)= %12.6f \n",ig,jg,kg,lg,v);
						
						if( i == j)  // (ii|kl)
						{
							if( k == l) // (ii|kk)
							{
								if( i == k)
								{  
									int in = LInd(ig,ig);
									if(ipr)PrintLog(" X D( %4d)  %12.6f \n", in+1,db[in]);
									fm[in] += db[in]*v;   // F(i,i) += (D(i,i) - DA(i,i)*(ii|ii)                           
								    
								}
								else
								{
									int in = LInd(ig,ig);
									int kn = LInd(kg,kg);
									if(ipr)PrintLog(" X D( %4d)  %12.6f \n", kn+1,(da[kn]+db[kn]));
									fm[in] += (da[kn]+db[kn])*v;  // F(i,i) += D(k,k)*(ii|kk)     
									

									if(ipr)PrintLog(" X D( %4d)  %12.6f \n", in+1,(da[in]+db[in]));
									fm[kn] += (da[in]+db[in])*v;  // F(k,k) += D(i,i)*(ii|kk)
									
									in = LInd(kg,ig);
									if(ipr)PrintLog(" X D( %4d)  %12.6f \n", in+1,-da[in]);
									fm[in] -= da[in]*v;  // F(k,i) -= DA(k,i)*(ii|kk) 
								}                        
							}
							else // (ii|kl)
							{
								int in = LInd(ig,ig);
								int kn = LInd(lg,kg);
								if(ipr)PrintLog(" X D( %4d)  %12.6f \n", in+1, 2.0*(da[kn]+db[kn]));
								fm[in] += 2.0*(da[kn]+db[kn])*v;  // F(i,i) += 2.0*D(l,k)*(ii|kl)
								if(ipr)PrintLog(" X D( %4d)  %12.6f \n", kn+1, (da[in]+db[in]));
								fm[kn] += (da[in]+db[in])*v;      // F(l,k) += D(i,i)*(ii|kl)
								in = LInd(kg,ig);
								kn = LInd(lg,ig);
								if(ipr)PrintLog(" X D( %4d)  %12.6f \n", kn+1,-da[kn]);
								fm[in] -= da[kn]*v;      // F(k,i) -= DA(l,i)* (ii|kl)
								if(ipr)PrintLog(" X D( %4d)  %12.6f \n", in+1,-da[in]);
								fm[kn] -= da[in]*v;      // F(l,i) -= DA(k,i)* (ii|kl)  
								if( i == k) // (ii|il)
								{
									if(ipr)PrintLog(" X D( %4d)  %12.6f \n", kn+1,-da[kn]);
									fm[in] -= da[kn]*v;   // F(i,i) -= DA(l,i)*(ii|il)                             
								}                   
							}
						}
						else  // i != j
						{
							if( k == l)
							{
								int in = LInd(ig,jg);                                                       
								int kn = LInd(kg,kg);
								fm[in] += (da[kn]+db[kn])*v;       // F(i,j) += D(k,k)*(ij|kk)   
								if(ipr)PrintLog(" X D( %4d)  %12.6f \n", in+1,(da[kn]+db[kn]));
								fm[kn] += 2.0*(da[in] + db[in])*v; // F(k,k) += 2*D(k,k)*(ij|kk)
								if(ipr)PrintLog(" X D( %4d)  %12.6f \n", kn+1,2.0*(da[in]+db[in]));
								in = LInd(kg,ig);
								kn = LInd(kg,jg);
								fm[in] -= da[kn]*v;  // F(k,i) -= DA(j,k)* (ij|kk)
								if(ipr)PrintLog(" X D( %4d)  %12.6f \n", in+1, -da[kn]);
								fm[kn] -= da[in]*v;  // F(k,j) -= DA(k,i)* (ij|kk)  
								if(ipr)PrintLog(" X D( %4d)  %12.6f \n", kn+1, -da[in]);
								if( k == j )
								{
									fm[kn] -= da[in]*v;  // F(k,k) -= DA(k,i)* (ik|kk)  
									if(ipr)PrintLog(" X D( %4d)  %12.6f \n", kn+1, -da[in]);
								}
							}
							else
							{
								int in = LInd(ig,jg);
								int kn = LInd(lg,kg);
								fm[in] += 2.0*(da[kn]+ db[kn])*v;    // F(i,j) += 2*D(k,l)*(ij|kl)  
								if(ipr)PrintLog(" X D( %4d)  %12.6f \n", in+1, 2.0*(da[kn]+db[kn]));
								if(in != kn)
								{ 
									fm[kn] += 2.0*(da[in] + db[in]) *v; // F(k,l) += 2*D(i,j)*(ij|kl) 
									if(ipr)PrintLog(" X D( %4d)  %12.6f \n", kn+1, 2.0*(da[in]+db[in]));
								}
								in = LInd(kg,ig);
								kn = LInd(lg,jg);
								fm[in] -= da[kn]*v;  // F(k,i) -= DA(l,j)*(ij|kl)     
                                if(ipr)PrintLog(" X D( %4d)  %12.6f \n", in+1, -da[kn]);
								fm[kn] -= da[in]*v;  // F(l,j) -= DA(k,i)*(ij|kl) 
								if(ipr)PrintLog(" X D( %4d)  %12.6f \n", kn+1, -da[in]);
								if( i == k)
								{		  
									if(j != l) 
									{
										fm[in] -= da[kn]*v;  // if( i == k && j != l) F(k,i) -= DA(l,j)*(ij|il)
							 	        if(ipr)PrintLog(" X D( %4d)  %12.6f \n", in+1, -da[kn]);
									}
								}
								else
								{
									if(j == l) 
									{
										fm[kn] -= da[in]*v; // if( i != k && j == l) F(l,j) -= DA(k,i)*(ij|kj)
								        if(ipr)PrintLog(" X D( %4d)  %12.6f \n", kn+1, -da[in]);
									}
								}
								in = LInd(lg,ig);
								kn = LInd(jg,kg);
								fm[in] -= da[kn]*v;  // F(l,i) -= DA(j,k)*(ij|kl) 
								if(ipr)PrintLog(" X D( %4d)  %12.6f \n", in+1, -da[kn]);
								if( in != kn)
								{
									fm[kn] -= da[in]*v; // if( l != j || i != k) F(j,k) -= F(l,i)*(ij|kl)
									if(ipr)PrintLog(" X D( %4d)  %12.6f \n", kn+1, -da[in]);
								}
								if( j == k)
								{
								    fm[kn] -= da[in]*v;    // if( l == j == i == k) F(j,k) -= F(l,i)*(ij|kl) 
								    if(ipr)PrintLog(" X D( %4d)  %12.6f \n", kn+1, -da[in]);
								}
							}														
						}
//                        if( fabs(fm[itr] - ftr) > 0.000001)
//						{
//			                PrintLog(" %3d %3d %3d %3d f[itr] change to = %12.6f \n",ig,jg,kg,lg,fm[itr]);
//							ftr = fm[itr];
//						}
					} // end of cycle on l
				}
		   }
		}
		for(ja = 0; ja < ia; ja++)
		{
			int il = 0;
			for( i = ifst_bf_at[ia]; i <= ilst_bf_at[ia]; i++)
			{
				il++;
				int jl = 0;
				for(j = ifst_bf_at[ja]; j <= ilst_bf_at[ja]; j++)
				{
				   jl++;
				   int in,kn;
                   in = LInd(i,j);
                   
				   double gg;
				   if( il <= 4 && jl <= 4) gg = gss[ LInd0(ia,ja)];
				   else if( il > 4  && jl <= 4) gg = gsd[ ja + natoms*ia];
				   else if( il <= 4 && jl > 4) gg = gsd[ ia + natoms*ja];
				   else gg = gdd[ LInd0(ia,ja)];

				   fm[in] -= da[in]*gg;
				   in = LInd(i,i);
				   kn = LInd(j,j);
				   fm[in] += (da[kn] + db[kn])*gg;
				   fm[kn] += (da[in] + db[in])*gg;
                        
//				   if( fabs(fm[itr] - ftr) > 0.000001)
//				   {
//			           PrintLog(" %3d %3d %3d %3d f[itr] change to = %12.6f \n",i,i,j,j,fm[itr]);
//					   ftr = fm[itr];
//				   }
				}
			}
		}
	}
    return TRUE;
}

double  HaQCMod::CalcNucRepEne(int natoms, const HaVec_int& ian, const HaMat_double& c, 
		                       const HaVec_double& core_ch, const double* gss)
{
	int ia,ja,k;
	double rene = 0.0;
	for(ia = 1; ia < natoms; ia++)
	{
		int elem1 = ian[ia];
		for(ja = 0; ja < ia; ja++)
		{
			int elem2 = ian[ja];
			if( elem1 < 1 && elem2 < 1) continue; // do not compute repulsion between point charges
			if( ndo_method == harlem::qc::ZINDO_S)
			{
                rene += core_ch[ia] * core_ch[ja] * gss[LInd0(ia,ja)];
			}
			else
			{
			    double dist = 0.0;
                for( k = 0; k < 3; k++)
				{
                   dist += (c.r0(k,ia) - c.r0(k,ja))*(c.r0(k,ia) - c.r0(k,ja));
				}
			    dist = sqrt(dist);
			    if(dist < 0.001) dist = 0.001;
                rene += core_ch[ia] * core_ch[ja]/dist;
			}
		}
	}
	return rene;
}


