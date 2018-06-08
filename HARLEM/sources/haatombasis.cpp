/*! \file haatombasis.cpp

   Gaussian Basis Set  object in HARLEM.

   \author Igor Kurnikov

   \date 1997-2003

*/
#include <mpi.h>

#include <float.h>
#include <math.h>

#include <boost/algorithm/string.hpp>

#include "haatombasis.h"
#include "hapseudopot.h"
#include "haqchem.h"
#include "hamolset.h"
#include "haatbasdb.h"
#include "gaufile.h"

#include "wx/wfstream.h"

#ifdef USE_IPACK  // IPACK headers
//#include "memchk.h"
#include "basis.h"

#endif

#include "tinyxml.h"

StrPtrMap HaBasisSet::ovlp_map_cache;

std::string HaBasisSet::GetID(const HaBasisSet* pbas)
{
	if( pbas == NULL) return "";
	
	char buf[20];
	sprintf(buf,"%x",pbas);
	return buf;
}

int HaBasisSet::RemoveCachedMatForBasis(const HaBasisSet* pbas)
{
	if(pbas == NULL) return FALSE;
	StrPtrMap::iterator itr;
	itr = ovlp_map_cache.begin();
	std::string id_bas = HaBasisSet::GetID(pbas);
	for(; itr != ovlp_map_cache.end(); )
	{
		std::string mat_id = (*itr).first;
		int ipos = mat_id.find(id_bas); 
		if( ipos < 0 ) 
		{
			itr++;
			continue;
		}
		
		HaMat_double* pmat = (HaMat_double*) (*itr).second;
		if( pmat != NULL) delete pmat;

  //>mikola 3/13/07
		//itr = ovlp_map_cache.erase(itr);
    ovlp_map_cache.erase(itr);
    itr = ovlp_map_cache.begin();
    //<mikola 3/13/07
	}
	return TRUE;
}

int HaBasisSet::ClearMatCache()
{
	StrPtrMap::iterator itr;
	itr = ovlp_map_cache.begin();
	for(; itr != ovlp_map_cache.end(); itr++)
	{
		HaMat_double* pmat = (HaMat_double*) (*itr).second;
		if( pmat != NULL) delete pmat;
	} 	
	ovlp_map_cache.clear();
	return TRUE;
}

HaMat_double* HaBasisSet::GetCachedOvlpMat(const HaBasisSet* pbas1, const HaBasisSet* pbas2)
{
	if( pbas1 == NULL || pbas2 == NULL) return NULL;

	std::string id_1 = HaBasisSet::GetID(pbas1);
	std::string id_2 = HaBasisSet::GetID(pbas2);

	std::string id_mat = id_1 + id_2;

	int ifound = ovlp_map_cache.count(id_mat);
	
	if( ifound == 0) return NULL;

	return (HaMat_double*) ovlp_map_cache[id_mat];
}

int
HaBasisSet::SaveInCacheOverlapMap(const HaBasisSet* pbas1, const HaBasisSet* pbas2, const HaMat_double& smat)
{
	if( pbas1 == NULL || pbas2 == NULL)
	{
		PrintLog("\nError in HaBasisSet::SaveInCacheOverlapMap() \n");
		PrintLog("NULL basis set pointer \n");
		return FALSE;
	}

	if( smat.num_rows() != pbas1->GetNBfunc() || smat.num_cols() != pbas2->GetNBfunc())
	{
		PrintLog("\nError in HaBasisSet::SaveInCacheOverlapMap() \n");
		PrintLog("The dimensions of the overlap matrix do not match the number of functions in basis sets \n");
		return FALSE;
	}

	std::string id_1 = HaBasisSet::GetID(pbas1);
	std::string id_2 = HaBasisSet::GetID(pbas2);
	std::string id_mat = id_1 + id_2;

	int ifound = ovlp_map_cache.count(id_mat);

	HaMat_double* p_smat;

	if( ifound == 0) 
	{
		p_smat = new HaMat_double();
		ovlp_map_cache[id_mat] = (void*) p_smat;
	}
	else
	{
		p_smat = (HaMat_double*) ovlp_map_cache[id_mat];
	}

	(*p_smat) = smat;

	return TRUE;
}



int ArrayOrb3D::ExpandInBas(HaMat_double& coef, HaBasisSet& bset)
{
	HaMat_double ss_bb;
	HaMat_double ss_ba;

    HaBasisSet::CalcOvlpMat(&bset,&bset,ss_bb);
	HaBasisSet::CalcOvlpMat(&bset,this,ss_ba);

	HaMat_double::solv_lin_syst_1(ss_bb,ss_ba);

	coef = ss_ba;
	return TRUE;
}

int ArrayOrb3D::GetTransfMat(HaMat_double& transf_mat, const HaMat_double& rot_mat) 
{
    int nb = GetNBfunc();
    transf_mat.newsize(nb,nb);
	transf_mat = 0.0;
	int i;
	for(i = 0; i < nb; i++)
	{
       transf_mat.r0(i,i) = 1.0;
	}
	return TRUE;
}

HaMat_double ArrayOrb3D::GetTransfMat(const HaMat_double& rot_mat )
{
	HaMat_double tr_mat;
	this->GetTransfMat( tr_mat, rot_mat);
	return tr_mat;
}

int ArrayOrb3D::SaveXML(FILE* file_out, int option) const
{
	if( file_out == NULL ) 
	{
		PrintLog("Error in ArrayOrb3D::SaveXML(FILE* file_out) \n");
		PrintLog(" file_out == NULL \n");
		return FALSE;
	}
	TiXmlDocument doc;
	TiXmlDeclaration * decl = new TiXmlDeclaration( "1.0", "", "" );
	doc.LinkEndChild( decl );

	TiXmlElement* root_element = new TiXmlElement("HARLEM_DATA");
		
	doc.LinkEndChild(root_element);

	this->AddXml( root_element );

	bool bres = doc.SaveFile( file_out );
	
	if( !bres ) return FALSE;

	return TRUE;
}

int ArrayOrb3D::LoadXmlFile(FILE* file_inp, const char* tag_name, int option)
{
	TiXmlDocument doc;
	bool bres = doc.LoadFile(file_inp); 

	const TiXmlElement* root_element;
	const TiXmlElement* data_element;

	if(bres)
	{
		root_element = doc.FirstChildElement();
		if( root_element == NULL )
		{
			PrintLog("Error in ArrayOrb3D::LoadXml() \n");
			PrintLog("No ROOT Element \n");
			return FALSE;
		}
		if( strlen(tag_name) == 0)
		{
			data_element = 	root_element->FirstChildElement();
		}
		else
		{
			data_element = 	root_element->FirstChildElement(tag_name);
		}
		if( data_element == NULL)
		{
			PrintLog("Error in ArrayOrb3D::LoadXml() \n");
			PrintLog("No Data Element with name %s \n",tag_name );
			return FALSE;
		}
		this->LoadXml(data_element);
	}
	else
	{
		PrintLog("Error in ArrayOrb3D::LoadXml() \n");
		PrintLog("Failed to Load File \n");
		return FALSE;
	}
	return TRUE;
}


ArrayOrb3D* ArrayOrb3D::CreateObjectWithType(const char* type)
{
	std::string class_type = type;
	
	ArrayOrb3D* parr = NULL;

	boost::to_upper(class_type);

	if( class_type == "GAUBASISSET")
	{
		parr = new GauBasisSet();
	}
	else if( class_type == "GAUATOMBASIS")
	{
		parr = new GauAtomBasis();
	}
	else if ( class_type == "GAUSHELL" )
	{
		parr = new GauShell();
	}

	return parr;
}

GauBasisSet::GauBasisSet()
{
	pure_fun_flag = 1;
}

GauBasisSet::~GauBasisSet()
{

}

int GauBasisSet::LoadToGaussianBCommon() const
{
#if defined(GAUSSVER) 
	LoadToGaussianBas(b_);
#endif
	return true;
}

int GauBasisSet::LoadToGaussianB2Common() const
{
#if defined(GAUSSVER) 
	LoadToGaussianBas(b2_);
#endif
	return true;
}

void GauBasisSet::Clear()
{
	at_bas_vec.clear();          
	atom_bas_idx.clear();                    
	fst_bas_fun_idx.clear();           
	bf_lbls.clear();  
	bset_name.clear();
}

std::string GauBasisSet::GetName() const
{
	return bset_name;
}

bool GauBasisSet::IsGeneric() const
{
	if( bset_name.empty() || bset_name == "GEN" ) return true;
	return false;
}

void GauBasisSet::SetGeneric()
{
	bset_name = "GEN";
}

int GauBasisSet::LoadToGaussianBas(b_type& gaub) const
//!
//! Note don't really determine determine atom order number 
//! Can be a problem if some atom do not have basis functions and atom numbers are used 
//!
{
	double* pc4= &gaub.c3[0]+MAXSHL;   // array of f orb coefficients (C4 in gaussian)
	int* padf= (int*)(&gaub.c3[0]+2*MAXSHL); // index array of D and F orb coefficients (ShlADF in Gaussian)

    int i;
	for(i=0; i< 3*MAXSHL; i++)
	{
		gaub.exx[i]=0.0; gaub.c1[i]=0.0; gaub.c2[i]=0.0; gaub.c3[i]=0.0;
	}

	for(i=0; i< MAXSHL; i++)
	{
		*(pc4+i)=0.0; *(padf+i)=1;
		gaub.x[i]=0.0; gaub.y[i]=0.0; gaub.z[i]=0.0;
		gaub.jan[i]=0;
		gaub.shella[i]=0; gaub.shelln[i]=0; gaub.shellt[i]=0; gaub.shellc[i]=0;
		gaub.aos[i]=0; gaub.aon[i]=0;
	}

	AtBasisType::const_iterator bsitr;
	int nshl=0;
	int mm_sp=1;
	int mm_df=1;
	gaub.maxtyp=0;
	gaub.nshell=0;

	int nao=1;

	const HaAtom* aptr_old = NULL;
	int iat_count = 0;
	for(bsitr= at_bas_vec.begin(); bsitr !=  at_bas_vec.end(); bsitr++)
	{
		const HaAtom* aptr= (*bsitr).GetAtHost();
		if(aptr != aptr_old)
		{
			iat_count++;
			aptr_old = aptr;
		}
		double x_cnt = aptr->GetX_Bohr();
		double y_cnt = aptr->GetY_Bohr();
		double z_cnt = aptr->GetZ_Bohr();
		
		ShellsType::const_iterator shitr;
		for(shitr=(*bsitr).Shells.begin(); shitr != (*bsitr).Shells.end(); shitr++)
		{
			nshl++;                                           // nshl - index of the current shell
			int ng=(*shitr).GetNumGauss();
			gaub.shelln[nshl-1]= (*shitr).GetNumGauss();
			gaub.shella[nshl-1]= mm_sp;

			gaub.x[nshl-1]= x_cnt;
			gaub.y[nshl-1]= y_cnt;
			gaub.z[nshl-1]= z_cnt;
           
            gaub.jan[nshl-1]= iat_count;
			gaub.aos[nshl-1]=nao;

			int l_sh = (*shitr).GetL();
			switch(l_sh)
			{
			case(0):	    
				gaub.shellt[nshl-1]= 0; gaub.shellc[nshl-1]=2; break;
//			case(SP_SHELL):	    
//				gaub.shellt[nshl-1]= 1; gaub.shellc[nshl-1]=2; break;
//			case(SPD_SHELL):    
//				gaub.shellt[nshl-1]= 2; gaub.shellc[nshl-1]=0; break;
			case(1):		
				gaub.shellt[nshl-1]= 1; gaub.shellc[nshl-1]=1; 
				gaub.aos[nshl-1]-= 1; // Gaussian weirdness starting AO number assuming 
				                    // lower angular momentum shells are here
				break;
			case(2):		
				gaub.shellt[nshl-1]= 2; gaub.shellc[nshl-1]=2; 
				gaub.aos[nshl-1]-= 4;
				break;
			case(3):		
				gaub.shellt[nshl-1]= 3; gaub.shellc[nshl-1]=2; 
				gaub.aos[nshl-1]-= 9;
				break;
//			case(G_SHELL):		gaub.shellt[nshell-1]= 4; gaub.shellc[nshell-1]=2; break;
			default:
				{
				cerr << "Unsupported Shell Momentum value in GauBasisSet::LoadToGaussianBCommon() " << l_sh << endl;
				return false;
				}
			}
			if(gaub.shellt[nshl-1] > gaub.maxtyp) gaub.maxtyp= gaub.shellt[nshl-1];
			
			const HaMat_double& cf_mat= (*shitr).GetCoef();
			if( cf_mat.num_cols() != ng)
			{
				
			    ErrorInMod( "GauBasisSet::LoadToGaussianBCommon()", 
					         " The number of columns in the matrix of shell coefficients not equal to the number of gaussian in the shell");
				return false;
            }
//  
// Set expansion coefficients correspond to normalized primitives
//
	        double pi2 = PI/2.0;
            double pi2_32 = pi2*sqrt(pi2);
   
			double fl = pi2_32;
			for(i = 1; i <= (2*l_sh -1); i += 2)
			{
				fl *= 0.25 * i;
			}

			for( i=1; i<= ng; i++)
			{
				gaub.exx[mm_sp-1]=cf_mat(1,i);
				gaub.c1[mm_sp-1]=0.0; gaub.c2[mm_sp-1]=0.0;
				
				double ss = fl/pow(cf_mat(1,i),(l_sh+1.5));

				if(l_sh == 0)
				{	
					gaub.c1[mm_sp-1]= cf_mat(2,i)/sqrt(ss);
				}

				if( l_sh == 1 )
					gaub.c2[mm_sp-1]= cf_mat(2,i)/sqrt(ss);
				
//				if( (shtype == SP_SHELL) || (shtype == SPD_SHELL) )
//					gaub.c2[mm_sp-1]= cf_mat(3,i);

				mm_sp++;

//				if( (shtype == D_SHELL) || (shtype == SPD_SHELL))
                if( l_sh == 2 || l_sh == 3)
				{
					if(i == 1) (*(padf+nshl-1))=mm_df;     // fill ShlADF - Position of contraction
					                                       // coefficients of the shell in C3

					gaub.c3[mm_df-1]=0.0; (*(pc4+mm_df-1))= 0.0;
					if(l_sh == 2) 
						gaub.c3[mm_df-1]= cf_mat(2,i)/sqrt(ss);
					if(l_sh == 3) 
						(*(pc4+mm_df-1))= cf_mat(2,i)/sqrt(ss);
//					if(shtype == SPD_SHELL) 
//						gaub.c3[mm_df-1]= cf_mat(4,i);
					mm_df++;
				}

			}
			nao+=(*shitr).GetNBfunc();
		}
	}
	gaub.nshell=nshl;
	logical ToRawP=TRUE;
//  logical ToRawP=FALSE;
//	renorm_(&ToRawP); // renormalize coefficients! 
	                  // Convert the contraction coefficients in /B/ to refer to normalized
                      // primitives (ToRawP true) or unnormalized primitives (ToRawP false).
	return true;
}


InternalBasis* GauBasisSet::CreateIPackBas()
//! InternalBasis object can be used to compute integrals using IPACK
{
    BasisSet basis;
	int ncnt = at_bas_vec.size();
	basis.loc.resize(ncnt);
	basis.pbasis.resize(ncnt);
	basis.name.resize(ncnt);

    list<PureBasisSet> pure_bset_list;

    pure_bset_list.clear();

	if(pure_fun_flag) // set Spherical or Cartesian Gaussian functions
		flaglist.list[Spherical] = true;
	else
		flaglist.list[Spherical] = false;

	int nb = this->GetNBfunc();
	HaVec_double fnorm(nb,1.0); // normalization coefficients for IPACK functions 

	int i,j;
	int icur = 0;
	for(i = 0; i < at_bas_vec.size(); i++)
	{
		pure_bset_list.push_back( PureBasisSet());
		PureBasisSet& pure_bset= pure_bset_list.back();
		
		for(j = 0; j < at_bas_vec[i].Shells.size(); j++)
		{
			GauShell& shl = at_bas_vec[i].Shells[j];
			int ig;
			int ng=shl.GetNumGauss();
			
			// Individual PureBasisFunction will be deleted when  PureBasisSet will be deleted
			PureBasisFunction* ptr_pure_bfun = new PureBasisFunction(shl.GetL()); 
			for( ig = 1; ig <= ng; ig++ )
			{
				ptr_pure_bfun->add(shl.GetExp(ig), shl.GetCoef(ig) );
			}
			pure_bset.add(ptr_pure_bfun);
			
            if( shl.GetL() == 2)
			{
				fnorm[icur] = 1.0/sqrt(12.0);
                fnorm[icur+3] = 1.0/2.0;
			}
			icur += shl.GetNBfunc();
		}
		HaAtom* aptr = at_bas_vec[i].GetAtHost();
		String str;
		Location at_loc(aptr->GetX_Bohr(), aptr->GetY_Bohr(), aptr->GetZ_Bohr());
		basis.add(&pure_bset, str, at_loc);
	}

	Location dum_loc;
	OrbInfo  orbinfo;

	InternalBasis* pibas = new InternalBasis(basis, dum_loc, ALL, orbinfo);
	
	if(pibas != NULL)
	{
		int nof1 = pibas->no_function();
		pibas->norm.reset(nof1); 
		for(int i=0; i< nof1;i++)
			pibas->norm[i] = fnorm[i];
	}

	return pibas;	
}

GauShell::GauShell()
{
   spherical = TRUE;
}

GauShell::GauShell(int new_l_ang, const int NGauss)
{
  if(new_l_ang < 0 || new_l_ang > 4)
  {
	  ErrorInMod("GauShell::GauShell()", 
		         "Invalid angular momentum value, set to zero");
	  l_ang = 0;
  }
  else
  {
	  l_ang = new_l_ang;
  }

  spherical = TRUE;
  SetNumGauss(NGauss);
}

GauShell::GauShell(int new_l_ang, const HaMat_double& new_cf_mat)
{
  if(new_l_ang < 0 || new_l_ang > 4)
  {
	  ErrorInMod("GauShell::GauShell()", 
		         "Invalid angular momentum value, set to zero");
	  l_ang = 0;
  }
  else
  {
	  l_ang = new_l_ang;
  }	

  spherical = TRUE;
  SetNumGauss(new_cf_mat.num_cols());
  SetCoef(new_cf_mat);
}


int
GauShell::GetNumGauss(void) const
  // set the number of Contracted gaussians 
{
 if(NumGauss > 0 && NumGauss < 100)
	 return(NumGauss);
 else
 {
	 cerr << "error in GauShell::GetNumGauss" << endl;
	 cerr << "invalid number of Elemental Gaussians: " << NumGauss << endl;
	 exit(1);
 }
 return -1;
}

bool
GauShell::SetNumGauss(const int NGauss) 
  // set the number of Contracted gaussians 
{
 if(NGauss > 0 && NGauss < 100)
 {
	 NumGauss=NGauss;
     coef.newsize(2,NGauss);
 }
 else
 {
	 cerr << "error in GauShell::SetNumGauss" << endl;
	 cerr << "invalid number of Elemental Gaussians: " << NGauss << endl;
	 return false;
 }
 return true;
}


int 
GauShell::GetNBfunc() const
{
	if( spherical || l_ang < 2)
		return(2*l_ang + 1);
	return (l_ang+2)*(l_ang+1)/2;
}

int 
GauShell::GetNBfuncCart() const
{
	return (l_ang+2)*(l_ang+1)/2;
}



const char* GauShell::GetShellFunSymbol(int ifun)
{
	int lm = GetL();
	
	if( lm == 0 )
		return "S";
	else if ( lm == 1)
	{
		if( ifun == 0) return "PX";
		else if (ifun == 1) return "PY";
		else if (ifun == 2) return "PZ";
		else
		{
			ErrorInMod("GauShell::GetShellFunSymbol()",
		    	"Invalid index of the function in the shell ");
			return "";
		}
	}
	else if ( lm == 2 )
	{
		if( spherical ) 
		{
			if( ifun == 0) return "D_Z2";  // 2ZZ - XX - YY
			else if (ifun == 1) return "D_XZ";
			else if (ifun == 2) return "D_YZ";
			else if (ifun == 3) return "D_XX_YY"; // XX - YY
			else if (ifun == 4) return "D_XY";
			else
			{
				ErrorInMod("GauShell::GetShellFunSymbol()",
					"Invalid index of the function in the shell ");
				return "";
			}
		}
		else
		{
			if( ifun == 0) return "D_XX";
			else if (ifun == 1) return "D_XY";
			else if (ifun == 2) return "D_XZ";
			else if (ifun == 3) return "D_YY";
			else if (ifun == 4) return "D_YZ";
			else if (ifun == 5) return "D_ZZ";
			else
			{
				ErrorInMod("GauShell::GetShellFunSymbol()",
					"Invalid index of the function in the shell ");
				return "";
			}
		}
    }
	return "G_XXX";
}

std::string GauShell::GetLabel(int idx)
{
	return GetShellFunSymbol(idx);
}

const char* GauShell::GetShellSymbol(void) const
{
 	switch (l_ang)
    {
    case 0:
		return "S";
    case 1:
		return "P";
    case 2:
		return "D";
    case 3:
		return "F";
    case 4:
		return "G";
    default:
		{
			cerr << "GauShell::GetShellSymbol unsupported Angular momentum " << l_ang << endl;
			return "U";
		}
    }
	return "S";
}  


bool
GauShell::DestroyCoef()
{
  NumGauss=0;
  coef.newsize(0,0);
  return true;
}


GauShell::~GauShell()
{
  DestroyCoef();
}

bool
GauShell::SetCoef(const double* NewCoef)
//! set gaussian contraction coefficients array
//! from array NewCoef(2,NumGauss) 
//! which contains sequentially 
//! contraction coefficients of the shell 
//! It is a resposibility of a  calling function to ensure NewCoef is valid
{
  int ii=0;
  int ng= GetNumGauss();
  for(int i=1; i <= ng; i++)
  { 
      coef(1,i)=NewCoef[ii++];
	  coef(2,i)=NewCoef[ii++];
  }
  return true;
}


bool
GauShell::SetCoef(const HaMat_double& new_cf_mat)
{
	if( new_cf_mat.num_cols() != NumGauss || new_cf_mat.num_rows() != 2 )
	{
		cerr << " Error in GauShell::SetCoef() " << endl;
		cerr << " Dimensions of the Coef matrix doesn't correspond to the number of Gaussians in the shell " << endl;
		return false;
	}
	coef= new_cf_mat;
	return true;
}

int 
GauShell::GetTransfMat(HaMat_double& transf_mat, const HaMat_double& rot_mat)
{
   int lorb = GetL();

   int nb = GetNBfunc();
   transf_mat.newsize(nb,nb);
   transf_mat = 0.0;
   if( nb != (2*lorb + 1)) 
   {
       PrintLog("Error in GauShell::GetTransfMat() \n");
	   PrintLog("number of basis functions not equal 2*l+1 \n");
	   return FALSE;
   }

   if( lorb == 0)
   {
      transf_mat.r0(0,0) = 1.0;
   }
   else if( lorb == 1)
   {
      transf_mat = rot_mat;
	  HaMat_double::mat_transpose(transf_mat);
   }
   else if( lorb == 2 )
   {
// From Ivanic J. and Ruedenberg K. J.Phys.Chem.(1996) v.100, 6342-6347

      double rxx = rot_mat.r0(0,0); double ryx = rot_mat.r0(1,0); double rzx = rot_mat.r0(2,0);
	  double rxy = rot_mat.r0(0,1); double ryy = rot_mat.r0(1,1); double rzy = rot_mat.r0(2,1);
      double rxz = rot_mat.r0(0,2); double ryz = rot_mat.r0(1,2); double rzz = rot_mat.r0(2,2);

      double sqrt3 = sqrt(3.0);

      transf_mat.r0(0,0) = (3*rzz*rzz - 1)/2;
      transf_mat.r0(0,1) = sqrt3 * rxz*rzz;
	  transf_mat.r0(0,2) = sqrt3 * ryz*rzz;
	  transf_mat.r0(0,3) = sqrt3 * (rxz*rxz - ryz*ryz)/2;
	  transf_mat.r0(0,4) = sqrt3 * rxz*ryz;

      transf_mat.r0(1,0) = sqrt3 * rzx * rzz;
      transf_mat.r0(1,1) = rxx * rzz + rxz * rzx;
	  transf_mat.r0(1,2) = ryx * rzz + ryz * rzx;
	  transf_mat.r0(1,3) = rxx * rxz - ryx * ryz;
	  transf_mat.r0(1,4) = rxx * ryz + rxz * ryx;

      transf_mat.r0(2,0) = sqrt3 * rzy * rzz;
      transf_mat.r0(2,1) = rxy * rzz + rxz * rzy;
	  transf_mat.r0(2,2) = ryy * rzz + ryz * rzy;
	  transf_mat.r0(2,3) = rxy * rxz - ryy * ryz;
	  transf_mat.r0(2,4) = rxy * ryz + rxz * ryy;

      transf_mat.r0(3,0) = sqrt3 * (rzx * rzx - rzy * rzy)/2 ;
      transf_mat.r0(3,1) = rxx * rzx - rxy * rzy;
	  transf_mat.r0(3,2) = ryx * rzx - ryy * rzy;
	  transf_mat.r0(3,3) = (rxx * rxx - rxy * rxy - ryx * ryx + ryy * ryy)/2;
	  transf_mat.r0(3,4) = rxx * ryx - rxy * ryy;

      transf_mat.r0(4,0) = sqrt3 * rzx * rzy;
      transf_mat.r0(4,1) = rxx * rzy + rxy * rzx;
	  transf_mat.r0(4,2) = ryx * rzy + ryy * rzx;
	  transf_mat.r0(4,3) = rxx * rxy - ryx * ryy;
	  transf_mat.r0(4,4) = rxx * ryy + rxy * ryx;
   }
   return TRUE;
}

int GauShell::TransferBetweenAtoms(PtrPtrMap& pt_corr_map)
{
	return TRUE; 
}
 
TiXmlElement* GauShell::AddXml(TiXmlElement* parent_element, const char* name, int option) const
{
	TiXmlElement* gau_shell_element;
	
	if( strlen(name) == 0)
	{
		gau_shell_element = new TiXmlElement("GauShell");
	}
	else
	{
		gau_shell_element = new TiXmlElement("name");
		gau_shell_element->SetAttribute("TYPE","GauShell");
	}
		
	parent_element->LinkEndChild(gau_shell_element);

	gau_shell_element->SetAttribute("l_ang", this->l_ang);
	gau_shell_element->SetAttribute("NumGauss", this->NumGauss);
	if( this->spherical)
	{
		gau_shell_element->SetAttribute("spherical",1);
	}
	else
	{
		gau_shell_element->SetAttribute("spherical",0);
	}

	TiXmlElement* mat_element = coef.AddXml(gau_shell_element,"coef");

	return gau_shell_element;
}

int GauShell::LoadXml(const TiXmlElement* xml_element,  int option )
{
	if( xml_element == NULL) return FALSE;

	int ires, itmp;

	ires = xml_element->QueryIntAttribute("l_ang",&itmp);
	
	if( ires == TIXML_SUCCESS )
	{
		this->l_ang = itmp;
	}

	ires = xml_element->QueryIntAttribute("NumGauss",&itmp);

	if(ires == TIXML_SUCCESS )
	{
		this->SetNumGauss(itmp);
	}

	ires = xml_element->QueryIntAttribute("spherical",&itmp);
	if( ires == TIXML_SUCCESS )
	{
		spherical = itmp;
	}

	const TiXmlElement* coef_element = xml_element->FirstChildElement("coef");
	if( coef_element != NULL)
	{
		coef.LoadXml(coef_element);
	}

	return TRUE;
}

void GauShell::SaveGaussianInp(std::ostream& os) const
{
	char buf[80];
	int ng = GetNumGauss();
	sprintf(buf," %4d 1.00 ",ng);
	os << GetShellSymbol() << buf << std::endl;

	for(int i=1; i <= ng; i++)
	{
		sprintf(buf,"%17.10e %17.10e ",coef.r1(1,i),coef.r1(2,i));
		os << buf << std::endl; 
	}
}


bool GauShell::operator == (const GauShell & rhs) const 
//! test if two shells are identical
{
	if(this->l_ang != rhs.l_ang)
	{
		return false;
	}
	else if ( this->GetNumGauss() != rhs.GetNumGauss() )
		return false;
	else  
	{
		int ng = this->GetNumGauss();
		
		for(int i=1; i <= ng; i++)
		{
			for( int j=1; j <= 2; j++)
			{
				if( this->coef(j,i) != rhs.coef(j,i) )
					return false;
			}
		}
		return true;
	}
	return false;
} 

bool 
GauShell::operator < (const GauShell & rhs) const
{
  if(*this == rhs)
    return false;
  else 
    return (this->GetNumGauss() < rhs.GetNumGauss());
} 




bool
GauAtomBasis::operator == (const GauAtomBasis & rhs) const
{
	return(this->BasName == rhs.BasName && this->AtomType == rhs.AtomType);

//	if(lhs.Shells.empty() || rhs.Shells.empty()) 
//		return false;
//	if( lhs.Shells.size() != lhs.Shells.size())
//		return false;
//	ShellsType::const_iterator itrl=lhs.Shells.begin();
//	ShellsType::const_iterator itrr=rhs.Shells.begin();
//	for(; itrl != lhs.Shells.end(); itrl++,itrr++)
//	{
//		if( (*itrl) != (*itrr) )
//			return false;
//	}
//	return true;


} 


bool
GauAtomBasis::operator < (const GauAtomBasis & rhs) const
{
	if(this->BasName != rhs.BasName)
		return( this->BasName < rhs.BasName);
    if(this->AtomType != rhs.AtomType)
		return( this->AtomType < rhs.AtomType);
	return false;
} 


double 
GauShell::EvalLinCombInPoint( double x, double y, double z, const double* cf) const
{
    int num_coef = this->GetNBfunc();
	double r2= x*x+y*y+z*z;
	double r = sqrt(r2);
	
	double fval=0.0;

	double frad= 0.0;
	double rnorm = 0.0; // accumulate normalization coeffcient in 
	// a multiexponential subshell
	
	
	if( l_ang == 0 )
	{
		for( int ic=1; ic <= GetNumGauss(); ic++)
		{
			if( coef(2,ic) <= 1e-8 )
				continue;
			rnorm = PI*sqrt(PI)/( 2.0*coef(1,ic)*sqrt(2.0*coef(1,ic)) );
			frad+= coef(2,ic)*exp(-coef(1,ic)*r2)/sqrt(rnorm);
		}
		
		fval = cf[0]*frad;
	}
	else if( l_ang == 1 )
	{
		const double cf_norm_p = sqrt(6.0/PI);
		for( int ic=1; ic <= GetNumGauss(); ic++)
		{
			if( coef(2,ic) <= 1e-8 )
				continue;
			frad+= coef(1,ic)*coef(2,ic)*exp(-coef(1,ic)*r2);
		}
		
		fval = cf_norm_p*( cf[0]*x + cf[1]*y + cf[2]*z )*frad; 
	}
	else if( l_ang == 2 )
	{
		const double cf_norm_d1= sqrt(30.0/PI);
		const double cf_norm_d2= sqrt(30.0/PI);
		const double cf_norm_d3= sqrt(30.0/PI);
		for( int ic=1; ic <= GetNumGauss(); ic++)
		{
			if( coef(2,ic) <= 1e-8 )
				continue;
			frad+= coef(1,ic)*sqrt(coef(1,ic))*coef(2,ic)*exp(-coef(1,ic)*r2);
		}
		fval =   cf_norm_d1* cf[0]*(r2 - 3.0*z*z);
		fval +=  cf_norm_d2* ( cf[1]*z*x + cf[2]*z*y + cf[4]*x*y);
		fval +=  cf_norm_d3* ( cf[3]*(x*x -y*y) );
		fval *= frad;
	}
	return fval;
}

double 
GauShell::GetExtent(int i, double tol) const
{
	if( tol < DBL_EPSILON) tol = 0.001;
	double exp_min = 1E-10;
	int nb = GetNBfunc();
	int j;
	for(j = 0; j < nb; j++)
	{
		double expv = coef.GetVal_idx0(0,i);
		if( expv < exp_min) exp_min = expv;
	}
	double r = log(tol)/exp_min;
	return r;
}

const double zero = 0.0;
const double f1 = 1.0;
const double f2 = 3.0;
const double f4 = 4.0;
const double small_num = 1e-8;

int
GauShell::Normalize()
{
   double pw = (3 + 2*l_ang)/4.0;
   int ng = GetNumGauss();
   int i,j;
   double ss = 0.0;
   for(i = 1; i <= ng; i++)
   {
	   double ci   = GetCoef(i);
	   double expi = GetExp(i);
	   for(j = 1; j <= ng; j++)
	   {
		   double cj   = GetCoef(j);
		   double expj = GetExp(j);
		   ss += ci*cj*pow(4.0*expi*expj/((expi+expj)*(expi+expj)),pw);
	   }
   }

   if(ss < small_num)
   {
	   ErrorInMod("GauShell::Normalize()",
				       " Unable to normalize shell");
	   return FALSE;
   }
   double ff = 1.0/sqrt(ss);
   for(i = 1; i <= ng; i++)
	   coef(2,i) *= ff;
   	  
   return TRUE;	
}

GauAtomBasis::GauAtomBasis(void)
{
	SetDefaultParams();
}

GauAtomBasis::GauAtomBasis(const std::string& NewBasName, const std::string& NewAtomType )
{
  SetDefaultParams();

  GauAtomBasis* pAtBas = bas_db.Extract(NewBasName.c_str(),NewAtomType.c_str());
  if(pAtBas == NULL) return;

  copy_from(*pAtBas);
}

GauAtomBasis::GauAtomBasis(const std::string& NewBasName, HaAtom* aptr )
{
  SetDefaultParams();
  
  std::string AtomLbl= aptr->GetStdSymbol();
  GauAtomBasis* pAtBas = bas_db.Extract(NewBasName.c_str(),AtomLbl.c_str());
  if(pAtBas == NULL) return;
  copy_from(*pAtBas);
  host_atom = aptr;
}

GauAtomBasis::GauAtomBasis(const GauAtomBasis & ref)
{
	SetDefaultParams();
	copy_from(ref);
}

GauAtomBasis::~GauAtomBasis()
{
	if( internal_atom_flag && host_atom != NULL)
	{
		delete host_atom;
	}
}

void
GauAtomBasis::SetDefaultParams()
{
	host_atom = NULL;
	ppot=NULL;
	internal_atom_flag = FALSE;
}


GauAtomBasis&
GauAtomBasis::copy_from(const GauAtomBasis & ref)
{
	Shells.clear();
	Shells.reserve(ref.Shells.size());

	ShellsType::const_iterator sitr;
	for(sitr = ref.Shells.begin(); sitr != ref.Shells.end(); sitr++)
		Shells.push_back(*sitr);

	BasName=ref.BasName;
	AtomType=ref.AtomType;
	host_atom=ref.host_atom;
	PseudoName= ref.PseudoName;
	ppot =ref.ppot;	
	return(*this);
}

int GauAtomBasis::SetForAtom(const char* bas_name, HaAtom* aptr )
{
	std::string AtomLbl = aptr->GetStdSymbol();
  GauAtomBasis* pAtBas = bas_db.Extract(bas_name, AtomLbl.c_str());
  if(pAtBas == NULL) return FALSE;
  copy_from(*pAtBas);
  host_atom = aptr;
  SetPseudoPotFromName();
  return TRUE;
}


const char* GauAtomBasis::GetBasName() const
{
  return BasName.c_str();
}

bool GauAtomBasis::SetBasName(const std::string & name) 
{
  int len = name.size();
  BasName.resize(len);
  for(int i=0; i < len; i++)
	BasName[i] = toupper(name[i]);
  return true;
}


bool GauAtomBasis::SetAtHost(HaAtom* new_host_atom)
{
 	host_atom=new_host_atom;
	return true;
}

int GauAtomBasis::TransferBetweenAtoms(PtrPtrMap& pt_corr_map)
{
	HaAtom* host_old = GetAtHost();
	HaAtom* host_new = (HaAtom*) pt_corr_map.GetVal(host_old );
	if( host_new == NULL) return FALSE;

	SetAtHost( host_new);

	return TRUE; 
}

TiXmlElement* GauAtomBasis::AddXml(TiXmlElement* parent_element, const char* name, int option) const
{
	if( parent_element == NULL) return NULL;

	TiXmlElement* at_bas_element;
	
	if( strlen(name) == 0 )
	{
		at_bas_element = new TiXmlElement("GauAtomBasis");
	}
	else
	{
		at_bas_element = new TiXmlElement(name);
		at_bas_element->SetAttribute("TYPE", "GauAtomBasis");
	}

	parent_element->LinkEndChild(at_bas_element);

	if(!BasName.empty()) at_bas_element->SetAttribute("BasName",  BasName.c_str());
	if(!AtomType.empty()) at_bas_element->SetAttribute("AtomType", AtomType.c_str());
	if(!PseudoName.empty()) at_bas_element->SetAttribute("PseudoName", PseudoName.c_str());

	if(host_atom != NULL) 
	{
//		at_bas_element->SetAttribute("HOST_ATOM_IDI",(int) host_atom);
		TiXmlElement* atom_element = host_atom->AddXml(at_bas_element,"host_atom");
	}
    
	int nsh = Shells.size();

	if( nsh > 0)
	{
		TiXmlElement* shells_element = new TiXmlElement("Shells");
		shells_element->SetAttribute("size", nsh);
		at_bas_element->LinkEndChild(shells_element);

		int i;
		for( i= 0 ; i < nsh; i++)
		{
			TiXmlElement* shell_element = Shells[i].AddXml(shells_element);
		}
	}

	return at_bas_element;
}

int GauAtomBasis::LoadXml(const TiXmlElement* xml_element,  int option )
{
	if( xml_element == NULL) return FALSE;

	int ires, itmp;
	std::string str;

	const char* sattr;

    sattr = xml_element->Attribute("BasName");
	if(sattr) BasName = sattr;
	sattr   = xml_element->Attribute("AtomType");
	if(sattr) AtomType = sattr;
    sattr = xml_element->Attribute("PseudoName");
	if(sattr) PseudoName = sattr;
	
	int nsh = 0;
    
	const TiXmlElement* shells_elem = xml_element->FirstChildElement("Shells");
	const TiXmlElement* child_elem;

	if( shells_elem != NULL)
	{
		ires = xml_element->QueryIntAttribute("size",&itmp);
		if(ires == TIXML_SUCCESS )
		{
			nsh = itmp;
		}
		else
		{		
			for ( child_elem = shells_elem->FirstChildElement(); child_elem != NULL; child_elem = child_elem->NextSiblingElement()) 
			{
				nsh++;
			}
		}
	}
	
	Shells.resize(nsh);
	int ish = 0;
    for ( child_elem = shells_elem->FirstChildElement(); child_elem != NULL; child_elem = child_elem->NextSiblingElement()) 
	{
		Shells[ish].LoadXml(child_elem);
		ish++;
	}
	
	const TiXmlElement* host_atom_elem = xml_element->FirstChildElement("host_atom");

	if( host_atom_elem != NULL )
	{
		HaAtom* aptr = new HaAtom();
		if( internal_atom_flag && this->host_atom != NULL) 
		{
			delete this->host_atom;
		}
		internal_atom_flag = TRUE;
        this->host_atom = aptr;
		aptr->LoadXml(host_atom_elem);
	}

	return TRUE;
}

const HaAtom* GauAtomBasis::GetAtHost() const
{
	return host_atom;
}

HaAtom* GauAtomBasis::GetAtHost() 
{
	return host_atom;
}

void GauAtomBasis::SaveGaussianInp(std::ostream& os) const
{
	wxString str;
	os << "-" << AtomType.c_str() << " " << std::endl;
	
	ShellsType::const_iterator itr;
	for(itr=Shells.begin(); itr != Shells.end(); itr++)
	{
		(*itr).SaveGaussianInp(os);
	}
}

std::string GauAtomBasis::GetAtomType() const
{
	return AtomType;
}

bool GauAtomBasis::SetAtomType(const std::string& atype)
{
  int len = atype.size();
  AtomType.resize(len);
  for(int i=0; i < len; i++)
	 AtomType[i] = toupper(atype[i]);
  return true;
}

int GauAtomBasis::GetNBfunc() const
{
	ShellsType::const_iterator i;
	int nb=0;

	if(Shells.empty())
		return 0;
	else
	{
		for(i=Shells.begin(); i != Shells.end(); i++)
			nb += (*i).GetNBfunc();
	}
	return nb;
}

int GauAtomBasis::GetNBfuncCart() const
{
	ShellsType::const_iterator i;
	int nb=0;

	if(Shells.empty())
		return 0;
	else
	{
		for(i=Shells.begin(); i != Shells.end(); i++)
			nb += (*i).GetNBfuncCart();
	}
	return nb;
}	
	
std::string GauAtomBasis::GetLabel(int idx)
{
	return "GAUSS_BAS_FUNC";
}


bool GauAtomBasis::AddShell(GauShell & shl) 
{
	Shells.push_back(shl);
	GauShell& slast = Shells.back();
	slast.Normalize();
	return true;
}

void GauAtomBasis::Clear()
{
   BasName="";
   AtomType="";
   PseudoName="";
   host_atom=NULL;
   ppot= NULL;

   Shells.clear();
}

static char bufl;

bool GauAtomBasis::SetFromGaussianInp(istream& is)
{
	std::string str;
	is >> str;
	if(str[0] == '-') 
	{
		str=str.substr(1);
	}
	SetAtomType(str);
	for(;;)
	{
		is >> str;
		if(is.eof()) break;
		if( str.substr(0,4) == "****") break;

		boost::to_upper(str);
		boost::trim(str);
		
		std::string sh_symb = str;

		int nsbsh;

		if(sh_symb == "S" || sh_symb == "P" || sh_symb == "D" || sh_symb == "F" ||  sh_symb == "G") 
			nsbsh = 1;
		else if(sh_symb == "SP")
			nsbsh = 2;
		else if(sh_symb == "SPD")
			nsbsh = 3;
		else
		{
			cerr << " Error in GauAtomBasis::SetFromGaussianInp() " << endl;
			cerr << " Invalid Shell Symbol  " << sh_symb << endl;
			return false;
		}

		int ng, i, j;
		double cf_scale;
		is >> ng;
		is >> cf_scale;
		
		HaMat_double cf(nsbsh+1,ng);
		for(i = 1; i <= ng; i++)
		{
			for(j = 1 ; j <= (nsbsh+1); j++)
			{
				is >> cf(j,i);
			}
		}
		
		int l_sh;
		for( int isub = 1; isub <= nsbsh; isub++)
		{
			HaMat_double cf2(2,ng);
			if( sh_symb == "S") l_sh = 0;
			if( sh_symb == "P") l_sh = 1;
			if( sh_symb == "D") l_sh = 2;
			if( sh_symb == "F") l_sh = 3;
			if( sh_symb == "G") l_sh = 4;
			if( sh_symb == "SP" || sh_symb == "SP")
			{
				l_sh = isub - 1; 	
			}
			for(i = 1; i <= ng; i++)
			{
				cf2(1,i) = cf(1,i);
				cf2(2,i) = cf(isub+1,i);
			}
			GauShell shl(l_sh,cf2);
			AddShell(shl);
		}
	}
	return true;
	
}

void GauAtomBasis::SetPseudoPotName(const std::string& new_pot_name) 
{
	ppot= NULL;
	PseudoName = new_pot_name;
}

void GauAtomBasis::SetPseudoPotPtr(const HaPseudoPot* new_ppot) 
{ 
	ppot= new_ppot; 
}

bool GauAtomBasis::SetPseudoPotFromName() 
//! Set PseudoPotential From DataBase using  Name of the PseudoPotential and the name of the atom
{
	if(PseudoName.empty()) return true;
	if(AtomType.empty())
	{
		cerr << " Error in GauAtomBasis::SetPseudoPotFromName() " << endl;
		cerr << " Pseudopotential name is set but Atom Type name is not " << endl;
		return false;
	}

	ppot= pseudo_db.Extract(PseudoName,AtomType);
	if(ppot == NULL)
	{
		cerr << " Error in GauAtomBasis::SetPseudoPotFromName() " << endl;
		cerr << " Cannot find in DB Pseudopotential " << PseudoName  << " For Atom " 
			 << AtomType << endl;
		return false;
	}
	return true;
}

bool GauAtomBasis::Print_info(ostream & sout, const int level) const
{
	sout << "Atomic Basis info: " << endl;
	sout << "Name = " << BasName << " AtomType= " 
		 << AtomType << endl;
	sout << "Number of shells: " << Shells.size() << endl;
    return true;
}

int GauAtomBasis::GetNumElectr() const
{
	if( host_atom == NULL) return 0;
	int ne = host_atom->GetElemNo();
	if(ppot != NULL)
	{
		ne -= ppot->GetNCoreEl();
	}
	std::string bas_name = GetBasName();
	if( bas_name == "MINB6G" )
	{
		if(  ne > 2 && ne <= 10)
		{
			ne -= 2;
		}
		else if( ne > 10 && ne <= 18)
		{
			ne -= 10;
		}
		else if(  ne > 18 && ne <= 36)
		{
			ne -= 18;
		}
		else if( ne > 36 && ne <= 54)
		{
			ne -= 36;
		}
		else if( ne > 54 && ne <= 86)
		{
			ne -= 54;
		}
		else if( ne > 86 )
		{
			ne -= 86;
		}
	}
	return ne;
}

int GauBasisSet::GetNBfunc() const
{
	vector<GauAtomBasis>::const_iterator bfitr;

	int nb = 0;
	for(bfitr = at_bas_vec.begin(); bfitr != at_bas_vec.end(); bfitr++)
	{
		nb += (*bfitr).GetNBfunc();
	}
	return nb;
}

int GauBasisSet::GetNBfuncCart() const
{
	vector<GauAtomBasis>::const_iterator bfitr;

	int nb = 0;
	for(bfitr = at_bas_vec.begin(); bfitr != at_bas_vec.end(); bfitr++)
	{
		nb += (*bfitr).GetNBfuncCart();
	}
	return nb;
}


int GauBasisSet::InitForAtoms(const char* bname, AtomContainer* pat_cont)
{
	at_bas_vec.clear();
	HaAtom* aptr;

	AtomIteratorGen aitr(pat_cont);
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		this->AddBasisToAtom(bname, aptr);
	}

	bset_name = bname;

	return TRUE;
}

std::string GauBasisSet::GetLabel(int idx)
{
	char buf[256];
	HaAtom* aptr = (HaAtom*) GetHostPt(idx);
	aptr->FillRef(buf);
	std::string lbl = buf;
	lbl += "@";

	if( idx < bf_lbls.size())
	{
		HaAtom* aptr = (HaAtom*) GetHostPt(idx);
		aptr->FillRef(buf);
		lbl += bf_lbls[idx];
	}
	else
	{
		lbl += "GAO";
	}
	return lbl;
}

static StrVec lbls_minb[10];

int GauBasisSet::InitForMolSet(const char* bname_str, HaMolSet* pmset)
{
//	PrintLog(" GauBasisSet::InitForMolSet() pt 1 \n");
	at_bas_vec.clear();
	HaAtom* aptr;
	AtomIteratorMolSet aitr(pmset);
	for(aptr= aitr.GetFirstAtom(); aptr; aptr= aitr.GetNextAtom())
	{
		this->AddBasisToAtom(bname_str, aptr);
	}
	std::string bname = bname_str;
	if( bname == "MINB6G")
	{
		int nb = this->GetNBfunc();
		bf_lbls.resize(nb);
		if( lbls_minb[0].size() == 0)
		{
			lbls_minb[0].resize(1);
			lbls_minb[0][0] = "1S";

			lbls_minb[1].resize(4);
			lbls_minb[1][0] = "2S";
			lbls_minb[1][1] = "2PX"; lbls_minb[1][2] = "2PY"; lbls_minb[1][3] = "2PZ";

			lbls_minb[2].resize(4);
			lbls_minb[2][0] = "3S";
			lbls_minb[2][1] = "3PX"; lbls_minb[2][2] = "3PY"; lbls_minb[2][3] = "3PZ";	

			lbls_minb[3].resize(4);
			lbls_minb[3][0] = "4S";
			lbls_minb[3][1] = "4PX"; lbls_minb[3][2] = "4PY"; lbls_minb[3][3] = "4PZ";	

			lbls_minb[4].resize(9);
			lbls_minb[4][0] = "4S";
			lbls_minb[4][1] = "4PX"; lbls_minb[4][2] = "4PY"; lbls_minb[4][3] = "4PZ";
			lbls_minb[4][4] = "3D0"; lbls_minb[4][5] = "3DP1"; lbls_minb[4][6] = "3DM1"; lbls_minb[4][7] = "3DP2"; lbls_minb[4][8] = "3DM2";

			lbls_minb[5].resize(9);
			lbls_minb[5][0] = "5S";
			lbls_minb[5][1] = "5PX"; lbls_minb[5][2] = "5PY"; lbls_minb[5][3] = "5PZ";
			lbls_minb[5][4] = "4D0"; lbls_minb[5][5] = "4DP1"; lbls_minb[5][6] = "4DM1"; lbls_minb[5][7] = "4DP2"; lbls_minb[5][8] = "4DM2";
		}
		int ig = 0;
		AtBasisType::iterator bitr = at_bas_vec.begin();
		for(; bitr != at_bas_vec.end(); bitr++)
		{
			aptr= (HaAtom*) (*bitr).GetAtHost();

			std::string bname ( (*bitr).GetBasName() );
			std::string atsymb= aptr->GetStdSymbol();

			int nf = (*bitr).GetNBfunc();

			int idx_lbl = 5;
			int elem = aptr->GetElemNo();
			if(elem <= 2) idx_lbl = 0;
			if(elem > 2 && elem <= 10) idx_lbl = 1;
			if(elem > 10 && elem <= 18) idx_lbl = 2;
			if(elem > 18 && elem <= 20) idx_lbl = 3;
			if(elem > 20 && elem <= 36) idx_lbl = 4;
			if(elem > 36 && elem <= 48) idx_lbl = 5;

			int il;
			for(il = 0; il < nf; il++)
			{
				bf_lbls[ig] = lbls_minb[idx_lbl][il];
				ig++;
			}
		}
	}

	bset_name = bname_str;

	return TRUE;
}

GauAtomBasis* GauBasisSet::AddBasisToAtom(const char* bas_name, HaAtom* aptr)
{
	if(aptr->IsDummy()) return NULL;
	GauAtomBasis at_bas;
	int ires = at_bas.SetForAtom(bas_name, aptr );
	if( ires == FALSE) return NULL;
	at_bas_vec.push_back(at_bas); 
	GauAtomBasis* pa_out = &(at_bas_vec.back());
	return pa_out;
}

int GauBasisSet::CalcOvlpMat(GauBasisSet* pbas1, GauBasisSet* pbas2, HaMat_double& ovlp_mat)
{
	if(pbas1 == NULL || pbas2 == NULL) return FALSE;
	
	int nb1 = pbas1->GetNBfunc();
	int nb2 = pbas2->GetNBfunc();

	ovlp_mat.newsize(nb1,nb2);

	int same_bas = TRUE;

	if(pbas1 != pbas2 ) same_bas = FALSE;

	if(HaQCMod::int_engine == QCIntEngineType::INT_ENGINE_GAUSS)
	{
		PrintLog("Compute Overlap Matrix using Gaussian Library \n");
#if GAUSSVER == 98
		int c0=0, c1=1;
		int iprint=0;
		double xx=0;
		HaMat_double coord;
		
		HaQCMod::LoadGauCom(*pbas1);
		if(!same_bas) pbas2->LoadToGaussianB2Common();

		int nbasis = nb1;
		
		int nbas6d = pbas1->GetNBfuncCart();
        int nbas6d_2 = pbas2->GetNBfuncCart();

//		getnb6_(&nbas6d);
		
		int ntt6d=nbas6d*(nbas6d+1)/2;
		int ntt= nbasis*(nbasis+1)/2;
		
		AtomBasIterator abitr;
		
		int natoms = pbas1->at_bas_vec.size();
		coord.newsize(3,natoms);
		
		int i = 0;
		for(abitr = pbas1->at_bas_vec.begin(); abitr != pbas1->at_bas_vec.end(); abitr++)
		{
			i++;
			HaAtom* aptr = (*abitr).GetAtHost();
			if( aptr == NULL) continue;
			coord(1,i) = aptr->GetX_Bohr();
			coord(2,i) = aptr->GetY_Bohr();
			coord(3,i) = aptr->GetZ_Bohr();
		}
		
		HaVec_double scr_val;
		scr_val.newsize( ntt6d*2);
//		long iptr = (long)scr_val.begin();
		HaMat_double valmat;
		valmat.set_ext_alloc(scr_val.begin(),ntt,2); 
		
		int icalc=1;  // Overlap and Kinetic Energy integrals
		if(!same_bas) icalc = 2; // Overlap between two basis sets
		int idgst=10; // return Operator Matricies
		
		double* pwork;
		int mdv=10*nbasis*nbasis;
		if(mdv < HaQCMod::max_gauss_mem) mdv= HaQCMod::max_gauss_mem;
		
		pwork=(double*)malloc(mdv*sizeof(double));
		if(pwork == NULL)
		{
			ErrorInMod(" GauBasisSet::CalcOvlpMat() ",
				"Error to allocate work array for DENBAS");
			return FALSE;
		}
		
		int* piwork=(int*)pwork;
		
		int isdim;
		if(same_bas)
		{
			int ipflag=0;
			logical allowp[50];
			for(i=0; i< 50; i++)
				allowp[i]= 0;
			
			denbas_(&io_.iout,&icalc,&idgst,
				&c0,&c0,&c0,&c0,
				&c0,&c0,&c0,&xx,&xx,
				&c0,&nbasis,&c1,&c0,
				coord.begin(),&c0,&c0,&xx, &natoms,
				&xx,&c0,
				&xx,valmat.v(),&xx,
				&c0,&c1,&c0,&xx,
				&iprint,&ipflag,&allowp[0],&c0,
				piwork,pwork,&mdv);
		}
		else
//		if(0)
		{
			int iprint = 1;
			isdim  = nbas6d;
			if(nbas6d_2 > nbas6d) isdim = nbas6d_2;
			int itrans = 0;
			int ipurd1 = !pbas1->pure_fun_flag;
			int ipurf1 = !pbas1->pure_fun_flag;
			int ipurd2 = !pbas2->pure_fun_flag;
			int ipurf2 = !pbas2->pure_fun_flag;


	    	ovlp_(&io_.iout,&iprint, valmat.begin(),&isdim, &itrans,
		          &ipurd1, &ipurf1, &ipurd2, &ipurf2,
		          pwork, &mdv);
		
		}
		free(pwork);
		{
			if(same_bas)
			{
				HaVec_double fmat_lin(ntt,&valmat(1,1)); 
				mat_symm_from_lin(ovlp_mat,fmat_lin,nbasis);
			}
			else
			{
				HaMat_double axx_ovlp;
				axx_ovlp.set_ext_alloc(&valmat(1,1),isdim,isdim); 
				ovlp_mat.newsize(nb1,nb2);
				
				for(i=1 ; i <= nb1; i++)
				{
					for(int j=1; j <= nb2; j++)
					{
						ovlp_mat(i,j) = axx_ovlp(i,j);	
					}
				}
			}
		}
#else
	PrintLog("Error in GauBasSet::CalcOvlpMat() \n");
	PrintLog("Gaussian Library is not initialized \n");
#endif
	}
	else if(HaQCMod::int_engine == QCIntEngineType::INT_ENGINE_IPACK) 
	{
        HaMat_double* p_save_mat = HaBasisSet::GetCachedOvlpMat(pbas1,pbas2);
		if( p_save_mat != NULL)
		{
			ovlp_mat = (*p_save_mat);
			return TRUE;
		}

//		PrintLog("Compute Overlap Matrix using IPACK Library \n");
    	HaQCMod::InitIPack();
		InternalBasis* ip_bas1 = pbas1->CreateIPackBas();
		InternalBasis* ip_bas2 = ip_bas1;
		if(!same_bas) ip_bas2 = pbas2->CreateIPackBas();

	    ARRAY<Mat> ovlp_arr(1);
        ovlp_arr[0].reset(nb1,nb2);
		ovlp_arr[0].set(0);

        compute_overlap(*ip_bas1, *ip_bas2, ovlp_arr, same_bas);
		normalize(*ip_bas1,*ip_bas2,ovlp_arr,0);
//		cout << ovlp_arr[0];
		ovlp_mat.set_from(ovlp_arr[0]);
		
		delete ip_bas1;
		if(!same_bas) delete ip_bas2;

		HaBasisSet::SaveInCacheOverlapMap(pbas1, pbas2, ovlp_mat);
	}
	return TRUE;
}

int
GauBasisSet::GetNumCnt() const
{
	return at_bas_vec.size();
}

GauAtomBasis& 
GauBasisSet::GetAtBasByIdx(int i)
{
	return at_bas_vec[i];
}

int 
GauBasisSet::RecompFstBasVec()
{
	int nb = at_bas_vec.size();
	
	fst_bas_fun_idx.resize(nb);

	if(nb == 0) return TRUE;

	int i; 
	int idx_f = 0;

	for(i = 0; i < nb; i++)
	{
       fst_bas_fun_idx[i] = idx_f;
	   idx_f+= at_bas_vec[i].GetNBfunc();
	}
	return TRUE;
}

int
GauBasisSet::MatchBasisSet(const GauBasisSet* basis_frag, IntIntMap& frag_bas_fun_map, HaVec_double* bas_pert_vec)
//!
//! frag_bas_fun_map - map of fragment atom basis functions to the current basis set functions  
//! bas_pert_vec - if not NULL - filled with estimated perterbation factors( 0 to 1.0) for fragment basis functions 
{
	frag_bas_fun_map.clear();

	int nb_frag = basis_frag->GetNBfunc();
    int nb = this->GetNBfunc();

	int natb_frag    = basis_frag->at_bas_vec.size(); // Number of atom basis in the basis of the fragment 
	int natb = this->at_bas_vec.size();               // Number of atom basis in the active basis of the main molecule 

	int iat,jat;
	int jfst = 0;

	set<const HaAtom*, less<const HaAtom*> > atoms_mol_non_matched; 
	HaVec_double at_bas_frag_pert(natb_frag);

	if( bas_pert_vec != NULL)
	{
		for(iat = 0; iat < natb; iat++)
		{
			GauAtomBasis& atbas = this->at_bas_vec[iat];
            const HaAtom* host_atom = atbas.GetAtHost();
			atoms_mol_non_matched.insert(host_atom);
		}		
	}

	for(jat=0; jat < natb_frag; jat++)
	{
		const GauAtomBasis& atbas_frag = basis_frag->at_bas_vec[jat];
		const HaAtom* host_atom_frag = atbas_frag.GetAtHost();

		at_bas_frag_pert[jat] = 1.0;		

		int ifst = 0;
		for(iat = 0; iat < natb; iat++)
		{
			GauAtomBasis& atbas = this->at_bas_vec[iat];
			int nb_at   = atbas.GetNBfunc();

            const HaAtom* host_atom = atbas.GetAtHost();
			if( host_atom_frag->IsClose(*host_atom))
			{
				atoms_mol_non_matched.erase(host_atom);
				at_bas_frag_pert[jat] = 0.0;

				if( atbas_frag.GetNBfunc() == nb_at )
				{
					PrintLog("Map fragment atom basis %d to atom basis %d of the host molecule \n", jat, iat); 
										
					int k;
					for( k = 0; k < nb_at ; k++) 
					{
						frag_bas_fun_map[jfst + k] = ifst + k;
					}

				}
				else
				{
					PrintLog("Error in GauBasisSet::MatchBasisSet() \n");
					PrintLog("Found match for the fragment atom but atom bases are different \n");
				}
			}
			ifst += nb_at;
		}
		jfst += atbas_frag.GetNBfunc();
	}
	
	int bf_match = frag_bas_fun_map.size();

	if(bf_match == 0)
	{
		PrintLog("\nWarning: GauBasisSet::MatchBasisSet() \n");
		PrintLog("No match found between atomic bases of the fragment and the host  \n");
		return FALSE;
	}

	if(bas_pert_vec != NULL)
	{
		vector< HaVec_int > at_bas_frag_close;
		at_bas_frag_close.resize(natb_frag);
		
		for( jat = 0; jat < natb_frag; jat++)
		{
			const HaAtom* host_atom = (basis_frag->at_bas_vec[jat]).GetAtHost();
			int jat2;
			for(jat2 = jat + 1; jat2 < natb_frag; jat2++)
			{
				const HaAtom* host_atom_2 = (basis_frag->at_bas_vec[jat2]).GetAtHost();
				if( Vec3D::CalcDistance(host_atom, host_atom_2) < 10.0) 
				{
					at_bas_frag_close[jat].append(jat2);
					at_bas_frag_close[jat2].append(jat);
				}
			}
			
			if( at_bas_frag_pert[jat] < 0.999)
			{
				set<const HaAtom*, less<const HaAtom*> >::iterator itr;
				
				double dist_min = 10.0;
				for( itr = atoms_mol_non_matched.begin(); itr != atoms_mol_non_matched.end(); itr++)
				{
					double dist =  Vec3D::CalcDistance(host_atom, (*itr));
					if( dist < dist_min)  dist_min = dist;
				}
				at_bas_frag_pert[jat] = exp( -1.0*dist_min); 
			}
		}
		
		int i;
		for(i=0; i < 20; i++)
		{
			int ipert_changed = 0;
			for( jat = 0; jat < natb_frag; jat++)
			{
				const HaAtom* host_atom = (basis_frag->at_bas_vec[jat]).GetAtHost();
				int nn = at_bas_frag_close[jat].size();
				int jn;
				
				for(jn = 0; jn < nn; jn++)
				{
					int jnb = at_bas_frag_close[jat][jn];
					const HaAtom* nb_atom = (basis_frag->at_bas_vec[jnb]).GetAtHost();
					double dist =  Vec3D::CalcDistance(host_atom, nb_atom);
					double pert = at_bas_frag_pert[jnb] * exp( -1.0*dist);
					if( pert >  at_bas_frag_pert[jat] )
					{
						at_bas_frag_pert[jat] = pert;
						ipert_changed  = 1;
					}
				}	
			}
			if( !ipert_changed) break;
		}
		
		bas_pert_vec->resize(nb_frag);
		int natb_frag    = basis_frag->at_bas_vec.size();
		int idx_fst = 0;

		for(jat= 0; jat < natb_frag; jat++)
		{
			int nf = basis_frag->at_bas_vec[jat].GetNBfunc();
			int jf;
			for( jf = 0; jf < nf; jf++)
			{
				(*bas_pert_vec)[idx_fst + jf] = at_bas_frag_pert[jat];
			}
			idx_fst += nf;
		}
	}

    return TRUE;
}

int GauBasisSet::GetNumAtBas() const
{
	return at_bas_vec.size();
}

Vec3D* GauBasisSet::GetHostPt(int i)
{
	int idx = GetAtBasIdxForOrb(i);
	GauAtomBasis& atb = GetAtBasByIdx(idx);
	return atb.GetHostPt(0);
}

const Vec3D* GauBasisSet::GetHostPt(int i) const
{
	if( at_bas_vec.size() == 0) return NULL;
	int idx = 0; 
	if(atom_bas_idx.size() > 0) 
		idx = atom_bas_idx[i];
	else
	{
		int naf = GetNumAtBas();
		int k;
		int nl =0;
		for( k = 0; k < naf; k++)
		{
			nl = nl + at_bas_vec[k].GetNBfunc();
			if(i < nl)
			{
				idx = k;
				break;
			}
		}
	}
	return at_bas_vec[idx].GetHostPt(0);
}

int 
GauBasisSet::GetTransfMat(HaMat_double& transf_mat, const HaMat_double& rot_mat)
{
	int nb = GetNBfunc();
    transf_mat.newsize(nb,nb);
	transf_mat = 0.0;
    int natb = GetNumAtBas();
	int i;
	int idx_cur = 0;
	for(i=0; i < natb; i++)
	{
        GauAtomBasis& atb = GetAtBasByIdx(i);
		int nsh = atb.Shells.size();
		int ish;
		for(ish = 0; ish < nsh; ish++)
		{
           GauShell& shl = atb.Shells[ish];
		   HaMat_double tmat;
		   int nf = shl.GetNBfunc();
		   shl.GetTransfMat(tmat,rot_mat);
		   int k,l;
		   for(k =0; k < nf; k++)
		   {
			  for(l = 0; l < nf; l++)
			  {
                  transf_mat.r0(idx_cur + k, idx_cur + l) = tmat.r0(k,l);
			  }
		   }
           idx_cur += nf;
		}
	}
	return TRUE;
}

int 
GauBasisSet::TransferBetweenAtoms(PtrPtrMap& pt_corr_map)
{
	int unknown_atom = FALSE;

	int natb = at_bas_vec.size();
	int i;

	for( i=0; i < natb; i++)
	{
	   HaAtom* host_old = at_bas_vec[i].GetAtHost();
	   HaAtom* host_new = (HaAtom*) pt_corr_map.GetVal( host_old );	
	   if( host_new == NULL )
	   {
		   unknown_atom = TRUE;
	   }
	   else
	   {
		   at_bas_vec[i].SetAtHost( host_new );
	   }
	}
	if( unknown_atom )
	{
		PrintLog("Error in GauBasisSet::TransferBetweenAtoms() \n");
		PrintLog("One Basis Function is centered on the atom unknown to the provide atom-atom correspondence map \n");
		return FALSE;
	}

	return TRUE; 
}

TiXmlElement*
GauBasisSet::AddXml(TiXmlElement* parent_element, const char* name, int option) const
{
	if( parent_element == NULL) return NULL;

	TiXmlElement* basis_set_element;
	
	if( strlen(name) == 0 )
	{
		basis_set_element = new TiXmlElement("GauBasisSet");
	}
	else
	{
		basis_set_element = new TiXmlElement(name);
		basis_set_element->SetAttribute("TYPE", "GauBasisSet");
	}
		
	parent_element->LinkEndChild(basis_set_element);
	
	int natb = at_bas_vec.size();
	int i;

	if( natb > 0)
	{
	   TiXmlElement* atb_vec_element = new TiXmlElement("at_bas_vec");
	   basis_set_element->LinkEndChild(atb_vec_element);
	   atb_vec_element->SetAttribute("size",natb);
	   for(i = 0; i < natb; i++)
	   {
		  at_bas_vec[i].AddXml(atb_vec_element);
	   }
	}

	return basis_set_element;
}

int GauBasisSet::LoadXml(const TiXmlElement* xml_element, int option )
{
	if( xml_element == NULL) return FALSE;

	const TiXmlElement* atb_vec_element = xml_element->FirstChildElement("at_bas_vec");

	int itmp, ires;

	if( atb_vec_element != NULL )
	{
		ires = atb_vec_element->QueryIntAttribute("size", &itmp);
		int nsize = itmp;
		if( ires == TIXML_SUCCESS )
		{
			at_bas_vec.resize(nsize);
		}
		else
		{
			nsize = 0;
			at_bas_vec.resize(nsize);
		}

		const TiXmlElement* child_elem;

		int i = 0;
		for(child_elem =  atb_vec_element->FirstChildElement(); child_elem != NULL;  child_elem = child_elem->NextSiblingElement())
		{
			if( i >= nsize )
			{
				at_bas_vec.push_back(GauAtomBasis());
				nsize++;
			}
			at_bas_vec[i].LoadXml(child_elem);
			i++;
		}
	}

	return TRUE;
}


int
GauBasisSet::GetAtBasIdxForOrb(int i)
{
	if(atom_bas_idx.size() == 0)
	{
		atom_bas_idx.newsize(GetNBfunc());
		int k;
		int j=0;
		int ml=0;
		
		int natf = GetNumAtBas();
		for(k=0; k < natf; k++)
		{
			GauAtomBasis& atb = at_bas_vec[k];
			ml += atb.GetNBfunc();
			while( j < ml)
			{
				atom_bas_idx[j] = k;
				j++;
			}
		}
	}
	return atom_bas_idx[i];
}

int GauBasisSet::GetNumElectr() const
{
	int natb = at_bas_vec.size();
	int i;
	int nel = 0;
	for( i=0 ; i < natb; i++)
	{
		nel += at_bas_vec[i].GetNumElectr();
	}
	return nel;
}

int GauBasisSet::GetCntCoord(HaMat_double& coord) const
{
	int nf = at_bas_vec.size();
	coord.newsize(3,nf);
	int i;
	for(i = 0; i < nf; i++)
	{
		const Vec3D* pt = at_bas_vec[i].GetHostPt(i);
		coord.SetVal_idx0(0,i,pt->GetX_Bohr());
        coord.SetVal_idx0(1,i,pt->GetY_Bohr());
		coord.SetVal_idx0(2,i,pt->GetZ_Bohr());
	}
	return TRUE;
}

int 
GauBasisSet::GetCntCoordArr(Vec3DValArray& crd_arr) const
{
	int nf = at_bas_vec.size();
	crd_arr.resize(nf);
	int i;
	for(i = 0; i < nf; i++)
	{
		const Vec3D* pt = at_bas_vec[i].GetHostPt(i);
		crd_arr[i].SetX(pt->GetX());
		crd_arr[i].SetY(pt->GetY());
		crd_arr[i].SetZ(pt->GetZ());
	}
	return TRUE;
}


