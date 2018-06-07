/*! \file halinalg.cpp 

   Classes and functions for linear algebra operations in HARLEM

   \author Igor Kurnikov

   \date 1998-2003 

*/
#define HALINALG_CPP

#include "haconst.h"
#include "g94_globals.h"
#include "g94_protos.h"
#include "halinalg.h"
#include "math.h"
#include "vtype.h"

#include "tinyxml.h"

size_t ij_indx0(size_t i, size_t j)
// Get position of the (i,j) matrix element 
// in the array presenting packed square matrix zero-based (C) 
{
		return( (i >= j) ? ( i*(i-1)/2 + j-1) 
	                	 : ( j*(j-1)/2 + i-1)  );
	
}

size_t ij_indx1(size_t i, size_t j)
// Get position of the (i,j) matrix element 
// in the array presenting packed square matrix 1-based (Fortran)
{
		return( (i >= j) ? ( i*(i-1)/2 + j) 
	                	 : ( j*(j-1)/2 + i)  );
	
}

extern "C" {
    // LU decomoposition of a general matrix
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

    // generate inverse of a matrix given its LU decomposition
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
}

int HaMat_double::mat_inverse(HaMat_double& aa)
{
//	PrintLog(" HaMat_double::mat_inverse() pt 1 \n");

	int n=aa.num_rows();

	if(n < 1)
	{
		PrintLog("Error in HaMat_double::mat_inverse() \n");
		PrintLog(" Matrix Dimensions < 1 \n");
		return FALSE;
	}
	if( aa.num_rows() != aa.num_cols() )
	{
		PrintLog("Error in HaMat_double::mat_inverse() \n");
		PrintLog("aa matrix is not square \n");
		return FALSE;
	}

    int LWORK = n*n;
  
	HaVec_int ipriv(n+1);
	HaVec_double work(LWORK);
    int INFO;

    dgetrf_(&n,&n,aa.v(),&n,ipriv.v(),&INFO);
    dgetri_(&n,aa.v(),&n,ipriv.v(),work.v(),&LWORK,&INFO);
//	PrintLog(" dgetri_()  info = %d \n", INFO );


// IGOR TMP	HaMat_double bb(n,n,0.0);
// IGOR TMP	int i;
// IGOR TMP	for(i = 0; i < n; i++)
// IGOR TMP	bb.r0(i,i) = 1.0;
// IGOR TMP	result = solv_lin_syst_1(aa,bb);
// IGOR TMP	aa = bb;


//	HaMat_int is(2,n);
//	HaVec_int iad1(n);
//	HaVec_int iad2(n);
//	HaVec_double d(n);  
//	int nn=n;
//	double det=0.0;
//	result= inv1_(aa.begin(),&nn,is.begin(),iad1.begin(),iad2.begin(),
//     		         d.begin(),&nn, &det);

	if(INFO == 0 ) return TRUE;
	
	ErrorInMod("HaMat_double::mat_inverse()",
		       " The matrix is singular ");

	return FALSE;
}

int 
HaMat_double::mat_transpose(HaMat_double& aa)
{
	size_t nn=aa.num_cols();
	if(aa.num_rows()== nn)
	{
		double xx;
        for(size_t i = 2; i <= nn; i++)
		{
			for(size_t j = 1; j < i; j++)
			{
				xx = aa(i,j);
				aa(i,j) = aa(j,i);
				aa(j,i) = xx;
			}
		}
	}
	else
	{	 
		HaMat_double bb=transpose(aa);
		aa=bb;
	}
	return 1;
}

int HaMat_double::mat_sdiag(const HaMat_double& aa,  HaMat_double& cc, HaVec_double& eig)
// diagonalize symmetrical matrix
{
	if(aa.num_cols() != aa.num_rows())
	{
		PrintLog(" Error in HaMat_double::mat_sdiag() \n");
		PrintLog(" aa matrix is not square \n");
	    return FALSE;
	}
	int nn=aa.num_cols();
	
	cc.newsize(nn,nn);
	eig.newsize(nn);

	cc = aa;
	//PrintLog("aa %2.3f %2.3f %2.3f\n",aa(1,1), aa(1,2), aa(1,3));
 //  	PrintLog("aa %2.3f %2.3f %2.3f\n",aa(2,1), aa(2,2), aa(2,3));
 //  	PrintLog("aa %2.3f %2.3f %2.3f\n",aa(3,1), aa(3,2), aa(3,3));

	char jobz[1],uplo[1];

	jobz[0] = 'V';
	uplo[0] = 'U';
    
	int info;
	double wk_opt;
	int lwork = -1;
// Call to determine an optimal size of work array:
	dsyev_(&jobz[0],&uplo[0],&nn, cc.begin(),&nn,eig.begin(),&wk_opt, &lwork, &info);
	lwork = wk_opt;
	HaVec_double work(lwork+1);
	dsyev_(&jobz[0],&uplo[0],&nn, cc.begin(),&nn,eig.begin(),work.begin(), &lwork, &info); 
    if (info < 0)
		PrintLog("INFO = %d, the i-th argument had an illegal value \n", info);
	if (info >0)
		PrintLog("INFO = %d, the algorithm failed to converge \n", info);

	//PrintLog("cc %2.3f %2.3f %2.3f\n",cc(1,1), cc(1,2), cc(1,3));
 //  	PrintLog("cc %2.3f %2.3f %2.3f\n",cc(2,1), cc(2,2), cc(2,3));
 //  	PrintLog("cc %2.3f %2.3f %2.3f\n",cc(3,1), cc(3,2), cc(3,3));

	return TRUE;

}

int HaMat_double::SqRoot(const int isign)
{
	int inv;
	int nn= num_cols();
	assert(nn > 0);
	assert(nn == num_rows());
	if(isign > 0) 
	{
		inv = 0;
	}
	else 
	{
		inv = 1;
	}

	HaMat_double b(nn,nn);
	HaMat_double c(nn,nn);

	HaVec_double aa(nn);
	HaVec_double bb(nn);

#if defined(GAUSSVER)

    rootmt_(this->begin(),b.begin(),aa.begin(),bb.begin(),&nn,&nn,&inv);

#else

	HaMat_double::mat_sdiag(*this, c, aa);

	int i,j;
	for(i = 0; i < nn; i++)
	{
		if(aa[i] <= 0 )
		{
			PrintLog(" Error in HaMat_double::SqRoot() \n");
			PrintLog(" Matrix has non-positive eigen values \n");
			return FALSE;
		}
		aa[i] = sqrt(aa[i]);
		if(isign < 0 ) aa[i] = 1.0/aa[i];
	}

	for(j=0; j < nn; j++)
	{
		for(i = 0; i < nn; i++)
		{
			b.r0(i,j) = aa[j]*c.r0(i,j);
		}
	}
    matmult_T2(*this,b,c);

#endif
	
	return True;
}

extern "C" DllExport
int write_double_array_chuncks( std::ostream& os, HaVec_double& dvec, int chunck_size, const std::string& form_str)
{
	char buf[256];
	int exp_format = 0;
	int i;
//	for(i= 0; i < strlen(form_str); i++)
//	{
//		if( form_str[i] == 'E')
//		{
//			exp_format = 1;
//			break;
//		}
//	}

	int nsize = dvec.size();
		if(nsize == 0)
		{
			os << endl;
			return TRUE;
		}

	int loc_idx = 0;

	for(i = 0; i <= nsize ; i++)
	{
		if( i >= nsize) 
		{
			if( loc_idx != 0) 
			{
				os << endl;
			}
			break;
		}
//#if defined(_MSC_VER) 
//		if( exp_format)  // microsoft printf(%xx.xE) format give three digit exponents   
//		{                // reduce exponents to two digits            
//			int n = sprintf(buf,form_str.c_str(),dvec[i]);
//			for(int j= n-3; j > 0; j--)
//			{
//				buf[j] = buf[j-1];
//			}
//			buf[0]= ' ';
//			os << buf;
//		}
//		else
//#endif
		{
			sprintf(buf,form_str.c_str(), dvec[i] );
			os << buf;
		}
		loc_idx++;
		if( loc_idx >=  chunck_size)
		{
			os << endl;
			loc_idx = 0;
		}
	}
	return TRUE;
}

extern "C" DllExport
int write_float_array_chuncks( FILE* fp, HaVec_double& fvec, int chunck_size, const std::string& form_str)
{
	if( fp == NULL)
		return FALSE;

	int exp_format = 0;
	int i;
	for(i= 0; i < form_str.size(); i++)
	{
		if( form_str[i] == 'E')
		{
			exp_format = 1;
			break;
		}
	}

	int nsize = fvec.size();
		if(nsize == 0)
		{
			fprintf(fp,"\n");
			return TRUE;
		}

	int loc_idx = 0;

	for(i = 0; i <= nsize ; i++)
	{
		if( i >= nsize) 
		{
			if( loc_idx != 0) 
			{
				fprintf(fp, "\n");
			}
			break;
		}
#if defined(_MSC_VER) 
		if( exp_format)  // microsoft printf(%xx.xE) format give three digit exponents   
		{                // reduce exponents to two digits
			char buf[128];            
			int n = sprintf(buf,form_str.c_str(),fvec[i]);
			for(int j= n-3; j > 0; j--)
			{
				buf[j] = buf[j-1];
			}
			buf[0]= ' ';
			fprintf(fp,"%s",buf);
		}
		else
#endif
			fprintf(fp, form_str.c_str(), fvec[i] );
		loc_idx++;
		if( loc_idx >=  chunck_size)
		{
			fprintf(fp,"\n");
			loc_idx = 0;
		}
	}
	return TRUE;
}

extern "C" DllExport
int write_int_array_chuncks( std::ostream& os, HaVec_int& ivec, int chunck_size, const std::string& form_str)
{
	char buf[128];

	int nsize = ivec.size();
	if(nsize == 0)
	{
		os << endl;
		return TRUE;
	}

	int loc_idx = 0;

	for(int i = 0; i <= nsize ; i++)
	{
		if( i >= nsize) 
		{
			if( loc_idx != 0) 
			{
				os << endl;
			}
			break;
		}
		sprintf(buf, form_str.c_str(), ivec[i] );
		os << buf;
		loc_idx++;
		if( loc_idx >=  chunck_size)
		{
			os << endl;
			loc_idx = 0;
		}
	}
	return TRUE;
}

int rot_vec(double a[], double n[], double ca, double sa)
{	 
	double atx = (n[0]*n[0] + (1-n[0]*n[0])*ca)*a[0]+
	             (n[0]*n[1]*(1-ca)  -  n[2]*sa)*a[1]+ 
				 (n[0]*n[2]*(1-ca)  +  n[1]*sa)*a[2];

	double aty = (n[0]*n[1]*(1-ca)  +  n[2]*sa)*a[0]+
		         (n[1]*n[1] + (1-n[1]*n[1])*ca)*a[1]+
                 (n[1]*n[2]*(1-ca)  -  n[0]*sa)*a[2];

	double atz = (n[0]*n[2]*(1-ca) - n[1]*sa)*a[0]+
		         (n[1]*n[2]*(1-ca) + n[0]*sa)*a[1]+
                 (n[2]*n[2] + (1-n[2]*n[2])*ca)*a[2];	

	a[0] = atx;
	a[1] = aty;
	a[2] = atz;

	return TRUE;

}

extern "C" double dot_double(double* a,double* b, int n)
{
   double dot = 0.0;
   int i;
   for( i = 0; i < n; i++)
   {
      dot += a[i]*b[i];
   }
   return dot;
}

int HaMat_double::solv_lin_syst_1(HaMat_double& a, HaMat_double& b)
//! a(N,N) - matrix of system coefficients (destroyed)
//! b(N,M) - right hand side matrix - on exit solution matrix
{
	int n = a.num_cols();

	if(n == 0)
	{
		ErrorInMod("HaMat_double::solv_lin_syst_1() ",
			       "Matrix A is empty");
		return FALSE;
	}

	if( a.num_cols() != a.num_rows())
	{
		ErrorInMod("HaMat_double::solv_lin_syst_1() ",
			       "Matrix A is not square");
		return FALSE;
	}

	if( a.num_cols() != b.num_rows())
	{
		ErrorInMod("HaMat_double::solv_lin_syst_1() ",
			       "Number of rows in matrix B not equal to number of cols in matrix A");
		return FALSE;
	}

	int lda = n;
	int ldb = n;
	int nrhs = b.num_cols();
	HaVec_int ipiv(n);
	int info;

    dgesv_( &n,  &nrhs , a.begin(), 
			&lda, ipiv.begin() , b.begin(), 
			&ldb, &info );

	if(info != 0 )
	{		
		ErrorInMod("HaMat_double::solv_lin_syst_1() ",
			       "Error solving system of linear equations");
		return FALSE;
	}	
	return TRUE;

}

int HaMat_double::set_from(Matrix& ipack_mat)
{
	int n1 = ipack_mat.size1();
	int n2 = ipack_mat.size2();

	this->newsize(n1, n2);

    int i,j;
	for(i=0; i < n2; i++)
	{
	   for(j=0; j < n1; j++)
	   {
		   this->SetVal(j+1,i+1,ipack_mat(j,i));
	   }
	}
	return TRUE;
}

int HaMat_double::SetFromStr(const char* str)
{
	istrstream is(str);
	int i;
	int n = num_rows() * num_cols();

	for( i = 0; i < n; i++)
	{
		is >> v_[i];
		if(is.fail())
		{
			ErrorInMod(" HaMat_double::SetFromStr()", 
				       " Failed to read a matrix element ");
			return FALSE;
		}
	}
	return TRUE;
}

int HaMat_double::DiagMat(HaMat_double& hmat, HaMat_double& ss, HaMat_double& eigv, HaVec_double& eig_ene)  
{
	int use_ovlp_mat = TRUE;
    int nb = hmat.num_cols();

	if( ss.num_cols() == 0)
	{
		use_ovlp_mat = FALSE;
	}
	else if(ss.num_rows() != nb ||  ss.num_cols() != nb )
	{
		PrintLog("Error in HaQCMod::DiagMat() \n");
		PrintLog(" Different dimensions of hamiltonian and overlap matricies \n");
		return FALSE;
	}

	HaMat_double ss_m12;

    eigv = hmat;
	eig_ene.newsize(nb);
		
    if(use_ovlp_mat)
	{
		HaMat_double ss_copy = ss;
		char jobz[1],uplo[1];
		
		jobz[0] = 'V';
		uplo[0] = 'U';
		
		int info;
		double wk_opt;
		int lwork = -1;
		int itype = 1; // type of generalized eigenvalue problem A*x=(lambda)*B*x
		
		// Call to determine an optimal size of work array:
		dsygv_(&itype,&jobz[0],&uplo[0],&nb, 
			eigv.begin(),&nb, ss_copy.begin(), &nb, 
			eig_ene.begin(),&wk_opt, &lwork, &info);
		lwork = wk_opt;
		HaVec_double work(lwork+1);
		dsygv_(&itype,&jobz[0],&uplo[0],&nb, 
			eigv.begin(),&nb, ss_copy.begin(), &nb, 
			eig_ene.begin(),work.begin(), &lwork, &info);
	}
	else
	{
		HaMat_double::mat_sdiag(hmat, eigv, eig_ene);
	}
//	if(use_ovlp_mat)
//	{
//		ss_m12 = ss;
//		ss_m12.SqRoot(-1);
//	}

//	HaMat_double scr;
//
//	if( use_ovlp_mat )
//	{
//		matmult(scr,ss_m12,hmat);
//		matmult(hmat,scr,ss_m12); // hmat = S^(-1/2) H S^(-1/2) 
//	}
//
//	HaMat_double::mat_sdiag(hmat, scr, eig_ene);
//	
//	if( use_ovlp_mat)
//	{
//		matmult(eigv,ss_m12,scr); // C = S^(-1/2) C'
//	}
//	else
//	{
//		eigv = scr;
//	}
	
	return true;
}

int HaMat_double::PrintSymmMat(ostream& sout, double* fmat,int nn, const char* frm_str)
{
	char buf[128];
	std::string format = frm_str;
	if( format.empty() ) format = "%12.6f";
	int i,j;
	int ij = 0;
	for( i = 0; i < nn; i++)
	{
		for(j = 0; j <=  i; j++)
		{
			sprintf(buf,format.c_str(),fmat[ij]);
			sout << buf;
			ij++;
		}
		sout << endl;
	}
	return TRUE;
}

TiXmlElement*
HaMat_double::AddXml( TiXmlElement* parent_element, const char* name,  int option) const
{
	if( parent_element == NULL) return NULL;

	TiXmlElement* mat_element;
	if(strlen(name) > 0 )
	{
		mat_element = new TiXmlElement(name);
		mat_element->SetAttribute("TYPE","HaMat_double");		
	}
	else
	{
		mat_element = new TiXmlElement("HaMat_double");
	}

	parent_element->LinkEndChild(mat_element);

	int nc =  this->num_cols();
	int nr =  this->num_rows();

	mat_element->SetAttribute("NCOL", nc);
	mat_element->SetAttribute("NROW", nr);

	int i,j;
		
//	char* buf = (char*) malloc( 20* this->num_cols()* this->num_rows());

	int ipos = 0;

	ostrstream str_buf;
	str_buf.precision(9);
	str_buf.width(16);

//	ipos += sprintf(buf+ipos,"\n");

	str_buf << "\n";

	for(j=0; j < nc; j++)
	{
		for(i = 0; i < nr; i++)
		{
//			ipos += sprintf(buf+ipos,"%16.9e ", this->r0(i,j));
			str_buf << this->r0(i,j) << " ";
		}
//		ipos += sprintf(buf+ipos,"\n");
		str_buf << "\n";
	}

	char* buf =  str_buf.str();
	int nsize = str_buf.pcount(); 
	buf[nsize] = 0;
//
//	ipos += sprintf(buf+ipos,"\n");
	TiXmlText* mat_text = new TiXmlText( buf );
	
	str_buf.freeze(0);

//	free(buf);

	mat_element->LinkEndChild( mat_text );

	return mat_element;
}

int HaMat_double::LoadXml( const TiXmlElement* mat_element, int option )
{
	if( mat_element == NULL ) 
	{
		PrintLog(" Error in HaMat_double::LoadXml() \n");
		PrintLog(" mat_element is NULL \n");
		newsize(0,0);
		return FALSE;

	}

	int nc;
	int nr;
	int ires;

	ires = mat_element->QueryIntAttribute("NCOL",&nc);

	if(ires == TIXML_NO_ATTRIBUTE ) 
	{
		PrintLog(" Error in HaMat_double::LoadXml() \n");
		PrintLog(" NCOL is not set \n");
		return FALSE;
	}

    ires = mat_element->QueryIntAttribute("NROW",&nr);

	if( ires == TIXML_NO_ATTRIBUTE ) 
	{
		PrintLog(" Error in HaMat_double::LoadXml() \n");
		PrintLog(" NROW is not set \n");
		return FALSE;
	}
	
	newsize(nr,nc);

	const char* elem_text = mat_element->GetText();

	istrstream istr(elem_text);

	int i;
	int j; 
	
	for(j=0; j < nc; j++)
	{
		for(i = 0; i < nr; i++)
		{
			double tmp;
			istr >> tmp;
			this->r0(i,j) = tmp;
		}
	}
	return TRUE;
}

HaSMat_double::HaSMat_double(AllocMode amode_new, double* ext_ptr )
{
   amode = amode_new;
   ndim = 0;
   v_ = NULL;
   if(amode == EXTERNAL_ALLOC)
   {
	   v_ = ext_ptr;
   }
}

HaSMat_double::HaSMat_double( int n, AllocMode amode_new, double* ext_ptr )
{
   amode = amode_new;
   ndim = 0;
   v_ = NULL;
   if(amode == EXTERNAL_ALLOC)
   {
	   v_ = ext_ptr;
   }
   if(n > 0)
   {
	  newsize(n);
   }
}


HaSMat_double::HaSMat_double(const HaSMat_double& dsmat)
{
	ndim = 0;
	v_ = NULL;
    CopyFrom(dsmat);
}

HaSMat_double::HaSMat_double(const HaMat_double& dmat)
{
	ndim = 0;
	v_ = NULL;
    CopyFrom(dmat);
}

HaSMat_double::~HaSMat_double()
{
    
}

int
HaSMat_double::newsize(int n)
{
	if( n < 0) return FALSE;
	ndim = n;
	if(amode == EXTERNAL_ALLOC)
	{

	}
	else
	{
	   vec_.newsize(n*(n+1));
	   v_ = vec_.begin();
	}
	return TRUE;
}



int 
HaSMat_double::CopyFrom(const HaSMat_double& dsmat)
{
   int nr = dsmat.num_rows();
   if(nr < 0) return FALSE;
   newsize(nr);
   if( nr == 0) return TRUE;
   memcpy(v_, dsmat.begin() , sizeof(double)*nr*(nr+1)/2);
   return TRUE;
}


int 
HaSMat_double::CopyFrom(const HaMat_double& dmat)
{
   int nr = dmat.num_rows();
   int nc = dmat.num_cols();
   if( nr != nc)
   {
       PrintLog("Error in HaSMat_double::CopyFrom() \n");
	   PrintLog(" nr = %d not equal nc = %d  \n",nr,nc);
	   return FALSE;
   }
   newsize(nr);
   int i,j;
   int ij = 0;
   for(i= 0; i < nr; i++)
   {
	  for(j = 0; j <= i; j++)
	  {
          v_[ij] = dmat.r0(i,j);
		  ij++;
	  }
   }
   return TRUE;
}

int 
HaSMat_double::CopyTo(HaMat_double& dmat) const
{
   if( dmat.num_rows() != ndim || dmat.num_cols() != ndim)
   {
       dmat.newsize(ndim,ndim);
   }
   int i,j;
   int ij = 0;
   for(i= 0; i < ndim; i++)
   {
	  for(j = 0; j <= i; j++)
	  {
          dmat.r0(i,j) = v_[ij];
          dmat.r0(j,i) = v_[ij];
		  ij++;
	  }
   }
   return TRUE;
}


MultiVarFunctor::MultiVarFunctor()
{

}

MultiVarFunctor::~MultiVarFunctor()
{

}

HaMinimizer::HaMinimizer()
{
	min_method = BFGS_MIN_METH;
	nitr = 0; 
	flast = 0.0;
}

HaMinimizer::~HaMinimizer()
{

}

int
HaMinimizer::GetNVar()
{
	return vvar.size();
}

void
HaMinimizer::SetNVar(int nvar)
{
	vvar.newsize(nvar);
}

int 
HaMinimizer::SetInitPoint(HaVec_double& init_pt) //!< Set Initial Point
{
	if(init_pt.size() != vvar.size())
	{
		PrintLog("Error in HaMinimizer::SetInitPoint() \n");
		PrintLog(" init val array size %5d != indep var array size %5d \n",
			     init_pt.size(),vvar.size());
		PrintLog(" Resize indep var array size \n");
	}
	vvar = init_pt;
	return TRUE;
}

extern "C"
void fmin_s(int* iptr,int* pn, double* x, double* f, double* g)
{
	HaMinimizer* pmin = (HaMinimizer*) iptr;
	int n = *pn;
	HaVec_double xv(n),grad(n);
	double val;
	int i;
	for(i = 0; i < n; i++)
	{
		xv[i] = x[i];
	}
	int ires = pmin->CalcValGrad(xv,val,grad);
	if(!ires)
	{
		PrintLog(" Error in fmin_s(): Error in evaluating function and gradients "); 
	}
	(*f) = val;
	for(i = 0; i < n; i++)
	{
		g[i] = grad[i];
	}
}

int 
HaMinimizer::Minimize()
{
	int nv = vvar.size();
	int i;
	if(scale.size() == 0)
	{
		scale.newsize( vvar.size());
	    for(i=0; i < nv; i++)
		{
		    double sh = MaxFun(0.1, fabs(vvar[i])/20.0);
		    scale[i] = sh;
		}
	}

	HaVec_double ww(nv*(nv+13)/2);
	HaVec_double grad(nv);
	double acc = 10E-5;
	
//	va13ad_(obj_ptr,ptf1,&nv,vvar.begin(),&fun_val,grad.begin(),scale.begin(),&acc,ww.begin());
	va13ad_((int*)this,&fmin_s,&nv,vvar.begin(),&fun_val,grad.begin(),scale.begin(),&acc,ww.begin());

	PrintLog("Min val of fun = %12.6e \n",fun_val);
//	PrintLog("Value of indep of var vector \n");
//	for(i=0; i < nv; i++)
//	{
//		PrintLog(" %d %12.6e \n",i+1,vvar[i]);
//	}

	return true;
}


int
HaMinimizer::CalcValGrad(HaVec_double& x, double& val, HaVec_double& grad)
{
	int n = x.size();
	val = 0.0;
	grad.newsize(n);
	grad = 0.0;
	return FALSE;
}

LanzPars::LanzPars()
{
	mat_order = 1;
	max_eigv_store = 1;
	num_eigv_search = 1;
	num_steps = 216;
	iret = 0;
	nfound = 0;
	debug_lvl = 1;
	problem_type = 0;

	inertia_check = 0;
	output_amount = 16 + 4 + 2;
	max_steps_shift = 45; // 45 recomended for range eigenvalues search
	                      // 35 - recomended for 
//	buckling_flag = 0;
	fc_mat_fmt = 0;
	m_mat_fmt = 0;
	loop_unroll_lvl = 6; // in DPREP said always 6 
	factor_type = 0;
	dyn_shift_off = 0;
	init_guess = FALSE;
	lead_y_idx = 1;
	num_delay_piv = 10; // in input files 8,10,0 - don't know what it is

	init_val = 0.0;
	accuracy = 1.0E-9;
	left_bound = -0.5;
	right_bound = 0.0;
	store_factor = 1.0; // in input files saw 1.0,1.1 













	


}

LanzPars::~LanzPars()
{

}
///////////////////////////////////////////////////////////////////////////////
HaHist::HaHist(double min,double max,double binsize)
{
	int nbins=((max-min)/binsize);
	nbins++;
	rleft.resize(nbins);
	ibins.resize(nbins);
	dbins.resize(nbins);
	
}
HaHist::~HaHist()
{
}
int HaHist::PutToIBins(double val)
{
	return EXIT_SUCCESS;
}















