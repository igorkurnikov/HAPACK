/*! \file  halinalg.h

  Classes to define Linear Algebra Objects like Vector, Matrix, Sparse matrix etc.

  \author Igor Kurnikov 
  \date 1999-2002

*/
#ifndef HALINALG_H
#define HALINALG_H

#define TNT_BOUNDS_CHECK


#include "haio.h" 
#include "hastring.h"
//#include "const.h"       

#include "vec.h"     //   TNT library
#include "fmat.h"    
#include "f2c.h"


#include "cg_sngl.h" // include files for IML++

typedef NumVector<short>   HaVec_short;   //!< Numerical vector of short
typedef NumVector<void*>   HaVec_ptr;     //!< Numerical vector of generic pointers (void*)

typedef Fortran_matrix<float>   HaMat_float;
typedef Fortran_matrix<int>     HaMat_int;

class Matrix; // Forward declaration for IPACK double-precision matrix
class TiXmlElement;


#if defined(SWIG)
%template(std_vector_double)  vector<double>;
%template(std_vector_int)  vector<int>;

	class HaVec_int_parent
	{
	public:
		HaVec_int();
		void newsize(int size);
		int* begin();
		int* end();
		int size();
		int GetVal(size_t i) const;
		void SetVal(size_t i, const int new_val);
		int GetVal_idx0(size_t i) const;
		void SetVal_idx0(size_t i, const int new_val);
	};

	class HaVec_float_parent
	{
	public:
		HaVec_float_parent();
		void newsize(int size);
		float* begin();
		float* end();
		int size();
		float GetVal(size_t i) const;
		void SetVal(size_t i, const float new_val);
		float GetVal_idx0(size_t i) const;
		void SetVal_idx0(size_t i, const float new_val);
	};	

    class HaVec_double_parent
	{
	public:
		HaVec_double_parent();
		void newsize(int size);
		
		double* begin();
		double* end();

		double& r0(size_t i);
		double& r1(size_t i);
		double r0(size_t i) const;
		double r1(size_t i) const;

		int size();
		double GetVal(size_t i) const;
		void SetVal(size_t i, const double new_val);
		double GetVal_idx0(size_t i) const;
		void SetVal_idx0(size_t i, const double new_val);
	};
#else
	typedef NumVector<double>  HaVec_double_parent;
 	typedef NumVector<float>   HaVec_float_parent;   
    typedef NumVector<int>     HaVec_int_parent;      
#endif


class HaVec_int: public HaVec_int_parent 
//! class for a numerical array of int
{
public:
	HaVec_int() {}
	HaVec_int( const HaVec_int& A): HaVec_int_parent( (HaVec_int_parent) A) {}
	HaVec_int(size_t N, const int& value = 0 ): HaVec_int_parent(N, value) {}
    HaVec_int(size_t N, const int* v): HaVec_int_parent(N, v) {}
	HaVec_int(size_t N, char *s): HaVec_int_parent(N, s) {}

	HaVec_int& operator=(const HaVec_int_parent &B) { HaVec_int_parent* pp = (HaVec_int_parent*)this; pp->operator=(B); return (*this); } 
	HaVec_int& operator=(const int& scalar) { set(scalar); return *this; } 
	HaVec_int& operator=(const std::vector<int> &B) 
	{ 
		resize(B.size()); 
		int i; 
		for(i=0; i < B.size(); i++) 
			(*this)[i] = B[i]; 
		return *this; 
	}
};

class HaVec_double: public HaVec_double_parent 
//! class for a numerical array of doubles
{
public:
	HaVec_double() {}
	HaVec_double( const HaVec_double& A): HaVec_double_parent( (HaVec_double_parent) A) {}
	HaVec_double(size_t N, const double& value = 0.0 ): HaVec_double_parent(N, value) {}
    HaVec_double(size_t N, const double* v): HaVec_double_parent(N, v) {}
	HaVec_double(size_t N, char *s): HaVec_double_parent(N, s) {}

	HaVec_double& operator=(const HaVec_double_parent &B) { HaVec_double_parent* pp = (HaVec_double_parent*)this; pp->operator=(B); return (*this); } 
	HaVec_double& operator=(const double& scalar) { set(scalar); return *this; } 
	HaVec_double& operator=(const std::vector<double> &B) 
	{ 
		resize(B.size()); 
		int i; 
		for(i=0; i < B.size(); i++) 
			(*this)[i] = B[i]; 
		return (*this); 
	}
//    HaVec_double& operator/(const HaVec_double &B) { HaVec_double_parent* pp = (HaVec_double_parent*)this; pp->operator/(B); return (*this); } 
};

class HaVec_float: public HaVec_float_parent 
//! class for a numerical array of float numbers 
{
public:
	HaVec_float() {}
	HaVec_float( const HaVec_float& A): HaVec_float_parent( (HaVec_float_parent) A) {}
	HaVec_float(size_t N, const float& value = 0.0 ): HaVec_float_parent(N, value) {}
    HaVec_float(size_t N, const float* v): HaVec_float_parent(N, v) {}
	HaVec_float(size_t N, char *s): HaVec_float_parent(N, s) {}

	HaVec_float& operator=(const HaVec_float_parent &B) { HaVec_float_parent* pp = (HaVec_float_parent*)this; pp->operator=(B); return (*this); } 
	HaVec_float& operator=(const float& scalar) { set(scalar); return *this; } 
	HaVec_float& operator=(const std::vector<float> &B) 
	{ 
		resize(B.size()); 
		int i; 
		for(i=0; i < B.size(); i++) 
			(*this)[i] = B[i]; 
		return (*this); 
	}
};


class HaHist 
//! class for histogramms not working yet
{
	public:
		HaHist(double min,double max,double binsize);
		~HaHist();
		
		HaVec_double rleft;
		HaVec_double rcenter;
		HaVec_int ibins;
		HaVec_double dbins;
		
		int PutToIBins(double val);
};
#if defined(SWIG)
	class HaMat_double_parent
	{
	public:
		void reserve(size_t M, size_t N);
		void newsize(size_t M, size_t N);
		void resize(size_t M, size_t N);

		double GetVal(size_t i, size_t j);
		void SetVal(size_t i, size_t j, const double new_val);
		double GetVal_idx0(size_t i, size_t j);
		void SetVal_idx0(size_t i, size_t j, const double new_val);
		double& r0(size_t i, size_t j);
		double& r1(size_t i, size_t j);
		double r0(size_t i, size_t j) const;
		double r1(size_t i, size_t j) const;

		int dim(int d) const;
		int num_rows() const;
		int num_cols() const;

		double* begin();
		double* end();
	};
#else
	typedef Fortran_matrix<double>   HaMat_double_parent ;	
#endif


class HaMat_double : public HaMat_double_parent
//! class for two-dimensional matrix of doubles
{
public:
	HaMat_double() {}
	HaMat_double(size_t M, size_t N, const double& value = 0.0): 
	                                       HaMat_double_parent(M,N,value) {} 

	HaMat_double(size_t M, size_t N, double* v): 
	                                       HaMat_double_parent(M,N,v) {}
										   
	HaMat_double(size_t M, size_t N, char* s): 
	                                       HaMat_double_parent(M,N,s) {}

	HaMat_double(const HaMat_double& A): 
										   HaMat_double_parent( (HaMat_double_parent) A) {}

	HaMat_double(const HaMat_double_parent& A): 
										   HaMat_double_parent(A) {}

	HaMat_double& operator=(double scalar) { set(scalar); return *this; } 

	int set_from(Matrix& ipack_mat); //!< Copy from IPACK 2D Matrix

	int SetFromStr(const char* str); //!< Set from string (expect list of M X N numberes

    static int mat_inverse(HaMat_double& aa);   //!< Inverse matrix
    static int mat_transpose(HaMat_double& aa); //!< Transpose matrix 
    static int mat_sdiag(const HaMat_double& aa, HaMat_double& cc, 
			             HaVec_double& eig);  //!< Diagonalize symmetric matrix

	int SqRoot(const int isign=1); //!< Square root and Inverse square of the matrix

	static int solv_lin_syst_1(HaMat_double& a, HaMat_double& b); //!< Solve system of linear equations 
	static int DiagMat(HaMat_double& hmat, HaMat_double& ss, HaMat_double& eigv, HaVec_double& eig_ene); //!< diagonalize matrix hmat (use overlap matrix ss if it is not empty)

	static int PrintSymmMat(ostream& sout, double* fmat,int nn, const char* frm_str = ""); //!< Print Symmetric matrix

    TiXmlElement* AddXml( TiXmlElement* parent_element,const char* name = "", int option = 0) const;
	virtual int LoadXml( const TiXmlElement* mat_element, int option = 0 ); 
};

#if defined(SWIG)
	int matmult(HaMat_double& C, const HaMat_double& A, const HaMat_double& B);
	int matmult_T1(HaMat_double& C, const HaMat_double& A, const HaMat_double& B);
	int matmult_T2(HaMat_double& C, const HaMat_double& A, const HaMat_double& B);
#endif


class HaSMat_double
//! Symmetric double precision real matrix
{
public:
	HaSMat_double( AllocMode amode_new = INTERNAL_ALLOC, double* ext_ptr = NULL );
	HaSMat_double( int n, AllocMode amode_new = INTERNAL_ALLOC, double* ext_ptr = NULL );
    HaSMat_double( const HaSMat_double& dsmat);
    HaSMat_double( const HaMat_double& dmat);
	virtual ~HaSMat_double();

	int newsize(int n ); 
	
	int CopyFrom(const HaSMat_double& dmat);
	int CopyFrom(const HaMat_double& dmat);
	int CopyTo(HaMat_double& dmat) const;
	
	int num_rows() const { return ndim; } 
	int num_cols() const { return ndim; } 
    
	double& r0(size_t i, size_t j) { return (i <= j) ? v_[j*(j+1)/2 + i] : v_[i*(i+1)/2 + j]; }
	const double& r0(size_t i, size_t j) const { return (i <= j) ? v_[j*(j+1)/2 + i] : v_[i*(i+1)/2 + j]; } 
	
	double* begin(){ return v_;}
	const double* begin() const { return v_;}

	int ndim;

	AllocMode amode;
	

private:
	double* v_;
	NumVector<double> vec_;
};

typedef map<std::string,HaVec_double,less<std::string> >  StrVecMap;
typedef vector<HaMat_double> HaMat_doubleArr;


class MultiVarFunctor
//! Abstract class for function of multiple variables - use as parent class to use with minimizers
{
  public:
      MultiVarFunctor(); 
      virtual ~MultiVarFunctor();
      virtual double operator() (const HaVec_double& xv) = 0;
	  virtual int CalcGrad(HaVec_double& grad,const HaVec_double& xv) = 0;
};

const int BFGS_MIN_METH = 0;

typedef void (*ptrMinFunc1)(int* iptr,int* n, double*x, double* f, double* g);

extern "C" { 
extern void va13ad_ ( int* iptr, ptrMinFunc1 pfunc, 
					  int* n, double* x, double* f ,double* g,double* scale,double* acc,double* w);
}

class HaMinimizer
//! Class to minimize a function of multiple variables
{
public:
	HaMinimizer();
	~HaMinimizer();

	int min_method;

	int GetNVar();              //!< Get the number of independent variables
	void SetNVar( int nvar);    //!< Set the number of independent variables

	int SetInitPoint(HaVec_double& init_pt); //!< Set Initial Point
	virtual int CalcValGrad(HaVec_double& x, double& val, HaVec_double& grad); //!< calc minimized fun value and gradient 
	                                                                               
	int Minimize(); //!< minimize function 

	HaVec_double vvar;  //!< value of indep variable vector
	HaVec_double scale; //!< vector values set to about 10% of expected changes of var in minimization
	double fun_val;     //!< minimized value of the function

	int nitr;
	double flast;
	HaVec_double glast;


};


// Matrix operations as global functions:

#ifdef HALINALG_CPP

size_t ij_indx0(size_t i, size_t j);
size_t ij_indx1(size_t i, size_t j);

extern "C" 
{
	int write_double_array_chuncks( std::ostream& os, HaVec_double& dvec, int chunck_size, const std::string& form_str);
	int write_float_array_chuncks( std::ostream& os, HaVec_float& fvec, int chunck_size,  const std::string& form_str);
	int write_int_array_chuncks( std::ostream& os, HaVec_int& ivec, int chunck_size,  const std::string& form_str);
}

#else

extern size_t ij_indx0(size_t i, size_t j);
extern size_t ij_indx1(size_t i, size_t j);

extern "C" 
{
	int write_double_array_chuncks( std::ostream& os, HaVec_double& dvec, int chunck_size, const std::string& form_str);
	int write_float_array_chuncks( std::ostream& os, HaVec_float& fvec, int chunck_size,  const std::string& form_str);
	int write_int_array_chuncks( std::ostream& os, HaVec_int& ivec, int chunck_size,  const std::string& form_str);
	
//  int rot_vec(double a[], double n[], double cosa, double sina); //!< rotate vector a around unit vector n by angle defined cosa and sina 
  double dot_double(double* a,double* b, int n); //!< dot product of two real vectors
}

#endif

class LanzPars
{
public:
// INPARI array:

   LanzPars();
   ~LanzPars();

   int mat_order;       //!< THE ORDER OF THE MATRIX (INPARI(1)) (OK)
   int max_eigv_store;  //!< THE MAXIMUM NUMBER OF EIGENVALUES THAT CAN BE STORED (OK)
   int num_eigv_search; //!< THE NUMBER OF EIGENVALUES BEING SOUGHT (OK)
   int num_steps;       //!< (ON INPUT) THE MAXIMUM NUMBER OF LANCZOS STEPS TO TAKE (OK)
                            //!< (ON OUTPUT) THE NUMBER OF LANCZOS STEPS TAKEN
   int iret;            //!< THE RETURN VALUE OF LANCZOS (0 MEANS NO ERROR) (INPARI(5))
   int nfound;          //!< (ON INPUT) THE NUMBER OF EIGENVALUES ALREADY FOUND (OK)
                            //!< (ON OUTPUT) THE NUMBER OF EIGENVALUES CORRECTLY FOUND
   int debug_lvl;       //!< THE LEVEL OF DEBUGGING OUTPUT (OK) (INPARI(7))

   int problem_type;    //!< TYPE OF EIGENPROBLEM (0=CLOSEST TO SIGMA, 1=SEARCH IN RANGE) (INPARI(8))

   int inertia_check;   //!< INERTIA CHECKING (0=NONE, 1=YES) (OK)
   int output_amount;   //!< AMOUNT OF OUTPUT (OK) 
                            //  PRINT EIGENVALUES AND INERTIAS =  64 + 4
                            //  PRINT EIGENVALUES, ESTIMATED ERRORS AND INERTIAS = 32 + 4
                            //  PRINT EIGENVALUES, ESTIMATED ERRORS, CALCULATED ERROR AND INERTIAS = 16 + 4
                            //  PRINT EIGENVALUES, ESTIMATED ERRORS,
                            //  CALCULATED ERROR, Y-ORTHOGONALITY AND INERTIAS = 16 + 4 + 2
                            //  PRINT EIGENVALUES, ESTIMATED ERRORS,
                            //  CALCULATED ERROR, Y-ORTHOGONALITY AND INERTIAS
                            //  AND THROW IN SOME FREQUENCIES 128 + 16 + 4 + 2
   int max_steps_shift; //!< THE MAX. NUMBER OF STEPS ON ONE SHIFT(LANCOS RUN) (INPARI(11))

//   int buckling_flag;     //!< VIBRATION (0) OR BUCKLING (1 W/ SHIFT) (2 W/O SHIFT)
   int search_in_boundary;  //!< same as problem type??

   int fc_mat_fmt;      //!< STORAGE FORMAT FOR K (0=SPARSE,1=FULL) (DON'T USE FULL) (INPARI(13)) (OK)
   int m_mat_fmt;       //!< STORAGE FORMAT FOR M (0=SPARSE,1=FULL) (DON'T USE FULL) (INPARI(14)) (OK)
   int loop_unroll_lvl; //!< LEVEL OF LOOP UNROLLING (1,4, OR 6, APPLIES ON TO POSITIVE (INPARI(15)) (OK)
                            //!< DEFINITE SOLVER, THE LEVEL OF
                            //!< UNROLLING FOR THE INDEFINITE
                            //!< SOLVER IS SET IN THE MAKEFILE
   int factor_type;     //!< Factorization type:                INPARI(16) - OK
                            //!< 0 = BUNCH-KAUFMAN IF INDEFINITE
                            //!< LDL IF P.D.
                            //!< 1 = ALWAYS USE LDL
                            //!< 2 = SAME AS 0 BUT SPARSE
                            //!< 3 = SAME AS 1 BUT SPARSE
   int dyn_shift_off;   //!< DYNAMIC SHIFTING TURNED OFF (0=ON, 1=OFF) INPARI(17) - OK
   int init_guess;      //!< INITIAL GUESS (1), NO GUESS (0)           INPARI(18) - OK
   int lead_y_idx;      //!< THE LEADING INDEX OF THE ARRAY Y
   int num_delay_piv;   //!< THE NUMBER OF DELAYED PIVOTS TO ALLOW SPACE FOR

// INPARR array:
   double init_val;     //!< value to search around (shift? )
   double accuracy;     //!< THE REQUIRED RELATIVE ACCURACY OF THE EIGENVALUES
   double left_bound;   //!< LEFT SIDE OF BOUNDARY TO SEARCH (IF ANY)
   double right_bound;  //!< RIGHT SIDE OF BOUNDARY TO SEARCH (IF ANY)
   double store_factor;  //!< THE STORAGE FACTOR FOR BUNCH-KAUFMAN

};

// LANZ - soubroutine to solve eigen vaalue problem with Lancos Algorithm
#if defined LANZ_SOLVER
extern "C" {
void lanz_(int* ldi,  //!< THE LEADING INDEX OF Y
		   int* inpari, double* inparr, //!< parameters of the run see above
		   double* otheta, //!< OLDTHETA(*) IS AN ARRAY OF CONVERGED EIGENVALUES
		   double* oldbj,  //!< OLDBJ(*) IS AN ARRAY OF ERROR BOUNDS ON CONVERGED EIGENVALUES
		   double* y,      //!< Y(LDI,*) IS THE ARRAY OF EIGENVECTORS USED IN A LANCZOS RUN
		   int* onumt, //!< ONUMT(*) IS A WORK VECTOR FOR MSHIFT
		   double* kaux,   //!< KAUX(*) IS THE STIFFNESS MATRIX
		   int* krp,   //!< KRP(*) IS AN AUXILARY INTEGER VECTOR USED WHEN USING KAUX
           int* kcp,   //!< KCP(*) IS AN AUXILARY INTEGER VECTOR USED WHEN USING KAUX
		   double* maux,   //!< MAUX IS THE MASS MATRIX
		   int* mrp,   //!< MRP IS AN AUXILARY INTEGER VECTOR USED WHEN USING MAUX
		   int* mcp,   //!< MCP IS AN AUXILARY INTEGER VECTOR USED WHEN USING MAUX
		   double*  r1j    //!< INITIAL GUESS VECTOR
		   );
}
#endif

// LAPACK subroutines prototypes
extern "C" {

void dgesv_( int* n,  int* nrhs , double* a, 
			 int* lda, int* ipiv , double* b, 
			 int* ldb, int* info );

void dsyev_( char* jobz,  char* uplo , int* n, double* a, 
			 int* lda, double* w, double* work, int* lwork, int* info);

void dsygv_( int* itype, char* jobz,  char* uplo , int* n, 
			 double* a, int* lda, double* b, int* ldb, 
			 double* w, double* work, int* lwork, int* info);

}




#endif /* !HALINALG_H */
