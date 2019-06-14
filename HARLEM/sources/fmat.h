/*! \file fmat.h

   Template Numerical Toolkit (TNT) for Linear Algebra

   BETA VERSION INCOMPLETE AND SUBJECT TO CHANGE
   Please see http://math.nist.gov/tnt for updates
  
   R. Pozo
   Mathematical and Computational Sciences Division
   National Institute of Standards and Technology

   Fortran-compatible matrix: column oriented, 1-based (i,j) indexing
*/

#ifndef FMAT_H
#define FMAT_H

#include "vec.h"
#include <stdlib.h>
#include <assert.h>

#ifdef TNT_USE_REGIONS
#include "region2d.h"
#endif

//<mikola have moved to vec.h
	//enum AllocMode {INTERNAL_ALLOC, EXTERNAL_ALLOC}; 
	// determine if Storage for v_ is allocated when matrix is 
	// created and is destroyed when matrix is destroyed or   
    // v points to externally allocated memory region assumed 
	// to be of the sufficient size (m_*n_)
//>mikola 


//! A simple 1-based, column oriented Matrix class
template <class T>
class Fortran_matrix 
{
  public:

    typedef    size_t   size_type;
    typedef         T   value_type;
    typedef         T   element_type;
    typedef         T*  pointer;
    typedef         T*  iterator;
    typedef         T&  reference;
    typedef const   T*  const_iterator;
    typedef const   T&  const_reference;

    size_t lbound() const { return 1;}
 
	AllocMode amode;

  protected:
    T* v_;               //!< these are adjusted to simulate 1-offset
    size_type m_;
    size_type n_;
    T** col_;           //!< these are adjusted to simulate 1-offset

	size_type nalloc_ ; //!< maximal number of elements in matrix fitting the allocated memory 
		                   
  public:


  public:

    T* v() { return v_; }
	const T* v() const { return v_;}
    T* begin() { return v_; }
    const T* begin() const { return v_;}

    T* end() { return v_ + m_*n_; }
    const T* end() const { return v_ + m_*n_; }


    // constructors

    Fortran_matrix() : v_(0), m_(0), n_(0), col_(0),
		               amode(INTERNAL_ALLOC), nalloc_(0) {}

    Fortran_matrix(const Fortran_matrix<T> & A)
    {
		amode=INTERNAL_ALLOC; 
		nalloc_ = 0;
		m_   = 0;
		n_   = 0;
		v_   = NULL;
		col_ = NULL;
        newsize(A.m_, A.n_);
        copy(A.v_);
    }

    Fortran_matrix(size_type M, size_type N, const T& value = T(0))
    {
		amode=INTERNAL_ALLOC; 
		nalloc_ = 0;
		m_   = 0;
		n_   = 0;
		v_   = NULL;
		col_ = NULL;

        newsize(M,N);
        set(value);
    }

	Fortran_matrix(size_type M, size_type N,  T* v)
    {
		amode=INTERNAL_ALLOC; 
		nalloc_ = 0;
		m_   = 0;
		n_   = 0;
		v_   = NULL;
		col_ = NULL;

        newsize(M,N);
        copy(v);
    }

//    Fortran_matrix(size_type M, size_type N, const T* v, 
//		const int New_amode=0)
//    {
//		Set_amode(New_amode);
//		if(amode == EXTERNAL_ALLOC)
//		{
//			MaxSize=M*N;
//			v_= v;
//		}
//       initialize(M,N);
//        if(amode == INTERNAL_ALLOC) copy(v);
//    }


    Fortran_matrix(size_type M, size_type N, char *s)
    {
		amode=INTERNAL_ALLOC; 
		nalloc_ = 0;
		m_   = 0;
		n_   = 0;
		v_   = NULL;
		col_ = NULL;

        newsize(M,N);
        istrstream ins(s);

        size_type i, j;

        for (i=0; i< M; i++)
            for (j=0; j< N; j++)
                ins >> (*this).r0(i,j);
    }

    // destructor
    ~Fortran_matrix()
    {
        destroy();
    }

	void set_ext_alloc(T* v, size_type M, size_type N, size_type N_alloc = 0)
	{
		if(v == NULL) 
		{
			throw("Error in Fortran_matrix<T>::set_ext_alloc() v == NULL \n");
		}

		destroy();
		amode = EXTERNAL_ALLOC;
		
		v_ = v;
		m_ = 0;
		n_ = 0;
		nalloc_ = N_alloc;
        
		newsize(M,N);
	}


    // assignments
    //
    Fortran_matrix<T>& operator=(const Fortran_matrix<T> &A)
    {
        if (v_ == A.v_)
            return *this;

        if (m_ == A.m_  && n_ == A.n_)      // no need to re-alloc
            copy(A.v_);

        else
        {
            newsize(A.m_, A.n_);
            copy(A.v_);
        }

        return *this;
    }
        
    Fortran_matrix<T>& operator=(const T& scalar)
    { 
        set(scalar); 
        return *this;
    }

	void reserve(size_type M, size_type N)
//!< Reserve space for N columns of size M
	{
		if(amode == EXTERNAL_ALLOC)
		{
			nalloc_ = M*N;
			return;
		}

		if(M*N <= nalloc_) return;
	
		T*  v_new   = (T*)  realloc(v_, M*N*sizeof(T));

		if( v_new == NULL)
		{
			v_ = NULL;
			m_ = 0;
			n_ = 0;
			nalloc_ = 0;

			ostrstream sstr; 
			sstr << endl << "Error in Fortran_matrix::reserve() " << endl;
			sstr << "Failed to allocate memory block of size " << M << " X " << N << "elements " << endl;
			throw sstr.str();
		}

		v_   = v_new; 
		nalloc_ = M*N;
	}

    
	void newsize(size_type M, size_type N)
//! Set new dimensions of the matrix M x N
    {
        if (num_rows() == M && num_cols() == N) return;

		if(nalloc_ < M*N) reserve(M,N);

		if( N > num_cols() )
		{
			T** col_new = (T**) realloc(col_, N*sizeof(T*));

			if( col_new == NULL)
			{
				destroy();
				ostrstream sstr;
				sstr << endl << "Error in Fortran_matrix::newsize() " << endl;
				sstr << "Failed to allocate space for columns poiters of matrix " << M << " X " << N << "elements " << endl;
				throw sstr.str();	
			}
			col_ = col_new;
		}

		T* p = v_;              
		for (size_type i=0; i<N; i++)
		{
            col_[i] = p;
            p += M ;
		}

		m_ = M;
		n_ = N;
    }

	inline void resize(size_type M, size_type N)
	{
		newsize(M,N);
	}	
	
	void destroy()
    {     
        /* do nothing, if no memory has been previously allocated */
        if (v_ == NULL) return ;

		if( col_)
		{
			free (col_);
			col_ = NULL;
		}

		if(amode == INTERNAL_ALLOC)
		{
			free (v_); // do not delete if externally allocated 
		}
		
		v_ = NULL;
		amode = INTERNAL_ALLOC;
		m_ = 0;
		n_ = 0;
		nalloc_ = 0;
    }

	void clear()
	{
		destroy();
	}
   
    void copy(const T*  v)
    {
        size_type N = m_ * n_;
        size_type i;

        for (i=0; i< N; i++)
            v_[i] = v[i];
    }

//! Set all values of the matrix equal 'val'
    void set(const T& val)
    {
        size_type N = m_ * n_;
        size_type i;

        for (i=0; i< N; i++)
            v_[i] = val;
    }

    size_type dim(size_type d) const 
    {
#ifdef TNT_BOUNDS_CHECK
        assert( d >= 1);
        assert( d <= 2);
#endif
        return (d==1) ? m_ : ((d==2) ? n_ : 0); 
    }

    size_type num_rows() const { return m_; }
    size_type num_cols() const { return n_; }

//! 1-based element access
    inline reference operator()(size_type i, size_type j)
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= m_) ;
        assert(1<=j);
        assert(j <= n_);
#endif
        return col_[j-1][i-1]; 
    }

    inline const_reference operator() (size_type i, size_type j) const
    {
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= m_) ;
        assert(1<=j);
        assert(j <= n_);
#endif
        return col_[j-1][i-1]; 
    }

//! 0-based linear element access
    inline reference operator[](size_type i)
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(i <  m_*n_) ;
#endif
        return v_[i]; 
    }

    inline const_reference operator[] (size_type i) const
    {
#ifdef TNT_BOUNDS_CHECK
        assert(i <  m_*n_) ;
#endif
        return v_[i]; 
    }

//! 0-based element access
	inline reference r0(size_type i, size_type j)
	{
#ifdef TNT_BOUNDS_CHECK
        assert(i <  m_);
		assert(j <  n_);
#endif
		return col_[j][i];
	}

	inline value_type r0(size_type i, size_type j) const
	{
#ifdef TNT_BOUNDS_CHECK
        assert(i <  m_);
		assert(j <  n_);
#endif
		return col_[j][i];
	}

//! 1-based element access
    inline reference r1(size_type i, size_type j)
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= m_) ;
        assert(1<=j);
        assert(j <= n_);
#endif
        return col_[j-1][i-1]; 
    }

	inline value_type r1(size_type i, size_type j) const    
    {
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= m_) ;
        assert(1<=j);
        assert(j <= n_);
#endif
        return col_[j-1][i-1]; 
	}

	inline value_type GetVal(size_type i, size_type j) const    
    {
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= m_) ;
        assert(1<=j);
        assert(j <= n_);
#endif
        return col_[j-1][i-1]; 
	}


	inline void SetVal(size_type i, size_type j, const value_type new_val)
	{
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= m_) ;
        assert(1<=j);
        assert(j <= n_);
#endif
        col_[j-1][i-1] = new_val; 
	}

	inline value_type GetVal_idx1(size_type i, size_type j) const    
    {
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= m_) ;
        assert(1<=j);
        assert(j <= n_);
#endif
        return col_[j-1][i-1]; 
	}


	inline void SetVal_idx1(size_type i, size_type j, const value_type new_val)
	{
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= m_) ;
        assert(1<=j);
        assert(j <= n_);
#endif
        col_[j-1][i-1] = new_val; 
	}

	inline value_type GetVal_idx0(size_type i, size_type j) const    
    {
#ifdef TNT_BOUNDS_CHECK
        assert(i < m_) ;
        assert(j < n_);
#endif
        return col_[j][i]; 
	}

	inline void SetVal_idx0(size_type i, size_type j, const value_type new_val)
	{
#ifdef TNT_BOUNDS_CHECK
        assert(i < m_) ;
        assert(j < n_);
#endif
        col_[j][i] = new_val; 
	}

	inline value_type GetVal_lidx(size_type i) const    
    {
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= m_*n_) ;
#endif
        return v_[i]; 
	}

	inline void SetVal_lidx(size_type i, const value_type new_val)
	{
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= m_*n_) ;
#endif
        v_[i] = new_val; 
	}

//friend istream& operator>><>(istream &s, Fortran_matrix<T> &A); 
        
#ifdef TNT_USE_REGIONS

    typedef Region2D<Fortran_matrix<T> > Region;
    typedef const_Region2D<Fortran_matrix<T>,T > const_Region;

    Region operator()(const Index1D &I, const Index1D &J)
    {
        return Region(*this, I,J);
    }

    const_Region operator()(const Index1D &I, const Index1D &J) const
    {
        return const_Region(*this, I,J);
    }

#endif


	int Print_format(ostream &sout, const char* format) const
	{
		
		char buf[120]; // set max 20 characters per number
		
		sout << m_ << " " << n_ << endl;
		
		for (size_type i=1; i<= m_; i++)
		{
			for (size_type j=1; j<= n_; j++)
			{
				sprintf(buf, format, (*this)(i,j) );
				sout << buf;
			}
			sout << endl;
		}
		
		return 0;
	}

    int GetSymmPart(pointer dsymm)
		// get Symmetrical part of the square matrix in low-diagonal form:
	{
		assert(dsymm != NULL);
		assert(m_ == n_);
		
		for (size_type i=1; i<= m_; i++)
		{
			for (size_type j=1; j<= i; j++)
			{
				(*dsymm)= ((*this)(i,j) + (*this)(j,i))/2;
				dsymm++;
			}
		}
		return 1;		
	}

    int GetASymmPart(pointer dasymm)
    // get AntiSymmetrical part of the square matrix in lower triangular form
	// (diagonal inlcuded)
	{
		assert(dasymm != NULL);
		assert(m_ == n_);
		
		for (size_type i=1; i<= m_; i++)
		{
			for (size_type j=1; j<= i; j++)
			{
				(*dasymm)= ((*this)(i,j) - (*this)(j,i))/2;
				dasymm++;
			}
		}
		return 1;		
	}

    int AddSymmPart(pointer dsymm)
	// add Symmetrical matrix in low triangular form 
	// to the current square matrix 
	{
		assert(dsymm != NULL);
		assert(m_ == n_);
		
		for (size_type i=1; i<= m_; i++)
		{
			for (size_type j=1; j<= i; j++)
			{
				(*this)(i,j)+= (*dsymm);
				if(i != j) (*this)(j,i)+= (*dsymm);
				dsymm++;
			}
		}
		return 1;		
	}

    int AddASymmPart(pointer dasymm)
	// add AntiSymmetrical matrix in low triangular form (diagonal included)
	// to the current square matrix 
	{
		assert(dasymm != NULL);
		assert(m_ == n_);
		
		for (size_type i=1; i<= m_; i++)
		{
			for (size_type j=1; j<= i; j++)
			{
				(*this)(i,j)-= (*dasymm);
				(*this)(j,i)+= (*dasymm);
				dasymm++;
			}
		}
		return 1;		
	}

	int Symmetrize()
	// A= 1/2(A+A^T)
	{
		assert(m_ == n_);
		if(m_ == 0)
			return 1;
		for(size_type i=1; i <= m_; i++)
		{
			for (size_type j=1; j< i; j++)
			{
				(*this)(i,j)= 0.5*((*this)(i,j)+ (*this)(j,i));
			}

		}
		return 1;
	}

	int ASymmetrize()
	// A= 1/2(A-A^T)
	{
		assert(m_ == n_);
		if(m_ == 0)
			return 1;
		for(size_type i=1; i <= m_; i++)
		{
			for (size_type j=1; j< i; j++)
			{
				(*this)(i,j)= 0.5*((*this)(i,j) - (*this)(j,i));
				(*this)(j,i)=-(*this)(i,j);
			}

		}
		return 1;
	}

};


/* ***************************  I/O  ********************************/

template <class T>
ostream& operator<<(ostream &s, const Fortran_matrix<T> &A)
{
    size_t M=A.num_rows();
    size_t N=A.num_cols();

    s << M << " " << N << endl;

    for (size_t i=1; i<=M; i++)
    {
        for (size_t j=1; j<=N; j++)
        {
            s << A(i,j) << " ";
        }
        s << endl;
    }


    return s;
}



template <class T>
istream& operator>>(istream &s, Fortran_matrix<T> &A)
{
    size_t M, N;

    s >> M >> N;

    if ( !(M == A.m_ && N == A.n_) )
    {
        A.newsize(M,N);
    }


    for (size_t i=1; i<=M; i++)
        for (size_t j=1; j<=N; j++)
        {
            s >>  A(i,j);
        }

    return s;
}

//*******************[ basic matrix algorithms ]***************************


template <class T>
Fortran_matrix<T> operator+(const Fortran_matrix<T> &A, 
    const Fortran_matrix<T> &B)
{
    size_t M = A.num_rows();
    size_t N = A.num_cols();

    assert(M==B.num_rows());
    assert(N==B.num_cols());

    Fortran_matrix<T> tmp(M,N);
    size_t i,j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
            tmp(i,j) = A(i,j) + B(i,j);

    return tmp;
}

template <class T>
Fortran_matrix<T> operator-(const Fortran_matrix<T> &A, 
    const Fortran_matrix<T> &B)
{
    size_t M = A.num_rows();
    size_t N = A.num_cols();

    assert(M==B.num_rows());
    assert(N==B.num_cols());

    Fortran_matrix<T> tmp(M,N);
    size_t i,j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
            tmp(i,j) = A(i,j) - B(i,j);

    return tmp;
}

// element-wise multiplication  (use matmult() below for matrix
// multiplication in the linear algebra sense.)
//
//
template <class T>
Fortran_matrix<T> mult_element(const Fortran_matrix<T> &A, 
    const Fortran_matrix<T> &B)
{
    size_t M = A.num_rows();
    size_t N = A.num_cols();

    assert(M==B.num_rows());
    assert(N==B.num_cols());

    Fortran_matrix<T> tmp(M,N);
    size_t i,j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
            tmp(i,j) = A(i,j) * B(i,j);

    return tmp;
}


template <class T>
Fortran_matrix<T> transpose(const Fortran_matrix<T> &A)
{
    size_t M = A.num_rows();
    size_t N = A.num_cols();

    Fortran_matrix<T> S(N,M);
    size_t i, j;

    for (i=1; i<=M; i++)
        for (j=1; j<=N; j++)
            S(j,i) = A(i,j);

    return S;
}


    
template <class T>
inline Fortran_matrix<T> matmult(const Fortran_matrix<T>  &A, 
    const Fortran_matrix<T> &B)
{

#ifdef TNT_BOUNDS_CHECK
    assert(A.num_cols() == B.num_rows());
#endif

    size_t M = A.num_rows();
    size_t N = A.num_cols();
    size_t K = B.num_cols();

    Fortran_matrix<T> tmp(M,K);
    T sum;

    for (size_t i=1; i<=M; i++)
    for (size_t k=1; k<=K; k++)
    {
        sum = 0;
        for (size_t j=1; j<=N; j++)
            sum = sum +  A(i,j) * B(j,k);

        tmp(i,k) = sum; 
    }

    return tmp;
}

template <class T>
inline Fortran_matrix<T> operator*(const Fortran_matrix<T> &A, 
    const Fortran_matrix<T> &B)
{
    return matmult(A,B);
}

template <class T>
inline int mat_diff(Fortran_matrix<T>& C, Fortran_matrix<T>  &A, 
    Fortran_matrix<T> &B)
// C= A-B
// C may coincide with A and B
{

    assert(A.num_cols() == B.num_cols());
	assert(A.num_rows() == B.num_rows());

    size_t M = A.num_rows();
    size_t N = A.num_cols();

    C.newsize(M,N);         // adjust shape of C, if necessary

    for(size_t i=1; i <= M; i++)
	{
		for(size_t j=1; j <= N; j++)
		{
			C(i,j)=A(i,j)-B(i,j);
		}
	}
	return 0;
}

template <class T>
inline int mat_add(Fortran_matrix<T>& C, Fortran_matrix<T>  &A, 
    Fortran_matrix<T> &B)
// C= A+B
// C may coincide with A and B
{

    assert(A.num_cols() == B.num_cols());
	assert(A.num_rows() == B.num_rows());

    size_t M = A.num_rows();
    size_t N = A.num_cols();

    C.newsize(M,N);         // adjust shape of C, if necessary

    for(size_t i=1; i <= M; i++)
	{
		for(size_t j=1; j <= N; j++)
		{
			C(i,j)=A(i,j)+B(i,j);
		}
	}
	return 0;
}



template <class T>
inline int mat_add_unit(Fortran_matrix<T>& C, Fortran_matrix<T>  &A, 
    const T& factor)
// C= A + factor*I , where I is a unit matrix
// C may coincide with A 
{
	size_t N=A.num_rows();
	assert(A.num_rows() == A.num_cols());
    C.newsize(N,N);    // adjust shape of C, if necessary
	C=A;
	for(size_t i=1; i<=N; i++)
	{
		C(i,i)=C(i,i)+factor;
	}
	return 0;
}

template <class T>
T SProd(const Fortran_matrix<T>  &A, 
    const Fortran_matrix<T> &B)
// Scalar Product (sum of element pair product of two matricies
{

    assert(A.num_cols() == B.num_cols());
	assert(A.num_rows() == B.num_rows());

    size_t M = A.num_rows();
    size_t N = A.num_cols();
	
	T ss;
    
    ss=0;
    for(size_t i=1; i <= M; i++)
	{
		for(size_t j=1; j <= N; j++)
		{
			ss+=A(i,j)*B(i,j);
		}
	}
	return ss;
}


template <class T>
inline int mat_scale(Fortran_matrix<T>& C, Fortran_matrix<T>  &A, 
    const T& factor)
// C= factor*A  , where I is a unit matrix
// C may coincide with A 
{
    size_t M = A.num_rows();
    size_t N = A.num_cols();

    size_t MN = M*N; 

    C.newsize(M,N);    // adjust shape of C, if necessary

    const T* a = A.begin();
    T* t = C.begin();
    T* tend = C.end();

    for (; t < tend; t++, a++)
        *t = *a * factor;

	return 0;
}

template <class T>
inline int mat_copy_scale(Fortran_matrix<T>& C, const Fortran_matrix<T>  &A, 
    const T& factor)
// C= factor*A  , where I is a unit matrix
// C may coincide with A 
{
    size_t M = A.num_rows();
    size_t N = A.num_cols();

    size_t MN = M*N; 

    C.newsize(M,N);    // adjust shape of C, if necessary

    const T* a = A.begin();
    T* t = C.begin();
    T* tend = C.end();

    for (; t < tend; t++, a++)
        *t = *a * factor;

	return 0;
}

template <class T>
inline int matmult(Fortran_matrix<T>& C, const Fortran_matrix<T>  &A, 
    const Fortran_matrix<T> &B)
{
    assert(A.num_cols() == B.num_rows());

    size_t M = A.num_rows();
    size_t N = A.num_cols();
    size_t K = B.num_cols();

    C.newsize(M,K);         // adjust shape of C, if necessary

    T sum; 

    const T* row_i;
    const T* col_k;

    for (size_t i=1; i<=M; i++)
    {
        for (size_t k=1; k<=K; k++)
        {
            row_i = &A(i,1);
            col_k = &B(1,k);
            sum = 0;
            for (size_t j=1; j<=N; j++)
            {
                sum +=  *row_i * *col_k;
                row_i += M;
                col_k ++;
            }
            C(i,k) = sum; 
        }
    }
    return 0;
}


template <class T>
NumVector<T> matmult(const Fortran_matrix<T>  &A, const NumVector<T> &x)
{

#ifdef TNT_BOUNDS_CHECK
    assert(A.num_cols() == x.dim());
#endif

    size_t M = A.num_rows();
    size_t N = A.num_cols();

    NumVector<T> tmp(M);
    T sum;

    for (size_t i=1; i<=M; i++)
    {
        sum = 0;
        for (size_t j=1; j<=N; j++)
            sum = sum +  A(i,j) * x(j);

        tmp(i) = sum; 
    }

    return tmp;
}

template <class T>
inline NumVector<T> operator*(const Fortran_matrix<T>  &A, const NumVector<T> &x)
{
    return matmult(A,x);
}

template <class T>
inline Fortran_matrix<T> operator*(const Fortran_matrix<T>  &A, const T &x)
{
    size_t M = A.num_rows();
    size_t N = A.num_cols();

    size_t MN = M*N; 

    Fortran_matrix<T> res(M,N);
    const T* a = A.begin();
    T* t = res.begin();
    T* tend = res.end();

    for (t=res.begin(); t < tend; t++, a++)
        *t = *a * x;

    return res;
} 

// my functions

template <class T>
inline int matmult_T1(Fortran_matrix<T>& C, const Fortran_matrix<T>  &A, 
    const Fortran_matrix<T> &B)
// C = A^T * B 
{
    assert(A.num_rows() == B.num_rows());

    size_t M = A.num_cols();
    size_t N = A.num_rows();
    size_t K = B.num_cols();

    C.newsize(M,K);         // adjust shape of C, if necessary

    T sum; 

    const T* col_i;
    const T* col_k;

    for (size_t i=1; i<=M; i++)
    {
        for (size_t k=1; k<=K; k++)
        {
            col_i = &A(1,i);
            col_k = &B(1,k);
            sum = 0;
            for (size_t j=1; j<=N; j++)
            {
                sum +=  *col_i * *col_k;
                col_i ++;
                col_k ++;
            }
        
            C(i,k) = sum; 
        }

    }

    return 0;
}

template <class T>
inline int matmult_T2(Fortran_matrix<T>& C, const Fortran_matrix<T>  &A, 
    const Fortran_matrix<T> &B)
// C = A * B^T
{

    assert(A.num_cols() == B.num_cols());

    size_t M = A.num_rows();
    size_t N = A.num_cols();
    size_t K = B.num_rows();

    C.newsize(M,K);         // adjust shape of C, if necessary


    T sum; 

    const T* row_i;
    const T* row_k;

    for (size_t i=1; i<=M; i++)
    {
        for (size_t k=1; k<=K; k++)
        {
            row_i = &A(i,1);
            row_k = &B(k,1);
            sum = 0;
            for (size_t j=1; j<=N; j++)
            {
                sum +=  *row_i * *row_k;
                row_i += M;
                row_k += K;
            }
        
            C(i,k) = sum; 
        }

    }

    return 0;
}

template<class NumMatrix, class NumVector>
int mat_symm_from_lin(NumMatrix & smat, const NumVector & slin, size_t n)
// convert Symmetrical Matrix from linear storage to Square form
// ( indexation of the matrix and the vector are assumed to be 1-based )
{
	smat.newsize(n,n);
	for(size_t i=1; i <= n ; i++)
	{
		for(size_t j=1; j <= n ; j++)
		{
			double tmp = ((i >= j) ? slin(i*(i-1)/2 + j) : slin(j*(j-1)/2 + i));
	        smat.r1(i,j)= tmp;
		}
	}
	return 1;
}

template<class NumMatrix, class NumVector>
int mat_asymm_from_lin(NumMatrix & smat, NumVector & slin, size_t n)
// convert AntiSymmetrical Matrix from linear storage to Square form
// ( indexation of the matrix and the vector are assumed to be 1-based )
{
	smat.newsize(n,n);
	for(size_t i=1; i <= n ; i++)
	{
		for(size_t j=1; j <= n ; j++)
		{
	      smat(i,j)=((i >= j) ? slin(i*(i-1)/2 + j) : -slin(j*(j-1)/2 + i));
		}
	}
	return 1;
}

#endif
// FMAT_H
