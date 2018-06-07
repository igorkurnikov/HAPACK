/*! \file vec.h

   Numerical vector classes 
   Derived from 
   Template Numerical Toolkit (TNT) for Linear Algebra
 
   BETA VERSION INCOMPLETE AND SUBJECT TO CHANGE
   Please see http://math.nist.gov/tnt for updates
 
   R. Pozo
   Mathematical and Computational Sciences Division
   National Institute of Standards and Technology
  
   Basic TNT  numerical vector (0-based [i] AND (i) indexing )
*/

#ifndef VEC_H
#define VEC_H

#include "haconst.h"
#include <stdlib.h>
#include <assert.h>

enum AllocMode {INTERNAL_ALLOC=0, EXTERNAL_ALLOC=1}; 
	// define type of storage for the vector 
    // v_ is allocated when matrix is 
	// created and is destroyed when vector is destroyed or   
    // v points to externally allocated memory region assumed 
	// to be of the sufficient size (m_*n_)

template <class T>
class NumVector
//! Numerical vector class - provides both 0- and 1- based index access
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

    size_type lbound() const { return 1;}
	
  protected:
    T* v_;                //!< pointer for allocated vector  
    T* vm1_;              //!< pointer adjustment (v_-1) for optimized 1-offset indexing
    size_type n_;         //!< size of the array 
    size_type nalloc_;    //!< maximal number of elements in the allocated memory 
	                      
    AllocMode amode;   //!< Allocation model (allocate memory or use supplied pointer to the memory block)
    
    void initialize(size_type N)
//! internal function for initialization internal allocaction array of size N
    {
        // adjust pointers so that they are 1-offset:
        // v_[] is the internal contiguous array, it is still 0-offset
        //
		amode = INTERNAL_ALLOC;
		v_   = NULL;
		vm1_ = NULL;
		n_      = 0;
		nalloc_ = 0;
		reserve(N);
		resize(N);
    }

  public:

    void reserve(size_type N)
	//!< Realocate the vector 
	{
		if( amode == EXTERNAL_ALLOC) 
		{
			nalloc_ = N;
			return;
		}
			
		T* v_new = (T*) realloc(v_, N*sizeof(T));
		if( v_new == NULL )
		{
			throw "Error in NumVector::reserve() \n";
		}

		v_ = v_new; 
		vm1_ = v_-1;
		nalloc_ = N;
	}


    void resize(size_type N)
    {
        if (n_ == N) return;

		if(N <= nalloc_ || (amode == EXTERNAL_ALLOC && nalloc_ == 0)  ) 
		{
			n_ = N;
			return;
		}

        reserve(N);
		n_ = N;
    }

    void newsize(size_type N)
    {
		resize(N);
    }

    void copy(const T*  v)
    {
        size_type N = n_;
        size_type i;

        for (i=0; i< N; i++)
            v_[i] = v[i];
    }

    void set(const T& val)
    {
        size_type N = n_;
        size_type i;

        for (i=0; i< N; i++)
            v_[i] = val;      
    }

	void push_back(const T& x)
	{
		n_++;
		if( n_ <= nalloc_)
		{
			vm1_[n_] = x;
		}
		else
		{
			size_type nalloc_new;
			if( nalloc_ < 4 )
			{
				nalloc_new = 4;
			}
			else
			{
				nalloc_new = nalloc_ + nalloc_/2;
			}
			reserve(nalloc_new);
			vm1_[n_] = x;
		}
	}

	inline void append( const T& x)
	{	
		push_back(x);
	}
    
	bool operator == (const NumVector<T>& rhs ) const	
	{
		if(n_ != rhs.n_) return false;
		size_type i;
		for(i = 0; i < n_; i++)
		{
			if(v_[i] != rhs.v_[i]) return false;
		}
		return true;
	}

	bool operator < (const NumVector<T>& rhs) const
	{
		if(n_ > rhs.n_) return false;
		if(n_ < rhs.n_) return true;
		size_type i;
		for(i = 0; i < n_; i++)
		{
			if(v_[i] > rhs.v_[i]) return false;
			if(v_[i] < rhs.v_[i]) return true;
		}
		return false;
	}

    void destroy()
    {
      if( amode == EXTERNAL_ALLOC) 
	  {
		  v_   = NULL;
		  vm1_ = NULL;
		  n_   = 0;
		  nalloc_ = 0;
		  amode = INTERNAL_ALLOC;
		  return;
	  }
		  
	  if( v_ != NULL )
	  {
	 	  free(v_);
	  }
	  
	  n_ = 0;
	  nalloc_ = 0;
      v_ = NULL;
      vm1_ = NULL;
    }

	void clear()
	{ 
		destroy();
	}


  public:

    // access

    iterator v() { return v_;}
    iterator begin() { return v_;}
    iterator end()   { return v_ + n_; }
	const_iterator v() const { return v_;}
    const_iterator begin() const { return v_;}
    const_iterator end() const  { return v_ + n_; }

    // destructor

    virtual ~NumVector() 
    {
        destroy();
    }

    // constructors

    NumVector() : v_(0), vm1_(0), n_(0), nalloc_(0), amode(INTERNAL_ALLOC)  {}

    NumVector(const NumVector<T> &A) : v_(0), vm1_(0), n_(0), nalloc_(0), amode(INTERNAL_ALLOC)
    {
        initialize(A.n_);
        copy(A.v_);
    }

    NumVector(size_type N, const T& value = T(0)) :  v_(0), vm1_(0), n_(0), nalloc_(0), amode(INTERNAL_ALLOC)
    {
        initialize(N);
        set(value);
    }

    NumVector(size_type N, const T* v) :  v_(0), vm1_(0), n_(0), nalloc_(0), amode(INTERNAL_ALLOC)
    {
        initialize(N);
        copy(v);
    }

    NumVector(size_type N, char *s) :  v_(0), vm1_(0), n_(0), nalloc_(0), amode(INTERNAL_ALLOC)
    {
       initialize(N);
       istrstream ins(s);

       size_type i;

       for (i=0; i<N; i++)
                ins >> v_[i];
    }

//>Mikola Jul 27,2006 QUESTIONABLE CHECK: should be changed
 //    NumVector(T* v,size_type N) :  v_(0), vm1_(0), n_(0), nalloc_(0), amode(INTERNAL_ALLOC)
 //    {
 //       assert(v != NULL);
 //       initialize_array_ctrl(N,v); 
 //   }

//    void initialize_array_ctrl(size_type N, T* v)
 //   {
        // adjust pointers so that they are 1-offset:
        // v_[] is the internal contiguous array, it is still 0-offset
      //
 //     assert(v_ == NULL);
 //     if(N==0) 
 //       return;
//        v_ = new T[N];
      //v_=(pointer)malloc(N*sizeof(T));
//      v_=v;
//      assert(v_  != NULL);
//      vm1_ = v_-1;
//      n_ = N;
//	  nalloc_ = N;
//    }

	void set_ext_alloc(T* v, size_type N, size_type N_alloc = 0)
	{
		if(v == NULL) 
		{
			throw("Error in NumVector<T>::set_ext_alloc() v == NULL \n");
		}
		if(v_ != 0) destroy();
		amode = EXTERNAL_ALLOC;
		v_ = v;
		vm1_ = v_-1;
		n_ = N;
		nalloc_ = N_alloc;
	}

    
    NumVector<T>& operator=(const NumVector<T> &A)
    {
        if (v_ == A.v_)
            return *this;

        if (n_ == A.n_) // no need to re-alloc
		{
			copy(A.v_);
		}
        else
        {
			resize(A.n_);
            copy(A.v_);
        }
        return *this;
    }
       
    NumVector<T>& operator=(const T& scalar)
    { 
        set(scalar);  
        return *this;
    }
	
	NumVector<T> operator*(const T& scalar)
	{
		NumVector<T> tmp(*this);
		tmp.scale(scalar);
		return tmp;
	}

	NumVector<T>& scale(const T& scalar)
    { 
		size_type i;
        for(i= 0; i < n_;i++)
			v_[i] *= scalar;
        return *this;
    }
	
	NumVector<T>& operator*=(const T& scalar)
    { 
        return scale(scalar);
    }
	
	NumVector<T>& operator+=(const NumVector<T> &A)
	{
		size_type i;
		for (i=0; i< n_; i++)
            (*this)[i] += A[i];

		return *this;
	}

    size_type dim() const
//!< number of elements in the array
    {
        return  n_; 
    }

    size_type size() const 
//!< number of elements in the array
	{
        return  n_; 
    }

	size_type capacity() const
//! return alocated memory of the array (in number of elements) 
	{
		return nalloc_;
	}

	inline value_type GetVal_idx0(size_type i) const
//! 0- based access to the vector elements
	{
#ifdef TNT_BOUNDS_CHECK
        assert(i < n_) ;
#endif
        return v_[i]; 
    }

	inline void SetVal_idx0(size_type i, const value_type new_val)
//! 0- based access to the vector elements
	{
#ifdef TNT_BOUNDS_CHECK
        assert(i < n_) ;
#endif
		v_[i] = new_val;
	}

	inline value_type GetVal_idx1(size_type i) const
//! 1- based access to the vector elements
	{
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= n_) ;
#endif
        return vm1_[i]; 
    }

	inline void SetVal_idx1(size_type i, const value_type new_val)
//! 1- based access to the vector elements
	{
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= n_) ;
#endif
		vm1_[i] = new_val;
	}

	inline value_type GetVal(size_type i) const
	{
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= n_) ;
#endif
        return vm1_[i]; 
    }

	inline void SetVal(size_type i, const value_type new_val)
	{
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= n_) ;
#endif
		vm1_[i] = new_val;
	}

    inline reference operator()(size_type i)
//! 1- based access to the vector elements
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= n_) ;
#endif
        return vm1_[i]; 
    }

    inline const_reference operator() (size_type i) const
 //! 1- based access to the vector elements
	{
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= n_) ;
#endif
        return vm1_[i]; 
    }

    inline reference r1(size_type i)
//! 1- based access to the vector elements
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= n_) ;
#endif
        return vm1_[i]; 
    }

    inline value_type r1(size_type i) const
 //! 1- based access to the vector elements
	{
#ifdef TNT_BOUNDS_CHECK
        assert(1<=i);
        assert(i <= n_) ;
#endif
        return vm1_[i]; 
    }


    inline reference operator[](size_type i)
//! 0- based access to the vector elements
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(i < n_) ;
#endif
        return v_[i]; 
    }

    inline const_reference operator[](size_type i) const
//! 0- based access to the vector elements
	{
#ifdef TNT_BOUNDS_CHECK
        assert(i < n_) ;
#endif
        return v_[i]; 
    }

	inline reference r0(size_type i)
//! 0- based access to the vector elements
    { 
#ifdef TNT_BOUNDS_CHECK
        assert(i < n_) ;
#endif
        return v_[i]; 
    }

    inline value_type r0(size_type i) const
//! 0- based access to the vector elements
	{
#ifdef TNT_BOUNDS_CHECK
        assert(i < n_) ;
#endif
        return v_[i]; 
    }

//    friend istream& operator>><>(istream &s, NumVector<T> &A);

};


/* ***************************  I/O  ********************************/

template <class T>
ostream& operator<<(ostream &s, const NumVector<T> &A)
{
    size_t N=A.dim();

    s <<  N << endl;

    for (size_t i=0; i<N; i++)
        s   << A[i] << " " << endl;
    s << endl;

    return s;
}

template <class T>
istream& operator>>(istream &s, NumVector<T> &A)
{
    size_t N;

    s >> N;

    if ( !(N == A.n_) )
    {
        A.destroy();
        A.initialize(N);
    }


    for (size_t i=0; i<N; i++)
            s >>  A[i];


    return s;
}

//*******************[ basic matrix algorithms ]***************************


template <class T>
NumVector<T> operator+(const NumVector<T> &A, 
    const NumVector<T> &B)
{
    size_t N = A.dim();

    if( N !=B.size()) throw "Error in NumVector<T> operator+() N != B.size() \n";

    NumVector<T> tmp(N);
    size_t i;

    for (i=0; i<N; i++)
        tmp[i] = A[i] + B[i];

    return tmp;
}

template <class T>
NumVector<T> operator-(const NumVector<T> &A, 
    const NumVector<T> &B)
{
    size_t N = A.size();

    if( N !=B.size()) throw "Error in NumVector<T> operator-() N != B.size() \n";

    NumVector<T> tmp(N);
    size_t i;

    for (i=0; i<N; i++)
         tmp[i] = A[i] - B[i];

    return tmp;
}

template <class T>
NumVector<T> operator*(const NumVector<T> &A, 
    const NumVector<T> &B)
{
    size_t N = A.size();

    if( N !=B.size()) throw "Error in NumVector<T> operator*() N != B.size() \n";

    NumVector<T> tmp(N);
    size_t i;

    for (i=0; i<N; i++)
            tmp[i] = A[i] * B[i];

    return tmp;
}

template <class T>
NumVector<T> operator/(const NumVector<T> &A, 
    const NumVector<T> &B)
{
    size_t N = A.size();

    if( N !=B.size()) throw "Error in NumVector<T> operator/() N != B.size() \n";

    NumVector<T> tmp(N);
    size_t i;

    for (i=0; i<N; i++)
            tmp[i] = A[i] / B[i];

    return tmp;
}

template <class T>
int operator == (const NumVector<T> &A, 
    const T & t)
{
    size_t N = A.size();
    size_t i;

    for (i=0; i<N; i++)
	{
       if( A[i] != t) return 0;
	}

    return 1;
}

template <class T>
T dot_product(const NumVector<T> &A, const NumVector<T> &B)
{
    size_t N = A.size();
    if( N !=B.size()) throw "Error in dot_product() N != B.size() \n";

    size_t i;
    T sum = 0;

    for (i=0; i<N; i++)
        sum += A[i] * B[i];

    return sum;
}

#endif /* if(!Defined(VEC_H))  */
 


