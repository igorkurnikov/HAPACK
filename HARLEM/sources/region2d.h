// Template Numerical Toolkit (TNT) for Linear Algebra
//
// BETA VERSION INCOMPLETE AND SUBJECT TO CHANGE
// Please see http://math.nist.gov/tnt for updates
//
// R. Pozo
// Mathematical and Computational Sciences Division
// National Institute of Standards and Technology

// 2D Regions for arrays and matrices

#ifndef REGION2D_H
#define REGION2D_H

#include "index.h"

// const_Region2D needs the T parameter -- cannot rely on 
// TNT_Array2D::element_type for definition later. (possible compiler
// g++ bug?)

template <class TNT_Array2D, class T > class const_Region2D;

template <class TNT_Array2D /*,  class T */>
class Region2D
{
    protected:

        TNT_Array2D &  A_;
        size_t offset_[2];       // 0-offset internally
        size_t dim_[2];

    public:
	typedef typename TNT_Array2D::value_type T;
	//        typedef TNT_Array2D::value_type T;
        typedef size_t   size_type;
        typedef         T   value_type;
        typedef         T   element_type;
        typedef         T*  pointer;
        typedef         T*  iterator;
        typedef         T&  reference;
        typedef const   T*  const_iterator;
        typedef const   T&  const_reference;

        TNT_Array2D & array() { return A_; }
        const TNT_Array2D & array()  const { return A_; }
        size_t lbound() const { return A_.lbound(); }
        size_t num_rows() const { return dim_[0]; }
        size_t num_cols() const { return dim_[1]; }
        size_t offset(size_t i) const                 // 1-offset
        {
#ifdef TNT_BOUNDS_CHECK
            assert( A_.lbound() <= i);
            assert( i<= dim_[0] + A_.lbound()-1);
#endif
            return offset_[i-A_.lbound()];
        }

        size_t dim(size_t i) const
        {
#ifdef TNT_BOUNDS_CHECK
            assert( A_.lbound() <= i);
            assert( i<= dim_[0] + A_.lbound()-1);
#endif
            return dim_[i-A_.lbound()];
        }



        Region2D(TNT_Array2D &A, size_t i1, size_t i2, size_t j1,
                size_t j2) : A_(A)
        {
#ifdef TNT_BOUNDS_CHECK
            assert( i1 <= i2 );
            assert( j1 <= j2);
            assert( A.lbound() <= i1);
            assert( i2<= A.dim(A.lbound()) + A.lbound()-1);
            assert( A.lbound() <= j1);
            assert( j2<= A.dim(A.lbound()+1) + A.lbound()-1 );
#endif


            offset_[0] = i1-A.lbound();
            offset_[1] = j1-A.lbound();
            dim_[0] = i2-i1+1;
            dim_[1] = j2-j1+1;
        }

        Region2D(TNT_Array2D &A, const Index1D &I, const Index1D &J) : A_(A)
        {
#ifdef TNT_BOUNDS_CHECK
            assert( I.lbound() <= I.ubound() );
            assert( J.lbound() <= J.ubound() );
            assert( A.lbound() <= I.lbound());
            assert( I.ubound()<= A.dim(A.lbound()) + A.lbound()-1);
            assert( A.lbound() <= J.lbound());
            assert( J.ubound() <= A.dim(A.lbound()+1) + A.lbound()-1 );
#endif

            offset_[0] = I.lbound()-A.lbound();
            offset_[1] = J.lbound()-A.lbound();
            dim_[0] = I.ubound() - I.lbound() + 1;
            dim_[1] = J.ubound() - J.lbound() + 1;
        }

        Region2D(Region2D<TNT_Array2D> &A, size_t i1, size_t i2,
            size_t j1, size_t j2) : A_(A.A_)
        {
#ifdef TNT_BOUNDS_CHECK
            assert( i1 <= i2 );
            assert( j1 <= j2);
            assert( A.lbound() <= i1);
            assert( i2<= A.dim(A.lbound()) + A.lbound()-1);
            assert( A.lbound() <= j1);
            assert( j2<= A.dim(A.lbound()+1) + A.lbound()-1 );
#endif
            offset_[0] = (i1 - A.lbound()) + A.offset_[0];
            offset_[1] = (j1 - A.lbound()) + A.offset_[1];
            dim_[0] = i2-i1 + 1;
            dim_[1] = j2-j1+1;
        }

        Region2D<TNT_Array2D> operator()(size_t i1, size_t i2,
                size_t j1, size_t j2)
        {
#ifdef TNT_BOUNDS_CHECK
            assert( i1 <= i2 );
            assert( j1 <= j2);
            assert( A_.lbound() <= i1);
            assert( i2<= dim_[0] + A_.lbound()-1);
            assert( A_.lbound() <= j1);
            assert( j2<= dim_[1] + A_.lbound()-1 );
#endif
            return Region2D<TNT_Array2D>(A_, 
                    i1+offset_[0], offset_[0] + i2, 
                    j1+offset_[1], offset_[1] + j2);
        }


        Region2D<TNT_Array2D> operator()(const Index1D &I,
                const Index1D &J)
        {
#ifdef TNT_BOUNDS_CHECK
            assert( I.lbound() <= I.ubound() );
            assert( J.lbound() <= J.ubound() );
            assert( A_.lbound() <= I.lbound());
            assert( I.ubound()<= dim_[0] + A_.lbound()-1);
            assert( A_.lbound() <= J.lbound());
            assert( J.ubound() <= dim_[1] + A_.lbound()-1 );
#endif

            return Region2D<TNT_Array2D>(A_, I.lbound()+offset_[0],
                offset_[0] + I.ubound(), offset_[1]+J.lbound(),
                offset_[1] + J.ubound());
        }

        inline T & operator()(size_t i, size_t j)
        {
#ifdef TNT_BOUNDS_CHECK
            assert( A_.lbound() <= i);
            assert( i<= dim_[0] + A_.lbound()-1);
            assert( A_.lbound() <= j);
            assert( j<= dim_[1] + A_.lbound()-1 );
#endif
            return A_(i+offset_[0], j+offset_[1]);
        }

        inline const T & operator() (size_t i, size_t j) const
        {
#ifdef TNT_BOUNDS_CHECK
            assert( A_.lbound() <= i);
            assert( i<= dim_[0] + A_.lbound()-1);
            assert( A_.lbound() <= j);
            assert( j<= dim_[1] + A_.lbound()-1 );
#endif
            return A_(i+offset_[0], j+offset_[1]);
        }


        Region2D<TNT_Array2D> & operator=(const Region2D<TNT_Array2D> &R)
        {
            size_t M = num_rows(); 
            size_t N = num_cols();

            // make sure both sides conform
            assert(M == R.num_rows());
            assert(N == R.num_cols());


            size_t start = R.lbound();
            size_t Mend =  start + M - 1;
            size_t Nend =  start + N - 1;
            for (size_t i=start; i<=Mend; i++)
              for (size_t j=start; j<=Nend; j++)
                (*this)(i,j) = R(i,j);

            return *this;
        }

        Region2D<TNT_Array2D> & operator=(const const_Region2D<TNT_Array2D,T> &R)
        {
            size_t M = num_rows(); 
            size_t N = num_cols();

            // make sure both sides conform
            assert(M == R.num_rows());
            assert(N == R.num_cols());


            size_t start = R.lbound();
            size_t Mend =  start + M - 1;
            size_t Nend =  start + N - 1;
            for (size_t i=start; i<=Mend; i++)
              for (size_t j=start; j<=Nend; j++)
                (*this)(i,j) = R(i,j);

            return *this;
        }

        Region2D<TNT_Array2D> & operator=(const TNT_Array2D &R)
        {
            size_t M = num_rows(); 
            size_t N = num_cols();

            // make sure both sides conform
            assert(M == R.num_rows());
            assert(N == R.num_cols());


            size_t start = R.lbound();
            size_t Mend =  start + M - 1;
            size_t Nend =  start + N - 1;
            for (size_t i=start; i<=Mend; i++)
              for (size_t j=start; j<=Nend; j++)
                (*this)(i,j) = R(i,j);

            return *this;
        }

        Region2D<TNT_Array2D> & operator=(const  T &scalar)
        {
            size_t start = lbound();
            size_t Mend = lbound() + num_rows() - 1;
            size_t Nend = lbound() + num_cols() - 1;


            for (size_t i=start; i<=Mend; i++)
              for (size_t j=start; j<=Nend; j++)
                (*this)(i,j) = scalar;

            return *this;
        }

};

//************************

template <class TNT_Array2D, class T>
class const_Region2D
{
    protected:

        const TNT_Array2D &  A_;
        size_t offset_[2];       // 0-offset internally
        size_t dim_[2];

    public:

        typedef         T   value_type;
        typedef         T   element_type;
        typedef const   T*  const_iterator;
        typedef const   T&  const_reference;

        const TNT_Array2D & array() const { return A_; }
        size_t lbound() const { return A_.lbound(); }
        size_t num_rows() const { return dim_[0]; }
        size_t num_cols() const { return dim_[1]; }
        size_t offset(size_t i) const                 // 1-offset
        {
#ifdef TNT_BOUNDS_CHECK
            assert( TNT_BASE_OFFSET <= i);
            assert( i<= dim_[0] + TNT_BASE_OFFSET-1);
#endif
            return offset_[i-TNT_BASE_OFFSET];
        }

        size_t dim(size_t i) const
        {
#ifdef TNT_BOUNDS_CHECK
            assert( TNT_BASE_OFFSET <= i);
            assert( i<= dim_[0] + TNT_BASE_OFFSET-1);
#endif
            return dim_[i-TNT_BASE_OFFSET];
        }


        const_Region2D(const TNT_Array2D &A, size_t i1, size_t i2, 
                size_t j1, size_t j2) : A_(A)
        {
#ifdef TNT_BOUNDS_CHECK
            assert( i1 <= i2 );
            assert( j1 <= j2);
            assert( TNT_BASE_OFFSET <= i1);
            assert( i2<= A.dim(TNT_BASE_OFFSET) + TNT_BASE_OFFSET-1);
            assert( TNT_BASE_OFFSET <= j1);
            assert( j2<= A.dim(TNT_BASE_OFFSET+1) + TNT_BASE_OFFSET-1 );
#endif

            offset_[0] = i1-TNT_BASE_OFFSET;
            offset_[1] = j1-TNT_BASE_OFFSET;
            dim_[0] = i2-i1+1;
            dim_[1] = j2-j1+1;
        }

        const_Region2D(const TNT_Array2D &A, const Index1D &I, const Index1D &J) 
                : A_(A)
        {
#ifdef TNT_BOUNDS_CHECK
            assert( I.lbound() <= I.ubound() );
            assert( J.lbound() <= J.ubound() );
            assert( TNT_BASE_OFFSET <= I.lbound());
            assert( I.ubound()<= A.dim(TNT_BASE_OFFSET) + TNT_BASE_OFFSET-1);
            assert( TNT_BASE_OFFSET <= J.lbound());
            assert( J.ubound() <= A.dim(TNT_BASE_OFFSET+1) + TNT_BASE_OFFSET-1 );
#endif

            offset_[0] = I.lbound()-TNT_BASE_OFFSET;
            offset_[1] = J.lbound()-TNT_BASE_OFFSET;
            dim_[0] = I.ubound() - I.lbound() + 1;
            dim_[1] = J.ubound() - J.lbound() + 1;
        }


        const_Region2D(const_Region2D<TNT_Array2D,T> &A, size_t i1, 
                size_t i2,
            size_t j1, size_t j2) : A_(A.A_)
        {
#ifdef TNT_BOUNDS_CHECK
            assert( i1 <= i2 );
            assert( j1 <= j2);
            assert( TNT_BASE_OFFSET <= i1);
            assert( i2<= A.dim(TNT_BASE_OFFSET) + TNT_BASE_OFFSET-1);
            assert( TNT_BASE_OFFSET <= j1);
            assert( j2<= A.dim(TNT_BASE_OFFSET+1) + TNT_BASE_OFFSET-1 );
#endif
            offset_[0] = (i1 - TNT_BASE_OFFSET) + A.offset_[0];
            offset_[1] = (j1 - TNT_BASE_OFFSET) + A.offset_[1];
            dim_[0] = i2-i1 + 1;
            dim_[1] = j2-j1+1;
        }

        const_Region2D<TNT_Array2D,T> operator()(size_t i1, size_t i2,
                size_t j1, size_t j2)
        {
#ifdef TNT_BOUNDS_CHECK
            assert( i1 <= i2 );
            assert( j1 <= j2);
            assert( TNT_BASE_OFFSET <= i1);
            assert( i2<= dim_[0] + TNT_BASE_OFFSET-1);
            assert( TNT_BASE_OFFSET <= j1);
            assert( j2<= dim_[0] + TNT_BASE_OFFSET-1 );
#endif
            return const_Region2D<TNT_Array2D,T>(A_, 
                    i1+offset_[0], offset_[0] + i2, 
                    j1+offset_[1], offset_[1] + j2);
        }


        const_Region2D<TNT_Array2D,T> operator()(const Index1D &I,
                const Index1D &J)
        {
#ifdef TNT_BOUNDS_CHECK
            assert( I.lbound() <= I.ubound() );
            assert( J.lbound() <= J.ubound() );
            assert( TNT_BASE_OFFSET <= I.lbound());
            assert( I.ubound()<= dim_[0] + TNT_BASE_OFFSET-1);
            assert( TNT_BASE_OFFSET <= J.lbound());
            assert( J.ubound() <= dim_[1] + TNT_BASE_OFFSET-1 );
#endif

            return const_Region2D<TNT_Array2D,T>(A_, I.lbound()+offset_[0],
                offset_[0] + I.ubound(), offset_[1]+J.lbound(),
                offset_[1] + J.ubound());
        }


        inline const T & operator() (size_t i, size_t j) const
        {
#ifdef TNT_BOUNDS_CHECK
            assert( TNT_BASE_OFFSET <= i);
            assert( i<= dim_[0] + TNT_BASE_OFFSET-1);
            assert( TNT_BASE_OFFSET <= j);
            assert( j<= dim_[1] + TNT_BASE_OFFSET-1 );
#endif
            return A_(i+offset_[0], j+offset_[1]);
        }

};


//  ************** ostream algorithms *******************************

template <class TNT_Array2D,class T>
std::ostream& operator<<(std::ostream &s, const const_Region2D<TNT_Array2D,T> &A)
{
    size_t start = A.lbound();
    size_t Mend=A.lbound()+ A.num_rows() - 1;
    size_t Nend=A.lbound() + A.num_cols() - 1;


    s << A.num_rows() << "  " << A.num_cols() << std::endl;
    for (size_t i=start; i<=Mend; i++)
    {
        for (size_t j=start; j<=Nend; j++)
        {
            s << A(i,j) << " ";
        }
        s << std::endl;
    }


    return s;
}

template <class TNT_Array2D>
std::ostream& operator<<(std::ostream &s, const Region2D<TNT_Array2D> &A)
{
    size_t start = A.lbound();
    size_t Mend=A.lbound()+ A.num_rows() - 1;
    size_t Nend=A.lbound() + A.num_cols() - 1;


    s << A.num_rows() << "  " << A.num_cols() << std::endl;
    for (size_t i=start; i<=Mend; i++)
    {
        for (size_t j=start; j<=Nend; j++)
        {
            s << A(i,j) << " ";
        }
        s << std::endl;
    }


    return s;

}

#endif
// REGION2D_H
