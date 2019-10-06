#ifndef MORTSRC_COMMON_MATRIX_HPP
#define MORTSRC_COMMON_MATRIX_HPP

#include <valarray>
#include "numvec.hpp"

namespace mort
{

    class matrix
    {
    public:

	matrix(int n1, int n2);

        matrix(const matrix& rhs);

        double& operator()( int i, int j );
        
        double operator()(int i, int j ) const;
        
        int size1() const 
        {
            return m_size1;
        }
        
        int size2() const
        {
            return m_size2;
        }
        
    private:

        int m_size1;
        int m_size2;
        std::valarray<double> m_data;
    };

    bool is_identity( const double* mat, int size );

    // jacobi decomposition of matrix
    void jacobi( matrix& m, numvec& eigvals, matrix& eigvecs);

} // namespace mort


#endif
