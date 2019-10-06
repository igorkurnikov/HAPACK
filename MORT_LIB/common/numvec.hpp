#ifndef MORT_OBJECT_NUMVEC_HPP
#define MORT_OBJECT_NUMVEC_HPP

#include <boost/numeric/ublas/vector.hpp>

namespace mort
{
    typedef boost::numeric::ublas::vector<double> numvec;
    typedef boost::numeric::ublas::zero_vector<double> zero_numvec;
    typedef boost::numeric::ublas::scalar_vector<double> scalar_numvec;

    numvec makevec(size_t, const double*);

    /// construct numvec from elements
    numvec makevec(const double& x, const double& y);

    numvec makevec(const double& x, const double& y, const double& z);

    numvec makevec(const double& x, const double& y, const double& z,
                   const double& s);

    numvec subvec(const numvec&, size_t b, size_t n);

    // maximum element of a numvec
    double max(const numvec& vec);

    /// norm of a vector 
    double norm(const numvec& vec);

    /// dot product of two vectors.
    double dotpd(const numvec& v1, const numvec& v2);

    /// cross product of two vector
    numvec cross(const numvec& v1, const numvec& v2);

    /// normalize a vector in place
    numvec& normalize(numvec& vec);

    /// return a normalized copy of a vector, the orginal one remain unchanged.
    numvec  normalcpy(const numvec& v);

} // namespace mort

#endif // MORT_OBJECT_NUMVEC_HPP
