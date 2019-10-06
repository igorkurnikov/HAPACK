#include <cmath>
#include <cassert>
#include "fortran.hpp"
#include "numvec.hpp"

namespace mort
{
    numvec makevec(size_t n, const double* x)
    {
        numvec result(n);
        std::copy(x, x + n, result.begin());
        return result;
    }

    numvec subvec(const numvec& src, size_t b, size_t n)
    {
        assert(n > 0);
        assert(b + n <= src.size());
        numvec result(n);
        std::copy(src.begin() + b, src.begin() + b + n, result.begin());
        return result;
    }

    numvec makevec(const double& x, const double& y)
    {
        numvec result(2);
        result[0] = x;
        result[1] = y;
        return result;
    }

    numvec makevec(const double& x, const double& y, const double& z)
    {
        numvec result(3);
        result[0] = x;
        result[1] = y;
        result[2] = z;
        return result;
    }

    numvec makevec(const double& x, const double& y, const double& z,
                   const double& s)
    {
        numvec result(4);
        result[0] = x;
        result[1] = y;
        result[2] = z;
        result[3] = s;
        return result;
    }

    double max(const numvec& v)
    {
        assert(v.size() > 0);

        double accu(v[0]);
        for (numvec::const_iterator i = v.begin() + 1; i != v.end(); ++i)
            accu = std::max(accu, *i);

        return accu;
    }

    double norm(const numvec& v)
    {
        assert(v.size() > 0);

        double accu(0);
        for (numvec::const_iterator i = v.begin(); i != v.end(); ++i)
            accu += (*i)*(*i);

        return std::sqrt(accu);
    }

    double dotpd(const numvec& a, const numvec& b)
    {
        assert(a.size() == b.size());

        double accu(0);
        for (size_t n(0); n < a.size(); ++n)
            accu += a[n]*b[n];

        return accu;
    }

    numvec cross(const numvec& v1, const numvec& v2)
    {
        return makevec(v1[1]*v2[2] - v2[1]*v1[2],
                       v1[2]*v2[0] - v2[2]*v1[0],
                       v1[0]*v2[1] - v2[0]*v1[1]);
    }

    numvec& normalize(numvec& vec)
    {
        return (vec /= norm(vec));
    }

    numvec normalcpy(const numvec& vec)
    {
        numvec copy(vec);
        return normalize(copy);
    }

} // namespace mort
