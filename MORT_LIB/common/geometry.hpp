#ifndef MORTSRC_COMMON_GEOMETRY_HPP
#define MORTSRC_COMMON_GEOMETRY_HPP

#include "numvec.hpp"
#include "matrix.hpp"

namespace mort
{
    /// distance square between 2 points
    double dis2( const numvec& v1, const numvec& v2 );

    /// distance between 2 points
    double dist( const numvec& v1, const numvec& v2 );

    /// angle between 3 points, v1-v2-v3
    double angl( const numvec& v1, const numvec& v2, const numvec& v3 );
 
    /// dihedral angle between 4 points (v1-v2-v3-v4)
    double tors( const numvec& v1, const numvec& v2, const numvec& v3, const numvec& v4 );

    /// in place rotation, ang around axis
    numvec& rotate( numvec& vec, const numvec& axs, double ang );

    /// in place rotation, on a matrix
    double* rotate( double* vec, const matrix& mat );

    /// rotate a vector and return a copy, do not change the original
    numvec rotcpy( const numvec& vec, numvec axs, double theta );

    /// superimpose the pts in srcpot in pts in dstpos, apply the transform to input
    numvec impose( std::vector< numvec > srcpos, const std::vector< numvec >& dstpos, const numvec& input );
    
    /// return the plane determined by the 3 points. 
    /// the result is a 4-element vector (n1, n2, n3, d)
    /// the plane is then n1*x + n2*y + n3*z = d
    numvec plane( const numvec& v1, const numvec& v2, const numvec& v3 );

    // check if a set of points is on a line
    bool same_line( const std::vector<numvec>& pts );
 
    // check if a point is on a plane p.  
    bool same_plane( const numvec& c, const numvec& p ); 

    // check if a set of points is on a plane, if so, p will be set to the plane
    bool same_plane( const std::vector<numvec>& crds );

    // check if a set of points is on a plane, if so, p will be set to the plane
    bool same_plane( const std::vector<numvec>& crds, numvec& p );

    numvec& rotate_x( numvec& src, double ang );

    numvec& rotate_y( numvec& src, double ang );

    numvec& rotate_z( numvec& src, double ang );

    // zmatrix routines: given a zmatrix, give the result.
    numvec zmatrix(const numvec& p1, double dist);
    
    numvec zmatrix(const numvec& p1, double dist, const numvec& p2, double angl);
    
    numvec zmatrix(const numvec& p1, double dist, const numvec& p2, double angl, const numvec& p3, double tors);
 
    // transformation based on a 4x4 transform matrix, could be either way
    void transform( double* v, const double* mat, bool inverse=false );    

    void translate( double* v, const numvec& sft );

} // namespace mort

#endif

