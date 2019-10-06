#ifndef MORTSRC_COMMON_CRDALGO_HPP
#define MORTSRC_COMMON_CRDALGO_HPP

#include <vector>
#include "numvec.hpp"

namespace mort
{
    using std::vector;

    class matrix;
    // crdalgo (coordinate algorithm) contains a set of algorithms operates on 3d coordinates
    // the coordinates were stored in a 1d vector of type vector<double>

    // return the geometrical center of the coordinates.
    numvec center( const vector<double>& crd );

    numvec center( const vector<double>& crd, const vector<int>& ids );

    // return the min and max of coordinates.
    //     return value: a 6 elements numvec, first 3 elements are min, last 3 are max.
    numvec minmax( const vector<double>& crd );

    // return the min and max of coordinates with consideration of vdw radius.
    //     return value: a 6 elements numvec, first 3 elements are min, last 3 are max.
    numvec minmax( const vector<double>& crd, const vector<double>& vdwr );

    // return the region center and extent
    //    center = (min+max)/2.0
    //    extent = (max-min) 
    numvec region( const vector<double>& crd );

    // return the region with vdwr
    numvec region( const vector<double>& crd, const vector<double>& vdwr );

    // return the inter moment.
    matrix moment( const vector<double>& crd );

    // return the buffer size of octrahedron. sometime user input is too small that part
    // of the solute is exposed to vacuum.
    double octbufsize( const vector<double>& crd, double uinput );

    void rotate( vector<double>& crd, const matrix& mat );

    void rotate( vector<double>& crd, const matrix& mat, const vector<int>& ids );

    // transform 
    void transform( vector<double>& crd, const double* mat );

    void transform( vector<double>& crd, const double* mat, const vector<int>& ids );

    // translation
    void translate( vector<double>& crd, const numvec& sft );

    void translate( vector<double>& crd, const numvec& sft, const vector<int>& ids );

    // translate so that center is on origin
    void centralize( vector<double>& crd );

    // translate so that region center is on origin.
    // return the extent of the region.
    numvec regionlize( vector<double>& crd );

    // translate so that region center is on origin (considering vdwr).
    // return the extent of the region.
    numvec regionlize( vector<double>& crd, const vector<double>& vdwr );

    // ewald rotation
    void rotatewald( vector<double>& crd );

    // longest rotation: rotate 
    void rotatelong( vector<double>& crd );

    

} // namespace mort

#endif


