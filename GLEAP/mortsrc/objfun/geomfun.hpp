#ifndef MORTSRC_OBJECT_GEOMFUN_HPP
#define MORTSRC_OBJECT_GEOMFUN_HPP

#include <common.hpp>

namespace mort
{
    class atom_t;

    class resd_t;

    class morf_t;

    class molecule_t;
 
    class database_t;

    double dis2( const atom_t& a1, const atom_t& a2 );

    double dist( const atom_t& a1, const atom_t& a2 );

    double angl( const atom_t& a1, const atom_t& a2, const atom_t& a3 );

    double tors( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4 );

    numvec center( const molecule_t& m );

    numvec center( const resd_t& mo );

    void translate( morf_t& mo, const numvec& sft );

    void translate( molecule_t& m, const numvec& sft );

    void transform( morf_t& mo, const double* mat );

    void transform( molecule_t& mol, const double* mat );

    void transform( database_t& mdb, const double* mat, bool invert=false );

    void rotate( molecule_t& m, const matrix& rot );

    void rotate( morf_t& mo, const matrix& rot );

}

#endif

