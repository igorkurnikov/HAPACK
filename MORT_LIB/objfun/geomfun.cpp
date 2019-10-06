#include <common.hpp>
#include <object.hpp>
#include "parmfun.hpp"
#include "geomfun.hpp"

namespace mort
{
    double dis2( const atom_t& a1, const atom_t& a2 )
    {
        return dis2( a1.pos(), a2.pos() );
    }

    double dist( const atom_t& a1, const atom_t& a2 )
    {
        return dist( a1.pos(), a2.pos() );
    }

    double angl( const atom_t& a1, const atom_t& a2, const atom_t& a3 )
    {
        return angl( a1.pos(), a2.pos(), a3.pos() );
    }

    double tors( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4 )
    {
        return tors( a1.pos(), a2.pos(), a3.pos(), a4.pos() );
    }

    numvec center( const molecule_t& m )
    {
        return center( get_vvec(m, ATOM, POSITION) );
    }

    numvec center( const resd_t& r )
    {
        const vector<double>& pos = get_vvec(r.getmol(), ATOM, POSITION);
        return center( pos, r.related_atom_ids() );
    }

    void translate( molecule_t& m, const numvec& sft )
    {
        translate( get_vvec(m, ATOM, POSITION), sft );
    }

    void translate( morf_t& mo, const numvec& sft )
    {
        vector<double>& crd = get_vvec(mo.getmol(), ATOM, POSITION);
        if( mo.cmpid()==ATOM )
        {
            int id = mo.absid();
            translate( &crd[3*id] , sft );
        }
        else if( mo.cmpid()==RESD )
        {
            translate( crd, sft, mo.related_atom_ids() );
        }
        else
        {
            throw std::runtime_error( "Error: cannot translate a " + unhash(mo.cmpid()) );
        }
    }

    void transform( molecule_t& m, const double* mat )
    {
        transform( get_vvec(m, ATOM, POSITION), mat );
    }

    void transform( morf_t& mo, const double* mat )
    {
        vector<double>& crd = get_vvec(mo.getmol(), ATOM, POSITION);
        if( mo.cmpid()==ATOM )
        {
            int id = mo.absid();
            transform( &crd[3*id] , mat );
        }
        else if( mo.cmpid()==RESD )
        {
            transform( crd, mat, mo.related_atom_ids() );
        }
        else
        {
            throw std::runtime_error( "Error: cannot transform a " + unhash(mo.cmpid()) );
        }
    }


    void transform( database_t& mdb, const double* mat )
    {
        database_t::iterator i = mdb.begin();
        for( ; i != mdb.end(); ++i )
        {
            molecule_ptr pm = dynamic_pointer_cast< molecule_t >(i->second);
            if( i->first[0] !='_' && pm !=NULL )
            {
                transform( *pm, mat );
            }
        }
    }

    void rotate( molecule_t& m, const matrix& mat )
    {
        rotate( get_vvec(m, ATOM, POSITION), mat );
    }

    void rotate( morf_t& mo, const matrix& mat )
    {
        vector<double>& crd = get_vvec(mo.getmol(), ATOM, POSITION);
        if( mo.cmpid()==ATOM )
        {
            int id = mo.absid();
            rotate( &crd[3*id] , mat );
        }
        else if( mo.cmpid()==RESD )
        {
            rotate( crd, mat, mo.related_atom_ids() );
        }
        else
        {
            throw std::runtime_error( "Error: cannot rotate a " + unhash(mo.cmpid()) );
        }
    }

} // namespace mort

