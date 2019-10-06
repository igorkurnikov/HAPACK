#include <iostream>
#include <stdexcept>
#include "atom.hpp"
#include "bond.hpp"
#include "angl.hpp"
#include "impr.hpp"
#include "mole.hpp"
#include "atmv.hpp"

namespace mort
{
    using std::runtime_error;

    atmvec sort( const atom_t& a1, const atom_t& a2, const atom_t& a3 )
    {
        if( a1 > a2 )
        {
            return sort( a2, a1, a3 );
        }

        assert( a1 <= a2 );

        atmvec as(3, a1);
        if( a1 < a3 )
        {
            as[0] = a1;
            as[1] = std::min(a2, a3);
            as[2] = std::max(a2, a3);
        }
        else
        {
            as[0] = a3;
            as[1] = a1;
            as[2] = a2;
        }

        return as;
    }


    static const hashid_t IMPR = OOPS;

    impr_t::impr_t( const morf_t& rhs )
        : morf_t( rhs )
    {
    }

    impr_t::impr_t( const molecule_t& m, int id )
        : morf_t( m, IMPR, id )
    {
    }

    impr_t impr_t::create( atom_t& a1, atom_t& a2, atom_t& a3, atom_t& a4, int p )
    {
        mole_t& m = a1.getmol();

        impr_t im = m.create( IMPR );
        im.set_i(PERIOD, p);


        im.connect( a1 );
        im.connect( a2 );
        im.connect( a3 );
        im.connect( a4 );

        bond_t b( a1.getmol(), -1 );
        if( !bond_t::get(a3, a4, b) )
        {
            throw runtime_error( "Error: while creating improper torsion, central bond: "+a3.name() + a4.name() 
                                 + " does not exist." );
        }

        b.connect( im );

        // the following connection is needed for amoeba force field.
        angl_t an(a1.getmol(), -1);
        if( angl_t::get(a1,a3,a2,an) )
        {
            an.connect( im );
        }

        if( angl_t::get(a1,a3,a4,an) )
        {
            an.connect( im );
        }

        if( angl_t::get(a2,a3,a4,an) )
        {
            an.connect( im );
        }

        a3.connect( im );
        return im;
    }

    impr_t impr_t::frcget( atom_t& a1, atom_t& a2, atom_t& a3, atom_t& a4, int p )
    {
        impr_t im(a1.getmol(),-1);
        return get(a1,a2,a3,a4,im) ? im : create(a1,a2,a3,a4,p);
    }

    impr_t impr_t::get( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4 )
    {
        impr_t im(a1.getmol(),-1);
        if( get(a1,a2,a3,a4,im) )
        {
            return im;
        }

        throw runtime_error( "Error: cannot get improper torsion among: " + a1.name() + "-" 
                             + a2.name() + "-" + a3.name() + "-" + a4.name() );
    }

    bool impr_t::get( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4, impr_t& im )
    {
        bond_t b = bond_t::get( a3, a4 );
        impriter_t i = b.impr_begin();
        impriter_t e = b.impr_end();
        for( ; i != e; ++i )
        {
            atom_range ar = i->atoms();
            if( a1==ar[0] && a2==ar[1] )
            {
                im = *i;
                return true;
            }
        }

        return false;
    }

    bool impr_t::has( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4 )
    {
        impr_t im(a1.getmol(),-1);
        return get(a1,a2,a3,a4,im);
    }

} // namespace mort

