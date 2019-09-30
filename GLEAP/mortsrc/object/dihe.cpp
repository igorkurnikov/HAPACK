#include <iostream>
#include <stdexcept>
#include "atom.hpp"
#include "bond.hpp"
#include "dihe.hpp"
#include "mole.hpp"

namespace mort
{
    using std::runtime_error;

    enum alignment_e {LEFT, RIGHT};

    typedef molecule_t mole_t;

    int side( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4 )
    {
        if( a2 < a3 ) return LEFT;
        if( a2 > a3 ) return RIGHT;
        return (a1 < a4) ? LEFT : RIGHT;
    }

    dihvec getdih( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4, int p )
    {
        dihvec ds;

        bond_t bc(a1.getmol(), -1);
        if( !bond_t::get(a2, a3, bc) )
        {
            return ds;
        }


        diheiter_t i = bc.dihe_begin();
        diheiter_t e = bc.dihe_end();
        for( ; i != e; ++i )
        {
            atom_range as = i->atoms();
            if( (a1==as[0]&&a2==as[1]&&a3==as[2]&&a4==as[3]) || 
                (a4==as[0]&&a3==as[1]&&a2==as[2]&&a1==as[3]) )
            {
                if( p==-1 )
                {
                    ds.push_back(*i);
                }
                else if( i->period()==p )
                {
                    ds.push_back(*i);
                    return ds;
                }
            }
        }

        return ds;
    }

    static const hashid_t DIHE =TORS;

    dihe_t::dihe_t( const morf_t& rhs )
        : morf_t( rhs )
    {
    }

    dihe_t::dihe_t( const mole_t& m, int id )
        : morf_t( m, TORS, id )
    {
    }

    int dihe_t::period() const
    {
        return get_i(PERIOD);
    }

    dihe_t dihe_t::create( atom_t& a1, atom_t& a2, atom_t& a3, atom_t& a4, int p )
    {
        mole_t& m = a3.getmol();
        bond_t bc( m, -1 );
        if( !bond_t::get(a2, a3, bc) )
        {
            throw runtime_error( "Error: while creating dihedral " + a1.name() + "-"
                                 + a2.name() + "-" + a3.name() + "-" + a4.name() +
                                 " central bond does not exist." );
        }

        dihe_t d = m.create( DIHE );
        d.set_i(PERIOD, p);

        bc.connect( d );
        a1.connect( d );
        a4.connect( d );
/*
        if( side(a1,a2,a3,a4)==LEFT )
        {
            d.connect( a1 );
            d.connect( a2 );
            d.connect( a3 );
            d.connect( a4 );
        }
        else
        {
            d.connect( a4 );
            d.connect( a3 );
            d.connect( a2 );
            d.connect( a1 );
        }
*/
        d.connect( a1 );
        d.connect( a2 );
        d.connect( a3 );
        d.connect( a4 );
      
        return d;
    }

    dihe_t dihe_t::frcget( atom_t& a1, atom_t& a2, atom_t& a3, atom_t& a4, int p )
    {
        dihe_t d(a1.getmol(),-1);
        return get(a1,a2,a3,a4,p,d) ? d : create(a1,a2,a3,a4,p);
    }

    dihe_t dihe_t::get( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4, int p )
    {
        dihe_t d(a1.getmol(),-1);
        if( get(a1,a2,a3,a4,p,d) )
        {
            return d;
        }

        throw runtime_error( "Error: cannot get dihedral torsion among: " + a1.name() + "-" 
                             + a2.name() + "-" + a3.name() + "-" + a4.name() );
    }

    dihvec dihe_t::get( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4 )
    {
        return getdih(a1,a2,a3,a4,-1);
    }

    bool dihe_t::get( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4, int p, dihe_t& d )
    {
        dihvec ds = getdih(a1,a2,a3,a4,p);
        if( ds.size()==0 )
        {
            return false;
        }

        d = ds[0];
        return true;
    }

    bool dihe_t::has( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4, int p )
    {
        dihe_t d( a1.getmol(), -1);
        return get(a1,a2,a3,a4,p,d);
    }

} // namespace mort

