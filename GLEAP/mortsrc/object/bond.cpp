#include <iostream>
#include <stdexcept>
#include "mole.hpp"
#include "atom.hpp"
#include "bond.hpp"
#include "resd.hpp"

namespace mort
{
    using std::runtime_error;

    bond_t::bond_t( const morf_t& rhs )
        : morf_t( rhs )
    {
        int cid = rhs.cmpid();
	if( cid!=BOND ) 
        {
            throw std::runtime_error( "Error: try to initialize bond with a " + unhash(cid) );
        }

    }

    bond_t::bond_t( const molecule_t& m, int id )
        : morf_t( m, BOND, id )
    {
    }

    bond_t bond_t::create( atom_t& a1, atom_t& a2 )
    {
        molecule_t& m = a1.getmol();
        if( &m != &(a2.getmol()) )
        {
            std::cout << "atom1: " << a1.get_s(NAME) << " of " << a1.getmol().get_s(NAME) << std::endl;
            std::cout << "atom2: " << a2.get_s(NAME) << " of " << a2.getmol().get_s(NAME) << std::endl;
            throw std::runtime_error( "Error: try to make bond between atom from two molecules" );
        }

        bond_t b = m.create(BOND);
	
	a1.connect( a2 );
        a2.connect( a1 );

        a1.connect( b );
	a2.connect( b );

	if( a1.relid() < a2.relid() )
        {
            b.connect( a1 );
	    b.connect( a2 );
        }
        else
        {
            b.connect( a2 );
            b.connect( a1 );
        }

        assert( a1.is_connected_to(a2) );
        assert( a2.is_connected_to(a1) );

        if( m.nresd()>0 )
        {
            // set head and tail parameter for residue
            resd_t r1 = a1.resd();
            resd_t r2 = a2.resd();

            if( r1==r2 )
            {
                r1.connect( b );
                b.connect( r1 );
            }
            else
            {
                r1.connect( b );
                r2.connect( b );
                b.connect( r1 );
                b.connect( r2 );
                r1.connect( r2 );
                r2.connect( r1 );
            }
        }

	return b;
    }

    bond_t bond_t::frcget( atom_t& a1, atom_t& a2 )
    {
        bond_t b( a1.getmol(), -1);
        return get(a1,a2,b) ? b : create(a1,a2);
    }

    bond_t bond_t::get( const atom_t& a1, const atom_t& a2 )
    {
        bond_t b( a1.getmol(), -1);
        if( get(a1,a2,b) )
        {
            return b;
        }

        throw runtime_error( "Error: cannot get bond between " + a1.name() + " and " + a2.name() );
    }

    bool bond_t::get( const atom_t& a1, const atom_t& a2, bond_t& b )
    {
        atomiter_t ai = a1.atom_begin();
        atomiter_t ae = a1.atom_end();
        bonditer_t bi = a1.bond_begin();
        for( ; ai!=ae; ++ai, ++bi )
        {
            if( *ai==a2 )
            {
                b = *bi;
                return true;
            }
        }
        return false;
    }

    bool bond_t::has( const atom_t& a1, const atom_t& a2 )
    {
        bond_t b(a1.getmol(),-1);
        return get(a1,a2,b);
    }

    bool bond_t::has( const resd_t& r1, const resd_t& r2 )
    {
        int rid2 = r2.get_i(ID);

        atomiter_t ai = r1.atom_begin();
        for( ; ai != r1.atom_end(); ++ai )
        {
            atomiter_t aj = ai->atom_begin();
            for( ; aj != ai->atom_end(); ++aj )
            {
                if( aj->get_i(RESID)==rid2 )
                {
                    return true;
                }
            }
        }

        return false;
    }

} // namespace mort

