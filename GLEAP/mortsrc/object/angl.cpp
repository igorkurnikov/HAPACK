#include <stdexcept>
#include "atom.hpp"
#include "bond.hpp"
#include "angl.hpp"
#include "mole.hpp"

namespace mort
{
    using std::runtime_error;

    angl_t::angl_t( const morf_t& rhs )
        : morf_t( rhs )
    {
    }

    angl_t::angl_t( const molecule_t& m, int id )
        : morf_t( m, ANGL, id )
    {
    }

    angl_t angl_t::create( atom_t& a1, atom_t& a2, atom_t& a3 )
    {
        molecule_t& m = a1.getmol();

        angl_t an = m.create( ANGL );
        atom_t af = std::min( a1, a3 );
        atom_t ab = std::max( a1, a3 );
        an.connect( a1 );
        an.connect( a2 );
        an.connect( a3 );
        a2.connect( an );
	return an;
    }

    angl_t angl_t::frcget( atom_t& a1, atom_t& a2, atom_t& a3)
    {
        angl_t an(a1.getmol(),-1);
        return get(a1,a2,a3,an) ? an : create(a1,a2,a3);
    }

    angl_t angl_t::get( const atom_t& a1, const atom_t& a2, const atom_t& a3 )
    {
        angl_t an(a1.getmol(),-1);
        if( get(a1,a2,a3,an) )
        {
            return an;
        }

        throw runtime_error( "Error: angl does not exist among " 
                             + a1.name() + "-" + a2.name() + "-" + a3.name() );
    }

    bool angl_t::get( const atom_t& a1, const atom_t& a2, const atom_t& a3, angl_t& an )
    {
        atom_t af = std::min( a1, a3 );
        atom_t ab = std::max( a1, a3 );

        angliter_t ai = a2.angl_begin();
        angliter_t ae = a2.angl_end();
        for( ; ai != ae; ++ai )
        {
            atom_range as = ai->atoms();
            if( as[0]==a1 && as[2]==a3 )
            {
                an = *ai;
                return true;
            }

            if( as[0]==a3 && as[2]==a1 )
            {
                an = *ai;
                return true;
            }
        }

        return false;
    }

    bool angl_t::has( const atom_t& a1, const atom_t& a2, const atom_t& a3 )
    {
        angl_t an(a1.getmol(),-1);
        return get(a1,a2,a3,an);
    }

    bool angl_t::has( const atom_t& a1, const atom_t& a3 )
    {
        atomiter_t ai = a1.atom_begin();
	for( ; ai != a1.atom_end(); ++ai )
	{
	    if( ai->is_connected_to(a3) )
	    {
	        return true;
	    }
	}

	return false;
    }


} // namespace mort

