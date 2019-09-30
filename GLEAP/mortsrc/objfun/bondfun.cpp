#include <cassert>
#include <sstream>
#include <stdexcept>
#include <common.hpp>
#include <object.hpp>
#include "bondfun.hpp"
#include "atomfun.hpp"
#include "geomfun.hpp"

#ifndef __GNUC__
#include <ctype.h>
#endif

namespace mort
{
    using std::runtime_error;

    bool has_prev( const morf_t& tors )
    {
        morf_t a1 = atom_1st( tors );
        morf_t a4 = atom_4th( tors );
        int id = tors.get_i(ID);

        diheiter_t ti = a1.dihe_begin();
        for( ; ti != a1.dihe_end(); ++ti )
        {
            if( ti->get_i(ID) >= id )
            {
                continue;
            }

            if( atom_1st(*ti)==a1 && atom_4th(*ti)==a4 )
            {
                return true;
            }
            
            if( atom_1st(*ti)==a4 && atom_4th(*ti)==a1 )
            {
                return true;
            }
        }
        
        return false;
    }

    morf_t create_tor2(atomvec_t& atms )
    {
        assert( atms.size() == 5 );

        angl_t angl = angl_t::get(atms[1], atms[2], atms[3]);
	morf_t tor2 = atms[0].getmol().create(TOR2);

	angl.connect( tor2 );
	tor2.connect( angl );

        atms[0].connect( tor2 );
	atms[4].connect( tor2 );

	tor2.connect( atms[0] );
	tor2.connect( atms[1] );
	tor2.connect( atms[2] );
	tor2.connect( atms[3] );
	tor2.connect( atms[4] );

	return tor2;
    }

    morf_t create_ptor(atomvec_t& atms, int period)
    {
        assert( atms.size() == 6 );
        morf_t ptor = atms[0].getmol().create(PTOR);
        ptor.set_i(PERIOD, period);
        ptor.connect( atms[0] );
        ptor.connect( atms[1] );
        ptor.connect( atms[2] );
        ptor.connect( atms[3] );
        ptor.connect( atms[4] );
        ptor.connect( atms[5] );
        return ptor;
    }

    bool bond_bydis( atom_t& a1, atom_t& a2, double cutoff )
    {
        if( a1.is_connected_to(a2) )
            return true;

        double d2 = dis2( a1, a2 );
        if( d2 < BUMPCUT2 )
        {
            std::cout << "Warning: atom " << a1.name() << " ";
            std::cout << " is too close to atom " << a2.name();
            std::cout << ", distance is only " << sqrt(d2) << " angstrom." << std::endl;
            std::cout << "a bond will not be created between them." << std::endl;
            return false;
        }

        double cut2 = cutoff*cutoff;
        if( a1.get_i(ELEMENT)==HYDROGEN || a2.get_i(ELEMENT)==HYDROGEN )
        {
            cut2 *= 0.49;
        }

        if( d2 < cut2 )
        {
            bond_t b = bond_t::create(a1, a2);
            b.set_i(ORDER, 1);
            return true;
        }

        return false;
    }

    void bond_bydis(atomiter_t begin, atomiter_t end, double cutoff)
    {
        for(atomiter_t i = begin; i != end; ++i )
	{
	    for(atomiter_t j=i+1; j != end; ++j )
	    {
	        if( i->is_connected_to(*j) )
		{
		    continue;
		}

                bond_bydis( *i, *j, cutoff );
       	    }
        }
    }


    void bond_bydis( molecule_t& m, double cutoff )
    {
        bond_bydis( m.atom_begin(), m.atom_end(), cutoff );
    }

    void bond_bydis(resd_t& r, double cutoff)
    {
        bond_bydis( r.atom_begin(), r.atom_end(), cutoff );
    }


} // namespace mort


