#include <common.hpp>
#include "status.hpp"
#include "atomexpr.hpp"

namespace mort
{
    
    namespace daylight
    {

        bondinfo_t::bondinfo_t()
        {
            reset();
        }
            
        void bondinfo_t::reset()
        {
            arom = false;
            order = 1;
            stereo = 1;
        }

        void bondinfo_t::apply(morf_t& bond)
        {
            bond.set_i(AROM,   arom);
            bond.set_i(ORDER,  order);
            bond.set_i(STEREO, stereo);
            reset();
        }

        status_t::status_t(const hashid_t& lan, molecule_t& mol )
            :curt_atom(mol, -1), 
             prev_atom(mol, -1),
             ring_atoms(10, atom_t(mol, -1) ) 
        {
            pmol      = &mol;
            bracket   = CLOSE;
            language  = lan;
        }

        void status_t::new_atom( )
        {
            prev_atom = curt_atom;
            curt_atom = atom_t::create( *pmol, "");

	    assert( curt_atom.absid() >= 0 );

            if( language == SMARTS )
            {
                curt_atom.set_a(PATTERN, shared_ptr< atomexpr_t >( new atomexpr_t() ) );
            }
        }

        void status_t::end_atom( )
        {
            if( prev_atom.absid() != -1 && bond_info.order != 0 )
            {
	        bond_t bond = bond_t::create(prev_atom, curt_atom);
                bond_info.apply(bond);
            }
        }

    } // namespace daylight
    
} // namespace mort

