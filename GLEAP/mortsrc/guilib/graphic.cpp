#include <stdexcept>
#include <GL/gl.h>
#include <object.hpp>
#include "ribbon.hpp"
#include "graphic.hpp"

namespace mort
{
    using std::logic_error;

    molecule_graphic::molecule_graphic( const molecule_t& mol, const hashid_t& style )
 	:m_molecule( mol ), m_name( "mol_" + mol.get_s(NAME) )
    {
        if( style == LINE )
        {
            m_style = shared_ptr< style_t >( new line_style( 1.5 ) );
        }
        else if( style == BALLSTICK )
        {
            m_style = shared_ptr< style_t >( new ballstick_style( 0.2, 0.1, ELEMENT_COLOR ) );
        }
        else
        {
            throw logic_error( "Error: unknown molecule render style." );
        }
    }

    molecule_graphic::~molecule_graphic()
    {
    }

    string molecule_graphic::name() const
    {
        return m_name;
    }

    void molecule_graphic::paint()
    {
        iterator_t atom = m_molecule.atom_begin();

        for( ; atom != m_molecule.atom_end(); ++atom )
        {
	    int visible;
            if( !atom->get_i(VISIBLE,visible) || visible > 0  )
            {
                m_style->render_atom( *atom );
            }
        }

        iterator_t bond = m_molecule.bond_begin();

        for( ; bond != m_molecule.bond_end(); ++bond )
        {
	    int visible;

            if( atom_1st(*bond).get_i(VISIBLE, visible) && visible == 0 )
	    {
	        continue;
	    }

            if( atom_2nd(*bond).get_i(VISIBLE, visible) && visible == 0 )
	    {
	        continue;
	    }

            m_style->render_bond( *bond );
        }
    }

    ribbon_graphic::ribbon_graphic( molecule_t& mol )
        : m_name( "ribbon_" + mol.get_s(NAME) )
    {
        // create_ribbon( mol );
        m_pmol = &mol;
    }

    ribbon_graphic::~ribbon_graphic()
    {
    }

    string ribbon_graphic::name( ) const
    {
         return m_name;
    }

    void ribbon_graphic::paint()
    {
        render_ribbon( *m_pmol );
    }

} // namespace mort


