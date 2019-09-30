#include "impose.hpp"
#include "zmatrix.hpp"

namespace amber
{
    void apply( resd_t& resd, const vector< string >& inter )
    {
        assert( resd.cmpid() == RESD );

        if( inter.size() == 3 )
        {
            atomvec_t atoms;
            atoms.push_back( atom_t::get(resd, inter[0]) );
            atoms.push_back( atom_t::get(resd, inter[1]) );

            impose_dist( atoms, atof( inter[2].c_str() ) );
        }
        else if( inter.size() == 5 )
        {
            atomvec_t atoms;
            atoms.push_back( atom_t::get(resd, inter[0]) );
            atoms.push_back( atom_t::get(resd, inter[1]) );
            impose_dist( atoms, atof( inter[3].c_str() ) );

            atoms.push_back( atom_t::get(resd, inter[2]) );
            impose_angl( atoms, atof( inter[4].c_str() ) );
        }
        else if( inter.size() == 7 )
        {
            atomvec_t atoms;
            atoms.push_back( atom_t::get(resd, inter[0]) );

            atoms.push_back( atom_t::get(resd, inter[1]) );
            impose_dist( atoms, atof( inter[4].c_str() ) );

            atoms.push_back( atom_t::get(resd, inter[2]) );
            impose_angl( atoms, atof( inter[5].c_str() ) );

            atoms.push_back( atom_t::get(resd, inter[3]) );
            impose_tors( atoms, atof( inter[6].c_str() ), false );
        }
        else
        {
            throw std::runtime_error( "Error: wrong inters" );
        }
    }
 
    zmatrix_command::zmatrix_command( )
        : command_i( "zmatrix" )
    {
    }

    zmatrix_command::zmatrix_command( const string& object, const vector< vector< string > >& inters )
        : m_object( object ), m_inters( inters )
    {
    }

    zmatrix_command::~zmatrix_command( )
    {
    }

    bool zmatrix_command::exec( )
    {
        int ndot = std::count( m_object.begin(), m_object.end(), '.' );

        if( ndot > 1 )
        {
            throw std::runtime_error( " zmatrix target should be molecule or residue" );
        }

        if( ndot == 0 )
        {
            molecule_ptr pmol = content().get_mol( m_object );
            assert( pmol != NULL );

            if( pmol->nresd() != 1 )
            {
                throw std::runtime_error( "Error: zmatrix should only be used on monomer" );
            }

            resditer_t resd = pmol->resd_begin();
            vector< vector<string> >::iterator i = m_inters.begin();
            vector< vector<string> >::iterator e = m_inters.end();
            for( ; i != e; ++i )
            {
                apply( *resd, *i );
            }
        }
        else
        {
            assert( ndot == 1 );
            resd_t resd = content().get_resd( m_object );
            vector< vector<string> >::iterator i = m_inters.begin();
            vector< vector<string> >::iterator e = m_inters.end();
            for( ; i != e; ++i )
            {
                apply(resd, *i );
            }
        }

	return true;
    }

    void zmatrix_command::undo( )
    {
        throw std::runtime_error( "Sorry: not implmented yet" );
    }

    const char* zmatrix_command::info( ) const
    {
        return "  usage: zmatrix object matrix ";
    }

    shared_ptr< command_i > zmatrix_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 3 )
        {
            throw std::runtime_error( "Error: wrong number of arguments" );
        }

        return shared_ptr< command_i >( new zmatrix_command( args[1], interpret2d( args[2] ) ) );
    }

}  // namespace amber

amber::zmatrix_command g_zmatrix_command;

         
