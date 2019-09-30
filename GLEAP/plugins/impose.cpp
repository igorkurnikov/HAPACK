#include "impose.hpp"

namespace amber
{
    atomvec_t select( const molecule_t& mol, const vector< string >& resds )
    {
        atomvec_t result;

        if( resds.size() == 2 )
        {
            int bgn = atoi( resds[0].c_str() );
            int end = atoi( resds[1].c_str() );
           
            if( end < bgn )
            {
                std::swap( bgn, end );
            }

            for( int i = bgn; i <= end; ++i )
            {
                result.push_atom( morf_t(mol, RESD, i-1) );
            }
        }
        else
        {
            for( int i = 0; i < (int)resds.size(); ++i )
            {
                int idx = atoi( resds[i].c_str() );
                result.push_atom( morf_t( mol, RESD, idx - 1 ) );
            }
        }
    
        return result;
    }

    impose_command::impose_command( )
        : command_i( "impose" )
    {
    }

    impose_command::impose_command( const string& unit, const vector< string >& resds, const vector< vector< string > >& inters )
        : m_unit( unit ), m_resds( resds ), m_inters( inters )
    {
    }

    impose_command::~impose_command( )
    {
    }

    bool impose_command::exec( )
    {
        molecule_ptr pmol = content().get_mol( m_unit );
        assert( pmol != NULL );
        atomvec_t atoms = select( *pmol, m_resds );

        for( int i=0; i < (int)m_inters.size(); ++i )
        {
            if( m_inters[i].size() == 3 )
            {
                impose_dist( atoms, m_inters[i] );
            }
            else if( m_inters[i].size() == 4 )
            {
                impose_angl( atoms, m_inters[i] );
            }
            else if( m_inters[i].size() == 5 )
            {
                impose_tors( atoms, m_inters[i] );
            }
            else
            {
                throw std::runtime_error( "Error: can't understand internal constraint " + m_inters[i][0] );
            }
        }

        return true;
/*
        for( int i=0; i < constraints.size(); ++i )
        {
            if( inters[i].size() == 3 )
            {
                exam_dist( atoms, inters[i] );
            }
            else if( inters[i].size() == 4 )
            {
                exam_angl( atoms, inters[i] );
            }
            else
            {
                assert( inters[i].size() == 5 );
                exam_tors( atoms, inters[i] );
            }
        }
*/
    }

    void impose_command::undo( )
    {
        throw std::runtime_error( "Sorry, not implemented yet" );
    }

    const char* impose_command::info( ) const 
    {
        return "  usage: impose unit seqlist internals ";
    }

    shared_ptr< command_i > impose_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 4 )
        {
            throw std::runtime_error( "Error: wrong number of arguments " );
        }

        return shared_ptr< command_i >( new impose_command( args[1], interpret1d( args[2] ), interpret2d( args[3] ) ) );
    }

} // namespace amber

amber::impose_command g_impose_command;

