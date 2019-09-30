#include <fstream>
#include <object.hpp>
#include <format.hpp>

#include "savepdb.hpp"

namespace amber
{
    using namespace mort;

    savepdb_command::savepdb_command( )
	:command_i( "savepdb" )
    {
    }

    savepdb_command::savepdb_command( const string& name, const string& file )
        :m_name( name ),m_file( file )
    {
    }

    savepdb_command::~savepdb_command()
    {
    }

    bool savepdb_command::exec( )
    {
        std::ofstream os( m_file.c_str() );

        if( os == NULL )
        {
            throw std::runtime_error( "Error: can not open file " + m_file + " for written" );
        }

        std::cout << "saving " << m_name << " to " << m_file << std::endl;

        molecule_ptr pmol = content().get_mol( m_name );
        assert( pmol != NULL );
        
        write_pdb( os, *pmol );

        return true;        
    }

    void savepdb_command::undo( )
    {
    }

    const char* savepdb_command::info( ) const
    {
        static const char usage[] = "Usage: savepdb molname filename";
        return usage;
    }
   
    shared_ptr< command_i > savepdb_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 3 )
        {
            throw std::logic_error( string( "Error: wrong arguments\n" ) + info() );
        }
    
        return shared_ptr< command_i >( new savepdb_command( args[1], args[2] ) );
    }

}

static amber::savepdb_command g_savepdb_command;


