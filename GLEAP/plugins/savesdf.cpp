#include <fstream>
#include <object.hpp>
#include <format.hpp>
#include "savesdf.hpp"

namespace amber
{
    using namespace mort;

    savesdf_command::savesdf_command( )
	:command_i( "savesdf" )
    {
    }

    savesdf_command::savesdf_command( const string& name, const string& file )
        :m_name( name ),m_file( file )
    {
    }

    savesdf_command::~savesdf_command()
    {
    }

    bool savesdf_command::exec( )
    {
        std::ofstream os( m_file.c_str() );

        if( os == NULL )
        {
            throw std::runtime_error( "Error: can not open file " + m_file + " for written" );
        }

        molecule_ptr pmol = content().get_mol( m_name );
        assert( pmol != NULL );

        write_sdf( os, *pmol );

        return true;        
    }

    void savesdf_command::undo( )
    {
    }

    const char* savesdf_command::info( ) const
    {
        static const char usage[] = "Usage: savesdf molname filename";
        return usage;
    }
   
    shared_ptr< command_i > savesdf_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 3 )
        {
            throw std::logic_error( string( "Error: wrong arguments\n" ) + info() );
        }
    
        return shared_ptr< command_i >( new savesdf_command( args[1], args[2] ) );
    }

}

static amber::savesdf_command g_savesdf_command;


