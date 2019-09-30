#include <fstream>
#include <object.hpp>
#include <format.hpp>
#include <guilib.hpp>
#include "savemol2.hpp"

namespace amber
{
    using namespace mort;

    savemol2_command::savemol2_command( )
	:command_i( "savemol2" )
    {
    }

    savemol2_command::savemol2_command( const string& name, const string& file )
        :m_name( name ),m_file( file )
    {
    }

    savemol2_command::~savemol2_command()
    {
    }

    bool savemol2_command::exec( )
    {
        std::ofstream os( m_file.c_str() );

        if( os == NULL )
        {
            throw std::runtime_error( "Error: can not open file " + m_file + " for written" );
        }

        molecule_ptr pmol = content().get_mol( m_name );
        if( pmol == NULL )
	{
	    throw std::runtime_error( "Error: molecule " + m_name + " does not exist." );
	}

        write_mol2( os, *pmol );

        return true;        
    }

    void savemol2_command::undo( )
    {
    }

    const char* savemol2_command::info( ) const
    {
        static const char usage[] = "Usage: savemol2 molname filename";
        return usage;
    }
   
    shared_ptr< command_i > savemol2_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 3 )
        {
            throw std::runtime_error( string( "Error: wrong arguments\n" ) + info() );
        }
    
        return shared_ptr< command_i >( new savemol2_command( args[1], args[2] ) );
    }

}

static amber::savemol2_command g_savemol2_command;


