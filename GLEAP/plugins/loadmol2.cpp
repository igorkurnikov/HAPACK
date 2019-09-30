#include <fstream>
#include <object.hpp>
#include <format.hpp>
#include "loadmol2.hpp"

namespace amber
{
    using namespace mort;

    loadmol2_command::loadmol2_command( )
	:command_i( "loadmol2" )
    {
    }

    loadmol2_command::loadmol2_command( const string& name, const string& file )
        :m_name( name ),m_file( file )
    {
    }

    loadmol2_command::~loadmol2_command()
    {
    }

    bool loadmol2_command::exec( )
    {
        std::ifstream is( m_file.c_str() );

        if( is == NULL )
        {
            throw std::runtime_error( "Error: can not open file " + m_file + " for read" );
        }

        shared_ptr< molecule_t > pmol( new molecule_t() );

        read_mol2( is, *pmol );

        pmol->set_s(NAME, m_name);
    
        content().set( m_name, pmol );

        return true;        
    }

    void loadmol2_command::undo( )
    {
        content().remove( m_name.c_str() );    
    }

    const char* loadmol2_command::info( ) const
    {
        static const char usage[] = "Usage: molname = loadmol2 filename || loadmol2 molname filename";
        return usage;
    }
   
    shared_ptr< command_i > loadmol2_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 3 )
        {
            throw std::runtime_error( string( "Error: wrong arguments\n" ) + info() );
        }
    
        return shared_ptr< command_i >( new loadmol2_command( args[1], args[2] ) );
    }

}

static amber::loadmol2_command g_loadmol2_command;


