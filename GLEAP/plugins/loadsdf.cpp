#include <fstream>
#include <object.hpp>
#include <format.hpp>

#include "loadsdf.hpp"

namespace amber
{
    using namespace mort;

    loadsdf_command::loadsdf_command( )
	:command_i( "loadsdf" )
    {
    }

    loadsdf_command::loadsdf_command( const string& name, const string& file )
        :m_name( name ),m_file( file )
    {
    }

    loadsdf_command::~loadsdf_command()
    {
    }

    bool loadsdf_command::exec( )
    {
        std::ifstream is( m_file.c_str() );

        if( is == NULL )
        {
            throw std::runtime_error( "Error: can not open file " + m_file + " for read" );
        }

        shared_ptr< molecule_t > pmol( new molecule_t() );

        read_sdf( is, *pmol );

        pmol->set_s(NAME, m_name);
    
        content().set( m_name, pmol );

        return true;        
    }

    void loadsdf_command::undo( )
    {
        content().remove( m_name.c_str() );    
    }

    const char* loadsdf_command::info( ) const
    {
        static const char usage[] = "Usage: molname = loadsdf filename || loadsdf molname filename";
        return usage;
    }
   
    shared_ptr< command_i > loadsdf_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 3 )
        {
            throw std::runtime_error( string( "Error: number of wrong arguments\n" ) + info() );
        }
    
        return shared_ptr< command_i >( new loadsdf_command( args[1], args[2] ) );
    }

}

static amber::loadsdf_command g_loadsdf_command;


