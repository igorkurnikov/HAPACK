#include <ltdl.h>
#include <sys/types.h>
#include <dirent.h>

#include <guilib/leaplog.hpp>
#include <guilib/control.hpp>
#include "loadlib.hpp"

namespace amber
{
    
    loadlib_command::loadlib_command( )
        : command_i( "loadlib" )
    {
    }
    
    loadlib_command::loadlib_command( const string& name )
        : m_name( name )
    {
    }
    
    loadlib_command::~loadlib_command( )
    {
    }

    bool loadlib_command::exec( )
    {
        string file = "lib" + m_name + ".la";

        int errors = lt_dlinit();
 
        if( errors != 0 )
        {
            throw logic_error( "Error: can't initlialize libltdl" );
        }

        lt_dlhandle module = lt_dlopen( file.c_str() );

        if( module == NULL )
        {
            throw logic_error( "Error: can't open library file " + file + " because " + string( lt_dlerror() ) );
        }

        std::cout << "library commands have been loaded from " << file << std::endl;
        
        return true;
    }
    
    void loadlib_command::undo( )
    {
        throw logic_error( "Sorry, not implemented" );
    }
    
    const char* loadlib_command::info( ) const
    {
        return "  usage: loadlib module ";
    }
    
    shared_ptr< command_i > loadlib_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 2 )
        {
            throw logic_error( "Error: wrong number of arguments" );
        }
        
        return shared_ptr< command_i >( new loadlib_command( args[1] ) );
    }
    
} // namespace amber

amber::loadlib_command g_loadlib_command;

