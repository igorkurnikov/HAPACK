#include <fstream>
#include <stdexcept>
#include "source.hpp"

using std::string;

using std::vector;

using std::ifstream;

namespace amber
{
    source_command::source_command( )
	:command_i( "source" )
    {
    }

    source_command::source_command( const string& file )
        :m_file( file )
    {
        m_number_cmd = 0;
    }

    source_command::~source_command()
    {
    }

    bool source_command::exec( )
    {
        ifstream stream( find_file( m_file  ).c_str() );

        string line;
    
        while( getline( stream, line ) )
        {
            console()->process( line );
        }
    
        return true;
    }

    void source_command::undo( )
    {
        throw std::runtime_error( "Sorry: not implemented yet." );
    }

    const char* source_command::info( ) const
    {
        return "  usage: source file";
    }

    shared_ptr< command_i > source_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 2 )
        {
            throw std::runtime_error( "Error: wrong number of arguments" );
        }
        
        return shared_ptr< command_i >( new source_command( args[1] ) );
    }

} // namespace amber
 
static amber::source_command g_source_command;
