#include "shell.hpp"

namespace amber
{
    shell_command::shell_command( )
       	:command_i( "shell" )
    {
    }
    
    shell_command::shell_command( const string& cmd )
        :m_command( cmd )
    {
    }

    shell_command::~shell_command()
    {
    }
    
    bool shell_command::exec()
    {
        if( empty( m_command ) )
        {
            system( "bash" );
        }
        else
        {
            system( m_command.c_str() );
        }

	return true;
    }
    
    void shell_command::undo()
    {
        leaplog_t::putline( "Error: can't undo shell command " + m_command );
    }
    
    shared_ptr< command_i > shell_command::clone( const vector< string >& args ) const
    {
        string cmd;

        for( int i=0; i < (int)args.size(); ++i )
        {
            cmd += args[i] + " ";
        }

        return shared_ptr< command_i >( new shell_command( cmd ) );
    }

} // namespace amber

amber::shell_command g_shell_command;
