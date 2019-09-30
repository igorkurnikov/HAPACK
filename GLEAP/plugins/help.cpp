#include "help.hpp"

namespace amber
{
    using std::map;

    help_command::help_command()
        :command_i( "help" )
    {
    }

    help_command::help_command( const string& name )
        :m_name( name )
    {
    }

    help_command::~help_command()
    {
    }

    bool help_command::exec()
    {
        if( m_name.empty() )
        {
            map< string, command_i* >::iterator i = control_t::begin();
            for( ; i != control_t::end(); ++i )
            {
                leaplog_t::putline( i->first );
            }
        }
        else
        {
            leaplog_t::putline( "sorry, not implemented yet!" );
        }

	return true;
    }

    void help_command::undo()
    {
    }

    shared_ptr< command_i > help_command::clone( const vector< string >& args ) const
    {
        string arg = ( args.size()>1 ) ? args[1] : "";

        return shared_ptr< command_i >( new help_command( arg ) );
    }

} // namespace amber

//amber::help_command g_help_command;

