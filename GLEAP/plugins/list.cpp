#include <object.hpp>
#include "list.hpp"

namespace amber
{
    list_command::list_command( )
	: command_i( "list" )
    {
    }

    list_command::list_command( const string& )
    {
    }

    list_command::~list_command()
    {
    }

    bool list_command::exec()
    {
        database_t::iterator i = content().begin();
        database_t::iterator e = content().end();
    
        string line;

        for( ; i != e; ++i )
        {
            line += i->first;
            line += " ";

            if( (int)line.length() >= MAX_LINE_WIDTH )
            {
                leaplog_t::putline( line );
                line.clear();
            }
        }

        leaplog_t::putline( line );
    
        return false;
    }

    void list_command::undo()
    {
    }

    shared_ptr< command_i > list_command::clone( const vector< string >& ) const
    {
        return shared_ptr< command_i >( new list_command( "type" ) );
    }

} // namespace amber

amber::list_command g_list_command;    
    

    
