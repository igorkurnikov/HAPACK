#include "command.hpp"
#include "control.hpp"

namespace mort
{
    using std::map;

    command_i::command_i()
    {
    }
    
    command_i::command_i( const string& name )
    {
        insert( name, this );
    }
 
    command_i::~command_i()
    {
    }

    map<string, command_i*>& command_i::dict()
    {
        static map<string, command_i*> g_commands;
        return g_commands;
    }

    void command_i::insert( const string& name, command_i* inst )
    {
        dict()[name] = inst;
    }

    command_i* command_i::find( const string& name )
    {
        map<string, command_i*>::iterator i = dict().find( name );
        if( i==dict().end() )
        {
            throw std::runtime_error( name + ": command not found!" );
        }

        return i->second;
    }

} // namespace amber

