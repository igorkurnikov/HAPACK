#ifndef GTKLEAP_COMMAND_SHELL_H
#define GTKLEAP_COMMAND_SHELL_H

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class shell_command : public command_i
    {
    public:

        shell_command();
        
        shell_command( const string& cmd );
        
        virtual ~shell_command();
        
        virtual bool exec();
        
        virtual void undo();
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;

    private:

        string m_command;
    };
    
} // namespace amber


#endif
