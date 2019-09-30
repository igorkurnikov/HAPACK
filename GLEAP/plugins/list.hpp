#ifndef GTKLEAP_COMMAND_LIST_H
#define GTKLEAP_COMMAND_LIST_H

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class list_command : public command_i
    {
    public:

        list_command();
    
        list_command( const string& type );
    
        virtual ~list_command();
    
        virtual bool exec();
    
        virtual void undo();
    
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
    
    };
    
} // namespace amber
    
#endif
