#ifndef MOLVIEW_COMMAND_LOADCRD_H
#define MOLVIEW_COMMAND_LOADCRD_H

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class loadcrd_command : public command_i
    {
    public:

        loadcrd_command( const string& action );

        loadcrd_command( const string& action, const string& name, const string& file );

        virtual ~loadcrd_command();

        virtual bool exec();
    
        virtual void undo();
    
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;

    private:

        string m_action;

        string m_name;
    
        string m_file;

    };
    
} // namespace amber

#endif

