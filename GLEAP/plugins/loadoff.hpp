#ifndef GTKLEAP_COMMAND_LOADOFF_H
#define GTKLEAP_COMMAND_LOADOFF_H

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class loadoff_command : public command_i
    {
    public:

        loadoff_command( );
    
        loadoff_command( const string& file, const string& name );
    
        virtual ~loadoff_command( );
    
        virtual bool exec( );
    
        virtual void undo( );
    
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
    
    private:

        string m_file;

        string m_name;
    
        vector< string > m_names;
    };
    
} // namespace amber

    
#endif
