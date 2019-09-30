#ifndef MOLVIEW_COMMAND_LOADMOL2_H
#define MOLVIEW_COMMAND_LOADMOL2_H

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class loadmol2_command : public command_i
    {
    public:

        loadmol2_command( );

        loadmol2_command( const string& name, const string& file );

        virtual ~loadmol2_command();

        virtual bool exec();
    
        virtual void undo();
 
        virtual const char* info() const;
   
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;

    private:

        string m_name;
    
        string m_file;

    };
    
} // namespace amber

#endif

