#ifndef MOLVIEW_COMMAND_SAVEMOL2_H
#define MOLVIEW_COMMAND_SAVEMOL2_H


#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class savemol2_command : public command_i
    {
    public:

        savemol2_command( );

        savemol2_command( const string& name, const string& file );

        virtual ~savemol2_command();

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

