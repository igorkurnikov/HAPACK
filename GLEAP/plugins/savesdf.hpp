#ifndef MOLVIEW_COMMAND_SAVESDF_H
#define MOLVIEW_COMMAND_SAVESDF_H


#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class savesdf_command : public command_i
    {
    public:

        savesdf_command( );

        savesdf_command( const string& name, const string& file );

        virtual ~savesdf_command();

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



