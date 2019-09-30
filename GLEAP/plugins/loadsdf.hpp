#ifndef MOLVIEW_COMMAND_LOADSDF_H
#define MOLVIEW_COMMAND_LOADSDF_H


#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class loadsdf_command : public command_i
    {
    public:

        loadsdf_command( );

        loadsdf_command( const string& name, const string& file );

        virtual ~loadsdf_command();

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



