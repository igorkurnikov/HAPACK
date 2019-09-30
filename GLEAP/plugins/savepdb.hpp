#ifndef MOLVIEW_COMMAND_SAVEPDB_H
#define MOLVIEW_COMMAND_SAVEPDB_H


#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class savepdb_command : public command_i
    {
    public:

        savepdb_command( );

        savepdb_command( const string& name, const string& file );

        virtual ~savepdb_command();

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



