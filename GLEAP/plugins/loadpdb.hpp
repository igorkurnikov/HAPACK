#ifndef MOLVIEW_COMMAND_LOADPDB_H
#define MOLVIEW_COMMAND_LOADPDB_H


#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class loadpdb_command : public command_i
    {
    public:

        loadpdb_command( const string& action );

        loadpdb_command( const string& action, const string& name, const string& file, const string& seq );

        virtual ~loadpdb_command();

        virtual bool exec();
    
        virtual void undo();
    
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;

    private:

        string m_action;
        
        string m_name;
    
        string m_file;

        string m_seq;
    };
    
} // namespace amber

#endif



