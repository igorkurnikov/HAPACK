#ifndef GLEAP_PLUGINS_DELETE_HPP
#define GLEAP_PLUGINS_DELETE_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class delete_command : public command_i
    {
    public:

        delete_command( const string& type );
        
        delete_command( const string& action, const string& oper1, const string& oper2 );
        
        virtual ~delete_command( );
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual command_ptr clone( const vector< string >& args ) const;

    private:

        string m_action;
        
        string m_oper1;
        
        string m_oper2;
    };
}

#endif
