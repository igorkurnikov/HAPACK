#ifndef GLEAP_PLUGINS_COMBINE_HPP
#define GLEAP_PLUGINS_COMBINE_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class merge_command : public command_i
    {
    public:

        merge_command( const string& action );
        
        merge_command( const string& action, const string& name, const string& list );
        
        virtual ~merge_command( );
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
        
    private:

        string m_action;

        string m_name;
        
        string m_list;
    };
    
} // namespace amber


#endif
