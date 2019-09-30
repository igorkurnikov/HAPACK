#ifndef GLEAP_PLUGINS_SET_HPP
#define GLEAP_PLUGINS_SET_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class set_command : public command_i
    {
    public:

        set_command( );
        
        set_command( const string& object, const string& parm, const string& value );
        
        virtual ~set_command( );
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;

    private:

        string m_object;
        
        string m_parm;
        
        string m_value;
    };

} // namespace amber

    
#endif
