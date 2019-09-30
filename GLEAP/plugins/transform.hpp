#ifndef GLEAP_PLUGINS_TRANSFORM_HPP
#define GLEAP_PLUGINS_TRANSFORM_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class transform_command : public command_i
    {
    public:

        transform_command( const string& action );
        
        transform_command( const string& action, const string& object, const string& offset );
        
        virtual ~transform_command( );
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual command_ptr clone( const vector< string >& args ) const;
        
    private:

        string m_action;

        string m_object;
        
        string m_offset;
    };
    
} // namespace amber

#endif

