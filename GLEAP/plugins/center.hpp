#ifndef GLEAP_PLUGINS_CENTER_HPP
#define GLEAP_PLUGINS_CENTER_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class center_command : public command_i
    {
    public:

        center_command( );
        
        center_command( const string& object );
        
        virtual ~center_command( );
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
        
    private:

        string m_object;
        
    };

} // namespace amber

#endif
