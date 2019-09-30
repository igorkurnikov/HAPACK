#ifndef GLEAP_PLUGINS_CHARGE_HPP
#define GLEAP_PLUGINS_CHARGE_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class charge_command : public command_i
    {
    public:

        charge_command( );
        
        charge_command( const string& object );
        
        virtual ~charge_command( );
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
        
    private:

        string m_object;
        
    };

} // namespace amber

#endif
