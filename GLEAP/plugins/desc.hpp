#ifndef GLEAP_PLUGINS_DESC_HPP
#define GLEAP_PLUGINS_DESC_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class desc_command : public command_i
    {
    public:

        desc_command( );
        
        desc_command( const string& object );
        
        virtual ~desc_command( );
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
        
    private:

        string m_object;
        
    };

} // namespace amber

#endif


