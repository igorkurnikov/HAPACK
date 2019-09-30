#ifndef GLEAP_LEAPSRC_LOADLIB_HPP
#define GLEAP_LEAPSRC_LOADLIB_HPP

#include <guilib/command.hpp>

using namespace mort;

namespace amber
{
    class loadlib_command : public command_i
    {
    public:

        loadlib_command( );
        
        loadlib_command( const string& name );
        
        virtual ~loadlib_command( );
        
        virtual bool exec( );
        
        virtual void undo( );

        virtual const char* info( ) const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
        
    private:

        string m_name;
    };

} // namespace amber

#endif
