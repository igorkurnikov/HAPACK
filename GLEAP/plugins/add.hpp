#ifndef GLEAP_PLUGINS_ADD_HPP
#define GLEAP_PLUGINS_ADD_HPP

#include <guilib.hpp>

namespace amber
{
    
    using namespace mort;

    class add_command : public command_i
    {
    public:

        add_command( );
        
        add_command( const string& dst, const string& src );
        
        ~add_command( );
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;

    private:

        string m_dst;
        
        string m_src;

    };

    class null_command : public command_i
    {
    public:

        null_command();

        null_command( const string& name );

        virtual ~null_command();
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
    };

} // namespace amber

#endif
