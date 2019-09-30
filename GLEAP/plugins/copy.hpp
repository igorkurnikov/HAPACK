#ifndef GLEAP_PLUGINS_COPY_HPP
#define GLEAP_PLUGINS_COPY_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class copy_command : public command_i
    {
    public:

        copy_command( );
        
        copy_command( const string& dst, const string& src );
        
        virtual ~copy_command( );
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
        
    private:

        string m_src;
        
        string m_dst;
        
    };

} // namespace amber

#endif
