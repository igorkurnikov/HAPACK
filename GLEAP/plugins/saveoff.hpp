#ifndef GLEAP_PLUGINS_SAVEOFF_HPP
#define GLEAP_PLUGINS_SAVEOFF_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class saveoff_command : public command_i
    {
    public:

        saveoff_command( );
        
        saveoff_command( const string& unit, const string& file );
        
        virtual ~saveoff_command( );
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
        
    private:

        string m_unit;
        
        string m_file;
    };

} // namespace amber


#endif
