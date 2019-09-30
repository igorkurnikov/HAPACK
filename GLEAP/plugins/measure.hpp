#ifndef GLEAP_PLUGINS_MEASURE_HPP
#define GLEAP_PLUGINS_MEASURE_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class measure_command : public command_i
    {
    public:

        measure_command( );
        
        measure_command( const vector< string >& atoms );
        
        virtual ~measure_command( );

        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
        
    private:

        vector< string > m_atoms;
    };
    
} // namespace amber

#endif
