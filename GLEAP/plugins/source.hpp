#ifndef MOLVIEW_COMMAND_SOURCE_H
#define MOLVIEW_COMMAND_SOURCE_H

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class source_command : public command_i
    {
    public:
    
        source_command( );
    
        source_command( const std::string& m_filename );

        virtual ~source_command();
    
        virtual bool exec( );
    
        virtual void undo( );

        virtual const char* info() const;

        virtual shared_ptr< command_i > clone( const std::vector< std::string >& args ) const;
    
    private:

        std::string m_file;
    
        int m_number_cmd;
    };
    
} // namespace amber

#endif

