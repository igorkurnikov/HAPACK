#ifndef GLEAP_COMMAND_ADDIONS_HPP
#define GLEAP_COMMAND_ADDIONS_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class addions_command : public command_i
    {
    public:

        addions_command();
        
        addions_command( const string& dest, const string& ion, int number );
        
        virtual ~addions_command();
        
        virtual bool exec();
        
        virtual void undo();
        
        virtual const char* info( ) const;

        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
    
    private:

        string m_dest;

        string m_ion;
        
        int m_number;
    };
    
} // namespace amber

#endif


