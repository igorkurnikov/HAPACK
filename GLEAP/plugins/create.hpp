#ifndef GLEAP_PLUGINS_CREATE_HPP
#define GLEAP_PLUGINS_CREATE_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class create_command : public command_i
    {
    public:

        create_command( const string& object );

        create_command( const string& object, const string& dest, const string& name );
        
        create_command( const string& object, const string& dest, const string& name, const string& type, const string& charge );

        virtual ~create_command();
        
        virtual bool exec();
        
        virtual void undo();
        
        virtual const char* info() const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
        
    private:

        string m_object;

        string m_dest;

        string m_name;
        
        string m_type;
        
        string m_chrg;
    };

}

#endif

