#ifndef GLEAP_PLUGINS_HELP_H
#define GLEAP_PLUGINS_HELP_H

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class help_command : public command_i
    {
    public:

        help_command();

        help_command( const string& name );

	virtual ~help_command( );

        virtual bool exec();

        virtual void undo();

	virtual shared_ptr< command_i > clone( const vector< string >& args ) const;

    private:

        std::string m_name;
    };

}  // namespace amber

#endif
   
