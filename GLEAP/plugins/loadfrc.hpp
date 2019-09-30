#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class loadfrc_command : public command_i
    {
    public:

        loadfrc_command( const string& name );
        
        loadfrc_command( const string& action, const string& file );
        
        virtual ~loadfrc_command();
        
        virtual bool exec();
        
        virtual void undo();
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
        
    private:

        string m_file;

        string m_action;
    };

} // namespace amber

        
