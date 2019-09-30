#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class saveprm_command : public command_i
    {
    public:

        saveprm_command( const string& name );
        
        saveprm_command( const string& action, const string& name, const string& top, const string& xyz );
        
        virtual ~saveprm_command();
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
        
    private:

        string m_name;
        
        string m_action;
        
        string m_topfile;
        
        string m_xyzfile;
    };
    
} // namespace amber

        
