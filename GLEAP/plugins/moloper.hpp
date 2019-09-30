#ifndef GLEAP_PLUGINS_MOLOPER_HPP
#define GLEAP_PLUGINS_MOLOPER_HPP

#include <common.hpp>
#include <guilib.hpp>

namespace amber
{
    using namespace mort;
    
    class moloper_command : public command_i
    {
    public:

        moloper_command( const string& action);
        
        moloper_command( const string& action, const string& mname );

        virtual ~moloper_command();

        virtual bool exec();
        
        virtual void undo();
        
        virtual command_ptr clone( const vector<string>& args ) const;
        
    private:

        string m_action;

        molecule_ptr m_pmol;
    };
    

} // namespace amber


#endif
