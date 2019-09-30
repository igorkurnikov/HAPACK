#ifndef GTKLEAP_COMMAND_LOADPREP_H
#define GTKLEAP_COMMAND_LOADPREP_H

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class loadprep_command : public command_i
    {
    public:

        loadprep_command();
    
        loadprep_command(const string& file, const string& prefix);
    
        virtual ~loadprep_command();
    
        virtual bool exec();
    
        virtual void undo();
    
        virtual command_ptr clone(const vector<string>& args) const;
    
    private:

        string m_file;
 
        string m_prefix;

        vector< string > m_names;
    };
    
} // namespace amber

    
#endif
