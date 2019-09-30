#ifndef GLEAP_PLUGINS_ENERGY_HPP
#define GLEAP_PLUGINS_ENERGY_HPP

#include <guilib.hpp>


namespace amber
{
    using namespace mort;

    class energy_command : public command_i
    {
    public:

        energy_command();
        
        energy_command(const string& name, const string& parm);
 
        energy_command(const string& name, const string& parm, const string& vec1, const string& vec2);

        virtual ~energy_command();
        
        virtual bool exec();
        
        virtual void undo();
        
        virtual const char* info() const;

        virtual command_ptr clone(const vector<string>& args) const;
        
    private:

        string m_name;

	string m_parm;

	string m_vec1;

	string m_vec2;
    };
    
} // namespace amber


#endif



        
