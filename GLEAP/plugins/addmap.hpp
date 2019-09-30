#ifndef GLEAP_COMMAND_ADDMAP_HPP
#define GLEAP_COMMAND_ADDMAP_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class addmap_command : public command_i
    {
    public:

        addmap_command(const string& type);

	addmap_command(const string& type, const string& nmap );

	virtual ~addmap_command();

	virtual bool exec();

	virtual void undo();

	virtual const char* info() const;

	virtual shared_ptr<command_i> clone(const vector<string>& args ) const;

    private:

        string m_type;

	vector< vector<string> > m_nmap;
    };

} // namespace amber

#endif
