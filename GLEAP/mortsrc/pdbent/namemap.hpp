#ifndef GLEAP_OBJECT_NAMEMAP_HPP
#define GLEAP_OBJECT_NAMEMAP_HPP

#include <vector>
#include <object.hpp>

namespace mort
{
    using std::vector;

    class morf_t;

    class namemap_t : public root_t
    {
    public:

        namemap_t();

        namemap_t(const namemap_t& rhs);

	virtual ~namemap_t();

	virtual void swap(namemap_t& rhs);

	void add_resd_map(const vector<string>& list);

	void add_atom_map(const vector<string>& list);

        string get_name( const morf_t& mo ) const;

    private:

	string get_resd_name(const string& rname, int aapos) const;

	string get_atom_name(const string& aname) const;

    private:

        map<string, string> m_resd_nterm;
	map<string, string> m_resd_cterm;
	map<string, string> m_resd;
	map<string, string> m_atom;
    };

} // namespace mort


#endif
