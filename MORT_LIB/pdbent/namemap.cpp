#include <common.hpp>
#include <object.hpp>
#include "namemap.hpp"

namespace mort
{

    namemap_t::namemap_t()
    {
    }

    namemap_t::namemap_t(const namemap_t& rhs)
        :root_t(rhs),m_resd_nterm(rhs.m_resd_nterm),m_resd_cterm(rhs.m_resd_cterm),
	 m_resd(rhs.m_resd),m_atom(rhs.m_atom)
    {
    }

    namemap_t::~namemap_t()
    {
    }

    void namemap_t::swap(namemap_t& rhs)
    {
        root_t::swap(rhs);
	m_resd_nterm.swap(rhs.m_resd_nterm);
	m_resd_cterm.swap(rhs.m_resd_cterm);
	m_resd.swap(rhs.m_resd);
	m_atom.swap(rhs.m_atom);
    }

    void namemap_t::add_resd_map(const vector<string>& list)
    {
	if( list.size()==3 )
	{
	    assert( list[0].size()==1 );
	    char p = list[0][0];
	    assert( p=='0' || p=='1' );

            int aapos = (p=='0')?NTERM:CTERM;
	    string name_a = strip_quota(list[1]);
	    string name_b = strip_quota(list[2]);

            if(aapos==NTERM)
	    {
	        m_resd_nterm[name_a]=name_b;
	    }
	    else
	    {
	        assert(aapos==CTERM);
		m_resd_cterm[name_a]=name_b;
	    }
	}
        else
	{
	    string name_a = strip_quota(list[0]);
	    string name_b = strip_quota(list[1]);
            m_resd[name_a]=name_b;
	}
    }

    void namemap_t::add_atom_map(const vector<string>& list)
    {
        string name_a = strip_quota(list[0]);
        string name_b = strip_quota(list[1]);
        m_atom[name_a] = name_b;
    }

    string namemap_t::get_name( const morf_t& mo ) const
    {
        if( mo.cmpid()==RESD )
        {
            return get_resd_name( mo.get_s(TYPE), mo.get_i(AAPOS) );
        }

        if( mo.cmpid()==ATOM )
        {
            return get_atom_name( mo.get_s(NAME) );
        }

        throw std::runtime_error( "Error: namemap cannot get name for " + unhash(mo.cmpid()) );
    }

    string namemap_t::get_resd_name(const string& rname, int aapos) const
    {
        const map<string,string>* pmap = NULL;
	if(aapos==NTERM)
	{
	    pmap = &m_resd_nterm;
	}
	else if(aapos==CTERM)
	{
            pmap = &m_resd_cterm;
	}
	else
	{
	    pmap = &m_resd;
	}

	map<string,string>::const_iterator i = pmap->find(rname);
	if( i == pmap->end() )
	{
	    return rname;
	}
	return i->second;
    }

    string namemap_t::get_atom_name(const string& aname) const
    {
        map<string,string>::const_iterator i = m_atom.find(aname);
	if( i==m_atom.end() )
	{
	    return aname;
	}

	return i->second;
    }

}

