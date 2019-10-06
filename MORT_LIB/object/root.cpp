#include <stdexcept>
#include <common/numvec.hpp>
#include <common/hashcode.hpp>
#include "root.hpp"
#include "mole.hpp"

namespace mort
{
    using std::logic_error;

    root_t::root_t()
    {
    }
    
    root_t::root_t( const root_t& rhs )
        :m_icontent(rhs.m_icontent), m_dcontent(rhs.m_dcontent),
	 m_scontent(rhs.m_scontent), m_vcontent(rhs.m_vcontent),
	 m_acontent(rhs.m_acontent)
    {
    }
    
    root_t& root_t::operator=( const root_t& rhs )
    {
        if( &rhs != this )
        {
            root_t temp( rhs );
            swap( temp );
        }
        
        return *this;
    }
    
    root_t::~root_t()
    {
    }
    
    void root_t::swap( root_t& rhs )
    {
        if( &rhs != this )
        {
            m_icontent.swap( rhs.m_icontent );
            m_dcontent.swap( rhs.m_dcontent );
	    m_scontent.swap( rhs.m_scontent );
	    m_vcontent.swap( rhs.m_vcontent );
	    m_acontent.swap( rhs.m_acontent );
        }
    }
    
    // setter functions
    void root_t::set_i(const hashid_t& parmid, const int& v) 
    { 
        m_icontent[parmid] = v; 
    }

    void root_t::set_d(const hashid_t& parmid, const double& v) 
    { 
        m_dcontent[parmid] = v; 
    }

    void root_t::set_s(const hashid_t& parmid, const string& v) 
    { 
        m_scontent[parmid] = v; 
    }

    void root_t::set_v(const hashid_t& parmid, const numvec& v) 
    {
        map<hashid_t,numvec>::iterator i = m_vcontent.find( parmid );
        if( i==m_vcontent.end() )
        {
            m_vcontent.insert( std::make_pair(parmid, v) );
        }
        else
        {
            i->second = v;
        }
    }

    void root_t::set_a(const hashid_t& parmid, const any& v) 
    { 
        m_acontent[parmid] = v; 
    }

    void root_t::set_i(const string& parmname, const int& value)
    {
        set_i( mort::hash(parmname), value );
    }

    void root_t::set_d(const string& parmname, const double& value)
    {
        set_d( mort::hash(parmname), value );
    }

    void root_t::set_s(const string& parmname, const string& value)
    {
        set_s( mort::hash(parmname), value );
    }

    void root_t::set_v(const string& parmname, const numvec& value)
    {
        set_v( mort::hash(parmname), value );
    }

    void root_t::set_a(const string& parmname, const any& value)
    {
        set_a( mort::hash(parmname), value );
    }

    // type I getter function
    bool root_t::get_i(const hashid_t& parmid, int& v) const
    {
        const int* pv = NULL;
	bool result = get_iptr(parmid, pv);
	if(result) v = *pv;
	return result;
    }

    bool root_t::get_d(const hashid_t& parmid, double& v) const
    {
        const double* pv = NULL;
	bool result = get_dptr(parmid, pv);
	if(result) v = *pv;
	return result;
    }

    bool root_t::get_s(const hashid_t& parmid, string& v) const
    {
        const string* pv = NULL;
	bool result = get_sptr(parmid, pv);
	if(result) v = *pv;
	return result;
    }

    bool root_t::get_v(const hashid_t& parmid, numvec& v) const
    {
        const numvec* pv = NULL;
	bool result = get_vptr(parmid, pv);
	if(result) 
        {
            v.resize( pv->size() );
            v = *pv;
        }
        
	return result;
    }

    bool root_t::get_a(const hashid_t& parmid, any& v) const
    {
        const any* pv = NULL;
	bool result = get_aptr(parmid, pv);
	if(result) v = *pv;
	return result;
    }

    bool root_t::get_i(const string& parmname, int& value) const
    {
        return get_i(mort::hash(parmname), value);
    }

    bool root_t::get_d(const string& parmname, double& value) const
    {
        return get_d(mort::hash(parmname), value);
    }

    bool root_t::get_s(const string& parmname, string& value) const
    {
        return get_s(mort::hash(parmname), value);
    }

    bool root_t::get_v(const string& parmname, numvec& value) const
    {
        return get_v(mort::hash(parmname), value);
    }

    bool root_t::get_a(const string& parmname, any& value) const
    {
        return get_a(mort::hash(parmname), value);
    }

    bool root_t::get_iptr(const hashid_t& parmid, int*& pv)
    {
        map<hashid_t,int>::iterator i = m_icontent.find(parmid);
	if( i == m_icontent.end() )  return false;
        pv = &(i->second);
	return true;
    }

    bool root_t::get_dptr(const hashid_t& parmid, double*& pv)
    {
        map<hashid_t,double>::iterator i = m_dcontent.find(parmid);
	if( i == m_dcontent.end() )  return false;
        pv = &(i->second);
	return true;
    }

    bool root_t::get_sptr(const hashid_t& parmid, string*& pv)
    {
        map<hashid_t,string>::iterator i = m_scontent.find(parmid);
	if( i == m_scontent.end() )  return false;
        pv = &(i->second);
	return true;
    }

    bool root_t::get_vptr(const hashid_t& parmid, numvec*& pv)
    {
        map<hashid_t,numvec>::iterator i = m_vcontent.find(parmid);
	if( i == m_vcontent.end() )  return false;
        pv = &(i->second);
	return true;
    }

    bool root_t::get_aptr(const hashid_t& parmid, any*& pv)
    {
        map<hashid_t,any>::iterator i = m_acontent.find(parmid);
	if( i == m_acontent.end() )  return false;
        pv = &(i->second);
	return true;
    }

    bool root_t::get_iptr(const string& parmname, int*& pv)
    {
        return get_iptr(mort::hash(parmname), pv);
    }

    bool root_t::get_dptr(const string& parmname, double*& pv)
    {
        return get_dptr(mort::hash(parmname), pv);
    }

    bool root_t::get_sptr(const string& parmname, string*& pv)
    {
        return get_sptr(mort::hash(parmname), pv);
    }

    bool root_t::get_vptr(const string& parmname, numvec*& pv)
    {
        return get_vptr(mort::hash(parmname), pv);
    }

    bool root_t::get_aptr(const string& parmname, any*& pv)
    {
        return get_aptr(mort::hash(parmname), pv);
    }
  
    bool root_t::get_iptr(const hashid_t& parmid, const int*& pv) const
    {
        map<hashid_t,int>::const_iterator i = m_icontent.find(parmid);
	if( i == m_icontent.end() )  return false;
        pv = &(i->second);
	return true;
    }

    bool root_t::get_dptr(const hashid_t& parmid, const double*& pv) const
    {
        map<hashid_t,double>::const_iterator i = m_dcontent.find(parmid);
	if( i == m_dcontent.end() )  return false;
        pv = &(i->second);
	return true;
    }

    bool root_t::get_sptr(const hashid_t& parmid, const string*& pv) const
    {
        map<hashid_t,string>::const_iterator i = m_scontent.find(parmid);
	if( i == m_scontent.end() )  return false;
        pv = &(i->second);
	return true;
    }

    bool root_t::get_vptr(const hashid_t& parmid, const numvec*& pv) const
    {
        map<hashid_t,numvec>::const_iterator i = m_vcontent.find(parmid);
	if( i == m_vcontent.end() )  return false;
        pv = &(i->second);
	return true;
    }

    bool root_t::get_aptr(const hashid_t& parmid, const any*& pv) const
    {
        map<hashid_t,any>::const_iterator i = m_acontent.find(parmid);
	if( i == m_acontent.end() )  return false;
        pv = &(i->second);
	return true;
    }

    bool root_t::get_iptr(const string& parmname, const int*& pv) const
    {
        return get_iptr(mort::hash(parmname), pv);
    }

    bool root_t::get_dptr(const string& parmname, const double*& pv) const
    {
        return get_dptr(mort::hash(parmname), pv);
    }

    bool root_t::get_sptr(const string& parmname, const string*& pv) const
    {
        return get_sptr(mort::hash(parmname), pv);
    }

    bool root_t::get_vptr(const string& parmname, const numvec*& pv) const
    {
        return get_vptr(mort::hash(parmname), pv);
    }

    bool root_t::get_aptr(const string& parmname, const any*& pv) const
    {
        return get_aptr(mort::hash(parmname), pv);
    }

    // type II getter function
    int& root_t::get_i(const hashid_t& parmid)
    {
        int* pv = NULL;

	if( ! get_iptr(parmid, pv) )
	    throw logic_error( "Error: undefined integer paramter " + unhash(parmid) + " of entity");

	return *pv;
    }

    double& root_t::get_d(const hashid_t& parmid)
    {
        double* pv = NULL;

	if( ! get_dptr(parmid, pv) )
	    throw logic_error( "Error: undefined double paramter " + unhash(parmid) );

	return *pv;
    }

    string& root_t::get_s(const hashid_t& parmid)
    {
        string* pv = NULL;

	if( get_sptr(parmid, pv) )
            return *pv;

        static string noname = "noname";
        if( parmid==NAME )
            return noname;

        throw logic_error( "Error: undefined string paramter " + unhash(parmid) + " of entity");
    }

    numvec& root_t::get_v(const hashid_t& parmid)
    {
        numvec* pv = NULL;
	if( ! get_vptr(parmid, pv) )
	    throw logic_error( "Error: undefined numvec paramter " + unhash(parmid) );

	return *pv;
    }

    any& root_t::get_a(const hashid_t& parmid)
    {
        any* pv = NULL;

	if( ! get_aptr(parmid, pv) )
	    throw logic_error( "Error: undefined any other typed paramter " + unhash(parmid) );

	return *pv;
    }

    int& root_t::get_i(const string& parmname)
    {
        return get_i( mort::hash(parmname) );
    }

    double& root_t::get_d(const string& parmname)
    {
        return get_d( mort::hash(parmname) );
    }

    string& root_t::get_s(const string& parmname)
    {
        return get_s( mort::hash(parmname) );
    }

    numvec& root_t::get_v(const string& parmname)
    {
        return get_v( mort::hash(parmname) );
    }

    any& root_t::get_a(const string& parmname)
    {
        return get_a( mort::hash(parmname) );
    }

    const int& root_t::get_i(const hashid_t& parmid) const
    {
        const int* pv = NULL;

	if( ! get_iptr(parmid, pv) )
	    throw logic_error( "Error: undefined integer paramter " + unhash(parmid) );

	return *pv;
    }

    const double& root_t::get_d(const hashid_t& parmid) const
    {
        const double* pv = NULL;

	if( ! get_dptr(parmid, pv) )
	    throw logic_error( "Error: undefined double paramter " + unhash(parmid) );

	return *pv;
    }

    const string& root_t::get_s(const hashid_t& parmid) const
    {
        const string* pv = NULL;

	if( get_sptr(parmid, pv) )
            return *pv;

        static const string noname="noname";
        if( parmid==NAME )
            return noname;

        throw logic_error( "Error: undefined string paramter " + unhash(parmid) + " of entity" );
    }

    const numvec& root_t::get_v(const hashid_t& parmid) const
    {
        const numvec* pv = NULL;
	if( ! get_vptr(parmid, pv) )
	    throw logic_error( "Error: undefined numvec paramter " + unhash(parmid) );

	return *pv;
    }

    const any& root_t::get_a(const hashid_t& parmid) const
    {
        const any* pv = NULL;

	if( ! get_aptr(parmid, pv) )
	    throw logic_error( "Error: undefined any other typed paramter " + unhash(parmid) );

	return *pv;
    }

    const int& root_t::get_i(const string& parmname) const
    {
        return get_i( mort::hash(parmname) );
    }

    const double& root_t::get_d(const string& parmname) const
    {
        return get_d( mort::hash(parmname) );
    }

    const string& root_t::get_s(const string& parmname) const
    {
        return get_s( mort::hash(parmname) );
    }

    const numvec& root_t::get_v(const string& parmname) const
    {
        return get_v( mort::hash(parmname) );
    }

    const any& root_t::get_a(const string& parmname) const
    {
        return get_a( mort::hash(parmname) );
    }

    bool root_t::has_i(const hashid_t& parmid) const
    {
        const int* pv=NULL;
        return get_iptr(parmid, pv);
    }
    
    bool root_t::has_d(const hashid_t& parmid) const
    {
        const double* pv=NULL;
	return get_dptr(parmid, pv);
    }

    bool root_t::has_s(const hashid_t& parmid) const
    {
        const string* pv=NULL;
	return get_sptr(parmid, pv);
    }

    bool root_t::has_v(const hashid_t& parmid) const
    {
        const numvec* pv=NULL;
	return get_vptr(parmid, pv);
    }

    bool root_t::has_a(const hashid_t& parmid) const
    {
        const any* pv=NULL;
	return get_aptr(parmid, pv);
    }

    bool root_t::has_i(const string& parmname) const
    {
        return has_i( mort::hash(parmname) );
    }

    bool root_t::has_d(const string& parmname) const
    {
        return has_d( mort::hash(parmname) );
    }

    bool root_t::has_s(const string& parmname) const
    {
        return has_s( mort::hash(parmname) );
    }

    bool root_t::has_v(const string& parmname) const
    {
        return has_v( mort::hash(parmname) );
    }

    bool root_t::has_a(const string& parmname) const
    {
        return has_a( mort::hash(parmname) );
    }

    map<hashid_t,int>::const_iterator root_t::ibegin() const
    {
        return m_icontent.begin();
    }

    map<hashid_t,int>::const_iterator root_t::iend() const
    {
        return m_icontent.end();
    }

    map<hashid_t,double>::const_iterator root_t::dbegin() const
    {
        return m_dcontent.begin();
    }

    map<hashid_t,double>::const_iterator root_t::dend() const
    {
        return m_dcontent.end();
    }

    map<hashid_t,string>::const_iterator root_t::sbegin() const
    {
        return m_scontent.begin();
    }

    map<hashid_t,string>::const_iterator root_t::send() const
    {
        return m_scontent.end();
    }

    map<hashid_t,numvec>::const_iterator root_t::vbegin() const
    {
        return m_vcontent.begin();
    }

    map<hashid_t,numvec>::const_iterator root_t::vend() const
    {
        return m_vcontent.end();
    }

    map<hashid_t,any>::const_iterator root_t::abegin() const
    {
        return m_acontent.begin();
    }

    map<hashid_t,any>::const_iterator root_t::aend() const
    {
        return m_acontent.end();
    }

    bool ismol( const entity_ptr& pe )    
    {
        molecule_ptr pm = dynamic_pointer_cast< molecule_t >( pe );
	return pm != NULL;
    }


} // namespace mort


