#include <sstream>
#include <common.hpp>
#include <object.hpp>
#include "parmfun.hpp"

namespace mort
{
    iparm_setter2::iparm_setter2(const hashid_t& parmid)
    {
        m_parmid = parmid;
    }

    void iparm_setter2::operator()( morf_t& p, int value) const
    {
        p.set_i(m_parmid, value);
    }

    dparm_setter2::dparm_setter2(const hashid_t& parmid)
    {
        m_parmid = parmid;
    }

    void dparm_setter2::operator()( morf_t& p, double value) const
    {
        p.set_d(m_parmid, value);
    }

    sparm_setter2::sparm_setter2(const hashid_t& parmid)
    {
        m_parmid = parmid;
    }

    void sparm_setter2::operator()( morf_t& p, const string& value) const
    {
        p.set_s(m_parmid, value);
    }

    vparm_setter2::vparm_setter2(const hashid_t& parmid)
    {
        m_parmid = parmid;
    }

    void vparm_setter2::operator()( morf_t& p, const numvec& value) const
    {
        p.set_v(m_parmid, value);
    }

    aparm_setter2::aparm_setter2(const hashid_t& parmid)
    {
        m_parmid = parmid;
    }

    void aparm_setter2::operator()( morf_t& p, const any& value) const
    {
        p.set_a(m_parmid, value);
    }

    iparm_setter1::iparm_setter1(const hashid_t& parmid, int value )
    {
        m_parmid = parmid;
	m_value  = value;
    }

    void iparm_setter1::operator()(morf_t& p) const
    {
        p.set_i(m_parmid, m_value);
    }

    dparm_setter1::dparm_setter1(const hashid_t& parmid, double value )
    {
        m_parmid = parmid;
	m_value  = value;
    }

    void dparm_setter1::operator()(morf_t& p) const
    {
        p.set_d(m_parmid, m_value);
    }

    sparm_setter1::sparm_setter1(const hashid_t& parmid, const string& value)
        : m_value(value)
    {
        m_parmid = parmid;
    }

    void sparm_setter1::operator()( morf_t& p ) const
    {
        p.set_s(m_parmid, m_value);
    }

    vparm_setter1::vparm_setter1(const hashid_t& parmid, const numvec& value )
        : m_parmid(parmid), m_value(value)
    {
    }

    void vparm_setter1::operator()(morf_t& p) const
    {
        p.set_v(m_parmid, m_value);
    }

    aparm_setter1::aparm_setter1(const hashid_t& parmid, const any& value )
        : m_parmid(parmid), m_value(value)
    {
    }

    void aparm_setter1::operator()(morf_t& p) const
    {
        p.set_a(m_parmid, m_value);
    }

    // all kinds of parm getters
    iparm_getter::iparm_getter(const hashid_t& parmid)
    {
        m_parmid = parmid;
    }

    int iparm_getter::operator()(const morf_t& p) const
    {
        return p.get_i(m_parmid);
    }

    dparm_getter::dparm_getter(const hashid_t& parmid)
    {
        m_parmid = parmid;
    }

    double dparm_getter::operator()(const morf_t& p) const
    {
        return p.get_d(m_parmid);
    }

    sparm_getter::sparm_getter(const hashid_t& parmid)
    {
        m_parmid = parmid;
    }

    string sparm_getter::operator()(const morf_t& p) const
    {
        return p.get_s(m_parmid);
    }

    vparm_getter::vparm_getter(const hashid_t& parmid)
    {
        m_parmid = parmid;
    }

    numvec vparm_getter::operator()(const morf_t& p) const
    {
        return p.get_v(m_parmid);
    }

    aparm_getter::aparm_getter(const hashid_t& parmid)
    {
        m_parmid = parmid;
    }

    any aparm_getter::operator()(const morf_t& p) const
    {
        return p.get_a(m_parmid);
    }


    iparm_cmper1::iparm_cmper1(const hashid_t& parmid,int value)
    {
        m_parmid = parmid;
	m_value  = value;
    }

    bool iparm_cmper1::operator()(const morf_t& p) const
    {
        return p.get_i(m_parmid)==m_value;
    }

    sparm_cmper1::sparm_cmper1(const hashid_t& parmid, const string& value)
    {
        m_parmid = parmid;
	m_value  = value;
    }

    bool sparm_cmper1::operator()(const morf_t& p) const
    {
        return p.get_s(m_parmid)==m_value;
    }


    iparm_cmper2::iparm_cmper2(const hashid_t& parmid)
    {
        m_parmid = parmid;
    }

    bool iparm_cmper2::operator()(const morf_t& p, int value) const
    {
        return p.get_i(m_parmid) == value;
    }

    void read_iparm(istream& is, morf_t& o, const hashid_t& parmid)
    {
        string v = next_word( is, '\n' );
	if( ! v.empty() )
	{
	    o.set_i(parmid, atoi(v.c_str()) );
	}
    }

    void read_dparm(istream& is, morf_t& o, const hashid_t& parmid)
    {
        string v = next_word( is, '\n' );
	if( ! v.empty() )
	{
	    o.set_d(parmid, atof(v.c_str()) );
	}
    }

    void read_sparm(istream& is, morf_t& o, const hashid_t& parmid)
    {
        string v = next_word( is, '\n' );
	if( ! v.empty() )
	{
	    o.set_s(parmid, v);
	}
    }

    void read_vparm(istream& is, morf_t& o, const hashid_t& parmid, int len)
    {
        numvec value(len);

	for( int i=0; i < len; ++i )
	{
	    is >> value[i];
	}

	o.set_v(parmid, value);
    }

    void write_iparm(ostream& os, const morf_t& o, const hashid_t& parmid, const string& f)
    {
        int i = o.get_i(parmid);
	os << format( f ) % i;
    }

    void write_dparm(ostream& os, const morf_t& o, const hashid_t& parmid, const string& f)
    {
        double d = o.get_d(parmid);
	os << format( f ) % d;
    }

    void write_sparm(ostream& os, const morf_t& o, const hashid_t& parmid, const string& f)
    {
        string s = o.get_s(parmid);
	os << format( f ) % s;
    }

    void write_vparm(ostream& os, const morf_t& o, const hashid_t& parmid, const string& f)
    {
        numvec v = o.get_v(parmid);
	for( int i=0; i < (int)v.size(); ++i )
	{
	    os << format( f ) % v[i];
	}
    }

    void fortran_write_iparm( fortran_t* fortran_format, morf_t& p, const hashid_t& parmid )
    {
        int ivalue;
	if( p.get_i(parmid, ivalue) )
	{
	    (*fortran_format)( p.get_i(parmid) );
	}
        else
	{
	    (*fortran_format)( p.get_s(parmid) );
	}

    }

    void fortran_write_dparm( fortran_t* fortran_format, morf_t& p, const hashid_t& parmid )
    {
        (*fortran_format)( p.get_d(parmid) );
    }

    void fortran_write_sparm( fortran_t* fortran_format, morf_t& p, const hashid_t& parmid )
    {
        (*fortran_format)( p.get_s(parmid) );
    }

    bool copy_iparm(const morf_t& src, const hashid_t& parmid, morf_t& dst)
    {
        int value;
	if( ! src.get_i( parmid, value) )
	{
	    return false;
	}

        dst.set_i(parmid, src.get_i(parmid) );
	return true;
    }

    bool copy_dparm(const morf_t& src, const hashid_t& parmid, morf_t& dst)
    {
        double value;
	if( ! src.get_d( parmid, value) )
	{
	    return false;
	}

        dst.set_d(parmid, src.get_d(parmid) );
	return true;
    }

    bool copy_sparm(const morf_t& src, const hashid_t& parmid, morf_t& dst)
    {
        string value;
	if( ! src.get_s( parmid, value) )
	{
	    return false;
	}

        dst.set_s(parmid, src.get_s(parmid) );
	return true;
    }

    bool copy_vparm(const morf_t& src, const hashid_t& parmid, morf_t& dst)
    {
        numvec value(0);
	if( ! src.get_v(parmid, value) )
	{
	    return false;
	}

        dst.set_v(parmid, src.get_v(parmid) );
	return true;
    }

    bool copy_aparm(const morf_t& src, const hashid_t& parmid, morf_t& dst)
    {
        any value;
	if( ! src.get_a(parmid, value) )
	{
	    return false;
	}

        dst.set_a(parmid, src.get_a(parmid) );
	return true;
    }

    void copy_allparms_except(const morf_t& src, morf_t& dst, const vector<hashid_t>& excepts)
    {
        const mcmpdata_t* compsrc = src.getcmp();
        iparm_iterator ii = compsrc->iparm_begin();
	for( ; ii != compsrc->iparm_end(); ++ii )
	{
            if( excepts.end()==std::find(excepts.begin(), excepts.end(), ii->first) )
	    {
	        copy_iparm(src, ii->first, dst);
	    }
        }

        dparm_iterator di = compsrc->dparm_begin();
	for( ; di != compsrc->dparm_end(); ++di )
	{
	    if( excepts.end()==std::find(excepts.begin(), excepts.end(), di->first) )
	    {
	        copy_dparm(src, di->first, dst);
	    }
        }

        sparm_iterator si = compsrc->sparm_begin();
	for( ; si != compsrc->sparm_end(); ++si )
	{
	    if( excepts.end()==std::find(excepts.begin(), excepts.end(), si->first) )
	    {
	        copy_sparm(src, si->first, dst);
	    }
        }

        vparm_iterator vi = compsrc->vparm_begin();
	for( ; vi != compsrc->vparm_end(); ++vi )
	{
	    if( excepts.end()==std::find(excepts.begin(), excepts.end(), vi->first) )
	    {
	        copy_vparm(src, vi->first, dst);
	    }
        }

        aparm_iterator ai = compsrc->aparm_begin();
	for( ; ai != compsrc->aparm_end(); ++ai )
	{
	    if( excepts.end()==std::find(excepts.begin(), excepts.end(), ai->first) )
	    {
	        copy_aparm(src, ai->first, dst);
	    }
        }
    }

    void copy_allparms_except(const morf_t& src, morf_t& dst, const hashid_t& except)
    {
        const mcmpdata_t* compsrc = src.getcmp();
        iparm_iterator ii = compsrc->iparm_begin();
	for( ; ii != compsrc->iparm_end(); ++ii )
	{
            if( except != ii->first )
	    {
	        copy_iparm(src, ii->first, dst);
	    }
        }

        dparm_iterator di = compsrc->dparm_begin();
	for( ; di != compsrc->dparm_end(); ++di )
	{
	    if( except != di->first )
	    {
	        copy_dparm(src, di->first, dst);
	    }
        }

        sparm_iterator si = compsrc->sparm_begin();
	for( ; si != compsrc->sparm_end(); ++si )
	{
	    if( except != si->first )
	    {
	        copy_sparm(src, si->first, dst);
	    }
        }

        vparm_iterator vi = compsrc->vparm_begin();
	for( ; vi != compsrc->vparm_end(); ++vi )
	{
	    if( except != vi->first )
	    {
	        copy_vparm(src, vi->first, dst);
	    }
        }

        aparm_iterator ai = compsrc->aparm_begin();
	for( ; ai != compsrc->aparm_end(); ++ai )
	{
	    if( except != ai->first )
	    {
	        copy_aparm(src, ai->first, dst);
	    }
        }
    }

    void copy_allparms(const morf_t& src, morf_t& dst)
    {
        const mcmpdata_t* compsrc = src.getcmp();
        iparm_iterator ii = compsrc->iparm_begin();
	for( ; ii != compsrc->iparm_end(); ++ii )
	{
            copy_iparm(src, ii->first, dst);
        }

        dparm_iterator di = compsrc->dparm_begin();
	for( ; di != compsrc->dparm_end(); ++di )
	{
            copy_dparm(src, di->first, dst);
        }

        sparm_iterator si = compsrc->sparm_begin();
	for( ; si != compsrc->sparm_end(); ++si )
	{
            copy_sparm(src, si->first, dst);
        }

        vparm_iterator vi = compsrc->vparm_begin();
	for( ; vi != compsrc->vparm_end(); ++vi )
	{
            copy_vparm(src, vi->first, dst);
        }

        aparm_iterator ai = compsrc->aparm_begin();
	for( ; ai != compsrc->aparm_end(); ++ai )
	{
            copy_aparm(src, ai->first, dst);
        }
    }

    void copy_allparms(const root_t& src, root_t& dst)
    {
        map<hashid_t, int>::const_iterator ii = src.ibegin();
	for( ; ii != src.iend(); ++ii )
	{
            dst.set_i( ii->first, ii->second );
        }

        map<hashid_t, double>::const_iterator di = src.dbegin();
	for( ; di != src.dend(); ++di )
	{
            dst.set_d( di->first, di->second );
        }

        map<hashid_t, string>::const_iterator si = src.sbegin();
	for( ; si != src.send(); ++si )
	{
            dst.set_s( si->first, si->second );
        }

        map<hashid_t, numvec>::const_iterator vi = src.vbegin();
	for( ; vi != src.vend(); ++vi )
	{
            dst.set_v( vi->first, vi->second );
        }

        map<hashid_t, any>::const_iterator ai = src.abegin();
	for( ; ai != src.aend(); ++ai )
	{
            dst.set_a( ai->first, ai->second );
        }
    }

    string parm2str(const morf_t& mo, const hashid_t& parmid)
    {
        std::ostringstream os;

        int iv;
	if( mo.get_i(parmid, iv) )
	{
	    os << iv;
	    return os.str();
	}

	double dv;
	if( mo.get_d(parmid, dv) )
	{
	    os << dv;
	    return os.str();
	}

	string sv;
	if( mo.get_s(parmid, sv) )
	{
	    return sv;
	}

	numvec vv(0);
	if( mo.get_v(parmid, vv) )
	{
            for( unsigned int i=0; i < vv.size(); ++i )
            {
	        os << vv[i] << " ";
            }
	    return os.str();
	}

	any av;
	if( mo.get_a(parmid, av) )
	{
	    throw std::logic_error( "Sorry: any2str not implemented yet" );
	}

	throw std::logic_error( "Error: undefined parameter " );
    }

    vector<int>& get_ivec(molecule_t& m, const hashid_t& cid, const hashid_t& parmid)
    {
        return m.getcmp(cid)->get_ivec(parmid);
    }

    vector<double>& get_dvec(molecule_t& m, const hashid_t& cid, const hashid_t& parmid)
    {
        return m.getcmp(cid)->get_dvec(parmid);
    }

    vector<string>& get_svec(molecule_t& m, const hashid_t& cid, const hashid_t& parmid)
    {
        return m.getcmp(cid)->get_svec(parmid);
    }

    vector<double>& get_vvec(molecule_t& m, const hashid_t& cid, const hashid_t& parmid)
    {
        return m.getcmp(cid)->get_vvec(parmid);
    }

    vector<int> const& get_ivec(const molecule_t& m, const hashid_t& cid, const hashid_t& parmid)
    {
        return m.getcmp(cid)->get_ivec(parmid);
    }

    vector<double> const& get_dvec(const molecule_t& m, const hashid_t& cid, const hashid_t& parmid)
    {
        return m.getcmp(cid)->get_dvec(parmid);
    }

    vector<string> const& get_svec(const molecule_t& m, const hashid_t& cid, const hashid_t& parmid)
    {
        return m.getcmp(cid)->get_svec(parmid);
    }

    vector<double> const& get_vvec(const molecule_t& m, const hashid_t& cid, const hashid_t& parmid)
    {
        return m.getcmp(cid)->get_vvec(parmid);
    }

    double charge(const molecule_t& m)
    {
        const vector<double>& chg = get_dvec( m, ATOM, PCHG );
        double sum = 0;
        for( int i=0; i < (int)chg.size(); ++i )
        {
            sum += chg[i];
        }
        return sum;
    }

    double charge(const morf_t& mo )
    {
        if( mo.cmpid()==ATOM )
        {
            return mo.get_d(PCHG);
        }
        else if( mo.cmpid()==RESD )
        {
            return sum( mo.atom_begin(), mo.atom_end(), 0.0, dparm_getter(PCHG) );
        }
        else
        {
            throw std::runtime_error( "Error: cannot get pchg of a " + unhash(mo.cmpid()) );
        }
    }    

    double get_vdwr(const string& type)
    {
        static molecule_ptr ff;

        try
        {
            ff = mortenv().get_mol("_amberffp");
	}
        catch( std::exception& )
        {
            //std::cout << "Warning: no amber force field has been loaded, VDW radii set to 1.5." << std::endl;
            return 1.5;
        }

        assert( ff !=NULL );

	atom_t pa( *ff, -1 );
	if( atom_t::get(*ff, type, pa) )
	{
            // in parm99.dat, some atom types has zero radius (HO and HW for example), which would be too small
            // change them to 1.0 here.
	    double r = pa.get_d(RSTAR);
            return (r<0.01) ? 1.0 : r;
	}

        // for extra points
	if( type[0]=='E' )
	{
	    return 0.2;
	}

        //std::cout << "Warning: unknown amber atom type " << atype << ", using radii 1.5" << std::endl;
	return 1.5;
    }


    void set_vdwr( const molecule_t& m, vector<double>& vdwr )
    {
        const vector<string>& type = get_svec( m, ATOM, TYPE );
        int natom = m.natom();
        for( int i=0; i < natom; ++i )
        {
            vdwr[i] = get_vdwr(type[i]);
        }
    }


    void set_vdwr( molecule_t& m )
    {
        m.atom_begin()->set_d(VDWR, 0.0);
        
        vector<double>& vdwr = get_dvec( m, ATOM, VDWR );
        
        set_vdwr(m, vdwr);
    }

    void set_vdwr( morf_t& mo )
    {
        if( mo.cmpid()==ATOM )
        {
            mo.set_d( VDWR, get_vdwr(mo.get_s(TYPE)) );
        }
        else if( mo.cmpid()==RESD )
        {
            atomiter_t ai = mo.atom_begin();
            for( ; ai != mo.atom_end(); ++ai )
            {
                set_vdwr( *ai );
            }
        }
        else
        {
            throw std::runtime_error( "Error: cannot set vdwr for a " + unhash(mo.cmpid()) );
        }
    }

} // namespace mort


