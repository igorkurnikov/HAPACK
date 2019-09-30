#include <sstream>
#include <stdexcept>
#include <algorithm>

#include <common.hpp>
#include "morf.hpp"
#include "mole.hpp"

namespace mort
{

    string str( const morf_t& mo )
    {
        std::ostringstream r;
        int natom = mo.natom();
        for( int i=0; i < natom; ++i )
        {
            r << mo.atoms()[i].get_s(NAME);
            if( i < natom-1 )
            {
                r << "-";
            }
        }

        int p;
        if( mo.get_i(PERIOD, p) )
        {
            r << " v" << p;
        }
        return r.str();
    }



    using std::runtime_error;

    void set_atom_elem_from_type(morf_t& a, const string& type ) 
    {
        int element = pertab_t::get_element( type.substr(0,2).c_str() );
        if( element == 0 )
        {
            element = pertab_t::get_element( toupper(type[0]) );
        }

        if( element != 0 )
        {
            a.set_i(ELEMENT, element );
	    assert( a.get_i(ELEMENT) == element );
        }
    }

    void set_bond_order_from_type(morf_t& b, const string& type)
    {
        assert( b.cmpid() == BOND );

        b.getcmp()->set_s(TYPE, b.absid(), type);        

        if( type == "am" || type == "ar" )
        {
            b.set_i(ORDER, 1);
        }
        else if( type.length() == 1 && isdigit( type[0] ) )
        {
            b.set_i(ORDER, type[0] - '0' );
        }
        else if( type == "du" || type == "nc" )
        {
            throw runtime_error( "dummy bond in mol2 file is not supported(yet)!" );
        }
        else
        {
            throw runtime_error( "unknown bond type "+type );
        }
    }

    string get_bond_type_from_order( const morf_t& b )
    {
        assert( b.cmpid() == BOND );
        std::ostringstream os;
        os << b.get_i(ORDER);
        return os.str();
    }

    string get_resd_type( const string& name )
    {
        return name.length()<3 ? name : name.substr(0,3);
    }

    string get_resd_type( const morf_t& r )
    {
        string name;
        if( r.get_s(NAME, name) )
        {
            return get_resd_type( name );
        }

        throw std::runtime_error( "type undefined for residues of molecule " + r.getmol().get_s(NAME) );
    }



    morf_t::morf_t( const molecule_t& mol, int cmpid, int absid )
        :m_molecule( const_cast< molecule_t* >( &mol ) ),
         m_cmpid( cmpid ),
         m_absid( absid )
    {
    }

    morf_t::morf_t( const morf_t& rhs )
        :m_molecule( rhs.m_molecule ),
         m_cmpid( rhs.m_cmpid ),
         m_absid( rhs.m_absid )
    {
    }

    morf_t::~morf_t()
    {
    }

    morf_t& morf_t::operator=( const morf_t& rhs )
    {
        m_molecule = rhs.m_molecule;
        m_absid = rhs.absid();
	return *this;
    }

    int morf_t::absid() const
    {
        return m_absid;
    }

    int morf_t::relid() const
    {
        const mcmpdata_t* pcmp = getcmp();
        vector<int>::iterator ipos = std::find( pcmp->absid_begin(), 
	                                   pcmp->absid_end(),
					   m_absid );
	return (ipos - pcmp->absid_begin());
    }

    hashid_t morf_t::cmpid() const
    {
        return m_cmpid;
    }
 
    molecule_t& morf_t::getmol()
    {
        return *m_molecule;
    }

    molecule_t const& morf_t::getmol() const
    {
        return *m_molecule;
    }

    mcmpdata_t* morf_t::getcmp()
    {
        return m_molecule->getcmp(m_cmpid);
    }

    mcmpdata_t const* morf_t::getcmp() const
    {
        return m_molecule->getcmp(m_cmpid);
    }

    mcmprela_t* morf_t::getadj( int cid)
    {
        return m_molecule->getadj(m_cmpid,cid);
    }

    mcmprela_t const* morf_t::getadj(int cid) const
    {
        return m_molecule->getadj(m_cmpid,cid);
    }

    // setter functions
    void morf_t::set_i(const hashid_t& parmid, const int& value)
    {
        if(parmid == ID)
	{
	    assert( value == absid()+1 );
	    return;
	}

        if(parmid == RESID && getmol().nresd() > 0)
	{
	    int resid = value;
	    morf_t resd( getmol(), RESD, resid-1 );
	    this->connect(resd);
	    resd.connect(*this);
	    return;
	}

	if(parmid == WEIGHT)
	{
	    int wgtdiff = int( value - pertab_t::get_weight(get_i(ELEMENT)) );
	    set_i( WGTDIFF, wgtdiff );
	    return;
	}

        getcmp()->set_i(parmid, m_absid, value);
    }

    void morf_t::set_d( const hashid_t& parmid, const double& value )
    {
        getcmp()->set_d( parmid, m_absid, value );
    }

    void morf_t::set_s( const hashid_t& parmid, const string& value )
    {
        if( parmid == RESNAME )
	{
	    if( nresd() > 0 )
	    {
	        resd_begin()->set_s( NAME, value );
                resd_begin()->set_s( TYPE, get_resd_type(value) );
            }
	    else
	    {
	        getmol().set_s( NAME, value );
	    }
	    return;
	}

	if( parmid==TYPE || parmid==SYMBOL )
	{
	    string type = value;

	    if( parmid == TYPE )
	    {
	        getcmp()->set_s( TYPE, m_absid, type );
	    }

	    if( m_cmpid == ATOM )
	    {
	        set_atom_elem_from_type( *this, type );
            }
            else if( m_cmpid == BOND )
	    {
	        set_bond_order_from_type( *this, type );
            }

	    return;
	}

        getcmp()->set_s( parmid, m_absid, value );
    }

    void morf_t::set_v( const hashid_t& parmid, const numvec& value )
    {
        getcmp()->set_v( parmid, m_absid, value );
    }

    void morf_t::set_a( const hashid_t& parmid, const any& value )
    {  
        typedef std::pair<int, int> ipair_t;
        if( parmid == ATOMPAIR )
	{
	    assert( m_cmpid == BOND );
	    ipair_t ip = any_cast< ipair_t >(value);
	    morf_t atm1( getmol(), ATOM, ip.first-1 );
	    morf_t atm2( getmol(), ATOM, ip.second-1);

            atm1.connect( atm2 );
	    atm2.connect( atm1 );

	    atm1.connect( *this );
	    atm2.connect( *this );

	    this->connect( atm1 );
	    this->connect( atm2 );
	    return;
        }
	    
        getcmp()->set_a( parmid, m_absid, value );
    }

    // functions to set parameter by name
    void morf_t::set_i( const string& parmname, const int& value ) 
    {
        this->set_i( mort::hash(parmname), value );
    }

    void morf_t::set_d( const string& parmname, const double& value )
    {
        this->set_d( mort::hash(parmname), value );
    }

    void morf_t::set_s( const string& parmname, const string& value )
    {
        this->set_s( mort::hash(parmname), value );
    }

    void morf_t::set_v( const string& parmname, const numvec& value )
    {
        this->set_v( mort::hash(parmname), value );
    }

    void morf_t::set_a( const string& parmname, const any& value )
    {
        this->set_a( mort::hash(parmname), value );
    }

    // type I getter function
    bool morf_t::get_i( const hashid_t& parmid, int& value ) const
    {
        if( parmid == ID )
	{
            value = getcmp()->isclean() ? absid()+1 : relid()+1;
            return true;
	}
        
	if( parmid == OFF3 )
	{
	    value = getcmp()->isclean() ? absid()*3 : relid()*3;
	    return true;
	}

	if( parmid == RESID )
	{
	    value = nresd()>0 ? resd_begin()->relid() + 1 : 1;
	    return true;
	}

	if( parmid == WGTDIFF )
	{
	    if( !getcmp()->get_i(WGTDIFF, absid(), value) )
	    {
	        value = 0;
	    }
	    return true;
	}   

        if( parmid == WEIGHT )
	{
	    int weight = int( pertab_t::get_weight(get_i(ELEMENT))+0.5 );
	    value = has_i(WGTDIFF)? weight+get_i(WGTDIFF) : weight;
	    return true;
	}

        return getcmp()->get_i( parmid, absid(), value );
    }

    bool morf_t::get_d( const hashid_t& parmid, double& value ) const
    {
        if( parmid == PCHG )
	{
	    if( ! getcmp()->get_d(PCHG, absid(), value) )
	    {
                value = 0.0;
	    }
            return true;
	}

	if( parmid == VDWR )
	{
            if( ! getcmp()->get_d(VDWR, absid(), value) )
	    {
	        value = pertab_t::get_rvdw( get_i(ELEMENT) );
            }
	    return true;
    	}

        if( parmid == WEIGHT )
	{
	    double weight = pertab_t::get_weight(get_i(ELEMENT));
	    value = has_i(WGTDIFF)? weight+get_i(WGTDIFF) : weight;
	    return true;
	}

        return getcmp()->get_d(parmid, absid(), value);
    }

    bool morf_t::get_s(const hashid_t& parmid, string& v) const
    {
        if( parmid == RESNAME )
	{
	    v = nresd()>0 ?resd_begin()->get_s(NAME):getmol().get_s(NAME);
            return true;
        }

	if( parmid == SYMBOL )
	{
	    assert( m_cmpid==ATOM || m_cmpid==BOND );
 
	    if( m_cmpid == ATOM )
	    {
	        v = pertab_t::get_symbol( get_i(ELEMENT) );
		return true;
	    }
	    else
	    {
	        v = get_bond_type_from_order( *this );
		return true;
            }
	}

	if( parmid == TYPE )
	{
	    if( getcmp()->get_s(TYPE, absid(), v) && v.length() > 0 )
	    {
	        return true;
	    }

	    if( m_cmpid == ATOM )
	    {
	        v = pertab_t::get_symbol( get_i(ELEMENT) );
		return true;
	    }
	    else if( m_cmpid == RESD )
	    {
	        v = get_resd_type(*this);
		return true;
	    }
	    else if( m_cmpid == BOND )
	    {
	        v = get_bond_type_from_order(*this);
		return true;
	    }
	    else
	    {
	        return false;
	    }
        }

	if( parmid == NAME )
	{
	    if( getcmp()->get_s(NAME, absid(),v) )
	    {
	        return true;
            }

	    if( m_cmpid == ATOM )
	    {
	        v = get_s( SYMBOL );
		return true;
	    }
	    else
	    {
	        return false;
	    }
	}

        return getcmp()->get_s( parmid, absid(), v );
    }

    bool morf_t::get_v(const hashid_t& parmid, numvec& v) const
    {
        if( parmid==COLOR )
	{
	    int elem = get_i(ELEMENT);
	    double red   = pertab_t::get_red( elem );
	    double blue  = pertab_t::get_blue( elem );
	    double green = pertab_t::get_green( elem );
            v = makevec(red, green, blue, 1.0);
	    return true;
	}

        return getcmp()->get_v(parmid, absid(), v);
    }

    bool morf_t::get_a(const hashid_t& parmid, any& v) const
    {
        typedef std::pair<int,int> ipair_t;

        if( parmid == ATOMPAIR )
	{
	    assert( m_cmpid == BOND );
            atom_range atms = atoms();
	    v = ipair_t( atms[0].get_i(ID), atms[1].get_i(ID) );
	    return true;
	}

        return getcmp()->get_a(parmid,absid(),v);
    }

    bool morf_t::get_i(const string& parmname, int& v) const
    {
        return get_i( mort::hash(parmname), v );
    }

    bool morf_t::get_d(const string& parmname, double& v) const
    {
        return get_d( mort::hash(parmname), v );
    }

    bool morf_t::get_s(const string& parmname, string& v) const
    {
        return get_s( mort::hash(parmname), v );
    }

    bool morf_t::get_v(const string& parmname, numvec& v) const
    {
        return get_v( mort::hash(parmname), v );
    }

    bool morf_t::get_a(const string& parmname, any& v) const
    {
        return get_a( mort::hash(parmname), v );
    }

    bool morf_t::get_iptr(const hashid_t& parmid, int*& pv)
    {
        return getcmp()->get_iptr(parmid, m_absid, pv);
    }

    bool morf_t::get_dptr(const hashid_t& parmid, double*& pv)
    {
        return getcmp()->get_dptr(parmid, m_absid, pv);
    }

    bool morf_t::get_sptr(const hashid_t& parmid, string*& pv)
    {
        return getcmp()->get_sptr(parmid, m_absid, pv);
    }

    bool morf_t::get_vptr(const hashid_t& parmid, double*& pv)
    {
        return getcmp()->get_vptr(parmid, m_absid, pv);
    }

    bool morf_t::get_aptr(const hashid_t& parmid, any*& pv)
    {
        return getcmp()->get_aptr(parmid, m_absid, pv);
    }

    bool morf_t::get_iptr(const hashid_t& parmid, const int*& pv) const
    {
        return getcmp()->get_iptr(parmid, m_absid, pv);
    }

    bool morf_t::get_dptr(const hashid_t& parmid, const double*& pv) const
    {
        return getcmp()->get_dptr(parmid, m_absid, pv);
    }

    bool morf_t::get_sptr(const hashid_t& parmid, const string*& pv) const
    {
        return getcmp()->get_sptr(parmid, m_absid, pv);
    }

    bool morf_t::get_vptr(const hashid_t& parmid, const double*& pv) const
    {
        return getcmp()->get_vptr(parmid, m_absid, pv);
    }

    bool morf_t::get_aptr(const hashid_t& parmid, const any*& pv) const
    {
        return getcmp()->get_aptr(parmid, m_absid, pv);
    }

    /// type II getter functions
    int morf_t::get_i(const hashid_t& parmid) const
    {
        int v;
        if( !get_i( parmid, v) )
        {
	    throw runtime_error( "Error: undefined integer parameter " + unhash(parmid) + " of " + unhash(m_cmpid) );
        }
        return v;
    }

    double morf_t::get_d(const hashid_t& parmid) const
    {
        double v;
        if( !get_d( parmid, v) ) 
	    throw runtime_error( "Error: undefined double parameter " + unhash(parmid) + " of " + unhash(m_cmpid) );
        return v;
    }

    string morf_t::get_s(const hashid_t& parmid) const
    {
        string v;
        if( !get_s( parmid, v) )
	    throw runtime_error( "Error: undefined string parameter " + v + " " + unhash(parmid) + " of " + unhash(m_cmpid) );
	return v;
    }

    numvec morf_t::get_v(const hashid_t& parmid) const
    {
        numvec v(0);
        if( !get_v( parmid, v) )
     	    throw runtime_error( "Error: undefined numvec parameter " + unhash(parmid) + " of " + unhash(m_cmpid) );
	return v;
    }

    any morf_t::get_a(const hashid_t& parmid) const
    {
        any v;
        if( !get_a( parmid, v) )
	    throw runtime_error( "Error: undefined any other typed parameter " + unhash(parmid) + " of " + unhash(m_cmpid) );
	return v;
    }

    int morf_t::get_i( const string& parmname ) const
    {
        return this->get_i( mort::hash(parmname) );
    }

    double morf_t::get_d( const string& parmname ) const
    {
        return this->get_d( mort::hash(parmname) );
    }

    string morf_t::get_s( const string& parmname ) const
    {
        return this->get_s( mort::hash(parmname) );
    }

    numvec morf_t::get_v( const string& parmname ) const
    {
        return this->get_v( mort::hash(parmname) );
    }

    any morf_t::get_a( const string& parmname ) const
    {
        return this->get_a( mort::hash(parmname) );
    }

    int* morf_t::get_iptr(const hashid_t& parmid)
    {
        return &(getcmp()->get_i(parmid, m_absid));
    }

    double* morf_t::get_dptr(const hashid_t& parmid)
    {
        return &(getcmp()->get_d(parmid, m_absid));
    }

    string* morf_t::get_sptr(const hashid_t& parmid)
    {
        return &(getcmp()->get_s(parmid, m_absid));
    }

    double* morf_t::get_vptr(const hashid_t& parmid)
    {
        return getcmp()->get_vptr(parmid, m_absid);
    }

    any* morf_t::get_aptr(const hashid_t& parmid)
    {
        return &(getcmp()->get_a(parmid, m_absid));
    }

    const int* morf_t::get_iptr(const hashid_t& parmid) const
    {
        return &(getcmp()->get_i(parmid, m_absid));
    }

    const double* morf_t::get_dptr(const hashid_t& parmid) const
    {
        return &(getcmp()->get_d(parmid, m_absid));
    }

    const string* morf_t::get_sptr(const hashid_t& parmid) const
    {
        return &(getcmp()->get_s(parmid, m_absid));
    }

    const double* morf_t::get_vptr(const hashid_t& parmid) const
    {
        return getcmp()->get_vptr(parmid, m_absid);
    }

    const any* morf_t::get_aptr(const hashid_t& parmid) const
    {
        return &(getcmp()->get_a(parmid, m_absid));
    }

    // series of member function to test existance of a paramter.
    //
    bool morf_t::has_i( const hashid_t& parmid ) const
    {
        int v;
        return get_i(parmid, v);
    }

    bool morf_t::has_d( const hashid_t& parmid ) const
    {
        double v;
        return get_d(parmid, v);
    }

    bool morf_t::has_s( const hashid_t& parmid ) const
    {
        string v;
        return get_s(parmid, v);
    }

    bool morf_t::has_v( const hashid_t& parmid ) const
    {
        numvec v(0);
        return get_v(parmid, v);
    }

    bool morf_t::has_a( const hashid_t& parmid ) const
    {
        any v;
        return get_a(parmid, v);
    }

    bool morf_t::has_i( const string& parmname ) const
    {
        return this->has_i( mort::hash(parmname) );
    }

    bool morf_t::has_d( const string& parmname ) const
    {
        return this->has_d( mort::hash(parmname) );
    }

    bool morf_t::has_s( const string& parmname ) const
    {
        return this->has_s( mort::hash(parmname) );
    }

    bool morf_t::has_v( const string& parmname ) const
    {
        return this->has_v( mort::hash(parmname) );
    }

    bool morf_t::has_a( const string& parmname ) const
    {
        return this->has_a( mort::hash(parmname) );
    }

    bool morf_t::has_iptr( const hashid_t& parmid ) const
    {
        return getcmp()->has_i( parmid, absid() );
    }

    bool morf_t::has_dptr( const hashid_t& parmid ) const
    {
        return getcmp()->has_d( parmid, absid() );
    }

    bool morf_t::has_sptr( const hashid_t& parmid ) const
    {
        return getcmp()->has_s( parmid, absid() );
    }

    bool morf_t::has_vptr( const hashid_t& parmid ) const
    {
        return getcmp()->has_v( parmid, absid() );
    }

    bool morf_t::has_aptr( const hashid_t& parmid ) const
    {
        return getcmp()->has_a( parmid, absid() );
    }

    // function for connections
    bool morf_t::connect(const morf_t& p)
    {
        mcmprela_t* adj = getmol().getadj( m_cmpid, p.cmpid() );
	return adj->add( absid(), p.absid() );
    }

    bool morf_t::is_connected_to( const morf_t& p ) const
    {
        const mcmprela_t* adj = getmol().getadj( m_cmpid, p.cmpid() );
        return adj->has( absid(), p.absid() );
    }

    bool morf_t::disconnect(const morf_t& p)
    {
        mcmprela_t* adj = getmol().getadj( m_cmpid, p.cmpid() );
	return adj->remove( absid(), p.absid() );
    }

    // series of member function to get related objects
    //
    mobjiter_t morf_t::mobj_begin( const hashid_t& id) const
    {
        const mcmprela_t* a = getadj(id);
        return mobjiter_t(morf_t(getmol(),id,-1),a->begin(absid()));
    }

    mobjiter_t morf_t::mobj_end  ( const hashid_t& id) const
    {
        const mcmprela_t* a = getadj(id);
	return mobjiter_t(morf_t(getmol(),id,-1),a->end  (absid()));
    }

    atomiter_t morf_t::atom_begin() const { return mobj_begin(ATOM); }
    atomiter_t morf_t::atom_end  () const { return mobj_end  (ATOM); }
    bonditer_t morf_t::bond_begin() const { return mobj_begin(BOND); }
    bonditer_t morf_t::bond_end  () const { return mobj_end  (BOND); }
    resditer_t morf_t::resd_begin() const { return mobj_begin(RESD); }
    resditer_t morf_t::resd_end  () const { return mobj_end  (RESD); }
    angliter_t morf_t::angl_begin() const { return mobj_begin(ANGL); }
    angliter_t morf_t::angl_end  () const { return mobj_end  (ANGL); }
    diheiter_t morf_t::dihe_begin() const { return mobj_begin(TORS); }
    diheiter_t morf_t::dihe_end  () const { return mobj_end  (TORS); }
    impriter_t morf_t::impr_begin() const { return mobj_begin(OOPS); }
    impriter_t morf_t::impr_end  () const { return mobj_end  (OOPS); }
    mobjiter_t morf_t::tor2_begin() const { return mobj_begin(TOR2); }
    mobjiter_t morf_t::tor2_end  () const { return mobj_end  (TOR2); }
    mobjiter_t morf_t::ptor_begin() const { return mobj_begin(PTOR); }
    mobjiter_t morf_t::ptor_end  () const { return mobj_end  (PTOR); }

    vector<int> const& morf_t::related_atom_ids() const
    {
        const mcmprela_t* a = getadj(ATOM);
        return a->getvec( absid() );
    }

    vector<int> const& morf_t::related_bond_ids() const
    {
        const mcmprela_t* a = getadj(BOND);
        return a->getvec( absid() );
    }

    vector<int> const& morf_t::related_resd_ids() const
    {
        const mcmprela_t* a = getadj(RESD);
        return a->getvec( absid() );
    }

    mobj_range morf_t::mobjs(const hashid_t& id) const { return mobj_range(mobj_begin(id),mobj_end(id)); }
    atom_range morf_t::atoms() const { return mobjs(ATOM); }
    bond_range morf_t::bonds() const { return mobjs(BOND); }
    angl_range morf_t::angls() const { return mobjs(ANGL); }
    dihe_range morf_t::dihes() const { return mobjs(TORS); }
    resd_range morf_t::resds() const { return mobjs(RESD); }
    mobj_range morf_t::tor2s() const { return mobjs(TOR2); }
    mobj_range morf_t::ptors() const { return mobjs(PTOR); }

    int morf_t::nmobj(const hashid_t& id) const 
    { 
        return getadj(id)->size( absid() ); 
    }
 
    void validate( const morf_t& lhs, const morf_t& rhs )
    {
        if( &lhs.getmol() != &rhs.getmol() )
        {
            throw runtime_error( "Error: compare morf_t of different molecule is meaningless" );
        }
        
        if( lhs.cmpid() != rhs.cmpid() )
        {
            throw runtime_error( "Error: compare morf_t of different component is meaningless" );
        }
    }

    bool operator==( const morf_t& lhs, const morf_t& rhs )
    {
        validate( lhs, rhs );
        return lhs.absid() == rhs.absid();
    }

    bool operator!=( const morf_t& lhs, const morf_t& rhs )
    {
        validate( lhs, rhs );
        return lhs.absid() != rhs.absid();
    }

    bool operator< ( const morf_t& lhs, const morf_t& rhs )
    {
        validate( lhs, rhs );
        return lhs.absid() < rhs.absid();
    }

    bool operator> ( const morf_t& lhs, const morf_t& rhs )
    {
        validate( lhs, rhs );
        return lhs.absid() > rhs.absid();
    }

    bool operator<=( const morf_t& lhs, const morf_t& rhs )
    {
        validate( lhs, rhs );
        return lhs.absid() <= rhs.absid();
    }

    bool operator>=( const morf_t& lhs, const morf_t& rhs )
    {
        validate( lhs, rhs );
        return lhs.absid() >= rhs.absid();
    }


} // namespace mort


