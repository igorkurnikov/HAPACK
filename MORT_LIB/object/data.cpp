#include <map>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <common.hpp>
#include "data.hpp"

namespace mort
{
    using namespace boost;
    using std::pair;
    using std::runtime_error;

    template<typename T>
    struct parm_resize
    {
        typedef void result_type;

        void operator()(pair<const hashid_t, vector<T> >& input, int newsize)
        {
            input.second.resize(newsize);
        }
    };

    template<typename T>
    bool get_mholder(map<hashid_t, vector<T> >& content,
            const hashid_t& parmid, vector<T>*& p)
    {
        typedef typename map<hashid_t, vector<T> >::iterator iter_t;
        iter_t i = content.find(parmid);
        if( i==content.end() ) return false;
        p = &(i->second);
        return true;
    }

    template< typename T>
    bool get_cholder( const map<hashid_t, vector<T> >& content, const hashid_t& parmid, const vector<T>*& p )
    {
        typedef typename map<hashid_t, vector<T> >::const_iterator iter_t;
        iter_t i = content.find(parmid);
        if( i==content.end() ) return false;
        p = &(i->second);
        return true;
    }

    template< typename T>
    vector<T>* frcget_holder( map<hashid_t, vector<T> > & content, const hashid_t& parmid, int size )
    {
        vector<T>* p = &( content[parmid] );

        if( (int)(p->size()) != size )
        {
            p->resize( size );
        }

        return p;
    }

    template< typename T >
    void compress_byid_n( vector<T>& src, const vector<int>& ids )
    {
        if( src.size()==ids.size() )
        return;

        assert( src.size() > ids.size() );

        vector<T> dst( ids.size() );

        for(int i=0; i < (int)ids.size(); ++i)
        {
            int id = ids[i];
            assert( id >=0 && id < (int)src.size() );
            dst[i] = src[id];
        }

        src.swap( dst );
    }

    template< typename T >
    void compress_content_n( map< hashid_t, vector<T> >& content, const vector<int>& ids )
    {
        typename map< hashid_t, vector<T> >::iterator i = content.begin();
        for(; i != content.end(); ++i )
        {
            //std::cout << "compressing " << unhash( i->first ) << " " << " src, dst: " << i->second.size() << " " << ids.size() << std::endl;
            compress_byid_n( i->second, ids );
        }
    }

    void compress_byid_v( vector<double>& src, const vector<int>& ids, int length )
    {
        assert( length>0 && src.size()>=length*ids.size() );

        vector< double > dst( length*ids.size() );

        for( int i=0; i < (int)ids.size(); ++i )
        {
            int id = ids[i];
            assert( id >=0 && id < (int)src.size()/length );
            for( int j=0; j < length; ++j )
            {
                dst[length*i+j] = src[length*id+j];
            }
        }

        src.swap( dst );
    }

    void compress_content_v( map< hashid_t, vector<double> >& content, const vector<int>& ids, const map<hashid_t, int>& vlength )
    {
        map<hashid_t, vector<double> >::iterator i = content.begin();
        for(; i != content.end(); ++i )
        {
            assert( vlength.count( i->first ) > 0 );
            compress_byid_v( i->second, ids, vlength.find( i->first )->second );
        }
    }

    void adjust_occupy( set<int>& src, const map<int, int>& old2new )
    {
        set<int> dst;

        set<int>::iterator i = src.begin();
        for(; i != src.end(); ++i )
        {
            int old = *i;
            assert( old2new.count(old) > 0 );
            dst.insert( old2new.find(old)->second );
        }

    }

    void adjust_occupy( map< hashid_t, set<int> >& occupy, const map<int, int>& old2new )
    {
        map< hashid_t, set<int> >::iterator i = occupy.begin();
        for(; i != occupy.end(); ++i )
        {
            adjust_occupy( i->second, old2new );
        }
    }

    mcmpdata_t::mcmpdata_t(int size)
    {
        m_clean = true;
        m_idcounter = 0;
        resize( size );
        assert( m_idcounter == size );
    }

    mcmpdata_t::mcmpdata_t( const mcmpdata_t& rhs )
    :m_icontent(rhs.m_icontent), m_dcontent(rhs.m_dcontent),
    m_scontent(rhs.m_scontent), m_vcontent(rhs.m_vcontent),
    m_acontent(rhs.m_acontent), m_vlength (rhs.m_vlength),
    m_absids(rhs.m_absids), m_idcounter(rhs.m_idcounter)
    {
    }

    mcmpdata_t::~mcmpdata_t()
    {
    }

    vector< int >::iterator mcmpdata_t::append( int n )
    {
        resize( size() + n );
        return m_absids.end()-n;
    }

    vector<int>::iterator mcmpdata_t::insert( int seq, int n )
    {
        vector<int>::iterator oldpos = m_absids.begin()+seq;

        m_absids.insert( oldpos, n, -1 );

        // it is possible that the memory has been reallocated
        // and moved to new address, thus it is necessary to create
        // nwe pos

        vector<int>::iterator newpos = m_absids.begin() + seq;

        // assign new absid
        vector<int>::iterator j = newpos;
        for(int i=0; i < n; ++i,++j)
        {
            assert( *j == -1 );
            *j = m_idcounter;
            ++m_idcounter;
        }

        resize_parm();

        m_clean = false;
        return newpos;
    }

    bool mcmpdata_t::remove( int id )
    {
        vector< int >::iterator i = std::find( m_absids.begin(), m_absids.end(), id );
        if( i == m_absids.end() ) return false;
        m_absids.erase( i );

        return true;
    }

    void mcmpdata_t::resize( int size )
    {
        resize_absids( size );
        resize_parm( );
    }

    int mcmpdata_t::size() const
    {
        return m_absids.size();
    }

    vector<int>::iterator mcmpdata_t::absid_begin() const
    {
        vector<int>* absids = const_cast< vector<int>* >( &m_absids );
        return absids->begin();
    }

    vector<int>::iterator mcmpdata_t::absid_end() const
    {
        vector<int>* absids = const_cast< vector<int>* >( &m_absids );
        return absids->end();
    }

    bool mcmpdata_t::isclean() const
    {
        return m_clean;
    }

    void mcmpdata_t::cleanup(map<int,int>& old2new)
    {
        compress_content_n( m_icontent, m_absids );
        compress_content_n( m_dcontent, m_absids );
        compress_content_n( m_scontent, m_absids );
        compress_content_v( m_vcontent, m_absids, m_vlength );
        compress_content_n( m_acontent, m_absids );

        for(int i=0; i < (int)m_absids.size(); ++i )
        {
            int old = m_absids[i];
            m_absids[i] = i;
            old2new[old] = m_absids[i];
        }
        m_idcounter = m_absids.size();

        adjust_occupy( m_ioccupy, old2new );
        adjust_occupy( m_doccupy, old2new );
        adjust_occupy( m_soccupy, old2new );
        adjust_occupy( m_voccupy, old2new );
        adjust_occupy( m_aoccupy, old2new );

    }

    void mcmpdata_t::set_i(const hashid_t& parmid, int absid, const int& value)
    {
        if( absid < 0 || absid >= m_idcounter )
        {
            throw runtime_error( "Error: trying to set integer parameter " + unhash(parmid) + " to object doesn't exist" );
        }

        vector<int>* p = frcget_holder(m_icontent, parmid, m_idcounter);
        p->at(absid) = value;
    }

    void mcmpdata_t::set_d(const hashid_t& parmid, int absid, const double& value)
    {
        assert( absid >= 0 && absid < m_idcounter );
        vector<double>* p = frcget_holder(m_dcontent, parmid, m_idcounter);
        p->at(absid) = value;
    }

    void mcmpdata_t::set_s( const hashid_t& parmid, int absid, const string& value )
    {
        assert( absid >= 0 && absid < m_idcounter );
        vector<string>* p = frcget_holder(m_scontent, parmid, m_idcounter);
        p->at(absid) = value;
    }

    void mcmpdata_t::set_v( const hashid_t& parmid, int absid, const numvec& value )
    {
        assert( absid >= 0 && absid < m_idcounter );
        vector<double>* p = frcget_holder(m_vcontent, parmid, m_idcounter*value.size() );

        if( p->size() != value.size()*m_idcounter )
        {
            throw std::runtime_error( "Error: you are trying to insert numvec with different size,"
                    "this is not the purpose of set_v, use set_a instead" );
        }

        int vlen = value.size();

        if( vlen==0 )
        {
            throw std::runtime_error( "Error: try to set vparm " + unhash(parmid) + " with a zero length vector." );
        }

        m_vlength[parmid] = vlen;
        for(int i=0; i < vlen; ++i )
        {
            (*p)[vlen*absid+i] = value[i];
        }

    }

    void mcmpdata_t::set_a( const hashid_t& parmid, int absid, const any& value )
    {
        assert( absid >= 0 && absid < m_idcounter );
        vector<any>* p = frcget_holder(m_acontent, parmid, m_idcounter);
        p->at(absid) = value;
    }

    // type I getter function, implementation one:
    bool mcmpdata_t::get_i(const hashid_t& parmid, int absid, int& v) const
    {
        const int* pv = NULL;
        bool result = get_iptr(parmid, absid, pv);
        if( result )
        {
            v = *pv;
        }
        return result;
    }

    bool mcmpdata_t::get_d(const hashid_t& parmid, int absid, double& v) const
    {
        const double* pv = NULL;
        bool result = get_dptr(parmid, absid, pv);
        if( result ) v = *pv;
        return result;
    }

    bool mcmpdata_t::get_s(const hashid_t& parmid, int absid, string& v) const
    {
        const string* pv = NULL;
        bool result = get_sptr(parmid, absid, pv);
        if( result ) v = *pv;
        return result;
    }

    bool mcmpdata_t::get_v(const hashid_t& parmid, int absid, numvec& v) const
    {
        assert( absid >=0 && absid < m_idcounter );
        const double* pv = NULL;
        bool result = get_vptr(parmid, absid, pv);
        if( result )
        {
            map<hashid_t,int>::const_iterator i = m_vlength.find(parmid);
            assert( i != m_vlength.end() );
            v = makevec(i->second, pv);
        }
        return result;
    }

    bool mcmpdata_t::get_a(const hashid_t& parmid, int absid, any& v) const
    {
        if( absid<0 || absid>=m_idcounter )
        {
            throw runtime_error( "Error: while retriving parameter " + unhash(parmid) + ", impossible abs ID");
        }
        const any* pv = NULL;
        bool result = get_aptr(parmid, absid, pv);
        if( result ) v = *pv;
        return result;
    }

    // type I getter function: implementation two
    bool mcmpdata_t::get_iptr(const hashid_t& parmid, int absid, int*& pv)
    {
        assert( absid >=0 && absid < m_idcounter );
        vector<int>* p = NULL;
        if( !get_mholder(m_icontent, parmid, p) ) return false;
        pv = &(p->at(absid));
        return true;
    }

    bool mcmpdata_t::get_dptr(const hashid_t& parmid, int absid, double*& pv)
    {
        assert( absid >=0 && absid < m_idcounter );
        vector<double>* p = NULL;
        if( !get_mholder(m_dcontent, parmid, p) ) return false;
        pv = &(p->at(absid));
        return true;
    }

    bool mcmpdata_t::get_sptr(const hashid_t& parmid, int absid, string*& pv)
    {
        assert( absid >=0 && absid < m_idcounter );
        vector<string>* p = NULL;
        if( !get_mholder(m_scontent, parmid, p) ) return false;
        pv = &(p->at(absid));
        return true;
    }

    bool mcmpdata_t::get_vptr(const hashid_t& parmid, int absid, double*& pv)
    {
        assert( absid >=0 && absid < m_idcounter );
        vector<double>* p = NULL;
        if( !get_mholder(m_vcontent, parmid, p) ) return false;
        map<hashid_t,int>::iterator i = m_vlength.find(parmid);
        assert( i != m_vlength.end() );
        int vlen = i->second;
        pv = &(p->at(vlen*absid));
        return true;
    }

    bool mcmpdata_t::get_aptr(const hashid_t& parmid, int absid, any*& pv)
    {
        assert( absid >=0 && absid < m_idcounter );
        vector<any>* p = NULL;
        if( !get_mholder(m_acontent, parmid, p) ) return false;
        pv = &(p->at(absid));
        return true;
    }

    bool mcmpdata_t::get_iptr(const hashid_t& parmid, int absid, const int*& pv) const
    {
        if( absid<0 || absid>=m_idcounter)
        {
            std::cout << "get " << unhash(parmid) << " of " << absid << std::endl;
        }

        assert( absid >=0 && absid < m_idcounter );
        const vector<int>* p = NULL;
        if( !get_cholder(m_icontent, parmid, p) ) return false;
        pv = &(p->at(absid));
        return true;
    }

    bool mcmpdata_t::get_dptr(const hashid_t& parmid, int absid, const double*& pv) const
    {
        assert( absid >=0 && absid < m_idcounter );
        const vector<double>* p = NULL;
        if( !get_cholder(m_dcontent, parmid, p) ) return false;
        pv = &(p->at(absid));
        return true;
    }

    bool mcmpdata_t::get_sptr(const hashid_t& parmid, int absid, const string*& pv) const
    {
        assert( absid >=0 && absid < m_idcounter );
        const vector<string>* p = NULL;
        if( !get_cholder(m_scontent, parmid, p) ) return false;
        pv = &(p->at(absid));
        return true;
    }

    bool mcmpdata_t::get_vptr(const hashid_t& parmid, int absid, const double*& pv) const
    {
        assert( absid >=0 && absid < m_idcounter );
        const vector<double>* p = NULL;
        if( !get_cholder(m_vcontent, parmid, p) ) return false;
        map<hashid_t,int>::const_iterator i = m_vlength.find(parmid);
        assert( i != m_vlength.end() );
        int vlen = i->second;
        pv = &(p->at(vlen*absid));
        return true;
    }

    bool mcmpdata_t::get_aptr(const hashid_t& parmid, int absid, const any*& pv) const
    {
        assert( absid >=0 && absid < m_idcounter );
        const vector<any>* p = NULL;
        if( !get_cholder(m_acontent, parmid, p) ) return false;
        pv = &(p->at(absid));
        return true;
    }

    // type II getter functions
    int& mcmpdata_t::get_i(const hashid_t& parmid, int absid)
    {
        assert( absid>=0 && absid<m_idcounter);
        int* pv = NULL;
        if( !get_iptr(parmid, absid, pv) )
        {
            throw runtime_error( "Error: undefinded integer parameter " + unhash(parmid) + " of component" );
        }
        return *pv;
    }

    double& mcmpdata_t::get_d(const hashid_t& parmid, int absid)
    {
        double* pv = NULL;
        if( !get_dptr(parmid, absid, pv) )
        {
            throw runtime_error( "Error: undefinded double parameter " + unhash(parmid) );
        }
        return *pv;
    }

    string& mcmpdata_t::get_s(const hashid_t& parmid, int absid)
    {
        string* pv = NULL;
        if( !get_sptr(parmid, absid, pv) )
        {
            throw runtime_error( "Error: undefinded double parameter " + unhash(parmid) );
        }
        return *pv;
    }

    any& mcmpdata_t::get_a(const hashid_t& parmid, int absid)
    {
        any* pv = NULL;
        if( !get_aptr(parmid, absid, pv) )
        throw runtime_error( "Error: undefinded any parameter " + unhash(parmid) );
        return *pv;
    }

    const int& mcmpdata_t::get_i(const hashid_t& parmid, int absid) const
    {
        const int* pv = NULL;
        if( !get_iptr(parmid, absid, pv) )
        throw runtime_error( "Error: undefinded integer parameter " + unhash(parmid) );
        return *pv;
    }

    const double& mcmpdata_t::get_d(const hashid_t& parmid, int absid) const
    {
        const double* pv = NULL;
        if( !get_dptr(parmid, absid, pv) )
        throw runtime_error( "Error: undefinded double parameter " + unhash(parmid) );
        return *pv;
    }

    const string& mcmpdata_t::get_s(const hashid_t& parmid, int absid) const
    {
        const string* pv = NULL;
        if( !get_sptr(parmid, absid, pv) )
        throw runtime_error( "Error: undefinded double parameter " + unhash(parmid) );
        return *pv;
    }

    numvec mcmpdata_t::get_v(const hashid_t& parmid, int absid) const
    {
        numvec value(0);

        if( !get_v(parmid, absid, value) )
        throw runtime_error( "Error: undefinded numvec parameter " + unhash(parmid) );

        return value;
    }

    double* mcmpdata_t::get_dptr(const hashid_t& parmid, int absid)
    {
        double* pd = NULL;
        if( !get_dptr(parmid, absid, pd) )
        throw runtime_error("Error: undefined double parameter " + unhash(parmid) );
        return pd;
    }

    double* mcmpdata_t::get_vptr(const hashid_t& parmid, int absid)
    {
        double* pv = NULL;
        if( !get_vptr(parmid, absid, pv) )
        throw runtime_error( "Error: undefinded numvec parameter " + unhash(parmid) );
        return pv;
    }

    string* mcmpdata_t::get_sptr(const hashid_t& parmid, int absid)
    {
        string* ps = NULL;
        if( !get_sptr(parmid, absid, ps) )
        throw runtime_error( "Error: undefinded numvec parameter " + unhash(parmid) );
        return ps;
    }

    double const* mcmpdata_t::get_dptr(const hashid_t& parmid, int absid) const
    {
        const double* pd = NULL;
        if( !get_dptr(parmid, absid, pd) )
        throw runtime_error("Error: undefined double parameter " + unhash(parmid) );
        return pd;
    }

    const double* mcmpdata_t::get_vptr(const hashid_t& parmid, int absid) const
    {
        const double* pv = NULL;
        if( !get_vptr(parmid, absid, pv) )
        {
            throw runtime_error( "Error: undefinded numvec parameter " + unhash(parmid) );
        }
        return pv;
    }

    string const* mcmpdata_t::get_sptr(const hashid_t& parmid, int absid) const
    {
        const string* ps = NULL;
        if( !get_sptr(parmid, absid, ps) )
        {
            throw runtime_error( "Error: undefinded numvec parameter " + unhash(parmid) );
        }
        return ps;
    }

    const any& mcmpdata_t::get_a(const hashid_t& parmid, int absid) const
    {
        const any* pv = NULL;
        if( !get_aptr(parmid, absid, pv) )
        {
            throw runtime_error( "Error: undefinded any parameter " + unhash(parmid) );
        }
        return *pv;
    }

    // get reference to holder
    vector<int>& mcmpdata_t::get_ivec(const hashid_t& parmid)
    {
        vector<int>* p = NULL;
        if( get_mholder(m_icontent, parmid, p) )
        {
            return *p;
        }
        throw runtime_error( "Error: iparm " + unhash(parmid) + " does not exist!" );
    }

    vector<double>& mcmpdata_t::get_dvec(const hashid_t& parmid)
    {
        vector<double>* p = NULL;
        if( get_mholder(m_dcontent, parmid, p) )
        {
            return *p;
        }
        throw runtime_error( "Error: dparm " + unhash(parmid) + " does not exist!" );
    }

    vector<double>& mcmpdata_t::get_vvec(const hashid_t& parmid)
    {
        vector<double>* p = NULL;
        if( get_mholder(m_vcontent, parmid, p) )
        {
            return *p;
        }
        throw runtime_error( "Error: vparm " + unhash(parmid) + " does not exist!" );
    }

    vector<string>& mcmpdata_t::get_svec(const hashid_t& parmid)
    {
        vector<string>* p = NULL;
        if( get_mholder(m_scontent, parmid, p) )
        {
            return *p;
        }
        throw runtime_error( "Error: sparm " + unhash(parmid) + " does not exist!" );
    }

    vector<int> const& mcmpdata_t::get_ivec(const hashid_t& parmid) const
    {
        const vector<int>* p = NULL;
        if( get_cholder(m_icontent, parmid, p) )
        {
            return *p;
        }
        throw runtime_error( "Error: iparm " + unhash(parmid) + " does not exist!" );
    }

    vector<double> const& mcmpdata_t::get_dvec(const hashid_t& parmid) const
    {
        const vector<double>* p = NULL;
        if( get_cholder(m_dcontent, parmid, p) )
        {
            return *p;
        }
        throw runtime_error( "Error: dparm " + unhash(parmid) + " does not exist!" );
    }

    vector<double> const& mcmpdata_t::get_vvec(const hashid_t& parmid) const
    {
        const vector<double>* p = NULL;
        if( get_cholder(m_vcontent, parmid, p) )
        {
            return *p;
        }
        throw runtime_error( "Error: vparm " + unhash(parmid) + " does not exist!" );
    }

    vector<string> const& mcmpdata_t::get_svec(const hashid_t& parmid) const
    {
        vector<string> const* p = NULL;
        if( get_cholder(m_scontent, parmid, p) )
        {
            return *p;
        }
        throw runtime_error( "Error: sparm " + unhash(parmid) + " does not exist!" );
    }

    // Testing functions
    bool mcmpdata_t::has_i( const hashid_t& parmid, int ) const
    {
        const vector<int>* p = NULL;
        return get_cholder(m_icontent, parmid, p);
    }

    bool mcmpdata_t::has_d( const hashid_t& parmid, int ) const
    {
        const vector<double>* p = NULL;
        return get_cholder(m_dcontent, parmid, p);
    }

    bool mcmpdata_t::has_s( const hashid_t& parmid, int ) const
    {
        const vector<string>* p = NULL;
        return get_cholder(m_scontent, parmid, p);
    }

    bool mcmpdata_t::has_v( const hashid_t& parmid, int ) const
    {
        const vector<double>* p = NULL;
        return get_cholder(m_vcontent, parmid, p);
    }

    bool mcmpdata_t::has_a( const hashid_t& parmid, int ) const
    {
        const vector<any>* p = NULL;
        return get_cholder(m_acontent, parmid, p);
    }

    iparm_iterator mcmpdata_t::iparm_begin() const
    {
        return m_icontent.begin();
    }

    iparm_iterator mcmpdata_t::iparm_end() const
    {
        return m_icontent.end();
    }

    dparm_iterator mcmpdata_t::dparm_begin() const
    {
        return m_dcontent.begin();
    }

    dparm_iterator mcmpdata_t::dparm_end() const
    {
        return m_dcontent.end();
    }

    sparm_iterator mcmpdata_t::sparm_begin() const
    {
        return m_scontent.begin();
    }

    sparm_iterator mcmpdata_t::sparm_end() const
    {
        return m_scontent.end();
    }

    vparm_iterator mcmpdata_t::vparm_begin() const
    {
        return m_vcontent.begin();
    }

    vparm_iterator mcmpdata_t::vparm_end() const
    {
        return m_vcontent.end();
    }

    aparm_iterator mcmpdata_t::aparm_begin() const
    {
        return m_acontent.begin();
    }

    aparm_iterator mcmpdata_t::aparm_end() const
    {
        return m_acontent.end();
    }

    void mcmpdata_t::resize_absids( int size )
    {
        int oldsize = m_absids.size();
        int diff = size - oldsize;

        if( diff > 0 )
        {
            m_absids.resize( size );

            for( int i=oldsize; i < size; ++i )
            {
                m_absids[i] = m_idcounter++;
            }
        }
    }

    void mcmpdata_t::resize_parm( )
    {
        if( m_idcounter > 0 )
        {
            map< hashid_t, vector<int> >::iterator ii = m_icontent.begin();
            map< hashid_t, vector<int> >::iterator ie = m_icontent.end();
            for(; ii != ie; ++ii )
            {
                ii->second.resize( m_idcounter );
            }

            map< hashid_t, vector<double> >::iterator di = m_dcontent.begin();
            map< hashid_t, vector<double> >::iterator de = m_dcontent.end();
            for(; di != de; ++di )
            {
                di->second.resize( m_idcounter );
            }

            map< hashid_t, vector<string> >::iterator si = m_scontent.begin();
            map< hashid_t, vector<string> >::iterator se = m_scontent.end();
            for(; si != se; ++si )
            {
                si->second.resize( m_idcounter );
            }

            map<hashid_t, vector<double> >::iterator vi = m_vcontent.begin();
            map<hashid_t, vector<double> >::iterator ve = m_vcontent.end();
            for(; vi != m_vcontent.end(); ++vi )
            {
                int vlen = m_vlength[ vi->first ];
                assert( vlen != 0 );
                vi->second.resize( vlen*m_idcounter );
            }

            map< hashid_t, vector<any> >::iterator ai = m_acontent.begin();
            map< hashid_t, vector<any> >::iterator ae = m_acontent.end();
            for(; ai != ae; ++ai )
            {
                ai->second.resize( m_idcounter );
            }

        }
    }

} // namespace mort


