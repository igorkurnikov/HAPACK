#include <iostream>
#include <stdexcept>
#include <common.hpp>
#include <object.hpp>
#include <objfun.hpp>

namespace mort
{
    static const long long BASE4 = 40*40*40*40;
    
    using std::runtime_error;

    using namespace boost;

    molecule_t::molecule_t()
    {
    }
    
    molecule_t::molecule_t( const molecule_t& rhs )
        : root_t( rhs ), 
	  m_components( rhs.m_components), 
	  m_adjacencys( rhs.m_adjacencys),
	  m_connecteds( rhs.m_connecteds)
    {
    }
    
    molecule_t& molecule_t::operator=( const molecule_t& rhs )
    {
        if( &rhs != this )
        {
            molecule_t temp( rhs );
            swap( temp );
        }
        
        return *this;
    }

    molecule_t::~molecule_t()
    {
    }
    
    void molecule_t::swap( molecule_t& rhs )
    {
        if( &rhs != this )
        {
            root_t::swap( rhs );
            m_components.swap( rhs.m_components );
            m_adjacencys.swap( rhs.m_adjacencys );
            m_connecteds.swap( rhs.m_connecteds );
        }
    }  
    
    ///////////////////////////////////////////////////////////////
    // accessing component
    //////////////////////////////////////////////////////////////
    mcmpdata_t const& blank_component( )
    {
        static const mcmpdata_t inst;
        return inst;
    }
    
    mcmpdata_t* molecule_t::getcmp( int cmpid )
    {
        assert( cmpid < BASE4 );
        mcmpdata_t* p = &(m_components[cmpid]);
        return p;
    }

    mcmpdata_t const* molecule_t::getcmp( int cid ) const
    {
        assert( cid < BASE4 );
        
        typedef map< hashid_t, mcmpdata_t>::const_iterator com_iterator;
        com_iterator i = m_components.find( cid );
        if( i == m_components.end() )
        {
	    return &blank_component();
        }
        return &(i->second);
    }

    /////////////////////////////////////////////////////////////////////////////
    // accessing adjacency
    ////////////////////////////////////////////////////////////////////////////
    const mcmprela_t& blank_adjacency()
    {
        static const mcmprela_t inst;
        return inst;
    }

    mcmprela_t* molecule_t::getadj(int cid1, int cid2)
    {
        assert( cid1 < BASE4 && cid2 < BASE4 );
        hashid_t aid = cid1 + BASE4*cid2;
        mcmprela_t* p = &m_adjacencys[aid];
	m_connecteds[cid1].insert(cid2);
        return p;
    }

    mcmprela_t const* molecule_t::getadj(int cid1, int cid2) const 
    {
        typedef map< hashid_t, mcmprela_t >::const_iterator adj_iter;
        assert( cid1 < BASE4 && cid2 < BASE4 );
        hashid_t aid = cid1 + BASE4*cid2;
        adj_iter i = m_adjacencys.find(aid);

        if( i == m_adjacencys.end() ) 
	{ 
	    return &blank_adjacency();
        }
        
        return &(i->second);
    }

    //////////////////////////////////////////////////////////////////////
    // set connection list
    //////////////////////////////////////////////////////////////////////
    bool molecule_t::get_nbrlist( int cid1, const set<int>*& p) const
    {
        map< hashid_t, set<int> >::const_iterator i = m_connecteds.find(cid1);
        
        if( i == m_connecteds.end() )
        {
            return false;
        }
        
        p = &(i->second);
        return true;
    }
    
    morf_t molecule_t::create(const hashid_t& cid)
    {
        vector<int>::iterator i = getcmp(cid)->append(1);
        return morf_t(*this, cid, *i );
    }

    morf_t molecule_t::insert(const hashid_t& cid, int lid )
    {
        assert( lid>=0 && lid <= getcmp(cid)->size() );

        vector<int>::iterator i = getcmp(cid)->insert(lid);
        assert( i - getcmp(cid)->absid_begin()==lid );

        return morf_t(*this, cid, *i);
    }

    void molecule_t::remove(const morf_t& mo)
    {
        int cmpid = mo.cmpid();
        int absid = mo.absid();

        getcmp(cmpid)->remove(absid);
        
        const set<int>* pnbr = NULL;
        assert( get_nbrlist( cmpid, pnbr) );

        set<int>::const_iterator i = pnbr->begin();
        for( ; i != pnbr->end(); ++i )
        {
            mcmprela_t* adj1 = getadj(cmpid, *i);
            mcmprela_t* adj2 = getadj(*i, cmpid);

            vector<int>& pvec = adj1->getvec( absid );
            vector<int>::iterator j = pvec.begin();
            for( ; j != pvec.end(); ++j )
            {
                adj2->remove( *j, absid);
            }
            
            adj1->clear(absid);
        }
    }
    
    void molecule_t::remove_atom(const morf_t& a)
    {
        assert( a.cmpid() == ATOM );
        
        vector<morf_t> rbonds( a.bond_begin(), a.bond_end() );
        for( int i=0; i < (int)rbonds.size(); ++i )
        {
            remove( rbonds[i] );
        }

        vector<morf_t> rangls( a.angl_begin(), a.angl_end() );
        for( int i=0; i < (int)rangls.size(); ++i )
        {
            remove( rangls[i] );
        }

        vector<morf_t> rtorss( a.dihe_begin(), a.dihe_end() );
        for( int i=0; i < (int)rtorss.size(); ++i )
        {
            remove( rtorss[i] );
        }
        
        vector<morf_t> roopss( a.impr_begin(), a.impr_end() );
        for( int i=0; i < (int)roopss.size(); ++i )
        {
            remove( roopss[i] );
        }

        vector<morf_t> rptors( a.ptor_begin(), a.ptor_end() );
        for( int i=0; i < (int)rptors.size(); ++i )
        {
            remove( rptors[i] );
        }

        remove(a);
    }

    void molecule_t::remove_bond(const morf_t& b)
    {
        mcmprela_t* b2a = getadj( BOND, ATOM );
	vector<int>& pas = b2a->getvec( b.absid() );
	assert( pas.size()==2 );
        int atom_1 = pas[0];
	int atom_2 = pas[1];


	mcmprela_t* a2a = getadj( ATOM, ATOM );
	a2a->remove( atom_1, atom_2 );
	a2a->remove( atom_2, atom_1 );
	remove( b );
    }
    
    void molecule_t::remove_resd(const morf_t& r)
    {
        assert( r.cmpid() == RESD );
        
        vector<morf_t> ratoms( r.atom_begin(), r.atom_end() );
        for( int i=0; i < (int)ratoms.size(); ++i )
        {
            remove_atom( ratoms[i] );
        }
        
        remove(r);
    }
    
    void molecule_t::clear()
    {
        m_components.clear();
	m_adjacencys.clear();
	m_connecteds.clear();
    }

    void print_vec( std::ostream& os, const vector<int>& v )
    {
        for( int i=0; i < (int)v.size(); ++i )
	{
	    os << v[i] << " ";
	}
    }

    void molecule_t::cleanup( )
    {
        map< hashid_t, map<int,int> > idicts;


        // cleanup each component. meanwhile make a map between old ID and new ID,
	// each component will have such a map, thus we use a map of map here.
        map< hashid_t, mcmpdata_t >::iterator ci = m_components.begin();
        for( ; ci != m_components.end(); ++ci )
        {
            ci->second.cleanup( idicts[ci->first] );
        }
        
        map< hashid_t, mcmprela_t >::iterator ai = m_adjacencys.begin();
        for( ; ai != m_adjacencys.end(); ++ai )
        {
            if( ai->second.empty() )
                continue;


            hashid_t cid_1 = ai->first % BASE4;
            hashid_t cid_2 = ai->first / BASE4;

            assert( idicts.count(cid_1) > 0 && idicts.count(cid_2) > 0 );
            
	    // idict_1 is the ID map of the component 1. idict_2 is the ID map of the component 2.
            map<int, int>& idict_1 = idicts[cid_1];
            map<int, int>& idict_2 = idicts[cid_2];

            mcmprela_t anew;

            // for each adjacency (ai->second), create a new adjacency, at the end use the
	    // new old replace the old one.
            int nobj = ai->second.size();
            for( int i=0; i < nobj; ++i )
            {
                vector<int>& pvec = ai->second.getvec( i );

                // if idmap do not a ID, this mo must already be deleted, its relation should already be deleted.
                if( idict_1.count(i) == 0 )
                {
                    assert( pvec.size() == 0 );
                    continue;
                }
                
                int newid_1 = idict_1[i];
                
                vector<int>::iterator dj = pvec.begin();
                for( ; dj != pvec.end(); ++dj )
                {
                    int oldid_2 = *dj;
                    assert( idict_2.count(oldid_2) > 0 );
                    int newid_2 = idict_2[oldid_2];
                    anew.add(newid_1, newid_2);
                }
            }
            
            ai->second.swap( anew );
        }
        
    }

} // namespace mort


