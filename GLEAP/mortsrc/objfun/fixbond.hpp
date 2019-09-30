#ifndef MORTSRC_OBJECT_FIXBO_HPP
#define MORTSRC_OBJECT_FIXBO_HPP

#include <set>

namespace mort
{
    using std::set;

    class morf_t;

    class molecule_t;

    class isomer_set
    {
    public:

        isomer_set( const vector<morf_t>& bnds ) : m_bnds( bnds ) {}

	void add_solution( const map<int, int>& solu )
	{
	    vector<int> tmpsolu( m_bnds.size()==0 );
            for( int i=0; i < (int)m_bnds.size(); ++i )
	    {
	        int bid = m_bnds[i].absid();
	        assert( solu.count(bid) > 0 );
		tmpsolu.push_back( solu.find(bid)->second );
	    }

	    for( int i=0; i < (int)m_solus.size(); ++i )
	    {
	        bool equal_i=true;
	        assert( m_solus[i].size()==tmpsolu.size() );
	        for( int j=0; j < (int)m_solus[i].size(); ++j )
		{
		    if( m_solus[i][j] != tmpsolu[j] )
		    {
		        equal_i=false;
			break;
		    }
		}

		if(equal_i)
		    return;
	    }

	    m_solus.push_back( tmpsolu );
	}

	void apply_best( )
	{
	    assert( m_solus.size() > 0 && m_bnds.size()==m_solus[0].size() );
	    for( int i=0; i < (int)m_bnds.size(); ++i )
	    {
                m_bnds[i].set_i( ORDER, m_solus[0][i] );
	    }
	}

	void apply( int j )
	{
	    int id = m_actives.size()>0 ? m_actives[j] : j;

	    assert( m_bnds.size()==m_solus[id].size() );

	    for( int i=0; i < (int)m_bnds.size(); ++i )
	    {
	        m_bnds[i].set_i( ORDER, m_solus[id][i] );
	    }
	}

	void choose( const vector<int>& list )
	{
	    m_actives.clear();
	    m_actives.insert( m_actives.begin(), list.begin(), list.end() );
	}
	       
	int size() const
	{
	    return m_actives.size()>0 ? m_actives.size() : m_solus.size();
	}

        const vector<morf_t>& bnds() const
	{
	    return m_bnds;
	}

	const vector<int>& operator[]( int id ) const
	{
	    return m_solus.at( id );
	}

	bool has_bnds( const vector<morf_t>& bnds ) const
	{
	    assert( bnds.size() > 0 );

	    if( bnds.size() > m_bnds.size() )
	        return false;

            if( std::find(m_bnds.begin(), m_bnds.end(), bnds[0])==m_bnds.end() )
	        return false;

	    for( int i=1; i < (int)bnds.size(); ++i )
	    {
	        if( std::find(m_bnds.begin(), m_bnds.end(), bnds[i])==m_bnds.end() )
		    throw std::runtime_error( "Error: part of bonds in isomer, part not" );
            }

	    return true;
	}

    private:

        vector< int > m_actives;
        vector< morf_t > m_bnds;
	vector< vector<int> > m_solus;
    };

    void fixbo_hard( morf_t& a, set<int>& fixed );

    void fixbo_bylen( morf_t& a, set<int>& fixed );

    void fixbo_conj( molecule_t& m, set<int>& fixed, vector< isomer_set >& isomers );

    void fixbo_hard_nnbr2( morf_t& a, set<int>& fixed );

    void fixbo_hard_nnbr3( morf_t& a, set<int>& fixed );

    void fixbo_acid( morf_t& a, set<int>& fixed );

    void fixbo_isolated( morf_t& b, set<int>& fixed );

    void fixbo_isomer( vector<morf_t>& atms, vector<morf_t>& bnds, set<int>& fixed, isomer_set& isomer );

    void setbo_allone( morf_t& a, set<int>& fixed );

    void divide_bonds( const vector<morf_t>& bnds, vector< vector<morf_t> >& atmgrps, vector< vector<morf_t> >& bndgrps );
 
    bool fixable( int elem );

    bool is_halogen( int elem );

    bool is_acid( morf_t& a );

    bool in_line( const morf_t& a );
 
    bool on_plane( const morf_t& a );

    bool on_plane( const vector<morf_t>& a );

    int nfixed_nbr( const morf_t& a, const set<int>& fixed );

    int ndouble_bnd( const morf_t& a, const set<int>& fixed );

    int order_from_length( const morf_t& bnd );

    numvec plane( const morf_t& a0, const morf_t& a1, const morf_t& a2 );

    void fixbo( molecule_t& m, vector< isomer_set >& isomers );
}

#endif
