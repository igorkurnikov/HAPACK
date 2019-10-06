#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <common.hpp>
#include "rela.hpp"

namespace mort
{

    mcmprela_t::mcmprela_t( )
    {
    }

    mcmprela_t::mcmprela_t( mcmprela_t const& rhs )
        :m_content( rhs.m_content )
    {
    }

    mcmprela_t::~mcmprela_t()
    {
    }

    void mcmprela_t::swap( mcmprela_t& rhs )
    {
        m_content.swap( rhs.m_content );
    }

    bool mcmprela_t::add( int a, int b )
    {
        assert( a >= 0 && b >= 0 );

        if( has(a, b) )
        {
            getvec(a).push_back(b);
            return false;
        }

        getvec(a).push_back(b);
        return true;
    }


    bool mcmprela_t::remove( int a, int b )
    {
        iterator i = find(a,b);

	if( i == end(a) )
	{
	    return false;
	}

        getvec(a).erase(i);
	return true;
    }

    void mcmprela_t::clear( int a )
    {
        getvec(a).clear();
    }

    bool mcmprela_t::has( int a, int b ) const
    {
        return find(a, b) != end( a );
    }

    mcmprela_t::iterator mcmprela_t::find( int a, int b ) const
    {
        return std::find( begin(a), end(a), b );
    }

    mcmprela_t::iterator mcmprela_t::begin( int a ) const
    {
        return const_cast< vector<int>& >( getvec(a) ).begin();
    }

    mcmprela_t::iterator mcmprela_t::end( int a ) const
    {
        return const_cast< vector<int>& >( getvec(a) ).end();
    }

    int mcmprela_t::size( int a ) const
    {
        return getvec(a).size();
    }

    bool mcmprela_t::empty() const
    {
        for( int i=0; i < (int)m_content.size(); ++i )
        {
            if( m_content[i].size() > 0 )
                return false;
        }

        return true;
    }


    vector<int>& mcmprela_t::getvec(int a)
    {
        assert( a >= 0 );

        if( a >= (int)m_content.size() )
        {
            m_content.resize( a + 1 );
        }

        return m_content[a];
    }

    vector<int> const& mcmprela_t::getvec( int a ) const
    {
        if( a < 0 )
        {
            throw std::runtime_error( "Error: in getvec, negative absid" );
        }

        assert( a >= 0 );
        if( a >= (int)m_content.size() )
        {
            static vector<int> null;
            return null;
        }

        return m_content[a];
    }

} // namespace mort

