#include <object.hpp>
#include "ring.hpp"

namespace mort
{
    using std::set;
    
    ring_finder::ring_finder( const morf_t& o, int type, int size )
        : finder_t( type == WANT_ALL ? ALLOW_MULTI_VISIT : NO_MULTI_VISIT )
    {
        assert( size == 0 || size > 2 );

        m_ringsize = size;

        //assert( size != 0  && type == WANT_ONE );

	if( o.cmpid()==ATOM )
	{
            init( o );
	}
	else
	{
            init( o.atom_begin(), o.atom_end() );
	}
    }

    ring_finder::~ring_finder()
    {
    }

    int ring_finder::predict( const atomvec_t& path)
    {
        if( m_ringsize == 0 ) 
        {
	    return path.is_ring() ? ENOUGH : KEEP_GOING;
        }

        assert( (int)path.size() <= m_ringsize );

        if( (int)path.size() == m_ringsize )
	{
	    return path.is_ring() ? ENOUGH : WRONG_DIRECTION;
	}

	return KEEP_GOING;
    }

    bool has_ring(const morf_t& o)
    {
        atomvec_t ring;
        return find_ring( o, ring );
    }

    bool find_ring( const morf_t& o, atomvec_t& ring )
    {
        ring.clear();
 
        ring_finder f( o, WANT_ONE, 0 );

	return f.find_one( ring );
    }

    bool has_ring(const morf_t& o, int size )
    {
        atomvec_t ring;
        return find_ring(o, size, ring);
    }
    
    bool find_ring(const morf_t& o, int size, atomvec_t& ring )
    {
        if( ! has_ring(o) )
        {
            return false;
        }
        
        ring_finder f(o, WANT_ALL, size);
        
        return f.find_one( ring );
    }

    void find_ring( const morf_t& a, vector<atomvec_t>& rings )
    {
        ring_finder f(a, 0, WANT_ALL);
     
	f.find_all( rings );
  
        return;
    }


    bool find_arom( const morf_t& a, atomvec_t& arom)
    {
        ring_finder f(a, 0, WANT_ALL);
     
        vector<atomvec_t> rings;

	f.find_all( rings );
 
        for( int i=0; i < (int)rings.size(); ++i )
	{
	    if( rings[i].is_arom() )
	    {
	        arom = rings[i];
		return true;
	    }
	}
        return false;
    }


    bool is_arom( const morf_t& atom )
    {
        int arom;
        if( atom.get_i(AROM,arom) )
        {
            return arom;
        }

        atomvec_t arom_ring;
        return find_arom( atom, arom_ring);
    }

} // namespace mort


