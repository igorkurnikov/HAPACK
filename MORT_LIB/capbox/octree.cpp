#include <set>
#include <object.hpp>
#include "shaper.hpp"
#include "octree.hpp"

#ifndef M_PI
#include <math.h>
#define M_PI 3.1415926535897932385
#endif

namespace mort
{
    
    namespace capbox
    {
        int intpow( int sub, int p )
	{
	    int r = 1;
	    for( int i=0; i < p; ++i )
	        r *= sub;
	    return r;
	}
        
        octree_t::octree_t( const numvec& corner, int ndepth, double resolution )
	    : m_corner(corner), m_center(corner), 
              m_ndepth(ndepth), m_resolution(resolution)
        {
            m_ngrid = intpow( 2, m_ndepth );
            m_length = m_ngrid * resolution;            
            m_radius = m_length * sqrt( 0.75 );
            m_center = m_corner + makevec(m_length/2.0, m_length/2.0, m_length/2.0);
            m_typeid = UNKNOWN_TYPE;
        }
        
        octree_t::~octree_t()
        {
        }
        
        int octree_t::make_shape( shaper_i& shapefun, double ion_radius )
        {
            // possible change for 2nd cut:
            //             excluded -> excluded, 
            //             included -> included or partial
            //             partial  -> partial  or excluded
            if( m_typeid==EXCLUDED )
            {
                // must be excluded
                return EXCLUDED;
            }
            
            // if m_ndepth is 0, the node should be treated as a point.
            // the result of shapefun can only be INCLUDED or EXCLUDED.
            // PARTIAL is impossible.
            if( m_ndepth==0 )
            {
                assert( m_typeid == UNKNOWN_TYPE || m_typeid==INCLUDED );
                m_typeid = shapefun.check_atm( m_corner, ion_radius );
                assert( m_typeid==INCLUDED || m_typeid==EXCLUDED );
                return m_typeid;
            }

            // if m_ndepth > 0, the node is treat as a sphere (with a radius),
            // thus the result of shapefun could be INCLUDED, EXCLUDED or partial
            int type = UNKNOWN_TYPE;
            if( m_typeid==UNKNOWN_TYPE || m_typeid==INCLUDED )
            {
                type = shapefun.check_box( m_center, makevec(m_length, m_length, m_length), ion_radius );

                if( type==EXCLUDED )
                {
                    // not in the shape, does not matter if 1st or 2nd, ruled out, do not go furthur
                    m_typeid = EXCLUDED;
                    return m_typeid;
                }
            
                if( type==INCLUDED )
                {
                    m_typeid = type;
                    return m_typeid;
                }

                assert( type==PARTIAL );

                if( m_typeid==UNKNOWN_TYPE )
                {
                    m_typeid = PARTIAL;
                }
            }

            assert( m_typeid==PARTIAL || m_typeid==INCLUDED );


            // whenever partial exist, we need to go to the next level.
	    if( m_subs[0]==NULL ) split();

            int nincluded = 0;
            int nexcluded = 0;
            
            for( int i=0; i < 8; ++i )
            {
                int subtype = m_subs[i]->make_shape( shapefun, ion_radius );
                if( subtype==INCLUDED )
                {
                    nincluded++;
                }
                
                if( subtype==EXCLUDED )
                {
                    nexcluded++;
                }
            }
            
            // each node is a cubic box, which is smaller then the sphere, thus
            // it is possible that the box is not really partial, change the result
            // accordingly. 
            if( nexcluded==8 )
            {
                m_typeid = EXCLUDED;
                return m_typeid;
            }
            
            if( nincluded==8 )
            {
                m_typeid = INCLUDED;
                return m_typeid;
            }
            
            // now we are sure the node must be partial
            //
            assert( m_typeid == PARTIAL || m_typeid==INCLUDED );
            m_typeid = PARTIAL;
            return m_typeid;
        }
       
	void octree_t::split()
	{
	    assert( m_typeid == INCLUDED || m_typeid==PARTIAL );
	    assert( m_ndepth >= 1 && m_ngrid >= 2 );

            for( int i=0; i < 8; ++i )
            {
                int ix = i / 4;
                int iy = ( i % 4 ) / 2;
                int iz = i % 2;

                numvec corner(3);
                corner[0] = m_corner[0] + ix * m_length / 2;
                corner[1] = m_corner[1] + iy * m_length / 2;
                corner[2] = m_corner[2] + iz * m_length / 2;
		assert( m_subs[i] == NULL );
                m_subs[i] = octree_ptr( new octree_t( corner, m_ndepth-1, m_resolution ) );

                // if parent is included (means it is the second shaper), children are included, 
		// otherwise stay as unknown_type (the init value)
		if( m_typeid==INCLUDED )
		{
                    m_subs[i]->m_typeid = INCLUDED;
                }

                if( m_subpol.size() == 0 )
		{
		    continue;
		}


                // split the value of electro potential
                int xstart = ix*m_ngrid/2;
		int ystart = iy*m_ngrid/2;
		int zstart = iz*m_ngrid/2;

                int subgrid = m_ngrid/2;
	        m_subs[i]->m_subpol.resize( subgrid*subgrid*subgrid );
		for( int subx=0; subx<subgrid; ++subx )
		{
		    for( int suby=0; suby<subgrid; ++suby )
		    {
		        for( int subz=0; subz<subgrid; ++subz )
			{
			    int icrd1 = subx*subgrid*subgrid + suby*subgrid + subz;
			    int icrd2 = xstart*m_ngrid*m_ngrid + ystart*m_ngrid + zstart;
			    m_subs[i]->m_subpol[icrd1] = m_subpol[icrd2];
			}
		    }
		}
            }
	}
                
        void octree_t::calculate( potfun_t potfun, std::pair< double, numvec >& min, std::pair< double, numvec >& max )
        {
	    if( m_typeid==EXCLUDED )
	    {
	        return;
	    }
 
            if( m_typeid==PARTIAL )
	    {
                for( int i=0; i < 8; i++ )
                {
                    if( m_subs[i] != NULL )
                    {
                        m_subs[i]->calculate( potfun, min, max );
                    }
                }

		return;
	    }
            
	    assert( m_typeid==INCLUDED );

            if( m_subpol.size()==0 ) m_subpol.resize( m_ngrid*m_ngrid*m_ngrid, 0.0 );

            double* psubpol = &m_subpol[0];
	    for( int i=0; i < m_ngrid; ++i)
	    {
	        double x = m_corner[0] + i*m_resolution;
	        for( int j=0; j < m_ngrid; ++j )
		{
		    double y = m_corner[1] + j*m_resolution;
		    for( int k=0; k < m_ngrid; ++k )
		    {
		        double z = m_corner[2] + k*m_resolution;
                        *psubpol = *psubpol + potfun( x, y, z );
                        if( *psubpol < min.first )
			{
			    min.first = *psubpol;
			    min.second[0] = x;
			    min.second[1] = y;
			    min.second[2] = z;
			}

                        if( max.first < *psubpol )
                        {
                            max.first = *psubpol;
                            max.second[0] = x;
			    max.second[1] = y;
			    max.second[2] = z;
			}

			psubpol++;
                    }
                }
            }


        }
        
    } // namespace capbox

} // namespace mort



