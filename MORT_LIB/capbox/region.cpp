#include "region.hpp"


namespace mort
{

    // For cap, box, and oct region, currently rad is ignored
    int box_region::check_sphere( const numvec& cnt, double ) const
    {
        assert( cnt.size()==3u );
        assert( m_min.size()==3u && m_max.size()==3u );
        if( cnt[0] < m_min[0] || cnt[0] > m_max[0] ) return OUTSIDE;
        if( cnt[1] < m_min[1] || cnt[1] > m_max[1] ) return OUTSIDE;
        if( cnt[2] < m_min[2] || cnt[2] > m_max[2] ) return OUTSIDE;
        return INCLUDE;
    }

    int box_region::check_rndbox( const numvec& cnt, const numvec& ext, double rad ) const
    {
        assert( cnt.size()==3u && ext.size()==3u );
        numvec corner1 = cnt - 0.5*ext;
        numvec corner2 = cnt + 0.5*ext;
        int r1 = check_sphere( corner1, rad );
        int r2 = check_sphere( corner2, rad );
        return (r1==r2) ? r1 : PARTIAL;
    }


    int cap_region::check_sphere( const numvec& cnt, double ) const
    {
        double dx = cnt[0] - m_cnt[0];
        double dy = cnt[1] - m_cnt[1];
        double dz = cnt[2] - m_cnt[2];
        return (dx*dx+dy*dy+dz*dz < m_rad*m_rad) ? INCLUDE : OUTSIDE;
    }

    int cap_region::check_rndbox( const numvec& cnt, const numvec& extent, double ) const
    {
        double boxrad = norm(extent)*0.5;
        double dis = dist( cnt,m_cnt );
        if( dis + boxrad < m_rad ) return INCLUDE;
        if( dis < boxrad + m_rad ) return PARTIAL;
        return OUTSIDE;
    }

    int oct_region::check_sphere( const numvec& cnt, double ) const
    {
        double dx = cnt[0] - m_cnt[0];
        double dy = cnt[1] - m_cnt[1];
        double dz = cnt[2] - m_cnt[2];

        if( std::abs(dx) > 0.5*m_ext[0] ) return OUTSIDE;
        if( std::abs(dy) > 0.5*m_ext[1] ) return OUTSIDE;
        if( std::abs(dz) > 0.5*m_ext[2] ) return OUTSIDE;

        double fx = std::abs( dx/m_ext[0] );
        double fy = std::abs( dy/m_ext[1] );
        double fz = std::abs( dz/m_ext[2] );
        return (fx+fy+fz < 0.75) ? INCLUDE : OUTSIDE;
    }

    int oct_region::check_rndbox( const numvec& cnt, const numvec& extent, double rad ) const
    {
        int result = UNKNOWN;
        for( int ic=0; ic < 8; ++ic )
        {
            int ix = ic%2;
            int iy = (ic/2)%2;
            int iz = ic/4;

            numvec corner(3);
            corner[0] = cnt[0] + (ix-0.5)*extent[0];
            corner[1] = cnt[1] + (iy-0.5)*extent[1];
            corner[2] = cnt[2] + (iz-0.5)*extent[2];

            if( ic==0 )
            {   
                result = check_sphere( corner, rad );
                continue;
            }

            if( result != check_sphere(corner, rad) )
            {
                return PARTIAL;
            }
        }

        return result;
    }

    ////////////////////////////////////////////////////////////////////////
    // 
    // implementation of out_solute
    //
    ////////////////////////////////////////////////////////////////////////

    int out_solute::check_sphere( const numvec& sphcnt, double sphrad ) const
    {
        double mind2 = DOUBLEMAX;
        return check_sphere(sphcnt, sphrad, mind2);
    }

    int out_solute::check_rndbox( const numvec& boxcnt, const numvec& boxext, double sphrad ) const
    {
        double mind2 = DOUBLEMAX;
        return check_rndbox(boxcnt, boxext, sphrad, mind2);
    }


    int out_solute::check_sphere( const numvec& sphcnt, double sphrad, double& mind2 ) const
    {
        for( int i=0; i < m_nsph; ++i )
        {
            double dx = sphcnt[0] - m_cnts[3*i];
            double dy = sphcnt[1] - m_cnts[3*i+1];
            double dz = sphcnt[2] - m_cnts[3*i+2];
            double d2 = dx*dx + dy*dy + dz*dz;

            double c = m_clsnss*(m_rads[i]+sphrad);
            if( d2 < c*c ) return OUTSIDE;
            if( d2 < mind2 ) mind2 = d2;
        }

        return INCLUDE;
    }


    int out_solute::check_rndbox( const numvec& boxcnt, const numvec& boxext, double sphrad, double& mind ) const
    {
        double boxrad = norm(boxext)*0.5;
        int npartial = 0;
        for( int i=0; i < m_nsph; ++i )
        {
            double dx = boxcnt[0] - m_cnts[3*i];
            double dy = boxcnt[1] - m_cnts[3*i+1];
            double dz = boxcnt[2] - m_cnts[3*i+2];
            double d2 = dx*dx + dy*dy + dz*dz;

            double d = sqrt(d2);
            double c = m_clsnss*(m_rads[i]+sphrad);
            if( d+boxrad < c ) return OUTSIDE;
            if( d-boxrad < c ) npartial++;

            if( d < mind ) mind = d;
        }

        return (npartial > 0) ? PARTIAL : INCLUDE;
    }    

    ///////////////////////////////////////////////////////////////
    // 
    //   implementation of shl_solute
    //
    ///////////////////////////////////////////////////////////////

    shl_solute::shl_solute( solute_i& solute, double clsnss, double shlext )
        : m_inner( solute, clsnss ), 
          m_center(3),
          m_extent(3)
    {
        numvec reg = region( solute.getcord(), solute.getvdwr() );
        m_center = makevec( reg[0], reg[1], reg[2] );
        m_extent = makevec( reg[3], reg[4], reg[5] ) + scalar_numvec(3, 2*shlext);

        m_maxrad = *std::max_element( solute.getvdwr().begin(), solute.getvdwr().end() );
        m_clsnss = clsnss;
        m_shlext = shlext;
    }

    int shl_solute::check_sphere( const numvec& sphcnt, double sphrad ) const
    {
        double mind2 = DOUBLEMAX;
        if( m_inner.check_sphere(sphcnt, sphrad, mind2)==OUTSIDE )
        {
            return OUTSIDE;
        }

        double rout = m_shlext; //m_clsnss*(m_maxrad + sphrad) + m_shlext;
        return (mind2 < rout*rout) ? INCLUDE : OUTSIDE;
    }

    int shl_solute::check_rndbox( const numvec& cnt, const numvec& ext, double rad ) const
    {
        // mind: minimum distance
        double mind = DOUBLEMAX;
        int type = m_inner.check_rndbox( cnt, ext, rad, mind);
        if( type != INCLUDE )
        {
            return type;
        }

        double rout = m_shlext; // m_clsnss*(m_maxrad + rad) + m_shlext;
        if( mind + rad < rout ) return INCLUDE;
        if( mind < rout + rad ) return PARTIAL;
        return OUTSIDE;
    }

    //////////////////////////////////////////////
    // implementation of and_region
    /////////////////////////////////////////////

    and_region::and_region( const region_i& lhs, const region_i& rhs )
        : m_lhs(&lhs), m_rhs(&rhs), m_cnt(3), m_ext(3)
    {
        funstack_t::push( "and_region::and_region" );

        numvec cnt1 = subvec(lhs.center(), 0, 3);
        numvec ext1 = lhs.extent();

        numvec cnt2 = subvec(rhs.center(), 0, 3);
        numvec ext2 = rhs.extent();

        numvec min1 = cnt1 - 0.5*ext1;
        numvec max1 = cnt1 + 0.5*ext1;
        numvec min2 = cnt2 - 0.5*ext2;
        numvec max2 = cnt2 + 0.5*ext2;

        numvec min(3), max(3);
        for( int i=0; i < 3; ++i )
        {
            min[i] = std::max( min1[i], min2[i] );
            max[i] = std::min( max1[i], max2[i] );
        }

        m_cnt = (min + max) * 0.5;
        m_ext = (max - min);
	funstack_t::pop(); 
    }

    int and_region::check_sphere( const numvec& sphcnt, double sphrad ) const
    {
        int r1 = m_lhs->check_sphere( sphcnt, sphrad );
        int r2 = m_rhs->check_sphere( sphcnt, sphrad );

        return (r1==OUTSIDE || r2==OUTSIDE) ? OUTSIDE : INCLUDE;
    }

    int and_region::check_rndbox( const numvec& cnt, const numvec& ext, double rad ) const
    {
        //int r1 = INCLUDE; 
        int r1 = m_lhs->check_rndbox( cnt, ext, rad );
        int r2 = m_rhs->check_rndbox( cnt, ext, rad );

        if( r1==OUTSIDE || r2==OUTSIDE )
            return OUTSIDE;

        if( r1==PARTIAL || r2==PARTIAL )
            return PARTIAL;

        assert( r1==INCLUDE && r2==INCLUDE );
        return INCLUDE;
    }

    
} // namespace mort

