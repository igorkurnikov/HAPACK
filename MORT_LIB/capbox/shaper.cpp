#include <algorithm>
#include <object.hpp>
#include "addions.hpp"
#include "shaper.hpp"

namespace mort
{
    namespace capbox
    {

        enum result_e { INNER, PART_INNER, SHELL, PART_SHELL, OUTER };

        solvent_shaper::solvent_shaper( ionee_i& solute, double closeness, double outer_radius )
            : m_natom( solute.solute_natom() ),
              m_coord( solute.getcord() ),
              m_pvdwr( solute.getvdwr() ),
              m_closeness( closeness ),
              m_outer_radius( outer_radius ),
              m_exclusion( solute.solute_natom(), 0 )
        {
        }

        void solvent_shaper::make_exclusion_list( double radius )
        {
            // the exclusion list cotains the ID of spheres who are totally
            // out of the previous node. If the previous node happen to be 
            // parent of the current node,  its exclusion list should be 
            // inherited since children is always inside the parent.

            // go back until find my father, whose radius is bigger than mine
            // or there might not be one.
            int old_depth = m_exclusion_radii.size();
            while( m_exclusion_radii.size()>0 && radius>=m_exclusion_radii.back() )
            {
                m_exclusion_radii.pop_back();
            }


            // clear the exclusion of my brothers and their descents. leave only 
            // the exclusion of my father and his parents 
            int new_depth = m_exclusion_radii.size();
            if( new_depth != old_depth )
            {
                assert( new_depth < old_depth );
                assert( (int)m_exclusion.size()==m_natom );

                for(int is=0; is < m_natom; ++is )
                {
                    if( m_exclusion[is] > new_depth ) 
                        m_exclusion[is] = 0;
                }
            }
        }


        // atom is a single point. It is either INCLUDED or EXCLUDED, cannot 
        // be PARTIAL.
        int solvent_shaper::check_atm( const numvec& pos, double rvdw )
        {
           
            make_exclusion_list( 0.0 );

            double mind2 = 1e20;
            for( int i=0; i < m_natom; ++i )
            {
                // if i is in my father's exclusion list, skip it
                if( m_exclusion[i]> 0 )
                {
                    continue;
                }

                double inner = (m_pvdwr[i] + rvdw) * m_closeness;
                double inner2 = inner * inner;

                double dx = pos[0] - m_coord[3*i  ];
                double dy = pos[1] - m_coord[3*i+1];
                double dz = pos[2] - m_coord[3*i+2];
                double d2 = dx*dx + dy*dy + dz*dz;
                
                if( d2 < inner2 )
                {
                    return EXCLUDED;
                }
                
                
                if( i==0 || d2 < mind2 )
                {
                    mind2 = d2;
                }
            } 

            if( m_outer_radius==0.0 ) return INCLUDED;

            return ( mind2 < m_outer_radius*m_outer_radius ) ? INCLUDED : EXCLUDED;
        }


        int solvent_shaper::check_box( const numvec& center, const numvec& size, const double rsolvent  )
        {
            double radius = sqrt( size[0]*size[0] + size[1]*size[1] + size[2]*size[2] ) / 2.0;
   
            make_exclusion_list( radius );
        

            vector<int> num(5, 0);

            vector<int> own_exclusion;
            for( int i=0; i < m_natom; ++i )
            {
                // if i is in my father's exclusion list, skip it
                if( m_exclusion[i]> 0 )
                {
                    continue;
                }


                int type = check_each( center, radius, rsolvent, i );
            
                // if I am sucked inside any sphere, I must be excluded
                if( type == INNER )
                {
                    return EXCLUDED;
                }
            

                num[type] += 1;


                if( m_outer_radius > 0.0 && type == OUTER )
                {
                    own_exclusion.push_back( i );
                }
                
                if( m_outer_radius == 0.0 && type ==SHELL )
                {
                    own_exclusion.push_back( i );
                }

            }

            if( own_exclusion.size()>0 )
            {
                m_exclusion_radii.push_back( radius );

                int ndepth = m_exclusion_radii.size();
            
                for( int i=0; i < (int)own_exclusion.size(); ++i )
                {
                    int is = own_exclusion[i];
                
                    // if a sphere is in my own exclusion list, it cannot be in my father's exclusion
                    // list, since I will skip it.
                    assert( m_exclusion[is]==0 );
                    m_exclusion[is] = ndepth;
                }

            }
                


            assert( num[INNER]==0 );

            if( num[PART_INNER] > 0 )
            {
                // intersect with an atom, must be partial.
                return PARTIAL;
            }

            if( num[SHELL] > 0 )
            {
                // inside some atoms' shell, must be included? taken from old leap code
                return INCLUDED;
            }
        
            // only possible: PART_SHELL, OUTER, if no part_shell, all outer
            // should be excluded
            return (num[PART_SHELL]==0) ? EXCLUDED : PARTIAL;
        }

        int solvent_shaper::check_each( const numvec& center, double radius, double rsolvent, int iatom ) const
        {
            double dx = center[0] - m_coord[3*iatom];
            double dy = center[1] - m_coord[3*iatom+1];
            double dz = center[2] - m_coord[3*iatom+2];
            double d = sqrt( dx*dx + dy*dy + dz*dz );

            double inner_radius = (m_pvdwr[iatom] + rsolvent)*m_closeness;


            assert( radius > 0.0 );
            
            if( d + radius < inner_radius ) return INNER;

            if( d - radius < inner_radius ) return PART_INNER;

            if( m_outer_radius==0.0 ) return SHELL;

            if( d + radius < m_outer_radius ) return SHELL;

            if( d - radius < m_outer_radius ) return PART_SHELL;

            return OUTER;
        }

        outsphere_shaper::outsphere_shaper( const numvec& sphere, double radius )
            :m_sphere( sphere ), m_radius(radius)
        {
        }

        int outsphere_shaper::check_atm( const numvec& pos, double rvdw )
        {
            double dx = pos[0] - m_sphere[0];
            double dy = pos[1] - m_sphere[1];
            double dz = pos[2] - m_sphere[2];
            double d2 = dx*dx + dy*dy + dz*dz;

            double tot_radius = m_radius + rvdw;
            return (d2 < tot_radius*tot_radius)? EXCLUDED : INCLUDED;
        }
        
            

        int outsphere_shaper::check_box( const numvec& center, const numvec& size, double rsolvent )
        {
            double radius = sqrt( size[0]*size[0] + size[1]*size[1] + size[2]*size[2] ) / 2.0;
            double dx = center[0] - m_sphere[0];
            double dy = center[1] - m_sphere[1];
            double dz = center[2] - m_sphere[2];
            double d = sqrt( dx*dx + dy*dy + dz*dz );

            double tot_radius = m_radius + rsolvent;

            assert( radius > 0.0 );


            if( d+radius < tot_radius )
                return EXCLUDED;

            if( d-radius < tot_radius )
                return PARTIAL;

            return INCLUDED;
        }



    } // namespace capbox

} // namespace mort

   
