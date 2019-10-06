#ifndef MORT_CAPBOX_OCTREE_HPP
#define MORT_CAPBOX_OCTREE_HPP

#include <object.hpp>

namespace mort
{
    
    namespace capbox
    {
        class shaper_i;

        typedef function< double ( double, double, double ) > potfun_t;

        class octree_t
        {
        public:

            octree_t( const numvec& corner, int ndepth, double resolution );
            
            virtual ~octree_t();

            void split( );

            int make_shape( shaper_i& shapefun, double ion_radius );

            void calculate( potfun_t potfun, std::pair< double, numvec >& min, std::pair< double, numvec >& max );

        private:

            shared_ptr< octree_t > m_subs[8];
            
            numvec m_corner;

            numvec m_center;

            int m_typeid;

	    int m_ndepth;

	    int m_ngrid;
           
            double m_length;

            double m_radius;

	    double m_resolution;
            
            vector<double> m_subpol;

        };

        typedef shared_ptr< octree_t > octree_ptr;
        
        struct elepot_t
        {
            elepot_t( const double* pcord, const double* pchrg, int natom )
		: m_spos(0)
            {
		m_natom = natom;
		m_pchg = pchrg;
		m_pcrd = pcord;
            }

            elepot_t( const numvec& pos, double chg )
                : m_spos( pos ), m_schg(chg)
            {
                m_natom = 0;
            }
            
            double operator()( double x, double y, double z )
            {
                if( m_natom==0 )
                {
                    assert( m_spos.size()==3u );
                    double dx = m_spos[0] - x;
                    double dy = m_spos[1] - y;
                    double dz = m_spos[2] - z;
                    double d2 = dx*dx + dy*dy + dz*dz;
                    return m_schg/d2;
                }
                

                double elep = 0.0;
                for( int i=0; i < m_natom; ++i )
                {
                    double dx = m_pcrd[3*i  ] - x;
		    double dy = m_pcrd[3*i+1] - y;
		    double dz = m_pcrd[3*i+2] - z;
                    double d2 = dx*dx + dy*dy + dz*dz;
                    elep = elep + m_pchg[i] / d2;
                }
                
                return elep;
            }
            
            int m_natom;
	    const double* m_pchg;
	    const double* m_pcrd;

            // this is for single atom: 's' stands for single.
            numvec m_spos;
            double m_schg;
        };

    } // namespace capbox

} // namespace mort

#endif

