#include <common.hpp>
#include <object.hpp>
#include "dist.hpp"

namespace mort
{
    
    namespace atmask
    {
 
        struct dist_atom
        {
            dist_atom( atomvec_t& core, int type, double cutoff )
            {
                m_core = &core;
                m_type = type;
                m_cutoff = cutoff;
            }
            
            bool for_atom(const morf_t& atom)
            {
                atomvec_t::iterator si = m_core->begin();
                for( ; si != m_core->end(); ++si )
                {
                    if( dis2(atom, *si) < m_cutoff*m_cutoff )
                    {
                        return m_type == '<';
                    }
                }
                
                return m_type == '>';
            }
            
            bool for_resd( const morf_t& resd )
            {
                atomiter_t ai = resd.atom_begin();
                for( ; ai != resd.atom_end(); ++ai )
                {
                    atomvec_t::iterator si = m_core->begin();
                    for( ; si != m_core->end(); ++si )
                    {
                        if( dis2( *ai, *si ) < m_cutoff*m_cutoff )
                        {
                            return m_type == '<';
                        }
                    }
                }
                
                return m_type == '>';
            }

	    bool operator()(const morf_t& p)
	    {
	        if( p.cmpid()==ATOM )
		{
		    return for_atom(p);
		}
		else 
		{
		    assert( p.cmpid()==RESD );
		    return for_resd(p);
		}
            }

            atomvec_t* m_core;
            
            int m_type;
            
            double m_cutoff;
        };

            dist_node::dist_node(bool reg)
                : node_i( reg )
            {
            }

            dist_node::~dist_node()
            {
            }
            
            atomvec_t dist_node::match( const molecule_t& mol ) const
            {
                atomvec_t result;
                
                atomvec_t core = m_core->match( mol );

                if( m_object == ATOM )
                {
                    copy_if( mol.atom_begin(), mol.atom_end(), atom_pusher(result), dist_atom(core, m_type, m_cutoff) );
                }
                else
                {
                    copy_if( mol.resd_begin(), mol.resd_end(), atom_pusher(result), dist_atom(core, m_type, m_cutoff) );
                }    

                return result;
            }
            
            bool dist_node::readable( const char* ptr ) const
            {
                return *ptr == '<' || *ptr == '>';
            }
            
            const char* dist_node::parse( const char* ptr, shared_ptr< node_i >& curt ) const
            {
                assert( curt != NULL );
                shared_ptr<dist_node> copy( new dist_node() );
                copy->m_core = curt;

                assert( ( *ptr == '<' || *ptr == '>' ) && curt != NULL );
                copy->m_type = *ptr;
                ptr++;

                assert( *ptr == '@' || *ptr == ':' );
                copy->m_object = (*ptr == '@') ? ATOM : RESD;
                ptr++;
                
                copy->m_cutoff = atof( ptr );
                ptr = skip_float( ptr );

                curt = copy;
                return ptr;                
            }
           
        
    } // namespace atmask
    
} // namespace mort


static const mort::atmask::dist_node g_dist_node(true);

    
