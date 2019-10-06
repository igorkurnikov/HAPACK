#include <object.hpp>
#include "resd.hpp"

namespace mort
{
    
    namespace atmask
    {

            resd_node::resd_node(bool reg)
                :node_i(reg)
            {
            }
            
            resd_node::~resd_node()
            {
            }

            atomvec_t resd_node::match( const molecule_t& mol ) const
            {
                atomvec_t result;
                copy_if( mol.resd_begin(), mol.resd_end(), atom_pusher(result), m_condition );
                return result;
            }
            
            bool resd_node::readable( const char* ptr ) const
            {
                return *ptr == ':';
            }
            
            const char* resd_node::parse( const char* ptr, shared_ptr< node_i >& curt ) const
            {
                assert( readable(ptr) && curt == NULL );
                ptr++;

                shared_ptr< resd_node > copy( new resd_node() );
                copy->m_condition.parm = TYPE;
                ptr = read_condition(ptr, copy->m_condition);
                curt = copy;
                return ptr;
            }

    } // namespace atmask

} // namespace mort

static const mort::atmask::resd_node g_resd_node( true );
