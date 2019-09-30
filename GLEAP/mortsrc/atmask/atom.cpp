#include <common.hpp>
#include "atom.hpp"
#include "logic.hpp"
#include "condition.hpp"

namespace mort
{
    
    namespace atmask
    {

            atom_node::atom_node(bool reg)
                :node_i( reg )
            {
            }

            atom_node::~atom_node()
            {
            }
            
            atomvec_t atom_node::match( const molecule_t& mol ) const
            {
                atomvec_t result;
                mort::copy_if( mol.atom_begin(), mol.atom_end(), std::back_inserter(result), m_condition );
                return result;
            }
            
            bool atom_node::readable( const char* ptr ) const
            {
                return *ptr == '@';
            }
            
            const char* atom_node::parse( const char* ptr, shared_ptr< node_i >& curt ) const
            {
                assert( *ptr == '@' );
                ptr++;

                shared_ptr< atom_node > copy( new atom_node() );

                if( *ptr == '%' ) 
                {
                    copy->m_condition.parm = TYPE;
                    ptr++;
                }

                ptr = read_condition( ptr, copy->m_condition );
                


                if( curt == NULL )
                {
                    curt = copy;
                }
                else
                {
                    curt = shared_ptr< node_i >( new logic_node( AND, curt, copy ) );
                }

                return ptr;
            }

    } // namespace atmask

} // namespace mort

static const mort::atmask::atom_node g_atom_node(true);


