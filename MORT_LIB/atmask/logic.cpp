#include "logic.hpp"

namespace mort
{
    
    namespace atmask
    {

        logic_node::logic_node( bool reg )
            :node_i( reg )
        {
        }

        logic_node::~logic_node()
        {
        }
            
        logic_node::logic_node( const hashid_t& type, shared_ptr< node_i > left, shared_ptr< node_i > right )
            :node_i( false )
        {
            assert( type == AND || type == OR );

            m_type = type;
            
            m_left = left;
            
            m_right = right;
        }
            
        atomvec_t logic_node::match( const molecule_t& mol ) const
        {
            assert( m_type == AND || m_type == OR || m_type == NOT );
                
            atomvec_t result;

            if( m_type == AND )
            {
                atomvec_t a = m_left ->match( mol );
                atomvec_t b = m_right->match( mol );
                set_intersection( a.begin(), a.end(), b.begin(), b.end(), back_inserter( result ) );
            }
            else if( m_type == OR )
            {
                atomvec_t a = m_left ->match( mol );
                atomvec_t b = m_right->match( mol );
                set_union( a.begin(), a.end(), b.begin(), b.end(), back_inserter( result ) );
            }
            else
            {
                atomvec_t a = m_left ->match( mol );
                set_difference( mol.atom_begin(), mol.atom_end(), a.begin(), a.end(), back_inserter(result) );
            }
            
            return result;
        }

        bool logic_node::readable( const char* ptr ) const
        {
            return *ptr == '&' || *ptr == '!' || *ptr == '|';
        }
            
        const char* logic_node::parse( const char* ptr, shared_ptr< node_i >& curt ) const
        {
            assert( *ptr == '&' || *ptr == '|' || *ptr == '!' );

            shared_ptr< logic_node > copy( new logic_node() );
                
            if( *ptr == '&' || *ptr == '|' )
            {
                assert( curt != NULL );

                copy->m_type = ( *ptr == '&' ) ? AND : OR;
                    
                ptr++;

                copy->m_left = curt;
                    
                ptr = read_node( ptr, copy->m_right );
                    
                curt = copy;

                return ptr;
            }
            else
            {
                assert( curt == NULL );
                    
                copy->m_type = NOT;
                
                ptr++;
                    
                ptr = read_node( ptr, copy->m_left );
                    
                curt = copy;
                
                return ptr;
            }
        }

        static const logic_node g_logic_node( true );

    } // namespace atmask
    
} // namespace mort
