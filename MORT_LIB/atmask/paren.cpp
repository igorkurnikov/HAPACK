#include "paren.hpp"

namespace mort
{
    
    namespace atmask
    {

            paren_node::paren_node( )
                : node_i( true )
            {
            }
          
            atomvec_t paren_node::match( const molecule_t& mol ) const
            {
                assert( false );
				atomvec_t av;
				return av;
            }
  
            bool paren_node::readable( const char* ptr ) const
            {
                return *ptr == '(';
            }
            
            const char* paren_node::parse( const char* ptr, shared_ptr< node_i >& curt ) const
            {
                assert( *ptr == '(' && curt == NULL );
                
                ptr++;
                
                ptr = read_until( ptr, ')', curt );
                
                assert( *ptr == ')' );
                
                ptr++;

                return ptr;
            }
            

        
    } // namespace atmask
    
} // namespace mort

static const mort::atmask::paren_node g_paren_node;
