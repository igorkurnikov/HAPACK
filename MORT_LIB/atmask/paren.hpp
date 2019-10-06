#include "mask.hpp"

namespace mort
{
    
    namespace atmask
    {

        class paren_node : public node_i
        {
        public:

            paren_node();
                      
            virtual atomvec_t match( const molecule_t& mol ) const;
  
            virtual bool readable( const char* ptr ) const;
            
            virtual const char* parse( const char* ptr, shared_ptr< node_i >& curt ) const;
                           
        };

    } // namespace atmask

} // namespace mort

