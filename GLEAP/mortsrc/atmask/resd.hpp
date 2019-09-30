#ifndef MORTSRC_ATMASK_RESD_HPP
#define MORTSRC_ATMASK_RESD_HPP

#include "mask.hpp"
#include "condition.hpp"

namespace mort
{
    
    namespace atmask
    {

        class resd_node : public node_i
        {
        public:

            resd_node(bool reg=false);
            
            virtual ~resd_node();

            virtual atomvec_t match( const molecule_t& mol ) const;
                        
            virtual bool readable( const char* ptr ) const;
                        
            virtual const char* parse( const char* ptr, shared_ptr< node_i >& curt ) const;
            
        private:
            
            condition_t m_condition;
        };

    } // namespace atmask

} // namespace mort

#endif

