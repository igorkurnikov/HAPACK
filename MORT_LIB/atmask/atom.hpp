#ifndef MORTSRC_ATMASK_ATOM_HPP
#define MORTSRC_ATMASK_ATOM_HPP

#include "mask.hpp"
#include "condition.hpp"

namespace mort
{
    namespace atmask
    {
        class atom_node : public node_i
        {
        public:

            atom_node(bool reg = false);
        
            virtual ~atom_node();
            
            virtual atomvec_t match( const molecule_t& mol ) const;
            
            virtual bool readable( const char* ptr ) const;
            
            virtual const char* parse( const char* ptr, shared_ptr< node_i >& curt ) const;

        private:

             condition_t m_condition;
        };
             
        
    } // namespace atmask

} // namespace mort

#endif

