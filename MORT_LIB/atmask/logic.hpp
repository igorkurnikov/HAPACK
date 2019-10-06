#ifndef MORT_ATMASK_LOGIC_HPP
#define MORT_ATMASK_LOGIC_HPP

#include "mask.hpp"

namespace mort
{
    
    namespace atmask
    {

        class logic_node : public node_i
        {
        public:

            logic_node( bool reg = false );

            logic_node( const hashid_t& type, shared_ptr< node_i > left, shared_ptr< node_i > right );

            virtual ~logic_node();
            
            virtual atomvec_t match( const molecule_t& mol ) const;

            virtual bool readable( const char* ptr ) const;
            
            virtual const char* parse( const char* ptr, shared_ptr< node_i >& curt ) const;
            
            hashid_t m_type;
            
            shared_ptr< node_i > m_left;
            
            shared_ptr< node_i > m_right;
        };

    } // namespace atmask
    
} // namespace mort

#endif
