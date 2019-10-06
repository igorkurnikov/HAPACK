#ifndef MORT_ATMASK_PATTERN_HPP
#define MORT_ATMASK_PATTERN_HPP

#include <object.hpp>

namespace mort
{
    namespace atmask
    {

        class node_i
        {
        public:

            node_i( bool reg );
            
            virtual ~node_i();

            virtual atomvec_t match( const molecule_t& mol ) const = 0;
            
            virtual bool readable( const char* ptr ) const = 0;

            virtual const char* parse( const char* ptr, shared_ptr< node_i >& curt ) const = 0;

        };

        const char* read_node( const char* ptr, shared_ptr< node_i >& curt );
 
        const char* read_until( const char* ptr, char end, shared_ptr< node_i >& curt );

        shared_ptr< node_i > read_mask( const string& mask );

    } // namespace atmaks

} // namespace mort

#endif
