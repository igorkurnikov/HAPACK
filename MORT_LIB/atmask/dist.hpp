#ifndef MORTSRC_ATMASK_DIST_HPP
#define MORTSRC_ATMASK_DIST_HPP


#include "mask.hpp"

namespace mort
{

    namespace atmask
    {

        class dist_node : public node_i
        {
        public:

            dist_node( bool reg = false );
            
            virtual ~dist_node();
            
            virtual atomvec_t match( const molecule_t& mol ) const;
            
            virtual bool readable( const char* ptr ) const;
                        
            virtual const char* parse( const char* ptr, shared_ptr< node_i >& curt ) const;

        private:

            shared_ptr< node_i > m_core;

            int m_type;
            
	    int m_object;

            double m_cutoff;
        };

    } // namspace atmask

} // namespace mort


#endif

