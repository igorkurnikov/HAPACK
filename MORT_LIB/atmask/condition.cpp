#include "condition.hpp"

namespace mort
{
    
    namespace atmask
    {
     
        const char* read_condition( const char* ptr, condition_t& condition )
        {
            string first = next_alnum( ptr );
            
            if( isdigit( first[0] ) ) condition.parm = ID;
            
            ptr = skip_alnum( ptr );
            
            if( *ptr == ',' )
            {
                condition.type = LIST;
                condition.list.push_back( first );
                
                while( *ptr == ',' )
                {
                    ptr++;
                    condition.list.push_back( next_alnum(ptr) );
                    ptr = skip_alnum( ptr );
                }
            }
            else if( *ptr == '-' )
            {
                ptr++;
                assert( isdigit( first[0] ) && condition.parm == ID );
                condition.type = RANGE;
                condition.lower = atoi( first.c_str() );
                condition.upper = atoi( next_alnum( ptr ).c_str() );
                ptr = skip_alnum( ptr );
            }
            else
            {
                condition.type = ONE;
                condition.word = first;
            }
            
            return ptr;
        }

    } // namespace atmask

} // namespace mort
