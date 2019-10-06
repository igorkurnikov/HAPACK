#ifndef MORT_ATMASK_CONDITION_HPP
#define MORT_ATMASK_CONDITION_HPP

#include <common.hpp>
#include <object.hpp>

namespace mort
{
    namespace atmask
    {
        enum type_e { ONE, LIST, RANGE };
        struct condition_t
        {
            condition_t()
            {
                parm = NAME;
            }
            
            bool match_one(const morf_t& obj, const string& value) const
            {
                if( parm==ID )
                {
                    assert( isdigit( value[0] ) );
                    return obj.get_i(ID) == atoi( value.c_str() );
                }

                assert(parm == NAME || parm == TYPE);
                return (parm == NAME) ? obj.get_s(NAME)==value : obj.get_s(TYPE)==value;
            }

            bool operator()(const morf_t& obj) const
            {
                if( type == ONE )
                {
                    return match_one(obj, word);
                }            

                if( type == LIST )
                {
                    atomvec_t result;

                    for( int i=0; i < (int)list.size(); ++i )
                    {
                        if( match_one(obj, list[i]) )
                        {
                            return true;
                        }
                    }
                    
                    return false;
                }
                
                assert( type == RANGE && parm == ID );

                int objid = obj.get_i(ID);
                return objid >= lower && objid <= upper;
            }

            int type;
            
            int parm;
            
            int upper;
            
            int lower;

            string word;
            
            vector< string > list;
        };

        const char* read_condition( const char* ptr, condition_t& condition );
            
    } // namespace atmask

} // namespace mort

#endif

