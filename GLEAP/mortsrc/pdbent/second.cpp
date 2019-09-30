#include <common.hpp>
#include <object.hpp>
#include "second.hpp"

namespace mort
{
    using std::find_if;

    namespace pdbent
    {
        
        void read_second( istream& is, second_s& second )
        {
            char line[MAX_LINE_WIDTH];
            is.getline( line, MAX_LINE_WIDTH );
            string keyw( line, line + 5 );

            if( keyw == "HELIX"  )
            {
                string start( line + 14, line + 25 );
                start.erase( std::remove( start.begin(), start.end(), ' ' ), start.end() );

                string end( line + 26, line + 37 ); 
                end.erase  ( std::remove( end.begin(),   end.end(),   ' ' ), end.end() );

                second.start = start;
                second.end   = end;
                second.type  = HELIX;
            }
            else
            {
                assert( keyw == "SHEET" );

                string start( line + 16, line + 26 );
                start.erase( std::remove( start.begin(), start.end(), ' ' ), start.end() );

                string end( line + 28, line + 38 );
                end.erase  ( std::remove( end.begin(),   end.end(),   ' ' ), end.end() );
                
                second.start = start;
                second.end   = end;
                second.type  = SHEET;
            }

        }
        
        void apply_second( const second_s& second, molecule_t& mol )
        {
            resditer_t start = find_if( mol.resd_begin(), mol.resd_end(), sparm_cmper1(NAME, second.start) );
            resditer_t end   = find_if( start+1, mol.resd_end(), sparm_cmper1(NAME, second.end) );

            assert( start != mol.resd_end() && end != mol.resd_end() );
            assert( second.type == HELIX || second.type == SHEET );

            for( ; start != end; ++start )
            {
                start->set_i(SECOND, second.type);
            }

            // start.set_i(SECOND, second.type);
        }
        
    } // namespace pdbent

} // namespace mort

