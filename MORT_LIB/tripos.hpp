#include <common.hpp>

namespace mort
{
    class molecule_t;

    namespace tripos
    {
        bool read_sect( istream& stream, molecule_t& mol );
        
        void write_atom( ostream& stream, const molecule_t& mol );
        
        void write_bond( ostream& stream, const molecule_t& mol );
        
        void write_resd( ostream& stream, const molecule_t& mol );
        
        void write_head( ostream& stream, const molecule_t& mol );
        
        void write_note( ostream& stream, const molecule_t& mol );    
    }

} // namespace mort


