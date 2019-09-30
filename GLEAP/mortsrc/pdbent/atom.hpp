#ifndef MORTSRC_PDBENT_ATOM_HPP
#define MORTSRC_PDBENT_ATOM_HPP

namespace mort
{

    namespace pdbent
    {
        enum pdbatomfield_e
        {
            ATOM_SEQ    = 6,
            ATOM_NAME   = 12,
            RESD_TYPE   = 17,
            RESD_CHAIN  = 21,
            RESD_SEQ    = 22,
            RESD_ICODE  = 26,
            ATOM_X      = 30,
            ATOM_Y      = 38,
            ATOM_Z      = 46
        };

       
        void read_atom( istream& is, resd_t& r );

        string read_atom_name( const string& line );

        void write_atom( ostream& os, const atom_t& a );

        void write_atom_name( ostream& os, const atom_t& a );

        void write_resd_info( ostream& os, const atom_t& a );
    
    } // namepace pdbent

} // namespace mort

#endif
