#ifndef MORTSRC_PDBENT_RESD_HPP
#define MORTSRC_PDBENT_RESD_HPP

#include <iosfwd>

namespace mort
{
    class molecule_t;

    namespace pdbent
    {

        resd_t read_resd( std::istream& is, molecule_t& m );

        bool bond_resd( resd_t& ri, resd_t& rj );

        void write_resd( std::ostream& os, const resd_t& resd );

        void setbgn( resd_t& r );

        void setend( resd_t& r );

        void setpos( resd_t& rprev, resd_t& rcurt );

    } // namespace pdbent

} // namespace mort

#endif


