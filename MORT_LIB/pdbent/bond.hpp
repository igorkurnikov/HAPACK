#ifndef MORTSRC_PDBENT_BOND_HPP
#define MORTSRC_PDBENT_BOND_HPP

#include <set>
#include <iosfwd>

namespace mort
{
    class molecule_t;

    namespace pdbent
    {

        void read_bond( std::istream& is, molecule_t& m );

        void mark_bond( const molecule_t& m, std::set<int>& ids );
    }


    class console_t;

    void disulfide( molecule_t& m, double cutoff, console_t* pcon );

} // namespace mort

#endif
