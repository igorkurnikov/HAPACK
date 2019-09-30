#ifndef MORTSRC_OBJECT_BOND_HPP
#define MORTSRC_OBJECT_BOND_HPP

#include "morf.hpp"

namespace mort
{

    class bond_t : public morf_t
    {
    public:

        bond_t( const morf_t& rhs );

        bond_t( const molecule_t& m, int id );

        int order() const { return get_i(ORDER); }

        static bond_t create( atom_t& a1, atom_t& a2 );

        static bond_t frcget( atom_t& a1, atom_t& a2 );

        static bond_t get( const atom_t& a1, const atom_t& a2 );

        static bool has( const atom_t& a1, const atom_t& a2 );

        static bool get( const atom_t& a1, const atom_t& a2, bond_t& b );

        static bool has( const resd_t& r1, const resd_t& r2 );
    };

} // namespace mort

#endif
