#ifndef MORT_OBJECT_BONDFUN_HPP
#define MORT_OBJECT_BONDFUN_HPP
#include <string>

namespace mort
{
    class atom_t;
    class bond_t;
    class resd_t;
    class atmvec;
    class molecule_t;

    bool has_prev(const morf_t& tors );

    // create bond
    morf_t create_tor2(atmvec& atms);
    morf_t create_ptor(atmvec& atms, int period );

    // make bond by distance. Two atoms will be connected by a bond if their
    // distance < cutoff.
    void bond_bydis(resd_t& resd, double cutoff);
    void bond_bydis(molecule_t& mol, double cutoff);
    bool bond_bydis(atom_t& a1, atom_t& a2, double cutoff);

} // namespace mort
 
#endif


