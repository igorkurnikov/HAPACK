#include <common.hpp>
#include <object.hpp>

namespace mort
{
    class molecule_t;

    /// \defgroup atmask Atmask: Atom selection using amber Mask
    ///
    /// using amber mask to selection atoms in a molecule

    /// \brief return a subset of atoms which is selected by the mask
    /// \ingroup atmask
    atomvec_t mask_atom(const molecule_t& mol, const string& mask);

    /// store the mask result in an array
    /// \ingroup atmask
    /// for each atom, if it fit for the input mask, maskvec[atom.getoid()] = 1
    ///  otherwise, maskvec[atom.getoid()] = 0
    void mask_atom( const molecule_t& mol, const string& mask, vector<int>& maskvec );

    atomvec_t smarts_mask_atom(const molecule_t& mol, const string& mask);

    void smarts_mask_atom(const molecule_t& mol, const string& mask, vector<int>& maskvec);
}
