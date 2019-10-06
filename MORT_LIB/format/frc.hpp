#ifndef MORTSRC_AMBFMT_AMBFRC_HPP
#define MORTSRC_AMBFMT_AMBFRC_HPP

#include <iosfwd>

namespace mort
{

    class molecule_t;

    /// \brief read force field paramter from an amber's paramter file, 
    ///
    /// the paramter is stored in a molecule 
    void read_frc( std::istream& is, molecule_t& fp );


    /// \brief read amoeba force field paramter from an tinker's paramter. 
    /// \param   is input stream
    /// \atomff  atomic part of force field parameters
    /// \poleff  dipole part of force field parameters
    /// the paramter is stored in two molecules, the dipole parameters are stored
    /// in poleff, and the rest part of it are stored in atomff, including vdw,
    /// bond, angle, torsion, out-of-plane bending, torsion-torsion, PI-torsion,
    /// urey-bradley, bond-angle interaction
    void read_amoeba_frc( std::istream& is, molecule_t& afp, molecule_t& pfp );
}

#endif

