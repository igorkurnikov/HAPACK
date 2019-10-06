#ifndef MORTSRC_FORMAT_HPP
#define MORTSRC_FORMAT_HPP

#include <iosfwd>
#include <common.hpp>
#include <format/ext.hpp>
#include <format/mol.hpp>
#include <format/mdb.hpp>
#include <format/top.hpp>

namespace mort
{
    class molecule_t;
    class database_t;
    /** \defgroup format Format: molecular IO functions
     *  This group contain functions related to molecule IO in differnet format
     *  @{
     */

    /// \brief read molecule from a mdl's SDF format stream.
    bool read_sdf(std::istream& is, molecule_t& mol);
    
    bool read_rdf(std::istream& is, molecule_t& mol);
    
    /// \brief read molecule from a brookhaven's PDB format stream.
    /// \param db pre-defined residue templates
    //
    bool read_pdb(std::istream& is, molecule_t& mol);

    /// \brief read molecule from a tripos's MOL2 format stream.
    bool read_mol2(std::istream& is, molecule_t& mol);
    
    /// \brief read molecule from a SMILES string.
    bool read_smiles(const char* ptr, molecule_t& mol);
    
    /// \brief read molecule from an amber's OFF format stream
    bool read_off(std::istream& is, molecule_t& mol);
    

    /// \brief read molecule from a SMARTS string.
    ///
    /// the string is interpreted to a pattern molecule which can be used with find_subst
    /// for molecular pattern recognition.
    bool read_smarts( const char* ptr, molecule_t& mol );

    /// \brief dump molecule to a stream in SDF format.
    void write_sdf( std::ostream& os, const molecule_t& mol );

    /// \brief dump molecule to a stream in SDF format.
    void write_mol2( std::ostream& os, const molecule_t& mol );
    
    /// \brief dump molecule to a stream in PDB format.
    bool write_pdb( std::ostream& os, const molecule_t& mol );
 
    /// \brief dump molecule to a stream in OFF format
    void write_off( std::ostream& os, const molecule_t& mol );
 
    /// \brief write molecule's prmtop file, used by sander for MD simulation
    /// \param os      output stream
    /// \param mol     target molecule
    /// \param ff      force field parameter
    /// \sa parametrize
    void write_amber_prmtop( std::ostream& os, molecule_t& m, const molecule_t& ff );
 
    /// \brief write molecule's prmtop file, used by sander for MD simulation
    /// \param os      output stream
    /// \param mol     target molecule
    /// \param params  all the force field paramters appeared in the molecule
    /// \param atomff  atomic part of amoeba force field parameters
    /// \param poleff  dipole part of amoeba force field parameters
    void write_amoeba_prmtop( std::ostream& os, molecule_t& m, const molecule_t& aff, const molecule_t& pff );
 
   

    /** @} */
}

#endif

