#ifndef MORTSRC_FORMAT_MOL_HPP
#define MORTSRC_FORMAT_MOL_HPP

#include <iosfwd>
#include <string>

namespace mort
{
    class molecule_t;

    class database_t;    

    /// \brief load molecule from input stream according to format code
    /// \param format integer code for different, could be SDF, PDB, MOL2, OFF
    /// throw an exception for unrecognized format.   
    void load_mol( std::istream& is, molecule_t& mol, const hashid_t& format );

    /// \brief save molecule to stream according to foramt
    /// \param format integer code for different, could be SDF, PDB, MOL2, OFF
    /// throw an exception for unrecognized format.   
    void save_mol( std::ostream& os, const molecule_t& mol, const hashid_t& format );

    /// \brief load molecule from file according to foramt
    /// \param format integer code for different, could be SDF, PDB, MOL2, OFF or UNKNOWN
    /// if format equals to UNKNOWN, it is determined by file name extension.
    /// throw an exception for unrecognized format.   
    void load_mol( const std::string& file, molecule_t& mol, const hashid_t& format=UNKNOWN );    

    /// \brief save molecule to file according to foramt
    /// \param format integer code for different, could be SDF, PDB, MOL2, OFF or UNKNOWN
    /// if format equals to UNKNOWN, it is determined by file name extension.
    /// throw an exception for unrecognized format.   
    void save_mol( const std::string& file, const molecule_t& mol, const hashid_t& format=UNKNOWN );
}

#endif

