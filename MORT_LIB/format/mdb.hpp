#ifndef MORTSRC_FORMAT_MDB_HPP
#define MORTSRC_FORMAT_MDB_HPP

#include <iosfwd>

namespace mort
{
    class molecule_t;

    class database_t;

    /// \brief load database from input stream according to format code
    /// \param format integer code for different, could be SDF, MOL2, OFF
    /// \param db   pre-defined residue templates for PDB format
    /// throw an exception for unrecognized format.   
    void load_mdb( std::istream& is, database_t& db, const hashid_t& format );

    /// \brief save database to stream according to foramt
    /// \param format integer code for different, could be SDF, MOL2, OFF
    /// throw an exception for unrecognized format.   
    void save_mdb( std::ostream& os, const database_t& db, const hashid_t& format );

    /// \brief load database from file according to foramt
    /// \param format integer code for different, could be SDF, MOL2, OFF or UNKNOWN
    /// \param db   pre-defined residue templates for PDB format
    /// if format equals to UNKNOWN, it is determined by file name extension.
    /// throw an exception for unrecognized format.   
    void load_mdb( const std::string& file, database_t& db, const hashid_t& format=UNKNOWN );
 
    /// \brief save database to file according to foramt
    /// \param format integer code for different, could be SDF, MOL2, OFF or UNKNOWN
    /// if format equals to UNKNOWN, it is determined by file name extension.
    /// throw an exception for unrecognized format.   
    void save_mdb( const std::string& file, const database_t& db, const hashid_t& format=UNKNOWN );


}



#endif
