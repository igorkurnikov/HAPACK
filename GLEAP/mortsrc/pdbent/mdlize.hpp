#ifndef MORTSRC_PDBENT_MDLIZE_HPP
#define MORTSRC_PDBENT_MDLIZE_HPP

#include <map>

namespace mort
{
    class atmvec;

    class namemap_t;

    class molecule_t;

    class database_t;

    void mdlize_mdb( molecule_t& m, const database_t& mdb );

    void mdlize_seq( molecule_t& m, const database_t& seq, const namemap_t& nmap );

    resd_t mdlize_resd( resd_t& r, const database_t& mdb, molecule_t& full, std::map<int,int>& idmap );

    resd_t mdlize_resd( resd_t& r, const molecule_t& mdl, const namemap_t& n, molecule_t& full, std::map<int,int>& idmap );

    void conect_resd( resd_t& r1, resd_t& r2 );

    void mdlize_auto( resd_t& r );

    numvec impose( const atom_t& a, const std::map<string, numvec>& pos );

    atmvec bfv( const atom_t& a, const std::map<string, numvec>& pos );
 
} // namespace mort

#endif
