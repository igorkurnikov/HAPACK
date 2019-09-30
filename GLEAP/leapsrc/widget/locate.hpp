#ifndef MORTSRC_OBJECT_LOCATE_HPP
#define MORTSRC_OBJECT_LOCATE_HPP

#include <common.hpp>



namespace mort
{
    class atomvec_t;

    class molecule_t;

    class database_t;

    bool locate_atom(const molecule_t& mol, const numvec& rect, atomvec_t& atms, int policy );

    bool locate_bond(const molecule_t& mol, const numvec& rect, bondvec_t& bnds, int policy );

    bool locate_atom( database_t& mdb, const numvec& rect, atomvec_t& atms, int policy );

    bool locate_bond( database_t& mdb, const numvec& rect, bondvec_t& bnds, int policy );
 
}

#endif
