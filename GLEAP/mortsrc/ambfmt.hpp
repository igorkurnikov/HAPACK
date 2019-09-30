#ifndef MORT_AMBFMT_HPP
#define MORT_AMBFMT_HPP

#include <object.hpp>

#include "ambfmt/prmtop.hpp"
#include "ambfmt/amboff.hpp"
#include "ambfmt/ambfrc.hpp"

namespace mort
{
    void assign_ffparm(molecule_t& mol, const molecule_t& ff);

    void parametrize(molecule_t& mol, const molecule_t& ff, parmset_t& ps);

    void read_amber_prep(istream& is, molecule_t& mol);

    void read_amber_prep(istream& is, database_t& mdb);
 
} // namespace mort
    
#endif

