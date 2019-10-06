#include <enefrc/nonbond.hpp>
#include <enefrc/ctrlparm.hpp>

namespace mort
{
    // bonding energy
    double eval_bond(const molecule_t& m);
    double eval_angl(const molecule_t& m);
    double eval_tors(const molecule_t& m);
    double eval_oops(const molecule_t& m);

    // nonbond energy, non-periodic 
    numvec nonbond_egb(const molecule_t& m, const ctrlparm_t& p);
    numvec nonbond_egb(const atomvec_t& vec1, const atomvec_t& vec2, const ctrlparm_t& p);

    // nonbond energy, periodic boundary condition
    numvec get_dir(const molecule_t& m, const ctrlparm_t& p);
    double get_rec(const molecule_t& m, const ctrlparm_t& p);

    
}

