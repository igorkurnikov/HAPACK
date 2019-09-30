#ifndef MORTSRC_OBJFWD_HPP
#define MORTSRC_OBJFWD_HPP

namespace mort
{
    template <typename valu_t> class iter_T;
    template <typename iter_t, typename valu_t > class range_T;

    class root_t;
    class morf_t;
    class atom_t;
    class bond_t;
    class angl_t;
    class dihe_t;
    class impr_t;
    class resd_t;

    class molecule_t;
    class mcmpdata_t;
    class mcmprela_t;
    class database_t;

    typedef molecule_t mole_t;
    typedef database_t comp_t;
    typedef iter_T< morf_t > mobjiter_t;
    typedef iter_T< atom_t > atomiter_t;
    typedef iter_T< bond_t > bonditer_t;
    typedef iter_T< angl_t > angliter_t;
    typedef iter_T< dihe_t > diheiter_t;
    typedef iter_T< impr_t > impriter_t;
    typedef iter_T< resd_t > resditer_t;

    typedef range_T< mobjiter_t, morf_t > mobj_range;
    typedef range_T< atomiter_t, atom_t > atom_range;
    typedef range_T< bonditer_t, bond_t > bond_range;
    typedef range_T< angliter_t, angl_t > angl_range;
    typedef range_T< diheiter_t, dihe_t > dihe_range;
    typedef range_T< impriter_t, impr_t > impr_range;
    typedef range_T< resditer_t, resd_t > resd_range;

} // namespace mort

#endif

