#include <object.hpp>

#include "parm.hpp"

namespace mort
{

    struct parmset_t;
    
    namespace prmtop
    {        

        void set_les(morf_t& obj );

        void les_scale(molecule_t& mol, parmset_t& ps);

        bool les_forbid(const morf_t& atm1, const morf_t& atm2);
        
        void les_exclude(const morf_t& atom, vector<int>& list, vector<int>& dist);
        
    } // namespace prmtop

} // namespace mort

        

