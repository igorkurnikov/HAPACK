#ifndef DDUTILS_MORT_DAYLIGHT_STATUS_H
#define DDUTILS_MORT_DAYLIGHT_STATUS_H

#include <object.hpp>

namespace mort
{
    namespace daylight
    {

        struct bondinfo_t
        {
            bondinfo_t();
            
            void apply(morf_t& bond);
            
            void reset( );
            
            int arom;

            int order;

            int stereo;

        };

        struct status_t
        {
            status_t(const hashid_t& lan, molecule_t& mol );
            
            void new_atom( );
            
            void end_atom( );
            
            hashid_t bracket;

            hashid_t language;

            bondinfo_t bond_info;            

            atom_t curt_atom;
            
            atom_t prev_atom;            
            
            vector<morf_t> stack;
            
            vector<atom_t> ring_atoms;

            molecule_t* pmol;
        };
        
    }  // namespace daylight
        
} // namespace mort

#endif
