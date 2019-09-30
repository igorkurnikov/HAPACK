#include <object.hpp>

namespace mort
{
    class molecule_t;
    
    namespace prmtop
    {        

        void atom2angl( const molecule_t* ff, const atom_t& ac, bool fixed_water_model ); 

        void bond2dihe( const molecule_t* ff, const bond_t& bc );

        void atom2impr( const molecule_t* ff, const atom_t& ac );

        void angl2tor2( const molecule_t* ff, const morf_t& angl );

        void bond2ptor( const morf_t& bond );
        
    } // namespace prmtop

} // namespace mort


