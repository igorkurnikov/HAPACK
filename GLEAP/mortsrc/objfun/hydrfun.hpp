#ifndef MORTSRC_OBJECT_HYDRFUN_HPP
#define MORTSRC_OBJECT_HYDRFUN_HPP

namespace mort
{

    // functions to add and delete hydrogens

    class atom_t;

    class atmvec;

    class molecule_t;

    class database_t;

    void addHs( atom_t& mo );
    
    void addHs( atmvec& avec );
    
    void addHs( molecule_t& mol );
    
    void addHs( database_t& mdb );

    void delHs( molecule_t& m );

}

#endif
