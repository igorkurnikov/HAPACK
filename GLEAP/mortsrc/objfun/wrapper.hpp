#ifndef MORTSRC_OBJECT_ACWRAPPER_HPP
#define MORTSRC_OBJECT_ACWRAPPER_HPP


// wrapper functions to external programs (antechamber, parmchk, etc)
namespace mort
{
    class molecule_t;

    // call antechamber to assign GAFF atom type and am1-bcc partical charge.
    void setpchg( molecule_t& m );
    
    // call parmchk (a program of antechamber) to add missing force field parameter.
    void parmchk( molecule_t& m, molecule_t& ff );

    void fixbond( molecule_t& m );
}

#endif


