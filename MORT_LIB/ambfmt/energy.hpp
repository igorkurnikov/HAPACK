#ifndef GLEAP_MORTSRC_AMBFMT_HPP
#define GLEAP_MORTSRC_AMBFMT_HPP


namespace mort
{
    double eval_bond( const molecule_t& mol );

    double eval_angl( const molecule_t& mol );

    double eval_oops( const molecule_t& mol );
    
    vector3d eval_tors( const molecule_t& mol );
    
    vector2d eval_nonbond( const molecule_t& mol, double cut );

    vector2d eval_direct( const molecule_t& mol, double cut );

    double eval_regewald( const molecule_t& mol, double cut );


}

#endif
