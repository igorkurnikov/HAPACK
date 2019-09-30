#ifndef MORTSRC_OBJECT_MERGE_HPP
#define MORTSRC_OBJECT_MERGE_HPP



namespace mort
{
    class resd_t;

    class molecule_t;
 
    resd_t merge( molecule_t& m, const molecule_t& s, int rid=-1 );

    resd_t merge( molecule_t& m, const resd_t& s,     int rid=-1 );
    
    molecule_t unmerge( const resd_t& r );

}

#endif

