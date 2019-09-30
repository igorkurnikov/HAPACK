#ifndef MORTSRC_FORMAT_TOP_HPP
#define MORTSRC_FORMAT_TOP_HPP

#include <iosfwd>
#include <string>

namespace mort
{
    class molecule_t;



    void load_frc( const std::string& file, molecule_t& ffp);
 
    void save_top( molecule_t& m, const std::string& file, const molecule_t& ffp );
 
}

#endif
