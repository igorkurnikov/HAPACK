#ifndef MORTSRC_AMBFMT_AMBOFF_HPP
#define MORTSRC_AMBFMT_AMBOFF_HPP

#include <iosfwd>

namespace mort
{
    class molecule_t;

    bool read_off( std::istream& is, molecule_t& m );

    void write_off(std::ostream& os, molecule_t& m );

}


#endif
