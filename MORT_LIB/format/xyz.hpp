
#ifndef MORTSRC_FORMAT_XYZ_HPP
#define MORTSRC_FORMAT_XYZ_HPP

namespace mort
{
    class molecule_t;

    void write_amber_xyz( ostream& os, const molecule_t& m );

    void write_tinkerxyz( ostream& os, const molecule_t& m );

} // namespace mort

#endif

