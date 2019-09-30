#ifndef MORTSRC_OBJECT_DIHE_HPP
#define MORTSRC_OBJECT_DIHE_HPP

#include <vector>
#include "morf.hpp"

namespace mort
{
    class dihe_t;

    typedef std::vector<dihe_t> dihvec;

    class dihe_t : public morf_t
    {
    public:

        dihe_t( const morf_t& rhs );

        dihe_t( const molecule_t& m, int id );

        int period() const;

        static dihe_t create( atom_t& a1, atom_t& a2, atom_t& a3, atom_t& a4, int p );

        static dihe_t frcget( atom_t& a1, atom_t& a2, atom_t& a3, atom_t& a4, int p );

        static dihe_t get( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4, int p );

        static dihvec get( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4 );

        static bool get( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4, int p, dihe_t& d );

        static bool has( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4, int p );

        static bool has( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4 );

    };

} // namespace mort

#endif
