#ifndef MORTSRC_OBJECT_ANGL_HPP
#define MORTSRC_OBJECT_ANGL_HPP

#include "morf.hpp"

namespace mort
{

    class angl_t : public morf_t
    {
    public:

        angl_t( const morf_t& rhs );

        angl_t( const molecule_t& m, int id );

        static angl_t create( atom_t& a1, atom_t& a2, atom_t& a3 );

        static angl_t frcget( atom_t& a1, atom_t& a2, atom_t& a3 );

        static angl_t get( const atom_t& a1, const atom_t& a2, const atom_t& a3 );

        static bool has( const atom_t& a1, const atom_t& a2, const atom_t& a3 );

        static bool get( const atom_t& a1, const atom_t& a2, const atom_t& a3, angl_t& an );

        static bool has( const atom_t& a1, const atom_t& a3 );
    };

} // namespace mort

#endif
