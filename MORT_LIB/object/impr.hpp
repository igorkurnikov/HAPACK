#ifndef MORTSRC_OBJECT_IMPR_HPP
#define MORTSRC_OBJECT_IMPR_HPP

#include "morf.hpp"

namespace mort
{
    class molecule_t;

    class impr_t : public morf_t
    {
    public:

        impr_t( const morf_t& rhs );

        impr_t( const molecule_t& m, int id );

        static impr_t create( atom_t& a1, atom_t& a2, atom_t& a3, atom_t& a4, int p );

        static impr_t frcget( atom_t& a1, atom_t& a2, atom_t& a3, atom_t& a4, int p );

        static impr_t get( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4 );

        static bool get( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4, impr_t& im );

        static bool has( const atom_t& a1, const atom_t& a2, const atom_t& a3, const atom_t& a4 );
    };

} // namespace mort

#endif

