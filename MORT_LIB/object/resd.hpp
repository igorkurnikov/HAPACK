#ifndef MORTSRC_OBJECT_RESD_HPP
#define MORTSRC_OBJECT_RESD_HPP

#include "morf.hpp"

namespace mort
{

    class resd_t : public morf_t
    {
    public:

        resd_t( const morf_t& rhs );

        resd_t( const mole_t& m, int id );

        resditer_t iter() const;

        atom_t create_atom( const string& name );

        static resd_t create( molecule_t& m, int lid=-1 );

        static resd_t get( const molecule_t& m, const string& name );

        string name() const  { return get_s(NAME); }

        int chain() const { return get_i(CHAIN); }

        int head() const { return get_i(HEAD); }

        int tail() const { return get_i(TAIL); }
    };


} // namespace mort

#endif

