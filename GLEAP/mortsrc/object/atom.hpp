#ifndef MORTSRC_OBJECT_ATOM_HPP
#define MORTSRC_OBJECT_ATOM_HPP

#include "morf.hpp"

namespace mort
{

    class atom_t : public morf_t
    {
    public:

        atom_t( const morf_t& rhs );

        atom_t( const molecule_t& m, int id );

        virtual ~atom_t() {}

        resd_t resd() const;

        numvec pos()  const { return get_v(POSITION); }

        string name() const { return get_s(NAME); }

        string type() const { return get_s(TYPE); }

        int element() const { return get_i(ELEMENT); }

        static atom_t create( molecule_t& m, const string& name );

        static atom_t create( resd_t& r, const string& name );

        static atom_t get( const mole_t& m, const string& name );

        static atom_t get( const resd_t& r, const string& name );

        static atom_t get( const resd_t& r, const hashid_t& pid, int pval );

        static bool get( const mole_t& m, const string& name, atom_t& a );

        static bool has( const mole_t& m, const string& name );

        static bool get( const resd_t& r, const string& name, atom_t& a );

        static bool has( const resd_t& r, const string& name );


    };

} // namespace mort

#endif
