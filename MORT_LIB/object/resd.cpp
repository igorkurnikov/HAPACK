#include <objfun.hpp>
#include "mole.hpp"
#include "resd.hpp"
#include "atom.hpp"
#include "bond.hpp"

namespace mort
{
    using std::runtime_error;

    resd_t::resd_t( const morf_t& rhs )
        : morf_t( rhs )
    {
        if( rhs.cmpid()!=RESD )
        {
            throw runtime_error( "Error: cannot initialize resd_t from other type of morf" );
        }
    }

    resd_t::resd_t( const molecule_t& m, int id )
        : morf_t( m, RESD, id )
    {
    }

    resd_t resd_t::create( molecule_t& m, int lid )
    {
        assert( lid <= m.nresd() );
        return lid < 0 ? m.create(RESD) : m.insert(RESD, lid);
    }

    resd_t resd_t::get( const molecule_t& m, const string& name )
    {
        resditer_t ri = std::find_if(m.resd_begin(), m.resd_end(), sparm_cmper1(NAME,name) );
        if( ri == m.resd_end() )
        {
            throw runtime_error("Error: cannot find resd " + name + " in molecule.");
        }

        return *ri;
    }

    atom_t resd_t::create_atom( const string& name )
    {
        molecule_t* m = &getmol();

        resditer_t n = iter()+1;
        resditer_t e = m->resd_end();
        while( n != e && n->natom()==0 )
        {
            ++n;
        }

        int pos = (n==e) ? m->natom() : n->atoms()[0].relid();
                    
        atom_t a = m->insert(ATOM, pos );
        a.set_s(NAME, name);

        a.connect( *this );
        connect( a );
        return a;
    }   

    resditer_t resd_t::iter() const
    {
        int lid = relid();
        return getmol().resd_begin() + lid;
    }


} // namespace mort

