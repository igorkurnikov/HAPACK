#include <stdexcept>
#include <objfun.hpp>
#include "mole.hpp"
#include "atom.hpp"

namespace mort
{
    using std::runtime_error;

    atom_t::atom_t( const morf_t& rhs )
        : morf_t( rhs )
    {
        int cid = rhs.cmpid();
        if( cid!=ATOM ) 
        {
            throw runtime_error( "Error: try to initialize atom with a " + unhash(cid) );
        }
    }

    atom_t::atom_t( const molecule_t& m, int id )
        : morf_t( m, ATOM, id )
    {
    }

    resd_t atom_t::resd() const
    {
        if( nresd()==0 )
        {
            throw runtime_error( "Error: no residue is associated with this atom." );
        }

        return *resd_begin();
    }

    atom_t atom_t::create( mole_t& m, const string& name )
    {
        if( m.nresd()==0 )
        {
            m.create(RESD);
        }

        resd_t r = m.resds()[m.nresd()-1];
        return create( r, name );
    }

    atom_t atom_t::create( resd_t& r, const string& name )
    {
        return r.create_atom(name);
    }

    atom_t atom_t::get( const molecule_t& m, const string& name )
    {
        atom_t a(m, -1);
        if( get( m, name, a) )
        {
            return a;
        }

        throw runtime_error( "Error: cannot find atom " + name + " in molecule." );
    }

    bool atom_t::get( const molecule_t& m, const string& name, atom_t& a )
    {
        atomiter_t ai = std::find_if( m.atom_begin(), m.atom_end(), sparm_cmper1(NAME,name) );
        if( ai == m.atom_end() )
        {
            return false;
        }

        a = *ai;
        return true;
    }   

    bool atom_t::has( const molecule_t& m, const string& name )
    {
        atom_t a(m,-1);
        return get(m, name, a);
    }
 
    atom_t atom_t::get( const resd_t& r, const string& name )
    {
        atom_t a(r.getmol(), -1);
        if( get( r, name, a) )
        {
            return a;
        }

        throw runtime_error( "Error: cannot find atom " + name + " in residue." );
    }

    bool atom_t::get( const resd_t& r, const string& name, atom_t& a )
    {
        atomiter_t ai = std::find_if( r.atom_begin(), r.atom_end(), sparm_cmper1(NAME,name) );
        if( ai == r.atom_end() )
        {
            return false;
        }

        a = *ai;
        return true;
    }   


} // namespace mort

