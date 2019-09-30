#include <stdexcept>
#include <common.hpp>
#include <objfun.hpp>
#include <pdbent.hpp>
#include "bond.hpp"
#include "comp.hpp"
#include "mole.hpp"

namespace mort
{
    using std::runtime_error;

    using boost::dynamic_pointer_cast;

    database_t::database_t()
    {
    }

    database_t::~database_t()
    {
    }

    bool database_t::has( const string& name ) const
    {
        const_iterator i = begin();
        for( ; i != end(); ++i )
        {
            if( i->first == name )
            {
                return true;
            }
        }

        return false;
    }

    entity_ptr database_t::get( const string& name ) const
    {
        const_iterator i = begin();
        for( ; i != end(); ++i )
        {
            if( i->first == name )
            {
                return i->second;
            }
        }
        
        return entity_ptr();
    }

    molecule_ptr database_t::get_mol(const string& name) const
    {
        entity_ptr pe = get( name );
        if( pe != NULL )
        {
            molecule_ptr pm = dynamic_pointer_cast<molecule_t>(pe);
            if( pm != NULL )
            { 
                return pm;
            }
            throw runtime_error( "Error: while getting molecule " + name + " found an entry but it is not a molecule" );
        }

        const_iterator i = begin();
	for( ; i != end(); ++i )
	{
	    database_ptr db = dynamic_pointer_cast<database_t>(i->second);
	    if( db != NULL )
	    {
                molecule_ptr pm = db->get_mol(name);
		if( pm != NULL )
                {
                    return pm;
                }
            }
        }
        throw runtime_error( "Error: cannot find molecule " + name + " in the database" );
    }    

    database_ptr database_t::get_mdb(const string& name) const
    {
        entity_ptr pe = get( name );
	if( pe == NULL )
	{
	    return database_ptr();
        }

	database_ptr db = dynamic_pointer_cast< database_t >( pe );
	if( db == NULL )
        {
            throw runtime_error( "Error: there is such an entry in db, but it is not a database" );
        }

        return db;
    }

    namemap_ptr database_t::get_nmap(const string& name) const
    {
        entity_ptr pe = get( name );
        if( pe == NULL )
        {
            return namemap_ptr();
        }

	namemap_ptr nmap = dynamic_pointer_cast< namemap_t >( pe );

	if( nmap == NULL )
        {
            throw runtime_error( "Error: there is such an entry in db, but it is not a namemap" );
        }

        return nmap;
    }

    atmvec_ptr database_t::get_avec(const string& name) const
    {
        entity_ptr pe = get( name );
        if( pe == NULL )
        {
            return atmvec_ptr();
        }
        
        atmvec_ptr avec = dynamic_pointer_cast<atmvec>(pe);
        if( avec == NULL )
        {
            throw runtime_error( "Error: there is such an entry in db, but it is not a atmvec" );
        }

        return avec;
    }

    bndvec_ptr database_t::get_bvec(const string& name) const
    {
        entity_ptr pe = get( name );
        if( pe == NULL )
        {
            return bndvec_ptr();
        }
        
        bndvec_ptr bvec = dynamic_pointer_cast<bndvec>(pe);
        if( bvec == NULL )
        {
            throw runtime_error( "Error: there is such an entry in db, but it is not a bndvec" );
        }

        return bvec;
    }        

    morf_t database_t::get_atom( const string& mask ) const
    {
        int dot = mask.rfind( '.' );
        
        if( dot == (int)string::npos )
        {
            throw runtime_error( "Error: not a good atom mask: " + mask );
        }
        
        resd_t r = get_resd( mask.substr( 0, dot ) );

        string name = mask.substr(dot+1, mask.length()-dot-1 );
        
        return atom_t::get(r, name);
    }

    morf_t database_t::get_resd( const string& mask ) const
    {
        int dot = mask.find( '.' );
        
        if( dot == (int)string::npos )
        {
            throw runtime_error( "Error: not a good atom mask: " + mask );
        }
        
        string molname = mask.substr( 0, dot );

        molecule_ptr pmol = get_mol( molname );

        if( pmol == NULL )
        {
            throw runtime_error("Error: can not find molecule " + molname + " in database" );
        }

        string r = mask.substr( dot + 1, mask.length() - dot - 1 );
        assert( r.length()>0 );

	if( isdigit(r[0]) )
	{
	    int rid = atoi( r.c_str() );
	    assert( rid > 0 );
            return resd_t( *pmol, rid-1 );
	}

	return resd_t::get( *pmol, r );
    }    

    void database_t::set( const string& name, const entity_ptr& pe )
    { 
        iterator i = begin();
        for( ; i != end(); ++i )
        {
            if( i->first == name )
            {
                i->second = pe;
            }
        }

        m_data.push_back( make_pair(name, pe) );
    }
 
    void database_t::add( const string& name, const entity_ptr& pe )
    {
        m_data.push_back( make_pair(name, pe) );
    }

    bool database_t::remove( const string& name )
    {
        iterator i = begin();
        for( ; i != end(); ++i )
        {
            if( i->first == name )
            {
                m_data.erase(i);
                return true;
            }
        }
            
        return false;
    }

    database_t& mortenv()
    {
        static database_t g_mortenv;
        return g_mortenv;
    }




} // namespace mort


