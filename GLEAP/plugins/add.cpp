#include <boost/algorithm/string.hpp>
#include <object.hpp>
#include "add.hpp"

namespace amber
{
    using namespace mort;
    using namespace boost;


    atom_t add_atom( resd_t& r, const root_t& a )
    {
	atom_t atom = r.create_atom( "" );
        if( a.has_s(NAME) ) atom.set_s(NAME, a.get_s(NAME) );
        if( a.has_s(TYPE) ) atom.set_s(TYPE, a.get_s(TYPE) );
        if( a.has_d(PCHG) ) atom.set_d(PCHG, a.get_d(PCHG) );
        if( a.has_i(ELEMENT) ) atom.set_i(ELEMENT, a.get_i(ELEMENT) );
        if( a.has_v(POSITION) ) atom.set_v(POSITION, a.get_v(POSITION) );
	return atom;
    }


    void add_resd( molecule_t& m, database_t& r )
    {
        resd_t resd = resd_t::create( m );

        if( r.has_s(NAME) ) resd.set_s(NAME, r.get_s(NAME) );

        database_t::iterator i = r.begin();
        for( ; i != r.end(); ++i )
        {
            entity_ptr pa = i->second;
            atom_t atom = add_atom( resd, *pa );
	    
            if( pa->has_s(NBRLIST) )
            {
                assert( pa->has_s(NBRTYPE) );
                
                vector< string > names, types;
                string line = pa->get_s(NBRLIST);
                split( names, line, is_any_of(" "), token_compress_on);
		line = pa->get_s(NBRTYPE);
                split( types, line, is_any_of(" "), token_compress_on);
                
                for( unsigned int i=0; i < names.size(); ++i )
                {
		    atom_t nbr(atom.getmol(), -1);
                    if( atom_t::get( resd, names[i], nbr) )
                    {
                        if( !bond_t::has( atom, nbr ) )
                        {
                            bond_t::create(atom, nbr).set_i(ORDER, atoi(types[i].c_str()) );
                        }
                    }
                }
            }
        }
    }





    add_command::add_command( )
        :command_i( "add" )
    {
    }
    
    add_command::add_command( const string& dst, const string& src )
        : m_dst( dst ), m_src( src )
    {
    }
    
    add_command::~add_command( )
    {
    }
    
    bool add_command::exec( )
    {
        int ndot = count( m_dst.begin(), m_dst.end(), '.' );

	if( ndot==0 )
	{
            entity_ptr pdst = content().get( m_dst );
	    entity_ptr psrc = content().get( m_src );

            // if dst is a database (a pesudo-residue), simple add src (should be a pesudo-atom) to it.
	    database_ptr presdst = dynamic_pointer_cast< database_t >( pdst );
	    if( presdst )
	    {
	        presdst->add( psrc->get_s(NAME), psrc );
	        return true;
	    }
        
            // if not a database, must be a molecule (sometimes called unit)
	    molecule_ptr pmoldst = dynamic_pointer_cast< molecule_t >( pdst );
	    if( pmoldst==NULL )
	    {
	        throw std::runtime_error( "Error: " + m_dst + " is neither a molecule or a reside, can not add things into it" );
	    }

	    database_ptr pressrc = dynamic_pointer_cast< database_t >( psrc );
	    if( pressrc==NULL )
	    {
	        throw std::runtime_error( "Error: " + m_src + " is not a residue. Only residues are allowed to be add to molecule." );
	    }

	    add_resd( *pmoldst, *pressrc );
	    return true;
        }

        if( ndot==2 )
	{
	    throw std::runtime_error( "Error: cannot add things to atom " + m_dst );
	}

	assert( ndot==1 );

        entity_ptr psrc = content().get( m_src );
	resd_t resdst = content().get_resd( m_dst );
	add_atom( resdst, *psrc );

        return true;
    }
    
    void add_command::undo( )
    {
        throw std::runtime_error( "sorry, not implemented yet" );
    }
    
    const char* add_command::info( ) const
    {
        return "  add list atom\n   or \n  add unit list";
    }
    
    shared_ptr< command_i > add_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 3 )
        {
            throw std::runtime_error( "wrong number of argument, try help add\n" );
        }
        
        return shared_ptr< command_i >( new add_command( args[1], args[2] ) );
    }

    null_command::null_command( )
    {
    }
    
    null_command::null_command( const string& name )
        : command_i( name )
    {
    }

    null_command::~null_command()
    {
    }
    
    bool null_command::exec( )
    {
        return true;
    }
    
    void null_command::undo( )
    {
    }
    
    const char* null_command::info( ) const
    {
        return " Warning: this command is depreted and is no longer supported, it is kept in system just for compatiability";
    }
    
    shared_ptr< command_i > null_command::clone( const vector< string >& ) const
    {
        return shared_ptr< command_i >( new null_command() );
    }
    
}

amber::add_command g_add_command;

amber::null_command g_addatomtypes_command( "addatomtypes" );

amber::null_command g_logfile_command( "logfile" );
