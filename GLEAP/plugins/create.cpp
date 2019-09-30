#include <object.hpp>
#include "create.hpp"

namespace amber
{
    create_command::create_command( const string& object )
        : command_i( "create" + object )
    {
    }
    
    create_command::create_command( const string& object, const string& dest, const string& name, const string& type, const string& charge )
        : m_object( object ), m_dest( dest ), m_name( name ), m_type( type ), m_chrg( charge )
    {
        if( m_object != "atom" )
        {
            throw std::runtime_error( "Error: only createAtom have three argument" );
        }
    }

    create_command::create_command( const string& object, const string& dest, const string& name )
        : m_object( object ), m_dest( dest ), m_name( name )
    {
        if( m_object == "atom" )
        {
            throw std::runtime_error( "Error: createAtom should have three argument" );
        }
    }

    create_command::~create_command()
    {
    }
    
    bool create_command::exec()
    {
        entity_ptr ptr;
        
        if( m_object == "atom" )
        {
            ptr = entity_ptr( new root_t() );
            ptr->set_s(NAME, m_name);
            ptr->set_s(TYPE, m_type);
            ptr->set_d(PCHG, atof( m_chrg.c_str() ) );
        }
        else if( m_object == "residue" )
        {
            ptr = database_ptr( new database_t() );
            ptr->set_s(NAME, m_name);
        }
        else
        {
            assert( m_object == "unit" );
            ptr = molecule_ptr( new molecule_t() );
            ptr->set_s(NAME, m_name);
        }
        
        content().set( m_dest, ptr );
	return true;
    }
    
    void create_command::undo( )
    {
        // content().remove( m_name );
        throw std::runtime_error( "Error: not implemented yet." );
    }
    
    const char* create_command::info( ) const
    {
        return "  a = createAtom name type charge or\n  r = createResidue name or\n  u = createUnit name\n";
    }

    shared_ptr< command_i > create_command::clone( const vector< string >& args ) const
    {
        string object = args[0].substr( 6, args[0].length() - 6 );
        
        if( object == "atom" )
        {
            if( args.size() != 5 )
            {
                throw std::runtime_error( "Error: wrong number of arguments" );
            }

            return shared_ptr< command_i >( new create_command( object, args[1], args[2], args[3], args[4] ) );
        }

        if( args.size() != 3 )
        {
            throw std::runtime_error( "Error: wrong number of arguments" );
        }

        return shared_ptr< command_i >( new create_command( object, args[1], args[2] ) );
    }
    
} // namespace amber

amber::create_command g_createAtom_command( "atom" );
amber::create_command g_createResd_command(  "residue" );
amber::create_command g_createMolecule_command( "unit" );


