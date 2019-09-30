#include <object.hpp>
#include "copy.hpp"

namespace amber
{
    
    copy_command::copy_command( )
        : command_i( "copy" )
    {
    }
    
    copy_command::copy_command( const string& dst, const string& src )
        : m_src(src), m_dst(dst) 
    {
    }
    
    copy_command::~copy_command( )
    {
    }
    
    bool copy_command::exec( )
    {
        int ndot = std::count( m_src.begin(), m_src.end(), '.' );
	if( ndot != 0 )
	{
	    throw std::runtime_error( "Error: cannot copy atom or resd in a molecule." );
	}

        entity_ptr psrc = content().get( m_src );

        database_ptr pmdb = dynamic_pointer_cast< database_t >( psrc );
        if( pmdb != NULL )
        {
            database_ptr pdst( new database_t(*pmdb) );
            content().set( m_dst, pdst );
            return true;
        }
        
        molecule_ptr pmol = dynamic_pointer_cast< molecule_t >( psrc );        
        if( pmol != NULL )
        {
            molecule_ptr pdst( new molecule_t(*pmol) );
            content().set( m_dst, pdst );
            return true;
        }

        entity_ptr pdst( new root_t(*psrc) );
        content().set( m_dst, pdst );
        return true;
    }
    
    void copy_command::undo( )
    {
        throw std::runtime_error( "Error: sorry, not implemented yet" );
    }
    
    const char* copy_command::info( ) const
    {
        return "  usage: newvar = copy var ";
    }
    
    shared_ptr< command_i > copy_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 3 )
        {
            throw std::runtime_error( "Error: wrong number of arguments" );
        }

        return shared_ptr< command_i >( new copy_command( args[1], args[2] ) );
    }
    
} // namespace amber
  
amber::copy_command g_copy_command;

