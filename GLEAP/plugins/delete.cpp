#include <object.hpp>
#include "delete.hpp"

namespace amber
{
    
    delete_command::delete_command( const string& action )
        : command_i( action )
    {
    }
    
    delete_command::delete_command( const string& action, const string& oper1, const string& oper2 )
        : m_action( action ), m_oper1( oper1 ), m_oper2( oper2 )
    {
    }
    
    delete_command::~delete_command()
    {
    }
    
    bool delete_command::exec( )
    {
        if( m_action=="remove" )
        {
	    int ndot = count( m_oper2.begin(), m_oper2.end(), '.' );

	    if( ndot==1 )
	    {   
	        morf_t r = content().get_resd( m_oper2 );
		molecule_t* pm = &(r.getmol());
		pm->remove_resd( r );
		pm->cleanup();
	    }
	    else
	    {
	        assert( ndot==2 );
                morf_t a = content().get_atom( m_oper2 );
                molecule_t* pm = &(a.getmol());
		pm->remove_atom( a );
		pm->cleanup();
	    }
        }
        else  
	{
            assert( m_action == "deletebond" );

            morf_t a1 = content().get_atom( m_oper1 );
            morf_t a2 = content().get_atom( m_oper2 );
            morf_t b0 = bond_t::get(a1, a2);

	    molecule_t* pm = &(b0.getmol());
            pm->remove( b0 );
	    pm->cleanup();
        }

        return true;
    }
    
    void delete_command::undo( )
    {
        throw std::runtime_error( "Error: sorry, not implemented." );
    }
    
    const char* delete_command::info( ) const
    {
        return "  usage: deleteBond atom1 atom2 | remove unit atom|resd";
    }
    
    command_ptr delete_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 3 )
        {
            throw std::runtime_error( "Error: wrong number of arguments" );
        }
        
        return command_ptr( new delete_command( args[0], args[1], args[2] ) );
    }
    
} // namespace amber

amber::delete_command g_deleteAtom_command("deleteAtom");
amber::delete_command g_deleteBond_command("deleteBond");
amber::delete_command g_deleteResd_command("deleteResd");


        
