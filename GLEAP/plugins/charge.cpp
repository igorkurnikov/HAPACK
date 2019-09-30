#include <sstream>
#include <object.hpp>
#include "charge.hpp"
#include <sstream>

namespace amber
{
    
    charge_command::charge_command( )
        : command_i( "charge" )
    {
    }
    
    charge_command::charge_command( const string& object )
        : m_object( object )
    {
    }
    
    charge_command::~charge_command( )
    {
    }
    
    bool charge_command::exec( )
    {
        int ndot = std::count( m_object.begin(), m_object.end(), '.' );

        if( ndot < 0 || ndot > 2 )
        {
            throw std::runtime_error( "Error: can not understand object " + m_object );
        }

        double result;

        if( ndot == 0 )
        {
            molecule_ptr pmol = content().get_mol( m_object );
            assert( pmol != NULL );
            result = charge( *pmol );
        }
        else if( ndot == 1 )
        {
            morf_t resd = content().get_resd( m_object );
            result = charge( resd );
        }
        else
        {
            morf_t atom = content().get_atom( m_object );            
            result = atom.get_d(PCHG);
        }
        

        std::ostringstream os;
        
        os << result;
        
        leaplog_t::putline( os.str() );

	return true;
    }

    void charge_command::undo( )
    {
    }
    
    const char* charge_command::info( ) const
    {
        return "  usage: charge object ";
    }
    
    shared_ptr< command_i > charge_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 2 )
        {
            throw std::runtime_error( "Error: wrong number of argument " );
        }
        
        return shared_ptr< command_i >( new charge_command( args[1] ) );
    }
    
} // namespace amber

amber::charge_command g_charge_command;
 
