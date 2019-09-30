#include <sstream>
#include <object.hpp>
#include "center.hpp"
#include <sstream>

namespace amber
{
    using namespace std;
    
    center_command::center_command( )
        : command_i( "center" )
    {
    }
    
    center_command::center_command( const string& object )
        : m_object( object )
    {
    }
    
    center_command::~center_command( )
    {
    }
    
    bool center_command::exec( )
    {
        int ndot = std::count( m_object.begin(), m_object.end(), '.' );

        if( ndot < 0 ||  ndot > 2 )
        {
            throw logic_error( "Error: can not understand object " + m_object );
        }

        numvec result(3);
        if( ndot == 0 )
        {
            molecule_ptr pmol = content().get_mol(m_object);
            if( pmol == NULL )
            {
                throw logic_error( "Error: can not find molecule " + m_object + " in database" );
            }
            
            result = center( *pmol );
        }
        else if( ndot == 1 )
        {
            morf_t resd = content().get_resd( m_object );
            result = center( resd );
        }
        else
        {
            morf_t atom = content().get_atom( m_object );
            result = atom.get_v(POSITION);
        }

        std::ostringstream os;
        os << string( "The center is at: " );
        os << result[0] << " " << result[1] << " " << result[2];
        leaplog_t::putline( os.str() );
	return true;
    }

    void center_command::undo( )
    {
    }
    
    const char* center_command::info( ) const
    {
        return "  usage: center object ";
    }
    
    shared_ptr< command_i > center_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 2 )
        {
            throw logic_error( "Error: wrong number of argument " );
        }
        
        return shared_ptr< command_i >( new center_command( args[1] ) );
    }
    
} // namespace amber

amber::center_command g_center_command;
 
