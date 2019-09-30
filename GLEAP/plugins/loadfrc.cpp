#include <fstream>
#include <object.hpp>
#include <format.hpp>
#include <ambfmt.hpp>
#include "loadfrc.hpp"

namespace amber
{
    using std::ifstream;

    loadfrc_command::loadfrc_command( const string& name )
	: command_i( name )
    {
    }
    
    loadfrc_command::loadfrc_command( const string& action, const string& file )
        : m_file( file ), m_action( action ) 
    {
    }
    
    loadfrc_command::~loadfrc_command()
    {
    }
    
    bool loadfrc_command::exec()
    {
        ifstream stream( find_file( m_file ).c_str() );
        
        if( m_action == "loadamberparams" )
        {
            shared_ptr< molecule_t > amberff;

            if( content().has( "_amberffp" ) )
            {
                amberff = dynamic_pointer_cast< molecule_t >( content().get( "_amberffp" ) );
            }
            else
            {
                amberff = shared_ptr< molecule_t >( new molecule_t() );
            }
            
            read_frc( stream, *amberff );
            
            if( ! content().has( "_amberffp" ) )
            {
                amberff->set_s(NAME, "_amberffp");
                content().set( "_amberffp", amberff );
            }
        }
        else if( m_action == "loadamoebaparams" )
        {
            shared_ptr< molecule_t > atomff( new molecule_t() );
            shared_ptr< molecule_t > poleff( new molecule_t() );
            
            read_amoeba_frc( stream, *atomff, *poleff );
            
            content().set( "_amoeba-atom", atomff );
            content().set( "_amoeba-pole", poleff );
        }
        else
        {
            assert( false );
        }

        return true;
    }
    
    void loadfrc_command::undo()
    {
    }
    
    shared_ptr< command_i > loadfrc_command::clone( const vector< string >& args ) const
    {
        if( args.size() < 2 || args.size() > 3 )
        {
            throw std::runtime_error( "Error: wrong number of arguments" );    
        }

        return shared_ptr< command_i >( new loadfrc_command( args.front(), args.back() ) );
    }

} // namespace amber

amber::loadfrc_command g_loadamberparams_command( "loadamberparams" );

amber::loadfrc_command g_loadamoebaparms_command( "loadamoebaparams" );
