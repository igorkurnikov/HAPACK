#include <fstream>
#include <object.hpp>
#include <format.hpp>
#include "loadcrd.hpp"

namespace amber
{
    void read_amoeba_crd( istream& is, molecule_t& mol )
    {
        string version;
        std::getline( is, version );

        // read title
        string flag, format, content;
        std::getline( is, flag );
        std::getline( is, format );
        std::getline( is, content );

        // read simulation time
        std::getline( is, flag );
        std::getline( is, format );
        std::getline( is, content );

        // read atom number
        std::getline( is, flag );
        std::getline( is, format );
        std::getline( is, content );
        int natom = atoi( content.c_str() );

        if( natom != mol.natom() )
        {
            throw std::runtime_error( "Error: wrong number of atom in the coordinate" );
        }

        string comment;
        std::getline( is, flag );
        std::getline( is, comment );
        std::getline( is, format );
 
        atomiter_t atom = mol.atom_begin();
        for( ; atom != mol.atom_end(); ++atom )
        {
            numvec pos(3);
            is >> pos[0] >> pos[1] >> pos[2];
            atom->set_v(POSITION, pos);
        }
    }

    loadcrd_command::loadcrd_command( const string& action )
        :command_i( action )
    { 
    }

    loadcrd_command::loadcrd_command( const string& action, const string& name, const string& file )
        :m_action( action ), m_name( name ), m_file( file )
    {
    }

    loadcrd_command::~loadcrd_command()
    {
    }

    bool loadcrd_command::exec( )
    {
        std::ifstream is( m_file.c_str() );

        if( ! is )
        {
            throw std::runtime_error( "Error: file " + m_file + " doesn't exist" );
        }

        if( m_action == "loadamoebacrd" )
        {
            molecule_ptr pmol = content().get_mol( m_name );
            assert( pmol != NULL );
            read_amoeba_crd( is, *pmol );
        }

        return true;        
    }

    void loadcrd_command::undo()
    {
    }

    shared_ptr< command_i > loadcrd_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 3 )
        {
            throw std::logic_error( "Error: wrong arguments " + args[1] + " " + args[2] );
        }
    
        return shared_ptr< command_i >( new loadcrd_command( args[0], args[1], args[2] ) );
    }

}

static amber::loadcrd_command g_loadamoebacrd_command( "loadamoebacrd" );

