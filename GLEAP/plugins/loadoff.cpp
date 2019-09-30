#include <fstream>
#include <object.hpp>
#include <format.hpp>

#include "loadoff.hpp"

using namespace mort;

using namespace boost;

namespace amber
{

    loadoff_command::loadoff_command( )
      	:command_i( "loadoff" )
    { 
    }

    loadoff_command::loadoff_command( const string& file, const string& name )
        :m_file( file ), m_name(name)
    {
    }

    loadoff_command::~loadoff_command()
    {
    }

    bool loadoff_command::exec()
    {
        std::ifstream stream( find_file( m_file ).c_str() );

        database_t* pmdb = NULL;
        if( m_name.empty() )
        {
            pmdb = &content();
        }
        else
        {
            database_ptr tmp = content().get_mdb(m_name);
            if( tmp== NULL )
            {
                tmp = database_ptr( new database_t() );
                content().set( m_name, tmp );
            }
            
            assert( tmp !=NULL );
            pmdb = tmp.get();
        }
    
        while( stream )
        {
            shared_ptr< molecule_t > pmol( new molecule_t() );
        
            read_off( stream, *pmol );

            m_names.push_back( pmol->get_s(NAME) );
        
            pmdb->set( pmol->get_s(NAME), pmol );
        }    

	return true;
    }

    void loadoff_command::undo()
    {
        throw std::runtime_error( "Sorry: not implemented yet." );
    }

    shared_ptr< command_i > loadoff_command::clone( const vector< string >& args ) const
    {
        if( args.size()!= 2 && args.size()!=3 )
        {
            throw std::runtime_error( "Error: wrong number of arguments" );
        }

        string name = (args.size()==2) ? "": args[2];
        return shared_ptr< command_i >( new loadoff_command(args[1], name) );
    }

} // namespace amber 

amber::loadoff_command g_loadoff_command;


