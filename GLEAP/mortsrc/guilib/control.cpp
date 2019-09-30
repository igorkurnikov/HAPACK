#include <sstream>
#include <boost/algorithm/string.hpp>
#include <object.hpp>

#include "grammar.hpp"
#include "command.hpp"
#include "control.hpp"
#include "leaplog.hpp"
#include "mainwin.hpp"

namespace mort
{
    using std::map;

    static vector<command_ptr>& history()
    {
        static vector< shared_ptr< command_i > > g_history;
        return g_history;
    }

    static map<string, command_i*> g_commands;

    control_t::control_t()
    {
    }

    control_t::~control_t()
    {
    }

    void control_t::insert( const string& name, command_i* ptr )
    {
        g_commands[name] = ptr;
    }

    entity_ptr recu_find( const database_t& mdb, const string& name )
    {
        database_t::const_iterator i = mdb.begin();
        for( ; i != mdb.end(); ++i )
        {
            if( i->first==name )
                return i->second;
        }
        
        database_t::const_iterator j = mdb.begin();
        for( ; j != mdb.end(); ++j )
        {
            database_ptr pdb = dynamic_pointer_cast< database_t >(j->second);
            if( pdb != NULL )
            {
                entity_ptr p = recu_find( *pdb, name );
                if( p != NULL )
                    return p;
            }
        }
        
        return entity_ptr();
    }



    void make_alias( const string& src, const string& dst )
    {
        entity_ptr p = recu_find( content(), src );
        
        if( p==NULL )
        {
            throw std::runtime_error( "Error: can not make alias " + dst + "->" + src + ", " + src + " does not exist." );
        }
        
        database_ptr palias = content().get_mdb( "_alias" );
        if( palias==NULL )
        {
            palias = database_ptr( new database_t() );
            content().set( "_alias", palias );
        }
        
        palias->set( dst, p );
    }
    
    bool control_t::run( const string& command )
    {

        assert( ! command.empty() );

        vector<string> args;

        try
        {
            interpret(command, args);

            if( command.find('=') != string::npos && args.size() == 2 )
            {
                if( ndim(args[0])==0 )
                {
                    make_alias( args[0], args[1] );                    
                }
                else
                {
                    vector<string> mols = interpret1d(args[0]);
                    database_ptr pdb = database_ptr( new database_t() );
                    for(int i=0; i < (int)mols.size(); ++i)
                    {
                        molecule_ptr pm = content().get_mol(mols[i]);
                        pdb->add( mols[i], pm );
                    }
                    content().set(args[1], pdb);
                }
            }
            else
            {
                args[0] = to_lower_copy(args[0]);

                string name = args[0];

                command_i* pcmd = command_i::find( name );

                assert( pcmd != NULL );

                command_ptr comm = pcmd->clone(args);
                
                comm->exec();

                history().push_back( comm );
            }

            return true;
        }
        catch( std::exception& e )
        {
            std::cout <<  e.what() << std::endl;
            funstack_t::print();
            return false;
        }

        return false;
    }

    map< string, command_i* >::iterator control_t::begin()
    {
        return g_commands.begin();
    }

    map< string, command_i* >::iterator control_t::end()
    {
        return g_commands.end();
    }

} // namespace mort



