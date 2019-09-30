#include <fstream>
#include <object.hpp>
#include <format.hpp>
#include <ambfmt.hpp>
#include "loadprep.hpp"

using namespace mort;

using namespace boost;

namespace amber
{

    loadprep_command::loadprep_command()
      	:command_i( "loadamberprep" )
    { 
    }

    loadprep_command::loadprep_command(const string& file, const string& prefix)
        :m_file( file ), m_prefix( prefix )
    {
    }

    loadprep_command::~loadprep_command()
    {
    }

    bool loadprep_command::exec()
    {
        std::string file = find_file(m_file);
        std::ifstream is( file.c_str() );

	database_t mdb;
        read_amber_prep( is, mdb );

	database_t::iterator mi = mdb.begin();
	database_t::iterator me = mdb.end();
	for( ; mi != me; ++mi )
	{
	    if( !m_prefix.empty() )
	    {
	    	string name = m_prefix + mi->first;
	        mi->second->set_s(NAME, name);
	        content().set( name, mi->second );
	    }
	    else
	    {
	        content().set( mi->first, mi->second );
	    }
	}

	return true;
    }

    void loadprep_command::undo()
    {
        throw std::runtime_error( "Sorry: not implemented yet." );
    }

    command_ptr loadprep_command::clone(const vector<string>& args) const
    {
	if( args.size()==2 )
	{
            return command_ptr( new loadprep_command(args[1], "") );
	}
	else if( args.size()==3 )
	{
            return command_ptr( new loadprep_command(args[1], args[2]) );
	}
        else
        {
            throw std::runtime_error( "Error: wrong number of arguments" );
        }
    }

} // namespace amber 

amber::loadprep_command g_loadprep_command;


