#include <fstream>
#include <boost/algorithm/string.hpp>
#include <object.hpp>
#include "mainwin.hpp"
#include "leaplog.hpp"

namespace mort
{
    shared_ptr< console_t >& console()
    {
        static shared_ptr< console_t > g_console;
        return g_console;
    }

    database_t& content()
    {
        return mortenv();
    }

    string find_file( const string& name )
    {
        std::ifstream is( name.c_str() );
	if( is )
	{
	    return name;
	}


        string path;
        if( !mortenv().get_s("path", path) )
        {
            throw std::runtime_error( "Error: mort path not set" );
        }

        vector<string> dirs;
        split( dirs, path, is_any_of(":") );

        // windows platform need special handling,
        // on windows filename is sometimes c:\...
        //
        vector<string> newdirs;
        unsigned int i=0;
        while( i < dirs.size() )
        {
            if( dirs[i].length()==1 )
            {
                // disk label
                newdirs.push_back( dirs[i] + ":\\" + dirs[i+1] );
                i += 2;
            }
            else
            {
                newdirs.push_back( dirs[i] );
                i += 1;
            }
        }

        dirs.swap( newdirs );
  

        
        for( int i=0; i<(int)dirs.size(); ++i )
        {
            string full_name = dirs[i] + "/" + name;

            std::ifstream is( full_name.c_str() );
            
            if(is) 
            {
                return full_name;
            }
        }
        throw std::runtime_error( "Error: can not find file " + name + " in all the search path." );
    }
    

} // namespace mort

