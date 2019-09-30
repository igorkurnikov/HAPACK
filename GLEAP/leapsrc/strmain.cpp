#include <fstream>
#include <common.hpp>
#include <object.hpp>
#include <guilib.hpp>
#include <cstdio>

#include "strbuff.hpp"

using namespace std;

using namespace mort;

using namespace amber;


int main( int argc, char** argv )
try
{
    const char* amberhome = getenv( "AMBERHOME" );
    if( amberhome == NULL )
    {
        cerr << "Environmental variable AMBERHOME must be set." << std::endl;
        return -1;
    }

    string path = "./:";
    path += amberhome + string( "/dat/leap/gleap:" );
    path += amberhome + string( "/dat/leap/prep:" );
    path += amberhome + string( "/dat/leap/parm:" );
    path += amberhome + string( "/dat/leap/lib:" );
    path += amberhome + string( "/dat/leap/cmd"  );
    mortenv().set_s( "path", path.c_str() );


    ifstream fleaprc;


    if( argc > 1 && argv[1] == string("-f") )
    {
        mortenv().set_s("echo",  "on" );
        mortenv().set_s("batch", "on" );
        assert( argc >= 3 );
        fleaprc.open( argv[2] );
        console() = console_ptr( new std_console(fleaprc, cout) );
        if( !console()->mainloop() )
	{
            return false;
        }
    }

    mortenv().set_s("echo",  "off");
    mortenv().set_s("batch", "off");
    console() = console_ptr( new std_console(cin, cout) );
    return console()->mainloop();
}
catch( std::exception& e )
{
    std::cout << e.what() << std::endl;
    funstack_t::print();
}


