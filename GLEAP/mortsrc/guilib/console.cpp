#include <object.hpp>
#include "console.hpp"
#include "control.hpp"
#include "leaplog.hpp"
#include "mainwin.hpp"
#include "grammar.hpp"

namespace mort
{

    static const char PROMPT[] = "[gtkleap]$ ";

    console_t::console_t( )
    {
    }

    console_t::~console_t()
    {
    }

    bool console_t::mainloop()
    {
        string line;

	while( getline(line, PROMPT) )
	{
	    if( !process(line) )
	    {
	        return false;
	    }
        }

	return true;
    }

    bool console_t::process( const string& line )
    {
        if( empty(line) )
        {
            return true;
        }

        string echo;
        if( mortenv().get_s("echo",echo)&&echo=="on" )
        {
            print( line + "\n" );
        }

        if( line.substr(0,4) == "quit" || line.substr(0,4) == "exit" )
        {
            return false;
        }

        if( line[0] == '#' )
        {
            return true;
        }

        m_pending += line + " ";

        if( !curly_closed(m_pending) )
        {
            return true;
        }

	string cmd = m_pending;

        m_pending = "";

        bool status = control_t::run(cmd);
        vector< string >::iterator i = leaplog_t::cur();
        for( ; i != leaplog_t::end(); ++i )
        {
            print( *i + "\n" );
        }

        leaplog_t::move_cursor( leaplog_t::end() - leaplog_t::cur() );

        string batch;
        if( !status && mortenv().get_s("batch", batch) && batch=="on"  )
        {
            return false;
        }
        return true;
    }

/*
    char std_console::getchr( const char* choise )
    {
        if( choise == NULL )
        {
            return m_is->get();
        }
        
        string chs = choise;

        char chr = m_is->get();
        while( chs.find(chr) == string::npos )
        {
            print( "Please choose from " + chs );
            chr = m_is->get();
        }

        return chr;
    }

    string std_console::getstr( const char* choise )
    {
        if( choise == NULL )
        {
            string str;
   
            std::getline( *m_is, str );

            return str;
        }
        
        vector<string> chs;
        split( chs, choise, is_any_of( " /" ), token_compress_on );
        
        string str;
        std::getline( *m_is, str );
        
        while( std::find( chs.begin(), chs.end(), str ) == chs.end() )
        {
            print( "Please choose from " + string(choise) );    
            std::getline( *m_is, str );
        }

        return str;
    }
 */

} // namespace mort

