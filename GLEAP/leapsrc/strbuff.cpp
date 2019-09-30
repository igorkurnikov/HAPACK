#include <cstdio>
//#include <boost/algorithm/string.hpp>
#include <readline/readline.h>
#include <readline/history.h>

#include <common.hpp>
#include <object.hpp>
#include "strbuff.hpp"

namespace amber
{
    using namespace std;
    using namespace mort;

    std_console::std_console( istream& is, ostream& os )
    {
        m_is = &is;
        m_os = &os;

	string batch;
	m_batch = mortenv().get_s("batch", batch) && (batch=="on");
    }
    
    std_console::~std_console( )
    {
    }
    
    void std_console::print( const string& str )
    {
        *m_os << str << std::flush;
    }

    bool std_console::getline( string& line, const char* prompt )
    {
        if( m_batch )
	{
	    if(prompt)
	        *m_os << prompt << std::flush;

	    return std::getline( *m_is, line );
	}
   
        assert( m_is == &std::cin );
        char* lptr = readline( prompt );
	if( lptr == NULL )
	    return false;
    
        if( !empty(lptr) )
            add_history( lptr );

        line = string(lptr);
	return true;
    }

    char std_console::getchar()
    {
        char emptys[] = "\n\t ";
        
        char o = m_is->get();
        while( std::find( emptys, emptys+3, o ) != emptys+3 )
        {
            o = m_is->get();
        }

        return o;
    }


} // namespace amber

