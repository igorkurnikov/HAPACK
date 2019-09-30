#include <sstream>
#include <boost/algorithm/string.hpp>
#include <common.hpp>
#include "grammar.hpp"

namespace mort
{
    using namespace boost;
 
    using std::runtime_error;

    using std::istringstream;

    void interpret( const string& command, vector<string>& args )
    {
        istringstream is( command );

        bool assign = false;
        while( is )
        {
            string arg =  parse( is );

            if( arg=="=" )
            {
                if( assign ) throw runtime_error( "Error: '=' appear more than one time" );

                if( args.size() != 1 ) throw runtime_error( "Error: no assignee or more than on assignee" );

                assign = true;
            }
            else if( !arg.empty() )
            {
                args.push_back(arg);
            }
        }

        if( assign )
        {
            if( args.size() < 2 ) 
            {
                throw runtime_error( "Error: no arguments for assignment" );
            }

            std::swap( args[0], args[1] );
        }
    }

    bool isblank_tab_newline_carriagereturn( char n )
    {
        return n==' ' || n=='\t' || n=='\n' || n=='\r';
    }

    string parse( istream& is )
    {
        char n = is.peek();
        while( is && isblank_tab_newline_carriagereturn(n) )
        {
            is.ignore();
            n = is.peek();
        }

        if( !is ) return "";

        if( n=='=' )
        {
            is.ignore();
            return "=";
        }

        if( n=='{' )
        {
            return parse_curly( is );
        }

        if( n=='"' )
        {
            return parse_quota( is );
        }

        if( n=='}' )
        {
            throw runtime_error( "Error: unpaired }" );
        }

        string reserve = "={}\" \n\t\r";
        string word;
  
        while( reserve.find(n)==string::npos )
        {
            if(is >> n)
            {
                word.append(1, n);
            }
            else
            {
                break;
            }

            n = is.peek();
        }
 
        return word;
    }

    string parse_curly( istream& is )
    {
        char c = is.get();
        assert( c == '{' );

        int ncurly = 1;
        string str = "{";

        while( ncurly > 0 )
        {
            c = is.get();

            if( !is ) throw std::runtime_error( "Error: unpaired {" );

            if( c=='{' ) ++ncurly;
            if( c=='}' ) --ncurly;
            str.append( 1, c );
        }

        return str;
    }

    string parse_quota( istream& is )
    {
        char c;
        is >> c;
        assert( c == '"' );

        string str;
        while( is >> std::noskipws >> c && c !='"' )
        {
            if(c !='"' ) str.append( 1, c );
        }

        if( c != '"' )
        {
            throw runtime_error( "Error: unpaired quota" );
        }

        return str;
    }

    int ndim( const string& list )
    {
        int ndimen = 0 ;
        int ncurly = 0;
        
        for( int i=0; i < (int)list.length(); i++ )
        {
            if( list[i] == '}' )
            {
                ncurly--;
            }

            if( ncurly < 0 )
            {
                throw std::runtime_error( "Error: bad list " + list );
            }

            if( list[i] == '{' )
            {
                ncurly++;
            }

            if( ndimen < ncurly )
            {
                ndimen = ncurly;
            }
        }
           
        if( ncurly != 0 )
        {
            throw runtime_error( "Error: bad list " + list );
        }

        return ndimen;
    }

    vector<string> interpret1d( const string& list )
    {
        assert( ndim(list)==1 && list[0]=='{' );

        vector<string> args;

        split( args, list, is_any_of(" {}"), token_compress_on );

        return vector<string>( args.begin() + 1, args.end() - 1 );
    }

    vector< vector<string> > interpret2d( const string& list )
    {
        assert( ndim(list)==2 && list[0]=='{' );

        vector< vector<string> > result;

        int start = list.find( '{' );   

        start = list.find( '{', start + 1);

        while( start != (int)string::npos )
        {
            int end = list.find( '}', start ) + 1;

            result.push_back( interpret1d( list.substr( start, end - start ) ) );

            start = list.find( '{', end + 1 ); 
        }

        return result;
    }

    bool closed( const string& list, char bgn, char end )
    {
        int ncurly=0;
        
        for( int i=0; i < (int)list.length(); ++i )
        {
            if( list[i] == bgn )
            {
                ncurly++;
            }
            else if( list[i] == end )
            {
                ncurly--;
            }

            if( ncurly < 0 )
            {
                throw runtime_error( "Error: bad list " + list );
            }
        }
        
        return ( ncurly == 0 );
    }

    bool curly_closed( const string& list )
    {
        return closed( list, '{', '}' );
    }

    bool quota_closed( const string& list )
    {
        return closed( list, '"', '"' );
    }







    
} // namespace mort

