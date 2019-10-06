#include <algorithm>
#include <cstring>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <numeric>
#include <boost/algorithm/string.hpp>
#include "constant.hpp"

namespace mort
{
    using namespace std;

    vector<string> split( const string& line, const string& dels )
    {
        vector<string> r;
        boost::split( r, line, boost::is_any_of(dels), boost::token_compress_on );

        vector<string>::iterator i = std::find( r.begin(), r.end(), "" );
        while( i != r.end() )
        {
            r.erase( i );
            i = std::find( r.begin(), r.end(), "" );
        }

        return r;
    }


    string next_word( istream& stream, char deli )
    {
        static const string empty = " \t\r\n";
        
        string word;
        
        char c = stream.peek();
        
        while( stream && c != deli && empty.find( c ) != string::npos )
        {
            stream.ignore();
            c = stream.peek();
        }
        
        while( stream && c != deli && empty.find( c ) == string::npos )
        {
            word.append( 1, c );
            stream.ignore();
	    while (c == -1 && stream) {
	    	stream.get();
	    }
            c = stream.peek();
        }

        if( !stream && word.length() >= 1 )
        {
            word.erase( word.length() - 1 );
        }

        return word;        
    }

    const char* skip_digit( const char* ptr )
    {
        while( isdigit( *ptr ) )
        {
            ptr++;
        }
        
        return ptr;
    }

    const char* skip_alnum( const char* ptr )
    {
        while( isalnum( *ptr ) )
        {
            ptr++;
        }
        
        return ptr;
    }
    
    const char* skip_float( const char* ptr )
    {
        ptr = skip_digit( ptr );
        
        if( *ptr != '.' )
        {
            return ptr;
        }

        ptr++;
        return skip_digit( ptr );
    }

    const char* skip_alpha( const char* ptr )
    {
        while( isalpha( *ptr ) )
            ++ptr;

        return ptr;
    }

    std::string next_digit( const char* ptr )
    {
        const char* end = skip_digit( ptr );
        
        return std::string( ptr, end );
    }
    
    std::string next_alnum( const char* ptr )
    {
        const char* end = skip_alnum( ptr );
        
        return std::string( ptr, end );
    }

    string strip_quota(const string& name )
    {
        string tmp(name);

        string::iterator quota_begin = std::remove( tmp.begin(), tmp.end(), '\"' );
            
        tmp.erase( quota_begin, tmp.end() );
            
        return tmp;
    }

    void replace( string& str, char a, char b )
    {
        for( int i=0; i < (int)str.length(); ++i )
        {
            if( str[i] == a )
            {
                str[i] = b;
            }
        }
    }
        
    string replace_copy( const string& str, char a, char b )
    {
        string copy( str );
            
        for( int i=0; i < (int)copy.length(); ++i )
        {
            if( copy[i] == a )
            {
                copy[i] = b;
            }
        }
            
        return copy;
    }

    bool empty( const string& str )
    {
        for( int i=0; i < (int)str.length(); ++i )
        {
            if( str[i] != ' ' && str[i] != '\n' && str[i] != '\t' && str[i] != '\r' )
            {
                return false;
            }
        }
        
        return true;
    }    

    int count_item( const string& line )
    {
        int i=0;
        string item;       
        istringstream ls( line );
     
        while( ls >> item )
        {
            i++;
        }
        
        return i;
    }


    string peek_keyw( istream& stream )
    {
        if( !stream )
        {
            return "    ";
        }

        char line[MAX_LINE_WIDTH];

        std::ios::pos_type pos = stream.tellg();
            
        stream.getline(line, MAX_LINE_WIDTH);
            
        stream.seekg(pos);

        if( strlen(line) > 4 )
        {
            return string(line, line+4);
        }
        else if( strlen(line) > 0 )
        {
            return line;
        }
  
        return "";
    } 



} // namespace mort
