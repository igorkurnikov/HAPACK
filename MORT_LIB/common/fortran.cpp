#include <cassert>
#include <iostream>
#include <stdexcept>
#include <cstdio>
#include "fortran.hpp"
#include "stralgo.hpp"
#include "constant.hpp"

namespace mort
{
    using namespace std;
    using boost::any_cast;

    void write_ivalue( ostream& os, const string& cfmt, const any& value )
    {
        char result[100];

        try
        {
            int v = any_cast<int>( value );
            sprintf( result, cfmt.c_str(), v );
            os << result;
        }
        catch( std::bad_cast& e )
        {
            throw runtime_error( "Error: during print format: " + cfmt + " : " + e.what() );
        }

    }

    void write_dvalue( ostream& os, const string& cfmt, const any& value )
    {
        char result[MAX_LINE_WIDTH];

        try
        {
            double v = any_cast<double>( value );
            sprintf( result, cfmt.c_str(), v );
            os << result;
        }
        catch( std::bad_cast& e )
        {
            throw runtime_error( "Error: during print format: " + cfmt + " : " + e.what() );
        }

    }

    void write_svalue( ostream& os, const string& cfmt, const any& value )
    {
        char result[100];
        try
        {
            string v = any_cast<string>( value );
            sprintf( result, cfmt.c_str(), v.c_str() );
            os << result;
        }
        catch( std::bad_cast& e )
        {
            throw runtime_error( "Error: during print format: " + cfmt + " : " + e.what() );
        }
    }

    format::format( const string& strfmt )
        : m_strfmt( strfmt )
    {
    }

    format::format( const format& rhs )
        : m_strfmt( rhs.m_strfmt ),
          m_values( rhs.m_values )
    {
    }

    format::~format()
    {
    }

    format format::operator%( int i )
    {
        format f2( *this );
        f2.m_values.push_back( i );
        return f2;
    }

    format format::operator%( size_t i )
    {
        format f2( *this );
        f2.m_values.push_back( int(i) );
        return f2;
    }


    format format::operator%( double d )
    {
        format f2( *this );
        f2.m_values.push_back( d );
        return f2;
    }

    format format::operator%( const string& s )
    {
        format f2( *this );
        f2.m_values.push_back( s );
        return f2;
    }

    void format::write( ostream& os, const string& cfmt, const any& v ) const
    {
        char type = cfmt[ cfmt.length()-1 ];

        if( type=='d' )
        {
            write_ivalue( os, cfmt, v );
        }
        else if( type=='f' || type=='e' )
        {
            write_dvalue( os, cfmt, v );
        }
        else if( type=='s' )
        {
            write_svalue( os, cfmt, v );
        }
        else
        {
            throw runtime_error( "Error: unknown format string: " + cfmt );
        }
    }
               

    void format::write( ostream& os ) const
    {
        size_t bgn(0), vid(0);
        size_t end = m_strfmt.find( '%' );

        while( end != string::npos )
        {
            // print content between format strings
            //
      
            os << m_strfmt.substr( bgn, end - bgn );
           
            bgn = end;
            while( end < m_strfmt.length() && 
                   m_strfmt[end] != 'f' && m_strfmt[end] != 'e' &&
                   m_strfmt[end] != 'd' && m_strfmt[end] != 's' )
            {
                end++;
            }

            if( end==m_strfmt.length() )
            {
                throw runtime_error( "Error: cannot understand format " + m_strfmt.substr(bgn, end-bgn) );
            }

            end++;
            string cfmt = m_strfmt.substr( bgn, end-bgn );
 
            if( vid==m_values.size() )
            {
                throw runtime_error( "Error: no argument was fed to format string: " + cfmt );
            }

            write( os, cfmt, m_values[vid] );
            vid++;

            bgn = end;
            end = m_strfmt.find( '%', bgn );
        }         
              
        os << m_strfmt.substr( bgn, m_strfmt.length()-bgn );

    } 

    ostream& operator<<( ostream& os, const format& f )
    {
        f.write( os );
        return os;
    }

    fortran_t::fortran_t( ostream& os, const string& format )
        :m_os( os ), m_format( format )
    {
        m_size = atoi( format.c_str() );
        m_form = skip_digit( format.c_str() );
        m_index = 0;
    }

    void fortran_t::begin()
    {
        m_index = 0;
    }

    void fortran_t::step()
    {
        m_index++;

        if( m_index % m_size == 0 )
        {
            m_os << endl;
        }
    }
        
    void fortran_t::end()
    {
        if( m_index==0 )
        {
            m_os << endl;
        }

        if( m_index % m_size != 0 )
        {
            m_os << endl;
        }
    }

    void fortran_t::operator()( const string& value )
    {
        assert( m_form[0] == 'a' || m_form[0] == 'I' );
        
        string prefix = ( m_form[0]=='a' ? "%-" : "%" );
  
        string cform = prefix + m_form.substr( 1, m_form.length() ) + "s";
        
        m_os << format( cform ) % value;

        step();
    }

    void fortran_t::operator()( double value )
    {
        assert( m_form[0] == 'E' || m_form[0]=='F');

        char f = m_form[0] - 'A' + 'a';
        
        string cform = "%" + m_form.substr( 1, m_form.length() );

        cform.append( 1, f );
        
        m_os << format( cform ) % value;
            
        step();
    }
    
    void fortran_t::operator()( int value )
    {
        assert( m_form[0] == 'I' || m_form[0] == 'E' || m_form[0] == 'F' );

        if( m_form[0] == 'I' )
        {
            string cform = "%" + m_form.substr( 1, m_form.length() ) + "d";
        
            m_os << format( cform ) % value;
                
            step();
        }
        else
        {
            double dval = value;
            operator()( dval );
        }
    }    
}

