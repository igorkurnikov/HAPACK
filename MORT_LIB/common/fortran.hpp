#ifndef MORT_SANDER_FORMAT_HPP
#define MORT_SANDER_FORMAT_HPP

#include <iosfwd>
#include <string>
#include <vector>
#include <boost/any.hpp>

namespace mort
{
    using std::string;
    using std::vector;
    using std::ostream;
    using boost::any;

    class format
    {
    public:

        format( const string& strfmt );

        format( const format& rhs );

        virtual ~format();

        format operator%( int i);

        format operator%( size_t i );

        format operator%( float f );

        format operator%( double d );

        format operator%( const string& s );

        void write( ostream& os ) const;

    private:

        void write( ostream& os, const string& cfmt, const any& value ) const;

        string m_strfmt;

        vector<any> m_values;
    };
   
    ostream& operator<<( ostream& os, const format& f );
           
    class fortran_t
    {
    public:

        fortran_t( ostream& os, const string& format );

        void begin();
        
        void step();
        
        void end();
        
        void operator()( int value );
        
        void operator()( double value );

        void operator()( const string& value );
        
    private:

        ostream& m_os;

        int m_size;

        int m_index;

        string m_form;

        string m_format;
        
    };

    template< typename T >
    void write_size( fortran_t* fortran, const vector< T >& vec )
    {
        (*fortran)( int(vec.size()) );
    }

    template< typename T >
    void write_list( fortran_t* fortran, const vector< T >& vec )
    {
        // for_each( vec.begin(), vec.end(), *fortran );
        for( int i=0; i < (int)vec.size(); ++i )
        {
            (*fortran)( vec[i] );
        }
    }
    
} // namespace mort

#endif

