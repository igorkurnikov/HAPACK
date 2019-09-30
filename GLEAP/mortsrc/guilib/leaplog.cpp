#include "leaplog.hpp"

namespace mort
{

    leaplog_t::leaplog_t()
    {
    }

    leaplog_t::~leaplog_t()
    {
    }

    void leaplog_t::putline( const string& line )
    {
    	std::cout << line << std::endl;
        //m_logs.push_back( line );
    }

    bool leaplog_t::getline( string& line )
    {
        if( m_cursor == (int)m_logs.size() )
        {
            return false;
        }
    
        line = m_logs[ m_cursor ];
    
        m_cursor++;
    
        return true;
    }

    vector< string >::iterator leaplog_t::cur()
    {
        return m_logs.begin() + m_cursor;
    }
    
    vector< string >::iterator leaplog_t::end()
    {
        return m_logs.end();
    }

    void leaplog_t::move_cursor( int size )
    {
        m_cursor += size;
    }
    

    vector< string > leaplog_t::m_logs;
    
    int leaplog_t::m_cursor;

} // namespace mort
