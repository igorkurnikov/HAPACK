#ifndef MOLVIEW_SYSLOG_H
#define MOLVIEW_SYSLOG_H

#include <common.hpp>

namespace mort
{
    /// \brief log system for leap
    /// \ingroup leaplib
    class leaplog_t
    {
    private:

        leaplog_t();
    
        virtual ~leaplog_t();
    
    public:

        /// add line to log
        static void putline( const string& line );
    
        /// get a line from log
        /// 
        /// \return  true  if log got successfully
        ///          false if no log in system
        static bool getline( string& line );

        static vector< string >::iterator cur();
        
        static vector< string >::iterator end();

        /// move current cursor 
        ///
        /// \param the distance for move, could be negative
        static void move_cursor( int size );
        
    private:

        static vector< string > m_logs;

        static int m_cursor;
    
    };
    
} // namespace mort

#endif

        
