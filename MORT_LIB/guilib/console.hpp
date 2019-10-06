#ifndef MOLVIEW_CONSOLE_H
#define MOLVIEW_CONSOLE_H

#include <common.hpp>
#include "history.hpp"

namespace mort
{

    /// \brief processing keyboard event
    /// \ingroup leaplib
    class console_t
    {
    public:

        /// constructor
        console_t( );
    
        /// deconstructor
        virtual ~console_t();

        bool mainloop();

        bool process( const string& input );

        virtual void print( const string& str ) = 0;

#ifdef getchar
#  undef getchar
#endif
        virtual char getchar( ) = 0;
        
        virtual bool getline( string& line, const char* prompt ) = 0;

      private:

        string m_pending;
    };

    typedef shared_ptr< console_t > console_ptr;

} // namespace mort

#endif
