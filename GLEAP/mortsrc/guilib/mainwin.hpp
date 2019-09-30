#ifndef GTKLEAP_MAINWIN_HPP
#define GTKLEAP_MAINWIN_HPP

#include <common.hpp>

namespace mort
{
    enum framework_e { GTK, TER, QT, MOTIF, MFC };

    class root_t;

    class drawing_i;
    
    class console_t;

    class database_t;
 
    shared_ptr< console_t >& console();

    database_t& content();

    string find_file( const string& name );

} // namespace mort
        
#endif

