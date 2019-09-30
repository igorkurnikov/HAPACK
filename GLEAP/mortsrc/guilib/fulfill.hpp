#ifndef GTKLEAP_FULFILL_HPP
#define GTKLEAP_FULFILL_HPP

#include <common.hpp>

namespace mort
{
    string fulfill( const string& cmd );
    
    string fulfill_cmd( const string& part );

    string fulfill_file( const string& part );
    
} // namespace mort

#endif

