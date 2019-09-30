#include <object.hpp>
#include <mortgl.hpp>
#include <moloper.hpp>
#include "addhydr.hpp"

namespace mortgtk
{
    using namespace mort;
    
    addhydr_t::addhydr_t()
        :basicdlg_t( "Add hydrogens" )
    {
        dlgfactory_t::add( "addhydr", this );
    }    

    addhydr_t::~addhydr_t()
    {
    }
    

    command_ptr addhydr_t::get_command()
    {
        return command_ptr( new amber::moloper_command( "addhydr", molname() ) );
    }
    
} // namespace mortgtk

    
