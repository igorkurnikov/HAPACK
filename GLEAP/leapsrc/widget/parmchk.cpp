#include <stdexcept>

#include <fstream>
#include <mortgl.hpp>
#include <format.hpp>
#include <moloper.hpp>
#include "parmchk.hpp"

namespace mortgtk
{
    using namespace mort;

    parmchk_t::parmchk_t()
        :basicdlg_t( "Running parmchk to load extra force field parameter" )
    {
        dlgfactory_t::add( "parmchk", this );
    }
    

    parmchk_t::~parmchk_t()
    {
    }
    

    command_ptr parmchk_t::get_command()
    {
        return command_ptr( new amber::moloper_command( "parmchk", molname()) );
    }
    
} // namespace mortgtk

    
