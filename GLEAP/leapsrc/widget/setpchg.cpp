#include <stdexcept>

#include <fstream>
#include <mortgl.hpp>
#include <format.hpp>
#include <moloper.hpp>
#include "setpchg.hpp"

namespace mortgtk
{
    using namespace mort;

    setpchg_t::setpchg_t()
        :basicdlg_t( "Set partial charge using antechamber" )
    {
        dlgfactory_t::add( "setpchg", this );
    }
    

    setpchg_t::~setpchg_t()
    {
    }
    

    command_ptr setpchg_t::get_command()
    {
        return command_ptr( new amber::moloper_command( "setpchg", molname()) );
    }
    
} // namespace mortgtk

    
