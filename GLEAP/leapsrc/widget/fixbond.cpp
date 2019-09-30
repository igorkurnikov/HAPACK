#include <mortgl.hpp>
#include <moloper.hpp>
#include "fixbond.hpp"

namespace mortgtk
{
    fixbond_t::fixbond_t()
        :basicdlg_t( "Fix bond order automatically" )
    {
        dlgfactory_t::add( "fixbond", this );
    }
    

    fixbond_t::~fixbond_t()
    {
    }

    command_ptr fixbond_t::get_command()
    {
        return command_ptr( new amber::moloper_command( "fixbond", molname() ) );
    }
    
    basicdlg_ptr fixbond_t::clone() const
    {
        return basicdlg_ptr( new fixbond_t() );
    }

} // namespace mortgtk

    
