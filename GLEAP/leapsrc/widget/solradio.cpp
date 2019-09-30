#include <stdexcept>
#include "solradio.hpp"

namespace mortgtk
{
    
    solradio_t::solradio_t( )
        : strradio_t( "Solvation type:" )
    {
        m_maskentry = gtk_entry_new( );
    }
    
    solradio_t::~solradio_t()
    {
    }
    

    void solradio_t::init( )
    {
        strradio_t::init( "box", "oct", "shell", "cap" );
        gtk_box_pack_start( GTK_BOX(vbox()), m_maskentry, true, true, 0 );
    }
    
    void solradio_t::onchange( GtkWidget* w )
    {
        strradio_t::onchange( w );
        
        string type = getstr();
        
        if( type=="cap" )
        {
            gtk_widget_show( m_maskentry );
        }
        else
        {
            gtk_widget_hide( m_maskentry );
        }
    }
    
    string solradio_t::getsol( )
    {
        string type = getstr();
        
        if( type=="cap" )
        {
            string mask = gtk_entry_get_text( GTK_ENTRY(m_maskentry) );
            
            if( mask.empty() )
            {
                throw std::runtime_error( "Error: must select a center for cap" );
            }
            
            return "solvatecap " + mask;
        }
        
        return "solvate" + type;
    }
    
} // namespace mortgtk
