#include <plugins/solvate.hpp>
#include "solvate.hpp"

namespace mortgtk
{
    
    solvate_t::solvate_t()
        : basicdlg_t( "Add solvation" ),
          m_solvent( "Select solvent:" )
    {
        database_ptr pmdb = content().get_mdb( "_solvents" );
        assert( pmdb != NULL );
        m_solvent.init( *pmdb );
        m_soltype.init( );

        gtk_box_pack_start( GTK_BOX(vbox()), (GtkWidget*)m_solvent, true, true, 0 );
        gtk_box_pack_start( GTK_BOX(vbox()), (GtkWidget*)m_soltype, true, true, 0 );

        GtkWidget* sizelabel = gtk_label_new( "Select size:" );
        m_solsize = gtk_entry_new();

        gtk_widget_show( sizelabel );
        gtk_widget_show( m_solsize );

        gtk_box_pack_start( GTK_BOX(vbox()), sizelabel, true, true, 0 );
        gtk_box_pack_start( GTK_BOX(vbox()), m_solsize, true, true, 0 );
        
        dlgfactory_t::add( "solvate", this );
    }

    solvate_t::~solvate_t()
    {
    }

    command_ptr solvate_t::get_command( )
    {
        string soltype = m_soltype.getsol();
        string solvent = m_solvent.getstr();
        string solsize = gtk_entry_get_text( GTK_ENTRY(m_solsize) );
        string solcent;
        

        string prefix = soltype.substr(0, 10);
        if( prefix=="solvatecap" )
        {
            solcent = soltype.substr(11, soltype.length()-11 );
            soltype = "solvatecap";
        }
        
        return command_ptr( new amber::solvate_command(soltype, molname(), solvent, solcent, solsize) );
    }
    
} // namespace mortgtk
