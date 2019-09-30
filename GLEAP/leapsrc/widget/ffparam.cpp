#include <plugins/saveprm.hpp>
#include "ffparam.hpp"

namespace mortgtk
{
    
    ffparam_t::ffparam_t()
        : basicdlg_t( "Save force field parameter..." ),
          m_fftype( "Select force field:" )
    {
        m_fftype.init( "amber", "amoeba" );
        
        GtkWidget* toplabel = gtk_label_new( "Topology file:" );
        GtkWidget* xyzlabel = gtk_label_new( "Position file:" );
        m_topfil = gtk_entry_new( );
        m_xyzfil = gtk_entry_new( );
        gtk_widget_show( toplabel );
        gtk_widget_show( xyzlabel );
        gtk_widget_show( m_topfil );
        gtk_widget_show( m_xyzfil );


        gtk_box_pack_start( GTK_BOX(vbox()), (GtkWidget*)m_fftype, true, true, 0 );
        gtk_box_pack_start( GTK_BOX(vbox()), toplabel, true, true, 0 );
        gtk_box_pack_start( GTK_BOX(vbox()), m_topfil, true, true, 0 );
        gtk_box_pack_start( GTK_BOX(vbox()), xyzlabel, true, true, 0 );
        gtk_box_pack_start( GTK_BOX(vbox()), m_xyzfil, true, true, 0 );        

        onchange(NULL);

        dlgfactory_t::add( "ffparam", this );
    }
    
    ffparam_t::~ffparam_t()
    {
    }
    

    void ffparam_t::onchange( GtkWidget* )
    {
        string mname = molname();
        string top = mname + ".top";
        string xyz = mname + ".xyz";
        gtk_entry_set_text( GTK_ENTRY(m_topfil), top.c_str() );
        gtk_entry_set_text( GTK_ENTRY(m_xyzfil), xyz.c_str() );
    }
    
    command_ptr ffparam_t::get_command()
    {
        string fftype = m_fftype.getstr();
        string topfil = gtk_entry_get_text( GTK_ENTRY(m_topfil) );
        string xyzfil = gtk_entry_get_text( GTK_ENTRY(m_xyzfil) );
        string action = "save" + fftype + "parm";
        
        return command_ptr( new amber::saveprm_command(action, molname(), topfil, xyzfil) );    
    }

} // namespace mortgtk


 
