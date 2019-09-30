#include <string>
#include "editbar.hpp"
#include "editbui.hpp"


editbar_t::editbar_t( )
    : m_buttons(NUM_TOOL)
{
    m_self = gtk_vbox_new(false, 0);
    gtk_widget_show( m_self );

    for(int i=0; i < NUM_TOOL; ++i )
    {
        char* amberhome = getenv( "AMBERHOME" );
        assert( amberhome != NULL );

        std::string ficon(amberhome);
        ficon += "/dat/leap/pixmaps/";
        ficon += edit_tool_entries[i].icon_file;

        GtkWidget* picon = gtk_image_new_from_file(ficon.c_str());
        									    
 	if(i == 0)
 	    m_buttons[i] = gtk_radio_button_new(NULL);
        else
            m_buttons[i] = gtk_radio_button_new_from_widget(GTK_RADIO_BUTTON(m_buttons[0]));

        gtk_toggle_button_set_mode(GTK_TOGGLE_BUTTON(m_buttons[i]), FALSE);

        gtk_widget_set_size_request(m_buttons[i], 34, 34);

        gtk_widget_show( m_buttons[i] );

        gtk_widget_show( picon );
	
	gtk_container_add(GTK_CONTAINER(m_buttons[i]), picon);

        g_signal_connect(G_OBJECT(m_buttons[i]), 
	                 "clicked",
                         G_CALLBACK(edit_tool_entries[i].callback),
                         this);

        gtk_box_pack_start( GTK_BOX(m_self), m_buttons[i], false, false, 0 );
    }

}

editbar_t::~editbar_t()
{
}


void editbar_t::show()
{
}


