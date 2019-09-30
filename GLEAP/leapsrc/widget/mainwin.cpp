#include <stdexcept>
#include <object.hpp>
#include <format.hpp>
#include "mainwin.hpp"
#include "mainwui.hpp"
#include "fixbond.hpp"
#include "addions.hpp"
#include "solvate.hpp"
#include "addhydr.hpp"
#include "setpchg.hpp"
#include "ffparam.hpp"
#include "shidedlg.hpp"
#include "labeldlg.hpp"



extern "C"
{
    
    void on_mainwin_delete_event(GtkWidget* widget, gpointer user_data)
    {
        gtk_main_quit();
    }
    

    void on_mainwin_destroy_signal(GtkWidget* widget, gpointer user_data)
    {
    }
}

int popup_message( const string& msg )
{
    GtkWidget* dlg = gtk_message_dialog_new( NULL, GTK_DIALOG_DESTROY_WITH_PARENT,
	                                         GTK_MESSAGE_ERROR,
			                         GTK_BUTTONS_CLOSE,
                                             msg.c_str() );
    int res = gtk_dialog_run( GTK_DIALOG(dlg) );
    gtk_widget_destroy(dlg);
    return res;
}



using namespace mort;

string remove_dir( const string& full )
{
#ifdef WIN32
    int pos = full.rfind( '\\' );
#else
    int pos = full.rfind( '/'  );
#endif

    if( pos == string::npos )
        return full;
    
    return full.substr( pos+1, full.size() );
}

string remove_ext( const string& full )
{
    int pos = full.rfind( '.' );
    
    if( pos == string::npos )
        return full;
    
    return full.substr( 0, pos );
}


string uniq_name( const string& path )
{
    string file = remove_dir( path );
    return remove_ext(file);
}

main_window::main_window()
{
    m_glarea = shared_ptr<gtk_glarea>( new gtk_glarea() );
    drawing() = m_glarea;

    m_pmainwin = NULL;
    m_pmanager = NULL;
    m_pmenubar = NULL;
    m_ptoolbar = NULL;
    m_pviewpanel = NULL;
    m_pstatusbar = NULL;
}

main_window::~main_window()
{
}

void main_window::init()
{
    assert( m_pmainwin == NULL );

    m_pmainwin = gtk_window_new(GTK_WINDOW_TOPLEVEL);

    // connect delete_event and destroy signal
    g_signal_connect(G_OBJECT(m_pmainwin), "delete_event", G_CALLBACK(on_mainwin_delete_event),   this);
    g_signal_connect(G_OBJECT(m_pmainwin), "destroy",      G_CALLBACK(on_mainwin_destroy_signal), this);
    
    // set property
    gtk_window_set_title(GTK_WINDOW(m_pmainwin), "gLEaP");
    gtk_window_set_position(GTK_WINDOW(m_pmainwin), GTK_WIN_POS_CENTER_ALWAYS);
    gtk_window_set_default_size(GTK_WINDOW(m_pmainwin), 800, 600);

    GtkWidget* pbox = gtk_vbox_new(true, DEFAULT_BOX_PADDING_SIZE);
    gtk_container_add( GTK_CONTAINER(m_pmainwin), pbox );
    gtk_box_set_homogeneous(GTK_BOX(pbox), false);
    gtk_widget_show( pbox );

    init_manager();
    init_menubar( GTK_BOX(pbox) );
    init_toolbar( GTK_BOX(pbox) );
    init_viewpanel( GTK_BOX(pbox) );
    init_statusbar( GTK_BOX(pbox) );
}


void main_window::init_manager( )
{
    GtkActionGroup* pactgrp = gtk_action_group_new("MainwinActionGroup");

    // notify gettext to translate menu label and tooltip
    // gtk_action_group_set_translation_domain(_M_action_group_ptr, PACKAGE);

    gtk_action_group_add_actions( pactgrp, MAINWIN_ACTION_ENTRYS, MAINWIN_ACTION_NENTRY, this);
    m_pmanager = gtk_ui_manager_new();
    gtk_ui_manager_insert_action_group(m_pmanager, pactgrp, 0);
	      
    g_assert(m_pmainwin);

    gtk_window_add_accel_group(GTK_WINDOW(m_pmainwin), gtk_ui_manager_get_accel_group(m_pmanager));
	      
    if(!gtk_ui_manager_add_ui_from_string(m_pmanager, MAINWIN_UIDESIGN, -1, NULL))
    {
        throw std::runtime_error("Error: can't load ui from string");
    }

}



void main_window::init_menubar(GtkBox* pbox)
{
    assert( m_pmanager );
    m_pmenubar = gtk_ui_manager_get_widget(m_pmanager, "/Menubar");

    assert( m_pmenubar );
    gtk_widget_show( m_pmenubar);

    gtk_box_pack_start( pbox, GTK_WIDGET(m_pmenubar), false, false, DEFAULT_BOX_PADDING_SIZE);
}

void main_window::init_toolbar(GtkBox* pbox)
{
    g_assert( m_pmanager );
    m_ptoolbar = gtk_ui_manager_get_widget(m_pmanager, "/Toolbar");

    g_assert( m_ptoolbar );
    gtk_widget_show( m_ptoolbar);

    gtk_box_pack_start( pbox, GTK_WIDGET(m_ptoolbar), false, false, DEFAULT_BOX_PADDING_SIZE);
}

void main_window::init_viewpanel( GtkBox* pbox )
{
    m_pviewpanel = gtk_frame_new(NULL);
    gtk_box_pack_start( pbox, m_pviewpanel, false, false, DEFAULT_BOX_PADDING_SIZE);
    gtk_widget_show( m_pviewpanel );

    GtkWidget* hbox = gtk_hbox_new(false, 0);
    gtk_widget_show( hbox );
    
    gtk_box_pack_start( GTK_BOX(hbox), (GtkWidget*)(m_editbar), false, false, DEFAULT_BOX_PADDING_SIZE );
    gtk_box_pack_start( GTK_BOX(hbox), (GtkWidget*)(*m_glarea), false, false, DEFAULT_BOX_PADDING_SIZE );
    gtk_container_add( GTK_CONTAINER(m_pviewpanel), hbox );

    gtk_widget_grab_focus( (GtkWidget*)(*m_glarea) );
}

void main_window::init_statusbar( GtkBox* pbox )
{
    g_assert( pbox );

    m_pstatusbar = gtk_statusbar_new();
    gtk_widget_show( m_pstatusbar );

    gtk_box_pack_start( pbox, m_pstatusbar, false, false, DEFAULT_BOX_PADDING_SIZE);
    m_status.init(GTK_STATUSBAR(m_pstatusbar), "Global Context");
}

void main_window::show()
{
    init();

    assert( m_pmainwin != NULL );

    gtk_widget_show( m_pmainwin );

    m_status.push("Ready");
}
    
void moveto_center(molecule_t& m)
{
    int natom = m.natom();
    double* crd = get_vptr( m, ATOM, POSITION );

    double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;
    for( int i=0; i < natom; ++i )
    {
        cx += crd[3*i  ];
	cy += crd[3*i+1];
	cz += crd[3*i+2];
    }

    cx /= natom;
    cy /= natom;
    cz /= natom;

    for( int i=0; i < natom; ++i )
    {
        crd[3*i  ] -= cx;
	crd[3*i+1] -= cy;
	crd[3*i+2] -= cz;
    }

}


void main_window::on_file_open_action_listener()
{
    GtkWidget* pdlg = gtk_file_chooser_dialog_new( "Open ...", GTK_WINDOW(m_pmainwin), 
                                                   GTK_FILE_CHOOSER_ACTION_OPEN,
                                                   GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                                   GTK_STOCK_OPEN,   GTK_RESPONSE_ACCEPT,
                                                   NULL);

    GtkFileFilter* filter = gtk_file_filter_new();
    gtk_file_filter_set_name( filter, "Molecular files" );
    gtk_file_filter_add_pattern( filter, "*.sdf" );
    gtk_file_filter_add_pattern( filter, "*.pdb" );
    gtk_file_filter_add_pattern( filter, "*.mol2" );
    gtk_file_chooser_add_filter( GTK_FILE_CHOOSER(pdlg), filter );
   
    if( gtk_dialog_run(GTK_DIALOG(pdlg))==GTK_RESPONSE_ACCEPT )
    {
        string fn = gtk_file_chooser_get_filename( GTK_FILE_CHOOSER(pdlg) );
        molecule_ptr pmol( new molecule_t() );

        try
        {
            load_mol( fn, *pmol );
        }
        catch( std::exception& e )
        {
            popup_message( e.what() );
        }

        moveto_center( *pmol );
        string name = uniq_name(fn);
        pmol->set_s(NAME, name );
        m_glarea->add( "molecule", *pmol, "line" );
        content().set( name, pmol );
        m_curt = pmol;
    }
    
    gtk_widget_destroy( pdlg );    
    
}

void main_window::on_file_close_action_listener()
{
    dlgfactory_t::run( "clsdata" );
}

void main_window::on_file_save_action_listener()
{
    GtkWidget* pdlg = gtk_file_chooser_dialog_new( "Open ...", GTK_WINDOW(m_pmainwin), 
                                                   GTK_FILE_CHOOSER_ACTION_SAVE,
                                                   GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                                                   GTK_STOCK_SAVE,   GTK_RESPONSE_ACCEPT,
                                                   NULL);

    strcombo_t molname( "Select molecule" );
    molname.init( content() );
    gtk_box_pack_start( GTK_BOX(GTK_DIALOG(pdlg)->vbox), (GtkWidget*)molname, true, true, 0 );
    

    if( gtk_dialog_run(GTK_DIALOG(pdlg))==GTK_RESPONSE_ACCEPT )
    {
        string fn = gtk_file_chooser_get_filename( GTK_FILE_CHOOSER(pdlg) );

        molecule_ptr pmol = content().get_mol( molname.getstr() );
        
        assert( pmol !=NULL );
        
        save_mol( fn, *pmol );
    }
    
    gtk_widget_destroy( pdlg );    
    
}

void main_window::on_file_quit_action_listener()
{
    gtk_main_quit();
}


void main_window::on_edit_fixbond_action_listener()
{
    dlgfactory_t::run( "fixbond" );
}

void main_window::on_edit_addhydr_action_listener()
{
    dlgfactory_t::run( "addhydr" );
}

void main_window::on_edit_setpchg_action_listener()
{
    dlgfactory_t::run( "setpchg" );
}

void main_window::on_edit_parmchk_action_listener()
{
    dlgfactory_t::run( "parmchk" );
}

void main_window::on_edit_fastbld_action_listener()
{
    dlgfactory_t::run( "fastbld" );
}

void main_window::on_edit_addions_action_listener()
{
    dlgfactory_t::run( "addions" );
}

void main_window::on_edit_solvate_action_listener()
{
    dlgfactory_t::run( "solvate" );
}

void main_window::on_edit_ffparam_action_listener()
{
    dlgfactory_t::run( "ffparam" );
}

void main_window::on_view_shide_action_listener()
{
    dlgfactory_t::run( "shide" );
}

void main_window::on_view_label_action_listener()
{
    dlgfactory_t::run( "label" );
}

void main_window::on_view_color_action_listener()
{
}

void main_window::on_view_ribbon_action_listener()
{
    dlgfactory_t::run( "ribbdlg" );
}

void main_window::on_view_surface_action_listener()
{
}
