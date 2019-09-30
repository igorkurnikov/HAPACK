#include <stdexcept>

#include <GL/gl.h>
#include <GL/glu.h>

#include <gtk/gtk.h>
#include <gtk/gtkgl.h>
#include <gdk/gdkkeysyms.h>
#include <object.hpp>


#include "glarea.hpp"
#include "mainwin.hpp"

#ifdef WIN32
#  define EXPORT __declspec(dllexport)
#else
#  define EXPORT
#endif

using namespace mort;

using mortgtk::gtk_glarea;

extern "C"
{

EXPORT void on_drawing_realize(GtkWidget *widget, gpointer user_data)
{
    gtk_glarea* parea = (gtk_glarea*)user_data;
    parea->glbegin();
    parea->init();
    parea->glend();
}

EXPORT gboolean on_drawing_configure(GtkWidget *widget, GdkEventConfigure *event, gpointer user_data)
{
    gtk_glarea* parea = (gtk_glarea*)user_data;
    parea->glbegin();
    parea->resize( widget->allocation.width, widget->allocation.height );
    parea->glend();
    return FALSE;
}

EXPORT gboolean on_drawing_expose( GtkWidget *widget, GdkEventExpose  *event, gpointer user_data)
{
    gtk_glarea* parea = (gtk_glarea*)user_data;
    parea->glbegin();
    parea->repaint();
    parea->glend();
    return FALSE;
}

EXPORT gboolean on_drawing_button_press( GtkWidget* widget, GdkEventButton* event, gpointer user_data )
{
    try
    {
        gtk_glarea* parea = (gtk_glarea*)user_data;
        gtk_widget_grab_focus( (GtkWidget*)(*parea) );

        if( event->button == 3 )
        {
            parea->popup_menu(event);
        }
        else
        {
            parea->mouse_press( event->button, event->state, event->x, event->y );
        }
    }
    catch( std::exception& e )
    {
        popup_message( e.what() );
    }
    
    return FALSE;
}

EXPORT gboolean on_drawing_button_release( GtkWidget* widget, GdkEventButton* event, gpointer user_data )
{
    try
    {
        gtk_glarea* parea = (gtk_glarea*)user_data;
        parea->glbegin();
        parea->mouse_release( event->button, event->state, event->x, event->y );
        parea->repaint();
        parea->glend();
    }
    catch( std::exception& e )
    {
        popup_message( e.what() );
    }

    return FALSE;
}
    
EXPORT gboolean on_drawing_motion_notify( GtkWidget* widget, GdkEventMotion* event, gpointer user_data )
{
    gint x, y;
    gtk_glarea* parea = (gtk_glarea*)user_data;

    try
    {
        parea->glbegin();
        parea->mouse_move( 1, event->state, event->x, event->y );
        parea->repaint();
        parea->glend();
        gtk_widget_get_pointer( widget, &x, &y );
    }
    catch( std::exception& e )
    {
        popup_message( e.what() );
    }
    
    return FALSE;
}

EXPORT gboolean on_drawing_key_press( GtkWidget* widget, GdkEventKey* event, gpointer user_data )
{
    gtk_glarea* parea = (gtk_glarea*)user_data;
    parea->glbegin();

    if( event->keyval==GDK_Delete )
    {
        atomvec_ptr pselected = content().get_avec( "selected" );
        if( pselected )
        {
            set< molecule_t* > pmols;

            for( int i=0; i < pselected->size(); ++i )
            {
                moref_t a = pselected->at(i);
                a.getmol().remove_atom( a );
                pmols.insert( &a.getmol() );
            }

            std::for_each( pmols.begin(), pmols.end(), boost::mem_fn(&molecule_t::cleanup) );  

            pselected->clear();
        }

    }
  
    parea->update();
    parea->repaint();
    parea->glend();

    return FALSE;
}

} 

namespace mortgtk
{
    using namespace mort;

    void unselect( vector<moref_t>& mos )
    {
        for( int i=0; i < (int)mos.size(); ++i )
        {
            mos[i].set_i( SELECTED, 0 );
        }
    }
    



    const gchar* menu_info =
        "<ui>"
        "  <popup>"
        "    <menu action='ChangeElementAction'>"
        "      <menuitem action='ChangeCAction' />"
        "      <menuitem action='ChangeHAction' />"
        "      <menuitem action='ChangeOAction' />"
        "      <menuitem action='ChangeNAction' />"
        "      <menuitem action='ChangeSAction' />"
        "      <menuitem action='ChangePAction' />"
        "      <menuitem action='ChangeFAction' />"
        "      <menuitem action='ChangeClAction' />"
        "    </menu>"
        "    <menu action='ChangeOrderAction'>"
        "      <menuitem action='ChangeSingleAction'/>"
        "      <menuitem action='ChangeDoubleAction' />"
        "      <menuitem action='ChangeTripleAction' />"
        "    </menu>"
        "    <menuitem action='PropertyAction' />"
        "    <menuitem action='StyleAction' />"
        "  </popup>"
        "</ui>";
    
    void S_on_change_c_action(GtkToggleAction* paction, gpointer userdata );
    void S_on_change_h_action(GtkToggleAction* paction, gpointer userdata );
    void S_on_change_o_action(GtkToggleAction* paction, gpointer userdata );
    void S_on_change_n_action(GtkToggleAction* paction, gpointer userdata );
    void S_on_change_s_action(GtkToggleAction* paction, gpointer userdata );
    void S_on_change_p_action(GtkToggleAction* paction, gpointer userdata );
    void S_on_change_f_action(GtkToggleAction* paction, gpointer userdata );
    void S_on_change_cl_action(GtkToggleAction* paction, gpointer userdata );
    void S_on_change_single_action(GtkToggleAction* paction, gpointer userdata );
    void S_on_change_double_action(GtkToggleAction* paction, gpointer userdata );
    void S_on_change_triple_action(GtkToggleAction* paction, gpointer userdata );
    void S_on_property_action(GtkToggleAction* paction, gpointer userdata );
    void S_on_style_action(GtkToggleAction* paction, gpointer userdata );
    

    GtkActionEntry menu_actions[] =
    {
        { "ChangeElementAction", NULL, "Change Element"},
        { "ChangeCAction", NULL, "C",
          NULL, "Set to carbon",
          G_CALLBACK(S_on_change_c_action)
        },
        { "ChangeHAction", NULL, "H",
          NULL, "Set to hydrogen",
          G_CALLBACK(S_on_change_h_action)
        },
        { "ChangeOAction", NULL, "O",
          NULL, "Set to oxygen",
          G_CALLBACK(S_on_change_o_action)
        },
        { "ChangeNAction", NULL, "N",
          NULL, "Set to nitrogen",
          G_CALLBACK(S_on_change_n_action)
        },
        { "ChangeSAction", NULL, "S",
          NULL, "Set to sulfur",
          G_CALLBACK(S_on_change_s_action)
        },
        { "ChangePAction", NULL, "P",
          NULL, "Set to sulfur",
          G_CALLBACK(S_on_change_p_action)
        },
        { "ChangeFAction", NULL, "F",
          NULL, "Set to flourine",
          G_CALLBACK(S_on_change_f_action)
        },
        { "ChangeClAction", NULL, "Cl",
          NULL, "Set to chlorine",
          G_CALLBACK(S_on_change_cl_action)
        },
        { "ChangeOrderAction", NULL, "Change Order"},
        { "ChangeSingleAction", NULL, "Single",
          NULL, "Set to single",
          G_CALLBACK(S_on_change_single_action)
        },
        { "ChangeDoubleAction", NULL, "Double",
          NULL, "Set to double",
          G_CALLBACK(S_on_change_double_action)
        },
        { "ChangeTripleAction", NULL, "Triple",
          NULL, "Set to triple",
          G_CALLBACK(S_on_change_triple_action)
        },
        { "PropertyAction", NULL, "Property",
          NULL, "Property",
          G_CALLBACK(S_on_property_action)
        },
        { "StyleAction", NULL, "Display Style...",
          NULL, "Set display style",
          G_CALLBACK(S_on_style_action)
        }
    };
    
    static const int menu_naction = G_N_ELEMENTS(menu_actions);

    gtk_glarea::gtk_glarea( )
    {
        m_self = gtk_drawing_area_new();
        gtk_widget_set_size_request( m_self, 800, 600 );
        gtk_widget_set_events( m_self, 
                               GDK_EXPOSURE_MASK|
                               GDK_POINTER_MOTION_HINT_MASK|
                               GDK_BUTTON_MOTION_MASK|
                               GDK_BUTTON_PRESS_MASK|
                               GDK_BUTTON_RELEASE_MASK|
                               GDK_KEY_PRESS_MASK);

        GdkGLConfigMode mode = GdkGLConfigMode( GDK_GL_MODE_RGB|GDK_GL_MODE_DEPTH|GDK_GL_MODE_DOUBLE );

        static GdkGLConfig* glconfig = gdk_gl_config_new_by_mode( mode );

        if (glconfig == NULL)
        {
            std::cerr << "Warning: Cannot find the double-buffered visual." << std::endl;
            std::cerr << "Warning: Trying single-buffered visual." << std::endl;

            /* Try single-buffered visual */
            glconfig = gdk_gl_config_new_by_mode ( (GdkGLConfigMode)(GDK_GL_MODE_RGB|
                                                                     GDK_GL_MODE_DEPTH) );
            if (glconfig == NULL)
            {
                std::cerr << "Error: No appropriate OpenGL-capable visual found." << std::endl;
                return;
            }
        }

        gtk_widget_set_gl_capability( m_self, glconfig, NULL, TRUE, GDK_GL_RGBA_TYPE );
 
        g_signal_connect_after( G_OBJECT( m_self ), "realize",
                                G_CALLBACK( on_drawing_realize ), this );
        g_signal_connect( G_OBJECT( m_self ), "configure_event",
                          G_CALLBACK( on_drawing_configure ), this );
        g_signal_connect( G_OBJECT( m_self ), "expose_event",
                          G_CALLBACK( on_drawing_expose ), this );
        g_signal_connect( G_OBJECT( m_self ), "button_press_event",
                          G_CALLBACK( on_drawing_button_press ), this );
        g_signal_connect( G_OBJECT( m_self ), "button_release_event",
                          G_CALLBACK( on_drawing_button_release ), this );
        g_signal_connect( G_OBJECT( m_self ), "motion_notify_event",
                          G_CALLBACK( on_drawing_motion_notify ), this );
        g_signal_connect( G_OBJECT( m_self ), "key_press_event",
                          G_CALLBACK( on_drawing_key_press ), this );

        GTK_WIDGET_SET_FLAGS( m_self, GTK_CAN_FOCUS );
        gtk_widget_show( m_self );

        m_actions = gtk_action_group_new("GlareaPopupMenuActionGroup");
        gtk_action_group_add_actions( m_actions, menu_actions, menu_naction, this );

        m_manager = gtk_ui_manager_new();
        gtk_ui_manager_insert_action_group( m_manager, m_actions, 0 );

        if( !gtk_ui_manager_add_ui_from_string(m_manager, menu_info, -1, NULL) )
        {
            throw std::runtime_error( "Error: can not parse ui from string" );
        }

        assert( m_manager != NULL );
        m_popmenu = gtk_ui_manager_get_widget( m_manager, "/popup" );
    }

    gtk_glarea::~gtk_glarea()
    {
    }

    void gtk_glarea::glbegin( )
    {
        GdkGLContext *glcontext = gtk_widget_get_gl_context( m_self );
        GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable( m_self );

        if ( !gdk_gl_drawable_gl_begin( gldrawable, glcontext ) )
        {
            throw logic_error( "can't begin OpenGL routine." );
        }
        
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    }

    void gtk_glarea::expose( )
    {

        GdkRectangle rect;
        rect.x = 0;
        rect.y = 0;
        rect.width = 640;
        rect.height = 480;

        GdkWindow* gdkwin = m_self->window;
        gdk_threads_enter();
        gdk_window_invalidate_rect( gdkwin, &rect, TRUE );
        gdk_flush();
        gdk_threads_leave();
    }

    void gtk_glarea::glend( )
    {        
        GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable( m_self );

        if(gdk_gl_drawable_is_double_buffered ( gldrawable ) )
        {
            gdk_gl_drawable_swap_buffers( gldrawable );
        }
        else
        {
            glFlush();
        }

        gdk_gl_drawable_gl_end( gldrawable);
    }

    void gtk_glarea::popup_menu(GdkEventButton* event)
    {
        assert( m_popmenu != NULL );

        double transmat[16];
        glGetDoublev( GL_MODELVIEW_MATRIX, transmat );
        if( !is_identity(transmat, 4) )
        {
            transform( content(), transmat );
            glMatrixMode( GL_MODELVIEW );
            glLoadIdentity( );
        }

        atomvec_ptr pav = content().get_avec( "selected" );
        bondvec_ptr pbv = content().get_bvec( "selected-bonds" );
        
        double dx = ortho().get_x( event->x );
        double dy = ortho().get_y( event->y );
        
        numvec rect( dx-0.3, dy-0.3, dx+0.3, dy+0.3 );

        atomvec_t a;
        locate_atom( content(), rect, a, FIND_ONE );

        bondvec_t b;
        locate_bond( content(), rect, b, FIND_ONE );

        if( a.size()||b.size()>0 )
        {
            if( pav )
            {
                unselect( *pav );
                pav->clear();
            }
            
            if( pbv )
            {
                unselect( *pbv );
                pbv->clear();
            }
        }

        if( a.size()>0 )
        {
            if( pav==NULL )
            {
                pav = atomvec_ptr( new atomvec_t() );
                content().set( "selected", pav );
            }
            
            a[0].set_i(SELECTED, 1);
            pav->push_back( a[0] );
        }
        

        if( b.size()>0 )
        {
            if( pbv==NULL )
            {
                pbv = bondvec_ptr( new bondvec_t() );
                content().set( "selected-bonds", pbv );
            }
            
            b[0].set_i(SELECTED, 1);
            pbv->push_back( b[0] );
        }

        
        if( a.size()>0 || b.size()>0 )
        {
            glbegin();
            update();
            repaint();
            glend();
        }
        
        gtk_menu_popup( GTK_MENU(m_popmenu), NULL, NULL, 
                        NULL, NULL, event->button, event->time );
        
    }

    void change_elem( gtk_glarea* parea, int elem )
    {
        atomvec_ptr ps = content().get_avec( "selected" );
        if( ps )
        {
            for( int i=0; i < ps->size(); ++i )
	    {
	        moref_t a = ps->at(i);
		a.set_i(SELECTED, 0);
		a.set_i(ELEMENT,  elem);
		a.set_s(NAME,     uniq_name(a) );
		a.set_s(TYPE,     pertab_t::get_symbol(elem) );
	    }

            ps->clear();
        }
            
        
        parea->glbegin();
        parea->update();
        parea->repaint();
        parea->glend();
    }


    void change_order( gtk_glarea* parea, int ord )
    {
        bondvec_ptr ps = content().get_bvec( "selected-bonds" );
        if( ps )
        {
            std::for_each( ps->begin(), ps->end(), iparm_setter1(ORDER, ord) );
            std::for_each( ps->begin(), ps->end(), iparm_setter1(SELECTED, 0) );
            ps->clear();
        }
        
        parea->glbegin();
        parea->update();
        parea->repaint();
        parea->glend();
    }

    void S_on_change_c_action(GtkToggleAction* paction, gpointer user_data )
    {
        gtk_glarea* parea = (gtk_glarea*)user_data;
        change_elem( parea, CARBON );
    }
    
    void S_on_change_h_action(GtkToggleAction* paction, gpointer user_data )
    {
        gtk_glarea* parea = (gtk_glarea*)user_data;
        change_elem( parea, HYDROGEN );
    }
    
    void S_on_change_o_action(GtkToggleAction* paction, gpointer user_data )
    {
        gtk_glarea* parea = (gtk_glarea*)user_data;
        change_elem( parea, OXYGEN );
    }
    
    void S_on_change_n_action(GtkToggleAction* paction, gpointer user_data )
    {
        gtk_glarea* parea = (gtk_glarea*)user_data;
        change_elem( parea, NITROGEN );
    }
    
    void S_on_change_s_action(GtkToggleAction* paction, gpointer user_data )
    {
        gtk_glarea* parea = (gtk_glarea*)user_data;
        change_elem( parea, SULFUR );
    }
    
    void S_on_change_p_action(GtkToggleAction* paction, gpointer user_data )
    {
        gtk_glarea* parea = (gtk_glarea*)user_data;
        change_elem( parea, PHOSPHORUS );
    }
    

    void S_on_change_f_action(GtkToggleAction* paction, gpointer user_data )
    {
        gtk_glarea* parea = (gtk_glarea*)user_data;
        change_elem( parea, FLUORINE );
    }
    
    void S_on_change_cl_action(GtkToggleAction* paction, gpointer user_data )
    {
        gtk_glarea* parea = (gtk_glarea*)user_data;
        change_elem( parea, CHLORINE );
    }
    
    void S_on_change_single_action(GtkToggleAction* paction, gpointer user_data )
    {
        gtk_glarea* parea = (gtk_glarea*)user_data;
        change_order( parea, 1 );
    }
    
    void S_on_change_double_action(GtkToggleAction* paction, gpointer user_data )
    {
        gtk_glarea* parea = (gtk_glarea*)user_data;
        change_order( parea, 2 );
    }
    
    void S_on_change_triple_action(GtkToggleAction* paction, gpointer user_data )
    {
        gtk_glarea* parea = (gtk_glarea*)user_data;
        change_order( parea, 3 );
    }
    
    void S_on_property_action(GtkToggleAction* paction, gpointer user_data )
    {
    }
    
    void S_on_style_action(GtkToggleAction* paction, gpointer user_data )
    {
    }
    

    
} // namespace amber

