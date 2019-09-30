#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#include <drawing.hpp>
#include <console.hpp>
#include <mainwin.hpp>

#ifdef WIN32
#  define EXPORT __declspec(dllexport)
#else
#  define EXPORT
#endif

using namespace mort;

extern "C"
{

EXPORT void on_new1_activate(GtkMenuItem *menuitem, gpointer user_data)
{
}

EXPORT void on_open1_activate(GtkMenuItem *menuitem, gpointer user_data)
{
}

EXPORT void on_save1_activate(GtkMenuItem *menuitem, gpointer user_data)
{
}

EXPORT void on_save_as1_activate(GtkMenuItem *menuitem, gpointer user_data)
{
}

EXPORT void on_cut1_activate(GtkMenuItem *menuitem, gpointer user_data)
{
}

EXPORT void on_copy1_activate(GtkMenuItem *menuitem, gpointer user_data)
{
}

EXPORT void on_paste1_activate(GtkMenuItem *menuitem, gpointer user_data)
{
}

EXPORT void on_delete1_activate(GtkMenuItem *menuitem, gpointer user_data)
{
}

EXPORT void on_about1_activate(GtkMenuItem *menuitem, gpointer user_data)
{
}

EXPORT void on_drawing_realize(GtkWidget *widget, gpointer user_data)
{
    drawing()->init();
}

EXPORT gboolean on_drawing_configure(GtkWidget *widget, GdkEventConfigure *event, gpointer user_data)
{
    drawing()->resize( widget->allocation.width, widget->allocation.height );
    return FALSE;
}

EXPORT gboolean on_drawing_expose( GtkWidget *widget, GdkEventExpose  *event, gpointer user_data)
{
    drawing()->repaint();
    return FALSE;
}

EXPORT gboolean on_drawing_button_press( GtkWidget* widget, GdkEventButton* event, gpointer user_data )
{
    drawing()->mouse_press( event->button, event->state, event->x, event->y );
    return FALSE;
}

EXPORT gboolean on_drawing_button_release( GtkWidget* widget, GdkEventButton* event, gpointer user_data )
{
    drawing()->mouse_release( event->button, event->state, event->x, event->y );
    return FALSE;
}
    
EXPORT gboolean on_drawing_motion_notify( GtkWidget* widget, GdkEventMotion* event, gpointer user_data )
{
    gint x, y;
    
    drawing()->mouse_move( 1, event->state, event->x, event->y );

    gtk_widget_get_pointer( widget, &x, &y );
    
    return FALSE;
}

EXPORT gboolean on_console_key_press_event( GtkWidget *widget, GdkEventKey* event, gpointer user_data )
{
    if( event->keyval == GDK_Return )
    {
        if( console()->on_enter( ) )
        {
            return true;
        }

        gtk_main_quit();
    }
    else if( event->keyval == GDK_Up )
    {
        console()->on_up();
        return true;
    }
    else if( event->keyval == GDK_Down )
    {
        console()->on_down();
        return true;
    }
    else if( event->keyval == GDK_Tab )
    {
        console()->on_tab();
        return true;
    }
    
    return false;
}

EXPORT void main_quit()
{
    gtk_main_quit();
}

} // extern "C"
