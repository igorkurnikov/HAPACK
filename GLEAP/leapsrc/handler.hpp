#include <gtk/gtk.h>

#ifdef WIN32
#  define EXPORT __declspec(dllexport)
#else
#  define EXPORT
#endif

extern "C"
{
    EXPORT void on_drawing_realize(GtkWidget *widget, gpointer user_data);

    EXPORT gboolean on_drawing_configure(GtkWidget *widget, GdkEventConfigure *event, gpointer user_data);

    EXPORT gboolean on_drawing_expose( GtkWidget *widget, GdkEventExpose  *event, gpointer user_data);
    
    EXPORT gboolean on_drawing_button_press( GtkWidget* widget, GdkEventButton* event, gpointer user_data );
    
    EXPORT gboolean on_drawing_button_release( GtkWidget* widget, GdkEventButton* event, gpointer user_data );
    
    EXPORT gboolean on_drawing_motion_notify( GtkWidget* widget, GdkEventMotion* event, gpointer user_data );
}

