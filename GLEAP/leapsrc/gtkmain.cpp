#include <object.hpp>
#include <drawing.hpp>
#include <mainwin.hpp>
#include <control.hpp>
#include <console.hpp>
#include <leaplog.hpp>

#include <GL/gl.h>
#include <gtk/gtkgl.h>

#include "handler.hpp"
#include "gtkarea.hpp"
#include "strbuff.hpp"

using namespace std;

using namespace mort;

using namespace amber;

void* console_thread( void* arg )
{
    while( console()->on_enter() );

    gtk_main_quit();

    return NULL;
}

void gl_init()
{
    GtkWidget* window = gtk_window_new( GTK_WINDOW_TOPLEVEL );
    GtkWidget* widget = gtk_drawing_area_new();

    gtk_window_set_title( GTK_WINDOW( window ), "gLEaP" );
    gtk_window_set_keep_above( GTK_WINDOW( window ), 1 );
    gtk_window_move( GTK_WINDOW( window ), 20, 20 );
    gtk_container_set_border_width( GTK_CONTAINER( window ), 0 );
    gtk_container_add( GTK_CONTAINER( window ), widget );
    gtk_widget_set_size_request( widget, 640, 480 );
    gtk_widget_set_events( widget, GDK_EXPOSURE_MASK|
                                   GDK_POINTER_MOTION_HINT_MASK|
                                   GDK_BUTTON_MOTION_MASK|
                                   GDK_BUTTON_PRESS_MASK|
                                   GDK_BUTTON_RELEASE_MASK);
 
    drawing() = shared_ptr< drawing_i >( new gtk_glarea( widget ) ); 

    g_signal_connect_after( G_OBJECT( widget ), "realize",
                            G_CALLBACK( on_drawing_realize ), NULL );
    g_signal_connect( G_OBJECT( widget ), "configure_event",
                      G_CALLBACK( on_drawing_configure ), NULL );
    g_signal_connect( G_OBJECT( widget ), "expose_event",
                      G_CALLBACK( on_drawing_expose ), NULL );
    g_signal_connect( G_OBJECT( widget ), "button_press_event",
                      G_CALLBACK( on_drawing_button_press ), NULL );
    g_signal_connect( G_OBJECT( widget ), "button_release_event",
                      G_CALLBACK( on_drawing_button_release ), NULL );
    g_signal_connect( G_OBJECT( widget ), "motion_notify_event",
                      G_CALLBACK( on_drawing_motion_notify ), NULL );
    g_signal_connect( G_OBJECT( window ), "delete_event",
                      G_CALLBACK( gtk_main_quit ), NULL );

    gtk_widget_show( widget );
    gtk_widget_show( window );
} 



int main( int argc, char** argv )
{
    const char* amberhome = getenv( "AMBERHOME" );
    
    if( amberhome == NULL )
    {
        cerr << "Environmental variable AMBERHOME must be set." << std::endl;
        return -1;
    }

    g_thread_init( NULL );

    gdk_threads_init();

    gtk_init( &argc, &argv );

    gtk_gl_init( &argc, &argv );

    gl_init();

    std::cout << "gl version   : " << glGetString( GL_VERSION ) << std::endl;
    std::cout << "gl vendor    : " << glGetString( GL_VENDOR )  << std::endl;
    std::cout << "gl renderer  : " << glGetString( GL_RENDERER ) << std::endl;
    std::cout << "gl extensions: " << glGetString( GL_EXTENSIONS ) << std::endl;



    shared_ptr< txtbuff_i > buff( new str_buff( cin, cout ) );

    console() = shared_ptr< console_t >( new console_t( buff ) );

    string path = "./:";

    path += amberhome + string( "/src/gleap/leapdat:" );

    path += amberhome + string( "/dat/leap/prep:" );
    
    path += amberhome + string( "/dat/leap/parm:" );

    path += amberhome + string( "/dat/leap/lib:" );

    path += amberhome + string( "/dat/leap/cmd"  );
    
    mortenv().set( "path", path.c_str() );

    GError* error = NULL;

    GThread* thread = g_thread_create( console_thread, NULL, FALSE, &error );

    if( thread == NULL )
    {
        cerr << "Error: can not create console thread" << std::endl;
        return -1;
    }

    gdk_threads_enter();

    gtk_main();

    gdk_threads_leave();

    return 0;
}

    
#include "../plugins/help.hpp"

amber::help_command g_help_command;


