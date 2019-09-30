#include <GL/gl.h>
#include <gtk/gtk.h>
#include <gtk/gtkgl.h>
#include <leaplog.hpp>

#include "gtkarea.hpp"

namespace amber
{
    using namespace mort;

    gtk_glarea::gtk_glarea( GtkWidget* widget )
    {
        m_widget = widget;

        GdkGLConfigMode mode = GdkGLConfigMode( GDK_GL_MODE_RGB|GDK_GL_MODE_DEPTH|GDK_GL_MODE_DOUBLE );

        static GdkGLConfig* glconfig = gdk_gl_config_new_by_mode( mode );

        if (glconfig == NULL)
        {
            leaplog_t::putline("*** Cannot find the double-buffered visual.\n");
            leaplog_t::putline("*** Trying single-buffered visual.\n");

            /* Try single-buffered visual */
            glconfig = gdk_gl_config_new_by_mode ( (GdkGLConfigMode)(GDK_GL_MODE_RGB   |
                                                                         GDK_GL_MODE_DEPTH) );
            if (glconfig == NULL)
            {
                leaplog_t::putline("*** No appropriate OpenGL-capable visual found.\n");
                return;
            }
        }

        gtk_widget_set_gl_capability( m_widget, glconfig, NULL, TRUE, GDK_GL_RGBA_TYPE );
    }

    gtk_glarea::~gtk_glarea()
    {
    }

    void gtk_glarea::glbegin( )
    {
        GdkGLContext *glcontext = gtk_widget_get_gl_context( m_widget );
        GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable( m_widget );

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

        GdkWindow* gdkwin = m_widget->window;
        gdk_threads_enter();
        gdk_window_invalidate_rect( gdkwin, &rect, TRUE );
        gdk_flush();
        gdk_threads_leave();
    }

    void gtk_glarea::glend( )
    {        
        GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable( m_widget );

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
    
} // namespace amber

