#include <GL/gl.h>
#include <GL/glu.h>

#include <map>
#include "mouse.hpp"

namespace mort
{
    using std::map;

    static map< string, const mouse_i* >& mouses()
    {
        static map< string, const mouse_i* > holder;
        return holder;
    }    

    mouse_i::mouse_i()
    {
    }
    
    mouse_i::mouse_i( const string& name )
    {
        add( name, this );
    }
    
    mouse_i::~mouse_i()
    {
    }    

    void mouse_i::add( const string& name, const mouse_i* ptr )
    {
        mouses().insert( make_pair( name, ptr ) );
    }
    
    shared_ptr< mouse_i > mouse_i::get( const string& name, double x, double y )
    {
        if( mouses().find( name ) == mouses().end() )
        {
            throw logic_error( "Error: can't find mouse handler named " + name );
        }
        
        return mouses().find( name )->second->clone( x, y );
    }


    scale_mouse::scale_mouse()
        :mouse_i( "scale" )
    {
    }
    
    scale_mouse::scale_mouse( double startx, double starty )
    {
        m_startx = startx;
        m_starty = starty;
        glGetDoublev( GL_MODELVIEW_MATRIX, m_modelview );
    }

    scale_mouse::~scale_mouse()
    {
    }

    void scale_mouse::move( int button, double x, double y )
    {
        double dy = ( y - m_starty ) * 0.1;

        double scale = dy > 0 ? 1.0 / ( 1.0 + dy ) : 1.0 - dy;

        double transform[] = {
            scale, 0.0, 0.0, 0.0,
            0.0, scale, 0.0, 0.0,
            0.0, 0.0, scale, 0.0,
            0.0, 0.0, 0.0, 1.0
        };

        glMatrixMode( GL_MODELVIEW );
    
        glLoadMatrixd( transform );
        glMultMatrixd( m_modelview );
    
    
//    glScaled( dy, 1.0, 0.0, 0.0 );
//    glScaled( dx, 0.0, 1.0, 0.0 );
    }

    void scale_mouse::release( int button, double x, double y )
    {
        move( button, x, y );
    }

    shared_ptr< mouse_i > scale_mouse::clone( double x, double y ) const
    {
        return shared_ptr< mouse_i >( new scale_mouse( x, y ) );
    }

    translate_mouse::translate_mouse()
        :mouse_i( "translate" )
    {
    }
    
    translate_mouse::translate_mouse( double startx, double starty )
    {
        m_startx = startx;
        m_starty = starty;
        glGetDoublev( GL_MODELVIEW_MATRIX, m_modelview );
    }

    translate_mouse::~translate_mouse()
    {
    
    }

    void translate_mouse::move( int button, double x, double y )
    {
        const double transform[] = {
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            x - m_startx,
            m_starty - y, 
            0.0, 1.0 
        };

        glLoadMatrixd( transform );
    
        glMultMatrixd( m_modelview );
    

//    double result[16];
    
//    result[0] = t[0] * m[0] + t[4] * m[1] + t[8] * m[2] + m12 * m[3];
//    result[1] = t[1] * m[

    


//    glTranslated( m_scale * ( x - m_startx ), 
//                  m_scale * ( m_starty - y ), 
//                  0.0 );

    }

    void translate_mouse::release( int button, double x, double y )
    {
        move( button, x, y );
    }

    shared_ptr< mouse_i > translate_mouse::clone( double x, double y ) const
    {
        return shared_ptr< mouse_i >( new translate_mouse( x, y ) );
    }
    

    rotate_mouse::rotate_mouse()
        :mouse_i( "rotate" )
    {
    }

    rotate_mouse::rotate_mouse( double startx, double starty )
    {
        m_startx = startx;
        m_starty = starty;
        glGetDoublev( GL_MODELVIEW_MATRIX, m_modelview );
    }

    rotate_mouse::~rotate_mouse()
    {
    }

    void rotate_mouse::move( int button, double x, double y )
    {
        double dx = ( x - m_startx ) / 10.0;
        double dy = ( y - m_starty ) / 10.0;


        const double transform[] = {
            cos(dx), sin(dx)*sin(dy), -sin(dx) * cos(dy), 0.0,
            0.0,     cos(dy),         sin(dy),            0.0,
            sin(dx), -cos(dx)*sin(dy), cos(dx) * cos(dy), 0.0,
            0.0,         0.0,             0.0,            1.0 
        };

        glMatrixMode( GL_MODELVIEW );
    
        glLoadMatrixd( transform );
        glMultMatrixd( m_modelview );
    
    
//    glRotated( dy, 1.0, 0.0, 0.0 );
//    glRotated( dx, 0.0, 1.0, 0.0 );
    }

    void rotate_mouse::release( int button, double x, double y )
    {
        move( button, x, y );
    }

    shared_ptr< mouse_i > rotate_mouse::clone( double x, double y ) const
    {
        return shared_ptr< mouse_i >( new rotate_mouse( x, y ) );
    }

    
} // namespace mort

mort::scale_mouse g_scale_mouse;

mort::rotate_mouse g_rotate_mouse;

mort::translate_mouse g_translate_mouse;

