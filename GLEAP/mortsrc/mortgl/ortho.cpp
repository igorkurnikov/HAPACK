#include <iostream>

#include <GL/gl.h>
#include <GL/glu.h>

#include "ortho.hpp"

namespace mort
{

    ortho_t::ortho_t( double s )
    {
        m_scale = s;
        m_front = 50;
        m_back  =-50;
    }

    ortho_t::~ortho_t()
    {
    }

    void ortho_t::resize( int iwidth, int iheight )
    {
        m_width = iwidth;
        m_height = iheight;

        glViewport(0, 0, iwidth, iheight );

        update();        
    }

    void ortho_t::set_scale( double scale )
    {
        m_scale = scale;
        update();
    }    

    double ortho_t::get_scale( ) const
    {
        return m_scale;
    }

    double ortho_t::get_x( double winx ) const
    {
        return m_scale * ( winx - m_width/2.0 );
    }

    double ortho_t::get_y( double winy ) const
    {
        return m_scale * ( m_height/2.0 - winy );
    }


    void ortho_t::update()
    {

        glMatrixMode(GL_PROJECTION);

        glLoadIdentity();

        double xborder = 0.5 * m_scale * m_width;
        double yborder = 0.5 * m_scale * m_height;

        glOrtho( -xborder, xborder, -yborder, yborder, -m_front, -m_back );

        glMatrixMode(GL_MODELVIEW);
    }


} // namespace mort

