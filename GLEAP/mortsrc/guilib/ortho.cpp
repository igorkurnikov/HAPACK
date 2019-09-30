

#include <GL/gl.h>
#include <GL/glu.h>

#include "ortho.hpp"

namespace mort
{

    ortho_t::ortho_t( double s )
    {
        scale = s;
        front = 50;
        back  =-50;
    }

    ortho_t::~ortho_t()
    {
    }

    void ortho_t::resize( int iwidth, int iheight )
    {
        width = iwidth;
        height = iheight;

        glViewport(0, 0, iwidth, iheight );

        glMatrixMode(GL_PROJECTION);

        glLoadIdentity();

        double xborder = 0.5 * scale * width;
        double yborder = 0.5 * scale * height;

        glOrtho( -xborder, xborder, -yborder, yborder, -front, -back );

        glMatrixMode(GL_MODELVIEW);
    }

    double ortho_t::get_x( double winx )
    {
        return scale * ( winx - width / 2.0 );
    }

    double ortho_t::get_y( double winy )
    {
        return scale * ( winy - height / 2.0 );
    }

} // namespace mort

