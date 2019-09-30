#ifndef GTKLEAP_ORTHO_HPP
#define GTKLEAP_ORTHO_HPP

namespace mort
{

    struct ortho_t
    {
        ortho_t( double s );
    
        ~ortho_t();
    
        void resize( int width, int height );

        double get_x( double winx );

        double get_y( double winy );

        double scale;

        double height;

        double width;
    
        double front;
    
        double back;
    };
    
} // namespace mort


#endif


