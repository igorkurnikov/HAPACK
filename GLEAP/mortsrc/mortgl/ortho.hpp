#ifndef GTKLEAP_ORTHO_HPP
#define GTKLEAP_ORTHO_HPP

namespace mort
{

    class ortho_t
    {
    public:
        
        ortho_t( double s );
    
        ~ortho_t();
    
        void resize( int width, int height );

        void set_scale( double ratio );
        
        double get_scale( ) const;

        double get_x( double winx ) const;

        double get_y( double winy ) const;

    private:

        void update();
        
    private:

        double m_scale;

        double m_height;

        double m_width;
    
        double m_front;
    
        double m_back;
    };
    
} // namespace mort


#endif


