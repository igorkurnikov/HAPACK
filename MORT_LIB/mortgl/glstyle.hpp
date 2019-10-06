#ifndef GLEAP_MORTSRC_GLSTYLE_HPP
#define GLEAP_MORTSRC_GLSTYLE_HPP

#include <GL/gl.h>
#include <GL/glu.h>
#include <object.hpp>

namespace mort
{

    class glstyle_t
    {
    public:

        /// initial
        glstyle_t( char strategy );

        virtual ~glstyle_t();

        char get_strategy() const;

        void render_atom( const moref_t& atom ) const;
        
        void render_bond( const moref_t& bond ) const;

        virtual void prepare( ) const = 0;

        virtual void display_atom( const numvec& pos_i, bool single ) const = 0;
        
        virtual void display_bond( const numvec& pos_i, const numvec& pos_j, int order ) const = 0;

      private:

        void render_single( const moref_t& bond ) const;
        
        void render_double( const moref_t& bond ) const;
        
        void display_bond( int element_i, const numvec& pos_i, int elmenet_j, const numvec& pos_j, int order ) const;

      private:

        char m_strategy;

    };

    typedef shared_ptr< glstyle_t > glstyle_ptr;

    static const char DISPLAY_ORDER = 0x01;
    
    static const char ELEMENT_COLOR = 0x02;
    
    static const char SPECIAL_RING_BOND = 0x04;

    class line_style : public glstyle_t
    {
    public:

        line_style( double width );
        
        virtual ~line_style();

        virtual void prepare() const;
        
        virtual void display_atom( const numvec& pos_i, bool single ) const;
        
        virtual void display_bond( const numvec& pos_i, const numvec& pos_j, int order ) const;

      private:

        double m_width;
    };

    class ballstick_style : public glstyle_t
    {
      public:

        ballstick_style( double radius, double thickness, char strategy );
        
        ~ballstick_style();

        virtual void prepare( ) const;
        
        virtual void display_atom( const numvec& pos_i, bool single ) const;
        
        virtual void display_bond( const numvec& pos_i, const numvec& pos_j, int order ) const;

      private:

        double get_thickness( int order ) const;
        
      private:

        double m_ball_radius;
        int m_ball_splice;
        int m_ball_vertex;
        
        double m_stick_thickness;
        int m_stick_splice;
        int m_stick_vertex;
        GLUquadric* m_quadric;
    };


} // namespace mort
 

#endif
