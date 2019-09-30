#ifndef MORT_GLUTIL_HPP
#define MORT_GLUTIL_HPP

#include <object.hpp>
#include "drawing.hpp"

namespace mort
{
    void set_material( int matid );

	class moref_t;

    /// molecular display style base class
    /// 
    class style_t
    {
      public:

        /// initial
        style_t( char strategy );

        virtual ~style_t();

        char get_strategy() const;

        void render_atom( const moref_t& atom ) const;
        
        void render_bond( const moref_t& bond ) const;

        virtual void display_atom( const numvec& pos_i, bool single ) const = 0;
        
        virtual void display_bond( const numvec& pos_i, const numvec& pos_j, int order ) const = 0;

      private:


        void render_single(const moref_t& bond) const;
        
        void render_double(const moref_t& bond) const;
        
        void display_bond(int element_i, const numvec& pos_i, int element_j, const numvec& pos_j, int order) const;

      private:

        char m_strategy;

    };

    void render(const molecule_t& mol, const style_t* style);

    static const char DISPLAY_ORDER = 0x01;
    
    static const char ELEMENT_COLOR = 0x02;
    
    static const char SPECIAL_RING_BOND = 0x04;

    class line_style : public style_t
    {
      public:

        line_style( double width );
        
        virtual ~line_style();
        
        virtual void display_atom( const numvec& pos_i, bool single ) const;
        
        virtual void display_bond( const numvec& pos_i, const numvec& pos_j, int order ) const;

      private:

        double m_width;
    };

    class ballstick_style : public style_t
    {
      public:

        ballstick_style( double radius, double thickness, char strategy );
        
        ~ballstick_style();
        
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

    };

    namespace GL
    {
        void light();
        
        void anti_alias( bool flag );
        
        void depth_cue( bool flag );

        void new_list( const string& name );
    
        bool has_list( const string& name );

        void end_list( );
    
        void call_list( const string& name );
    
    } // namespace GL
 
    class mouse_i;

    struct ortho_t;

    class gldrawing_t : public drawing_i
    {
    public:

        using drawing_i::add;

        /// constructor
        gldrawing_t();

        /// deconstructor
        virtual ~gldrawing_t();


        virtual void add(const string& type, molecule_t& mol);

        /// handler for mouse pressed
        virtual void mouse_press( int button, int state, double x, double y );

        /// handler for mouse move
        virtual void mouse_move ( int button, int state, double x, double y );

        /// handler for mouse released
        virtual void mouse_release( int button, int state, double x, double y );

        /// OpenGL initializaiton
        virtual void init();
 
        /// paint all graphics
        virtual void repaint();
    
        /// adjust size of the display
        virtual void resize( int width, int height );

    private:
        shared_ptr< mouse_i > mp_mouse;
        shared_ptr< ortho_t > mp_ortho;
        quantity_e m_quantity;
    };


} // namespace mort

#endif

