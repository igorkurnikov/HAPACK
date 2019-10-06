#ifndef MOLVIEW_MOUSE_H
#define MOLVIEW_MOUSE_H

#include <common.hpp>

namespace mort
{

    class mouse_i
    {
    public: 

        mouse_i();
        
        mouse_i( const string& name );
        
        virtual ~mouse_i();

        virtual void move( int button, double x, double y ) = 0;
        
        virtual void release( int button, double x, double y ) = 0;

        virtual shared_ptr< mouse_i > clone( double x, double y ) const = 0;
        
        static void add( const string& name, const mouse_i* ptr );
        
        static shared_ptr< mouse_i > get( const string& name, double x, double y );
        
    };

    class scale_mouse : public mouse_i
    {
    public:

        scale_mouse();

        scale_mouse( double startx, double starty );
    
        virtual ~scale_mouse( );
    
        virtual void move( int button, double x, double y );
    
        virtual void release( int button, double x, double y );
    
        virtual shared_ptr< mouse_i > clone( double x, double y ) const;
        
    private:

        double m_modelview[16];
    
        double m_startx;
    
        double m_starty;
    };

    class translate_mouse : public mouse_i
    {
    public:

        translate_mouse();

        translate_mouse( double x, double y );

        virtual ~translate_mouse();
    
        virtual void move( int button, double x, double y );
    
        virtual void release( int button, double x, double y );

        virtual shared_ptr< mouse_i > clone( double x, double y ) const;        
    
    private:

        double m_modelview[16];

        double m_startx;
    
        double m_starty;
    
    };

    class rotate_mouse : public mouse_i
    {
    public:

        rotate_mouse();

        rotate_mouse( double startx, double starty );
    
        virtual ~rotate_mouse( );
    
        virtual void move( int button, double x, double y );
    
        virtual void release( int button, double x, double y );
    
        virtual shared_ptr< mouse_i > clone( double x, double y ) const;
        
    private:

        double m_modelview[16];
    
        double m_startx;
    
        double m_starty;
    };

    
} // namespace mort


#endif

