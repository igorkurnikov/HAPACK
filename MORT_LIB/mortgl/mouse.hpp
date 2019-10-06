#ifndef MOLVIEW_MOUSER_H
#define MOLVIEW_MOUSER_H

#include <common.hpp>
#include <object.hpp>

namespace mort
{

    class mouser_i
    {
    public: 

        mouser_i();
        
        mouser_i( const string& name );
        
        virtual ~mouser_i();

        virtual void move( int button, double x, double y ) = 0;
        
        virtual void release( int button, double x, double y ) = 0;

        virtual shared_ptr< mouser_i > clone( double x, double y ) const = 0;
        
        static void add( const string& name, const mouser_i* ptr );
        
        static shared_ptr< mouser_i > get( const string& name, double x, double y );
        
    };

    typedef shared_ptr< mouser_i > mouser_ptr;

    class scale_mouser : public mouser_i
    {
    public:

        scale_mouser();

        scale_mouser( double startx, double starty );
    
        virtual ~scale_mouser( );
    
        virtual void move( int button, double x, double y );
    
        virtual void release( int button, double x, double y );
    
        virtual mouser_ptr clone( double x, double y ) const;
        
    private:

        double m_modelview[16];
    
        double m_startx;
    
        double m_starty;

        double m_scale;
    };

    class translate_mouser : public mouser_i
    {
    public:

        translate_mouser();

        translate_mouser( double x, double y );

        virtual ~translate_mouser();
    
        virtual void move( int button, double x, double y );
    
        virtual void release( int button, double x, double y );

        virtual mouser_ptr clone( double x, double y ) const;        
    
    private:

        double m_modelview[16];

        double m_startx;
    
        double m_starty;
    
    };

    class rotate_mouser : public mouser_i
    {
    public:

        rotate_mouser();

        rotate_mouser( double startx, double starty );
    
        virtual ~rotate_mouser( );
    
        virtual void move( int button, double x, double y );
    
        virtual void release( int button, double x, double y );
    
        virtual mouser_ptr clone( double x, double y ) const;
        
    private:

        double m_modelview[16];
    
        double m_startx;
    
        double m_starty;
    };

    class select_mouser : public mouser_i
    {
    public:

        select_mouser( );

        select_mouser( double startx, double starty );
        
        virtual ~select_mouser( );

        virtual void move( int button, double x, double y );
        
        virtual void release( int button, double x, double y );
        
        virtual mouser_ptr clone( double x, double y ) const;
        
    private:

        atomvec_ptr m_pselected;
        double m_startx;
        double m_starty;
    };
    
        
    class drawbond_mouser : public mouser_i
    {
    public:

        drawbond_mouser( );

        drawbond_mouser( double startx, double starty );
        
        virtual ~drawbond_mouser( );

        virtual void move( int button, double x, double y );
        
        virtual void release( int button, double x, double y );
        
        virtual mouser_ptr clone( double x, double y ) const;
        
    private:

        void init_frag( molecule_t& m );

        void init_2d( double startx, double starty );
        
        void init_3d( double startx, double starty );
        
    private:

        double m_startx;
        double m_starty;
        double m_curth;

        atomvec_t m_atm1;
        atomvec_t m_atm2;
        bondvec_t m_bond;
        
        
        int m_dim;
        molecule_t m_frag;
        molecule_t* m_pmol;
    };


    
} // namespace mort


#endif

