#include <GL/gl.h>
#include <GL/glu.h>
#include <object.hpp>
#include <common.hpp>
#include "glutils.hpp"

namespace mort
{
    template< typename T > struct unequal
    {
        unequal( const T& x )
            :m_x( &x )
        {
        }
        
        bool operator()( const T& y ) const
        {
            return y != *m_x;
        }
        
        const T* m_x;
    };  

    pointer_t atom_ref( const pointer_t& ai, const pointer_t& aj, char strategy )
    {     
        atomvec_t ring;
        if( strategy & SPECIAL_RING_BOND )  
        {
            assert( find_ring(get_bond(ai, aj), ring) );
            return ring.back();
        }

        iterator_t nbr = std::find_if( ai.related_atom_begin(), ai.related_atom_end(), unequal< pointer_t >( aj ) );

        BOOST_ASSERT( nbr != ai.related_atom_end() );

        return *nbr;
    }
    
    double get_distance( const pointer_t& bond )
    {
        static double dis[] = { 0.0, 0.2, 0.4, 0.2 };

        int order = bond.get_i(ORDER);

    
        BOOST_ASSERT( order > 0 && order <= 4 );
        
        return dis[ order - 1 ];
    } 
   
    style_t::style_t( char strategy )
    {
        m_strategy = strategy;
    }
    
    style_t::~style_t( )
    {
    }

    void style_t::render_atom( const pointer_t& atom ) const
    {   
        bool single = ( atom.related_atom_number() == 0 );
    
        set_material( atom.get_i(ELEMENT) );
    
        display_atom( atom.get_v(POSITION), single );
    }

    void style_t::render_bond( const pointer_t& bond ) const
    {
        if( m_strategy & DISPLAY_ORDER )
        {
            int order = bond.get_i(ORDER);
            
            if( order == 1 )
            {
                render_single( bond );
            }
            else if( order == 2 )
            {
                render_double( bond );
            }
            else if( order == 3 )
            {
                render_single( bond );
                render_double( bond );
            }
            else if( order == 4 )
            {
                render_double( bond );
            }
        }
        else
        {
            render_single( bond );
        }
    }

    void style_t::render_single( const pointer_t& bond ) const
    {

        pointer_t ai = atom_1st( bond );
        pointer_t aj = atom_2nd( bond );
        
        numvec pos_i = ai.get_v(POSITION);
        numvec pos_j = aj.get_v(POSITION);
        
        int element_i = ai.get_i(ELEMENT);
        int element_j = aj.get_i(ELEMENT);
        
        display_bond( element_i, pos_i, element_j, pos_j, bond.get_i(ORDER) );
    }

    void style_t::render_double( const pointer_t& bond ) const
    {
        double dis = get_distance( bond );

        bool symetrical = true;

        pointer_t ai = atom_1st( bond );
        pointer_t aj = atom_2nd( bond );

        if( ai.related_atom_number() < aj.related_atom_number() )
        {
            std::swap( ai, aj );
        }        

        pointer_t ak = atom_ref( ai, aj, m_strategy );
        
        if( (m_strategy & SPECIAL_RING_BOND) && has_ring(bond) )
        {
            symetrical = false;
        }
        
        double shift = symetrical ? dis / 2.0 : dis;
        double theta = symetrical ? 90.0 : 60.0;
        
        numvec vi = ai.get_v(POSITION);
        numvec vj = aj.get_v(POSITION);
        numvec vk = ak.get_v(POSITION);
        
        numvec vki = vk - vi;
        numvec vji = vj - vi;

        numvec axs = normalize_copy( cross(vki,vji) );

        numvec vni = rotate_copy(vji*shift, axs,  theta);
        numvec vnj = rotate_copy(vji*shift, axs, -theta);
        
        numvec beg_1 = vi + vni;
        numvec end_1 = vj - vnj;

        int element_i = ai.get_i(ELEMENT);
        int element_j = aj.get_i(ELEMENT);
        int order = bond.get_i(ORDER);

        if( order == 4 )
        {
            glLineStipple( 1, 0x00FF );
            glEnable( GL_LINE_STIPPLE );
        }

        display_bond( element_i, beg_1, element_j, end_1, order );
        
        if( order == 4 )
        {
            glDisable( GL_LINE_STIPPLE );
        }

        numvec beg_2 = symetrical ? vi - vni : vi;
        numvec end_2 = symetrical ? vj + vnj : vj;
        
       display_bond( element_i, beg_2, element_j, end_2, order );
    }

    char style_t::get_strategy() const
    {
        return m_strategy;
    }

    void set_material( int element )
    {
        const double ambient_ratio = 0.4, diffuse_ratio=0.5, specular_ratio = 0.2;
        double red = pertab_t::get_red( element );
        double blue = pertab_t::get_blue( element );
        double green = pertab_t::get_green( element );
        float shineness=50.0;
        
        float ambient[4], diffuse[4], specular[4];
        
        ambient[0] = ambient_ratio*red;
        ambient[1] = ambient_ratio*green;
        ambient[2] = ambient_ratio*blue; 
        ambient[3]=1.0;

        diffuse[0] = diffuse_ratio*red; 
        diffuse[1] = diffuse_ratio*green; 
        diffuse[2] = diffuse_ratio*blue; 
        diffuse[3]=1.0;
	
        specular[0]=specular_ratio*red; 
        specular[1]=specular_ratio*green; 
        specular[2]=specular_ratio*blue; 
        specular[3]=1.0;


	glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT, ambient);
	glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE, diffuse);
        glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR, specular);
        glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&shineness);
    }

    void style_t::display_bond( int mat_i, const numvec& pos_i, int mat_j, const numvec& pos_j, int order ) const
    {
        if( ( m_strategy & ELEMENT_COLOR ) && mat_i != mat_j )
        {
            numvec vm = 0.5 * ( pos_i + pos_j );
            
            set_material( mat_i );
            display_bond( pos_i, vm, order );
            
            set_material( mat_j );
            display_bond( pos_j, vm, order);
        }
        else
        {
            set_material( mat_i );
            display_bond( pos_i, pos_j, order );
        }
    }

    line_style::line_style( double width )
        : style_t( ELEMENT_COLOR | DISPLAY_ORDER | SPECIAL_RING_BOND )
    {
        m_width = width;
    }

    line_style::~line_style()
    {
    }

    void line_style::display_atom( const numvec& pos_i, bool single ) const
    {
        if( single )
        {
            glBegin( GL_LINES );
            glVertex3d( pos_i[0] - 0.2, pos_i[1], pos_i[2] );
            glVertex3d( pos_i[0] + 0.2, pos_i[1], pos_i[2] );
            glVertex3d( pos_i[0], pos_i[1] - 0.2, pos_i[2] );
            glVertex3d( pos_i[0], pos_i[1] + 0.2, pos_i[2] );
            glVertex3d( pos_i[0], pos_i[1], pos_i[2] - 0.2 );
            glVertex3d( pos_i[0], pos_i[1], pos_i[2] + 0.2 );
            glEnd();
        }
    }
    
    
    void line_style::display_bond( const numvec& pos_i, const numvec& pos_j, int ) const
    {
        glLineWidth( m_width );

        glBegin( GL_LINES );

        glVertex3d( pos_i[0], pos_i[1], pos_i[2] );
        glVertex3d( pos_j[0], pos_j[1], pos_j[2] );

        glEnd();
    }

    ballstick_style::ballstick_style( double radius, double thickness, char strategy )
        : style_t( strategy )
    {
        m_ball_radius = radius;
        m_ball_splice = 10;
        m_ball_vertex = 5;
        
        m_stick_thickness = thickness;
        m_stick_splice = 20;
        m_stick_vertex = 5;
    }

    ballstick_style::~ballstick_style()
    {
    }

    double ballstick_style::get_thickness( int order ) const
    {
        if( get_strategy() & DISPLAY_ORDER )
        {
            static double scales[] = { 1.0, 1.0 / 3.0, 1.0 / 4.0 };
            
            return m_stick_thickness * scales[ order - 1 ];
        }
        
        return m_stick_thickness;
    }
    
    void ballstick_style::display_atom( const numvec& pos_i, bool ) const
    {
     
        glPushMatrix();

        {        
            glTranslated( pos_i[0], pos_i[1], pos_i[2] );

            GLUquadricObj * q = gluNewQuadric();

            gluSphere( q, m_ball_radius, m_ball_splice, m_ball_vertex );

            gluDeleteQuadric(q);
        }
        
	glPopMatrix();
    }


    void ballstick_style::display_bond( const numvec& pos_i, const numvec& pos_j, int order ) const
    {

        double thickness = get_thickness( order ); 

        numvec vij = pos_j - pos_i;

	double height = normal( vij );

        double rot_angle = 180.0 / M_PI * acos( vij[2] / height );
        
	glPushMatrix();

        glTranslated( pos_i[0], pos_i[1], pos_i[2] );

        glRotated( rot_angle, -vij[1], vij[0], 0.0 );

        GLUquadricObj * q = gluNewQuadric();

        gluCylinder( q, thickness, thickness, height, m_stick_splice, m_stick_vertex );

        gluDeleteQuadric( q );

	glPopMatrix();

    }

    namespace GL
    {        
        void light( )
        {        
            GLfloat position[] = { 0.0, 0.0, 10.0, 0.0 };
            GLfloat ambient[] = { 0.5, 0.5, 0.5, 1.0 };
            GLfloat diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
            GLfloat local_view = 0.0;

  
            glLightfv(GL_LIGHT0, GL_POSITION, position);
            glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse );
            glLightfv(GL_LIGHT0, GL_AMBIENT, ambient );

            glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, &local_view);
    

            glFrontFace(GL_CW);
            glEnable(GL_LIGHTING);
            glEnable(GL_LIGHT0);
        
            glEnable(GL_AUTO_NORMAL);
            glEnable(GL_NORMALIZE);

        }    

        void depth_cue( bool flag )
        {
            if( flag )
            {
                GLfloat color[4] = { 0.0, 0.0, 0.0, 1.0 };

                glEnable( GL_FOG );            
                glEnable( GL_DEPTH_TEST );
                glHint( GL_FOG_HINT, GL_NICEST );
                glFogf( GL_FOG_MODE, GL_LINEAR );
                glFogf( GL_FOG_START, 10.f );
                glFogf( GL_FOG_END,   -40.0f );
                glFogf(GL_FOG_DENSITY, 0.4f);
                glFogfv( GL_FOG_COLOR, color);
                glDepthFunc(GL_LESS);
            }
            else
            {
                glDisable( GL_FOG );
            }
        }

        void anti_alias( bool flag )
        {
            if( flag )
            {
                glEnable(GL_BLEND);
                glEnable(GL_LINE_SMOOTH);
                glEnable(GL_POLYGON_SMOOTH);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
            }
            else
            {
                glDisable( GL_BLEND );
                glDisable( GL_LINE_SMOOTH );
            }
        }

        void new_list( const string& name )
        {
            glNewList( hash( name.c_str() ), GL_COMPILE );
        }
    
        void call_list( const string& name )
        {
            glCallList( hash( name.c_str() ) );
        }

        bool has_list( const string& name )
        {
            return glIsList( hash( name.c_str() ) );
        }
    
        void end_list( )
        {
            glEndList();
        }

    } // namespace GL


} // namespace mort


