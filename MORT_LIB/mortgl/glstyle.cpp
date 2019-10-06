#include "glutils.hpp"
#include "glstyle.hpp"



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

    moref_t atom_ref( const moref_t& ai, const moref_t& aj, char strategy )
    {     
        atomvec_t ring;
//        if( strategy & SPECIAL_RING_BOND )  
//        {
            if( find_ring(get_bond(ai, aj), ring) )
                return ring.back();
//        }

        moiter_t nbr = std::find_if( ai.related_atom_begin(), ai.related_atom_end(), unequal< moref_t >( aj ) );

        BOOST_ASSERT( nbr != ai.related_atom_end() );

        return *nbr;
    }
    
    double get_distance( const moref_t& bond )
    {
        static double dis[] = { 0.0, 0.2, 0.4, 0.2 };

        int order = bond.get_i(ORDER);

    
        BOOST_ASSERT( order > 0 && order <= 4 );
        
        return dis[ order - 1 ];
    } 
   
    glstyle_t::glstyle_t( char strategy )
    {
        m_strategy = strategy;
    }
    
    glstyle_t::~glstyle_t( )
    {
    }

    void glstyle_t::render_atom( const moref_t& atom ) const
    {   
        bool single = ( atom.related_atom_number() == 0 );

        int selected=0;
        if( atom.get_i(SELECTED,selected) && selected>0 )
        {
            GL::set_material( int(SELECTED) );
        }
        else
        {
            GL::set_material( atom.get_i(ELEMENT) );
        }
    
        display_atom( atom.get_v(POSITION), single );
    }

    void glstyle_t::render_bond( const moref_t& bond ) const
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

    void glstyle_t::render_single( const moref_t& bond ) const
    {

        moref_t ai = atom_1st( bond );
        moref_t aj = atom_2nd( bond );
        
        numvec pos_i = ai.get_v(POSITION);
        numvec pos_j = aj.get_v(POSITION);
        
        int material_i = ai.get_i(ELEMENT);
        int material_j = aj.get_i(ELEMENT);

        int selected=0;

        if( ai.get_i(SELECTED, selected) && selected>0 )
            material_i = int(SELECTED);
        
        if( aj.get_i(SELECTED, selected) && selected>0 )
            material_j = int(SELECTED);

        if( bond.get_i(SELECTED, selected) && selected>0 )
        {
            material_i = int(SELECTED);
            material_j = int(SELECTED);
        }
        
        display_bond( material_i, pos_i, material_j, pos_j, bond.get_i(ORDER) );
    }

    void glstyle_t::render_double( const moref_t& bond ) const
    {
        double dis = get_distance( bond );

        bool symetrical = true;

        moref_t ai = atom_1st( bond );
        moref_t aj = atom_2nd( bond );

        if( ai.related_atom_number() < aj.related_atom_number() )
        {
            std::swap( ai, aj );
        }

        numvec vi = ai.get_v(POSITION);
        numvec vj = aj.get_v(POSITION);
        numvec vk(3);


        if( ai.related_atom_number()>1 )
        {
            moref_t ak = atom_ref( ai, aj, m_strategy );
            vk = ak.get_v(POSITION);
        }
        else
        {
            vk = numvec( vi[0]+1.0, vi[1], vi[2] );
        }
        
        
        if( (m_strategy & SPECIAL_RING_BOND) && has_ring(bond) )
        {
            symetrical = false;
        }
        
        double shift = symetrical ? dis / 2.0 : dis;
        double theta = symetrical ? 90.0 : 60.0;
        
        numvec vki = vk - vi;
        numvec vji = vj - vi;

        numvec axs = normalize_copy( cross(vki,vji) );

        numvec vni = rotate_copy(vji*shift, axs, -theta);
        numvec vnj = rotate_copy(vji*shift, axs,  theta);
        
        numvec beg_1 = vi + vni;
        numvec end_1 = vj - vnj;

        int material_i = ai.get_i(ELEMENT);
        int material_j = aj.get_i(ELEMENT);
        int selected=0;

        if( ai.get_i(SELECTED, selected) && selected>0 )
            material_i = int(SELECTED);
        
        if( aj.get_i(SELECTED, selected) && selected>0 )
            material_j = int(SELECTED);

        if( bond.get_i(SELECTED, selected) && selected>0 )
        {
            material_i = int(SELECTED);
            material_j = int(SELECTED);
        }

        int order = bond.get_i(ORDER);

        if( order == 4 )
        {
            glLineStipple( 1, 0x00FF );
            glEnable( GL_LINE_STIPPLE );
        }

        display_bond( material_i, beg_1, material_j, end_1, order );
        
        if( order == 4 )
        {
            glDisable( GL_LINE_STIPPLE );
        }

        numvec beg_2 = symetrical ? vi - vni : vi;
        numvec end_2 = symetrical ? vj + vnj : vj;
        
       display_bond( material_i, beg_2, material_j, end_2, order );
    }

    char glstyle_t::get_strategy() const
    {
        return m_strategy;
    }

    void glstyle_t::display_bond( int mat_i, const numvec& pos_i, int mat_j, const numvec& pos_j, int order ) const
    {
        if( ( m_strategy & ELEMENT_COLOR ) && mat_i != mat_j )
        {
            numvec vm = 0.5 * ( pos_i + pos_j );
            
            GL::set_material( mat_i );
            display_bond( pos_i, vm, order );
            
            GL::set_material( mat_j );
            display_bond( pos_j, vm, order);
        }
        else
        {
            GL::set_material( mat_i );
            display_bond( pos_i, pos_j, order );
        }
    }

    line_style::line_style( double width )
        : glstyle_t( ELEMENT_COLOR | DISPLAY_ORDER | SPECIAL_RING_BOND )
    {
        m_width = width;
    }

    line_style::~line_style()
    {
    }

    void line_style::prepare( ) const
    {
        GL::anti_alias( true );
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
        : glstyle_t( strategy | DISPLAY_ORDER )
    {
        m_ball_radius = radius;
        m_ball_splice = 10;
        m_ball_vertex =  5;
        
        m_stick_thickness = thickness;
        m_stick_splice = 10;
        m_stick_vertex =  5;

        m_quadric = gluNewQuadric();
    }

    ballstick_style::~ballstick_style()
    {
        gluDeleteQuadric( m_quadric );
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

    void ballstick_style::prepare() const
    {
        GL::anti_alias( false );
    }
    
    void ballstick_style::display_atom( const numvec& pos_i, bool ) const
    {
        glPushMatrix();
        {        
            glTranslated( pos_i[0], pos_i[1], pos_i[2] );
            gluSphere( m_quadric, m_ball_radius, m_ball_splice, m_ball_vertex );
        }
	glPopMatrix();
    }


    void ballstick_style::display_bond( const numvec& pos_i, const numvec& pos_j, int order ) const
    {
        numvec vij = pos_j - pos_i;
        double thickness = get_thickness( order );
	double height = normal( vij );
        double rot_angle = 180.0 / M_PI * acos( vij[2] / height );
        
	glPushMatrix();
        {
            glTranslated( pos_i[0], pos_i[1], pos_i[2] );
            glRotated( rot_angle, -vij[1], vij[0], 0.0 );
            gluCylinder( m_quadric, thickness, thickness, height, m_stick_splice, m_stick_vertex );
        }
	glPopMatrix();
    }

} // namespace mort
 
