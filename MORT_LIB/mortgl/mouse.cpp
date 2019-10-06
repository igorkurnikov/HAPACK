#include <GL/gl.h>
#include <GL/glu.h>

#include <map>
#include <guilib.hpp>

#include "ortho.hpp"
#include "mouse.hpp"
#include "glutils.hpp"
#include "drawing.hpp"

namespace mort
{

    void print_transmat( const double* mat )
    {
	std::cout << "transmat:" << std::endl;
        for(int i=0; i < 4; ++i)
	{
	    for( int j=0; j < 4; ++j ) 
	    {
                std::cout << format("%10.3f") % mat[4*i+j] << " ";
	    }
	    std::cout << "\n";
	}
    }
    

    using std::map;

    static map< string, const mouser_i* >& mousers()
    {
        static map< string, const mouser_i* > holder;
        return holder;
    }    

    mouser_i::mouser_i()
    {
    }
    
    mouser_i::mouser_i( const string& name )
    {
        add( name, this );
    }
    
    mouser_i::~mouser_i()
    {
    }    

    void mouser_i::add( const string& name, const mouser_i* ptr )
    {
        mousers().insert( make_pair( name, ptr ) );
    }
    
    shared_ptr< mouser_i > mouser_i::get( const string& name, double x, double y )
    {
        if( mousers().find( name ) == mousers().end() )
        {
            throw logic_error( "Error: can't find mouser handler named " + name );
        }
        
        return mousers().find( name )->second->clone( x, y );
    }


    scale_mouser::scale_mouser()
        :mouser_i( "scale" )
    {
    }
    
    scale_mouser::scale_mouser( double startx, double starty )
    {
        m_startx = startx;
        m_starty = starty;
        glGetDoublev( GL_MODELVIEW_MATRIX, m_modelview );
        m_scale = drawing()->ortho().get_scale();
    }

    scale_mouser::~scale_mouser()
    {
    }

    void scale_mouser::move( int button, double x, double y )
    {
        double dy = ( m_starty - y ) * 0.1;

        double ratio = dy > 0 ? 1.0 / ( 1.0 + dy ) : 1.0 - dy;
        
        drawing()->ortho().set_scale( m_scale/ratio );
    }

    void scale_mouser::release( int button, double x, double y )
    {
        move( button, x, y );
    }

    shared_ptr< mouser_i > scale_mouser::clone( double x, double y ) const
    {
        return shared_ptr< mouser_i >( new scale_mouser( x, y ) );
    }

    translate_mouser::translate_mouser()
        :mouser_i( "translate" )
    {
    }
    
    translate_mouser::translate_mouser( double startx, double starty )
    {
        m_startx = startx;
        m_starty = starty;
        glGetDoublev( GL_MODELVIEW_MATRIX, m_modelview );
    }

    translate_mouser::~translate_mouser()
    {
    
    }

    void translate_mouser::move( int button, double x, double y )
    {
        const double transform[] = {
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            x - m_startx,
            y - m_starty, 
            0.0, 1.0 
        };

        glLoadMatrixd( transform );
        glMultMatrixd( m_modelview );
   

//    double result[16];
    
//    result[0] = t[0] * m[0] + t[4] * m[1] + t[8] * m[2] + m12 * m[3];
//    result[1] = t[1] * m[

    


//    glTranslated( m_scale * ( x - m_startx ), 
//                  m_scale * ( m_starty - y ), 
//                  0.0 );

    }

    void translate_mouser::release( int button, double x, double y )
    {
        move( button, x, y );
    }

    mouser_ptr translate_mouser::clone( double x, double y ) const
    {
        return mouser_ptr( new translate_mouser( x, y ) );
    }
    

    rotate_mouser::rotate_mouser()
        :mouser_i( "rotate" )
    {
    }

    rotate_mouser::rotate_mouser( double startx, double starty )
    {
        m_startx = startx;
        m_starty = starty;
        glGetDoublev( GL_MODELVIEW_MATRIX, m_modelview );
    }

    rotate_mouser::~rotate_mouser()
    {
    }

    void rotate_mouser::move( int button, double x, double y )
    {
        double dx = ( x - m_startx ) / 10.0;
        double dy = ( m_starty - y ) / 10.0;


        const double transform[] = {
            cos(dx), sin(dx)*sin(dy), -sin(dx) * cos(dy), 0.0,
            0.0,     cos(dy),         sin(dy),            0.0,
            sin(dx), -cos(dx)*sin(dy), cos(dx) * cos(dy), 0.0,
            0.0,         0.0,             0.0,            1.0 
        };

        glMatrixMode( GL_MODELVIEW );
    
        glLoadMatrixd( transform );
        glMultMatrixd( m_modelview );
    
    }

    void rotate_mouser::release( int button, double x, double y )
    {
        move( button, x, y );
    }

    mouser_ptr rotate_mouser::clone( double x, double y ) const
    {
        return mouser_ptr( new rotate_mouser( x, y ) );
    }

    select_mouser::select_mouser( )
        : mouser_i( "select" )
    {
    }

    select_mouser::select_mouser( double startx, double starty )
    {
        m_startx = startx;
        m_starty = starty;

        m_pselected = content().get_avec( "selected" );
        if( m_pselected )
        {
            for( int i=0; i < m_pselected->size(); ++i )
            {
                m_pselected->at(i).set_i(SELECTED, 0);
            }
            m_pselected->clear();
        }        
        else
        {
            m_pselected = atomvec_ptr( new atomvec_t() );
            content().set( "selected", m_pselected );
        }

        double transmat[16];
        glGetDoublev( GL_MODELVIEW_MATRIX, transmat );
        if( !is_identity(transmat, 4) )
        {
            transform( content(), transmat );
            glMatrixMode( GL_MODELVIEW );
            glLoadIdentity( );
        }
    }

    select_mouser::~select_mouser()
    {
    }

    void select_mouser::move( int button, double x, double y )
    {
        assert( m_pselected != NULL );

        double transmat[16];
        glGetDoublev( GL_MODELVIEW_MATRIX, transmat );
        assert( is_identity(transmat, 4) );

        for( int i=0; i < m_pselected->size(); ++i )
        {
            m_pselected->at(i).set_i(SELECTED, 0);
        }

        m_pselected->clear();

        GL::set_material( int(SELECTED) );
        glBegin( GL_LINE_LOOP );
        glVertex3d( m_startx, m_starty, 0.0 );
        glVertex3d( m_startx, y, 0.0 );
        glVertex3d( x, y, 0.0 );
        glVertex3d( x, m_starty, 0.0 );
        glVertex3d( m_startx, m_starty, 0.0 );
        glEnd();

        assert( m_pselected->size()==0 );
        locate_atom( content(), numvec(m_startx, m_starty, x, y), *m_pselected, FIND_ALL );

        for( int i=0; i < m_pselected->size(); ++i )
        {
            m_pselected->at(i).set_i( SELECTED, 1 );
        }

        drawing()->update();
    }

    void select_mouser::release( int button, double x, double y )
    {
    }

    mouser_ptr select_mouser::clone( double x, double y ) const
    {
        return mouser_ptr( new select_mouser( x, y ) );
    }


    void incbo( moref_t& b )
    {
        int order = b.get_i( ORDER );
        b.set_i( ORDER, 1+order%3 );
    }
    
    void unselect( molecule_t& m )
    {
        int natom = m.natom();
        for( int i=0; i < natom; ++i )
        {
            m.atoms()[i].set_i( SELECTED, 0 );
        }
    }

    drawbond_mouser::drawbond_mouser()
        : mouser_i( "drawbond" )
    {
    }

    drawbond_mouser::drawbond_mouser( double startx, double starty )
    {
        double transmat[16];
        glGetDoublev( GL_MODELVIEW_MATRIX, transmat );
        if( !is_identity(transmat, 4) )
        {
            transform( content(), transmat );
            glMatrixMode( GL_MODELVIEW );
            glLoadIdentity( );
        }
  
        numvec rect( startx-0.3, starty-0.3, startx+0.3, starty+0.3);
        if( locate_bond(content(), rect, m_bond, FIND_ONE) )
        {
            assert( m_bond.size()==1 );
            incbo( m_bond[0] );
            return;
        }

        assert( m_atm1.size()==0 );
        locate_atom( content(), rect, m_atm1, FIND_ONE );
        if( m_atm1.size()>0 && m_atm1[0].getmol().get_s(NAME) != "skeleton" )
        {
            init_3d( startx, starty );
        }
        else
        {
            init_2d( startx, starty );
            drawing()->update( "molecule_fragment" );
        }

        drawing()->update( "molecule_" + m_pmol->get_s(NAME) );
        
    }

    drawbond_mouser::~drawbond_mouser()
    {
    }
    
    void drawbond_mouser::move( int button, double x, double y )
    {
        if( m_bond.size()>0 )
        {
            drawing()->update( "molecule_" + m_bond[0].getmol().get_s(NAME) );
            return;
        }

        double transmat[16];
        glGetDoublev( GL_MODELVIEW_MATRIX, transmat );
        assert( is_identity(transmat, 4) );

        static const double BOND_LENGTH = 1.5;
        static const double THTSTEP = M_PI/6.0;

        if( m_atm2.size()>0 )
        {
            assert( m_atm2[0] != m_atm1[0] );
            m_atm2[0].set_i( SELECTED, 0 );
            m_atm2.clear();
        }

        numvec rect1( x-0.3, y-0.3, x+0.3, y+0.3);
        assert( m_atm2.size()==0 );
        if( locate_atom(content(), rect1, m_atm2, FIND_ONE) && m_atm1[0]==m_atm2[0] )
        {
            m_atm2.clear();
        }

        if( m_atm2.size()==0 || m_atm2[0]==m_atm1[0] )
        {
            double dx = x - m_startx;
            double dy = y - m_starty;
            double th = dx>0.0 ? atan(dy/dx) : M_PI+atan(dy/dx);
            th = int(th/THTSTEP+1)*THTSTEP;

            double nx = m_startx + BOND_LENGTH*cos(th);
            double ny = m_starty + BOND_LENGTH*sin(th);
            m_frag.atoms()[1].set_v(POSITION, numvec(nx, ny,0.0) );
        
            numvec rect2( nx-0.3, ny-0.3, nx+0.3, ny+0.3);
            assert( m_atm2.size()==0 );
            if( locate_atom(content(), rect2, m_atm2, FIND_ONE) && m_atm2[0]==m_atm1[0] )
            {
                m_atm2.clear();
            }
        }
        
        if( m_atm2.size()>0 )
        {
            assert( m_atm2.size()==1 );
            m_atm2[0].set_i( SELECTED, 1 );
            numvec p = m_atm2[0].get_v(POSITION);
            if(m_dim==2) m_frag.atoms()[1].set_v(POSITION, p);
        }
        
        drawing()->update( "molecule_" + m_pmol->get_s(NAME) );
        drawing()->update( "molecule_fragment" );
    }
    
    void drawbond_mouser::release( int button, double x, double y )
    {
        if( m_bond.size()>0 )
        {
            drawing()->update( "molecule_" + m_bond[0].getmol().get_s(NAME) );
            return;
        }

        if( m_atm1.size() >0 )
        {
            m_atm1[0].set_i( SELECTED, 0 );
        }

        if( m_atm2.size() >0 )
        {
            m_atm2[0].set_i( SELECTED, 0 );
        }
        
        unselect( m_frag );

        if( m_atm1.size()>0 && m_atm2.size()>0 )
        {
            if( &m_atm1[0].getmol() != &m_atm2[0].getmol() )
            {
                throw logic_error("Error: can not create bond between two molelcules, merge them together first!");
            }
            
            moref_t b( m_atm1[0].getmol(), BOND, -1 );
            if( get_bond(m_atm1[0], m_atm2[0], b) )
            {
                incbo( b );
            }
            else
            {
                b = create_bond( m_atm1[0], m_atm2[0] );
                b.set_i( ORDER, 1 );
            }

        }
        else
        {
            assert( m_atm2.size()==0 );
        
            if( m_dim==2 )
            {
                merge_bypos( *m_pmol, m_frag );
            }
            else
            {
                assert( m_dim==3 );

                if( m_atm1[0].get_i(ELEMENT) != CARBON )
                {
                    m_atm1[0].set_i(ELEMENT, CARBON);
                    m_atm1[0].set_s(NAME,    uniq_name(m_atm1[0]) );
                    m_atm1[0].set_s(TYPE,    "C");
                }
        
                addHs( m_atm1[0] );
            }
        }
        
        drawing()->update( "molecule_" + m_pmol->get_s(NAME) );
        drawing()->remove( "molecule_fragment" );
    }
    
    mouser_ptr drawbond_mouser::clone( double x, double y ) const
    {
        return mouser_ptr( new drawbond_mouser(x, y) );
    }

    void drawbond_mouser::init_frag( molecule_t& frag )
    {
        moref_t atm1 = frag.create_atom();
        moref_t atm2 = frag.create_atom();
        moref_t bond = create_bond( atm1, atm2 );

        atm1.set_v( POSITION, numvec(0.0, 0.0, 0.0) );
        atm2.set_v( POSITION, numvec(1.5, 0.0, 0.0) );
        
        atm1.set_i( ELEMENT, CARBON );
        atm2.set_i( ELEMENT, CARBON );

        atm1.set_s( NAME, "C1" );
        atm2.set_s( NAME, "C2" );
        
        atm1.set_s( TYPE, "C" );
        atm2.set_s( TYPE, "C" );

        atm1.set_i( SELECTED, 1 );
        atm2.set_i( SELECTED, 1 );
        
        bond.set_i( ORDER, 1 );
        
        m_frag.set_s( NAME, "fragment" );
    }

    void drawbond_mouser::init_2d( double startx, double starty )
    {
        double transmat[16];
        glGetDoublev( GL_MODELVIEW_MATRIX, transmat );
        assert( is_identity(transmat, 4) );
        
        m_dim = 2;
        m_pmol = content().get_mol( "skeleton" ).get();
        if( m_pmol == NULL )
        {
            molecule_ptr pskt = molecule_ptr( new molecule_t() );
            pskt->set_s( NAME, "skeleton" );
            content().set( "skeleton", pskt );
            drawing()->add( "molecule", *pskt, "ballstick" );
            m_pmol = pskt.get();
        }
        
        if( m_atm1.size()==0 )
        {
            m_atm1.push_back( m_pmol->create_atom() );
            m_atm1[0].set_v( POSITION, numvec(startx, starty, 0.0) );
            m_atm1[0].set_i( ELEMENT,  CARBON );
            m_atm1[0].set_s( NAME,     uniq_name(m_atm1[0]) );
            m_atm1[0].set_s( TYPE,     "C" );
        }
        else
        {
            assert( &m_atm1[0].getmol() == m_pmol );
            numvec p = m_atm1[0].get_v(POSITION);
            startx = p[0];
            starty = p[1];
        }
        
        m_atm1[0].set_i(SELECTED, 1);
        
    
        init_frag( m_frag );
        translate( m_frag, startx, starty );
        
        drawing()->add( "molecule", m_frag, "line" );
        m_startx = startx;
        m_starty = starty;
        m_curth  = 0.0;
    }

    void drawbond_mouser::init_3d( double startx, double starty )
    {
        assert( m_atm1.size()>0 );
        m_dim = 3;
        m_pmol = &(m_atm1[0].getmol());

        init_frag( m_frag );
        translate( m_frag, startx, starty );
        
        drawing()->add( "molecule", m_frag, "line" );
        m_startx = startx;
        m_starty = starty;
        m_curth  = 0.0;
    }
    

} // namespace mort

mort::scale_mouser     g_scale_mouser;

mort::rotate_mouser    g_rotate_mouser;

mort::select_mouser    g_select_mouser;

mort::drawbond_mouser  g_drawbond_mouser;

mort::translate_mouser g_translate_mouser;

