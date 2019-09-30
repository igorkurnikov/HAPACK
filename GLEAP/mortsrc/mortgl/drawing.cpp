#ifdef WIN32
#include <windef.h>
#include <wingdi.h>
#endif 

#include <GL/gl.h>
#include <GL/glu.h>
#include "ortho.hpp"
#include "mouse.hpp"
#include "drawing.hpp"
#include "graphic.hpp"

namespace mort
{

    drawing_t::drawing_t()
        :mp_ortho( new ortho_t(1.0/32.0) )
    {
        m_quantity = QUANTITY_HIGH;
        m_mouser = "select";
    }

    drawing_t::~drawing_t()
    {
    }

    void drawing_t::add(const string& type, molecule_t& m, const string& argument)
    {
        addone( graphic_t::construct(m, type, argument) );
    }



    void drawing_t::addone(const shared_ptr<graphic_t>& graphic)
    {
        string gname = graphic->name();
        
        for(int i=0; i < m_graphics.size(); ++i )
        {
            if( m_graphics[i]->name() == gname )
            {
                m_graphics[i] = graphic;
                return;
            }
        }
        
        m_graphics.push_back( graphic );
    }

    bool drawing_t::has( const string& gname ) const
    {
        vector< shared_ptr<graphic_t> >::const_iterator i = m_graphics.begin();

        for( ; i != m_graphics.end(); ++i )
        {
            if( (*i)->name() == gname )
            {
                return true;
            }
        }

        return false;
    }

    shared_ptr<graphic_t> drawing_t::get(const string& gname) const
    {
        vector< shared_ptr<graphic_t> >::const_iterator i = m_graphics.begin();
	for( ; i != m_graphics.end(); ++i )
	{
	    if( (*i)->name() == gname )
	    {
	        return *i;
	    }
	}

	return shared_ptr<graphic_t>();
    }

    void drawing_t::remove( const string& gname )
    {
        vector< shared_ptr<graphic_t> >::iterator i = m_graphics.begin();

        for( ; i != m_graphics.end(); ++i )
        {
            string curt = (*i)->name();
            if(curt.find(gname) != string::npos)
            {
                m_graphics.erase( i );
                return;
            }
        }
    }

    void drawing_t::init()
    {
        glClearColor( 0.0, 0.0, 0.0, 1.0 );

        glClearDepth( 1.0 );

        glDepthFunc(GL_LEQUAL);
        glEnable(GL_DEPTH_TEST);

        GL::light( );

        if( m_quantity == QUANTITY_HIGH )
        {
            GL::depth_cue( true );
        }
        else
        {
            GL::depth_cue( false );
        }

#ifdef WIN32
        HDC    hdc = wglGetCurrentDC();
        HGLRC  hglrc = wglGetCurrentContext();
        wglUseFontBitmaps( hdc, 0, 255, 1000);
#endif

    }

    void drawing_t::resize( int width, int height )
    {
        mp_ortho->resize( width, height );
    }

    void drawing_t::update( )
    {
        for( int i=0; i < (int)m_graphics.size(); ++i )
        {
	    string name = m_graphics[i]->name();
            update( name );
        }
    }
    

    void drawing_t::update(const string& name)
    {
        graphic_ptr graf = get( name );
	assert( graf != NULL );

	GL::new_list(name);
        graf->paint();
	glEndList();

    }

    void drawing_t::repaint()
    {
        for( int i=0; i < (int)m_graphics.size(); ++i )
        {
	    string name = m_graphics[i]->name();
	    if( !GL::has_list(name) )
	    {
	        GL::new_list(name);
                m_graphics[i]->paint();
		glEndList();
	    }

	    GL::call_list(name);
        }
    }

    void drawing_t::mouse_press( int button, int state, double x, double y )
    {
        mp_mouser.reset();
       
        if( button == 1 )
        {
            mp_mouser = mouser_i::get( m_mouser.c_str(), mp_ortho->get_x(x), mp_ortho->get_y(y) );
        }
        else if(button == 3) 
        {
	}
	else
	{
	    throw logic_error("Error: unknown button");
	}
    }

    void drawing_t::mouse_move( int button, int , double x, double y )
    {
        if( mp_mouser != NULL )
        {
            mp_mouser->move( button, mp_ortho->get_x( x ), mp_ortho->get_y( y ) );
        }
    }

    void drawing_t::mouse_release( int button, int , double x, double y )
    {
        if( mp_mouser != NULL )
        {
            mp_mouser->release( button, mp_ortho->get_x( x ), mp_ortho->get_y( y ) );
        }

	mp_mouser.reset();
    }

    void drawing_t::mouse_settype( const string& mouser )
    {
        m_mouser = mouser;
    }

    shared_ptr<drawing_t>& drawing()
    {
        static shared_ptr<drawing_t> inst;
	return inst;
    }

} // namespace mort
 
