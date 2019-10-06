#include <stdexcept>
#include <object.hpp>
#include <GL/glx.h>

#include "ribbon.hpp"
#include "graphic.hpp"
#include "glstyle.hpp"

namespace mort
{
    using std::logic_error;

    graphic_t::graphic_t(const molecule_t& m)
        : m_pmol(&m)
    {
    }

    void graphic_t::setall(const hashid_t& parmid, const hashid_t& value)
    {
        moiter_t ai = m_pmol->atom_begin();
	for( ; ai != m_pmol->atom_end(); ++ai )
	{
	    int absid = ai->absid();
	    m_params[parmid][absid] = value;
	}
    }

    void graphic_t::set(const atomvec_t& bound, const hashid_t& parmid, const hashid_t& value)
    {
        for( int i=0; i < bound.size(); ++i )
	{
	    int absid = bound[i].absid();
	    m_params[parmid][absid] = value;
	}
    }

    bool graphic_t::get(int absid, const hashid_t& parmid, hashid_t& v) const
    {
        map< hashid_t, map<int, hashid_t> >::const_iterator i = m_params.find(parmid);
	if( i == m_params.end() )
	{
	    return false;
	}

	map<int, hashid_t>::const_iterator j = i->second.find(absid);
	if( j == i->second.end() )
	{
	    return false;
	}

	v = j->second;
	return true;
    }

    shared_ptr<graphic_t> graphic_t::construct(const molecule_t& m, const string& name, const string& argument)
    {
        vector<string> args;
        interpret(argument, args);

        map<string, const graphic_t*>::iterator i = factory().find(name);
	if( i == factory().end() )
	{
	    throw logic_error("Error: undefined graphic type: " + name);
        }

	return i->second->clone(m, args);
    }

    map<string, const graphic_t*>& graphic_t::factory()
    {
        static molecule_graphic mg;
	static label_graphic    lg;
	static ribbon_graphic   rg;
        static map<string, const graphic_t*> g_factory;

	if( g_factory.size()==0 )
	{
            g_factory["molecule"] = &mg;
            g_factory["label"] = &lg;
	    g_factory["ribbon"] = &rg;
	}

        return g_factory;
    }

    molecule_graphic::molecule_graphic()
    {
        m_pmol = NULL;
    }

    molecule_graphic::molecule_graphic(const molecule_t& m, const hashid_t& style)
 	:graphic_t(m), m_pmol(&m), m_name( "molecule_" + m.get_s(NAME) )
    {
        if( style == LINE )
        {
            m_style = glstyle_ptr( new line_style(1.5) );
        }
        else if( style == BALLSTICK )
        {
            m_style = glstyle_ptr( new ballstick_style(0.4, 0.2, ELEMENT_COLOR) );
        }
        else if( style == STICK )
        {
            m_style = glstyle_ptr( new ballstick_style(0.2, 0.2, ELEMENT_COLOR) );
        }
        else
        {
            throw logic_error( "Error: unknown molecule render style." );
        }
    }

    molecule_graphic::~molecule_graphic()
    {
    }

    string molecule_graphic::name() const
    {
        return m_name;
    }

    void molecule_graphic::paint()
    {
        m_style->prepare();
        
        moiter_t ai = m_pmol->atom_begin();
        for( ; ai != m_pmol->atom_end(); ++ai )
        {
	    hashid_t visible;

            if( get( ai->absid(), VISIBLE, visible) && visible==OFF )
            {
	        continue;
	    }

            m_style->render_atom(*ai);
        }

        moiter_t bi = m_pmol->bond_begin();
        for( ; bi != m_pmol->bond_end(); ++bi )
        {
	
	    int absid_1 = atom_1st(*bi).absid();
	    int absid_2 = atom_2nd(*bi).absid();

	    hashid_t visible1;
	    hashid_t visible2;

            bool exist1 = get(absid_1, VISIBLE, visible1);
	    bool exist2 = get(absid_2, VISIBLE, visible2);

            if(exist1 && visible1==OFF)
	    {
	        continue;
	    }

	    if(exist2 && visible2==OFF )
	    {
	        continue;
	    }
 
            m_style->render_bond(*bi); 
        }
 

    }

    void molecule_graphic::update()
    {
    }

    shared_ptr<graphic_t> molecule_graphic::clone(const molecule_t& m, const vector<string>& args) const
    {
        assert( args.size()==1 );
	return shared_ptr<graphic_t>( new molecule_graphic(m, mort::hash(args[0]) ) );
    }

    label_graphic::label_graphic()
    {
    }

    label_graphic::label_graphic(const molecule_t& m, const hashid_t& level, const hashid_t& parm)
        : graphic_t(m), m_pmol(&m), m_level(level), m_parm(parm)
    {
        m_base = 1000;
        m_name = "label_" + m_pmol->get_s(NAME);
        Display* dpy = glXGetCurrentDisplay();
        string fontname = "-*-times-medium-i-normal-*-14-100-*-*-*-*-*-*";
	XFontStruct* xfont = XLoadQueryFont(dpy, fontname.c_str() );
	assert(xfont != NULL);
	glXUseXFont( xfont->fid, 0, 128, m_base);
    }

    label_graphic::~label_graphic()
    {
    }

    string label_graphic::name() const
    {
        return m_name;
    }

    void label_graphic::paint()
    {
        GL::set_material(int(SELECTED));

        glListBase(m_base);
	moiter_t ai = m_pmol->atom_begin();
        for( ; ai != m_pmol->atom_end(); ++ai)
	{
	    int absid = ai->absid();
            hashid_t visible;
            if( get(absid, VISIBLE, visible) && visible==OFF )
	    {
	        continue;
	    }

	    numvec pos = ai->get_v(POSITION);
	    string value = parm2str( *ai, m_parm);
	    glRasterPos3d( pos[0], pos[1], pos[2] );
            glCallLists( value.length(), GL_BYTE, value.c_str() );
	}
    }

    void label_graphic::update()
    {
    }

    shared_ptr<graphic_t> label_graphic::clone(const molecule_t& m, const vector<string>& args) const
    {
        assert( args.size()==2 );
	return shared_ptr<graphic_t>( new label_graphic(m, mort::hash(args[0]), mort::hash(args[1])) );
    }

    ribbon_graphic::ribbon_graphic()
    {
    }

    ribbon_graphic::ribbon_graphic(const molecule_t& m)
        : graphic_t(m), m_name( "ribbon_" + m.get_s(NAME) )
    {
        m_pmol = &m;
    }

    ribbon_graphic::~ribbon_graphic()
    {
    }

    string ribbon_graphic::name( ) const
    {
         return m_name;
    }

    void ribbon_graphic::paint()
    {
        assert( m_pmol != NULL );
        render_ribbon(*m_pmol);
    }

    void ribbon_graphic::update()
    {
    }

    shared_ptr<graphic_t> ribbon_graphic::clone(const molecule_t& m, const vector<string>& args) const
    {
        return shared_ptr<graphic_t>( new ribbon_graphic(m) );
    }


} // namespace mort


