#include <stdexcept>
#include <atmask.hpp>
#include "mskcombo.hpp"


namespace mortgtk
{
    map<string, string>& aliases()
    {
        static map<string, string> inst;
	if( inst.size()==0 )
	{
	    inst["all"]       = "smarts.... [*]";
	    inst["hydrogens"] = "smarts.... [#1]";
            inst["backbone" ] = "ambmask... @C,O,N,CA";
	}

	return inst;
    }

    atomvec_t unified_mask_atom(const molecule_t& m, const string& mask )
    {
        string realmask = aliases().count(mask)>0? aliases()[mask] : mask;

        std::istringstream is(realmask);
        string type, value;
	is >> type >> value;
	
	// std::cout << "type : " << type  << std::endl;
	// std::cout << "value: " << value << std::endl;
        assert( type=="smarts...." || type=="ambmask..." );

	return type=="smarts...." ? smarts_mask_atom( m, value )
            : mask_atom( m, value );
    }

    mskcombo_t::mskcombo_t( const string& title )
        : strcombo_t( title )
    {
        m_maskentry = gtk_entry_new();
        m_entry_on = false;
    }
    
    mskcombo_t::~mskcombo_t( )
    {
    }



    void mskcombo_t::init( )
    {

        vector<string> strs;
        strs.push_back( "all" );
	strs.push_back( "hydrogens" );
	strs.push_back( "backbone" );
        strs.push_back( "ambmask..." );
	strs.push_back( "smarts...." );
        init( strs );

	gtk_box_pack_start( GTK_BOX(vbox()), m_maskentry, true, true, 0 );
    }

    void mskcombo_t::onchange( GtkWidget* w )
    {

        strcombo_t::onchange( w );


        string str = getstr();

        
	if( str=="ambmask..." || str=="smarts...." )
	{
	    if(!m_entry_on)
            {
                gtk_widget_show( m_maskentry );
                m_entry_on = true;
            }
            
	}
	else
	{
            if( m_entry_on)
            {
                gtk_widget_hide( m_maskentry );
                m_entry_on = false;
            }
	}

    }    

    string mskcombo_t::getmsk( )
    {
        string str = getstr();

	if( str=="ambmask..." || str =="smarts...." )
	{
            string mask = gtk_entry_get_text( GTK_ENTRY(m_maskentry) );
            if( mask.empty() )
	    {
	        throw std::runtime_error( "Warning: no mask is given" );
	    }

            string result = str + " " + mask;

            return result;
	}

        return str;
    }
    

} // namespace mortgtk
 
