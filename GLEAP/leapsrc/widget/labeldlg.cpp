#include <sstream>
#include <stdexcept>
#include <object.hpp>
#include <atmask.hpp>
#include "mainwin.hpp"
#include "labeldlg.hpp"

namespace mortgtk
{
    using std::runtime_error;

    using namespace mort;

    class label_command : public command_i
    {
    public:

        label_command()
	{
	}

        label_command( const string& mname, const string& level, const string& param )
	    : m_mname(mname), m_level(level), m_param(param) 
	{
	}

        virtual bool exec()
	{
            molecule_ptr pmol = content().get_mol( m_mname );
            
            assert( pmol != NULL );
            
            drawing()->add( "label", *pmol, m_level + " " + m_param );
	}

        virtual void undo()
	{
	}

	virtual command_ptr clone(const vector<string>& args) const
	{
            assert( args.size()==4 );
            
	    return command_ptr( new label_command(args[1], args[2], args[3]) );
	}

        string m_mname;
	string m_level;
	string m_param;
    };

    labeldlg_t::labeldlg_t()
        : basicdlg_t( "Show/hide part of molecule..." ),
          m_level( "Select level:" ),
          m_param( "Select param:" )
    {
        m_level.init( "atom", "resd", "molecule" );
        m_param.init( "name", "type", "pchg" );

	gtk_box_pack_start( GTK_BOX(vbox()), (GtkWidget*)m_level, true, true, 0);
	gtk_box_pack_start( GTK_BOX(vbox()), (GtkWidget*)m_param, true, true, 0);

	GtkWidget* nulllabel = gtk_label_new( " " );
	gtk_box_pack_start( GTK_BOX(vbox()), nulllabel, true, true, 0 );
	gtk_widget_show( nulllabel );

        dlgfactory_t::add( "label", this );        
    }

    labeldlg_t::~labeldlg_t()
    {
    }

    command_ptr labeldlg_t::get_command()
    {
	if( molname().empty() )
	{
	    throw runtime_error( "Error: no molecule selected" );
	}

        string param = m_param.getstr();
        string level = m_level.getstr();

        if( param.empty() )
	{
	    throw runtime_error( "Error: no param selected" );
	}

        assert( !level.empty() );

	return command_ptr( new label_command( molname(), level, param) );
    }

} // namespace mortgtk

