#include <sstream>
#include <stdexcept>
#include <object.hpp>
#include <atmask.hpp>
#include "mainwin.hpp"
#include "shidedlg.hpp"

namespace mortgtk
{
    using std::runtime_error;

    using namespace mort;

    class show_command : public command_i
    {
    public:

        show_command()
	{
	}

        show_command( const string& mole, const string& oper, const string& graf, const string& mask )
	    : m_mole(mole), m_oper(oper), m_graf(graf), m_mask(mask) 
	{
	}

        virtual bool exec()
	{
            string grafname = m_graf + "_" + m_mole;
            graphic_ptr graf = drawing()->get( grafname );
            if( graf==NULL )
            {
                return true;
            }

	    assert( m_oper=="on" || m_oper=="off" || m_oper=="only" );

	    if( m_oper=="only" )
	    {
	        graf->setall( VISIBLE, OFF );
	    }

            molecule_ptr pmol = content().get_mol( m_mole );

	    atomvec_t atms = unified_mask_atom( *pmol, m_mask );
            // std::cout << atms.size() << " selected" << std::endl;

	    if( m_oper=="on" || m_oper=="only" )
	    {
	        graf->set( atms, VISIBLE, ON );
	    }
	    else
	    {
	        graf->set( atms, VISIBLE, OFF );
	    }

            drawing()->update( grafname );
	    return true;
	}

        virtual void undo()
	{
	}

	virtual command_ptr clone(const vector<string>& args) const
	{
	    return command_ptr( new show_command() );
	}

        string m_mole;
        string m_graf;
	string m_mask;
	string m_oper;
    };

    shidedlg_t::shidedlg_t()
        : basicdlg_t( "Show/hide part of molecule..." ),
          m_operate( "Select operator:" ),
          m_graphic( "Select graphic:" ),
          m_atomask( "Atom selection:" )
    {

        m_operate.init( "on", "off", "only" );
        m_graphic.init( "molecule", "label" );
        m_atomask.init( );
        
        gtk_box_pack_start( GTK_BOX(vbox()), (GtkWidget*)m_operate, true, true, 0);
	gtk_box_pack_start( GTK_BOX(vbox()), (GtkWidget*)m_graphic, true, true, 0);
	gtk_box_pack_start( GTK_BOX(vbox()), (GtkWidget*)m_atomask, true, true, 0);


	GtkWidget* nulllabel = gtk_label_new( " " );
	gtk_box_pack_start( GTK_BOX(vbox()), nulllabel, true, true, 0 );
	gtk_widget_show( nulllabel );

        dlgfactory_t::add( "shide", this );
    }

    shidedlg_t::~shidedlg_t()
    {
    }

    command_ptr shidedlg_t::get_command()
    {
        string atomask = m_atomask.getmsk();
        string operate = m_operate.getstr();
        string graphic = m_graphic.getstr();

        assert( !atomask.empty() && !molname().empty() );
        
        molecule_ptr pmol =  content().get_mol(molname());
        
        assert( pmol );

        return command_ptr( new show_command(molname(), operate, graphic, atomask) );
    }

} // namespace mortgtk

