#include <sstream>
#include <stdexcept>
#include <object.hpp>
#include <atmask.hpp>
#include "mainwin.hpp"
#include "showdlg.hpp"

namespace mortgtk
{
    using std::runtime_error;

    using namespace mort;

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

    class show_command : public command_i
    {
    public:

        show_command()
	{
	}

        show_command( drawing_t* area, molecule_t* pmol, const string& oper, const string& mask )
	    : m_mask( mask ),
	      m_oper( oper )
	{
	    m_area = area;
	    m_pmol = pmol;
	    m_graf = area->get( "molecule_" + pmol->get_s(NAME) );
	}

        virtual bool exec()
	{
	    assert( m_oper=="on" || m_oper=="off" || m_oper=="only" );

	    if( m_oper=="only" )
	    {
	        m_graf->setall( VISIBLE, OFF );
	    }

	    atomvec_t atms = unified_mask_atom( *m_pmol, m_mask );
            // std::cout << atms.size() << " selected" << std::endl;

	    if( m_oper=="on" || m_oper=="only" )
	    {
	        m_graf->set( atms, VISIBLE, ON );
	    }
	    else
	    {
	        m_graf->set( atms, VISIBLE, OFF );
	    }

            m_area->update( "molecule_" + m_pmol->get_s(NAME) );
	    return true;
	}

        virtual void undo()
	{
	}

	virtual command_ptr clone(const vector<string>& args) const
	{
	    return command_ptr( new show_command() );
	}

        drawing_t* m_area;
	molecule_t* m_pmol;
	graphic_ptr m_graf;

	string m_mask;
	string m_oper;
    };

    bool ismol( entity_ptr& pe )
    {
        molecule_ptr pm = dynamic_pointer_cast< molecule_t >( pe );
	return pm != NULL;
    }

    void on_operate_change(GtkWidget* toggled, gpointer user_data )
    {
        showdlg_t* dlg = (showdlg_t*)user_data;
	dlg->on_operate_change(toggled);
    }

    void on_atomask_change(GtkWidget* sel, gpointer user_data )
    {
        showdlg_t* dlg = (showdlg_t*)user_data;
	dlg->on_atomask_change(sel);
    }

    void on_molname_change(GtkWidget* sel, gpointer user_data )
    {
        showdlg_t* dlg = (showdlg_t*)user_data;
	dlg->on_molname_change(sel);
    }

    showdlg_t::showdlg_t(main_window& mainwin)
    {
        m_self = (GtkDialog*)gtk_dialog_new_with_buttons( "Show/hide part of molecule...", NULL,
                                                         GTK_DIALOG_MODAL,
	 	                                         GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                                         GTK_STOCK_OK,     GTK_RESPONSE_OK,
					                 NULL );

        m_area = drawing().get();
	m_data = &content();
        m_operate = "on";

        GtkWidget* molelabel = gtk_label_new( "\nWhich molecule to show/hide?" );
	gtk_box_pack_start( GTK_BOX(m_self->vbox), molelabel, true, true, 0);
        gtk_widget_show( molelabel );

        GtkWidget* molebutton = gtk_combo_box_new_text();
	gtk_box_pack_start( GTK_BOX(m_self->vbox), molebutton, true, true, 0);
	gtk_widget_show( molebutton );
        database_t::iterator i = m_data->begin();
	for( ; i != m_data->end(); ++i )
	{
	    string name = i->first;
	    if( name[0] != '_' && ismol(i->second) )
	        gtk_combo_box_append_text( GTK_COMBO_BOX(molebutton), name.c_str() );
	}
	g_signal_connect( G_OBJECT(molebutton), "changed", G_CALLBACK(::on_molname_change), this );

        GtkWidget* separator1 = gtk_hseparator_new();
	gtk_box_pack_start( GTK_BOX(m_self->vbox), separator1, true, true, 0);
	gtk_widget_show( separator1 );

        GtkWidget* operlabel = gtk_label_new( "\nWhat operation to perform?" );
	gtk_box_pack_start( GTK_BOX(m_self->vbox), operlabel, true, true, 0);
        gtk_widget_show( operlabel );

        m_oper_on   = gtk_radio_button_new_with_label( NULL, "on" );
        m_oper_off  = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(m_oper_on), "off" );
        m_oper_only = gtk_radio_button_new_with_label_from_widget( GTK_RADIO_BUTTON(m_oper_on), "only" );

	gtk_box_pack_start( GTK_BOX(m_self->vbox), m_oper_on,   true, true, 0);
	gtk_box_pack_start( GTK_BOX(m_self->vbox), m_oper_off,  true, true, 0);
	gtk_box_pack_start( GTK_BOX(m_self->vbox), m_oper_only, true, true, 0);

	g_signal_connect( G_OBJECT(m_oper_on),   "toggled", G_CALLBACK(::on_operate_change), this );
	g_signal_connect( G_OBJECT(m_oper_off),  "toggled", G_CALLBACK(::on_operate_change), this );
	g_signal_connect( G_OBJECT(m_oper_only), "toggled", G_CALLBACK(::on_operate_change), this );

	gtk_widget_show( m_oper_on );
	gtk_widget_show( m_oper_off );
	gtk_widget_show( m_oper_only );

        GtkWidget* separator2 = gtk_hseparator_new();
	gtk_box_pack_start( GTK_BOX(m_self->vbox), separator2, true, true, 0);
	gtk_widget_show( separator2 );

        GtkWidget* selabel = gtk_label_new( "\nWhich part of molecule?" );
	gtk_box_pack_start( GTK_BOX(m_self->vbox), selabel, true, true, 0);
        gtk_widget_show( selabel );

        GtkWidget* selwidget = gtk_combo_box_new_text();
	gtk_combo_box_append_text( GTK_COMBO_BOX(selwidget), "all" );
	gtk_combo_box_append_text( GTK_COMBO_BOX(selwidget), "hydrogens" );
	gtk_combo_box_append_text( GTK_COMBO_BOX(selwidget), "backbone" );
        gtk_combo_box_append_text( GTK_COMBO_BOX(selwidget), "ambmask..." );
	gtk_combo_box_append_text( GTK_COMBO_BOX(selwidget), "smarts...." );
	g_signal_connect( G_OBJECT(selwidget), "changed", G_CALLBACK(::on_atomask_change), this );
	gtk_box_pack_start( GTK_BOX(m_self->vbox), selwidget, true, true, 0 );
	gtk_widget_show( selwidget );

        m_maskentry = gtk_entry_new();
	gtk_box_pack_start( GTK_BOX(m_self->vbox), m_maskentry, true, true, 0 );

	GtkWidget* nulllabel = gtk_label_new( " " );
	gtk_box_pack_start( GTK_BOX(m_self->vbox), nulllabel, true, true, 0 );
	gtk_widget_show( nulllabel );

    }

    showdlg_t::~showdlg_t()
    {
        gtk_widget_destroy( GTK_WIDGET(m_self) );
    }

    void showdlg_t::on_operate_change(GtkWidget* toggled)
    {
        if( toggled==m_oper_on )
	{
	    m_operate = "on";
	}
	else if( toggled==m_oper_off )
	{
	    m_operate = "off";
	}
	else
	{
	    assert( toggled==m_oper_only );
	    m_operate = "only";
	}

    }

    void showdlg_t::on_atomask_change(GtkWidget* widget)
    {
	GtkTreeIter iter;
        GtkTreeModel* model = gtk_combo_box_get_model( GTK_COMBO_BOX(widget) );

        if( !gtk_combo_box_get_active_iter(GTK_COMBO_BOX(widget), &iter) )
	    return;

	gchar* seltype;
	gtk_tree_model_get(model, &iter, 0, &seltype, -1);
	m_atomask = seltype;
	g_free(seltype);

	if( m_atomask=="ambmask..." || m_atomask=="smarts...." )
	{
	    gtk_widget_show( m_maskentry );
	}
	else
	{
	    gtk_widget_hide( m_maskentry );
	}

    }

    void showdlg_t::on_molname_change(GtkWidget* widget)
    {
        GtkTreeIter iter;
	GtkTreeModel* model = gtk_combo_box_get_model( GTK_COMBO_BOX(widget) );

        if( !gtk_combo_box_get_active_iter(GTK_COMBO_BOX(widget), &iter) )
	    return;

	gchar* molname;
	gtk_tree_model_get(model, &iter, 0, &molname, -1);
	m_molname = molname;
	g_free(molname);
    }

    command_ptr showdlg_t::get_command()
    {
	assert( !m_operate.empty() );

        if( m_atomask.empty() )
	{
	    throw runtime_error( "Error: no atom selected" );
	}

	if( m_molname.empty() )
	{
	    throw runtime_error( "Error: no molecule selected" );
	}

	if( m_atomask=="ambmask..." || m_atomask=="smarts...." )
	{
            string mask = gtk_entry_get_text( GTK_ENTRY(m_maskentry) );
            if( mask.empty() )
	    {
	        throw runtime_error( "Warning: no mask is given" );
	    }

            m_atomask.append( 1, ' ' );
	    m_atomask += mask;
	}

        molecule_ptr pmol =  m_data->get_mol(m_molname);
        if( pmol == NULL )
	{
	    throw runtime_error( "Error: no such molecule " + m_molname );
	}

	return command_ptr( new show_command(m_area, pmol.get(), m_operate, m_atomask) );
    }

    int showdlg_t::run( )
    {
        return gtk_dialog_run( m_self );
    }

} // namespace mortgtk

