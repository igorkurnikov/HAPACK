#include <stdexcept>
#include "mainwin.hpp"
#include "basicdlg.hpp"


namespace mortgtk
{


    basicdlg_t::basicdlg_t( const char* title)
        : m_molname( "Select molecule:" )
    {
        m_self = (GtkDialog*)gtk_dialog_new_with_buttons( title, NULL,
                                                          GTK_DIALOG_MODAL,
                                                          GTK_STOCK_CANCEL, GTK_RESPONSE_REJECT,
                                                          GTK_STOCK_OK,     GTK_RESPONSE_OK,
                                                          NULL );

        m_molname.init( content() );

        m_molname.listen( this );

        gtk_box_pack_start( GTK_BOX(m_self->vbox), (GtkWidget*)m_molname, true, true, 0 );
    }
    
    basicdlg_t::~basicdlg_t()
    {
        gtk_widget_destroy( GTK_WIDGET(m_self) );
    }

    string basicdlg_t::molname() const
    {
        return m_molname.getstr();
    }

    dlgfactory_t::dlgfactory_t()
    {
    }
    
    dlgfactory_t::~dlgfactory_t()
    {
    }

    dlgfactory_t& dlgfactory_t::inst()
    {
        static dlgfactory_t in;
        return in;
    }
    
    void dlgfactory_t::add( const string& name, basicdlg_t* impl )
    {
        if( inst().m_data.count(name)==0 )
            inst().m_data[name] = impl;
    }
    
    void dlgfactory_t::run( const string& name )
    {
        if( inst().m_data.count(name)==0 )
        {
            popup_message( "Error: can not find dialog with type " + name );
            return;
        }
        
        
        basicdlg_ptr pdlg = inst().m_data[name]->clone();
    
        try
        {
            if( pdlg->run()==GTK_RESPONSE_OK )
            {
                command_ptr pcmd = pdlg->get_command();
                assert( pcmd != NULL );
                pcmd->exec();
            }

        }
        catch( std::exception& e )
        {
            popup_message( e.what() );
        }

        drawing()->update();

    }
    

    
} // namespace mortgtk

