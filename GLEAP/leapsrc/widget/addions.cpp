#include <plugins/addions.hpp>
#include "addions.hpp"

namespace mortgtk
{
    
    addions_t::addions_t()
        : basicdlg_t( "Add ions to molecule..." ),
          m_ioname( "Select an ion:" ),
          m_ionumb( "Number of ion:" )
    {
        database_ptr pmdb = content().get_mdb( "_ions" );
        assert( pmdb != NULL );
        m_ioname.init( *pmdb );

        vector< string > strs( 1, "to neutral" );
        for( int i=0; i < 9; ++i )
        {
            string tmp = "1";
            tmp[0] += i;
            strs.push_back( tmp );
        }
        m_ionumb.init( strs );
        
        gtk_box_pack_start( GTK_BOX(vbox()), (GtkWidget*)m_ioname, true, true, 0 );
        gtk_box_pack_start( GTK_BOX(vbox()), (GtkWidget*)m_ionumb, true, true, 0 );

        dlgfactory_t::add( "addions", this );
    }
    
    addions_t::~addions_t()
    {
    }
    
    command_ptr addions_t::get_command()
    {
        string ioname = m_ioname.getstr();
        string ionumb = m_ionumb.getstr();
        
        int nion = (ionumb=="to neutral") ? 0 : atoi( ionumb.c_str() );
        
        return command_ptr( new amber::addions_command( molname(), ioname, nion ) );
    }
    
}
