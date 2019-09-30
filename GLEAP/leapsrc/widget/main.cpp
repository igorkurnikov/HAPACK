#include <GL/gl.h>
#include <fstream>
#include <iostream>
#include <gtk/gtk.h>
#include <gtk/gtkgl.h>
#include <format.hpp>
#include <guilib.hpp>
#include "mainwin.hpp"
#include "clsdata.hpp"
#include "fixbond.hpp"
#include "addhydr.hpp"
#include "setpchg.hpp"
#include "parmchk.hpp"
#include "fastbld.hpp"
#include "ribbdlg.hpp"
#include "addions.hpp"
#include "solvate.hpp"
#include "ffparam.hpp"
#include "gconsole.hpp"
#include "shidedlg.hpp"
#include "labeldlg.hpp"

void init_lib( const string& amberhome )
{
    string path = "./:";
    path += amberhome + string( "/dat/leap/gleap:" );
    path += amberhome + string( "/dat/leap/prep:" );
    path += amberhome + string( "/dat/leap/parm:" );
    path += amberhome + string( "/dat/leap/lib:" );
    path += amberhome + string( "/dat/leap/cmd:"  );
    mortenv().set_s( "path", path.c_str() );

    console()->process( "source leaprc.gleap" );

    assert( content().has( "_ions" ) );
    assert( content().has( "_solvents" ) );
}

void init_dlg( )
{
    
    static mortgtk::fixbond_t g_fixbond_dlg;
    static mortgtk::addhydr_t g_addhydr_dlg;
    static mortgtk::setpchg_t g_setpchg_dlg;
    static mortgtk::parmchk_t g_parmchk_dlg;
    static mortgtk::clsdata_t g_clsdata_dlg;
    static mortgtk::fastbld_t g_fastbld_dlg;
    static mortgtk::addions_t g_addions_dlg;
    static mortgtk::solvate_t g_solvate_dlg;
    static mortgtk::ffparam_t g_ffparam_dlg;
    static mortgtk::ribbdlg_t g_ribbdlg_dlg;

    static mortgtk::shidedlg_t g_shide_dlg;
    static mortgtk::labeldlg_t g_label_dlg;

}


int main( int argc, char** argv )
{
    const char* amberhome = getenv( "AMBERHOME" );
    
    if( amberhome == NULL )
    {
        std::cerr << "Environmental variable AMBERHOME must be set." << std::endl;
        return -1;
    }

    // g_thread_init( NULL );
    
    // gdk_threads_init();

    gtk_init( &argc, &argv );
    gtk_gl_init( &argc, &argv );

    console() = console_ptr( new gtk_console() );

    init_lib( amberhome );
    init_dlg( );


    main_window mainwin;
    mainwin.show();
    gtk_main();
    return 0;
}
