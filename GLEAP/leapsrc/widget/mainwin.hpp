#ifndef GLEAP_WIDGET_MAINWIN_HPP
#define GLEAP_WIDGET_MAINWIN_HPP

#include <gtk/gtk.h>
#include <object.hpp>
#include "status.hpp"
#include "glarea.hpp"
#include "macros.hpp"
#include "editbar.hpp"

using namespace mort;
using namespace mortgtk;

class main_window
{
public:

    main_window();

    virtual ~main_window();

    void show();

    void init();

    void init_manager();

    void init_menubar(GtkBox* pbox);

    void init_toolbar(GtkBox* pbox);

    void init_viewpanel(GtkBox* pbox);

    void init_statusbar(GtkBox* pbox);

#define DECLARE_MENU_SIGNAL_LISTENER(menu)	\
DECLARE_CLASS_ACTION_LISTENER(menu)		\
DECLARE_CLASS_MENU_HINTS_LISTENER(menu)

    DECLARE_MENU_SIGNAL_LISTENER(file_open)
    DECLARE_MENU_SIGNAL_LISTENER(file_close)
    DECLARE_MENU_SIGNAL_LISTENER(file_save)
    DECLARE_MENU_SIGNAL_LISTENER(file_quit)

    DECLARE_MENU_SIGNAL_LISTENER(edit_fixbond)
    DECLARE_MENU_SIGNAL_LISTENER(edit_addhydr)
    DECLARE_MENU_SIGNAL_LISTENER(edit_setpchg)
    DECLARE_MENU_SIGNAL_LISTENER(edit_parmchk)

    DECLARE_MENU_SIGNAL_LISTENER(edit_fastbld)
    DECLARE_MENU_SIGNAL_LISTENER(edit_addions)
    DECLARE_MENU_SIGNAL_LISTENER(edit_solvate)
    DECLARE_MENU_SIGNAL_LISTENER(edit_ffparam)

    DECLARE_MENU_SIGNAL_LISTENER(view_shide)
    DECLARE_MENU_SIGNAL_LISTENER(view_label)
    DECLARE_MENU_SIGNAL_LISTENER(view_color)
    DECLARE_MENU_SIGNAL_LISTENER(view_ribbon)
    DECLARE_MENU_SIGNAL_LISTENER(view_surface)

#undef DECLARE_MENU_SIGNAL_LISTENER

private:

    gtk_status m_status;

    GtkWidget* m_pmainwin;
    GtkWidget* m_pmenubar;
    GtkWidget* m_ptoolbar;
    GtkWidget* m_pviewpanel;
    GtkWidget* m_pstatusbar;

    GtkUIManager* m_pmanager;
    
    molecule_ptr m_curt;

    editbar_t m_editbar;

    shared_ptr<gtk_glarea> m_glarea;

};

int popup_message( const string& msg );




#endif
