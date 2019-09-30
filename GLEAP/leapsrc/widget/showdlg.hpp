#ifndef GLEAP_LEAPSRC_WIDGET_SHOWDLG_HPP
#define GLEAP_LEAPSRC_WIDGET_SHOWDLG_HPP

#include <guilib.hpp>
#include "glarea.hpp"


class main_window;

namespace mortgtk
{

    using namespace mort;

    class showdlg_t
    {
    public:

        showdlg_t(main_window& mainwin);

        virtual ~showdlg_t();

	int run();

	command_ptr get_command();

        void on_operate_change(GtkWidget* toggled);
	
	void on_atomask_change(GtkWidget* widget);

	void on_molname_change(GtkWidget* widget);

    private:

        string m_operate;
        string m_atomask;
        string m_molname;

        GtkDialog* m_self;
	GtkWidget* m_oper_on;
	GtkWidget* m_oper_off;
	GtkWidget* m_oper_only;
	GtkWidget* m_maskentry;

        database_t* m_data;
        drawing_t* m_area;
    };

} // namespace mortgtk

#endif
