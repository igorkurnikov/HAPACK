#include "status.hpp"


namespace mortgtk
{


    gtk_status::gtk_status(GtkStatusbar* pbar, const string& context_name )
    {
        init( pbar, context_name );
    }

    gtk_status::gtk_status( )
    {
        m_pbar = NULL;
        m_context_id = 0;
    }

    void gtk_status::init( GtkStatusbar* pbar, const string& context_name )
    {
        g_assert(pbar);
        m_pbar = pbar;
        m_context_name = context_name;
        m_context_id   = gtk_statusbar_get_context_id( pbar, m_context_name.c_str());
    }

    guint gtk_status::push(const string& msg)
    {
        g_assert(m_pbar);
        return gtk_statusbar_push(m_pbar, m_context_id, msg.c_str());
    }

    void gtk_status::pop()
    {
        g_assert(m_pbar);
        gtk_statusbar_pop(m_pbar, m_context_id);
    }

    void gtk_status::remove(guint msgid)
    {
        g_assert(m_pbar);
        gtk_statusbar_remove(m_pbar, m_context_id, msgid);
    }


} // namespace mortgtk

 
