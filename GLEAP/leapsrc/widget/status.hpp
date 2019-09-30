#ifndef GLEAP_LEAPSRC_WIDGET_STATUS_HPP
#define GLEAP_LEAPSRC_WIDGET_STATUS_HPP


#include <string>
#include <gtk/gtkstatusbar.h>

namespace mortgtk
{
    using std::string;

    class gtk_status
    {
    public:
        // ctor/dtor
        explicit 
        gtk_status(GtkStatusbar* __ptr, const string& __context_name);
  
        gtk_status();
  
        virtual ~gtk_status() {}
  
        void init(GtkStatusbar* ptr, const string& context_name);
  
        guint push(const string& msg);
  
        guint getid() { return m_context_id; }

        void  pop();

        void  remove(guint msgid);
  
        string getdesc() { return m_context_name; }

    private:
  
        GtkStatusbar* m_pbar;
        string m_context_name;
        guint  m_context_id;
    };

} // namespace mortgtk
  

#endif

