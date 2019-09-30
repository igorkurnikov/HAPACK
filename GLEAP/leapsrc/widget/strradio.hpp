#ifndef GLEAP_LEAPSRC_WIDGET_STRRADIO_HPP
#define GLEAP_LEAPSRC_WIDGET_STRRADIO_HPP

#include <vector>
#include <string>

#include <gtk/gtk.h>


namespace mortgtk
{
    using std::vector;
    
    using std::string;
    
    class strradio_t
    {
    public:

        strradio_t( const string& title );

        virtual ~strradio_t();

        void init( const string& s1, const string& s2 );
        
        void init( const string& s1, const string& s2, const string& s3 );

        void init( const string& s1, const string& s2, const string& s3, const string& s4 );

        void init( const vector<string>& strs );

        virtual void onchange(GtkWidget* toggled);

        string getstr();

        operator GtkWidget*() 
        {
            return m_self;
        }

        GtkWidget* vbox()
        {
            return m_vbox;
        }

    private:

        void init( );
        
    private:

        GtkWidget* m_self;

        GtkWidget* m_vbox;

        string m_title;

        string m_result;

        vector<string> m_strings;
        
        vector<GtkWidget*> m_widgets;
    };

} // namespace mortgtk



#endif
