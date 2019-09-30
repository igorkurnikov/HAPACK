#ifndef GLEAP_LEAPSRC_WIDGET_STRCOMBO_HPP
#define GLEAP_LEAPSRC_WIDGET_STRCOMBO_HPP

#include <vector>
#include <string>
#include <gtk/gtk.h>
#include <object.hpp>

namespace mortgtk
{
    using namespace mort;
    
    using std::vector;
    
    using std::string;

    class listener_i
    {
    public:

        virtual void onchange( GtkWidget* w ) = 0;
    };
    
    class strcombo_t
    {
    public:

        strcombo_t( const string& title );
        
        virtual ~strcombo_t( );

        void init( const string& s1, const string& s2, const string& s3 );
        
        void init( const vector<string>& ss );
        
        void init( const database_t& mdb );
        
        virtual void onchange(GtkWidget* w);
        
        operator GtkWidget*()
        {
            return m_self;
        }

        GtkWidget* vbox( )
        {
            return m_vbox;
        }

        string getstr() const;

        void listen( listener_i* l )
        {
            m_listeners.push_back( l );
        }
        
    private:

        vector<string> m_strings;

        string m_title;
        
        string m_result;

        GtkWidget* m_self;

        GtkWidget* m_vbox;
      
        vector< listener_i* > m_listeners;
  
    };

} // namespace mortgtk

        
#endif
