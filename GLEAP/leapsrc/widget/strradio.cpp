#include <cassert>
#include "strradio.hpp"


namespace mortgtk
{
    void S_on_strradio_change(GtkWidget* w, gpointer user_data)
    {
        strradio_t* pr = (strradio_t*)user_data;
        pr->onchange( w );
    }
    
    strradio_t::strradio_t(const string& title)
        : m_title( title )
    {
    }
    

    strradio_t::~strradio_t()
    {
    }

    void strradio_t::init( const string& s1, const string& s2 )
    {
        assert( m_strings.size()==0 );
        m_strings.push_back( s1 );
        m_strings.push_back( s2 );
        init();
    }

    void strradio_t::init( const string& s1, const string& s2, const string& s3 )
    {
        assert( m_strings.size()==0 );
        m_strings.push_back( s1 );
        m_strings.push_back( s2 );
        m_strings.push_back( s3 );
        init();
    }

    void strradio_t::init( const string& s1, const string& s2, const string& s3, const string& s4 )
    {
        assert( m_strings.size()==0 );
        m_strings.push_back( s1 );
        m_strings.push_back( s2 );
        m_strings.push_back( s3 );
        m_strings.push_back( s4 );
        init();
    }

    void strradio_t::init( )
    {
        assert( m_strings.size()>0 );
        
        m_self = gtk_frame_new( m_title.c_str() );
        gtk_widget_show( m_self );


        m_vbox = gtk_vbox_new(true, 0);
        gtk_widget_show( m_vbox );
        gtk_container_add( GTK_CONTAINER(m_self), m_vbox );

        m_widgets.resize( m_strings.size() );
        for( int i=0; i < m_strings.size(); ++i )
        {
            GtkRadioButton* parent = (i==0)? NULL : GTK_RADIO_BUTTON(m_widgets[0]);
            m_widgets[i] = gtk_radio_button_new_with_label_from_widget( parent, m_strings[i].c_str() );
            gtk_box_pack_start( GTK_BOX(m_vbox), m_widgets[i],   true, true, 0);
            gtk_widget_show( m_widgets[i] );
            g_signal_connect( G_OBJECT(m_widgets[i]),   "toggled", G_CALLBACK(S_on_strradio_change), this );
        }

        m_result = m_strings[0];
    }
    
    void strradio_t::onchange(GtkWidget* toggled)
    {
        int i = std::find(m_widgets.begin(), m_widgets.end(), toggled) - m_widgets.begin();
        assert( i < m_strings.size() );
        m_result = m_strings[i];
    }

    string strradio_t::getstr()
    {
        return m_result;
    }
    

} // namespace mortgtk



    

        

