#include "strcombo.hpp"


namespace mortgtk
{
    void S_on_strcombo_change( GtkWidget* w, gpointer user_data )
    {
        strcombo_t* sc = (strcombo_t*)user_data;
        sc->onchange( w );
    };


    
    strcombo_t::strcombo_t(const string& title)
        :m_title( title )
    {
    }
    

    strcombo_t::~strcombo_t()
    {
    }

    void strcombo_t::init( const database_t& mdb )
    {
        vector<string> strs;
        
        database_t::const_iterator i = mdb.begin();
        for( ; i != mdb.end(); ++i )
        {
            string str = i->first;
            
            if( str[0]!='_' && ismol(i->second) )
            {
                strs.push_back( str );
            }
        }
        
        init( strs );
    }

    void strcombo_t::init( const string& s1, const string& s2, const string& s3 )
    {
        vector<string> strs;
        strs.push_back( s1 );
        strs.push_back( s2 );
        strs.push_back( s3 );
        init( strs );
    }
    
    void strcombo_t::init( const vector<string>& strs )
    {
        m_self = gtk_frame_new( m_title.c_str() );
        gtk_widget_show( m_self );
        
        m_vbox = gtk_vbox_new( true, 0 );
        gtk_container_add( GTK_CONTAINER(m_self), m_vbox );
	gtk_widget_show( m_vbox );

        GtkWidget* button = gtk_combo_box_new_text();
        gtk_box_pack_start( GTK_BOX(m_vbox), button, true, true, 0 );
        gtk_widget_show( button );

        for( int i=0; i < strs.size(); ++i )
        {
            gtk_combo_box_append_text( GTK_COMBO_BOX(button), strs[i].c_str() );
        }        
        g_signal_connect( G_OBJECT(button), "changed", G_CALLBACK(S_on_strcombo_change), this );
        
 

        GtkTreeModel* model = gtk_combo_box_get_model( GTK_COMBO_BOX(button) );

        GtkTreeIter iter;

        if( gtk_tree_model_get_iter_first( model, &iter ) )
        {
            gchar* text;
            gtk_tree_model_get( model, &iter, 0, &text, -1 );
            gtk_combo_box_set_active_iter( GTK_COMBO_BOX(button), &iter);
        }
    }
    
    void strcombo_t::onchange(GtkWidget* w)
    {
    	GtkTreeIter iter;
	if( gtk_combo_box_get_active_iter( GTK_COMBO_BOX(w), &iter) )
	{
	    GtkTreeModel* model = gtk_combo_box_get_model( GTK_COMBO_BOX(w) );
            gchar* select;
	    gtk_tree_model_get( model, &iter, 0, &select, -1 );
            m_result = select;
            g_free(select);
        }

        for( int i=0; i < m_listeners.size(); ++i )
        {
            m_listeners[i]->onchange( w );
        }

    }
    
    std::string strcombo_t::getstr( ) const
    {
        return m_result;
    }
    

} // namespace mortgtk

 
