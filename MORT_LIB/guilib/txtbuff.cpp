#include <pango/pango-font.h>
#include "leaplog.hpp"
#include "txtbuff.hpp"

namespace amber
{
    gtk_txtbuff::gtk_txtbuff( GtkWidget* view )
    {
        m_view = view;
        
        m_buff = gtk_text_view_get_buffer( GTK_TEXT_VIEW( view ) );
    }

    gtk_txtbuff::~gtk_txtbuff()
    {
    }

    void gtk_txtbuff::markon( postype_e pos, const string& name, bool flag )
    {
        assert( pos == END );
        
        GtkTextIter end;

        gtk_text_buffer_get_end_iter( m_buff, &end );

        gtk_text_buffer_create_mark( m_buff, name.c_str(), &end, flag );
    }
    

    void gtk_txtbuff::insert( postype_e pos, const string& str, const string& suffix )
    {
        assert( pos == END );
        
        GtkTextIter end;

        string line = str + suffix;

        gtk_text_buffer_get_end_iter( m_buff, &end );

        gtk_text_buffer_insert( m_buff, &end, line.c_str(), line.length() );

        gtk_text_view_scroll_mark_onscreen( GTK_TEXT_VIEW( m_view ), gtk_text_buffer_get_insert( m_buff ) );
        
    }

    void gtk_txtbuff::remove( const string& name, postype_e pos )
    {
        assert( pos == END );

        GtkTextIter begin, end;
        
        GtkTextMark* mark = gtk_text_buffer_get_mark( m_buff, name.c_str() );
        
        assert( mark != NULL );

        gtk_text_buffer_get_iter_at_mark( m_buff, &begin, mark );
        
        gtk_text_buffer_get_end_iter( m_buff, &end );
        
        gtk_text_buffer_delete( m_buff, &begin, &end );
    }
    
    string gtk_txtbuff::copy( const string& name, postype_e pos )
    {
        assert( pos == END );
            
        GtkTextIter begin, end;

        GtkTextMark* mark = gtk_text_buffer_get_mark( m_buff, name.c_str() );
        
        assert( mark != NULL );
    
        gtk_text_buffer_get_iter_at_mark( m_buff, &begin, mark );

        gtk_text_buffer_get_end_iter( m_buff, &end );

        return gtk_text_buffer_get_text( m_buff, &begin, &end, false );
    }
    
} // namespace amber


    
