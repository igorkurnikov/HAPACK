#ifndef AMBER_GLEAP_GTKAREA_HPP
#define AMBER_GLEAP_GTKAREA_HPP

#include <drawing.hpp>

namespace amber
{
    
    class gtk_glarea : public mort::drawing_i
    {
    public:

        gtk_glarea( GtkWidget* widget );
        
        ~gtk_glarea();

        virtual void glbegin();
        
        virtual void glend();

        virtual void expose( );
       
    private:

        GtkWidget* m_widget;
    };
    
}

#endif
