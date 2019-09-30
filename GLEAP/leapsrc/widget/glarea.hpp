#ifndef AMBER_GLEAP_GTKAREA_HPP
#define AMBER_GLEAP_GTKAREA_HPP

#include <mortgl.hpp>

namespace mortgtk
{
    class gtk_glarea : public mort::drawing_t
    {
    public:

        gtk_glarea( );
        
        virtual ~gtk_glarea();

        virtual void glbegin();
        
        virtual void glend();

        virtual void expose();

        operator GtkWidget* () 
        {
            return m_self;
        }

        void popup_menu( GdkEventButton* event);

    private:

        GtkWidget* m_self;

        GtkWidget* m_popmenu;
        
        GtkUIManager* m_manager;
        
        GtkActionGroup* m_actions;
        
    };



}

#endif
