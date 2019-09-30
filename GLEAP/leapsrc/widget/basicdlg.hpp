#ifndef GLEAP_LEAPSRC_BASICDLG_HPP
#define GLEAP_LEAPSRC_BASICDLG_HPP


#include <gtk/gtk.h>

#include <guilib.hpp>

#include "strcombo.hpp"

namespace mortgtk
{

    using namespace mort;

    
    class basicdlg_t : public listener_i
    {
    public:

        basicdlg_t( const char* title);

        virtual ~basicdlg_t();
        
        int run() 
        {
            return gtk_dialog_run(m_self);
        }
        
        GtkWidget* vbox() 
        { 
            return m_self->vbox; 
        }

        string molname() const;

        virtual void onchange( GtkWidget* w ) 
        {
        }

        virtual command_ptr get_command()=0;

        virtual shared_ptr<basicdlg_t> clone() const=0;
        
    private:

        GtkDialog* m_self;

        strcombo_t m_molname;
        
    };

    typedef shared_ptr<basicdlg_t> basicdlg_ptr;

    class dlgfactory_t
    {
    private:

        dlgfactory_t();
        
        virtual ~dlgfactory_t();

    public:

        static dlgfactory_t& inst();
        
        static void add( const string& name, basicdlg_t* impl );
        
        static void run( const string& name );

    private:

        map<string, basicdlg_t*> m_data;
    };
    







} // namespace mortgtk

#endif



    
