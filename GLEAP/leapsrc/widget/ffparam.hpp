#ifndef GLEAP_LEAPSRC_WIDGET_FFPARAM_HPP
#define GLEAP_LEAPSRC_WIDGET_FFPARAM_HPP

#include "basicdlg.hpp"
#include "strradio.hpp"

namespace mortgtk
{
    
    class ffparam_t : public basicdlg_t
    {
    public:

        ffparam_t( );
        
        virtual ~ffparam_t( );
        
        virtual void onchange( GtkWidget* w );
        
        virtual command_ptr get_command();

        virtual basicdlg_ptr clone() const
        {
            return basicdlg_ptr( new ffparam_t() );
        }

    private:

        strradio_t m_fftype;
        
        GtkWidget* m_topfil;
        
        GtkWidget* m_xyzfil;
    };
    

} // namespace mortgtk



#endif
