#ifndef GLEAP_LEAPSRC_WIDGET_SOLRADIO_HPP
#define GLEAP_LEAPSRC_WIDGET_SOLRADIO_HPP

#include "strradio.hpp"

namespace mortgtk
{
    
    class solradio_t : public strradio_t
    {
    public:

        solradio_t( );
        
        virtual ~solradio_t( );
        
        void init();
       
        virtual void onchange(GtkWidget* toggled);
 
        string getsol( );
        
    private:

        GtkWidget* m_maskentry;
    };
    

} // namespace mortgtk

        

#endif
