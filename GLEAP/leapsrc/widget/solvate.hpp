#ifndef GLEAP_LEAPSRC_WIDGET_SOLVATE_HPP
#define GLEAP_LEAPSRC_WIDGET_SOLVATE_HPP

#include "solradio.hpp"
#include "strcombo.hpp"
#include "basicdlg.hpp"

namespace mortgtk
{
    
    class solvate_t : public basicdlg_t
    {
    public:

        solvate_t();
        
        virtual ~solvate_t();
        
        virtual command_ptr get_command();
        
        virtual basicdlg_ptr clone() const
        {
            return basicdlg_ptr( new solvate_t() );
        }

    private:

        solradio_t m_soltype;
        strcombo_t m_solvent;
        GtkWidget* m_solsize;
    };
    

} // end 

        


#endif
