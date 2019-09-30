#ifndef GLEAP_LEAPSRC_WIDGET_ADDIONS_HPP
#define GLEAP_LEAPSRC_WIDGET_ADDIONS_HPP

#include "basicdlg.hpp"

namespace mortgtk
{

    class addions_t : public basicdlg_t
    {
    public:

        addions_t();
        
        virtual ~addions_t();

        virtual command_ptr get_command( );
        
        virtual basicdlg_ptr clone( ) const
        {
            return basicdlg_ptr( new addions_t() );
        }

    private:

        strcombo_t m_ioname;
        
        strcombo_t m_ionumb;
    };
    
} // namespace mortgtk

#endif
