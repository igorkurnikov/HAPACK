#ifndef GLEAP_LEAPSRC_PARMCHK_HPP
#define GLEAP_LEAPSRC_PARMCHK_HPP

#include "basicdlg.hpp"

namespace mortgtk
{
    
    class parmchk_t : public basicdlg_t
    {
    public:

        parmchk_t();
        
        virtual ~parmchk_t();
        
        command_ptr get_command();
        
        virtual basicdlg_ptr clone() const
        {
            return basicdlg_ptr( new parmchk_t() );
        }
        
    };
    
        
} // namespace mortgtk




#endif
