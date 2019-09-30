#ifndef GLEAP_LEAPSRC_ADDHYDR_HPP
#define GLEAP_LEAPSRC_ADDHYDR_HPP

#include "basicdlg.hpp"

namespace mortgtk
{
    
    class addhydr_t : public basicdlg_t
    {
    public:

        addhydr_t();
        
        virtual ~addhydr_t();
        
        virtual command_ptr get_command();
        
        virtual basicdlg_ptr clone() const
        {
            return basicdlg_ptr( new addhydr_t() );
        }
    };
    
        
} // namespace mortgtk




#endif
