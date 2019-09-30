#ifndef GLEAP_LEAPSRC_FIXBOND_HPP
#define GLEAP_LEAPSRC_FIXBOND_HPP

#include "basicdlg.hpp"

namespace mortgtk
{
    
    class fixbond_t : public basicdlg_t
    {
    public:

        fixbond_t();
        
        virtual ~fixbond_t();
        
        virtual command_ptr get_command();

        virtual basicdlg_ptr clone() const;
    };
    
        
} // namespace mortgtk




#endif
