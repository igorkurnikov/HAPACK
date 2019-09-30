#ifndef GLEAP_LEAPSRC_SETPCHG_HPP
#define GLEAP_LEAPSRC_SETPCHG_HPP

#include "basicdlg.hpp"

namespace mortgtk
{
    
    class setpchg_t : public basicdlg_t
    {
    public:

        setpchg_t();
        
        virtual ~setpchg_t();
        
        command_ptr get_command();
        
        virtual basicdlg_ptr clone() const
        {
            return basicdlg_ptr( new setpchg_t() );
        }
        
    };
    
        
} // namespace mortgtk




#endif
