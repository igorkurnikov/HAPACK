#ifndef GLEAP_LEAPSRC_WIDGET_FASTBLD_HPP
#define GLEAP_LEAPSRC_WIDGET_FASTBLD_HPP


#include "basicdlg.hpp"

namespace mortgtk
{
    
    class fastbld_t : public basicdlg_t
    {
    public:

        fastbld_t();
        
        virtual ~fastbld_t();
        
        virtual command_ptr get_command();
        
        virtual basicdlg_ptr clone() const
        {
            return basicdlg_ptr( new fastbld_t() );
        }
        
    private:

    };
    

} // namespace mortgtk


        











#endif
