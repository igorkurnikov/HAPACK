#ifndef GLEAP_LEAPSRC_CLSDATA_HPP
#define GLEAP_LEAPSRC_CLSDATA_HPP

#include "basicdlg.hpp"

namespace mortgtk
{

    class clsdata_t : public basicdlg_t
    {
    public:

        clsdata_t();

        virtual ~clsdata_t();

        virtual command_ptr get_command();

        virtual basicdlg_ptr clone() const;
    };

        
} // namespace mortgtk





#endif
