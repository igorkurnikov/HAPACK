#ifndef GLEAP_LEAPSRC_RIBBDLG_HPP
#define GLEAP_LEAPSRC_RIBBDLG_HPP

#include "basicdlg.hpp"

namespace mortgtk
{

    class ribbdlg_t : public basicdlg_t
    {
    public:

        ribbdlg_t();

        virtual ~ribbdlg_t();

        virtual command_ptr get_command();

        virtual basicdlg_ptr clone() const;
    };

        
} // namespace mortgtk





#endif
