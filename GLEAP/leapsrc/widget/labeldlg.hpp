#ifndef GLEAP_LEAPSRC_WIDGET_LABELDLG_HPP
#define GLEAP_LEAPSRC_WIDGET_LABELDLG_HPP

#include <guilib.hpp>
#include "basicdlg.hpp"
#include "strradio.hpp"
#include "strcombo.hpp"

namespace mortgtk
{

    using namespace mort;

    class labeldlg_t : public basicdlg_t
    {
    public:

        labeldlg_t();

        virtual ~labeldlg_t();

	virtual command_ptr get_command();

        virtual basicdlg_ptr clone() const
        {
            return basicdlg_ptr( new labeldlg_t() );
        }

    private:

        strradio_t m_level;
        strcombo_t m_param;


    };

} // namespace mortgtk

#endif
