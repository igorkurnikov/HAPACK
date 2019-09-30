#ifndef GLEAP_LEAPSRC_WIDGET_SHIDEDLG_HPP
#define GLEAP_LEAPSRC_WIDGET_SHIDEDLG_HPP

#include <guilib.hpp>
#include "basicdlg.hpp"
#include "strcombo.hpp"
#include "strradio.hpp"
#include "mskcombo.hpp"

namespace mortgtk
{

    using namespace mort;

    class shidedlg_t : public basicdlg_t
    {
    public:

        shidedlg_t();

        virtual ~shidedlg_t();

	virtual command_ptr get_command();

        virtual basicdlg_ptr clone() const
        {
            return basicdlg_ptr( new shidedlg_t() );
        }

    private:

        strradio_t m_operate;
        strradio_t m_graphic;
        mskcombo_t m_atomask;

    };

} // namespace mortgtk

#endif
