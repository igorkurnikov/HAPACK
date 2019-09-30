#ifndef GLEAP_LEAPSRC_WIDGET_MSKCOMBO_HPP
#define GLEAP_LEAPSRC_WIDGET_MSKCOMBO_HPP

#include "strcombo.hpp"

namespace mortgtk
{
    
    class mskcombo_t : public strcombo_t
    {
    public:

        mskcombo_t( const string& title );
        
        virtual ~mskcombo_t();

        void init( );        

        virtual void onchange( GtkWidget* w );
        
        string getmsk();

        using strcombo_t::init;

    private:

        GtkWidget* m_maskentry;

        bool m_entry_on;
        
    };
    
    atomvec_t unified_mask_atom(const molecule_t& m, const string& mask );
    

} // namespace mortgtk
        
        
#endif
