#ifndef MORT_OBJECT_SUBST_HPP
#define MORT_OBJECT_SUBST_HPP

#include "finder.hpp"

namespace mort
{
    class subst_finder : public finder_t
    {
    public:

        subst_finder(const morf_t& atom, const molecule_t& sub);
        
        virtual ~subst_finder();
        
        virtual int predict(const atomvec_t& path);

        virtual void list_next(const atomvec_t& path, atomvec_t& next);

        void find_next(atomiter_t begin, atomiter_t end, const morf_t& s_atom, const atmvec& path, atmvec& next );
 
    private:

        bool match(const morf_t& s_atom, const morf_t& m_atom );

        const molecule_t& m_sub;
    };

    bool has_subst(const morf_t& atom, const molecule_t& sub);

    bool find_subst(const morf_t& atom, const molecule_t& sub, atomvec_t& atoms);

}// namespace mort

#endif


