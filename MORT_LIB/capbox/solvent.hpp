#ifndef MORTSRC_CAPBOX_SOLVENT_HPP
#define MORTSRC_CAPBOX_SOLVENT_HPP

#include <common.hpp>
#include <object.hpp>

namespace mort
{
    class molecule_t;

    class solvent_t
    {
    public:

        solvent_t( const molecule_t& m );

        int natom() const { return m_pmol->natom(); }

        int nresd() const { return m_pmol->nresd(); }
         
        const molecule_t& getfull() const 
        {
            return *m_pmol;
        }

        const molecule_t& getresd(int rid) const;

        const vector<double>& getvdwr() const
        {
            return m_vdwr;
        }

        const vector<double>& getcord() const
        {
            return get_vvec( *m_pmol, ATOM, POSITION );
        }

        const numvec center() const { return mort::center(*m_pmol); }

        const numvec extent() const { return subvec(m_pmol->get_v(BOX), 0, 3); }

    private:

        const molecule_t* m_pmol;

        vector<double> m_vdwr;

        mutable molecule_t m_resd;
    };

} // namespace mort

#endif

