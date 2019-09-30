#include "solvent.hpp"


namespace mort
{

    solvent_t::solvent_t( const molecule_t& m )
        : m_vdwr( m.natom() ),
          m_resd( unmerge(m.resds()[0]) )
    {
        m_pmol = &m;
        set_vdwr( m, m_vdwr );
    }

    const molecule_t& solvent_t::getresd(int rid) const
    {
        int rsize = natom() / nresd();
        int begin = rid*rsize*3;

        vector<double> const& full = get_vvec( *m_pmol, ATOM, POSITION );
        vector<double>      & part = get_vvec(  m_resd, ATOM, POSITION );
        for( int i=0; i < rsize*3; ++i )
        {
            part[i] = full[begin+i];
        }
    
        return m_resd;
    }
} 
