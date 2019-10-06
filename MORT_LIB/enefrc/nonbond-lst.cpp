#include <object.hpp>
#include <ambfmt.hpp>
#include "nonbond.hpp"

namespace mort
{

    nonbond_t::nonbond_t(const molecule_t& m)
    	: list( m.natom() ), type( m.natom() ),
	  dist( m.natom() ), dis2( m.natom() ), m_box(3)
    {
        m_pmol = &m;
    }

    void nonbond_t::list_for_egb(double cut)
    {
        const vector<double>& pos = get_vvec(*m_pmol, ATOM, POSITION);

        double cut2 = cut*cut;

        excl_t excl;
        exclude( *m_pmol, excl, 3 );

        atomiter_t ai = m_pmol->atom_begin();
	for( ; ai != m_pmol->atom_end(); ++ai)
	{
	    int i = ai->absid();
            vector<int>::iterator excl_i = excl.list[i].begin();
        
	    atomiter_t aj = ai+1;
	    for( ; aj != m_pmol->atom_end(); ++aj)
	    {
	        int j = aj->absid();
                double dx = pos[3*i  ] - pos[3*j  ];
		double dy = pos[3*i+1] - pos[3*j+1];
		double dz = pos[3*i+2] - pos[3*j+2];
		double d2 = dx*dx + dy*dy + dz*dz;
		if( d2 < cut2 )
		{
		    list[i].push_back(j);
		    dis2[i].push_back(d2);
		    dist[i].push_back(sqrt(d2));
		    if( excl_i != excl.list[i].end() && aj->get_i(ID)==*excl_i )
		    {
		        ++excl_i;
			type[i].push_back(0);
		    }
		    else
		    {
		        type[i].push_back(1);
		    }
	        }
	    }
	}
    }

    void nonbond_t::list_nbrs_i(int i, const double* xi, const cell_t& cell, double cut2, int jstart)
    {
        for(int j=jstart; j < cell.size(); ++j)
        {
            int aj = cell.atomids[j];
            const double* xj = &cell.imgcrds[0] + 3*j;
            
            double d2 = img_dis2(xi, xj, m_box);
            if( d2 > 0.0 && d2 < cut2 )
            {
                list[i].push_back(aj);
		dis2[i].push_back(d2);
                dist[i].push_back(sqrt(d2));
            }
        }
    }


    void nonbond_t::list_nbrs(const cell_t& cell_a, const cell_t& cell_b, double cut2)
    {
        bool same_cell = (&cell_a==&cell_b);

        const int* ai = &cell_a.atomids[0];
        const double* xi = &cell_a.imgcrds[0];
        for(int i=0;i < cell_a.size();++i, ++ai, xi+=3)
        {    
            int jstart = same_cell? (i+1) : 0;
            list_nbrs_i(*ai, xi, cell_b, cut2, jstart);
        }
    }


    void nonbond_t::list_for_pbc(double cut)
    {
        int natom = m_pmol->natom();
        m_box = m_pmol->get_v(BOX);
	m_grid.init(m_box, cut+2.0);

	const vector<double>& pos = get_vvec(*m_pmol, ATOM, POSITION);
	vector<double> frac( 3*m_pmol->natom() );
        get_frac(3*natom, &pos[0], m_box, &frac[0]);

	m_grid.assign_atom(frac.size(), &frac[0]);

	excl_t excl;
	exclude(*m_pmol, excl, 3);

        double cut2 = cut*cut;

        for(int i=0; i < (int)m_grid.size(); ++i)
        {
            for(int j=0; j < (int)m_grid[i].nbrcells.size(); ++j)
            {
                int k = m_grid[i].nbrcells[j];
                assert( k>=0 && k < (int)m_grid.size() );
                list_nbrs(m_grid[i], m_grid[k], cut2);
            }
        }
    }

} // namespace mort

