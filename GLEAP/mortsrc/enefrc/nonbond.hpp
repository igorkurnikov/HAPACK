#include <vector>
#include "nonbond-pbc.hpp"

namespace mort
{
    using std::vector;

    class molecule_t;

    struct nonbond_t
    {
        nonbond_t(const molecule_t& m);
        void list_for_egb(double cut);
        void list_for_pbc(double cut);

	vector< vector<int> > list;
	vector< vector<int> > type;
	vector< vector<double> > dist;
	vector< vector<double> > dis2;

    private:

        void list_nbrs(const cell_t& a, const cell_t& b, double cut2);
        void list_nbrs_i(int i, const double* xi, const cell_t& cell, double cut2, int jstart);

        grid_t m_grid;
        numvec m_box;
        const molecule_t* m_pmol;

    };

} // namespace mort
