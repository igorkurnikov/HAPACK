#include <common.hpp>
#include <object.hpp>

namespace mort
{

    namespace tripos
    {

        using namespace std;
        
        void read_resd(istream& is, molecule_t& mol)
        {
	    resditer_t ri = mol.resd_begin();
	    for( ; ri != mol.resd_end(); ++ri )
	    {
	        //read_iparm( is, *ri, ID );
		//read_sparm( is, *ri, NAME );
                is.ignore( MAX_LINE_WIDTH, '\n' );
	    }
        }

        void write_resd(ostream& os, const molecule_t& mol)
        {
            if( mol.nresd() == 0 )
            {
                return;
            }

            os << "@<TRIPOS>SUBSTRUCTURE" << endl;

            resditer_t resd = mol.resd_begin();
            for( ; resd != mol.resd_end(); ++resd )
            {
                write_iparm(os, *resd, ID, "%8d ");
                write_sparm(os, *resd, TYPE, "%8s ");
                os << std::endl;
            }
        }

    } // namespace tripos
    
} // namespace mort



