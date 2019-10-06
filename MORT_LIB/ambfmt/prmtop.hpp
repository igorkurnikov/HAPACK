#ifndef MORTSRC_FRCFLD_PRMTOP_HPP
#define MORTSRC_FRCFLD_PRMTOP_HPP

#include <object.hpp>
#include <common.hpp>
#include "energee.hpp"
#include "exclude.hpp"

namespace mort
{
    namespace prmtop
    {
        void write_les( ostream& os, const molecule_t& m );
        void write_string( ostream& os, const string& flag, const string& str );
        void write_sarray( ostream& os, const string& flag, const string& str, int nitem );
        void write_iarray( ostream& os, const string& flag, const int* ia, int nitem );
        void write_darray( ostream& os, const string& flag, const double* da, int nitem );
        void write_prmtop( ostream& os, const nabparm_t& prm );
	
		namespace amoeba
		{
			void setup_adjust( const molecule_t& mol, excl_t& excl, const molecule_t& poleff );
			void create_atomic_frame( morf_t& atom, const molecule_t& poleff, vector< vector<int> >& chirials, vector< vector<int> >& regulars );
		}
 
    } // namespace prmtop

    void write_amber_prmtop( ostream& os, molecule_t& m, const molecule_t& ffp );

} // namespace mort

#endif

