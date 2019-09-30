#ifndef MHTL_SYBYL_BOND_HPP
#define MHTL_SYBYL_BOND_HPP

#include <map>
#include <common.hpp>
#include <object.hpp>

namespace mort
{
    namespace tripos
    {
        using std::pair;
        using std::find_if;
	using std::make_pair;
	typedef pair< int, int > ipair_t;

        void read_bond( istream& is, molecule_t& mol )
        {
	    bonditer_t bi = mol.bond_begin();
            bonditer_t be = mol.bond_end();
	    for( ; bi != be; ++bi )
	    {
	        read_iparm( is, *bi, ID );
		int first,second;
		is >> first >> second;

                atom_t a1 = *find_if( mol.atom_begin(), mol.atom_end(),  iparm_cmper1(SEQUENCE, first) );
                atom_t a2 = *find_if( mol.atom_begin(), mol.atom_end(),  iparm_cmper1(SEQUENCE, second) );

                resditer_t r1 = a1.resd_begin();
                resditer_t r2 = a2.resd_begin();

                int aid1 = a1.get_i(ID);
                int aid2 = a2.get_i(ID);

                int rid1 = r1->relid();
                int rid2 = r2->relid();
                if( rid1 != rid2 )
                {
                    if( rid1 < rid2 )
                    {
                        r1->set_i( TAIL, aid1 );
                        r2->set_i( HEAD, aid2 );
                    }
                    else
                    {
                        r1->set_i( HEAD, aid1 );
                        r2->set_i( TAIL, aid2 );
                    }
                }

                a1.connect( a2 );
                a2.connect( a1 );
               
                a1.connect( *bi );
                a2.connect( *bi );

                bi->connect( a1 );
                bi->connect( a2 );

                read_sparm( is, *bi, TYPE );
                is.ignore(MAX_LINE_WIDTH, '\n');
            } 
        } 


        void write_bond( ostream& os, const molecule_t& mol )
        {
            os << "@<TRIPOS>BOND" << std::endl;

            bonditer_t bond = mol.bond_begin();
            for( ; bond != mol.bond_end(); ++bond )
            {
                //write_iparm(os, *bond, ID, "%8d ");
                os << format("%8d ") % (bond->relid()+1);

                ipair_t ap = any_cast< ipair_t >( bond->get_a(ATOMPAIR) );
                os << format( "%8d " ) % std::min( ap.first, ap.second );
		os << format( "%8d " ) % std::max( ap.first, ap.second );

                write_iparm(os, *bond, ORDER, "%8d ");
                os << std::endl;
            }
        }
    
    } // namespace tripos
    
} // namespace mort

#endif
