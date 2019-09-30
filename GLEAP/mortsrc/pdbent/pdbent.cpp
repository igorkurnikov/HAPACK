#include <common.hpp>
#include <object.hpp>
#include <pdbent.hpp>
#include "second.hpp"
#include "atom.hpp"
#include "bond.hpp"
#include "resd.hpp"

namespace mort
{
    using pdbent::second_s;

    bool read_pdb( istream& is, molecule_t& mol )
    {
		int icount = 0;
        molecule_t m;
        resd_t rprev( m, -1 );
        vector<second_s> seconds;

        string keyw = peek_keyw( is );
        while( is && ! keyw.empty() )
        {
			icount++;
//			printf(" icount = %d %s \n",icount,keyw.c_str());  // IGOR_TMP
            if( keyw == "ATOM" || keyw == "HETA" )
            {
                resd_t rcurt = pdbent::read_resd( is, m );
                pdbent::setpos( rprev, rcurt );
                rprev = rcurt;
            }
            else if( keyw.find("TER")==0 || keyw.find("END")==0 )
            {
                pdbent::setend( rprev );
                is.ignore( MAX_LINE_WIDTH, '\n' );
            }    
            else if( keyw == "HELI" || keyw == "SHEE" )
            {
                seconds.push_back( pdbent::second_s() );
                pdbent::read_second( is, seconds.back() );
            }
            else if( keyw == "CONE" )
            {
                pdbent::read_bond( is, m );
            }
            else
            {
                is.ignore(MAX_LINE_WIDTH, '\n');
            }

            keyw = peek_keyw( is );    
        }    

        pdbent::setend( rprev );

        std::for_each( m.resd_begin(), m.resd_end(), iparm_setter1(SECOND, LOOP) );
        for( int i=0; i < (int)seconds.size(); ++i )
        {
            apply_second( seconds[i], m );
        }

        mol.swap(m);
        return true;
    }

    using std::map;
    using std::set;

    bool write_pdb( ostream& os, const molecule_t& mol )
    {
        funstack_t::push( "write_pdb" );

        if( mol.nresd() > 0 )
        {
            resditer_t i = mol.resd_begin();
            resditer_t e = mol.resd_end();
            for( ; i != mol.resd_end(); ++i )
            {
                pdbent::write_resd( os, *i );
                resditer_t j = i+1;
                if( j != mol.resd_end() && !bond_t::has(*i, *j) )
                    os << "TER" << std::endl;
            }
        }
        else
        {
            atomiter_t ai = mol.atom_begin();
            atomiter_t ae = mol.atom_end();
            for( ; ai != ae; ++ai )
            {
                pdbent::write_atom( os, *ai );
            }
        }

       
        set<int> atoms;
        pdbent::mark_bond( mol, atoms );

        set<int>::iterator i = atoms.begin();
        for( ; i != atoms.end(); ++i )
        {
            os << "CONECT" << format( "%5d" ) % (*i);

            std::set<int> nbrs;
            atom_t a = *( mol.atom_begin()+(*i-1) );
            atomiter_t ni = a.atom_begin();
            for( ; ni != a.atom_end(); ++ni )
            {
                nbrs.insert( ni->get_i(ID) );
            }

            std::set<int>::iterator j = nbrs.begin();
            for( ; j != nbrs.end(); ++j )
            {
                os << format( "%5d" ) % (*j);
            }
            os << std::endl;
        }

        os << "END" << std::endl;

        funstack_t::pop();
        return true;
    }    

} // namespace mort

