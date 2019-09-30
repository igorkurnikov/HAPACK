#include <set>
#include <map>
#include <object.hpp>

#ifndef __GNUC__
#include <string.h>
#endif

namespace mort
{
    using namespace boost;
    using std::for_each;

    void set_chain( morf_t& atom )
    {
        static const char marks[] = { "MESB3456X" };

        int level = atom.get_i(LEVEL);
        
        if( level > int(strlen(marks)-1) && level != 999 )
        {
            level = int( strlen( marks ) - 1 );
        }
        
        string str;

        if( level != 999 )
        {
            atom.set_s(TREE, str.append(1, marks[level]) );
        }
        else
        {
            atom.set_s(TREE, "BLA" );
        }
    }

    void visit( morf_t& atom, morf_t& parent, std::map< int, int >& levels, std::set< int >& visited )
    {
        atomiter_t nbr = atom.atom_begin();
        
        for( ; nbr != atom.atom_end(); ++nbr )
        {
            int index = bond_t::get(*nbr, atom).absid();

            if( visited.count( index ) > 0 )
            {
                continue;
            }
            
            if( levels.count( nbr->absid() ) == 0 && *nbr->resd_begin() == *atom.resd_begin() )
            {
                levels[ nbr->absid() ] = nbr->natom();
            }
            else if( *nbr != parent )
            {
                visited.insert( bond_t::get(*nbr, atom).absid() );
                levels[nbr->absid() ]--;
                levels[atom.absid() ]--;
            }
        }
        
        nbr = atom.atom_begin();
        
        for( ; nbr != atom.atom_end(); ++nbr )
        {
            int index = bond_t::get(*nbr, atom).absid();
          
            if( *nbr != parent && visited.count( index ) == 0 )
            {
                visit(*nbr, atom, levels, visited);
            }

            visited.insert( index );
        }
    }
    
    void mark_resd( resd_t& resd )
    {
        bool skip = false;
       
        int aapos;
        if( resd.get_i(AAPOS,aapos) && aapos != (int)OTHER )
        {
            skip = true;
        }

        int head;
        if( resd.get_i(HEAD, head) && head==0 )
        {
            skip = true;
        }

        int tail;
        if( resd.get_i(TAIL, tail) && tail==0 )
        {
            skip = true;
        }

        atom_t conn_1(resd.getmol(), -1);
        if( ! atom_t::get( resd, "N", conn_1) )
        {
            skip = true;
        }

        atom_t conn_2(resd.getmol(), -1);
        if( ! atom_t::get( resd, "C", conn_2) )
        {
            skip = true;
        }

        if( !skip )
        {
            atomvec_t path = find_path(conn_1, conn_2);

            std::set< int > bonds;
            std::map< int, int > levels;

            visit( conn_1, conn_1, levels, bonds );
            
            for( int i=0; i < (int)path.size(); ++i )
            {
                levels[ path.at( i ).absid() ] = 0;
            }

            atomiter_t ai = resd.atom_begin();
            
            for( ; ai != resd.atom_end(); ++ai )
            {
                if( levels.count( ai->absid() ) > 0 )
                {
                    ai->set_i(LEVEL, levels[ ai->absid() ]);
                }
                else
                {
                    ai->set_i(LEVEL, 999);
                }
            }
        }
        else
        {
            atomiter_t ai = resd.atom_begin();
            atomiter_t ae = resd.atom_end();
            for( ; ai != ae; ++ai )
            {
                ai->set_i(LEVEL, 999);
            }
        }
        
        for_each( resd.atom_begin(), resd.atom_end(), set_chain );
    }

    void mark_chain( molecule_t& mol )
    {
        if( mol.nresd()==0 )
        {
            for_each( mol.atom_begin(), mol.atom_end(), sparm_setter1( TREE, "BLA" ) );
        }
        else
        {
            for_each( mol.resd_begin(), mol.resd_end(), std::ptr_fun(mark_resd) );
        }
    }
        
    
} // namespace mort


