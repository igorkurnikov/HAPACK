#include <iostream>
#include <object.hpp>
#include "merge.hpp"
#include "atomfun.hpp"
#include "bondfun.hpp"
#include "parmfun.hpp"

namespace mort
{

    resd_t merge( molecule_t& m, const resd_t& s, int lid, std::map<int,int>& idmap );
 
    // merge two molecules into one. return the newly created residue.
    resd_t merge( molecule_t& m, const molecule_t& s, int lid)
    {
	if( m.nresd()==0 && s.nresd()==1 )
	{
            string tmp;
	    if( s.get_s(NAME,tmp) ) m.set_s( NAME, tmp );
            if( s.get_s(TYPE,tmp) ) m.set_s( TYPE, tmp );
        }

        std::map<int, int> idmap;
        resd_t r( m, -1 );
        for( int i=0; i < s.nresd(); ++i )
        {
            int curlid = (lid==-1) ? lid : lid + i;
            resd_t tmp = merge( m, s.resds()[i], curlid, idmap );
            if( i==0  )
            {
                r = tmp;
            }
        }

	return r;
    }

    resd_t merge( molecule_t& m, const resd_t& s, int lid )
    {
        std::map<int, int> idmap;
        return merge( m, s, lid, idmap );
    }

    resd_t merge( molecule_t& m, const resd_t& s, int lid, std::map<int,int>& idmap )
    {
        assert( s.cmpid() == RESD );

        resd_t r = resd_t::create( m, lid );
        copy_allparms( s, r );

        atomiter_t ai = s.atom_begin();
        for( ; ai != s.atom_end(); ++ai )
        {
            atom_t a2 = r.create_atom( ai->get_s(NAME) );
            copy_allparms( *ai, a2 );
            idmap[ai->absid()] = a2.absid();
        }

        
        ai = s.atom_begin();
        for( ; ai != s.atom_end(); ++ai )
        {
            atomiter_t aj = ai->atom_begin();
            bonditer_t bj = ai->bond_begin();
            for( ; aj != ai->atom_end(); ++aj, ++bj )
            {
                int id1 = ai->absid();
                int id2 = aj->absid();
                assert( idmap.count(id1)>0 );
                if( idmap.count(id2)==0 )
                {
                    continue;
                }

                atom_t a1( m, idmap[id1] );
                atom_t a2( m, idmap[id2] );
                if( bond_t::has(a1,a2) ) 
                {
                    continue;
                }

                bond_t b = bond_t::create( a1, a2 );
                copy_allparms( *bj, b ); 
            }

        }

        int head;
        if( s.get_i(HEAD, head) && head > 0 )
        {
            r.set_i(HEAD, idmap[head-1]+1 );
        }

        int tail;
        if( s.get_i(TAIL, tail) && tail > 0 )
        {
            r.set_i(TAIL, idmap[tail-1]+1 );
        }
        
        return r;
    }
 

    // extract a residue to a single molecule. 
    molecule_t unmerge( const resd_t& r )
    {
        molecule_t m;
        merge( m, r );
 
        string name;
        if( r.get_s(NAME, name) )
        {
            m.set_s(NAME, name);
        }

        return m;
    }


}
