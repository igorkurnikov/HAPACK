#include <sstream>
#include <algorithm>
#include <stdexcept>

#include <common.hpp>
#include <object.hpp>
#include "atomfun.hpp"
#include "parmfun.hpp"

namespace mort
{
    using std::count_if;

    using namespace boost;

    int get_atom_pie( const morf_t& a )
    {
        assert( a.cmpid() == ATOM );
        // an: atomic number
        int an = a.get_i(ELEMENT);
        
        if( an==8 || an==16 )
        {
            return 2;
        }
        
        if( an==7 && a.nbond()==3 )
        {
            return 2;
        }
        
        bonditer_t bi = std::find_if( a.bond_begin(), a.bond_end(), 
                                    iparm_cmper1(ORDER, 2) );
 
        if( bi == a.bond_end() )
        {
	    return 0;
	}

        morf_t nbr = atom_oth(*bi,a);
        
        if( nbr.get_i(ELEMENT) == 8 || nbr.get_i(ELEMENT) == 16 )
        {
            return 0;
        }
        
        if( nbr.get_i(ELEMENT) ==7 && an==6 )
        {
            return 0;
        }

        if( nbr.get_i(ELEMENT)==6 && an==7 )
        {
            return 2;
        }
    
        return 1;
    }


    int get_atom_hybrid( const morf_t& a )
    {
        int hybrid = 0;
        if( a.get_i(HYBRID,hybrid) && hybrid != 0 )
        {
            return hybrid;
        }

        int ndouble = count_double_bond(a);
        int ntriple = count_triple_bond(a);

        if( ntriple==1 || ndouble==2 )
        {
            return SP1;
        }
        
        if( ndouble==1 )
        {
            return SP2;
        }
        
        int an = a.get_i(ELEMENT);
        
        if( an==NITROGEN || an==OXYGEN || an==SULFUR )
        {
            if( count_if( a.atom_begin(), a.atom_end(), count_double_bond ) )
            {
                return SP2;
            }            
        }
        
        return SP3;
    }

    int get_atom_degree( const morf_t& a )
    {
        return a.natom() - get_atom_hydrogen( a );
    }

    int get_atom_valence( const morf_t& a )
    {
        return sum( a.bond_begin(), a.bond_end(), 0, 
                    iparm_getter(ORDER) );
    }
    
    int get_atom_hydrogen( const morf_t& a )
    {
        return count_if( a.atom_begin(), a.atom_end(), 
                         iparm_cmper1(ELEMENT,HYDROGEN) );
    }

    int count_bond(const morf_t& a, int order)
    {
        assert( a.cmpid() == ATOM );
        return count_if( a.bond_begin(), a.bond_end(),
                         iparm_cmper1(ORDER,order) );
    }

    int count_single_bond( const morf_t& atom )
    {
        return count_bond( atom, SINGLE_BOND );
    }    

    int count_double_bond( const morf_t& atom )
    {
        return count_bond( atom, DOUBLE_BOND );
    }    

    int count_triple_bond( const morf_t& atom )
    {
        return count_bond( atom, TRIPLE_BOND );
    }

    bool is_sp1( const morf_t& atom )
    {
        return get_atom_hybrid( atom )==SP1;
    }

    bool is_sp2( const morf_t& atom )
    {
        return get_atom_hybrid( atom )==SP2;
    }

    bool is_sp3( const morf_t& atom )
    {
        return get_atom_hybrid( atom )==SP3;
    }

    bool is_amide_carbon( const morf_t& a )
    {
        if( a.get_i(ELEMENT) != CARBON )
	    return false;

	if( a.natom() != 3 )
	    return false;

        bool alone_oxygen= false;
	bool nitrogen = false;

	atomiter_t ai = a.atom_begin();
	for( ; ai != a.atom_end(); ++ai )
	{
	    if( ai->get_i(ELEMENT)==NITROGEN )
	        nitrogen = true;

	    if( ai->get_i(ELEMENT)==OXYGEN && ai->natom()==1 )
	        alone_oxygen = true;
	}

        return alone_oxygen && nitrogen;
    }


    morf_t atom_1st(const morf_t& o)
    {
        if( o.natom()==0 )
	{
	    throw std::runtime_error( "Error: " + unhash(o.cmpid()) + " has zero related atom." );
	}

        return *o.atom_begin();
    }

    morf_t atom_2nd(const morf_t& o)
    {
        assert( o.natom() > 1 );
        return *(o.atom_begin()+1);
    }

    morf_t atom_3rd(const morf_t& o)
    {
        assert( o.natom() > 2 );
        return *(o.atom_begin()+2);
    }

    morf_t atom_4th(const morf_t& o)
    {
        assert( o.natom() > 3 );
        return *(o.atom_begin()+3);
    }

    morf_t atom_oth(const morf_t& b, const morf_t& a)
    {
        assert( b.cmpid()==BOND && a.cmpid()==ATOM );

	atomiter_t first = b.atom_begin();

	return ( *first == a ) ? *(first+1): *first;
    }
 
    string uniq_name( const morf_t& a )
    {
        assert( a.cmpid()==ATOM );

        std::ostringstream os;
        os << a.get_s(SYMBOL) << a.get_i(ID);
        return os.str();
    }

    // section for create_atom
    bool is_last_resd( morf_t& r )
    {
        assert( r.cmpid()==RESD );    
        return r.relid()==r.getmol().nresd()-1;
    }

    // for last residue: insert point is natom, 
    // otherwise it is the relative ID of the first atom of the next residue.
    int insert_point( morf_t& r )
    {
        if( is_last_resd(r) )
        {
            return r.getmol().natom();
        }

        molecule_t* pmol = &r.getmol();
        int resdpos = r.relid();
        resditer_t rnext = pmol->resd_begin() + resdpos + 1;
        assert( rnext != pmol->resd_end() );

        return rnext->atom_begin()->relid();
    }
 
/*
    morf_t create_atom(morf_t& r, const string& name)
    {
        molecule_t* pmol = &r.getmol();
        morf_t a = is_last_resd(r) ? pmol->create_atom() : pmol->insert_atom( insert_point(r) );
        a.set_s(NAME, name);

        a.connect( r );
        r.connect( a );
        return a;
    }

    morf_t frcget_atom(morf_t& r, const string& name)
    {
        morf_t a(r.getmol(), ATOM, -1);
	return get_atom(r, name, a) ? a : create_atom( r, name );
    }

    bool get_atom(const morf_t& r, const string& name, morf_t& a)
    {
        assert( r.cmpid()==RESD );
        return r.get_related_atom( sparm_cmper1(NAME, name), a);
    }

    morf_t get_atom(const morf_t& r, const string& name)
    {
        assert( r.cmpid()==RESD );
        return r.get_related_atom( sparm_cmper1(NAME, name) );
    }
*/


} // namespace mort


