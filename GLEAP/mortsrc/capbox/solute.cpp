#include <common.hpp>
#include "solute.hpp"
 
namespace mort
{

    numvec regionlize( solute_i& s, bool vdwr )
    {
        funstack_t::push( "regionlize_solute" );
        numvec extent = vdwr ? regionlize( s.getcord(), s.getvdwr() ) : regionlize( s.getcord() );
        s.finish();
        funstack_t::pop();
        return extent;
    }

    double octbufsize( solute_i& s, double desired )
    {
        return octbufsize( s.getcord(), desired );
    }

    void centralize( solute_i& s )
    {
        centralize( s.getcord() );
        s.finish();
    }

    void rotatewald( solute_i& s )
    {
        rotatewald( s.getcord() );
        s.finish();
    }

    void rotatelong( solute_i& s )
    {
        rotatelong( s.getcord() );
        s.finish();
    }

    int get_solute_nresd( const molecule_t& m )
    {
        int solute_nresd;
        if( m.get_i(SOLUTE_NRESD, solute_nresd) )
        {
            return solute_nresd;
        }

        solute_nresd = m.nresd();
        resditer_t ri = m.resd_end() - 1;
		
		if( m.nresd() < 10 )
		{
			solute_nresd = m.nresd();
			return solute_nresd;
		}

        string svt = ri->get_s(TYPE);

        while( solute_nresd > 1 && ri->get_s(TYPE)=="WAT" || ri->get_s(TYPE) == "HOH") 
//		while( ri->get_s(TYPE) == svt )
        {
            --ri;
            --solute_nresd;
        }

        if( m.nresd() - solute_nresd < 5 )
        {
            solute_nresd = m.nresd();
        }

        return solute_nresd;
    }

    int get_solute_natom( const molecule_t& m )
    {
        int natom = 0;
        int nresd = get_solute_nresd(m);
        for( int i=0; i < nresd; ++i )
        {
            natom += m.resds()[i].natom();
        }

        return natom;
    }

    int is_sequential( const molecule_t& m )
    {
        int curt = -1;
        atomiter_t ai = m.atom_begin();
        atomiter_t ae = m.atom_end();
        for( ; ai != ae; ++ai )
        {
            if( curt==-1) 
            {
                curt = ai->resd_begin()->relid();
                continue;
            }

            int rid = ai->resd_begin()->relid();
            if( rid != curt )
            {
                if( rid==curt+1 )
                {
                    curt++;
                }
                else
                {
                    return false;
                }
            }
        }
        
        return true;
    } 

} // namespace mort

