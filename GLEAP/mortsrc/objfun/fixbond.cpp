#include <fstream>
#include <sstream>
#include <stdexcept>

#include <common.hpp>
#include <object.hpp>
#include "atomfun.hpp"
#include "bondfun.hpp"
#include "parmfun.hpp"
#include "geomfun.hpp"
#include "fixbond.hpp"

namespace mort
{
    using std::runtime_error;

    ostream& flog()
    {
        static std::ofstream os( "fixbo.log" );
	return os;
    }

    void fixbond( molecule_t& m )
    {
        vector< isomer_set > isomers;
	fixbo( m, isomers );
    }

    void fixbo( molecule_t& m, vector< isomer_set >& isomers )
    {
        flog() << "fixbo for molecule " << m.get_s(NAME) << std::endl;

        std::set<int> fixed;

	flog() << "    Step 1. fixbo using hard rules" << std::endl;
        atomiter_t ai = m.atom_begin();
        for( ; ai != m.atom_end(); ++ai )
        {                
            fixbo_hard( *ai, fixed ); 
        }

	flog() << "    Step 1. fixbo using bond length rules" << std::endl;
        ai = m.atom_begin();
        for( ; ai != m.atom_end(); ++ai )
        {                
            fixbo_bylen( *ai, fixed ); 
        }

	fixbo_conj( m, fixed, isomers );

    }

    void fixbo_hard( morf_t& a, set<int>& fixed )
    {
        int nnbr = a.natom();
	if( nnbr== nfixed_nbr(a, fixed) )
	{
	    return;
        }

        flog() << "        fixing atom " << a.get_s(NAME) <<  " by hard rules, nnbr is " << nnbr << std::endl;
        int elem = a.get_i(ELEMENT);
        if( ! fixable(elem) )
        {
            throw runtime_error( "Error: can not handle element like this" );
        }

        if( nnbr==1 )
        {
	    morf_t nbr = a.atoms()[0];
            morf_t bnd = a.bonds()[0];

            if( elem == HYDROGEN || is_halogen(elem) )
            {
                bnd.set_i( ORDER, 1 );
                fixed.insert( bnd.absid() );
		flog() << "            all bonds (only one) set to single for hydrogen or halgen" << std::endl;
                return;
            }

	    if( elem==SULFUR && nbr.get_i(ELEMENT)==PHOSPHORUS )
	    {
	        bnd.set_i( ORDER, 2 );
		fixed.insert( bnd.absid() );
		flog() << "            P-S bond set to double when S is alone" << std::endl;
                return;
            }


	    flog() << "            can not be fixed at this step." << std::endl;
	    return;
        }

	if( nnbr==2 )
	{
	    fixbo_hard_nnbr2( a, fixed );
	    return;
	}

	if( nnbr==3 )
	{
	    fixbo_hard_nnbr3( a, fixed );
	    return;
	}

        if( nnbr==4 )
	{
	    if( elem==CARBON || elem==NITROGEN ) 
            {
                setbo_allone( a, fixed );
                flog() << "            all bonds set to single for carbon and nitrogen and nnbr equal 4" << std::endl;
            }
            else
            {
                assert(elem==PHOSPHORUS || elem==SULFUR || elem==CHLORINE || elem==BROMINE || elem==IODINE);
                fixbo_acid( a, fixed );
            }
            return;
        }
    }

    int order_from_angl( const morf_t& atml, const morf_t& atmc, const morf_t& atmr, const set<int>& fixed )
    {
	bond_t bnd0 = bond_t::get( atml, atmc );
	bond_t bnd1 = bond_t::get( atmc, atmr );

        if( ndouble_bnd(atmc, fixed) > 0 )
	{ 
	    return 0;
	}

	if( ndouble_bnd(atml, fixed) > 0 )
	{
	    return 4;
	}

	if( ndouble_bnd(atmr, fixed) > 0 )
	{
	    return 5;
	}

        int order0 = bnd0.get_i(ORDER);
	int order1 = bnd1.get_i(ORDER);

	if( order0==5 && order1==6 )
	{
            return 2;
	}

	if( order0==6 && order1==5 )
	{
	    return 1;
	}

	if( order0==4 && order1 >4 )
	{
	    return order1==6 ? 2 : 1;
	}

	if( order0 >4 && order1==4 )
	{
	    return order0==5 ? 2 : 1;
	}


        int nnbrl = atml.natom();
	int nnbrc = atmc.natom();
	int nnbrr = atmr.natom();

	if( nnbrl > 1 )
	{
	    return on_plane(bnd0) ? 2 : 4;
	}

	if( nnbrr > 1 )
	{
	    return on_plane(bnd1) ? 1 : 5;
	}

        assert( nnbrl==1 && nnbrr==1 );

	if( nnbrc==3 )
	{
	    return on_plane(atmc) ? 2 : 0;
	}

        return 2;
    }

    void fixbo_angl( vector< morf_t >& bnds, set<int>& fixed )
    {
        assert( bnds.size()==2 );

	vector< morf_t > atms;

	morf_t atm1 = atom_1st( bnds[0] );
        morf_t atm2 = atom_2nd( bnds[0] );
	morf_t atm3 = atom_1st( bnds[1] );
	morf_t atm4 = atom_2nd( bnds[1] );

	if( atm1==atm3 )
	{
	    atms.push_back( atm2 );
	    atms.push_back( atm1 );
	    atms.push_back( atm4 );
	}
	else if( atm1==atm4 )
	{
	    atms.push_back( atm2 );
            atms.push_back( atm1 );
	    atms.push_back( atm3 );
	}
	else if( atm2==atm3 )
	{
	    atms.push_back( atm1 );
	    atms.push_back( atm2 );
	    atms.push_back( atm4 );
	}
	else
	{
	    assert( atm2==atm4 );
	    atms.push_back( atm1 );
	    atms.push_back( atm2 );
	    atms.push_back( atm3 );
	}

        int code = order_from_angl( atms[0], atms[1], atms[2], fixed );
	if( code < 3 )
	{
	    int order1 = 1 + code/2;
	    int order2 = 1 + code%2;

	    bnds[0].set_i( ORDER, order1 );
	    bnds[1].set_i( ORDER, order2 );
            fixed.insert( bnds[0].absid() );
	    fixed.insert( bnds[1].absid() );
	    return;
	}
	
	assert( code==4 || code==5 );
	int id = (code==4) ? 0 : 1;

	bnds[id].set_i( ORDER, 1 );
	fixed.insert( bnds[id].absid() );
	fixbo_isolated( bnds[1-id], fixed );
    }


    void fixbo_conj( molecule_t& m, set<int>& fixed, vector< isomer_set >& isomers )
    {
        flog() << "    fixing conjugated system" << std::endl;
        vector< morf_t > unfixed;
	bonditer_t bi = m.bond_begin();
	for( ; bi != m.bond_end(); ++bi )
	{
            if( bi->get_i(ORDER) <4 )
	        continue;

	    assert( fixed.count( bi->absid() )==0 );

            morf_t atm1 = atom_1st( *bi );
	    morf_t atm2 = atom_2nd( *bi );

	    int nnbr1 = atm1.natom();
	    int nnbr2 = atm2.natom();

	    if( bi->get_i(ORDER)==6 && nnbr1 >=2 && nnbr2 >=2 && !on_plane(*bi) )
	    {
	        flog() << "        bond " << atm1.get_s(NAME) << "-" << atm2.get_s(NAME) << " set to single since now all related atoms are on plane." << std::endl;
		bi->set_i( ORDER, 1 );
		fixed.insert( bi->absid() );
            }
            else
	    { 
	        unfixed.push_back( *bi );
	    }
	}

        vector< vector<morf_t> > atmgrps;
	vector< vector<morf_t> > bndgrps;
        divide_bonds( unfixed, atmgrps, bndgrps );

        flog() << "    totally " << unfixed.size() << " bonds unfixed, divided into " << atmgrps.size() << " groups." << std::endl;

        assert( atmgrps.size()==bndgrps.size() );
	for( int i=0; i < (int)atmgrps.size(); ++i )
	{
	    flog() << "    group " << i << " has " << bndgrps[i].size() << " bonds." << std::endl;

	    if( bndgrps[i].size()==1 )
	    {
	        fixbo_isolated( bndgrps[i][0], fixed );
		flog() << "        bond " << atmgrps[i][0].get_s(NAME) << "-" << atmgrps[i][1].get_s(NAME) << " set to " << bndgrps[i][0].get_i(ORDER) << std::endl;
		continue;
	    }

            if( bndgrps[i].size()==2 )
	    {
	        fixbo_angl( bndgrps[i], fixed );
		continue;
	    }

	    assert( atmgrps[i].size() > 3 );

/*
	    if( on_plane(atmgrps[i]) )
	    {
 */
	        isomer_set isomer( bndgrps[i] );
                fixbo_isomer( atmgrps[i], bndgrps[i], fixed, isomer );
                if( isomer.size() > 1 )
	        {
	            isomers.push_back( isomer );
	        }
/*
            }
            else
	    {
	        fixbo_other( atmgrps[i], bndgrps[i], fixed, isomers );
	    }
 */
	}
 
    }

    
	    





    morf_t other_nbr( const morf_t& a, const morf_t& n )
    {
        assert( a.natom()==2 );

	if( a.atoms()[0] == n )
	    return a.atoms()[1];
	
	return a.atoms()[0];
    }

    double longer_diff( const morf_t& a, const morf_t& bnd )
    {
        double d = dist( atom_1st(bnd), atom_2nd(bnd) );
	double maxd = 0.0;

	for( int i=0; i < a.natom(); ++i )
	{
	    morf_t nbr = a.atoms()[i];
	    double d2 = dist( a, nbr );
	    if( d2 > maxd )
	        maxd = d2 ;
	    flog() << "    " << nbr.get_s(NAME) << "    " << d2 << std::endl;
	}

	return maxd - d;
    }

    void fixbo_isolated( morf_t& bnd, set<int>& fixed )
    {
        morf_t atm1 = atom_1st( bnd );
	morf_t atm2 = atom_2nd( bnd );

        flog() << "        fixing isolated bond " << atm1.get_s(NAME) << "-" << atm2.get_s(NAME) << std::endl;
	int ndouble1 = ndouble_bnd( atm1, fixed );
	int ndouble2 = ndouble_bnd( atm2, fixed );

	if( ndouble1==0 && ndouble2==0 )
	{
	    if( on_plane(bnd) )
	    {
	        flog() << "       bond " << atm1.get_s(NAME) << "-" << atm2.get_s(NAME);
	        flog() << " set to double since all related atoms are on same plane." << std::endl;
	        bnd.set_i( ORDER, 2 );
		fixed.insert( bnd.absid() );
		return;
	    }


            if( bnd.get_i(ORDER)==5 )
	    {
	    /*
                flog() << " not all related atom on same plane." << std::endl;
                double diff1 = longer_diff( atm1, bnd );
	        double diff2 = longer_diff( atm2, bnd );

                flog() << " diff of atm " << atm1.get_s(NAME) << " " << diff1 << std::endl;
	        flog() << " diff of atm " << atm2.get_s(NAME) << " " << diff2 << std::endl;

	        if( diff1 > 0.1 && diff2 > 0.1 )
	        {
	            flog() << " set to double since length difference is big." << std::endl;
	     */
	            bnd.set_i( ORDER, 2 );
		    fixed.insert( bnd.absid() );
		    return;
	//	}
            }

	    flog() << "       bond " << atm1.get_s(NAME) << "-" << atm2.get_s(NAME);
	    flog() << " set to single since not all related atoms are on same plane." << std::endl;
	    bnd.set_i( ORDER, 1 );
            fixed.insert( bnd.absid() );
            return;
        }

	assert( ndouble1==1 || ndouble2==1 );
	if( ndouble1 != ndouble2 )
	{
	    flog() << "    Warning: bond " << atm1.get_s(NAME) << "-" << atm2.get_s(NAME);
	    flog() << ", one has double bond, one has not, can not form a conjugate system.";
	    flog() << " set to single anyway." << std::endl;
        }
        else
	{
	    flog() << "    bond " << atm1.get_s(NAME) << "-" << atm2.get_s(NAME);
	    flog() << " set to single" << std::endl;
	}

        bnd.set_i( ORDER, 1 );
        fixed.insert( bnd.absid() );
    }

    int order_of_conj( const morf_t& bnd, const morf_t& atm, set<int>& fixed, map<int, int>& solu )
    {
        assert( fixed.count(bnd.absid())==0 && solu.count(bnd.absid())==0 );

        int nnbr = atm.natom();
	int isingle = 0;
        int idouble = 0;
        for( int i=0; i < nnbr; ++i )
	{
	    morf_t b = atm.bonds()[i];
	    if( fixed.count(b.absid()) > 0 )
	    {
	        int order = b.get_i(ORDER);
		if( order==1 )
		    isingle++;
		else
		    idouble++;
	    }

	    if( solu.count(b.absid()) > 0 )
	    {
		int order = solu[b.absid()];
                if( order==1 )
		    isingle++;
		else
		    idouble++;
	    }
	}

        if( idouble > 0 )
	    return 1;

	if( nnbr==3 && isingle==2 && bnd.get_i(ORDER)!=6)
	{
	    assert( atm.get_i(ELEMENT)==CARBON );
            return 2;
	}

	if( nnbr==2 && isingle==1 && atm.get_i(ELEMENT)==CARBON && bnd.get_i(ORDER)!=6 )
	{
	    return 2;
	}

	return 4;
    }

    int fixbo_isomer_each( vector<morf_t>& atms, vector<morf_t>& bnds, set<int>& fixed, map<int, int>& solu )
    {
        int nfixed = solu.size();

	do
	{
	    nfixed = solu.size();
            for( int i=0; i < (int)bnds.size(); ++i )
	    {
	        if( solu.count(bnds[i].absid()) > 0 )
	            continue;

                morf_t atm1 = atom_1st( bnds[i] );
		morf_t atm2 = atom_2nd( bnds[i] );

                int order1 = order_of_conj( bnds[i], atm1, fixed, solu );
		int order2 = order_of_conj( bnds[i], atm2, fixed, solu );

		if( order1==4 && order2==4 )
		    continue;

		int bid = bnds[i].absid();

		if( order1==4 )
		{
		    assert( order2==1 || order2==2 );
		    solu[bid] = order2;
		    flog() << "        try setting bond " << atm1.get_s(NAME) << "-" << atm2.get_s(NAME) << " to " << order2 << " according to conj rule." << std::endl;
		    continue;
                }

		if( order2==4 )
		{
		    assert( order1==1 || order1==2 );
		    solu[bid] = order1;
		    flog() << "        try setting bond " << atm1.get_s(NAME) << "-" << atm2.get_s(NAME) << " to " << order1 << " according to conj rule." << std::endl;
		    continue;
		}

		if( order1 != order2 )
		{
                    flog() << "        confliction encountered while trying setting bond " << atm1.get_s(NAME) << "-" << atm2.get_s(NAME) << std::endl;
		    return WRONG_DIRECTION;
		}

		flog() << "        try setting bond " << atm1.get_s(NAME) << "-" << atm2.get_s(NAME) << " to " << order1 << " according to conj rule." << std::endl;
		solu[bid] = order1;
	    }    
        }
	while( (int)solu.size() > nfixed );

        return solu.size()==bnds.size() ? ENOUGH : KEEP_GOING;    
    }

    int fixbo_isomer_conj( vector<morf_t>& atms, vector<morf_t>& bnds, set<int>& fixed, map<int, int>& solu, isomer_set& solus )
    {
        assert( solu.size() < bnds.size() );

        vector< morf_t > atms_left;
        vector< morf_t > bnds_left;

        for( int i=0; i < (int)bnds.size(); ++i )
        {
            if( solu.count(bnds[i].absid())>0 )
                continue;

            bnds_left.push_back( bnds[i] );
		 
            morf_t atm1 = atom_1st( bnds[i] );
            morf_t atm2 = atom_2nd( bnds[i] );
	    if( std::find( atms_left.begin(), atms_left.end(), atm1)== atms_left.end() )
	        atms_left.push_back( atm1 );

	    if( std::find( atms_left.begin(), atms_left.end(), atm2)== atms_left.end() )
	        atms_left.push_back( atm2 );
        }

        assert( atms_left.size() > 0 );

        vector< map<int,int> > attempts;
        for( int i=0; i < (int)bnds_left.size(); ++i )
        {
            morf_t bnd = bnds_left[i];
	    morf_t atm1 = atom_1st( bnd );
	    morf_t atm2 = atom_2nd( bnd );

            map<int, int> attempt( solu );
            attempt[ bnd.absid() ] = 2;
             
	    flog() << "        try setting bond " << atm1.get_s(NAME) << "-" << atm2.get_s(NAME) << " to double, start a journey" << std::endl;
            int result = fixbo_isomer_each( atms, bnds, fixed, attempt );
	    
            if( result == WRONG_DIRECTION )
	    {
	        flog() << "        attemption failed" << std::endl << std::endl;
                continue;
            }

            if( result == ENOUGH )
	    {
	        flog() << "        attemption is enough to fix the whole system." << std::endl;
		assert( attempt.size()==bnds.size() );
	        solus.add_solution( attempt );
		continue;
	    }

            assert( result==KEEP_GOING );

	    flog() << "        attemption is success, but not enough to fix the whole system, need another wild guess." << std::endl;
            attempts.push_back( attempt );
        }

	for( int i=0; i < (int)attempts.size(); ++i )
	{
	    flog() << "        resuming No. " << i << " attempts." << std::endl;
	    fixbo_isomer_conj( atms, bnds, fixed, attempts[i], solus );
        }

	return solus.size()==0? WRONG_DIRECTION : ENOUGH;
    }     

    void fixbo_isomer( vector<morf_t>& atms, vector<morf_t>& bnds, set<int>& fixed, isomer_set& solus )
    {
        map<int,int>  hard_solu;
        int result1 = fixbo_isomer_each( atms, bnds, fixed, hard_solu );
        if( result1 == ENOUGH )
	{
	    solus.add_solution( hard_solu );
	    solus.apply_best();
	    return;
	}

	if( result1 == WRONG_DIRECTION )
	{
	    string errmsg = "Error: can not fix conjugate system, initial attempt failed";
	    flog() << "        " << errmsg << std::endl;
	    throw runtime_error( errmsg );
	}

        fixbo_isomer_conj( atms, bnds, fixed, hard_solu, solus );
	if( solus.size() > 0 )
        {
	    solus.apply_best();
	}
        else
	{
	    flog() << "            can not fix isomer groups" << std::endl;
	}
    }


    int find_atom( const vector< vector<morf_t> >& atmgrps, const morf_t& atm )
    {
        for( int i=0; i < (int)atmgrps.size(); ++i )
        {
	    if( std::find(atmgrps[i].begin(), atmgrps[i].end(), atm) != atmgrps[i].end() )
	            return i;
	}

	return atmgrps.size();
    }

    void divide_bonds( const vector<morf_t>& bnds, vector< vector<morf_t> >& atmgrps, vector< vector<morf_t> >& bndgrps )
    {
        for( int i=0; i < (int)bnds.size(); ++i )
        {
            morf_t atm1 = atom_1st( bnds[i] );
	    morf_t atm2 = atom_2nd( bnds[i] );

            int id_1 = find_atom( atmgrps, atm1 );
            int id_2 = find_atom( atmgrps, atm2 );

            if( id_1 == (int)atmgrps.size() && id_2 == (int)atmgrps.size() )
	    {
	        vector< morf_t > tmpatms;
                vector< morf_t > tmpbnds;
	        tmpatms.push_back( atm1 );
                tmpatms.push_back( atm2 );
                tmpbnds.push_back( bnds[i] );
                atmgrps.push_back( tmpatms );
	        bndgrps.push_back( tmpbnds );
	        continue;
            }

            if( id_1 == (int)atmgrps.size() && id_2 != (int)atmgrps.size() )
            {
                atmgrps[id_2].push_back( atm1 );
	        bndgrps[id_2].push_back( bnds[i] );
	        continue;
	    }

	    if( id_1 != (int)atmgrps.size() && id_2 == (int)atmgrps.size() )
	    {
	        atmgrps[id_1].push_back( atm2 );
	        bndgrps[id_1].push_back( bnds[i] );
	        continue;
            }

            assert( id_1 < (int)atmgrps.size() && id_2 < (int)atmgrps.size() );

            if( id_1 == id_2 )
	    {
	        bndgrps[id_1].push_back( bnds[i] );
                assert( atmgrps[id_1].size() <= bndgrps[id_1].size() );
	        continue;
            }

            int minid = std::min( id_1, id_2 );
            int maxid = std::max( id_1, id_2 );

	    atmgrps[minid].insert( atmgrps[minid].end(), atmgrps[maxid].begin(), atmgrps[maxid].end() );
	    bndgrps[minid].insert( bndgrps[minid].end(), bndgrps[maxid].begin(), bndgrps[maxid].end() );
            bndgrps[minid].push_back( bnds[i] );

	    atmgrps.erase( atmgrps.begin()+maxid );
	    bndgrps.erase( bndgrps.begin()+maxid );
        }
    }




    void fixbo_hard_nnbr2( morf_t& a, set<int>& fixed )
    {
        assert( a.natom()==2 );

        int elem = a.get_i(ELEMENT);
        if( elem==OXYGEN || elem==SULFUR ) 
        {
	    a.set_i( HYBRID, SP3 );
            setbo_allone( a, fixed );
            flog() << "            all bonds set to single for oxygen and nnbr equals 2" << std::endl;
            return;
        }

        assert( elem==CARBON || elem==NITROGEN || elem==PHOSPHORUS );

        morf_t nbr0 = a.atoms()[0];
        morf_t nbr1 = a.atoms()[1];
        morf_t bnd0 = a.bonds()[0];
        morf_t bnd1 = a.bonds()[1];
        if( nbr0.natom() < nbr1.natom() )
        {
            std::swap( nbr0, nbr1 );
            std::swap( bnd0, bnd1 );
        }

        int nnbr0 = nbr0.natom();
     
        if( nnbr0 == 1 )
	{
	    // two nbrs are both alone. can only determine by later;
	    return;
	}

        double ang = angl( nbr0, a, nbr1 );
        flog() << "            angle " << nbr0.get_s(NAME) << "-" << a.get_s(NAME) << "-" << nbr1.get_s(NAME) << ": " << ang << std::endl;
	if( ang > 175.0 && ang < 185.0 )
	{
            // sp1 system, here we assume no continuous triple bond systems like C#C-C#C, which might not be true
	    // should be fixed in the future.
	    //
	    assert( elem == CARBON && nnbr0 > 1 );

	    flog() << "            So this is a sp1 atom, test if nbr " << nbr0.get_s(NAME) << " in line? ";
            if( nnbr0==2 && in_line(nbr0) )
            {
	        flog() << " yes! make triple bond with it " << std::endl;
                bnd0.set_i( ORDER, 3 );
                bnd1.set_i( ORDER, 1 );
            }
            else
            {
                int elem1 = nbr1.get_i(ELEMENT);
		if( elem1==CARBON || elem1==NITROGEN )
                {
                    flog() << " no!  make single bond with it, make triple with the other one " << std::endl;
                    bnd0.set_i( ORDER, 1 );
                    bnd1.set_i( ORDER, 3 );
		}
		else
		{
                    flog() << " no! but the other is oxygen or sulfur, make two double bonds " << std::endl;
		    bnd0.set_i( ORDER, 2 );
		    bnd1.set_i( ORDER, 2 );
		}
            }

            flog() << "            bond order set for SP1 system" << std::endl;
            fixed.insert( bnd0.absid() );
            fixed.insert( bnd1.absid() );
            return;
        }

	flog() << "            Can not be fixed at this step. leave it later." << std::endl;

    }

    void fixbo_hard_nnbr3( morf_t& a, set<int>& fixed )
    {
        int elem = a.get_i( ELEMENT );
        if( elem==NITROGEN || elem==PHOSPHORUS )
        {
	    flog() << "            N or P";
            if( is_acid(a) )
            {
	        flog() << " is acid, set order accordingly" << std::endl;
                fixbo_acid( a, fixed ); 
            }
            else
            {
                flog() << " is not acid, set all bonds to single." << std::endl;
                setbo_allone(a, fixed );
            }
            
            return;    
        }
        
        if( elem==SULFUR || elem==CHLORINE || elem==BROMINE || elem==IODINE ) 
        {
            if( is_acid(a) )
	    {
                fixbo_acid( a, fixed );
	    }
	    else
	    {
	        flog() << "Warning: can not understand such group: (S,Cl,Br,I) nnbr==3 but not acid" << std::endl;
	    }

            return;
        }

	flog() << "            can not be fixed at this time, leave it later." << std::endl;
    }

    void fixbo_bylen( morf_t& a, set<int>& fixed )
    {
        flog() << "        fixing atom " << a.get_s(NAME) << " by length rules." << std::endl;
        int idouble = 0;
        int nnbr = a.natom();
        for( int i=0; i < nnbr; ++i )
	{
	    morf_t atm = a.atoms()[i];
	    morf_t bnd = a.bonds()[i];
	    int order = order_from_length( bnd );
	
            if( fixed.count(bnd.absid())>0 )
	    {
	        if( bnd.get_i(ORDER) !=order && order <4 )
                {
    	            flog() << "Warning: confliction found for bond " << a.get_s(NAME) << "-" << atm.get_s(NAME) << std::endl;
		    flog() << "was set to " << bnd.get_i(ORDER) << " before, but should be " << order << " according to its length." << std::endl;
		}

		if( bnd.get_i(ORDER)==2 ) idouble++;
		continue;
            }

	    if( order >=4 )
	    {
		flog() << "            bond " << a.get_s(NAME) << "-" << atm.get_s(NAME) << "(len " << dist(a, atm) << ")" << " set to " << order << " , means not sure,";
		if( order==5 )
		    flog() << " has better chance to be double";
		if( order==6 )
		    flog() << " has better chance to be single";
		flog() << std::endl;
                bnd.set_i(ORDER, order );
	    }
	    else
	    {
                if( order==2 && nnbr==3 && !on_plane(a) )
		{
		    flog() << "Warning: confliction found for bond " << a.get_s(NAME) << "-" << atm.get_s(NAME) << std::endl;
		    flog() << "should be double bond according to the length, but not all nbrs of " << a.get_s(NAME) << " are on the same plane" << std::endl;
                }

		if( order==2 ) idouble++;
                
		flog() << "            bond " << a.get_s(NAME) << "-" << atm.get_s(NAME) << " set to " << order << " for its length " << dist(a, atm) << std::endl;
		bnd.set_i(ORDER, order);
		fixed.insert( bnd.absid() );
	    }
	}

	if( idouble <= 1)
	    return;

        if( a.get_i(ELEMENT)==SULFUR && nnbr==4 )
	    return;

        if( a.get_i(ELEMENT)==NITROGEN && nnbr==3 )
	    return;

	if( a.get_i(ELEMENT)==CARBON && nnbr==2 && in_line(a) )
	    return;

        throw runtime_error( "Error: more than one double bond made to atom " + a.get_s(NAME) +
                      ", there must be something wrong, are you sure your structure is reasonable?" );

    }

    void setbo_allone( morf_t& a, std::set<int>& fixed )
    {
        bonditer_t bi = a.bond_begin();
        for( ; bi != a.bond_end(); ++bi )
        {
            bi->set_i(ORDER, 1);
            fixed.insert( bi->absid() );
        }

        int elem = a.get_i(ELEMENT);
	int nnbr = a.natom();

        if( elem== NITROGEN && nnbr==4 )
        {
            a.set_i(FCHG, 1);
        }
            
        return;
    }

    bool fixable(int elem) 
    {
        static const int ELEMS[] = { HYDROGEN, CARBON, NITROGEN, OXYGEN, FLUORINE, PHOSPHORUS, SULFUR, CHLORINE, BROMINE, IODINE};
        static const int NELEM   = sizeof(ELEMS) / sizeof(int);    
        
        return std::find(ELEMS, ELEMS+NELEM, elem) != ELEMS+NELEM;
    }

    bool is_halogen(int elem)
    {
        return elem==FLUORINE || elem==CHLORINE || elem ==BROMINE || elem==IODINE;
    }

    bool in_line( const morf_t& a )
    {
        assert( a.natom()==2 );

        vector<numvec> pts( 3, numvec(3) );
        pts[0] = a.get_v(POSITION);
        pts[1] = a.atoms()[0].get_v(POSITION);
        pts[2] = a.atoms()[1].get_v(POSITION);
        return same_line( pts );
    }

    int nfixed_nbr( const morf_t& a, const set<int>& fixed )
    {
        int ifixed = 0;

        bonditer_t bi = a.bond_begin();
        for( ; bi != a.bond_end(); ++bi )
        {
            if( fixed.count(bi->absid())>0 )
            {
                ifixed++;
            }
        }

        return ifixed;
    }

    int ndouble_bnd( const morf_t& a, const set<int>& fixed )
    {
        int idouble = 0;

        bonditer_t bi = a.bond_begin();
        for( ; bi != a.bond_end(); ++bi )
        {
            if( fixed.count(bi->absid())>0 && bi->get_i(ORDER)==2 )
            {
                idouble++;
            }
        }

        return idouble;
    }



    void fixbo_acid( morf_t& a, std::set<int>& fixed )
    {
        int elem1 = a.get_i(ELEMENT);
        int nnbr1 = a.natom();
        
        int ndouble = 0;
        if( nnbr1 == 4 )
        {
            assert( elem1 == SULFUR || elem1 == PHOSPHORUS || elem1 == CHLORINE || elem1==BROMINE || elem1==IODINE );
            ndouble = elem1 - 14;
        }
        else if( nnbr1 == 3 )
        {
            if( elem1 == NITROGEN || elem1 == PHOSPHORUS || elem1 == CHLORINE || elem1==BROMINE || elem1==IODINE)
            {
                ndouble = 2;
            }
            else 
            {
                assert( elem1 == SULFUR || elem1==CARBON );
                ndouble = 1;
            }
        }
	else
	{
	    ndouble = 1;
	}

        int idouble=0;
        for( int i=0; i < a.natom(); ++i )
	{
	    morf_t atm = a.atoms()[i];
	    morf_t bnd = a.bonds()[i];
            int elem = atm.get_i(ELEMENT);

	    if( (elem==OXYGEN||elem==SULFUR) && atm.natom()==1 )
	    {
	        idouble++;
	        bnd.set_i( ORDER, 2 );
            }
	    else
	    {
	        bnd.set_i( ORDER, 1 );
	    }

	    fixed.insert( bnd.absid() );
	}

        if( idouble != ndouble )
	{
	    flog() << "Warning: acid group has more or less oxygens then expected" << std::endl;
	}

    }

    bool is_acid( morf_t& a )
    {
        for( int i=0; i < a.natom(); ++i )
	{
	    morf_t atm = a.atoms()[i];
	    if( atm.get_i(ELEMENT)==OXYGEN && atm.natom()==1 )
	        return true;
	}

	return false;
    }

    bool on_plane( const morf_t& o )
    {
        if( o.cmpid()==ATOM )
	{
	    morf_t a = o;
            assert( a.natom() == 3 );

            morf_t nbr0 = a.atoms()[0];
	    morf_t nbr1 = a.atoms()[1];
            morf_t nbr2 = a.atoms()[2];

            numvec v = a.get_v(POSITION);
            numvec crd0 = nbr0.get_v(POSITION);
	    numvec crd1 = nbr1.get_v(POSITION);
	    numvec crd2 = nbr2.get_v(POSITION);

	    numvec p01 = plane( v, crd0, crd1 );
            if( same_plane(crd2, p01) ) return true;

            numvec p12 = plane( v, crd1, crd2 );
	    if( same_plane(crd0, p12) ) return true;
 
	    numvec p20 = plane( v, crd2, crd0 );
            if( same_plane(crd1, p20) ) return true;

	    return false;
	}

	if( o.cmpid()==BOND )
	{
            morf_t b = o;
            morf_t atm1 = atom_1st( b );
	    morf_t atm2 = atom_2nd( b );

            int nnbr1 = atm1.natom();
	    int nnbr2 = atm2.natom();

	    vector<morf_t> atms;
	    for( int i=0; i < nnbr1; ++i )
	    {
	        morf_t tmp = atm1.atoms()[i];
		if( tmp != atm2 )
		    atms.push_back( tmp );
            }

	    atms.push_back( atm1 );
            atms.push_back( atm2 );

	    for( int i=0; i < nnbr2; ++i )
	    {
	        morf_t tmp = atm2.atoms()[i];
		if( tmp != atm1 )
		    atms.push_back( tmp );
	    }

	    if( atms.size() <=3 )
	    {
	        return true;
            }

	    return on_plane(atms);
	}
  
        throw runtime_error( "Error: this on_plane only works on atoms and bonds" );
       
    }

    bool on_plane( const vector<morf_t>& atms )
    {
        assert( atms.size() > 3 );

        for( int i=0; i < (int)atms.size(); ++i )
	{
	    morf_t a0 = atms[i];
	    morf_t a1 = atms[ (i+1)%atms.size() ];
	    morf_t a2 = atms[ (i+2)%atms.size() ];
            morf_t a3 = atms[ (i+3)%atms.size() ];

	    numvec p = plane( a0, a1, a2 );

	    if( !same_plane(a3.get_v(POSITION), p) )
                return false;
	}

	return true;
    }

    numvec plane( const morf_t& a0, const morf_t& a1, const morf_t& a2 )
    {
        return plane( a0.get_v(POSITION), a1.get_v(POSITION), a2.get_v(POSITION) );
    }

    static const char LENTABLE[] = "C C  1.33 1.49 "
                                   "C N  1.30 1.43 "
				   "C O  1.28 1.41 "
				   "C S  1.72 1.80 "
                                   "N N  1.29 1.38 "
                                   "N O  1.24 1.39 ";
    struct lenrule_t
    {
        int elem1;
	int elem2;
	double cutoff_double;
	double cutoff_single;
    };

    int order_from_length( const morf_t& bnd )
    {
        static vector<lenrule_t> rules;

	if( rules.size()==0 )
	{
	    std::istringstream is( LENTABLE );
	    char sym1[10], sym2[10];
	    double cut1, cut2;
	    while( is >> sym1 >> sym2 >> cut1 >> cut2 )
	    {
	        lenrule_t r;
		r.elem1 = pertab_t::get_element( sym1 );
                r.elem2 = pertab_t::get_element( sym2 );
		r.cutoff_double = cut1;
		r.cutoff_single = cut2;
	        rules.push_back( r );
	    }

	}

        morf_t atm1 = atom_1st( bnd );
	morf_t atm2 = atom_2nd( bnd );

        int elem1 = atm1.get_i(ELEMENT);
	int elem2 = atm2.get_i(ELEMENT);

        if( elem1 > elem2 )
	    std::swap( elem1, elem2 );

        if( elem1==HYDROGEN )
	    return 1;

	if( is_halogen(elem2) )
	    return 1;

        double d = dist( atm1, atm2 );
        for( int i=0; i < (int)rules.size(); ++i )
	{
	    if( rules[i].elem1 != elem1 )
	        continue;
	
	    if( rules[i].elem2 != elem2 )
	        continue;

	    if( d < rules[i].cutoff_double )
	        return 2;

	    if( d < rules[i].cutoff_double + 0.05 )
	        return 5;
	
	    if( d < rules[i].cutoff_single - 0.04 )
	        return 4;

	    if( d < rules[i].cutoff_single )
	        return 6;

	    return 1;
	}

        // flog() << "Error: can not find len rule for bond " << atm1.get_s(NAME) << "-" << atm2.get_s(NAME) << std::endl;
	return 4;
    }




} // namespace mort
