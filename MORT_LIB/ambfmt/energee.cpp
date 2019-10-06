#include <capbox.hpp>
#include "energee.hpp"

namespace mort
{
    

    void parametrize( molecule_t& mol, const molecule_t& ff, parmset_t& ps );

    int max_resd_size( const molecule_t& m )
    {
        int r = 0;
        resditer_t ri = m.resd_begin();
        for( ; ri != m.resd_end(); ++ri )
        {
            r = std::max( r, ri->natom() );
        }
        return r;
    }

    bool has_H( const morf_t& mo )
    {
        return mo.atom_end()!=std::find_if(mo.atom_begin(),mo.atom_end(),iparm_cmper1(ELEMENT,1));
    }


    void list_conn( mobjiter_t bgn, mobjiter_t end, conninfo_t& h, conninfo_t& o )
    {
        mobjiter_t i = bgn;
        for( ; i != end; ++i )
        {
            if( has_H(*i) )
            {
                h.push( *i );
            }
            else
            {
                o.push( *i );
            }
        }
    }

    void make_vdwp(const vector<double>& e, const vector<double>& r, vector<int>& id, vector<double>& c1, vector<double>& c2 )
    {
        assert( r.size()==e.size() );

        int ntype = r.size();
        int ntyp2 = ntype*ntype;
        int nttyp = ntype*(ntype+1)/2;

        id.resize(ntyp2);
        c1.resize(nttyp);
        c2.resize(nttyp);

        int s=0;
        for( int i=0; i < ntype; ++i )
        {
            for( int j=0; j <= i; ++j )
            {
                double r6 = std::pow ( r[i]+r[j], 6);
                double d = sqrt( e[i]*e[j] );
                c2[s] = 2.0*d*r6;
                c1[s] = d*r6*r6;

                s++;
                id[i*ntype+j] = s;
                id[j*ntype+i] = s;
            }
        }
    }

    char* pack( string& dst, const vector<string>& src, int width )
    {
        assert( dst.length()==0 );

        for( unsigned int i=0; i < src.size(); ++i )
        {
            if( (int)src[i].length() >= width )
            {
                dst += src[i].substr(0, 4);
            }
            else
            {
                dst += src[i];
                dst.append( width-src[i].length(), ' ' );
            }
        }

        assert( dst.length()==width*src.size() );
        return const_cast<char*>( dst.c_str() );
    }

    int* get_iptr( molecule_t& m, const hashid_t& cid, const hashid_t& parmid )
    {
        vector<int>& v = get_ivec( m, cid, parmid );
        return &v[0];
    }

    double* get_dptr( molecule_t& m, const hashid_t& cid, const hashid_t& parmid )
    {
        vector<double>& v = get_dvec( m, cid, parmid );
        return &v[0];
    }

    double* make_charge( const molecule_t& m, vector<double>& charge )
    {
        const vector<double>& tmp = get_dvec( m, ATOM, PCHG );
        charge.resize( m.natom() );
        for( int i=0; i < m.natom(); ++i )
        {
            charge[i] = tmp[i]*18.2223;
        }
        return &charge[0];
    }

    double* make_solty( int natom, vector<double>& solty )
    {
        solty.resize( natom );
        fill( solty.begin(), solty.end(), 0.0 );
        return &solty[0];
    }

    int* make_joint( int natom, vector<int>& joint )
    {
        joint.resize( natom );
        fill( joint.begin(), joint.end(), 0 );
        return &joint[0];
    }

    int* make_rotat( int natom, vector<int>& rotat )
    {
        rotat.resize( natom );
        fill( rotat.begin(), rotat.end(), 0 );
        return &rotat[0];
    }

    int* resd_begin( const molecule_t& m, vector<int>& resdp )
    {
        resdp.clear();
        if( m.nresd()==0 )
        {
            resdp.push_back( 0 );
        }
        else
        {
            resditer_t ri = m.resd_begin();
            for( ; ri != m.resd_end(); ++ri )
            {
                resdp.push_back( ri->atom_begin()->get_i(ID) );
            }

        }

        resdp.push_back( m.natom() );
        return &resdp[0];
    }

    int* list_usize( const molecule_t& m, vector<int>& usize )
    {
        usize.clear();
        if( m.nresd()==0 )
        {
            usize.push_back( m.natom() );
        }
        else
        {
            resditer_t ri = m.resd_begin();
            int ai = 0;
            while( ri != m.resd_end() )
            {
                resditer_t rj = std::find_if( ri + 1, m.resd_end(), iparm_cmper1(HEAD,0) );
                int aj=0;
                if( rj == m.resd_end() )
                {
                    usize.push_back( m.natom() - ai );
                }
                else
                {
                    aj = rj->atom_begin()->relid();
                    usize.push_back( aj - ai );
                }
            
                ri = rj;
                ai = aj;
            }
        }

        return &usize[0];
    }

 
    double* make_gvdw( const nabparm_t& prm, vector<double>& gvdw )
    {
        static const double SIGMAW   = 3.15365;
        static const double EPSILONW = 0.155;
        static const double RHOW     = 0.33428;
        static const double PI       = 3.141592650;

        vector<double> atype( prm.Ntypes);
        double sigmaw3 = SIGMAW * SIGMAW * SIGMAW;
        for (int i = 0; i < prm.Ntypes; i++) 
        {
            int iaci = prm.Cno[prm.Ntypes*i + i] - 1;
            if (prm.Cn1[iaci] == 0.0 || prm.Cn2[iaci] == 0.0) 
            {
                atype[i] = 0.0;
            } 
            else 
            {
                double sigma_iw6  = sigmaw3 * sqrt(prm.Cn1[iaci] / prm.Cn2[iaci]);
                double epsilon_iw = 0.5 * sqrt(EPSILONW/prm.Cn1[iaci]) * prm.Cn2[iaci];
                atype[i] = -16. * PI * RHOW * epsilon_iw * sigma_iw6 / 3.;
#if 0
         fprintf(stderr, "%5d  %15.8f  %15.8f  %15.8f\n", iaci, epsilon_iw,
                 pow(sigma_iw6, 1. / 6.), atype[i]);
#endif
            }
        }

        /*
         *
         *  Now, use these to fill in an array indexed by atom number:
         */

        gvdw.resize( prm.Natom );

        for(int i = 0; i < prm.Natom; i++) 
        {
            gvdw[i] = atype[prm.Iac[i] - 1];
        }

#if 0
   fprintf(stderr, "Gvdw values:\n");
   for (i = 0; i < prm->Natom; i++) {
      fprintf(stderr, "%5d  %15.8f\n", i + 1, prm->Gvdw[i]);
   }
#endif
        return &gvdw[0];            
    }



    /*
    * -------CONSTRUCT A 1-4 LIST -------
    */
    void make_14( const nabparm_t& prm, vector<int>& size, vector<int>& list )
    {
        vector< vector<int> > lst2;
        lst2.resize( prm.Natom );


        for (int i = 0; i < prm.Nphih; i++) 
        {
            int iat = prm.DihHAt1[i] / 3;
            int kat = prm.DihHAt3[i] / 3;
            int lat = prm.DihHAt4[i] / 3;
            if (kat >= 0 && lat >= 0) {
                int isml = std::min( iat, lat );
                int ibig = std::max( iat, lat );
                lst2[isml].push_back( ibig );
            }
        }

        for (int i = 0; i < prm.Mphia; i++) 
        {
            int iat = prm.DihAt1[i] / 3;
            int kat = prm.DihAt3[i] / 3;
            int lat = prm.DihAt4[i] / 3;
            if (kat >= 0 && lat >= 0) 
            {
                int isml = std::min( iat, lat );
                int ibig = std::max( iat, lat );
                lst2[isml].push_back( ibig );
            }
        }

        assert( lst2.back().size()==0 );
        for( int i=0; i < prm.Natom; ++i )
        {
            size.push_back( lst2[i].size() );
            for( unsigned int j=0; j < lst2[i].size(); ++j )
            {
                list.push_back( lst2[i][j] );
            }
        }


#ifdef PRINT_14PAIRS
        if (get_mytaskid() == 0) 
        {
            fprintf(nabout, "npairs:\n");
            for (k = 0; k < prm->Natom; k++) 
            {
                fprintf(nabout, "%4d", prm->N14pairs[k]);
                if ((k + 1) % 20 == 0)
                    fprintf(nabout, "\n");
            }

            fprintf(nabout, "\npairlist:\n");
            for (k = 0; k < idum; k++) 
            {
                fprintf(nabout, "%4d", prm->N14pairlist[k]);
                if ((k + 1) % 20 == 0)
                    fprintf(nabout, "\n");
            }
            fprintf(nabout, "\n");
        }
#endif
    }



    void energee_t::assignparm( const molecule_t& ffp, int ftype )
    {
        int lestype=0, solute = 0;

        // creating shortcuts to data member
        molecule_t& m = *m_pmol;
        parmset_t& ps = m_parmset;
        nabparm_t& prm = m_nabparm;

		int mol_sol_par = m.get_i(SOLUTE);

        // assign paramters
        parametrize( m, ffp, ps );

        // make exclusion list
        if( ftype==AMBER )
        {
            exclude( m, m_excl, 3 );
        }

        m_excl.finish( m_exsize, m_exlist );

        // make list of bond, angl and tors
        list_conn( m.bond_begin(), m.bond_end(), m_bondh, m_bondo );
        list_conn( m.angl_begin(), m.angl_end(), m_anglh, m_anglo );
        list_conn( m.dihe_begin(), m.dihe_end(), m_torsh, m_torso );
        list_conn( m.impr_begin(), m.impr_end(), m_torsh, m_torso );

        make_vdwp( ps.vdw[0], ps.vdw[1], m_vdwidx, m_vdwcn1, m_vdwcn2 );


        // standard amber prmtop header
        prm.Natom  	= m.natom();
        prm.Ntypes 	= ps.vdw[0].size();
        prm.Nbonh  	= m_bondh.size();
        prm.Mbona  	= m_bondo.size();
        prm.Ntheth 	= m_anglh.size();
        prm.Mtheta 	= m_anglo.size();
        prm.Nphih  	= m_torsh.size();
        prm.Mphia  	= m_torso.size();
        prm.Nhparm 	= 0; // ?
        prm.Nparm  	= m.get_i(LESTYPE,lestype)? lestype:0 ;
        prm.Nnb    	= m_exlist.size(); // size of nonbond exclusion list
        prm.Nres   	= std::max(1, m.nresd());
        prm.Nbona  	= m_bondo.size(); // use to be constrain bond number, now same as Mbona
        prm.Ntheta 	= m_anglo.size();
        prm.Nphia  	= m_torso.size();
        
        prm.Numbnd 	= ps.bond[0].size();
        prm.Numang 	= ps.angl[0].size();
        prm.Nptra  	= ps.tors[0].size();

        prm.Natyp  	= ps.vdw[0].size();
        prm.Nphb   	= 0; // 10-12 HBond, now obseleted

        prm.IfBox  	= m.get_i(SOLUTE,solute) && (solute==BOX);
        prm.Nmxrs 	= max_resd_size( m );
        prm.IfCap  	= (solute==CAP);
        prm.Numextra 	= std::count_if(m.atom_begin(), m.atom_end(), sparm_cmper1(NAME,"EP") );

        
        // some variable for memory allocation
        prm.Nat3   	= prm.Natom*3;
        prm.Ntype2d	= prm.Ntypes*prm.Ntypes;
        prm.Nttyp  	= prm.Ntypes*(prm.Ntypes+1)/2;
        
	// start making pointers
        prm.AtomNames	= pack( m_atomname, get_svec(m,ATOM,NAME), 4 );
        prm.Charges   	= make_charge( m, m_charge );
        prm.Masses	= get_dptr(m, ATOM, MASS);
        prm.Iac 	= get_iptr(m, ATOM, TYPEID);
        prm.Iblo 	= &m_exsize[0];
        prm.Cno 	= &m_vdwidx[0];

        if( m.nresd()==0 )
        {
            prm.ResNames= (char*)"UNK";
        }
        else
        {
            prm.ResNames= pack( m_resdname, get_svec(m,RESD,TYPE), 4 );
        }

        // must be nresd + 1.
        prm.Ipres 	= resd_begin( m, m_resdp);

        prm.Rk  	= &ps.bond[0][0];
        prm.Req 	= &ps.bond[1][0];

        prm.Tk 		= &ps.angl[0][0];
        prm.Teq 	= &ps.angl[1][0];

        if( !ps.tors[0].empty() ) prm.Pk   	= &ps.tors[0][0];
        if( !ps.tors[1].empty() ) prm.Pn   	= &ps.tors[1][0];
        if( !ps.tors[2].empty() ) prm.Phase	= &ps.tors[2][0];
//        prm.Scee	= &ps.tors[3][0];  // IGOR TEMP
//        prm.Scnb	= &ps.tors[4][0];  // IGOR TEMP

        prm.Solty 	= make_solty( m.natom(), m_solty );

        prm.Cn1 	= &m_vdwcn1[0];
        prm.Cn2 	= &m_vdwcn2[0];
        
	if( m_bondh.size() > 0 )
        {
            prm.BondHAt1 	= &m_bondh.id3[0][0];
            prm.BondHAt2 	= &m_bondh.id3[1][0];
            prm.BondHNum 	= &m_bondh.typ[0];
        }

        if( m_bondo.size() > 0 )
        {
            prm.BondAt1  	= &m_bondo.id3[0][0];
            prm.BondAt2  	= &m_bondo.id3[1][0];
            prm.BondNum  	= &m_bondo.typ[0];
        }

        if( m_anglh.size() > 0 )
        {
            prm.AngleHAt1 	= &m_anglh.id3[0][0];
            prm.AngleHAt2 	= &m_anglh.id3[1][0];
            prm.AngleHAt3 	= &m_anglh.id3[2][0];
            prm.AngleHNum 	= &m_anglh.typ[0];
        }

        if( m_anglo.size() > 0 )
        {
            prm.AngleAt1  	= &m_anglo.id3[0][0];
            prm.AngleAt2  	= &m_anglo.id3[1][0];
            prm.AngleAt3  	= &m_anglo.id3[2][0];
            prm.AngleNum  	= &m_anglo.typ[0];
        }

        if( m_torsh.size() > 0 )
        {
            prm.DihHAt1 	= &m_torsh.id3[0][0];
            prm.DihHAt2 	= &m_torsh.id3[1][0];
            prm.DihHAt3 	= &m_torsh.id3[2][0];
            prm.DihHAt4 	= &m_torsh.id3[3][0];
            prm.DihHNum 	= &m_torsh.typ[0];
        }

        if( m_torso.size() > 0 )
        {
            prm.DihAt1  	= &m_torso.id3[0][0];
            prm.DihAt2  	= &m_torso.id3[1][0];
            prm.DihAt3  	= &m_torso.id3[2][0];
            prm.DihAt4  	= &m_torso.id3[3][0];
            prm.DihNum  	= &m_torso.typ[0];
        }

	if( m_exlist.size() > 0 )
        {
            prm.ExclAt  	= &m_exlist[0];
        }

        prm.AtomSym  	= pack( m_atomtype, get_svec(m, ATOM, TYPE), 4 );
        prm.AtomTree	= pack( m_atomtree, get_svec(m, ATOM, TREE), 4 );
        prm.TreeJoin 	= make_joint( m.natom(), m_joint );
        prm.AtomRes  	= make_rotat( m.natom(), m_rotat );

        if( prm.IfBox )
        {
            int nresd 	= m.nresd();
            int nsolv 	= m.nresd() - get_solute_nresd(m);
            int nunit 	= std::count_if( m.resd_begin(), m.resd_end(), iparm_cmper1(HEAD,0) );
            prm.Iptres 	= nresd - nsolv;
            prm.Nspm   	= nunit;
            prm.Nspsol	= nunit - nsolv + 1;
            numvec b 	= m.get_v(BOX);
            prm.Box[0] 	= b[0];
            prm.Box[1] 	= b[1];
            prm.Box[2] 	= b[2];
            prm.Box[3]	= b[3];

            if( prm.Box[3] > 89.5 && prm.Box[3] < 90.5 )
            {
                prm.IfBox = 1;
            }
            else if( prm.Box[3] > 109.0 && prm.Box[3] < 110.0 )
            {
                prm.IfBox = 2;
            }
            else
            {
                throw std::runtime_error( "Error: unknown box angle, it is neither 90.0 nor 109.5." );
            }

        }
        else
        {
            prm.Iptres	= 0;
            prm.Nspm 	= 1;
            prm.Box[0] 	= prm.Box[1] = prm.Box[2] = 999999999.;
        }


        prm.Boundary 	= list_usize( m, m_usize );

        if( prm.Iptres )
            prm.Ipatm 	= prm.Ipres[prm.Iptres]-1;

        if( prm.IfCap )
        {
            prm.Natcap 	= get_solute_natom( m );
            numvec cap 	= m.get_v(CAP);
            prm.Cutcap 	= cap[3];
            prm.Xcap 	= cap[0];
            prm.Ycap 	= cap[1];
            prm.Zcap 	= cap[2];
        }

        const double BOFFSET = 0.09;
        prm.Rborn 	= get_dptr( m, ATOM, BRNR );
        prm.Fs    	= get_dptr( m, ATOM, GBFS );
        prm.Fsmax 	= 0.0;
        for( int i=0; i < prm.Natom; ++i )
        {
            double tmp 	= prm.Fs[i] * (prm.Rborn[i]-BOFFSET);
            prm.Fsmax 	= std::max( tmp, prm.Fsmax );
        }


        // post processing, 
        if( !prm.IfBox )
        {
            prm.Gvdw = make_gvdw( prm, m_gvdw );
        }

        make_14( prm, m_14size, m_14list );
        prm.N14pairs = &m_14size[0];
        if( !m_14list.empty() )prm.N14pairlist = &m_14list[0];


    }

} // namespace mort



