#include <stdexcept>
#include <common.hpp>
#include <object.hpp>
#include "mdlize.hpp"
#include "namemap.hpp"

namespace mort
{
    using std::runtime_error;

    // mdlize (modelize) contain routines to build up a molecule
    // according to a set of predefined models or (templates).
    // Actions includs: 1. add missing atoms and position them
    // properly (usually by superimpose its neighbours); 
    // 2. copy the parameters from the models: atom name, type,
    // charge, etc.


   
    // modelize a molecule according to a database
    void mdlize_mdb( molecule_t& m, const database_t& mdb )
    {
        funstack_t::push( "mdlize_mdb" );
        
        if( !mdb.has("_namemap") )
        {
            std::cerr << "Warning: namemap does not exist when modelize resd" << std::endl;
        }

        map<int,int> idmap;
        molecule_t full;

        resd_t rprev(full, -1);
        resditer_t ri = m.resd_begin();
        for( ; ri != m.resd_end(); ++ri )
        {
            resd_t rcurt = mdlize_resd( *ri, mdb, full, idmap );
            conect_resd( rprev, rcurt );
            rprev = rcurt;
        }

        m.swap( full );
        funstack_t::pop();
    }

    void mdlize_seq( molecule_t& m, const database_t& seq, const namemap_t& nmap )
    {
        if( m.nresd() != seq.size() )
        {
            throw runtime_error( "Error: while building by sequence, nresd != seq.size" );
        }

        map<int,int> idmap;
        molecule_t full;

        resd_t rprev(full, -1);
        database_t::const_iterator ei = seq.begin();
        for( resditer_t ri=m.resd_begin(); ri != m.resd_end(); ++ri, ++ei )
        {
            molecule_ptr pm = dynamic_pointer_cast< molecule_t >( ei->second );
            assert( pm != NULL );
            resd_t rcurt = mdlize_resd( *ri, *pm, nmap, full, idmap );
            conect_resd( rprev, rcurt );
            rprev = rcurt;
        }

        m.swap( full );
    }

    resd_t mdlize_resd( resd_t& r, const database_t& mdb, molecule_t& full, map<int,int>& idmap )
    {
        funstack_t::push( "mdlize_resd_bymdb" );

        // namemap is used to get the proper form of a residue .
        // Each type of amino acid has at least 3 forms: normal form
        // N-terminal form and C-terminal form. Some residues have 
        // more forms. For example histaine on different prontate 
        // state could be HIA, HID and HIE.
        // Meanwhile, nucleic acid residues usually has different nick
        // names.
        // namemap handle all these complexities.
        namemap_ptr n = mdb.get_nmap( "_namemap" );
        if( n==NULL) 
        {
            n  = namemap_ptr( new namemap_t() );
        }

        resd_t rout(full, -1);
        try
        {
            string mdlname = n->get_name(r);

            molecule_ptr mdl = mdb.get_mol( mdlname );
            rout = mdlize_resd( r, *mdl, *n, full, idmap );

            if( mdlname.length() < 4 )
            {
                rout.set_s( TYPE, mdlname );
            }
            else
            {
                rout.set_s( TYPE, mdlname.substr(mdlname.length()-3, 3) );
            }

        }
        catch( std::exception&  )
        {
            string fastbld( "off" );
            mortenv().get_s( "fastbld", fastbld );

            if( fastbld=="on" )
            {
                // this is a new feature of sleap. when a mdl is not availiable
                // it can assign bond order and parameter automatically.
                std::cerr << "Warning: cannot find template for residue " + r.get_s(TYPE) << " in our library." << std::endl;
                std::cerr << "         I am going to construct the template by myself." << std::endl;
 
                molecule_t frg = unmerge( r );
                fixbond( frg );
                addHs  ( frg );


                std::cerr << "         OK. I have fixed bond order and added hydrogens." << std::endl;
		std::cerr << "         Now I am calling antechamber to assign atom type and partial charge." << std::endl;
                setpchg( frg );
 
 
                std::cerr << "         OK. Partial charge and atom types has been assigned." << std::endl;
		std::cerr << "         Please run parmchk before saving topology file." << std::endl;

                rout = merge( full, frg );
            }
            else
            {
                std::cerr << "Warning: cannot find template for residue " + r.get_s(TYPE) << " in our library." << std::endl;
                std::cerr << "         You will not be able to save prmtop for this molecule." << std::endl;
                rout = merge( full, r );
            }
        }

        funstack_t::pop();
        return rout;
    }

    resd_t mdlize_resd( resd_t& r, const molecule_t& m, const namemap_t& nmap, molecule_t& full, map<int,int>& idmap )
    {
        funstack_t::push( "mdlize_resd_bymol" );

        map< string, numvec > pos;
        map< string, int    > seq;
        map< string, int    > ids;

        // backup information of the molecule.
        atomiter_t ai = r.atom_begin();
        for( ; ai != r.atom_end(); ++ai )
        {
            string name = nmap.get_name( *ai );
            if( !atom_t::has(m, name) )
            {
                std::cerr << "Warning: there is an atom "  << ai->get_s(NAME) << " in residue " + r.name();
                std::cerr << " which does not exist in the model " << m.get_s(NAME) << "!" << std::endl;
                std::cerr << "         This atom will not be created." << std::endl;
            }

            pos[ name ] = ai->get_v(POSITION);
            seq[ name ] = ai->get_i(SEQUENCE);
            ids[ name ] = ai->absid();
        }

        if( r.natom() > m.natom() )
        {
            std::cerr << "residue " << r.name() << " has " << r.natom() << " atoms" << std::endl;
            for( int i=0; i < r.natom(); ++i )
            {
                std::cerr << "       " << r.atoms()[i].name() << std::endl;
            }

            std::cerr << "model   " << m.get_s(NAME) << " has " << m.natom() << " atoms." << std::endl;
            for( int i=0; i < m.natom(); ++i )
            {
                std::cerr << "       " << m.atoms()[i].name() << std::endl; 
            }
            throw runtime_error( "Error: there are more atom in residue than in template." );
        }

        resd_t r2 = merge( full, m );

        vector<numvec> crd;
        ai = r2.atom_begin();
        for( ; ai != r2.atom_end(); ++ai )
        {
            string name = ai->name();
      
            map<string,numvec>::iterator i = pos.find(name);
            if( i!=pos.end() )
            {
                crd.push_back( i->second );
                idmap[ ids[name] ] = ai->absid();
                assert( crd.back().size()==3u );
            }
            else
            {
                crd.push_back( impose(*ai, pos) );
                pos[name] = crd.back();
                assert( crd.back().size()==3u );
            }
        }

        // copy backup informations back.
        copy_sparm( r, NAME, r2 );
        copy_sparm( r, TYPE, r2 );
        copy_iparm( r, AAPOS, r2 );
        copy_iparm( r, CHAIN, r2 );
        if( r.get_i(AAPOS)==CTERM || r.get_i(AAPOS)==ALONE )
        {
            r2.set_i(TAIL, 0);
        }

        if( r.get_i(AAPOS)==NTERM || r.get_i(AAPOS)==ALONE )
        {
            r2.set_i(HEAD, 0);
        }


        ai = r2.atom_begin();
        for( int i=0; ai != r2.atom_end(); ++ai, ++i )
        {
            ai->set_v( POSITION, crd[i] );
            ai->set_i( SEQUENCE, seq[ai->name()] );
        }

        // make inter residue bond
        ai = r.atom_begin();
        for( ; ai != r.atom_end(); ++ai )
        {
            atomiter_t aj = ai->atom_begin();
            for( ; aj != ai->atom_end(); ++aj )
            {
                if( idmap.find(aj->absid())==idmap.end() )
                    continue;

                bond_t b1 = bond_t::get( *ai, *aj );
                atom_t a1( full, idmap[ai->absid()] );
                atom_t a2( full, idmap[aj->absid()] );
                if( bond_t::has(a1, a2) ) continue;
                bond_t b2 = bond_t::create( a1, a2 );
                copy_allparms( b1, b2 );
            }
        }

        funstack_t::pop();
        return r2;
    }

    void conect_resd( resd_t& rprev, resd_t& rcurt )
    {
        if( rprev.absid()==-1 )
        {
            return;
        }
        
        int id1 = rprev.tail();
        int id2 = rcurt.head();
        if( id1!=0 && id2!=0 )
        {
            atom_t a1( rprev.getmol(), id1-1 );
            atom_t a2( rprev.getmol(), id2-1 );

			if( !mort::bond_t::has(a1,a2) )
			{
				bond_t b = bond_t::create( a2, a1 );
				b.set_i(ORDER, 1);
			}
        }
    }


    numvec impose( const atom_t& a, const map<string, numvec>& pos )
    {
        atomvec_t nbrs = bfv( a, pos );


        vector<numvec> spos( nbrs.size()-1, numvec(3) );
        vector<numvec> dpos( nbrs.size()-1, numvec(3) );
        for( unsigned int i=1; i < nbrs.size(); ++i )
        {
            string name = nbrs[i].get_s(NAME);
            spos[i-1] = nbrs[i].get_v(POSITION);
            dpos[i-1] = pos.find(name)->second;
        }

        return impose( spos, dpos, a.get_v(POSITION) );
    }

    // breath first visit.
    atomvec_t bfv( const atom_t& a, const map<string, numvec>& pos )
    {
        atomvec_t nbrs(1, a);

        int bgn = 0;
        int end = nbrs.size();
        while( nbrs.size() < 4 && bgn < end)
        {
            for(int i=bgn; i < end; ++i )
            {
                atomiter_t ni = nbrs[i].atom_begin();
                for( ; ni != nbrs[i].atom_end(); ++ni )
                {
                    // skip visited
                    if( nbrs.find(*ni) != nbrs.end() )
                    {
                        continue;
                    }

                    // keep those positions are known
                    string name = ni->get_s(NAME);
                    if( pos.find(name)!=pos.end() )
                    {
                        nbrs.push_back(*ni);
                        if( nbrs.size()==4 ) return nbrs;
                    }
                }
            }

            bgn = end;
            end = nbrs.size();
        }

        return nbrs;
    }



} // namespace mort


