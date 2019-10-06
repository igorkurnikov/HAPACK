
#include <boost/algorithm/string.hpp>

#include <common.hpp>
#include <object.hpp>

#include "les.hpp"
#include "parm.hpp"
#include "extent.hpp"
#include "parmset.hpp"

namespace mort
{
    using std::find_if;
    using std::count_if;
    using std::for_each;

    string str( const morf_t& mo );

    namespace prmtop
    {

        void parm_atom( const molecule_t* pff, atom_t& a )
        {
            string type = a.type();
 	    atom_t patm(*pff, -1);
	    if( ! atom_t::get(*pff, type, patm) )
            {
                throw std::runtime_error( "Error: unknown atom type: " + type );
            }

            copy_dparm(patm, MASS,  a );
            copy_dparm(patm, POLAR, a );
            copy_dparm(patm, DEPTH, a );
            copy_dparm(patm, RSTAR, a );
            copy_dparm(patm, VBUFF, a );
            copy_iparm(patm, ELEMENT, a );
            a.set_i(TYPEID, patm.get_i(ID) );
        }

        void parm_gbsa(atom_t& a, const hashid_t& gbtype)
        {
            static const int elements[] = { HYDROGEN, CARBON, NITROGEN, OXYGEN, FLUORINE, 
                                            SILICON, PHOSPHORUS, SULFUR, CHLORINE };
        
            static const double radius[] = { 1.20, 1.70, 1.55, 1.50, 1.50,
                                             2.10, 1.85, 1.80, 1.70, 1.50 };
            
            static const double screen[] = { 0.85, 0.72, 0.79, 0.85, 0.88,
                                             0.80, 0.86, 0.96, 0.80, 0.80 };

            int elem = a.element();
            const int* ptr = std::find(elements, elements+sizeof(elements)/sizeof(int), elem);
            a.set_d(BRNR, radius[ptr-elements] );
            a.set_d(GBFS, screen[ptr-elements] );

            if( elem == HYDROGEN && gbtype != BONDI && a.natom() > 0 )
            {
                int nbre = a.atom_begin()->element();        
                if( gbtype == AMBER6 || gbtype == MBONDI )
                {
                    if( nbre == 6  ) a.set_d(BRNR, 1.3);
                    if( nbre == 8  ) a.set_d(BRNR, 0.8);
                    if( nbre == 16 ) a.set_d(BRNR, 0.8);
                    if( nbre == 7 && gbtype == MBONDI ) a.set_d(BRNR, 1.3);
                }
                else
                {
                    assert( gbtype == MBONDI2 );
                    if( nbre == 7 ) a.set_d(BRNR, 1.3);
                }
            }
        }

        void parm_bond( const molecule_t* pff, bond_t& b )
        {
            atom_t a1 = b.atoms()[0];
            atom_t a2 = b.atoms()[1];
            atom_t pa1(*pff, a1.get_i(TYPEID)-1);
	    atom_t pa2(*pff, a2.get_i(TYPEID)-1);
            bond_t pb (*pff, -1);
            if( !bond_t::get(pa1, pa2, pb) )
            {
                std::cout << "bond parm not found for ";
                std::cout << a1.get_s(RESNAME) << ":" << a1.get_s(NAME) << "-";
                std::cout << a2.get_s(RESNAME) << ":" << a2.get_s(NAME) << std::endl;

                throw std::runtime_error( "Error: can't find parameter for bond " 
                                   + a1.type() + "-" 
                                   + a2.type() );
            } 

            copy_allparms( pb, b );
            numvec opbend(2);
            if( pb.get_v(OPBEND, opbend) )
            {
                if( pb.atoms()[0].get_i(ID)==b.atoms()[0].get_i(TYPEID) )
                {
                    b.set_v(OPBEND, opbend);
                }
                else
                {
                    b.set_v(OPBEND, makevec(opbend[1], opbend[0]) );
                }
            }
        }

        void parm_angl( const molecule_t* pff, atmvec& tmp )
        {
            assert( tmp.size() == 3 );
	    atmvec as(tmp);
            
            angl_t an = angl_t::create( as[0], as[1], as[2] );
            set_les( an );

            atom_t pa1( *pff, as[0].get_i(TYPEID)-1 );
            atom_t pa2( *pff, as[1].get_i(TYPEID)-1 );
	    atom_t pa3( *pff, as[2].get_i(TYPEID)-1 );

	    angl_t pan(*pff, -1);
	    if( !angl_t::get( pa1, pa2, pa3, pan) )
	    {
                morf_t a1 = an.atoms()[0];
                morf_t a2 = an.atoms()[1];
                morf_t a3 = an.atoms()[2];
                std::cout << "angl parm not found for ";
                std::cout << a1.get_s(RESNAME) << ":" << a1.get_s(NAME) << "-";
                std::cout << a2.get_s(RESNAME) << ":" << a2.get_s(NAME) << "-";
                std::cout << a3.get_s(RESNAME) << ":" << a3.get_s(NAME) << std::endl;

                throw std::runtime_error( "Error: can't find parameter for angl "
                                       + as[0].get_s(TYPE) + "-" 
                                       + as[1].get_s(TYPE) + "-"
                                       + as[2].get_s(TYPE) );
            }

	    copy_dparm(pan, FORCE, an);
	    
            numvec addition(2);
            if( pan.get_v(ADDITION, addition) )
            {
                int nhang = count_if( an.atom_begin(), an.atom_end(), iparm_cmper1(ELEMENT, HYDROGEN) );
                int nhnbr = count_if( an.atoms()[1].atom_begin(), an.atoms()[1].atom_end(), iparm_cmper1(ELEMENT, HYDROGEN) );
                int nh = nhnbr - nhang;

                if( nh == 1 && addition[0] != 0.0 )
                {
                    an.set_d(EQUIL, addition[0]);
                }
                else if( nh == 2 && addition[1] != 0.0 )
                {
                    an.set_d(EQUIL, addition[1]);
                }
                else
                {
                    copy_dparm(pan, EQUIL, an);
                }
            }
            else
            {
                copy_dparm(pan,EQUIL,an);
            }
            
            copy_vparm(pan, UREY, an);

            numvec sbparm(3);
            if( pan.atoms()[1].get_v(STRBND, sbparm) )
            {
                int nh = count_if( an.atom_begin(), an.atom_end(), iparm_cmper1(ELEMENT, HYDROGEN) );
                an.set_v(STRBND, makevec(sbparm[nh], 0.0, 0.0) );
            }
        }


        void compress( vector<dihe_t>& pdihs )
        {
            assert( pdihs.size() > 0 );

            string dummy;
            bool verbose = mortenv().get_s("ffverbose", dummy) && dummy=="on";

            if( pdihs.size()==1 )
            {
                // only one matching parm
                return;
            }

            vector<dihe_t> genes;
            vector<dihe_t> specs;
            for( unsigned int i=0; i < pdihs.size(); ++i )
            {
                int s = score( pdihs[i] );
                if( s==4 )
                {
                    specs.push_back( pdihs[i] );
                }
                else
                {
                    genes.push_back( pdihs[i] );
                }
            }

            if( genes.size()==0 || specs.size()==0 )
                return;

            for( unsigned int i=0; i < genes.size(); ++i )
            {
                if( genes[i].get_d(FORCE) < 1e-6 )
                    continue;

                bool found=false;
                int p = genes[i].get_i(PERIOD);
                for( unsigned int j=0; j < specs.size(); ++j )
                {
                    if( specs[j].get_i(PERIOD)==p )
                    {
                        found = true;
                        break;
                    }
                }

                if( !found && verbose )
                {
                    std::cout << "NOTE: generic torsion term " << str(genes[i]) << " will be overwritten by\n";
                    std::cout << "NOTE:        specific term " << str(specs[i]) << std::endl;
                }
            }

            pdihs.swap( specs );
        }


        void parm_dihe( const molecule_t* pff, atmvec& as )
        {
            assert( as[1] != as[2] );
            assert( as[3] != as[0] && as[3] != as[1] );

            atmvec pas(4, atom_t(*pff, -1) );
	    pas[0] = atom_t(*pff, as[0].get_i(TYPEID)-1);
            pas[1] = atom_t(*pff, as[1].get_i(TYPEID)-1);
            pas[2] = atom_t(*pff, as[2].get_i(TYPEID)-1);
            pas[3] = atom_t(*pff, as[3].get_i(TYPEID)-1); 

            bond_t pbc = bond_t::get(pas[1], pas[2]);
            
            vector< dihe_t > pdihs;   
            mort::copy_if( pbc.dihe_begin(), pbc.dihe_end(), back_inserter(pdihs), sequence_Match2(pas) );
            
            if( pdihs.size() == 0 )
            {
                throw std::runtime_error( "Error: can't find parameter for tors " 
                                   + as[0].get_s(TYPE) + "-" 
                                   + as[1].get_s(TYPE) + "-"
                                   + as[2].get_s(TYPE) + "-"
                                   + as[3].get_s(TYPE) );
            }

            std::sort( pdihs.begin(), pdihs.end(), by_score() );

            compress( pdihs );


            for( unsigned int i=0; i < pdihs.size(); ++i )
            {
                int period = pdihs[i].period();

                if( dihe_t::has(as[0],as[1],as[2],as[3],period) )
                {
                    continue;
                }

                if( pdihs[i].get_d(FORCE) > 1e-6 || score(pdihs[i]) == 4 || pdihs.size() == 1 )
                {
                    dihe_t di = dihe_t::create( as[0], as[1], as[2], as[3], period );
                    copy_dparm(pdihs[i], FORCE, di);
                    copy_dparm(pdihs[i], EQUIL, di);
                    copy_dparm(pdihs[i], SCEE,  di);
                    copy_dparm(pdihs[i], SCNB,  di);
                    set_les( di );
                }                
            }
        }

        void parm_tor2( const molecule_t* pff, atomvec_t& atoms )
        {
            assert( atoms[1]!=atoms[2] && atoms[3]!=atoms[0] && atoms[3]!=atoms[1] );

            mobjiter_t parm = find_if( pff->tor2_begin(), pff->tor2_end(), sequence_Match(atoms) );
            
            if( parm==pff->tor2_end() )
            {
                return;
            }
            
            morf_t tor2 = create_tor2( atoms );
            tor2.set_i(TYPEID, parm->absid()+1 );
        }

        void adjust_order( const morf_t& parm, atomvec_t& atoms )
        {
            if( atoms[0] > atoms[3] && atoms[0].get_s(TYPE) == atoms[3].get_s(TYPE) )
            {
                std::swap( atoms[0], atoms[3] );                
            }

            if( atoms[1] > atoms[3] && atoms[1].get_s(TYPE) == atoms[3].get_s(TYPE) )
            {
                std::swap( atoms[1], atoms[3] );
            }

            if( atoms[0] > atoms[1] && atoms[0].get_s(TYPE) == atoms[1].get_s(TYPE) )
            { 
                std::swap( atoms[0], atoms[1] );
            }

            if( atoms[0] > atoms[1] && atom_1st(parm).get_s(NAME)=="X" && atom_2nd(parm).get_s(NAME)=="X" )
            {
                std::swap( atoms[0], atoms[1] );
            }
        }

        void parm_amoeba_impr( atomvec_t& as )
        {
            if( as[0] > as[1] )
            {
                return;
            }

            if( as[2].natom( ) != 3 )
            {
                return;
            }

            bond_t cent = bond_t::get( as[2], as[3] );

            numvec opbend(2);
	    if( ! cent.get_v(OPBEND, opbend) )
	    {
	        throw std::runtime_error( "bond do not have opbend" );
	    }
            
            double force = (atom_1st(cent) == as[2]) ? opbend[0] : opbend[1];
            impr_t::create(as[0],as[1],as[2],as[3],2).set_v(OPBEND, makevec(force, 0.0) );
        }

        void parm_sander_impr( const molecule_t* pff, atomvec_t& as, int& max_score )
        {
            atom_t a3(*pff, as[2].get_i(TYPEID)-1);
            atom_t a4(*pff, as[3].get_i(TYPEID)-1);
            bond_t cent = bond_t::get(a3, a4);

            vector<impr_t> pimps;
            mort::copy_if( cent.impr_begin(), cent.impr_end(), back_inserter(pimps), sequence_Match(as, MONO_DIRECTION) );

            if( pimps.size() == 0 )
            {
                return;
            }            

            std::sort( pimps.begin(), pimps.end(), by_score() );

            if( max_score > score(pimps[0]) )
            {
                return;
            }
                
            int period = pimps[0].get_i(PERIOD);

            max_score = score( pimps[0] );
                
            adjust_order( pimps[0], as );

            impr_t im = impr_t::frcget(as[0], as[1], as[2], as[3], period);
            copy_dparm(pimps[0], FORCE, im);
            copy_dparm(pimps[0], EQUIL, im);
            copy_dparm(pimps[0], SCEE,  im);
            copy_dparm(pimps[0], SCNB,  im);
            set_les( im );
        }
        
    } // namespace prmtop

    using namespace prmtop;

    bool fixed_water()
    {
        if( ! mortenv().has_s( "FlexibleWater" ) )
        {
            return true;
        }
            
        return mortenv().get_s( "FlexibleWater" ) != "on";
    }

    hashid_t get_gbtype( )
    {
        if( ! mortenv().has_s( "PBradii" ) )
        {
            return MBONDI;
        }

        string parm = to_upper_copy( string( mortenv().get_s( "PBradii" ) ) );
            
        if( parm == "MBONDI" )
        {
            return MBONDI;
        }
        else if( parm == "AMBER6" )
        {
            return AMBER6;
        }
        else if( parm == "MBONDI2" )
        {
            return MBONDI2;
        }
        else if( parm == "BONDI" )
        {
            return BONDI;
        }
        else
        {
            throw std::runtime_error( "Error: unknown PBradii " + parm );
        }
    }

    void assign_ffparm(molecule_t& mol, const molecule_t& ff)
    {
        hashid_t gbtype = get_gbtype();
        atomiter_t ai = mol.atom_begin();
        atomiter_t ae = mol.atom_end();
        for( ; ai != ae; ++ai )
        {
//			std::cout << " atom: " << (*ai).get_s(RESNAME) << ":" << (*ai).get_s(NAME) << std::endl;
            parm_atom( &ff, *ai );
            parm_gbsa( *ai, gbtype );
        }

        bonditer_t bi = mol.bond_begin();
        bonditer_t be = mol.bond_end();
        for( ; bi != be; ++bi )
        {
            parm_bond( &ff, *bi );
        }

        bool fw = fixed_water();
        ai = mol.atom_begin();
        ae = mol.atom_end();
        for( ; ai != ae; ++ai )
        {
            atom2angl( &ff, *ai, fw );
        }

        bi = mol.bond_begin();
        be = mol.bond_end();
        for( ; bi != be; ++bi )
        {
            bond2dihe( &ff, *bi );
        }

        ai = mol.atom_begin();
        ae = mol.atom_end();
        for( ; ai != ae; ++ai )
        {
            atom2impr( &ff, *ai );
        }

        string name = ff.get_s(NAME);        
 
        if(name.find("amoeba") != string::npos)
        {
            angliter_t ani = mol.angl_begin();
            angliter_t ane = mol.angl_end();
            for( ; ani != ane; ++ani )
            {
                angl2tor2( &ff, *ani );
            }

            bi = mol.bond_begin();
            be = mol.bond_end();
            for( ; bi != be; ++bi )
            {
                bond2ptor( *bi );
            }
        }
    }
        
    void parametrize( molecule_t& mol, const molecule_t& ff, parmset_t& ps )
    {
        assign_ffparm(mol, ff);

        int lestype;
        if( mol.get_i(LESTYPE, lestype) && lestype == 1 )
        {
            les_scale(mol, ps);
        }

        for_each( mol.atom_begin(), mol.atom_end(), record_t(TYPEID, ps.vdw,  list_parm(DEPTH, RSTAR) ) );
        for_each( mol.bond_begin(), mol.bond_end(), record_t(TYPEID, ps.bond, list_parm(FORCE, EQUIL) ) );
        for_each( mol.angl_begin(), mol.angl_end(), record_t(TYPEID, ps.angl, list_parm(FORCE, EQUIL) ) );
        
        mark_chain( mol );        

        string name = ff.get_s(NAME);        
        
        if( name.find( "amber" ) != string::npos )
        {
            for_each( mol.dihe_begin(), mol.dihe_end(), record_t(TYPEID, ps.tors, list_parm(FORCE,PERIOD,EQUIL,SCEE,SCNB)));
            for_each( mol.impr_begin(), mol.impr_end(), record_t(TYPEID, ps.oops, list_parm(FORCE,PERIOD,EQUIL,SCEE,SCNB)));
            pack_tors_oops_parm( ps.tors, ps.oops, mol );
        }
        else
        {
            assert( name.find( "amoeba" ) != string::npos );
            for_each( mol.dihe_begin(), mol.dihe_end(), record_t(TYPEID, ps.tors, list_tors() ) );
            for_each( mol.impr_begin(), mol.impr_end(), record_t(TYPEID, ps.oops, list_oops() ) );
            for_each( mol.angl_begin(), mol.angl_end(), record_t(UREYID, ps.urey, list_urey() ) );
            for_each( mol.angl_begin(), mol.angl_end(), record_t(STRBNDID, ps.strbnd, list_strbnd() ) );
            for_each( mol.ptor_begin(), mol.ptor_end(), record_t(TYPEID, ps.ptor, list_parm(FORCE, PERIOD, EQUIL) ) );
        }

    }


} // namespace mort

