#include <boost/algorithm/string.hpp>
#include <iomanip>
#include <sstream>

#include <object.hpp>
#include <common.hpp>

#include "parm.hpp"
#include "prmtop.hpp"
#include "exclude.hpp"

namespace mort
{
    using std::endl;
    using std::find_if;
    using std::for_each;
    using std::count_if;

    namespace prmtop
    {
        void write_float_array( ostream& os, const vector<double>& param )
        {
            fortran_t format( os, "5E16.8" );

            for( int i=0; i < (int)param.size(); ++i )
            {
                format( param[i] );
            }
            
            format.end();
           
        }



        template< typename select_t, typename makeid_t >
        void write_object( fortran_t& format, mobjiter_t begin, mobjiter_t end, select_t select, makeid_t makeid )
        {
            vector< int > ids;
                
            for( ; begin != end; ++begin )
            {
                ids.clear();
            
                if( select( *begin ) )
                {
                    makeid( *begin, ids );
                
                    for( int i=0; i < (int)ids.size(); ++i )
                    {
                        format( ids[i] );
                    }                
                }
            }        
        }


        void write_ffparm( ostream& os, const vector< vector< double > >& params, int nparam, const string& title )
        {
            os << "%FLAG " << title << "_FORCE_CONSTANT" << endl;
            os << "%FORMAT(5E16.8)" << endl;
            
            write_float_array( os, params[0] );

            if( nparam > 1 )
            {
                os << "%FLAG " << title;

                if(nparam==2)
                {
                    os << "_EQUIL_VALUE" << endl;
                }
                else if(nparam==3)
                {
                    os << "_PERIODICITY" << endl;
                }
                else
                {
                    assert(nparam==4);
                    os << "_ANGLE_EQUIL_VALUE" << endl;
                }
                
                os << "%FORMAT(5E16.8)" << endl;
                write_float_array( os, params[1] );
            }
            
            if( nparam > 2 )
            {
                os << "%FLAG " << title;
                os << ( (nparam==3) ? "_PHASE" : "_BOND1_EQUIL_VALUE" );
                os << std::endl;
                os << "%FORMAT(5E16.8)" << endl;
                write_float_array( os, params[2] );
            }
            
            if( nparam > 3 )
            {
                assert( nparam==4 );
                os << "%FLAG " << title << "_BOND2_EQUIL_VALUE" << endl;
                os << "%FORMAT(5E16.8)" << endl;
                write_float_array( os, params[3] );
            }
        }
 
        namespace amoeba
        {

            void write_head( ostream& os )
            {
                os << "%FLAG AMOEBA_FORCEFIELD" << std::endl;
                os << "%COMMENT This indicates that this parm file is specific to amoeba" << std::endl;
                os << "%COMMENT This must be present if do_amoeba(in mdin) is 1" << std::endl;
                os << "%COMMENT This must NOT be present if do_amoeba is 0" << std::endl;
                os << "%FORMAT(i5)" << std::endl;
                os << "    1" <<  std::endl;
            }
        
            void write_pole_type( ostream& os, const molecule_t& mol )
            {
                os << "%FLAG AMOEBA_ATOM_TYPE_INDEX" << std::endl;
                os << "%COMMENT   dimention = (" << mol.natom() << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;

                fortran_t format( os, "10I8" );
                atomiter_t ai = mol.atom_begin();
                atomiter_t ae = mol.atom_end();
                for( ; ai != ae; ++ai )
                {
                    fortran_write_iparm( &format, *ai, POLTYPE );
                }
                format.end();
            }
        
            void write_atom_element( ostream& os, const molecule_t& mol )
            {
                os << "%FLAG AMOEBA_ATOMIC_NUMBER" << std::endl;
                os << "%COMMENT   dimention = (" << mol.natom() << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;

                fortran_t format( os, "10I8" );
                atomiter_t ai = mol.atom_begin();
                atomiter_t ae = mol.atom_end();
                for( ; ai != ae; ++ai )
                {
                    fortran_write_iparm( &format, *ai, ELEMENT );
                }
                format.end();
            }

            void write_atom_type( ostream& os, const molecule_t& mol )
            {
                os << "%FLAG AMOEBA_ATOM_CLASS_INDEX" << std::endl;
                os << "%COMMENT   dimention = (" << mol.natom() << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;

                fortran_t format( os, "10I8" );
                atomiter_t ai = mol.atom_begin();
                atomiter_t ae = mol.atom_end();
                for( ; ai != ae; ++ai )
                {
                    fortran_write_sparm( &format, *ai, TYPE );
                }
                format.end();
            }

            struct all_true
            {
                template< typename T >
                bool operator()( const T& )
                {
                    return true;
                }
            };

            struct list_atom_and_type
            {
	        list_atom_and_type(const hashid_t& idparm, const hashid_t& typeparm=TYPEID )
		{
		    m_idparm = idparm;
		    m_typeparm = typeparm;
		}

                void operator()( const morf_t p, vector<int>& ids )
                {
                    assert( ids.size() == 0 );
            
                    atomiter_t atom = p.atom_begin();
                    for( ; atom != p.atom_end(); ++atom )
                    {
                        ids.push_back( atom->get_i(m_idparm) );
                    }

                    ids.push_back( p.get_i(m_typeparm) );
                }

		int m_idparm;
		int m_typeparm;
            };

            void write_bond( ostream& os, const molecule_t& mol )
            {
                os << "%FLAG AMOEBA_REGULAR_BOND_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw( 8 ) << mol.nbond() << std::endl;
        
                os << "%FLAG AMOEBA_REGULAR_BOND_LIST" << std::endl;
                os << "%COMMENT dimension = (3," << mol.nbond() << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;
        
                fortran_t format( os, "10I8" );
                write_object( format, mol.bond_begin(), mol.bond_end(), all_true(), list_atom_and_type(ID) );
                format.end();
            }

            void write_bond_parm( ostream& os, const molecule_t& ff, const vector< vector<double> >& bondparm )
            {
                os << "%FLAG AMOEBA_REGULAR_BOND_NUM_PARAMS" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << bondparm[0].size() << std::endl;

                write_ffparm( os, bondparm, 2, "AMOEBA_REGULAR_BOND" );
            
                os << "%FLAG AMOEBA_REGULAR_BOND_FTAB_DEGREE" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << 4 << std::endl;
        
                os << "%FLAG AMOEBA_REGULAR_BOND_FTAB_COEFFS" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;

                fortran_t format( os, "5E16.8" );
                format( 0.0 );
                format( 0.0 );
                format( 1.0 );
                format( ff.get_d("bond-cubic") );
                format( ff.get_d("bond-quartic") );
                format.end();
            }

            struct has_urey
            {
                bool operator()( const morf_t& angl )
                {
                    numvec urey = angl.get_v(UREY);
                    return urey[0] != 0.0;
                }
            };
        
            struct non_urey
            {
                bool operator()( const morf_t& angl )
                {
                    numvec urey = angl.get_v(UREY);
                    return urey[0] == 0.0;
                }
            };

            struct angl_id13
            {
                void operator()( const morf_t& angl, vector< int >& ids )
                {
                    ids.clear();
                    ids.push_back( atom_1st(angl).get_i(ID) );
                    ids.push_back( atom_3rd(angl).get_i(ID) );
                    ids.push_back( angl.get_i(UREYID) );
                }
            };

            void write_urey( ostream& os, const molecule_t& mol )
            {
                int nurey = count_if( mol.angl_begin(), mol.angl_end(), has_urey() );
        
                os << "%FLAG AMOEBA_UREY_BRADLEY_BOND_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << nurey << std::endl;
        
                os << "%FLAG AMOEBA_UREY_BRADLEY_BOND_LIST" << std::endl;
                os << "%COMMENT   dimension = (3," << nurey << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;

                fortran_t format( os, "10I8" );
                write_object( format, mol.angl_begin(), mol.angl_end(), has_urey(), angl_id13() );
                format.end();
            }

            void write_urey_parm( ostream& os, const vector< vector<double> >& ureyparm )
            {
                os << "%FLAG AMOEBA_UREY_BRADLEY_BOND_NUM_PARAMS" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << ureyparm[0].size() << std::endl;

                write_ffparm( os, ureyparm, 2, "AMOEBA_UREY_BRADLEY_BOND" );
        
                os << "%FLAG AMOEBA_UREY_BRADLEY_BOND_FTAB_DEGREE" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << 2 << std::endl;
        
                os << "%FLAG AMOEBA_UREY_BRADLEY_BOND_FTAB_COEFFS" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;
                fortran_t format( os, "5E16.8" );
                format( 0.0 );
                format( 0.0 );
                format( 1.0 );
                format.end();
            }

            struct reg_angl
            {
                bool operator()( const morf_t& angl ) const
                {
                    return angl.noops() == 0;
                }
            };

            void write_reg_angl( ostream& os, const molecule_t& mol )
            {
                int nreg = count_if( mol.angl_begin(), mol.angl_end(), reg_angl() );

                os << "%FLAG AMOEBA_REGULAR_ANGLE_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw( 8 ) << nreg << std::endl;
        
                os << "%FLAG AMOEBA_REGULAR_ANGLE_LIST" << std::endl;
                os << "%COMMENT dimension = (4," << nreg << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;
        
                fortran_t format( os, "10I8" );
                write_object( format, mol.angl_begin(), mol.angl_end(), reg_angl(), list_atom_and_type(ID) );
		format.end();
            }
        
            struct trig_id
            {
                void operator()( const morf_t& oops, vector< int >& ids ) const
                {
                    ids.push_back( atom_1st(oops).get_i(ID) );
                    ids.push_back( atom_3rd(oops).get_i(ID) );
                    ids.push_back( atom_2nd(oops).get_i(ID) );
                    ids.push_back( atom_4th(oops).get_i(ID) );
                    ids.push_back( angl_t::get( atom_1st(oops), atom_3rd(oops), atom_2nd(oops) ).get_i(TYPEID) );
                }
            };

            void write_tri_angl( ostream& os, const molecule_t& mol )
            {
                os << "%FLAG AMOEBA_TRIGONAL_ANGLE_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw( 8 ) << mol.noops() << std::endl;
        
                os << "%FLAG AMOEBA_TRIGONAL_ANGLE_LIST" << std::endl;
                os << "%COMMENT dimension = (5," << mol.noops() << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;
        
                fortran_t format( os, "10I8" );
                write_object( format, mol.impr_begin(), mol.impr_end(), all_true(), trig_id() );
                format.end();
            }
 
            void write_angl_parm( ostream& os, const molecule_t& ff, const vector< vector< double > >& anglparm, const string& prefix )
            {
                os << "%FLAG AMOEBA_" << prefix  << "_ANGLE_NUM_PARAMS" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << anglparm[0].size() << std::endl;
 
                write_ffparm( os, anglparm, 2, "AMOEBA_" + prefix +"_ANGLE" );

                os << "%FLAG AMOEBA_" << prefix << "_ANGLE_FTAB_DEGREE" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << 6 << std::endl;
        
                os << "%FLAG AMOEBA_" << prefix << "_ANGLE_FTAB_COEFFS" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;
                fortran_t format( os, "5E16.8" );
                format( 0 );
                format( 0 );
                format( 1 );
                format( ff.get_d("angle-cubic")  );
                format( ff.get_d("angle-quartic") );
                format( ff.get_d("angle-pentic") );
                format( ff.get_d("angle-sextic") );
                format.end();
            }

            struct oops_id
            {
                 void operator()( const morf_t& oops, vector<int>& parm ) const
                 {
                     parm.clear();
                     parm.push_back( atom_1st(oops).get_i(ID) );
                     parm.push_back( atom_3rd(oops).get_i(ID) );
                     parm.push_back( atom_2nd(oops).get_i(ID) );
                     parm.push_back( atom_4th(oops).get_i(ID) );
                     parm.push_back( oops.get_i(TYPEID) );
                 }
            };

            void write_oops( ostream& os, const molecule_t& mol )
            {
                os << "%FLAG AMOEBA_OPBEND_ANGLE_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw( 8 ) << mol.noops() << std::endl;
        
                os << "%FLAG AMOEBA_OPBEND_ANGLE_LIST" << std::endl;
                os << "%COMMENT dimension = (5," << mol.noops() << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;
            
                fortran_t format( os, "10I8" );
                write_object( format, mol.impr_begin(), mol.impr_end(), all_true(), oops_id() );
                format.end();
            }    
        
            void write_oops_parm( ostream& os, const molecule_t& ff, const vector< vector< double > >& params )
            {
                os << "%FLAG AMOEBA_OPBEND_ANGLE_NUM_PARAMS" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << params[0].size() << std::endl;

                write_ffparm( os, params, 1, "AMOEBA_OPBEND_ANGLE" );
            
                os << "%FLAG AMOEBA_OPBEND_ANGLE_FTAB_DEGREE" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << 6 << std::endl;
        
                os << "%FLAG AMOEBA_OPBEND_ANGLE_FTAB_COEFFS" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;
                fortran_t format( os, "5E16.8" );
                format( 0 );
                format( 0 );
                format( 1 );
                format( ff.get_d("angle-cubic") );
                format( ff.get_d("angle-quartic") );
                format( ff.get_d("angle-pentic") );
                format( ff.get_d("angle-sextic") );
                format.end();
            }
        
            struct force_not_zero
            {
                bool operator()( const morf_t& tors ) const
                {
                    double force = tors.get_d(FORCE);
                    return force > 1e-6 || force < -1e-6;
                }
            };

            void write_tors( ostream& os, const molecule_t& mol )
            {
                int ntors = count_if( mol.dihe_begin(), mol.dihe_end(), force_not_zero() );

                os << "%FLAG AMOEBA_TORSION_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << ntors << std::endl;
        
                os << "%FLAG AMOEBA_TORSION_LIST" << std::endl;
                os << "%COMMENT   dimension = (5," << ntors << std::endl;
                os << "%FORMAT(10I8)" << std::endl;
        
                fortran_t format( os, "10I8" );
                write_object( format, mol.dihe_begin(), mol.dihe_end(), force_not_zero(), list_atom_and_type(ID) );
                format.end();
            }
        
            void write_tors_parm( ostream& os, const vector< vector< double > >& torsparm )
            {
                os << "%FLAG AMOEBA_TORSION_NUM_PARAMS" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << torsparm[0].size() << std::endl;
            
                write_ffparm( os, torsparm, 3, "AMOEBA_TORSION" );
            }
        
            void write_pitors( ostream& os, const molecule_t& mol )
            {
                os << "%FLAG AMOEBA_PI_TORSION_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << mol.nptor() << std::endl;
        
                os << "%FLAG AMOEBA_PI_TORSION_LIST" << std::endl;
                os << "%COMMENT   dimension = (7," << mol.nptor() << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;
        
                fortran_t format( os, "10I8" );
                write_object( format, mol.ptor_begin(), mol.ptor_end(), all_true(), list_atom_and_type(ID) );   
                format.end();
            }
        
            void write_pitors_parm( ostream& os, const vector< vector< double > >& ptorparm )
            {
                os << "%FLAG AMOEBA_PI_TORSION_NUM_PARAMS" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << ptorparm[0].size() << std::endl;
            
                write_ffparm( os, ptorparm, 3, "AMOEBA_PI_TORSION" );
            }

            struct has_strbnd
            {
                bool operator()( const morf_t& angl )
                {
                    numvec strbnd = angl.get_v(STRBND);
                    return strbnd[0]!=0.0;
                }
            };

            struct non_strbnd
            {
                bool operator()( const morf_t& angl )
                {
                    numvec strbnd = angl.get_v(STRBND);
                    return strbnd[0] == 0.0;
                }
            };
        
            void write_strbnd( ostream& os, const molecule_t& mol )
            {
                int nstrbnd = count_if( mol.angl_begin(), mol.angl_end(), has_strbnd() );
        
                os << "%FLAG AMOEBA_STRETCH_BEND_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << nstrbnd << std::endl;
            
                os << "%FLAG AMOEBA_STRETCH_BEND_LIST" << std::endl;
                os << "%COMMENT   dimension = (4," << nstrbnd << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;
        
                fortran_t format( os, "10I8" );
                write_object( format, mol.angl_begin(), mol.angl_end(), has_strbnd(), list_atom_and_type(ID, STRBNDID) );
                format.end();
            }
        
            void write_strbnd_parm( ostream& os, const vector< vector< double > >& strbndparm )
            {
                os << "%FLAG AMOEBA_STRETCH_BEND_NUM_PARAMS" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << strbndparm[0].size() << std::endl;
                write_ffparm( os, strbndparm, 4, "AMOEBA_STRETCH_BEND" );
            }
        
            void write_tor2( ostream& os, const molecule_t& mol )
            {
                os << "%FLAG AMOEBA_TORSION_TORSION_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << mol.ntor2() << std::endl;

                os << "%FLAG AMOEBA_TORSION_TORSION_LIST" << std::endl;
                os << "%COMMENT   dimension = (6," << 147 << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;
        
                fortran_t format( os, "10I8" );
                write_object( format, mol.tor2_begin(), mol.tor2_end(), all_true(), list_atom_and_type(ID) );
                format.end();
            }
        
                
            void write_tor2_parm( ostream& os, const molecule_t& atomff )
            {
                os << "%FLAG AMOEBA_TORSION_TORSION_NUM_PARAMS" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << atomff.ntor2() << std::endl;
                
                string prefix = "AMOEBA_TORSION_TORSION_TORTOR_TABLE_0";
                
                mobjiter_t parm = atomff.tor2_begin();
                for( ; parm != atomff.tor2_end(); ++parm )
                {
                    os << "%FLAG " << prefix << parm->absid() + 1 << "_DIMS" << std::endl;
                    os << "%COMMENT dimension = (2)" << std::endl;
                    os << "%FORMAT(2I8)" << std::endl;
                    os << std::setw(8) << 25;
                    os << std::setw(8) << 25 << std::endl;
                    
                    os << "%FLAG " << prefix << parm->absid() + 1 << "_ANGLE1" << std::endl;
                    os << "%COMMENT   dimension = (25)" << std::endl;
                    os << "%FORMAT(5E16.8)" << std::endl;
                    
                    fortran_t format( os, "5E16.8" );

                    double start = -180.0;
                    for( int j=0; j < 25; j++ )
                    {
                        format( start + j*15.0 );
                    }
                    
                    format.end();
                    
                    os << "%FLAG " << prefix << parm->absid() + 1 << "_ANGLE2" << std::endl;
                    os << "%COMMENT   dimension = (25)" << std::endl;
                    os << "%FORMAT(5E16.8)" << std::endl;

                    format.begin();
                    
                    for( int j=0; j < 25; j++ )
                    {
                        format( start + j*15.0 );
                    }
                    
                    format.end();

                    typedef shared_ptr< vector<double> > dvec_ptr;

                    dvec_ptr func = any_cast< dvec_ptr >( parm->get_a(TOR2FUNC) );                    

                    assert( func->size() == 3750 );

                    os << "%FLAG " << prefix << parm->absid() + 1 << "_FUNC" << std::endl;
                    os << "%COMMENT   dimension = (25,25)" << std::endl;
                    os << "%FORMAT(5E16.8)" << std::endl;

                    format.begin();                    

                    for( int j=0; j < 625; ++j )
                    {
                        format( func->at(6*j+2) );
                    }
                    format.end();
                    
                    os << "%FLAG " << prefix << parm->absid() + 1 << "_DFUNC_DANGLE1" << std::endl;
                    os << "%COMMENT   dimension = (25,25)" << std::endl;
                    os << "%FORMAT(5E16.8)" << std::endl;

                    format.begin();                    

                    for( int j=0; j < 625; ++j )
                    {
                        format( func->at(6*j+3) );
                    }
                    
                    format.end();

                    os << "%FLAG " << prefix << parm->absid() + 1 << "_DFUNC_DANGLE2" << std::endl;
                    os << "%COMMENT   dimension = (25,25)" << std::endl;
                    os << "%FORMAT(5E16.8)" << std::endl;

                    format.begin();                    

                    for( int j=0; j < 625; ++j )
                    {
                        format( func->at(6*j+4) );
                    }
                    
                    format.end();
                    
                    os << "%FLAG " << prefix << parm->absid() + 1 << "_D2FUNC_DANGLE1_DANGLE2" << std::endl;
                    os << "%COMMENT   dimension = (25,25)" << std::endl;
                    os << "%FORMAT(5E16.8)" << std::endl;

                    format.begin();                    

                    for( int j=0; j < 625; ++j )
                    {
                        format( func->at(6*j+5) );
                    }
                    
                    format.end();
                }
            }
        
            void write_vdw_type( ostream& os, const molecule_t& mol )
            {
                os << "%FLAG AMOEBA_VDW_ATOM_TYPES_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << mol.natom() << std::endl;
            
                os << "%FLAG AMOEBA_VDW_ATOM_TYPES_LIST" << std::endl;
                os << "%COMMENT   dimension = (1," << mol.natom() << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;
            
                fortran_t format( os, "10I8" );

                atomiter_t ai = mol.atom_begin();
                atomiter_t ae = mol.atom_end();
                for( ; ai != ae; ++ai )
                {
                    fortran_write_iparm( &format, *ai, TYPEID );
                }

                format.end();
            }
        
            void write_vdw_parent( ostream& os, const molecule_t& mol )
            {
                os << "%FLAG AMOEBA_VDW_ATOM_PARENT_LIST" << std::endl;
                os << "%COMMENT   dimension = (1," << mol.natom() << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;
        
                fortran_t i_format( os, "10I8" );
                atomiter_t atom = mol.atom_begin();
                for( ; atom != mol.atom_end(); ++atom )
                {
                    if( atom->get_i(ELEMENT) == HYDROGEN )
                    {
                        i_format( atom->atom_begin()->get_i(ID) );
                    }
                    else
                    {
                        i_format( atom->get_i(ID) );
                    }
                }

                i_format.end();
               
                os << "%FLAG AMOEBA_VDW_PARENT_COORD_WEIGHT_LIST" << std::endl;
                os << "%COMMENT   dimension = (1," << mol.natom() << ")" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;

                fortran_t d_format( os, "5E16.8" );

                atomiter_t ai = mol.atom_begin();
                atomiter_t ae = mol.atom_end();
                for( ; ai != ae; ++ai )
                {
                    fortran_write_dparm( &d_format, *ai, VBUFF );
                }

                d_format.end();
            }

        
            void  write_vdw_parm( ostream& os, const vector< vector<double> >& vdwparm )
            {
                assert( vdwparm[0].size() == vdwparm[1].size() );
            
                vector< double > rstars;
                vector< double > depths;

                int nvdw = vdwparm[0].size();
            
                for( int i=0; i < nvdw; ++i )
                {
                    for( int j=0; j < nvdw; ++j )
                    {
                        double ri = vdwparm[1][i];
                        double rj = vdwparm[1][j];
                        rstars.push_back( (ri*ri*ri + rj*rj*rj)/(ri*ri + rj*rj) );

                        double ei = vdwparm[0][i];
                        double ej = vdwparm[0][j];
                        depths.push_back( (4*ei*ej)/pow( sqrt(ei)+sqrt(ej), 2 ) );
                    }
                }
            
                os << "%FLAG AMOEBA_VDW_BUFFER_DELTA" << std::endl;
                os << "%FORMAT(E16.8)" << std::endl;
                os << "  0.70000000E-01" << std::endl;

                os << "%FLAG AMOEBA_VDW_BUFFER_GAMMA" << std::endl;
                os << "%FORMAT(E16.8)" << std::endl;
                os << "  0.12000000E+00" << std::endl;
        
                os << "%FLAG AMOEBA_VDW_PARAMS_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << nvdw << std::endl;
        
                os << "%FLAG AMOEBA_VDW_MIXED_RADII_LIST" << std::endl;
                os << "%COMMENT   dimension = (" << nvdw << "," << nvdw << ")" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;

                fortran_t format( os, "5E16.8" );

                for( int i=0; i < (int)rstars.size(); ++i )
                {
                    format( rstars[i] );
                }            

                format.end();
        
                os << "%FLAG AMOEBA_VDW_MIXED_EPSILONS_LIST" << std::endl;
                os << "%COMMENT   dimension = (" << nvdw << "," << nvdw << ")" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;

                format.begin();

                for( int i=0; i < (int)depths.size(); ++i )
                {
                    format( depths[i] );
                }

                format.end();            
            }

            void write_mpole( ostream& os, const molecule_t& mol )
            {
                os << "%FLAG AMOEBA_LOCAL_FRAME_MULTIPOLES_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << mol.natom() << std::endl;
        
                os << "%FLAG AMOEBA_LOCAL_FRAME_MULTIPOLES_LIST" << std::endl;
                os << "%COMMENT dimension = (10," << mol.natom() << ")" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;
        
                fortran_t format( os, "5E16.8" );

                atomiter_t atom = mol.atom_begin();
                for( ; atom != mol.atom_end(); ++atom )
                {
                    numvec pole = atom->get_v(POLE);
                    
                    format( atom->get_d(PCHG) );
                    
                    format( pole[0] );
                    format( pole[1] );
                    format( pole[2] );
                    format( pole[3] * 0.5 );
                    format( pole[5] * 0.5 );
                    format( pole[8] * 0.5 );
                    format( pole[4] );
                    format( pole[6] );
                    format( pole[7] );
                }
            }

            bool create_chirial_frame(morf_t& atom, const morf_t& pole, vector< vector<int> >& regulars, vector< vector<int> >& chirials )
            {
                impriter_t j = pole.impr_begin();
		
                for( ; j != pole.impr_end(); ++j )
                {
		    int poltypeid_1 = atom_1st(*j).get_i(ID);
		    int poltypeid_2 = atom_2nd(*j).get_i(ID);
		    int poltypeid_4 = atom_4th(*j).get_i(ID);

                    atomiter_t atm1 = find_if( atom.atom_begin(), atom.atom_end(), iparm_cmper1(POLTYPEID, poltypeid_1 ) );
                    if( atm1 == atom.atom_end() )
                    {
                        continue;
                    }
                    
                    atomiter_t atm2 = find_if( atom.atom_begin(), atom.atom_end(), iparm_cmper1(POLTYPEID, poltypeid_2) );
                    if( atm2 == atom.atom_end() )
                    {
                        continue;
                    }
                    
                    atomiter_t atm4 = find_if( atom.atom_begin(), atom.atom_end(), iparm_cmper1(POLTYPEID, poltypeid_4) );
                    if( atm4 == atom.atom_end() )
                    {
                        continue;
                    }
                    
                    atom.set_d(PCHG, j->get_d(PCHG) );
                    atom.set_v(POLE, j->get_v(POLE) );

                    vector<int> regular(5);

		    regular[0] = atom.get_i(ID);
		    regular[1] = 1;
		    regular[2] = atom.get_i(ID);
		    regular[3] = atm1->get_i(ID);
		    regular[4] = 1;


                    regulars.push_back( regular );
                                
        	    regular[0] = atom.get_i(ID);
		    regular[1] = 2;
		    regular[2] = atom.get_i(ID);
		    regular[3] = atm2->get_i(ID);
		    regular[4] = 1;

                    regulars.push_back( regular );
                    
		    vector<int> chirial(3);
                    chirial[0] = atom.get_i(ID);
		    chirial[1] = atm4->get_i(ID);
		    chirial[2] = 1;

                    chirials.push_back( chirial );
                    
                    return true;
                }

                return false;
            }
            
            bool create_regular_frame( morf_t& atom, const morf_t& pole, vector< vector<int> >& regulars )
            {
                angliter_t j = pole.angl_begin();
                for( ; j != pole.angl_end(); ++j )
                {
                    int poltypeid_1 = atom_1st(*j).get_i(ID);
		    int poltypeid_3 = atom_3rd(*j).get_i(ID);

                    atomiter_t atm1 = find_if( atom.atom_begin(), atom.atom_end(), iparm_cmper1(POLTYPEID, poltypeid_1) );
                    if( atm1 == atom.atom_end() )
                    {
                        continue;
                    }

                    atomiter_t begin= (poltypeid_1==poltypeid_3) ? atm1+1 : atom.atom_begin();
                    atomiter_t atm3 = find_if( begin, atom.atom_end(), iparm_cmper1(POLTYPEID, poltypeid_3) );
                    
                    if( atm3 == atom.atom_end() )
                    {
                        atm3 = find_if( atm1->atom_begin(), atm1->atom_end(), iparm_cmper1(POLTYPEID, poltypeid_3) );
						if( atm3 == atm1->atom_end() )
                        {
                            continue;
                        }
                        if( *atm3 == atom )
                        {
                            atm3 = find_if( atm3+1, atm1->atom_end(), iparm_cmper1(POLTYPEID, poltypeid_3) );
						}     
                    }
                    
                    atom.set_d(PCHG, j->get_d(PCHG) );
                    atom.set_v(POLE, j->get_v(POLE) );

                    if( j->get_i(TYPEID) == (int)BISECTOR )
                    {
		        vector<int> regular(5);
                        regular[0] = atom.get_i(ID);
			regular[1] = 1;
                        regular[2] = atom.get_i(ID);
			regular[3] = atm1->get_i(ID);
                        regular[4] = 2; 
                        regulars.push_back( regular );

                        regular[0] = atom.get_i(ID);
			regular[1] = 1;
			regular[2] = atom.get_i(ID);
			regular[3] = atm3->get_i(ID);
                        regular[4] = 2; 
                        regulars.push_back( regular );

                        regular[0] = atom.get_i(ID);
			regular[1] = 2;
                        regular[2] = atom.get_i(ID);
			regular[3] = atm3->get_i(ID),
                        regular[4] = 1; 
                        regulars.push_back( regular );
                    }
                    else
                    {
                        vector<int> regular(5);
                        regular[0] = atom.get_i(ID);
			regular[1] = 1;
                        regular[2] = atom.get_i(ID);
			regular[3] = atm1->get_i(ID);
                        regular[4] = 1; 
                        regulars.push_back( regular );

                        regular[0] = atom.get_i(ID);
			regular[1] = 2;
                        regular[2] = atom.get_i(ID);
			regular[3] = atm3->get_i(ID);
                        regular[4] = 1; 
                        regulars.push_back( regular );
                    }

                    return true;
                }

                return false;
            }

            morf_t find_pole( const molecule_t& poleff, const morf_t& atom )
            {
                string polt = atom.get_s(POLTYPE);
                return morf_t( poleff, ATOM, atoi(polt.c_str())-1 );
            }
            

            void create_atomic_frame( morf_t& atom, const molecule_t& poleff, vector< vector<int> >& chirials, vector< vector<int> >& regulars )
            {
                morf_t pole = find_pole( poleff, atom );
                
                if( pole.nangl() > 0 )
                {
                    if( ! create_regular_frame( atom, pole, regulars ) )
                    {
                        throw std::runtime_error( "Error: can not create regular frame for atom " + atom.get_s(POLTYPE) );
                    }
                }
                else if( pole.noops() > 0 )
                {
                    if( ! create_chirial_frame( atom, pole, regulars, chirials ) )
                    {
                        throw std::runtime_error( "Error: can not create chirial frame for atom " + atom.get_s(POLTYPE) );
                    }
                }
            }

            void create_frame( molecule_t& mol, const molecule_t& poleff, vector< vector<int> >& chirials, vector< vector<int> >& regulars )
            {
                atomiter_t atom = mol.atom_begin();
                for( ; atom != mol.atom_end(); ++atom )
                {
                    create_atomic_frame( *atom, poleff, chirials, regulars );
                }
            }

            void write_chirial_frame( ostream& os, const vector< vector<int> >& chirials )
            {
                os << "%FLAG AMOEBA_CHIRAL_FRAME_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << chirials.size() << std::endl;
                
                os << "%FLAG AMOEBA_CHIRAL_FRAME_LIST" << std::endl;
                os << "%COMMENT   dimension = (3," << chirials.size() << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;

                fortran_t format( os, "10I8" );
                for( int i=0; i < (int)chirials.size(); ++i )
                {
                    format( chirials[i][0] );
                    format( chirials[i][1] );
                    format( chirials[i][2] );
                }
                
                format.end();
            }
                
            void write_regular_frame( ostream& os, const vector< vector<int> >& regulars )
            {
                os << "%FLAG AMOEBA_FRAME_DEF_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << regulars.size() << std::endl;
        
                os << "%FLAG AMOEBA_FRAME_DEF_LIST" << std::endl;
                os << "%COMMENT   dimension = (5," << regulars.size() << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;

                fortran_t format( os, "10I8" );
                
                for( int i=0; i < (int)regulars.size(); ++i )
                {
                    format( regulars[i][0] );
                    format( regulars[i][1] );
                    format( regulars[i][2] );
                    format( regulars[i][3] );
                    format( regulars[i][4] );
                }

                format.end();
            }

            void write_frame( std::ostream& os, const vector< vector<int> >& chirials, const vector< vector<int> >& regulars )
            {
                write_chirial_frame( os, chirials );
                
                write_regular_frame( os, regulars );
            }


            void traverse( const morf_t& root, const morf_t& curt, atomvec_t& group, const molecule_t& ffpol )
	    {
                group.push_back( curt );
            
	        string curt_name = curt.get_s(POLTYPE);

                atom_t pol_curt = atom_t::get( ffpol, curt_name );
        
                atomiter_t next = curt.atom_begin();
        
                for( ; next != curt.atom_end(); ++next )
                {
                    if( *next == root )
                    {
                        continue;
                    }
            
                    if( group.count( *next ) > 0 )
                    {
                        continue;
                    }
            
                    atom_t pol_next = atom_t::get( ffpol, next->get_s(POLTYPE) );
                    
                    if( pol_curt.is_connected_to(pol_next) )
                    {
                        traverse( curt, *next, group, ffpol );
                    }
                }
            }
    
            void divide_pgroup( const molecule_t& mol, vector< atomvec_t >& pgroups, const molecule_t& ffpol )
            {
                atomiter_t atom = mol.atom_begin();
        
                for( ; atom != mol.atom_end(); ++atom )
                {
                    bool newgroup = true;
            
                    for( int i=0; i < (int)pgroups.size(); ++i )
                    {
                        if( pgroups[i].count( *atom ) > 0 )
                        {
                            newgroup = false;
                        }
                    }
            
                    if( newgroup )
                    {
                        pgroups.push_back( atomvec_t() );
                
                        traverse( morf_t(mol, ATOM, -1), *atom, pgroups.back(), ffpol );
                    }   
                }   
            }

            void setup_atom_adjust( int atom, atomvec_t& group, vector< int >& list, vector< int >& dist )
            {
                std::set< int > ids;
                    
                for( int i=0; i< (int)group.size(); ++i )
                {
                    if( group[i].get_i(ID) > atom )
                    {
                        ids.insert( group[i].get_i(ID) );
                    }
                }

                for( int i=0; i < (int)list.size(); ++i )
                {
                    ids.insert( list[i] );
                }

                vector< int > mylist;
                vector< int > mydist;

                std::set< int >::iterator i = ids.begin();
                for( ; i != ids.end(); ++i )
                {
                    mylist.push_back( *i );

                    int pos = std::find( list.begin(), list.end(), *i ) - list.begin();                        

                    if( pos == (int)list.size() )
                    {
                        mydist.push_back( 9 );
                    }
                    else if( group.count( morf_t( group[0].getmol(), ATOM, *i-1 ) ) > 0 )
                    {
                        mydist.push_back( dist[ pos ] + 4 );
                    }
                    else
                    {
                        mydist.push_back( dist[ pos ] );
                    }                    
                }

                list.swap( mylist );
                dist.swap( mydist );
            }

            void setup_adjust( const molecule_t& mol, excl_t& excl, const molecule_t& poleff )
            {
                vector< atomvec_t > pgroups;
                divide_pgroup( mol, pgroups, poleff );
        
                atomiter_t atom = mol.atom_begin();
                for( ; atom != mol.atom_end(); ++atom )
                {
                    int aid = atom->get_i(ID);
            
                    vector<atomvec_t>::iterator pgi = pgroups.begin();
                    vector<atomvec_t>::iterator pge = pgroups.end();

                    for( ; pgi != pge && pgi->count(*atom)==0; ++pgi ); 
                    
                    assert( pgi != pge );

                    setup_atom_adjust( aid, *pgi, excl.list[ aid-1 ], excl.dist[ aid-1 ] );
                }
            }

            void write_adjust( ostream& os, const molecule_t& , excl_t& excl, const molecule_t& )
            {
                os << "%FLAG AMOEBA_ADJUST_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << excl.full_size() << std::endl;
            
                os << "%FLAG AMOEBA_ADJUST_LIST" << std::endl;
                os << "%COMMENT dimension = (3," <<  excl.full_size() << ")" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;

                fortran_t format( os, "10I8" );

                for( int i=0; i < (int)excl.list.size(); ++i )
                {
                    assert( excl.list[i].size() == excl.dist[i].size() );
                
                    for( int j=0; j < (int)excl.list[i].size(); ++j )
                    {
                        format( i+1);
                        format( excl.list[i][j] );
                        format( excl.dist[i][j] );
                    }
                }

                format.end();

            }
    
            void write_adjust_mask( ostream& os )
            {
                os << "%FLAG AMOEBA_ADJUST_VDW_WEIGHTS_LIST" << std::endl;
                os << "%COMMENT   dimension = (9)" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;

                fortran_t format( os, "5E16.8" );
                format( 0.0 ); format( 0.0 ); format( 1.0 ); 
                format( 1.0 ); format( 0.0 ); format( 0.0 );
                format( 1.0 ); format( 1.0 ); format( 1.0 );
                format.end();

                os << "%FLAG AMOEBA_ADJUST_MPOLE_WEIGHTS_LIST" << std::endl;
                os << "%COMMENT   dimension = (9)" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;

                format.begin();
                format( 0.0 ); format( 0.0 ); format( 0.4 ); 
                format( 0.8 ); format( 0.0 ); format( 0.0 );
                format( 0.4 ); format( 0.8 ); format( 1.0 );
                format.end();

                os << "%FLAG AMOEBA_ADJUST_DIRECT_WEIGHTS_LIST" << std::endl;
                os << "%COMMENT   dimension = (9)" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;

                format.begin();
                format( 1.0 ); format( 1.0 ); format( 1.0 ); 
                format( 1.0 ); format( 0.0 ); format( 0.0 );
                format( 0.0 ); format( 0.0 ); format( 0.0 );
                format.end();

                os << "%FLAG AMOEBA_ADJUST_POLAR_WEIGHTS_LIST" << std::endl;
                os << "%COMMENT   dimension = (9)" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;

                format.begin();
                format( 0.0 ); format( 0.0 ); format( 1.0 ); 
                format( 1.0 ); format( 0.0 ); format( 0.0 );
                format( 0.5 ); format( 1.0 ); format( 1.0 );
                format.end();

                os << "%FLAG AMOEBA_ADJUST_MUTUAL_WEIGHTS_LIST" << std::endl;
                os << "%COMMENT   dimension = (9)" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;

                format.begin();
                format( 1.0 ); format( 1.0 ); format( 1.0 ); 
                format( 1.0 ); format( 1.0 ); format( 1.0 );
                format( 1.0 ); format( 1.0 ); format( 1.0 );
                format.end();
            }


            void write_polar( ostream& os, const molecule_t& mol, const molecule_t& poleff )
            {
                os << "%FLAG AMOEBA_POLARIZABILITY_NUM_LIST" << std::endl;
                os << "%FORMAT(I8)" << std::endl;
                os << std::setw(8) << mol.natom() << std::endl;

                os << "%FLAG AMOEBA_POLARIZABILITY_LIST" << std::endl;
                os << "%COMMENT   dimension = (" << mol.natom() << ")" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;
                
                fortran_t format( os, "5E16.8" );
        
                atomiter_t atom = mol.atom_begin();
                for( ; atom != mol.atom_end(); ++atom )
                {
                    morf_t pole = find_pole( poleff, *atom );
                    format( pole.get_d(POLAR) );
                }
                
                format.end();
            }

        
        } // namespace amoeba

    } // namespace prmtop
 

    using namespace prmtop;

    void write_amoeba_prmtop( ostream& os, molecule_t& mol, const molecule_t& aff, const molecule_t& pff )
    {
        energee_t e( mol );

        excl_t& excl = e.m_excl;
        exclude( mol, excl, 4 );
        amoeba::setup_adjust( mol, excl, pff );

        e.assignparm( aff, AMOEBA );
        parmset_t& params = e.m_parmset;
        write_prmtop( os, e.getnabparm() );
        
        amoeba::write_head( os );
        amoeba::write_pole_type( os, mol );
        amoeba::write_atom_element( os, mol );
        amoeba::write_atom_type( os, mol );
        amoeba::write_bond( os, mol );
        amoeba::write_bond_parm( os, aff, params.bond );
        amoeba::write_urey( os, mol );
        amoeba::write_urey_parm( os, params.urey );
        amoeba::write_reg_angl( os, mol );
        amoeba::write_angl_parm( os, aff, params.angl, "REGULAR" );
        amoeba::write_tri_angl( os, mol );
        amoeba::write_angl_parm( os, aff, params.angl, "TRIGONAL" );
        amoeba::write_oops( os, mol );
        amoeba::write_oops_parm( os, aff, params.oops );
        amoeba::write_tors( os, mol );
        amoeba::write_tors_parm( os, params.tors );
        amoeba::write_pitors( os, mol );
        amoeba::write_pitors_parm( os, params.ptor );
        amoeba::write_strbnd( os, mol );
        amoeba::write_strbnd_parm( os, params.strbnd );
        amoeba::write_tor2( os, mol );
        amoeba::write_tor2_parm( os, aff );
        amoeba::write_vdw_type( os, mol );
        amoeba::write_vdw_parent( os, mol );
        amoeba::write_vdw_parm( os, params.vdw );

        vector< vector<int> > chirials;
        vector< vector<int> > regulars;

        amoeba::create_frame( mol, pff, chirials, regulars );
        amoeba::write_mpole( os, mol );
        amoeba::write_frame( os, chirials, regulars );
        amoeba::write_adjust( os, mol, excl, pff );
        amoeba::write_adjust_mask( os );
        amoeba::write_polar( os, mol, pff );    
    } 
   
} // namespace mort

