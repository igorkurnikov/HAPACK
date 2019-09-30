#include <iomanip>
#include <common.hpp>
#include <object.hpp>
#include <capbox.hpp>
#include "prmtop.hpp"
#include "energee.hpp"

namespace mort
{
    using namespace std;

    namespace prmtop
    {
        void write_string( ostream& os, const string& flag, const string& str )
        {
            os << "%FLAG " << flag << endl;
            os << "%FORMAT(1a80)" << endl;
            os << str << endl;
        }

        void write_sarray( ostream& os, const string& flag, const string& str, int nitem )
        {
            os << "%FLAG " << flag << endl;
            os << "%FORMAT(20a4)" << endl;

            int ITEM_WIDTH = 4;
            int LINE_WIDTH = 80;
            int nline = ITEM_WIDTH*nitem/LINE_WIDTH;
            for( int i=0; i < nline; ++i )
            {
                int bgn = i*LINE_WIDTH;
                os << str.substr( bgn, LINE_WIDTH ) << std::endl;
            }

            int nleft = ITEM_WIDTH*nitem%LINE_WIDTH;
            if( nleft != 0 )
            {
                os << str.substr( nline*LINE_WIDTH, nleft );
                os << string( LINE_WIDTH-nleft, ' ' ) << std::endl;
            }
        }

        void write_iarray( ostream& os, const string& flag, const int* ia, int nitem )
        {
            os << "%FLAG " << flag << endl;
            os << "%FORMAT(10I8)" << endl;

            fortran_t format( os, "10I8" );
            for( int i=0; i < nitem; ++i )
            {
                format( ia[i] );
            }
            format.end();
        }

        void write_darray( ostream& os, const string& flag, const double* da, int nitem )
        {
            os << "%FLAG " << flag << endl;
            os << "%FORMAT(5E16.8)" << endl;

            fortran_t format( os, "5E16.8" );
            for( int i=0; i < nitem; ++i )
            {
                format( da[i] );
            }

            if( nitem==0 ) os << endl;
            format.end();
        }



        void write_head( ostream& os, const nabparm_t& prm )
        {
            // time_t result = time( NULL );
        
            os << "%VERSION VERSION_STAMP = V0001.000 ";
            os << "DATE = 05/22/06  12:10:21" << endl;
        
            os << "%FLAG TITLE" << endl;
            os << "%FORMAT(20a4)" << endl;
            os << endl;

            os << "%FLAG POINTERS" << endl;
            os << "%FORMAT(10I8)"  << endl;

            os << format( "%8d" ) % prm.Natom; // NTOTAT
            os << format( "%8d" ) % prm.Ntypes; // NTYPES
            os << format( "%8d" ) % prm.Nbonh;
            os << format( "%8d" ) % prm.Mbona;
            os << format( "%8d" ) % prm.Ntheth;
            os << format( "%8d" ) % prm.Mtheta;
            os << format( "%8d" ) % prm.Nphih;
            os << format( "%8d" ) % prm.Mphia;
            os << format( "%8d" ) % prm.Nhparm; // JHPARM
            os << format( "%8d" ) % prm.Nparm;
            os << endl; // JPARM
            
            os << format( "%8d" ) % prm.Nnb;
            os << format( "%8d" ) % prm.Nres;
            os << format( "%8d" ) % prm.Nbona; 
            os << format( "%8d" ) % prm.Ntheta; 
            os << format( "%8d" ) % prm.Nphia;
            os << format( "%8d" ) % prm.Numbnd;
            os << format( "%8d" ) % prm.Numang;
            os << format( "%8d" ) % prm.Nptra; 
            os << format( "%8d" ) % prm.Natyp;
            os << format( "%8d" ) % prm.Nphb; 
            os << endl; // ps.hbond.size( ) << endl;
        
            os << format( "%8d" ) % 0; // perturb, no longer supported
            os << format( "%8d" ) % 0;
            os << format( "%8d" ) % 0;
            os << format( "%8d" ) % 0;
            os << format( "%8d" ) % 0;
            os << format( "%8d" ) % 0;
            os << format( "%8d" ) % 0;
            os << format( "%8d" ) % prm.IfBox;
            os << format( "%8d" ) % prm.Nmxrs; 
            os << format( "%8d" ) % prm.IfCap;
            os << endl;
 
            os << format( "%8d" ) % prm.Numextra;
            //os << format( "%8d" ) % 0; //prm.lescopy;
            os << endl;
        }

        void write_bond( ostream& os, const nabparm_t& prm )
        {
            os << "%FLAG BONDS_INC_HYDROGEN" << endl;
            os << "%FORMAT(10I8)" << endl;

            fortran_t format( os, "10I8" );
            for( int i=0; i < prm.Nbonh; ++i )
            {
                format( prm.BondHAt1[i] );
                format( prm.BondHAt2[i] );
                format( prm.BondHNum[i] );
            }
            format.end();
    
            os << "%FLAG BONDS_WITHOUT_HYDROGEN" << endl;
            os << "%FORMAT(10I8)" << endl;

            format.begin();
            for( int i=0; i < prm.Nbona; ++i )
            {
                format( prm.BondAt1[i] );
                format( prm.BondAt2[i] );
                format( prm.BondNum[i] );
            }
            format.end();
        }

        void write_angl( ostream& os, const nabparm_t& prm )
        {
            os << "%FLAG ANGLES_INC_HYDROGEN" << endl;
            os << "%FORMAT(10I8)" << endl;

            fortran_t format( os, "10I8" );
            for( int i=0; i < prm.Ntheth; ++i )
            {
                format( prm.AngleHAt1[i] );
                format( prm.AngleHAt2[i] );
                format( prm.AngleHAt3[i] );
                format( prm.AngleHNum[i] );
            }
            format.end();
            
            os << "%FLAG ANGLES_WITHOUT_HYDROGEN" << endl;
            os << "%FORMAT(10I8)" << endl;

            format.begin();
            for( int i=0; i < prm.Ntheta; ++i )
            {
                format( prm.AngleAt1[i] );
                format( prm.AngleAt2[i] );
                format( prm.AngleAt3[i] );
                format( prm.AngleNum[i] );
            }
            format.end();
        }
    
        void write_tors( ostream& os, const nabparm_t& prm )
        {
            os << "%FLAG DIHEDRALS_INC_HYDROGEN" << endl;
            os << "%FORMAT(10I8)" << endl;

            fortran_t format( os, "10I8" );
            for( int i=0; i < prm.Nphih; ++i )
            {
                format( prm.DihHAt1[i] );
                format( prm.DihHAt2[i] );
                format( prm.DihHAt3[i] );
                format( prm.DihHAt4[i] );
                format( prm.DihHNum[i] );
            }
            format.end();

            format.begin();
            os << "%FLAG DIHEDRALS_WITHOUT_HYDROGEN" << endl;
            os << "%FORMAT(10I8)" << endl;
            format.begin();
            for( int i=0; i < prm.Nphia; ++i )
            {
                format( prm.DihAt1[i] );
                format( prm.DihAt2[i] );
                format( prm.DihAt3[i] );
                format( prm.DihAt4[i] );
                format( prm.DihNum[i] );
            }
            format.end();
        }
    
        static const char GBPARM_NAME[][80] = 
        {
            "Bondi radii (bondi)", 
            "amber6 modified Bondi radii (amber6)", 
            "modified Bondi radii (mbondi)",
            "H(N)-modified Bondi radii (mbondi2)"
        };

        void write_polar( ostream& os,  const molecule_t& m )
        {
            const vector<double> polar = get_dvec( m, ATOM, POLAR );
            write_darray( os, "POLARIZABILITY", &polar[0], m.natom() );
        }

        void write_solvent_inf( ostream& os, const nabparm_t& prm )
        {
            if( prm.IfBox )
            {
                os << "%FLAG SOLVENT_POINTERS" << std::endl;
                os << "%FORMAT(3I8)" << std::endl;
                os << std::setw( 8 ) << prm.Iptres;
                os << std::setw( 8 ) << prm.Nspm;
                os << std::setw( 8 ) << prm.Nspsol << std::endl;       
                write_iarray( os, "ATOMS_PER_MOLECULE", prm.Boundary, prm.Nspm );

                os << "%FLAG BOX_DIMENSIONS" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;
                fortran_t format( os, "5E16.8" );

                format( prm.Box[3] );
                format( prm.Box[0] );
                format( prm.Box[1] );
                format( prm.Box[2] );
                format.end();
            }
            else
            {
                assert( prm.IfCap );

                os << "%FLAG CAP_INFO" << std::endl;
                os << "%FORMAT(10I8)" << std::endl;
                os << std::setw(8) << prm.Natcap << std::endl;
            
                os << "%FLAG CAP_INFO2" << std::endl;
                os << "%FORMAT(5E16.8)" << std::endl;
                fortran_t format( os, "5E16.8" );
                format( prm.Cutcap );
                format( prm.Xcap );
                format( prm.Ycap );
                format( prm.Zcap );
                format.end();
            }        
        }
    

        void write_prmtop( ostream& os, const nabparm_t& prm )
        {
            write_head( os, prm );

            write_sarray( os, "ATOM_NAME", 		prm.AtomNames, 	prm.Natom );
            write_darray( os, "CHARGE",    		prm.Charges, 	prm.Natom );
            write_darray( os, "MASS",      		prm.Masses, 	prm.Natom );
            write_iarray( os, "ATOM_TYPE_INDEX", 	prm.Iac, 	prm.Natom );
            write_iarray( os, "NUMBER_EXCLUDED_ATOMS", 	prm.Iblo, 	prm.Natom );
            write_iarray( os, "NONBONDED_PARM_INDEX", 	prm.Cno,  	prm.Ntype2d );
            write_sarray( os, "RESIDUE_LABEL",        	prm.ResNames, 	prm.Nres );
            write_iarray( os, "RESIDUE_POINTER",      	prm.Ipres,    	prm.Nres );
            write_darray( os, "BOND_FORCE_CONSTANT",  	prm.Rk,       	prm.Numbnd );
            write_darray( os, "BOND_EQUIL_VALUE",     	prm.Req,      	prm.Numbnd );
            write_darray( os, "ANGLE_FORCE_CONSTANT", 	prm.Tk,       	prm.Numang );
            write_darray( os, "ANGLE_EQUIL_VALUE",    	prm.Teq,      	prm.Numang );
            write_darray( os, "DIHEDRAL_FORCE_CONSTANT",prm.Pk,       	prm.Nptra );
            write_darray( os, "DIHEDRAL_PERIODICITY", 	prm.Pn,       	prm.Nptra );
            write_darray( os, "DIHEDRAL_PHASE",      	prm.Phase,    	prm.Nptra );
            
            string tmp;
            if( mortenv().get_s("write14scale", tmp) && tmp=="on" )
            {
                write_darray( os, "SCEE_SCALE_FACTOR",  prm.Scee,	prm.Nptra );
		write_darray( os, "SCNB_SCALE_FACTOR",	prm.Scnb,	prm.Nptra );
            }


            write_darray( os, "SOLTY",                	prm.Solty,    	prm.Ntypes );
            write_darray( os, "LENNARD_JONES_ACOEF",  	prm.Cn1,      	prm.Nttyp );
            write_darray( os, "LENNARD_JONES_BCOEF",  	prm.Cn2,      	prm.Nttyp );

            write_bond( os, prm );
            write_angl( os, prm );
            write_tors( os, prm );
        
            write_iarray( os, "EXCLUDED_ATOMS_LIST",  	prm.ExclAt,   prm.Nnb  ); 
            write_darray( os, "HBOND_ACOEF", 		NULL, 0 );
            write_darray( os, "HBOND_BCOEF", 		NULL, 0 );
            write_darray( os, "HBCUT"      , 		NULL, 0 );
        
            write_sarray( os, "AMBER_ATOM_TYPE",	prm.AtomSym,	prm.Natom );
            write_sarray( os, "TREE_CHAIN_CLASSIFICATION",prm.AtomTree, prm.Natom );
            write_iarray( os, "JOIN_ARRAY",         	prm.TreeJoin, prm.Natom );
            write_iarray( os, "IROTAT",             	prm.AtomRes,  prm.Natom );

            if( prm.IfBox || prm.IfCap )
            {
                write_solvent_inf( os, prm );
            }

            write_string( os, "RADIUS_SET", GBPARM_NAME[2] );
            write_darray( os, "RADII",  prm.Rborn, prm.Natom );
            write_darray( os, "SCREEN", prm.Fs,    prm.Natom );  
 
        }             

    } // namespace prmtop
 
    using namespace prmtop;

    void write_amber_prmtop( ostream& os, molecule_t& m, const molecule_t& ffp )
    {
        energee_t e( m );

        e.assignparm( ffp );

        write_prmtop( os, e.getnabparm() );

 
        if( m.has_i(POLAR) )
        {
            write_polar( os, m );
        }

        int lestype;
        if( m.get_i(LESTYPE,lestype) && lestype > 0 )
        {
            write_les( os, m );
        }

    }

   
} // namespace mort

   
