#ifndef MHTL_SYBYL_ATOM_HPP
#define MHTL_SYBYL_ATOM_HPP

#include <iostream>
#include <object.hpp>

namespace mort
{
    using namespace std;

    namespace tripos
    {

        struct atominfo_s
        {
            atominfo_s() 
                : position(3) 
            {}

            atominfo_s( const atominfo_s& rhs ) 
                : name( rhs.name ), type( rhs.type ), resname( rhs.resname ),
                  position( rhs.position )
            {
                pchg = rhs.pchg;
                id = rhs.id;
                resid = rhs.resid;
            }

            atominfo_s& operator=( const atominfo_s& rhs )
            {
                name = rhs.name;
                type = rhs.type;
                resname = rhs.resname;
                position = rhs.position;
                pchg = rhs.pchg;
                id = rhs.id;
                resid = rhs.resid;
                return *this;
            }


            string name;
            string type;
            string resname;
            numvec position;
            double pchg;
            int id;
            int resid;
        };

        bool is_one_of( const string& resname, const char names[][10], int ntype )
        {
            for( int i=0; i < ntype; ++i )
            {
                int len = strlen( names[i] );
                if( resname.substr(0, len)==names[i] )
                    return true;
            }

            return false;
        }

        bool is_ion( const string& resname )
        {
            const char IONS[][10] = { "Br-", "Cl-", "Cs+", "F-", "I-", "K+", "Li+", "Mg+", "Na+", "Rb+" };
            int nion = 10;
            return is_one_of( resname, IONS, nion );
        }

        bool is_solvent( const string& resname )
        {
            const char SOLVENTS[][10] = { "CL3", "WAT", "HOH", "MOH", "NMA" };
            int nsolvent = 5;
            return is_one_of( resname, SOLVENTS, nsolvent );
        }

        void add_resds( vector<atominfo_s>& atominfos, vector< vector<atominfo_s*> >& resd_atoms, int nresd, vector<int>& rids )
        {
            for( int i=0; i < (int)rids.size(); ++i )
            {
                int id = rids[i];
                for( int j=0; j < (int)resd_atoms[id].size(); ++j )
                {
                    atominfos.push_back( *resd_atoms[id][j] );
                    atominfos.back().resid = nresd + i + 1;
                }
            }
        }

        // this function does the following thing:
        // 1. make sure atoms from same residue are packed together
        // 2. make the sequence of residues are arranged as:
        //     study object --> counter ions --> solvents
        //
        void rearrange( vector<atominfo_s>& atominfos )
        {
            vector<int> ions;
            vector<int> solute;
            vector<int> solvents;

            vector< vector<atominfo_s*> > resd_atoms;
            int natom = atominfos.size();
            for( int i=0; i < natom; ++i )
            {
                int resid = atominfos[i].resid;

                bool found = false;
                for( int j=0; j < (int)resd_atoms.size(); ++j )
                {
                    assert( resd_atoms[j].size() > 0 );
                    if( resd_atoms[j][0]->resid == resid )
                    {
                        found = true;
                        resd_atoms[j].push_back( &atominfos[i] );
                        break;
                    }
                }

            
                if( found ) continue;

                // new residue
                resd_atoms.push_back( vector<atominfo_s*>(1, &atominfos[i]) );
                
                string resname = atominfos[i].resname;
                if( is_ion(resname) )
                {
                    ions.push_back( resd_atoms.size()-1 );
                }
                else if( is_solvent(resname) )
                {
                    solvents.push_back( resd_atoms.size()-1 );
                }
                else
                {
                    solute.push_back( resd_atoms.size()-1 );
                }

            }

            
            vector< atominfo_s > replacements;
            add_resds( replacements, resd_atoms, 0, solute );
            add_resds( replacements, resd_atoms, solute.size(), ions );
            add_resds( replacements, resd_atoms, solute.size()+ions.size(), solvents );
            atominfos.swap( replacements );
        }         


        void read_atom( istream& is, molecule_t& mol )
	{
            vector<atominfo_s> atominfos( mol.natom() );
            for( int i=0; i < mol.natom(); ++i )
            {
                is >> atominfos[i].id;
                is >> atominfos[i].name;
                is >> atominfos[i].position[0];
                is >> atominfos[i].position[1];
                is >> atominfos[i].position[2];
                is >> atominfos[i].type;
                is >> atominfos[i].resid;
                is >> atominfos[i].resname;
                is >> atominfos[i].pchg;
                is.ignore( MAX_LINE_WIDTH, '\n' );   
            }

            string re;
            if( mortenv().get_s("rearrangeResidue", re) && re=="on" )
            {
                rearrange( atominfos );
            }

	    atomiter_t ai = mol.atom_begin();
	    atomiter_t ae = mol.atom_end();
	    for( int i=0; ai != ae; ++ai, ++i )
	    {

                ai->set_i( SEQUENCE, atominfos[i].id );
                ai->set_s( NAME, atominfos[i].name );
		ai->set_v( POSITION, atominfos[i].position );
		ai->set_s( TYPE, atominfos[i].type );
		ai->set_i( RESID, atominfos[i].resid );
	        ai->set_s( RESNAME, atominfos[i].resname );
		ai->set_d( PCHG, atominfos[i].pchg );

	    }
        }

        void write_atom( ostream& os, const molecule_t& mol )
        {
            os << "@<TRIPOS>ATOM" << endl;

            atomiter_t ai = mol.atom_begin();
	    atomiter_t ae = mol.atom_end();
            for( ; ai != ae; ++ai )
            {
                write_iparm(os, *ai, ID, "%8d ");
                write_sparm(os, *ai, NAME, "%-8s ");
                write_vparm(os, *ai, POSITION, "%9.3f ");
                write_sparm(os, *ai, TYPE, "%-8s ");
                write_iparm(os, *ai, RESID, "%8d ");
                write_sparm(os, *ai, RESNAME, "%-8s ");
                write_dparm(os, *ai, PCHG, "%9.3f ");
                os << endl;
            }
        }
        
    } // namespace tripos
    
} // namespace mort

#endif

