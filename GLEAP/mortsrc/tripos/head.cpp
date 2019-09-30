#include <sstream>
#include <object.hpp>
#include <common.hpp>

namespace mort
{
    using namespace std;
    
    namespace tripos
    {
        void read_size( istream& stream, molecule_t& mol, int compid )
        {
            string word = next_word( stream );
        
            if( !word.empty() )
            {
                int size = atoi( word.c_str() );
            
                if( size > 0 )
                {
                    mol.getcmp(compid)->resize( size );
                }
            }
        }

        void read_head( istream& stream, molecule_t& mol )
        {
            // cout << "reading mole ";
            string name, moltype, chgtype, status, comment, temp;
        
            // line 0
            getline( stream, name ); 
            mol.set_s(NAME, name);
        
            // line 1
            read_size( stream, mol, ATOM );
            read_size( stream, mol, BOND );
            read_size( stream, mol, RESD );        
            getline( stream, temp );

            // line 2
            getline( stream, moltype ); 

            // line 3
            getline( stream, chgtype );

            // line 4 : optional
            if( stream.peek() != '@' )
            {
                getline( stream, status );
            }
        
            // line 5: optional
            if( stream.peek() != '@' )
            {
                getline(stream, comment);
                mol.set_s( "comment", comment);
            }
        }

        void write_note( ostream& stream, const molecule_t& mol )
        {
            /*
            stream << "##################################################################" << endl;
            stream << "#" << endl;
	    string molname;
            stream << "#   Name       : " << ( mol.get_s(NAME, molname) ? molname : "untitled" ) << endl;
            stream << "#" << endl; 
            stream << "#   Created by : mort" << endl;
            stream << "#";
            stream << "##################################################################" << endl;
            */

            int solvent_shape;
            numvec solvent_info(4);
            if( mol.get_i(SOLUTE, solvent_shape) )
            {
                if(solvent_shape==CAP)
                {
                    mol.get_v(CAP, solvent_info);
                }
                else
                {
                    mol.get_v(BOX, solvent_info);
                }

                stream << "# ";
                stream << unhash(solvent_shape) << " ";
                stream << solvent_info[0] << " ";
                stream << solvent_info[1] << " ";
                stream << solvent_info[2] << " ";
                stream << solvent_info[3];
                stream << std::endl;
            }


        }

        void write_head( ostream& stream, const molecule_t& mol )
        {
            stream << "@<TRIPOS>MOLECULE\n";

            string molname;
            stream << ( mol.get_s(NAME, molname) ? molname : "untitled" ) << endl;

            stream << format( "%8d " ) % mol.natom();
            
            stream << format( "%8d " ) % mol.nbond();
            
            stream << format( "%8d " ) % mol.nresd() << endl;
            
            stream << ( mol.nresd() < 2 ? "SMALL" : "POLYMER" ) << endl;
            
            stream << "USER" << endl;

            stream << endl;

            stream << endl;
        }

    } // namespace tripos

} // namespace mort

