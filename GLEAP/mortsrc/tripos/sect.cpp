#include <istream>
#include <object.hpp>

namespace mort
{
    namespace tripos
    {
        using namespace std;

        void read_head( istream& stream, molecule_t& mol );

        void read_atom( istream& stream, molecule_t& mol );
        
        void read_bond( istream& stream, molecule_t& mol );
        
        void read_resd( istream& stream, molecule_t& mol );
        
        bool read_sect( istream& stream, molecule_t& mol )
        {
            std::ios::pos_type pos = stream.tellg();
            string keyw;
            getline( stream, keyw );

            // std::cout << "reading section: " << keyw << std::endl;
            
            if( keyw.find( "<TRIPOS>ATOM" ) != string::npos )
            {
                read_atom( stream, mol );
            }
            else if( keyw.find( "<TRIPOS>BOND" )!= string::npos)
            {
                read_bond( stream, mol );
            }
            else if( keyw.find( "<TRIPOS>SUBSTRUCTURE" ) != string::npos )
            {
                read_resd( stream, mol );
            }
            else if( keyw.find( "<TRIPOS>MOLECULE" ) != string::npos )
            {
                if( mol.natom() > 0 )
                {
                    stream.seekg( pos );
                    return false;
                }
                else
                {
                    read_head( stream, mol );
                }
            }
            else
            {
                throw logic_error( "unknown sybyl keyword " + keyw );
            }
        
            return true;
        }

    } // namespace sybyl
    
} // namespace mort
