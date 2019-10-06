#include <fstream>
#include <object.hpp>
#include "xyz.hpp"

namespace mort
{

    void write_amber_xyz( ostream& os, const molecule_t& mol )
    {
        string molname;
        os << ( mol.get_s(NAME, molname) ? molname : "untitled" ) << std::endl;
        
        os << format( "%6d" ) % mol.natom() << std::endl;

        atomiter_t atom = mol.atom_begin();
        for( ; atom != mol.atom_end(); ++atom )
        {
            numvec pos = atom->get_v(POSITION);
            
            os << format( "%12.7f" ) % pos[0];
            os << format( "%12.7f" ) % pos[1];
            os << format( "%12.7f" ) % pos[2];
 
 
            if( atom->absid() % 2 == 1 )
            {
                os <<std::endl;
            }
  
        }
    
    
        if( mol.natom() % 2 == 1 )
        {
            // os << "        ";
            os <<std::endl;
        }
  

        numvec box(4);
        if( mol.get_v(BOX, box) )
        {
            os << format( "%12.7f" ) % box[0];
            os << format( "%12.7f" ) % box[1];
            os << format( "%12.7f" ) % box[2];
            os << format( "%12.7f" ) % box[3];
            os << format( "%12.7f" ) % box[3];
            os << format( "%12.7f" ) % box[3];
            os << "        " <<std::endl;
        }
        
    }

    void write_tinkerxyz( ostream& os, const molecule_t& mol )
    {
        string molname;
        os << format( "%6d" ) % mol.natom() << "  ";
        os << ( mol.get_s(NAME, molname) ? molname : "Untitled" ) << std::endl;

        atomiter_t atom = mol.atom_begin();
        for( ; atom != mol.atom_end(); ++atom )
        {
            numvec pos = atom->get_v(POSITION);
            os << format( "%6d"    ) % atom->get_i(ID)    << "  ";
            os << format( "%-3s"   ) % atom->get_s(SYMBOL) << " ";
            os << format( "%11.6f" ) % pos[0] << " ";
            os << format( "%11.6f" ) % pos[1] << " ";
            os << format( "%11.6f" ) % pos[2] << " ";
            os << format( "%5s"    ) % atom->get_s(POLTYPE) << " ";

            atomiter_t nbr = atom->atom_begin();
            for( ; nbr != atom->atom_end(); ++nbr )
            {
                os << format( "%5d" ) % nbr->get_i(ID) << " ";
            }

            os << std::endl;
        }
    }
}
