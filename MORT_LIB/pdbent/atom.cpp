#include <common.hpp>
#include <object.hpp>
#include "atom.hpp"

namespace mort
{
    namespace pdbent
    {

        void read_atom( istream& is, resd_t& r )
        {
            char line[MAX_LINE_WIDTH];
            is.getline(line, MAX_LINE_WIDTH);

            string name = read_atom_name( line + ATOM_NAME );
            int atomseq  = atoi( line + ATOM_SEQ );
            int element = pertab_t::get_element( name.c_str() );

            double crdx = atof( line + ATOM_X );
            double crdy = atof( line + ATOM_Y );
            double crdz = atof( line + ATOM_Z );

	    atom_t a = atom_t::create(r, name);
            a.set_i( ELEMENT,  element );
	    a.set_i( SEQUENCE, atomseq );
            a.set_v( POSITION, makevec(crdx,crdy,crdz) );
        }

  
        // in pdb format, atom name takes 4 characters 4 characters long.
        // If shorter than 4, the first character is blank, such as " CA ",
        // " CG1".
        // 4 character name should be reorganized such that the last one is
        // on head, as "3HB1", however some softwares write the name as it is, 
        // i.e. "HB13". Thus, if the
        // heading character is digit, it should be put to the back.
        string read_atom_name(const string& input)
        {
            int bgn = 0;
            while( bgn < 3 && input[bgn]==' ' ) bgn+=1;

            int end = std::min( (int)input.find(' ', bgn+1), 4 );
            if( bgn > 0 )
            {
                return input.substr(bgn, end-bgn);
            }

            if( isdigit(input[0]) )
            {
                return input.substr(1, end-1).append(1, input[0]);
            }
        
            return input.substr(0, end);
        }
    
        void write_atom( ostream& os, const atom_t& atom )
        {
            funstack_t::push( "pdbent::write_atom" );

            os << "ATOM  ";

            write_iparm( os, atom, ID, "%5d" );
            
            write_atom_name( os, atom );
            
            write_resd_info( os, atom );
            
            write_vparm( os, atom, POSITION, "%8.3f");
            
            os << "  1.00  0.00" << std::endl;

            funstack_t::pop();
        }
   
        void write_atom_name( ostream& os, const atom_t& atom )
        {
            string name = atom.get_s(NAME);
            
            if( name.length() < 4 )
            {
                os << format( "  %-3s " ) % name;
            }
            else 
            {
                os << " " << name[3];
                os << format( "%-3s " ) % name.substr( 0, 3 );
            }
        }
        
        void write_resd_info( ostream& os, const atom_t& atom )
        {
            string name = atom.get_s(RESNAME);
            
            os << format( "%3s " ) % ( name.length() > 3 ? name.substr( 0, 3 ) : name );
 
            if( name.length() > 3 && !isdigit( name[3] ) )
            {
                os << name[3];
            }
            else
            {
                os << ' ';
            }

            os << format( "%4d    " ) % (atom.resd().relid()+1);
        }
       
    } // namespace pdbent

} // namespace mort



