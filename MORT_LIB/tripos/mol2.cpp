#include <object.hpp>
#include <tripos.hpp>

namespace mort
{
    bool read_mol2( istream& stream, molecule_t& mol )
    {        
        molecule_t temp;

        char tag;

        bool end = false;

        while( !end && stream >> tag )
        {
            if( tag == '#' )
            {
                string head = next_word( stream );

                if( !empty(head) && (head=="cap" || head=="box" || head=="oct") )
                {
                    string solvent_shape = head;
                    temp.set_i( SOLUTE, hash(solvent_shape) );
                    numvec solvent_info( 4 );
                    for( int i=0; i < 4; ++i )
                    {
                        string s = next_word( stream );
                        solvent_info[i] = atof( s.c_str() );
                    }
 
                    if( solvent_shape=="cap" )
                    {
                        temp.set_v( CAP, solvent_info );
                    }
                    else
                    {
                        temp.set_v( BOX, solvent_info );
                    }
                }               

                stream.ignore(MAX_LINE_WIDTH, '\n' );
            }
            else if( tag == '@' )
            {
                end = ! tripos::read_sect( stream, temp );
            }
            else
            {
                // std::cout << "tag: " << tag << std::endl;
                throw std::runtime_error( "unknown tag!\n" );
            }            
        }

        stream.putback( '@' );

        mol.swap( temp );
		return true;
    }

    void write_mol2( ostream& stream, const molecule_t& mol )
    {
        tripos::write_note( stream, mol );

        tripos::write_head( stream, mol );

        tripos::write_atom( stream, mol );
    
        tripos::write_bond( stream, mol );
    
        tripos::write_resd( stream, mol );
    
    }
    
} // namespace mort
 
