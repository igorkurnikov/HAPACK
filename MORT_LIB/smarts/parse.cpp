#include <map>
#include <object.hpp>
#include "reader.hpp"
#include "status.hpp"
#include "atomexpr.hpp"

namespace mort
{
    namespace daylight
    {
        
        void parse( const char* exp, molecule_t& mol, const reader_set& readers, int lan )
        {
            status_t state( lan, mol );

            while( *exp != '\0' )
            {
                reader_t reader = readers.get( *exp, state.bracket );

                try
                {
                    exp = reader( exp, state, mol ); 
                }
                catch( std::exception& e )
                {
                    std::cout << "can't find function for " << exp << std::endl;
                }
                
            }
        }

        const char* read_alpha( const char* ptr, int& element, int& arom )
        {
            int len = std::min( strlen( ptr ), (size_t)2 );

            element = pertab_t::get_element( string( ptr, ptr + len ).c_str() );
            arom = false;

            const char* next = ptr + len;
        
            if( element == 0 && islower( *ptr ) )
            {
                element = pertab_t::get_element( toupper( *ptr ) );
                arom = true;
                next = ptr + 1;
            }
            
            if( element == 0 )
            {
                element = pertab_t::get_element( *ptr );
                arom = false;
                next = ptr + 1;
            }

            BOOST_ASSERT( element != 0 );

            return next;
        }        

        using namespace atomexpr;

	typedef shared_ptr< atomexpr_t > atomexpr_ptr;

        const char* parse_alpha( const char* ptr, status_t& status, molecule_t& )
        {            
            BOOST_ASSERT( isalpha( *ptr ) );
            

            if( status.bracket == CLOSE )
            {
                status.new_atom( );
            }
        
            int element;
            int arom;

            ptr = read_alpha( ptr, element, arom );
        
            if( status.language == SMARTS )
            {
                atomexpr_ptr pattern = any_cast< atomexpr_ptr>( status.curt_atom.get_a(PATTERN) );

                pattern->insert( shared_ptr< node_i >( new iparm_node(ELEMENT,element) ) );
                
                pattern->insert( shared_ptr< node_i >( new iparm_node(AROM,arom) ) );
            }
            else
            {
                assert( status.language == SMILES );

                status.curt_atom.set_i(ELEMENT, element);

                status.curt_atom.set_i(AROM, arom);
            }

            if( status.bracket == CLOSE )
            {
                status.end_atom( );
            }
        
            return ptr;
        }
        
        const char* parse_rinfo( const char* ptr, status_t& status, molecule_t& )
        {
            assert( (*ptr == 'r' || *ptr == 'R' ) && status.language == SMARTS );

            int arom = ( *ptr == 'r' );

            int size = atoi( ptr + 1 );
            
            atomexpr_ptr tmp = any_cast< atomexpr_ptr >( status.curt_atom.get_a(PATTERN) );
	    tmp->insert( shared_ptr< node_i >( new ring_node( arom, size ) ) );

            return skip_digit( ptr + 1 );
        }   

        const char* parse_asterisk( const char* ptr, status_t& status, molecule_t& )
        {
            BOOST_ASSERT( *ptr == '*' && status.language == SMARTS );

            atomexpr_ptr tmp = any_cast< atomexpr_ptr >( status.curt_atom.get_a(PATTERN) );
	    tmp->insert( shared_ptr< node_i >( new true_node() ) );

            return ptr + 1; 
        }

        const char* parse_sharp( const char* ptr, status_t& status, molecule_t& )
        {
            BOOST_ASSERT( *ptr == '#' && status.language == SMARTS );

            int element = atoi( ptr + 1 );

            atomexpr_ptr tmp = any_cast< atomexpr_ptr >( status.curt_atom.get_a(PATTERN) );
	    tmp->insert( shared_ptr< node_i >(new iparm_node(ELEMENT, element)) );
            
            return skip_digit( ptr + 1 );
        }

        const char* parse_aromatic( const char* ptr, status_t& status, molecule_t& )
        {
            assert( *ptr == 'a' && status.language == SMARTS );
                
            atomexpr_ptr tmp = any_cast< atomexpr_ptr >( status.curt_atom.get_a(PATTERN) );
	    tmp->insert( shared_ptr< node_i >(new iparm_node(AROM,true)) );
            
            return ptr + 1;
        }
        
        const char* parse_bond( const char* ptr, status_t& status, molecule_t& )
        {
            bondinfo_t& next_bond  = status.bond_info;

            next_bond.reset();

            switch( *ptr )
            {
            case '.':
                next_bond.order = 0;
                break;
            case '-':
                break;
            case '=':
                next_bond.order = 2;
                break;
            case '#':
                next_bond.order = 3;
                break;
            case ':':
                next_bond.arom = true;
                break;
            case '/':
                next_bond.stereo = 1;
                break;
            case '\\':
                next_bond.stereo = 6;
                break;
            case '~':
                next_bond.order = 99;
                break;
            }

            return ptr+1;
        }


        const char* parse_brack( const char* ptr, status_t& status, molecule_t& )
        {
            if( status.bracket == CLOSE )
            {
                assert( *ptr == '[' );

                status.bracket = OPEN;
                status.new_atom( );
            }
            else
            {
                assert( *ptr == ']' );
            
                status.bracket = CLOSE;
                status.end_atom( );
            }
            
            return ptr + 1;
        }        

        const char* parse_charge( const char* ptr, status_t& status, molecule_t& )
        {
            assert( *ptr == '+' || *ptr == '-' );

            int sign = ( *ptr == '+' ) ? 1 : -1;
            
            ptr++;

            int quantity = 1;

            if( isdigit( *ptr ) )
            {
                quantity = atoi( ptr );
                ptr = skip_digit( ptr );
            }
            else
            {
                while( *ptr == '+' || *ptr == '-' )
                {
                    quantity +=( *ptr == '+' ? 1 : -1 );                    
                    ptr++;
                }
            }

            int charge = sign * quantity;
       
            if( status.language == SMARTS )
            {   
                atomexpr_ptr tmp = any_cast< atomexpr_ptr >( status.curt_atom.get_a(PATTERN) );
		tmp->insert( shared_ptr< node_i >(new iparm_node(FCHG, charge)) );
            }
            else
            {
                status.curt_atom.set_i(FCHG, charge);
            }
 
            return ptr;
        }

        env_e get_env( char c )
        {
            switch( c )
            {
            case '^':
                return HYBRID;
            case 'D':
                return DEGREE;
            case 'H':
                return HYDROGEN;
            case 'V':
                return VALENCE;
            case 'X':
                return NEBR;
            }

	    throw std::runtime_error( string("Error: unrecognizable env character").append(1,c) );
        }

        const char* parse_env( const char* ptr, status_t& status, molecule_t& )
        {
            env_e env = get_env( *ptr );

            int count = 0;

            if( isdigit( *( ptr + 1 ) ) )
            {
                count = atoi( ptr + 1 );
            }
            else
            { 
                count = (env == HYDROGEN) ? 1 : 0;
            }

            atomexpr_ptr tmp = any_cast< atomexpr_ptr >( status.curt_atom.get_a(PATTERN) );
	    tmp->insert( shared_ptr< node_i >(  new env_node( env, count ) ) );
            
            return skip_digit( ptr + 1 );
        }

        std::pair< int, int > get_operation( char c )
        {
            static const char symbols[] = "!,;";

            static const int logics[] = { NOT, OR, AND };
            
            static const int levels[] = { SURPRISE, COMMA, SEMICOLON };            

            assert( std::count( symbols, symbols + strlen( symbols ), c ) );

            int id = std::find( symbols, symbols + strlen( symbols ), c ) - symbols;
            
            return std::make_pair( logics[id], levels[id] ) ;
        }
        
        const char* parse_logic( const char* ptr, status_t& status, molecule_t& )
        {
            std::pair<int, int> oper( get_operation( *ptr ) );
            
	    atomexpr_ptr tmp = any_cast< atomexpr_ptr >( status.curt_atom.get_a(PATTERN) );
	    tmp->change_to( oper.second, oper.first );

            return ptr + 1;
        }

        const char* parse_paren( const char* ptr, status_t& status, molecule_t& mol )
        {
            assert( *ptr == '(' || *ptr == ')' );

            vector< morf_t >& atoms = status.stack;

            if( *ptr == '(' )
            {
                if( atoms.size() > 0 && atoms.back() == status.curt_atom )
                {
                    atoms.push_back( morf_t(mol, ATOM, -1) );
                }
                else
                {
                    atoms.push_back( status.curt_atom );
                }
            }
            else
            {
                assert( atoms.size() > 0 );
                
                if( atoms.back().absid() != -1 )
                {
                    status.curt_atom = atoms.back();
                }

                atoms.pop_back();
            }

            return ptr + 1;
        }

        const char* parse_ring( const char* ptr, status_t& status, molecule_t& )
        {
            vector<atom_t>& ring_atoms = status.ring_atoms;
        
            int rid = *ptr - '0';

            if( ring_atoms[ rid ].absid() != -1 )
            {
                bond_t bond = bond_t::create( status.curt_atom, ring_atoms[ rid ] );
                bond.set_i(ORDER, 1);
                bond.set_i(AROM, false);
            }
            else
            {
                ring_atoms[ rid ] = status.curt_atom;
            }

            return ptr + 1;
        }

        const char* parse_weight( const char* ptr, status_t& status, molecule_t& )
        {
            int weight = atoi( ptr );
            
            if( status.language == SMARTS )
            {
                atomexpr_ptr tmp = any_cast< atomexpr_ptr >( status.curt_atom.get_a(PATTERN) );
	        tmp->insert( shared_ptr< node_i >(new iparm_node(WEIGHT, weight)) );
            }
            else
            {
                assert( status.language == SMILES );

                status.curt_atom.set_i(WEIGHT, weight);
            }

            return skip_digit( ptr );
        }

    } // namespace daylight
    
} // namespace mort



    
