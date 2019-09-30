#include "reader.hpp"
#include "parse.hpp"

namespace mort
{
    
    namespace daylight
    {
        
       
        reader_set::reader_set( int lan )
        {
            if( lan == SMILES )
            {
                init_smiles( *this );
            }
            else
            {
                assert( lan == SMARTS );
                init_smarts( *this );
            }    
        }
        
        reader_set::~reader_set()
        {
        }
        
        void reader_set::add( char c, int brack, reader_t reader )
        {
            int id = (brack == OPEN) ? c : c + 128;
            m_readers[ id ] = reader;
        }

        void reader_set::add( const char* ptr, int brack, reader_t reader )
        {
            for( ; *ptr; ptr++ )
            {
                add( *ptr, brack, reader );
            }
        }
        
        reader_t reader_set::get( char c, int brack ) const
        {
            int id = (brack == OPEN) ? c : c + 128;
            return m_readers[ id ];
        }
        
        void init_smiles( reader_set& readers )
        {
            readers.add( "0123456789", CLOSE, parse_ring );
            readers.add( ".-=#:/\\~",  CLOSE, parse_bond );
            readers.add( "[", CLOSE,  parse_brack );
            readers.add( "]", OPEN,   parse_brack );
            readers.add( "()", CLOSE, parse_paren );
        
            readers.add( "ABCEFGIJKLMNOPQSTUWYZcnos", OPEN,  parse_alpha );
            readers.add( "ABCEFGIJKLMNOPQSTUWYZcnos", CLOSE, parse_alpha );
            readers.add( "123456789", OPEN, parse_weight );
            readers.add( "+-", OPEN, parse_charge);
        }    

        void init_smarts( reader_set& readers )
        {
            init_smiles( readers );

            readers.add( "*", OPEN, parse_asterisk );
            readers.add( "#", OPEN, parse_sharp );
            readers.add( "a", OPEN, parse_aromatic );
            readers.add( "rR", OPEN, parse_rinfo );
            readers.add( "^DHVX", OPEN, parse_env );
            readers.add( "!,;", OPEN, parse_logic );
        }

    } // namespace daylight
    
} // namespace mort



    
