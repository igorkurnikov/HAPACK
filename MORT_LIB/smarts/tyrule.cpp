#include <fstream>
#include <sstream>
#include "tyrule.hpp"

namespace mort
{
    item_t::item_t( int size )
        : std::vector< std::string >( size )
    {
    }

    void item_t::read_pattern( int id )
    {
        read_smarts( at(id).c_str(), pattern );
    }

    typing_rule::typing_rule( const char* filename )
    {
        std::ifstream stream( filename );

        if( !stream )
        {
            throw std::logic_error( "can't open file " + std::string( filename ) );
        }
        
        while( stream )
        {
            char tag;
            stream >> tag;
            
            if( tag == '#' )
            {
                stream.ignore( 160, '\n' );
            }
            else
            {
                if( !stream.eof() )
                {
                    stream.putback( tag );
                    m_title.empty() ? read_title( stream ) : read_item( stream );
                }
            }
        }
    }

    typing_rule::~typing_rule()
    {
    }
    
    void typing_rule::read_title( std::istream& is )
    {
        string word, line;
	std::getline( is, line );

	std::istringstream ls( line );
        while( ls >> word )
        {
            m_title[ mort::hash(word) ] = m_title.size() - 1;
        }

    }
    
    void typing_rule::read_item( istream& stream )
    {
        BOOST_ASSERT( m_title.count( mort::hash("pattern") ) );

        push_back( item_t( m_title.size() ) );

        for( int i=0; i < (int)back().size(); i++ )
        {
            stream >> back()[i];
        }

        back().read_pattern( m_title[ mort::hash( "pattern" ) ] );
        stream.ignore( MAX_LINE_WIDTH, '\n' );
    }

    const char* typing_rule::get(const morf_t& atom, const char* name) const
    {
        BOOST_ASSERT( m_title.count(mort::hash(name)) );

        size_t id = m_title.find(mort::hash(name))->second;

        const_reverse_iterator i = rbegin();
        for( ; i != rend(); i++ )
        {
            if( has_subst( atom, i->pattern ) )
            {
                return i->at( id ).c_str();
            }
        }

        throw std::logic_error( "can't determine atom type" );  
    }
    
    const char* typing_rule::get( int i, const char* name ) const
    {
        BOOST_ASSERT( m_title.count(mort::hash(name)) );

        int id = m_title.find(mort::hash(name))->second;
 
        return at( i ).at( id ).c_str();
    }

} // namespace mort

