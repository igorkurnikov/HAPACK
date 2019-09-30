#ifndef AMBER_GLEAP_STRBUFF_HPP
#define AMBER_GLEAP_STRBUFF_HPP

#include <string>
#include <guilib.hpp>

namespace amber
{
    using std::string;

    class std_console : public mort::console_t
    {
    public:

        std_console( std::istream& is, std::ostream& os );
        
        virtual ~std_console();
        
        virtual void print( const std::string& str );

        virtual char getchar( );

        virtual bool getline( string& line, const char* prompt=NULL );
 
    private:
        
        std::istream* m_is;
   
        std::ostream* m_os;

	bool m_batch;
    };
    
} // namespace amber

#endif
