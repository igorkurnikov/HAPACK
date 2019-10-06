#ifndef MORTSRC_COMMON_FUNSTACK_HPP
#define MORTSRC_COMMON_FUNSTACK_HPP

#include <vector>
#include <string>


namespace mort
{

    // funstack (function stack) is a tool for debugging the code.
    // When coming to a function (or member function), push the function and argument,
    // when error encountered, the stack trace can be print which is valuable information
    // for locating the error.

    // the current version does not consider multi-threading cases
    class funstack_t
    {
    private:

        typedef std::vector< std::vector<std::string> > data_t;
        
        // ensure singleton
        funstack_t() {}

        virtual ~funstack_t() {}

        static data_t& inst()
        {
            static data_t data;
            return data;
        }
            
    public:

        static void print( );

        static void print( std::ostream& os );

        static void push( const std::string& name );

        static void push( const std::string& name, const std::string& arg1 );

        static void pop( );
    };

} // namespace mort

#endif
