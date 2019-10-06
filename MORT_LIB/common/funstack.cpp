#include <iostream>
#include "funstack.hpp"

namespace mort
{
    using std::string;

    void funstack_t::print( )
    {
        print( std::cout );
    }

    void funstack_t::print( std::ostream& os )
    {
        for(unsigned int i=0; i < inst().size(); ++i )
        {
            for(unsigned int j=0; j < i; ++j )
                os << "    ";

            os << "calling \"" << inst()[i][0] << "( ";
            for(unsigned int j=1; j < inst()[i].size(); ++j )
            {
                os << " " << inst()[i][1];
                if( j!=inst()[i].size()-1 )
                    os << ",";
            }
            os << ")\"" << std::endl;
        }
    }

    void funstack_t::push( const string& name )
    {
        inst().push_back( std::vector<string>(1, name) );
    }

    void funstack_t::push( const string& name, const string& arg1 )
    {
        inst().push_back( std::vector<string>(1, name) );
        inst().back().push_back( arg1 );
    }

    void funstack_t::pop( )
    {
        inst().pop_back();
    }

} // namespace mort

