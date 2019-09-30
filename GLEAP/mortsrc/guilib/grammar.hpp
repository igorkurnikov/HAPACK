#ifndef MORTSRC_GUILIB_GRAMMAR_HPP
#define MORTSRC_GUILIB_GRAMMAR_HPP
#include <iosfwd>
#include <vector>
#include <string>

namespace mort
{
    using std::vector;

    using std::string;

    int ndim( const string& list );

    void interpret( const string& command, vector<string>& args );

    vector<string> interpret1d(const string& list);

    vector< vector<string> > interpret2d(const string& list);

    bool curly_closed(const string& list);

    bool quota_closed(const string& list);

    string parse( std::istream& is );

    string parse_curly( std::istream& is );

    string parse_quota( std::istream& is );
}

#endif
