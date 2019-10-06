#include <stdexcept>
#include <boost/algorithm/string.hpp>
#include "ctrlparm.hpp"

#ifndef __GNUC__
#include <stdlib.h>
#endif

namespace mort
{
    using namespace boost;

    ctrlparm_t::ctrlparm_t(const string& str)
    {
        rgbmax = 25.0;
	offset = 0.09;

        vector<string> terms;
	split(terms, str, is_any_of(" ,:"), token_compress_on);
	parse( terms );
    }

    void ctrlparm_t::parse(const vector<string>& terms)
    {
        assert( terms.size()%2==0 );
	int nterm = terms.size()/2;
	for(int i=0; i < nterm; ++i)
	{
	    if( terms[2*i] == "igb" )
	    {
	        igb = atoi( terms[2*i+1].c_str() );
	    }
	    else if( terms[2*i] == "cut" )
	    {
	        cut = atof( terms[2*i+1].c_str() );
	    }
	    else
	    {
	        throw std::logic_error("Error: unknown ctrl parm " + terms[2*i]);
	    }
	}
    }

} // namespace mort


