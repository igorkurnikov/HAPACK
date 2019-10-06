#include <cassert>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include "hashcode.hpp"

#ifndef __GNUC__
#include <ctype.h>
#endif

namespace mort
{
    hashid_t hash(const string& str)
    {
        assert( sizeof(hashid_t)==8 );

        hashid_t sum=0;

	int len = std::min( (int)str.length(), 12);

        for(int i=len-1; i>=0; --i)
	{
	    int c = (int)str[i];
      
            if( islower(c) )
	    {
	        c = c - 'a';
	    }
	    else if( isupper(c) )
	    {
	        c = tolower(c) - 'a';
	    }
	    else if( isdigit(c) )
	    {
	        c -= '0';
		c += 26;
	    }
            else if( c=='_' )
            {
                c = 36;
            }
            else
	    {
	        continue;
	    }

	    sum = 40*sum + c;
	}

	return sum;
    }

    string unhash(const hashid_t& id)
    {
        hashid_t code = id;
        string str;
        while(code != 0LL)
        {
	    hashid_t next = code/40LL;
	    hashid_t c = code - 40LL*next;
            assert(c < 40LL);
            if( c < 26LL )
	    {
	        c = c + 'a';
	    }
	    else if( c < 36LL )
	    {
	        c = c - 26LL + '0';
	    }
            else
            {
                assert( c==36LL );
                c = '_';
            }

            str.append(1, (char)c);
	
            code = next;
	}

	return str;
    }

} // namespace mort


