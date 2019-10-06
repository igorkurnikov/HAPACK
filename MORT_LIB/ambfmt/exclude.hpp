#ifndef MORTSRC_AMBERFMT_EXCLUDE_HPP
#define MORTSRC_AMBERFMT_EXCLUDE_HPP

#include <vector>

namespace mort
{
    using std::vector;

    struct excl_t
    {
        int full_size( ) const
        {
            int sum=0;
                
            for( int i=0; i < (int)list.size(); i++ )
            {
                sum += list[i].size();
            }

            return sum;
        }
        
        int size( int i ) const
        {
            return list[i].size();
        }

        int found( int i, int j )
        {
            assert( i >=0 && i < (int)list.size() );

            if( list[i].end() != std::find( list[i].begin(), list[i].end(), j+1 ) )
            {
                return true;
            }

            if( list[j].end() != std::find( list[j].begin(), list[j].end(), i+1 ) )
            {
                return true;
            }

            return false;
        }

        void finish( vector<int>& exsize, vector<int>& exlist )
        {
            for( unsigned int i=0; i < list.size(); ++i )
            {
                exsize.push_back( list[i].size() );
                for( unsigned int j=0; j < list[i].size(); ++j )
                {
                    exlist.push_back( list[i][j] );
                }
            }
        }

        
        vector< vector< int > > list;
        vector< vector< int > > dist;
    };

     void exclude( const molecule_t& mol, excl_t& excl, int level );

   
} // namespace mort

#endif

