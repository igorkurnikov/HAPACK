#include <object.hpp>

#include "les.hpp"
#include "exclude.hpp"

namespace mort
{
    namespace prmtop
    {

        void build_list(const morf_t& atom, vector<int>& list, vector<int>& dist, int max_level)
        {
            atomvec_t nbrs;            
            nbrs.push_back( atom );
            
            int level = 1;
            int start = 0;

            while( level <= max_level )
            {
                int end = nbrs.size();

                
                for( int i=start; i < end; i++ )
                {
                    atomiter_t nbr = nbrs[i].atom_begin();
                    for( ; nbr != nbrs[i].atom_end(); ++nbr )
                    {
                        if( std::find( nbrs.begin(), nbrs.end(), *nbr) != nbrs.end() )
                        {
                            continue;
                        }
                        
                        nbrs.push_back( *nbr );
                            
                        if(*nbr > atom)
                        {
                            list.push_back( nbr->get_i(ID) );
                            dist.push_back( level );
                        }
                    }
                }

                start = end;                
                level++;
            }

            les_exclude( atom, list, dist );
        }
        

        void sort_key( vector< int >& key, vector< int >& value )
        {
            assert( key.size() == value.size() );
            
            for( int i=0; i < (int)key.size(); ++i )
            {
                for( int j=i+1; j < (int)key.size(); ++j )
                {
                    if( key[i] > key[j] )
                    {
                        int tmp;
                        tmp = key[i];
                        key[i] = key[j];
                        key[j] = tmp;
                        
                        tmp = value[i];
                        value[i] = value[j];
                        value[j] = tmp;
                    }
                }
            }
        }

        struct exclude_atom
        {
            exclude_atom( excl_t& excl, int level )
                : m_excl( &excl )
            {
                m_level = level;
            }
        
            void operator()(const morf_t& atom)
            {
                int id = atom.absid();

                build_list( atom, m_excl->list[ id ], m_excl->dist[id], m_level );

                sort_key( m_excl->list[id], m_excl->dist[id] );

                if( m_level == 3 && m_excl->list[ id ].size() == 0 )
                {
                    m_excl->list[ id ].push_back( 0 );
                    m_excl->dist[ id ].push_back( 0 );
                }
            }

            excl_t* m_excl;

            int m_level;
        };
        
    } // namespace prmtop
    
    void exclude( const molecule_t& mol, excl_t& excl, int level )
    {
        excl.list.resize( mol.natom() );
            
        excl.dist.resize( mol.natom() );
            
        std::for_each( mol.atom_begin(), mol.atom_end(), prmtop::exclude_atom( excl, level ) );
    }
    
} // namespace mort


