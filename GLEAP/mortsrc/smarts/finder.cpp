#include <object.hpp>
#include "finder.hpp"

namespace mort
{
    finder_t::finder_t(int strategy)
    {
        m_allow_multi_visit = (strategy==ALLOW_MULTI_VISIT);

    }

    finder_t::~finder_t()
    {
    }

    void finder_t::init(const morf_t& begin)
    {
        m_path.push_back(begin);
	m_visited.insert( begin.absid() );

        init();
    }

    void finder_t::init(const atomvec_t& path)
    {
        m_path.insert( m_path.end(), path.begin(), path.end() );
	init();
    }

    void finder_t::init(const atomiter_t& begin, const atomiter_t& end)
    {
        m_path.insert(m_path.end(), begin, end);
	init();
    }

    void finder_t::init()
    {
    /*
        for( int i=0; i < (int)m_path.size(); ++i )
	{
	    m_visited.insert( m_path[i].absid() );
	}
     */
        assert( m_path.size() > 0 );
        m_visited.insert( m_path.back().absid() );
	if( predict(m_path)==ENOUGH )
	    return;

	atomvec_t next;
	list_next(m_path, next);
	m_nexts.push_back(next);
    }
	
    bool finder_t::dead()
    {
        for( int i=0; i < (int)m_nexts.size(); ++i )
	{
	    if(m_nexts[i].size() > 0)
	    {
	        return false;
	    }
	}

	return true;
    }

    void finder_t::back_one()
    {
        assert( ! m_path.empty() && ! m_nexts.empty() );
        m_path.pop_back();
        m_nexts.pop_back();
    }

    morf_t finder_t::next_atom()
    {
        assert( m_nexts.size() > 0 );

        while( m_nexts.back().empty() )
        {
            back_one();
        }

        assert( m_nexts.size() > 0 && m_nexts.back().size() > 0);

        morf_t atom = m_nexts.back().back();
        m_nexts.back().pop_back();
        return atom;
    }

    int finder_t::step( )
    {
        m_path.push_back( next_atom() );
        m_visited.insert( m_path.back().absid() );
        int result = predict(m_path);
	atomvec_t next;

        // std::cout << "visiting: " << m_path.back().get_i(ID) << " result: " << result << std::endl;


	if( result == KEEP_GOING )
	{
	    list_next(m_path, next);
        }

	m_nexts.push_back(next);
	return result;
    }

    bool finder_t::find_one( atomvec_t& path )
    {

        if( predict(m_path) == ENOUGH )
	{
	    path = m_path;
	    return true;
	}

        while( !dead() )
	{
	    int result = step();
	    if( result ==  ENOUGH )
	    {
	        path = m_path;
	        return true;
	    }
	    else if( result == WRONG_DIRECTION )
	    {
	        back_one();
            }
        }

        return false;

    }

    bool finder_t::find_all(vector<atomvec_t>& paths)
    {
        atomvec_t path;
	while( find_one(path) )
	{
	    paths.push_back(path);
	    back_one();
	    path.clear();
	}
	return paths.size()>0;
    }

    void finder_t::list_next(const atomvec_t& path, atomvec_t& next)
    {
        assert( path.size() > 0 );

        morf_t curt = path.back();
	
	atomiter_t i = curt.atom_begin();
	for( ; i != curt.atom_end(); ++i )
	{
	    if( path.count( *i ) > 0 )
	    {
	        continue;
	    }

	    if( !m_allow_multi_visit && m_visited.count( i->absid() ) > 0 )
	    {
	        continue;
	    }

	    next.push_back(*i);
	}

/*
        std::cout  << "for atom " << curt.get_i(ID) << " will visit ";
	for( int j=0; j < next.size(); ++j )
	{
	    std::cout << next[j].get_i(ID) << " ";
	}

	std::cout << std::endl;
 */

    }

    void finder_t::list_visited(const molecule_t& mol, atomvec_t& visited) const
    {
        std::set<int>::const_iterator i = m_visited.begin();
	for( ; i != m_visited.end(); ++i )
	{
	    visited.push_back( morf_t(mol, ATOM, *i) );
	}
    }
    
} // namespace mort

