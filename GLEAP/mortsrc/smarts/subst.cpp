#include <object.hpp>
#include "subst.hpp"
#include "atomexpr.hpp"


namespace mort
{
    typedef boost::shared_ptr< atomexpr_t > atomexpr_ptr;

    using boost::any_cast;

    bool match_atom( const morf_t& satom, const morf_t& matom )
    {
        any pattern;
        if( satom.get_a(PATTERN, pattern) )
        {
            return any_cast<atomexpr_ptr>(pattern)->match( matom );
        }
        
        if( matom.get_a(PATTERN, pattern) )
        {
            return any_cast<atomexpr_ptr>(pattern)->match( satom );
        }

        return satom.get_i(ELEMENT) == matom.get_i(ELEMENT);
    }
    
    bool match_bond( const morf_t& sbond, const morf_t& mbond )
    {
        if( sbond.get_i(ORDER) == 99 || mbond.get_i(ORDER) == 99 )
        {
            return true;
        }

        if( sbond.has_i(AROM) && sbond.get_i(AROM) )
        {
            return atom_1st(mbond).get_i(AROM) && atom_2nd(mbond).get_i(AROM);
        }

        return sbond.get_i(ORDER) == mbond.get_i(ORDER);
    }

    subst_finder::subst_finder(const morf_t& atom, const molecule_t& sub) 
        :  finder_t(ALLOW_MULTI_VISIT), m_sub(sub)
    { 
        init(atom);
    }

    subst_finder::~subst_finder() 
    {
    }

    int subst_finder::predict(const atomvec_t& path)
    {
        if( (int)path.size() == m_sub.natom() )
	{
	    return ENOUGH;
	}

        return KEEP_GOING;
    }

    void subst_finder::find_next(atomiter_t begin, atomiter_t end, const morf_t& s_atom, const atmvec& path, atmvec& next )
    {
        atomiter_t mi = begin; 
        for( ; mi != end; ++mi )
        {
            if( path.count(*mi) ==0 && match( s_atom, *mi ) )
     	    {
	        next.push_back( *mi );
	    }
	}
    }

    void subst_finder::list_next(const atomvec_t& path, atomvec_t& next)
    {
        assert( m_sub.natom() > (int)m_path.size() );
        atomiter_t s_atom = m_sub.atom_begin() + m_path.size();

        int id = -1;
        
        atomiter_t s_nbr = s_atom->atom_begin();
        for( ; s_nbr != s_atom->atom_end(); ++s_nbr )
        {
            id = std::find(m_sub.atom_begin(), s_atom, *s_nbr) - m_sub.atom_begin();

            if( id != (int)m_path.size() )
            {
                break;
            }
        }      
        
        if( s_nbr == s_atom->atom_end() )
        {
            find_next( m_path[0].getmol().atom_begin(), m_path[0].getmol().atom_end(), 
                       *s_atom, path, next );
        }
        else
        {
            morf_t& m_nbr = m_path[ id ];
            find_next( m_nbr.atom_begin(), m_nbr.atom_end(), 
                       *s_atom, path, next );
        }
    }

    bool subst_finder::match( const morf_t& s_atom, const morf_t& m_atom )
    {
        if( ! match_atom( s_atom, m_atom ) )
        {
            return false;
        }

        atomiter_t s_nbr = s_atom.atom_begin();
        atomiter_t s_end = s_atom.atom_end();
        for( ; s_nbr != s_end; ++s_nbr )
        {
            if( *s_nbr < s_atom )
            {
                atom_t& m_nbr = m_path[ s_nbr->absid() ];
                
                if( ! bond_t::has( m_nbr, m_atom ) )
                {
                    return false;
                }
                
                bond_t s_bond = bond_t::get( s_atom, *s_nbr );
                bond_t m_bond = bond_t::get( m_atom, m_nbr );
                
                if( ! match_bond( s_bond, m_bond ) )
                {
                    return false;
                }
            }
        }

        return true;
    }

    bool find_subst( const morf_t& atom, const molecule_t& sub, atomvec_t& path)
    {
        if( !match_atom( *sub.atom_begin(), atom ) ) 
        {
            return false;
        }

        // std::cout << "first check passed!\n";
        subst_finder f( atom, sub );
    
        return f.find_one(path);
    }

    bool has_subst( const morf_t& atom, const molecule_t& sub )
    {
        if( atom.getmol().natom() < sub.natom() )
	    return false;

        atomvec_t path;
        return find_subst( atom, sub, path);
    }

}  // namespace mort


