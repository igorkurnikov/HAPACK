#ifndef MORTSRC_OBJECT_ITER_HPP
#define MORTSRC_OBJECT_ITER_HPP

#include <vector>

namespace mort
{
    using std::vector;

    class molecule_t;

    template <typename valu_t > 
    class iter_T
    {
    public:

        /// \brief indicating an random access iterator
        typedef std::random_access_iterator_tag iterator_category;

        typedef std::ptrdiff_t difference_type;
 
        /// \brief the self type
        typedef valu_t value_type;

        typedef valu_t& reference;

        typedef valu_t* pointer;

        typedef iter_T<valu_t> iter_t;

        iter_T( const morf_t& mo, vector<int>::iterator iter )
            : m_ptee(mo),
              m_iter(iter)
        {
        }

        template< typename T>
        iter_T( const iter_T<T>& rhs )
            : m_ptee(rhs.m_ptee),
              m_iter(rhs.m_iter)
        {
        }

        virtual ~iter_T() {}

        iter_t& operator=( const iter_t& rhs )
        {
            m_iter = rhs.m_iter;
            m_ptee = rhs.m_ptee;
            return *this;
        }

        valu_t& operator*()
        {
            m_ptee.m_absid = *m_iter;
            return m_ptee;
        }

        valu_t const& operator*() const
        {
            m_ptee.m_absid = *m_iter;
            return m_ptee;
        }


        valu_t* operator->()  
        { 
            m_ptee.m_absid = *m_iter; 
            return &m_ptee; 
        }
        
        morf_t const* operator->() const 
        { 
            m_ptee.m_absid = *m_iter; 
            return &m_ptee; 
        }

        iter_t& operator++(void) { ++m_iter; return *this; }

        iter_t& operator--(void) { --m_iter; return *this; }

	iter_t operator++(int) { iter_t tmp(*this); ++m_iter; return tmp; }

	iter_t operator--(int) { iter_t tmp(*this); --m_iter; return tmp; }

        iter_t& operator+=( ptrdiff_t dif ) { m_iter += dif; return *this; }

	iter_t& operator-=( ptrdiff_t dif ) { m_iter -= dif; return *this; }

        friend ptrdiff_t operator-( const iter_t& lhs, const iter_t& rhs )
	{
	    return (lhs.m_iter - rhs.m_iter);
	}

	friend bool operator<( const iter_t& lhs, const iter_t& rhs )
	{
	    return (lhs.m_iter < rhs.m_iter);
	}

	friend bool operator>( const iter_t& lhs, const iter_t& rhs )
	{
	    return (lhs.m_iter > rhs.m_iter);
	}

	friend bool operator<=( const iter_t& lhs, const iter_t& rhs )
	{
	    return (lhs.m_iter <= rhs.m_iter);
	}

	friend bool operator>=( const iter_t& lhs, const iter_t& rhs )
	{
	    return (lhs.m_iter >= rhs.m_iter);
	}

	friend bool operator==( const iter_t& lhs, const iter_t& rhs )
	{
	    return (lhs.m_iter == rhs.m_iter);
	}

	friend bool operator!=( const iter_t& lhs, const iter_t& rhs )
	{
	    return (lhs.m_iter != rhs.m_iter);
	}

        friend iter_t operator+( const iter_t& lhs, ptrdiff_t dif )
        {
            iter_t tmp(lhs);
            tmp += dif;
            return tmp;
        }

        friend iter_t operator-( const iter_t& lhs, ptrdiff_t dif )
        {
            iter_t tmp(lhs);
            tmp -= dif;
            return tmp;
        }

    public:

        valu_t m_ptee;

        vector<int>::iterator m_iter;
    };

} // namespace mort

#endif
