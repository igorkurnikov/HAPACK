#ifndef GLEAP_MORTSRC_OBJECT_RANG_HPP
#define GLEAP_MORTSRC_OBJECT_RANG_HPP


namespace mort
{

    /// \brief morange_t provide a convinent way of pointer's random access.
    /// \ingroup objdef
    /// it works just like a normal array.
    template <typename iter_t, typename valu_t> 
    class range_T
    {
    public:

        /// \brief constructor
        /// \params bgn the bgn iterator 
        /// \params end   the ending iterator
        range_T( const iter_t& bgn, const iter_t& end )
            : m_bgn( bgn ), m_end ( end ) 
        {
        }

        template<typename T1, typename T2>
        range_T( const range_T<T1, T2>& rhs )
            : m_bgn(rhs.begin()), m_end(rhs.end())
        {
        }

        iter_t begin() const { return m_bgn; }
        
        iter_t end()   const { return m_end; }
        
        /// return number of objects in this range
        int size() const { return m_end - m_bgn; }
        
        /// return object at certain position
        valu_t  at( int id ) const
        {
	    return *(m_bgn+id);
        }
        
        /// return object at certain position, 
        valu_t operator[]( int id ) const
        {
            return *(m_bgn+id);
        }
        
    private:

        iter_t m_bgn;
        iter_t m_end;
    };

} // namespace mort


#endif

