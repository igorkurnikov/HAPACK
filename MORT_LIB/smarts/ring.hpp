#ifndef MORT_SUBSET_RING_HPP
#define MORT_SUBSET_RING_HPP

#include "finder.hpp"

namespace mort
{
    class morf_t;
    class atmvec;

    /// \defgroup ringfun Fucntions for ring detection
    /// \ingroup objfun
    /// This group contains functions related to ring detection
    /// @{
    ///

    enum ringfinder_type { WANT_ALL, WANT_ONE };

    class ring_finder : public finder_t
    {
    public:

        ring_finder(const morf_t& o, int type=WANT_ALL, int size=0);

        virtual ~ring_finder();
        
        virtual int predict(const atmvec& path);

    private:
        
        int m_ringsize;
        
    };
    
    /// \brief test if an atom is in ring
    bool has_ring(const morf_t& p);

    /// \brief test if an atom is in a certain sized ring.
    bool has_ring(const morf_t& a, int size);

    /// \brief find any ring contain the given atom
    bool find_ring(const morf_t& p, atmvec& ring);
    
    /// \brief find a certain sized ring contain the given atom
    /// \param size ring size, must > 0.
    /// returns empty atomvec.
    bool find_ring(const morf_t& a, int size, atmvec& ring);
 
    void find_ring( const morf_t& a, vector<atmvec>& rings );
   
    
    // \brief find an aromatic ring contain the given atom
    ///
    /// return NULL if such a ring doesn't exist.
    bool find_arom(const morf_t& a, atomvec_t& ring);
    
    bool is_arom(const morf_t& a );

    /** @} */ // end of group ringfun

}

#endif
