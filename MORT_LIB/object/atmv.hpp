#ifndef MORT_OBJECT_ATOMVEC_HPP
#define MORT_OBJECT_ATOMVEC_HPP

#include <vector>
#include "root.hpp"
#include "atom.hpp"
#include "bond.hpp"

namespace mort
{
    using std::vector;
    /// \brief atmvec is an atom container, uses vector to store atoms.
    /// \ingroup objdef
    class atmvec : public std::vector<atom_t>, public root_t
    {
    public:

        atmvec();

        atmvec( const morf_t& a0 );

        atmvec( int size, const atom_t& a0 );

        atmvec( const atom_range& range );

        atmvec( const atom_t& a0, const atom_t& a1 );
        
        atmvec( const atom_t& a0, const atom_t& a1, const atom_t& a2 );

        atmvec( const atom_t& a0, const atom_t& a1, 
		const atom_t& a2, const atom_t& a3 );

        atmvec( const atom_t& a0, const atom_t& a1, const atom_t& a2, 
	        const atom_t& a3, const atom_t& a4 );

        atmvec( const atom_t& a0, const atom_t& a1, const atom_t& a2, 
	        const atom_t& a3, const atom_t& a4, const atom_t& a5 );
        
        virtual ~atmvec();

        atmvec& operator=( const atmvec& rhs );
    
        void push_atom( const morf_t& p );
        
        int count( const atom_t& atom ) const;

	bool is_ring() const;

        bool is_arom() const;

        std::vector<atom_t>::const_iterator find(const atom_t& a) const
        {
            return std::find( begin(), end(), a);
        }
    };
    

    struct atom_pusher
    {
    protected:
        atmvec* subset;

    public:
        typedef atmvec container_type;

        explicit
        atom_pusher( atmvec& __x) : subset( &__x ) { }

        atom_pusher&
        operator=( const morf_t& o )
	{
	    if( o.cmpid() == ATOM )
	    {
                subset->push_back( o );
	    }
	    else
	    {
	        subset->push_atom( o );
	    }
            return *this;
        }

        atom_pusher&
        operator*()
        { return *this; }

        atom_pusher&
        operator++()
        { return *this; }

        atom_pusher
        operator++(int)
        { return *this; }
    };

    /// \brief print the name of atoms in a subset 
    std::ostream& operator<<( std::ostream& stream, const atmvec& path );

    atmvec find_path( const atom_t& begin, const atom_t& end );

    typedef shared_ptr<atmvec> atmvec_ptr;
        
    typedef atmvec atomvec_t;

    class bndvec : public std::vector<bond_t>, public root_t
    {
    public:

        bndvec() {}
        
        virtual ~bndvec() {}
    };    

    typedef shared_ptr<bndvec> bndvec_ptr;


} // namespace mort


#endif

