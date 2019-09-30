#ifndef MORTSRC_CAPBOX_SOLUTE_HPP
#define MORTSRC_CAPBOX_SOLUTE_HPP

#include <vector>
#include <common.hpp>
#include <object.hpp>

namespace mort
{
    class molecule_t;

    class solute_i
    {
    public:

        virtual ~solute_i() {}

        // Note: it returns the natom of solute. Not the natom of
        // the whole system. the later can be obtained by size of vdwr
        virtual int natom() const = 0;

        virtual std::vector<double>& getcord() = 0;

        virtual std::vector<double>& getvdwr() = 0;

        // Note: when adding new atom, set vdwr property
        virtual void insert_resd( const molecule_t& svt, const numvec& shift, int idx ) = 0;

        // sync the coordinate after making change to it.
        virtual void finish() = 0;
    };

    class molecule_t;

    class mort_solute : public solute_i
    {
    public:

        mort_solute( molecule_t& m )
        {
            m_pmol = &m;
            m.set_i( SOLUTE_NATOM, m.natom() );
            m.set_i( SOLUTE_NRESD, m.nresd() );
            set_vdwr( m );
        }

        virtual ~mort_solute() {}

        virtual int natom() const { return m_pmol->get_i(SOLUTE_NATOM); }

        virtual std::vector<double>& getcord()
        {
            return get_vvec( *m_pmol, ATOM, POSITION );
        }

        virtual std::vector<double>& getvdwr()
        {
            return get_dvec( *m_pmol, ATOM, VDWR );
        }

        virtual void insert_resd( const molecule_t& svt, const numvec& sft, int rid)
        {
            morf_t r = merge( *m_pmol, svt, rid );
            translate( r, sft );
            set_vdwr( r );
        }
 
        virtual void finish() {}

    private:

        molecule_t* m_pmol;
    };

    // translate solute so that region center is origin
    //        params:
    //            solute 
    //            vdwr: Van der Waals radius is considered
    //        returns: the extent of the region
    numvec regionlize( solute_i& solute, bool vdwr );

    // calculate the buffer size of octrahedron solvent.
    // there is a minimum required buffer size for octrahedron
    // otherwise part of the solute might will be exposed to vacuum
    //        params:
    //            desired: user input
    //        return value: the proper buffer size, must be >= dersired.
    double octbufsize( solute_i& solute, double desired );

    // translate solute so that geometrical center is on origin.
    void centralize( solute_i& solute );

    // rotate the solute so that the longest distance between atoms are
    // on axis. used by octrahedron
    void rotatelong( solute_i& solute );

    // ewald rotation
    void rotatewald( solute_i& solute );

    int get_solute_nresd( const molecule_t& m );

    int get_solute_natom( const molecule_t& m );

    int is_sequential( const molecule_t& m );

} // namespace mort

#endif


