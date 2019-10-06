#ifndef MORTSRC_AMBFMT_ENERGEE_HPP
#define MORTSRC_AMBFMT_ENERGEE_HPP

#include <vector>
#include <common.hpp>
#include <object.hpp>
#include "sff.h"
#include "parmset.hpp"
#include "exclude.hpp"

namespace mort
{
    using std::vector;

    typedef PARMSTRUCT_T nabparm_t;

    struct conninfo_t
    {
        vector< vector<int> > id3;
        vector<int> typ;

        int size()
        {
            for( unsigned int i=0; i < id3.size(); ++i )
            {
                assert( id3[i].size()==typ.size() );
            }
            return typ.size();
        }

        void push_dihe( dihe_t tors )
        {
            assert( tors.natom()==4 );

            atom_t a1 = atom_1st( tors );
            atom_t a2 = atom_2nd( tors );
            atom_t a3 = atom_3rd( tors );
            atom_t a4 = atom_4th( tors );
        
            if( a3.absid() == 0 || a4.absid() == 0 )
            {
                std::swap( a1, a4 );
                std::swap( a2, a3 );
            }

            if( id3.size()==0 )
            {
                id3.resize(4);
            }

            id3[0].push_back( a1.get_i(OFF3) );
            id3[1].push_back( a2.get_i(OFF3) );
            if( bond_t::has( a1, a4 ) || angl_t::has( a1, a4 ) || has_prev(tors)  )
            {
                id3[2].push_back( a3.get_i(OFF3) * (-1) );
            }
            else
            {
                id3[2].push_back( a3.get_i(OFF3) );
            }
        
            id3[3].push_back( a4.get_i(OFF3) );
            typ.push_back( tors.get_i(TYPEID) );

            size();
        }

        void push_impr( impr_t im )
        {
            if( id3.size()==0 )
            {
                id3.resize(4);
            }

            atom_range as = im.atoms();
            id3[0].push_back( as[0].get_i(OFF3) );
            id3[1].push_back( as[1].get_i(OFF3) );
            id3[2].push_back( as[2].get_i(OFF3) * (-1) );
            id3[3].push_back( as[3].get_i(OFF3) * (-1) );
            typ.push_back( im.get_i(TYPEID) );
            size();
        }

        void push( morf_t& mo )
        {
            if( mo.cmpid()==TORS )
            {
                push_dihe(mo);
            }
            else if( mo.cmpid()==OOPS )
            {
                push_impr(mo);
            }
            else
            {
                if( id3.size()==0 )
                {
                    id3.resize( mo.natom() );
                }

                atomiter_t ai = mo.atom_begin();
                for(int i=0; ai != mo.atom_end(); ++ai, ++i)
                {
                    id3[i].push_back( ai->relid()*3 );
                }

                typ.push_back( mo.get_i(TYPEID) );
            }
        }
    };



    class energee_t
    {
    public:

        energee_t( molecule_t& m ) { m_pmol = &m; }

        virtual ~energee_t() {}

        void assignparm( const molecule_t& ffp, int ftype=AMBER );

        nabparm_t& getnabparm() { return m_nabparm; }

        const nabparm_t& getnabparm() const { return m_nabparm; }


        parmset_t m_parmset;
        nabparm_t m_nabparm;
        excl_t    m_excl;

    private:

        molecule_t* m_pmol;
        conninfo_t m_bondh;
        conninfo_t m_bondo;
        conninfo_t m_anglh;
        conninfo_t m_anglo;
        conninfo_t m_torsh;
        conninfo_t m_torso;

        vector<int> m_joint;
        vector<int> m_rotat;
        vector<int> m_resdp;
        vector<int> m_usize;

        vector<int> m_exsize;
        vector<int> m_exlist;
        vector<int> m_vdwidx;
        vector<int> m_14size;
        vector<int> m_14list;
        

        vector<double> m_gvdw;
        vector<double> m_solty;
        vector<double> m_vdwcn1;
        vector<double> m_vdwcn2;
        vector<double> m_charge;

        string m_atomname;
        string m_resdname;
        string m_atomtype;
        string m_atomtree;

    };


} // namespace mort


#endif
