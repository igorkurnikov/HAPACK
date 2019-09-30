#ifndef MORTSRC_CAPBOX_ADDIONS_HPP
#define MORTSRC_CAPBOX_ADDIONS_HPP

#include <object.hpp>
#include "solute.hpp"

namespace mort
{

    class ionee_i
    {
    public: 

        virtual ~ionee_i() {}

        virtual int natom()=0;
        virtual int nresd()=0;
        virtual int solute_nresd() = 0;
        virtual int solute_natom() = 0;

        virtual vector<double>& getcord() = 0;
        virtual vector<double>& getchrg() = 0;
        virtual vector<double>& getvdwr() = 0;
	virtual vector<string>& gettype() = 0;

        virtual void insert_resd( const molecule_t& ion, const numvec& pos, int idx) = 0;
        virtual void remove_resd( int rid ) = 0;
    };


    class mort_ionee : public ionee_i
    {
    public:

        mort_ionee( molecule_t& m ) 
        {
            m.cleanup();

            m_pmol = &m;

            m_natom = m.natom();
            m_nresd = m.nresd();

            m_solute_nresd = get_solute_nresd( m );
            if( m_solute_nresd == m_nresd )
            {
                m_solute_natom = m_natom;
            }
            else
            {
                int solvent_size = (m.resd_end()-1)->natom();
                int solvent_numb = m_nresd - m_solute_nresd;
                m_solute_natom = m_natom - solvent_size*solvent_numb;
            }

            set_vdwr( m );

            if( is_sequential(m) )
            {
                m_pcord = &get_vvec(m, ATOM, POSITION);
                m_pchrg = &get_dvec(m, ATOM, PCHG);
                m_pvdwr = &get_dvec(m, ATOM, VDWR);
		m_ptype = &get_svec(m, ATOM, TYPE);
            }
            else
            {
                m_cord.resize( 3*m_natom );
                m_chrg.resize( m_natom );
                m_vdwr.resize( m_natom );
		m_type.resize( m_natom );

                resditer_t ri = m_pmol->resd_begin();
                resditer_t re = m_pmol->resd_end();
                int i=0;
                for( ; ri != re; ++ri )
                {
                    atomiter_t ai = ri->atom_begin();
                    atomiter_t ae = ri->atom_end();
                    for( ; ai != ae; ++ai )
                    {
                        numvec pos = ai->get_v(POSITION);
                        m_cord[3*i  ] = pos[0];
                        m_cord[3*i+1] = pos[1];
                        m_cord[3*i+2] = pos[2];
                        m_vdwr[i] = ai->get_d(VDWR);
                        m_chrg[i] = ai->get_d(PCHG);
			m_type[i] = ai->get_s(TYPE);
                        ++i;
                    }
                }

                if( i != m_natom )
                {
                    throw std::runtime_error( "Error: there are atoms do not belong to any residue" );
                }

                m_pcord = &m_cord;
                m_pchrg = &m_chrg;
                m_pvdwr = &m_vdwr;
		m_ptype = &m_type;
            }
        }

        virtual ~mort_ionee()
        {
        }

        vector<double>& getcord()
        {
            return *m_pcord;
        }

        vector<double>& getchrg()
        {
            return *m_pchrg;
        }

        vector<double>& getvdwr()
        {
            return *m_pvdwr;
        }

        vector<string>& gettype()
        {
            return *m_ptype;
        }

        int natom()
        {
            return m_natom;
        }

        int nresd()
        {
            return m_nresd;
        }

        int solute_nresd()
        {
            return m_solute_nresd;
        }

        int solute_natom()
        {
            return m_solute_natom;
        }

        void insert_resd( const molecule_t& ion, const numvec& sft, int rid )
        {
            morf_t r = merge( *m_pmol, ion, rid );
            translate( r, sft );
        }

        void remove_resd( int rid )
        {
            morf_t r(*m_pmol, RESD, rid);
            m_pmol->remove_resd( r );
        }


    private:

        molecule_t* m_pmol;

        int m_natom;
        int m_nresd;
        int m_solute_natom;
        int m_solute_nresd;

        vector<double>* m_pcord;
        vector<double>* m_pchrg;
        vector<double>* m_pvdwr;
	vector<string>* m_ptype;

        vector<double> m_cord;
        vector<double> m_chrg;
        vector<double> m_vdwr;
	vector<string> m_type;
    };

    void addions_core( ionee_i& ionee, const molecule_t& ion, int nion, double shlext, double resolution );
 
    void addions( molecule_t& m, const molecule_t& ion, int nion, double shlext, double res );
}

 



#endif
