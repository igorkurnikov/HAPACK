#ifndef GLEAP_MORTSRC_AMBFMT_PARM_HPP
#define GLEAP_MORTSRC_AMBFMT_PARM_HPP


#include <vector>
#include <common.hpp>

namespace mort
{
    class molecule_t;

    namespace prmtop
    {
        enum direction_e { MONO_DIRECTION, BIDIRECTION };

        inline int score( const morf_t& t )
        {
            return 4 - std::count_if( t.atom_begin(), t.atom_end(), sparm_cmper1(NAME, "X") );
        }

        struct by_score
        {
            bool operator()( const morf_t& ti, const morf_t& tj ) const
            {
                return score(ti) > score(tj);
            }
        };

        struct sequence_Match2
        {
            sequence_Match2( const atomvec_t& parmseq, int direct = BIDIRECTION )
                : m_parmseq(parmseq)
            {
	        const molecule_t& ff( parmseq[0].getmol() );

	        atom_t xatom(ff, -1);
		if( atom_t::get( ff, "X", xatom ) )
		{
		    m_xatomid = xatom.absid();
		}
		else
		{
		    m_xatomid = -1;
		}

                m_direction = direct;
            }

            template< typename T >
            bool run( atomiter_t begin, const atomiter_t& end, T input )
            {
                while( begin != end )
                {
                    int id1 = begin->absid();
                    int id2 = input->absid();

                    if( id1 != m_xatomid && id1 != id2 )
                    {
                        return false;
                    }

                    ++begin;
                    ++input;
                }
                
                return true;
            }
            
            bool operator()( const morf_t& parm )
            {
	        if( run( parm.atom_begin(), parm.atom_end(), m_parmseq.begin() ) )
                {
                    return true;
                }
                
                if( m_direction == MONO_DIRECTION )
                {
                    return false;
                }

                return run( parm.atom_begin(), parm.atom_end(), m_parmseq.rbegin() );
            }

            atomvec_t m_parmseq;

	    int m_xatomid;

            int m_direction;
        };
 
        struct sequence_Match
        {
            sequence_Match( const atomvec_t& seq, int direct = BIDIRECTION )
                : m_seq(seq)
            {
                m_direction = direct;
            }

            template< typename T >
            bool run( atomiter_t begin, const atomiter_t& end, T input )
            {
                while( begin != end )
                {
                    string name = begin->get_s(NAME);
                    string type = input->get_s(TYPE);

                    if( name[0] != 'X' && name != type )
                    {
                        return false;
                    }

                    ++begin;
                    ++input;
                }
                
                return true;
            }
            
            bool operator()( const morf_t& parm )
            {
	        if( run( parm.atom_begin(), parm.atom_end(), m_seq.begin() ) )
                {
                    return true;
                }
                
                if( m_direction == MONO_DIRECTION )
                {
                    return false;
                }

                return run( parm.atom_begin(), parm.atom_end(), m_seq.rbegin() );
            }

            atomvec_t m_seq;

            int m_direction;
        };
        

        void parm_atom( const molecule_t* pff, atom_t& a );
        
        void parm_bond( const molecule_t* pff, bond_t& b );

        void parm_angl( const molecule_t* pff, atmvec& as );

        void parm_dihe( const molecule_t* pff, atmvec& as );

        void parm_tor2( const molecule_t* pff, atmvec& as );

        void parm_amoeba_impr( atmvec& as );

        void parm_sander_impr( const molecule_t* pff, atmvec& as, int& max_score );        

        void parm_gbsa(atom_t& atom, const hashid_t& parmid);

    } // namespace prmtop    

    void mark_chain( molecule_t& mol );


} // namespace mort


#endif

