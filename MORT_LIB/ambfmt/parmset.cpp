#include "parmset.hpp"

namespace mort
{
    
    namespace prmtop
    {
        record_t::record_t(const hashid_t& idtype, vector< vector<double> >& parm, parmlister_t lister )
            : m_parm( parm ), m_lister(lister)
        {
	    m_idtype = idtype;
        }
            
        void record_t::operator()(morf_t& obj)
        {
            vector<double> parm;
            m_lister(obj, parm);

            if( parm.size()==0 ) obj.set_i(m_idtype, -1);                

            int i=0;
            for(; i < (int)m_parm[0].size(); ++i)
            {
                bool equil = true;
                int j=0;
                while( j < (int)parm.size() && equil )
                {                        
                    if( parm[j]>m_parm[j][i]+1e-6 || parm[j]<m_parm[j][i]-1e-6 )
                        equil = false;
                    j++;
                }
                    
                if( equil ) break;
            }
                
            if( i == (int)m_parm[0].size() )
            {
                for( int j=0; j < (int)parm.size(); ++j )
                {
                    m_parm[j].push_back( parm[j] );
                }

                
                obj.set_i(m_idtype, m_parm[0].size() );
            }
            else
            {
                obj.set_i(m_idtype, i + 1 );
            }
        }

	list_parm::list_parm(const hashid_t& pid0, const hashid_t& pid1)
            : m_pids(2)
        {
	    m_pids[0] = pid0;
	    m_pids[1] = pid1;
	}

	list_parm::list_parm(const hashid_t& pid0, const hashid_t& pid1, const hashid_t& pid2)
            : m_pids(3)
        {
	    m_pids[0] = pid0;
	    m_pids[1] = pid1;
            m_pids[2] = pid2;
	}

	list_parm::list_parm(const hashid_t& pid0, const hashid_t& pid1, const hashid_t& pid2,
            const hashid_t& pid3, const hashid_t& pid4)
            : m_pids(5)
        {
	    m_pids[0] = pid0;
	    m_pids[1] = pid1;
            m_pids[2] = pid2;
            m_pids[3] = pid3;
            m_pids[4] = pid4;
	}

        double get_d_or_i(const morf_t& obj, const hashid_t& parmid)
	{
	     double dvalue;
	     if( obj.get_d(parmid, dvalue) )
	     {
	         return dvalue;
             }
	     
	     int ivalue;
	     if( obj.get_i(parmid, ivalue) )
	     {
	         return ivalue;
             }
	     
             throw std::runtime_error( "Error: can not get double or integer parameter " + unhash(parmid) );
	}


        void list_parm::operator()( const morf_t& obj, vector<double>& parm ) const
        {
            parm.resize( m_pids.size() );
            for( unsigned int i=0; i < m_pids.size(); ++i )
            {
                parm[i] = get_d_or_i( obj, m_pids[i] );
            }
        }

        void list_urey::operator()( const morf_t& angl, vector< double >& parm ) const
        {
            parm.clear();
            
            numvec urey = angl.get_v(UREY);
            if( urey[0] != 0.0 )
            {
                parm.push_back( urey[0] );
                parm.push_back( urey[1] );
            }
        }

        void list_oops::operator()( const morf_t& oops, vector< double >& parm ) const
        {
            parm.clear();
                
            numvec opbend = oops.get_v(OPBEND);    
            if( opbend[0] != 0.0 )
            {
                parm.push_back( opbend[0] );
            }
        }

        void list_tors::operator()( const morf_t& tors, vector< double >& parm ) const
        {
            parm.clear();
      
            if( tors.get_d(FORCE) != 0.0 )
            {
                parm.push_back( tors.get_d(FORCE) );
                parm.push_back( tors.get_i(PERIOD) );
                parm.push_back( tors.get_d(EQUIL) );
            }
        }

        void list_strbnd::operator()( const morf_t& angl, vector< double >& parm ) const
        {
            parm.clear();
            numvec strbnd = angl.get_v(STRBND);

            if( strbnd[0] != 0.0 )
            {
                parm.push_back( strbnd[0] );
                parm.push_back( angl.get_d(EQUIL) );
                parm.push_back( bond_t::get(atom_1st(angl), atom_2nd(angl)).get_d(EQUIL) );
                parm.push_back( bond_t::get(atom_3rd(angl), atom_2nd(angl)).get_d(EQUIL) );
            }
        }        

        void pack_tors_oops_parm( vector< vector< double > >& torsparm, const vector< vector< double > >& oopsparm, molecule_t& mol )
        {
            assert( torsparm.size() == 5 && oopsparm.size() == 5 );

            int baseid = torsparm[0].size();
        
            impriter_t oops = mol.impr_begin();
            for( ; oops != mol.impr_end(); ++oops )
            {
                oops->set_i(TYPEID, oops->get_i(TYPEID)+baseid );
            }

            for( unsigned int i=0; i < torsparm.size(); ++i )
            {
                torsparm[i].insert( torsparm[i].end(), oopsparm[i].begin(), oopsparm[i].end() );
            }
        }
        
    } // namespace prmtop
    
} // namespace mort



    
