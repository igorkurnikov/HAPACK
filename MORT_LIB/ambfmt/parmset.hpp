#include <common.hpp>
#include <object.hpp>

namespace mort
{
    using std::vector;

    struct parmset_t
    {
        vector< vector< double > > vdw;
        vector< vector< double > > bond;
        vector< vector< double > > angl;
        vector< vector< double > > tors;
        vector< vector< double > > oops;
        vector< vector< double > > urey;
        vector< vector< double > > ptor;
        vector< vector< double > > strbnd;    
        vector< vector< double > > opbend;
        vector< int > tor2id;

        parmset_t( ) 
            : vdw(2), bond(2), angl(2),
              tors(5),oops(5), urey(2),
              ptor(3),strbnd(4),opbend(2)
        {
        }
    };


    namespace prmtop
    {
        typedef function< void (const morf_t&, vector<double>& ) > parmlister_t;

        struct record_t
        {
            record_t(const hashid_t& idtype, vector< vector<double> >& parm, parmlister_t lister );
	    void operator()(morf_t& o);
            
            vector< vector<double> >& m_parm;
            parmlister_t m_lister;
	    int m_idtype;
        };

        struct list_parm
        {
            list_parm(const hashid_t& pid0, const hashid_t& pid1);

            list_parm(const hashid_t& pid0, const hashid_t& pid1, const hashid_t& pid2);

            list_parm(const hashid_t& pid0, const hashid_t& pid1, const hashid_t& pid2, const hashid_t& pid3, const hashid_t& pid4);

            void operator()( const morf_t& obj, vector<double>& parm ) const;
 
            vector<hashid_t> m_pids;
        };

        struct list_urey
        {
            void operator()( const morf_t& angl, vector<double>& parm ) const;
        };    

        struct list_oops
        {
            void operator()( const morf_t& oops, vector< double >& parm ) const;
        };
        
        struct list_tors
        {
            void operator()( const morf_t& tors, vector< double >& parm ) const;
        };

        struct list_strbnd
        {
            void operator()( const morf_t& angl, vector< double >& parm ) const;
        };
        
        void pack_tors_oops_parm( vector< vector<double> >& torsparm, const vector< vector<double> >& oopsparm, molecule_t& mol );
   
    } // namespace prmtop
    
} // namespace mort

