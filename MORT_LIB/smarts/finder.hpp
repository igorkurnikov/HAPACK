#ifndef MORT_SUBSET_STATE_HPP
#define MORT_SUBSET_STATE_HPP

#include <object.hpp>

namespace mort
{
    class morf_t;

    class finder_t 
    {
    public:

        finder_t(int strategy);

	virtual ~finder_t();

        void init( );

        void init(const morf_t& p);
 
        void init(const atomvec_t& path);

        void init(const atomiter_t& begin, const atomiter_t& end);

        bool dead();

        void back_one();

	int step();

        morf_t next_atom();
 
        bool find_one(atomvec_t& path);

	bool find_all(vector<atomvec_t>& paths);

        virtual int predict(const atomvec_t& path) = 0;
 
        virtual void list_next(const atomvec_t& path, atomvec_t& candidate);

        void list_visited(const molecule_t& mol, atomvec_t& visited) const;
   protected:

        atomvec_t m_path;
        
        vector<atomvec_t> m_nexts;

	bool m_allow_multi_visit;

        std::set<int> m_visited;
    };

} // namespace mort

#endif


