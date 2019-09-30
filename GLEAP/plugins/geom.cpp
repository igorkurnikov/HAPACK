#include <common.hpp>
#include <object.hpp>
#include <smarts.hpp>
#include "impose.hpp"

namespace amber
{
    using namespace mort;

    void rotate( morf_t& atom, const numvec& orig, const numvec& axis, double angl )
    {
        assert( atom.cmpid()==ATOM );

        numvec pos = atom.get_v(POSITION);
        pos -= orig;
        mort::rotate( pos, axis, angl );
        pos += orig;
        atom.set_v(POSITION, pos);
    }


    class geom_finder : public finder_t
    {
    public:

        geom_finder( const morf_t& start, const atomvec_t& range, const vector< string >& names, int size )
            : finder_t(ALLOW_MULTI_VISIT), m_range( &range ), m_names( names ), m_size( size )
        {
	    init(start);
        }

        int predict( const atomvec_t& path )
        {
            assert( (int)path.size() <= m_size );

            morf_t atom = path.back();
            if( std::find(m_range->begin(), m_range->end(), atom) == m_range->end() )
            {
                return WRONG_DIRECTION;
            }

            if( atom.get_s(NAME) != m_names[ path.size()-1 ] )
            {
                return WRONG_DIRECTION;
            }

            if( (int)path.size() == m_size )
            {
                return ENOUGH;
            }

            return KEEP_GOING;
        }

    private:

        const  atomvec_t* m_range;
        vector< string > m_names;
        int m_size;
    };

    void find_geom( vector<atomvec_t>& geoms, const atomvec_t& range, const vector<string>& names, int size )
    {
        for( int i=0; i < (int)range.size(); ++i )
        {
            if( range[i].get_s(NAME) != names[0] )
	    {
	        continue;
	    }

            vector< atomvec_t > paths;
            geom_finder f(range[i], range, names, size);

            f.find_all(paths);
            for( int i=0; i < (int)paths.size(); ++i )
	    {
                if( paths[i].size() > 0 )
		{ 
		    geoms.push_back( paths[i] );
		}
	    }
        }
    }

    class related_finder : public finder_t
    {
    public:

        related_finder( const morf_t& root, const morf_t& start )
            : finder_t(NO_MULTI_VISIT), m_start( start )
        {
	    atomvec_t path;
	    path.push_back(root);
	    path.push_back(start);
	    init( path );
        }

        int predict( const atomvec_t& list )
        {
	    morf_t atom = list.back();

            if( atom == list[0] )
            {
                throw std::logic_error( "Error: can't form a clean related atom set, ring found." );
            }

            return KEEP_GOING;
        }    

    private:

        morf_t m_start;
    };

    void find_related( const morf_t& atom, const morf_t& root, atomvec_t& list )
    {
        atomvec_t null;
        related_finder f(root, atom);
        f.find_one(null);

	f.list_visited(atom.getmol(), list);
    /*
        std::cout << "visited ";
	for( int i=0; i < list.size(); ++i )
	{
	    std::cout << list[i].get_i(ID) << " ";
	}

	std::cout << std::endl;
    */
    }

    void impose_dist( const atomvec_t& atoms, double dist )
    {
        assert( atoms.size() == 2 );

        // std::cout << "imposing dist between: " << atoms[0].get_i(ID) << " " << atoms[1].get_i(ID) << std::endl;

	atomvec_t set_0, set_1;
        find_related( atoms[0], atoms[1], set_0 );
        find_related( atoms[1], atoms[0], set_1 );

        numvec v0 = atoms[0].get_v(POSITION);
        numvec v1 = atoms[1].get_v(POSITION);

        double curt = mort::dist( v0, v1 );

        atomvec_t* target;
        numvec  offset(3);

        if( set_0.size() <= set_1.size() )
        {
            target = &set_0;
            offset = (v1-v0) * ( 1 - dist / curt );
        }
        else
        {
            target = &set_1;
            offset = (v0-v1) * ( 1 - dist / curt );
        }

        atomvec_t::iterator ai = target->begin();
        atomvec_t::iterator ae = target->end();
        for( ; ai != ae; ++ai )
        {
            translate( *ai, offset );
        }

    /*
	std::cout << "moving: " << std::endl;
        for( int i=0; i < target->size(); ++i)
	{
	    std::cout << target->at(i).get_i(ID) << " to ";
	    numvec pos = target->at(i).get_v(POSITION);
	    std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	}
        std::cout << std::endl;
    */
    }


    void impose_angl( const numvec& v0, const numvec& v1, const numvec& v2, double angl, atomvec_t& set )
    {
        numvec axis = cross( v0-v1, v2-v1 );
        numvec orig = v1;

        double curt = mort::angl( v0, v1, v2 );
        double offset = angl - curt;
 
        atomvec_t::iterator ai = set.begin();
        atomvec_t::iterator ae = set.end();
        for( ; ai != ae; ++ai )
        {
            //std::cout << "rotating: " << ai->get_s(NAME) << std::endl;
            rotate( *ai, orig, axis, offset );
        }


     /*
	std::cout << "rotating: " << std::endl;
        for( int i=0; i < target->size(); ++i)
	{
	    std::cout << target->at(i).get_i(ID) << " to ";
	    numvec pos = target->at(i).get_v(POSITION);
	    std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	}
        std::cout << std::endl;
     */
    }


    void impose_angl( const morf_t& a0, const morf_t& a1, const morf_t& a2, double angl, atomvec_t& set )
    {
        numvec v0 = a0.get_v(POSITION);
        numvec v1 = a1.get_v(POSITION);
        numvec v2 = a2.get_v(POSITION);
        impose_angl( v0, v1, v2, angl, set );
    }


    void impose_angl( const atomvec_t& atoms, double angl )
    {
        assert( atoms.size()==3 );
        //std::cout << "imposing angl between: ";
        //std::cout << atoms[0].get_s(NAME) << "-";
        //std::cout << atoms[1].get_s(NAME) << "-";
        //std::cout << atoms[2].get_s(NAME) << std::endl;


	atomvec_t set_0, set_2;

        find_related( atoms[0], atoms[1], set_0 );
        find_related( atoms[2], atoms[1], set_2 );

        numvec rot(4);
        if( set_0.size() <= set_2.size() )
        {
            impose_angl( atoms[2], atoms[1], atoms[0], angl, set_0 );
        }
        else
        {
            impose_angl( atoms[0], atoms[1], atoms[2], angl, set_2 );
        }
    }

    void impose_tors( const numvec& v0, const numvec& v1, const numvec& v2, const numvec& v3, double tors, atomvec_t& set )
    {
        numvec axis = normalcpy( v1 - v2 );
        numvec orig = v1;
        double curt = mort::tors( v0, v1, v2, v3 );
        double offset = curt - tors;
 
        atomvec_t::iterator ai = set.begin();
        atomvec_t::iterator ae = set.end();
        for( ; ai != ae; ++ai )
        {
            //std::cout << "rotating: " << ai->get_s(NAME) << std::endl;
            rotate( *ai, orig, axis, offset );
        }

     /*
	std::cout << "rotating: " << std::endl;
        for( int i=0; i < target->size(); ++i)
	{
	    std::cout << target->at(i).get_i(ID) << " to ";
	    numvec pos = target->at(i).get_v(POSITION);
	    std::cout << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
	}
        std::cout << std::endl;
     */
    }

    void impose_tors( const morf_t& a0, const morf_t& a1, const morf_t& a2, const morf_t& a3, double tors, atomvec_t& set)
    {

        //std::cout << "imposing tors between: ";
        //std::cout << a0.get_s(NAME) << "-";
        //std::cout << a1.get_s(NAME) << "-";
        //std::cout << a2.get_s(NAME) << "-";
        //std::cout << a3.get_s(NAME) << std::endl;

        numvec v0 = a0.get_v(POSITION);
        numvec v1 = a1.get_v(POSITION);
        numvec v2 = a2.get_v(POSITION);
        numvec v3 = a3.get_v(POSITION);

        impose_tors( v0, v1, v2, v3, tors, set );
    }

    void impose_tors( const atomvec_t& atoms, double tors, bool middle_start )
    {
        assert( atoms.size() == 4 );


	atomvec_t set_1, set_2;

        if( middle_start )
        {
            find_related( atoms[1], atoms[2], set_1 );
            find_related( atoms[2], atoms[1], set_2 );
        }
        else
        {
            find_related( atoms[0], atoms[1], set_1 );
            find_related( atoms[3], atoms[2], set_2 );
        }

        if( set_1.size() <= set_2.size() )
        {
            impose_tors( atoms[3], atoms[2], atoms[1], atoms[0], tors, set_1 );
        }
        else
        {
            impose_tors( atoms[0], atoms[1], atoms[2], atoms[3], tors, set_2 );
        }
    }

    void impose_dist( atomvec_t& atoms, const vector<string>& inter )
    {
        vector< atomvec_t > geoms;
        find_geom( geoms, atoms, inter, 2 );

        for( int i=0; i < (int)geoms.size(); ++i )
        {
            impose_dist( geoms[i], atof( inter[2].c_str() ) );
        }
    }

    void impose_angl( atomvec_t& atoms, const vector<string>& inter )
    {
        vector< atomvec_t > geoms;

        find_geom( geoms, atoms, inter, 3 );
        
        for( int i=0; i < (int)geoms.size(); ++i )
        {
            impose_angl( geoms[i], atof( inter[3].c_str() ) );
        }
    }

    void impose_tors( atomvec_t& atoms, const vector<string>& inter )
    {
        vector< atomvec_t > geoms;

        find_geom( geoms, atoms, inter, 4 );

        for( int i=0; i < (int)geoms.size(); ++i )
        {
            impose_tors( geoms[i], atof( inter[4].c_str() ), true );
        }
    }

}  // namespace mort   



