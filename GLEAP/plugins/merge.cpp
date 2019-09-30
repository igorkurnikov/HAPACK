#include <stdexcept>

#include <object.hpp>
#include <boost/algorithm/string.hpp>
#include "merge.hpp"
#include "impose.hpp"

namespace amber
{

    struct less_nbrs
    {
        bool operator()( const morf_t& lhs, const morf_t& rhs )
	{
	    return lhs.natom() < rhs.natom();
	}
    };

    morf_t get_other_nbr( const morf_t& a1, const morf_t& a2 )
    {
        assert( a1.natom()>=2 );

	morf_t n1 = a1.atoms()[0];
	morf_t n2 = a1.atoms()[1];

	return n1==a2 ? n2 : n1;
    }

    morf_t get_nbr( const morf_t& a )
    {
        assert( a.natom() > 0 );
	    
        atomiter_t ai = a.atom_begin();
        atomiter_t ae = a.atom_end();
        for( ; ai != ae; ++ai )
        {
            if( ai->get_i(ELEMENT) != HYDROGEN )
                return *ai;
        }

        return a.atoms()[0];
    }

    numvec get_ang_ref( const morf_t& a )
    {
        assert( a.natom() > 0 );

        numvec orig = a.get_v(POSITION);

        std::vector<numvec> refs;
        atomiter_t ai = a.atom_begin();
        atomiter_t ae = a.atom_end();
        for( ; ai != ae; ++ai )
        {
            numvec pos = ai->get_v(POSITION);
            refs.push_back( normalcpy(pos - orig) );
        }

        assert( refs.size()==2 || refs.size()==3 );

        numvec tmp = rotcpy( refs[0], refs[1], 120.0 );
        if( refs.size()==2 || dis2(tmp,refs[2]) > 0.25)
        {
            return orig+tmp;
        }

        return orig+rotcpy(refs[0], refs[1], 240.0);

    }

    numvec get_tor_ref( const morf_t& a )
    {
        assert( a.natom() > 0 );

        std::vector<numvec> refs;
        atomiter_t ai = a.atom_begin();
        atomiter_t ae = a.atom_end();
        for( ; ai != ae; ++ai )
        {
            numvec pos = ai->get_v(POSITION);
            if( ai->get_i(ELEMENT)!=HYDROGEN )
            {
                refs.push_back( ai->get_v(POSITION) );
            }
        }

        if( refs.size()==0 )
            return a.atoms()[0].get_v(POSITION);

        if( refs.size()==1 )
            return refs[0];

        numvec orig = a.get_v(POSITION);
        numvec v0 = normalcpy( refs[0] - orig );
        numvec v1 = normalcpy( refs[1] - orig );
        numvec vm = (v0 + v1)*0.5;
        return orig + vm;
    }




    atomvec_t get_merge_pts( morf_t& prev, morf_t& tail )
    {

        atomvec_t atoms;

	atoms.push_back(tail);  

        if( tail.natom()==0 )
	    return atoms;

	if( tail.natom()==1 )
	{
	    morf_t atm_2 = tail.atoms()[0];
            atoms.push_back( atm_2 );

	    if( atm_2.natom()==1 )
	        return atoms;

	    morf_t atm_3 = get_other_nbr( atm_2, tail );
	    atoms.push_back( atm_3 );
	    return atoms;
	}
	    
        assert( tail.natom()>=2 );
	
	morf_t atm_2 = tail.atoms()[0];
	morf_t atm_3 = tail.atoms()[1];
	if( atm_2.get_s(NAME)=="CA" )
	{
	    atoms.push_back( atm_2 );
	    atoms.push_back( atm_3 );
	}
	else
	{
	    atoms.push_back( atm_3 );
	    atoms.push_back( atm_2 );
	}

        return atoms;
    }

    numvec get_merge_pos( const atomvec_t& atoms )
    {
        static const double LENGTH = 1.7;
	static const double DEGRAD = 120.0;
	static const double THETA  = DEGRAD*M_PI/180.0;

        vector<numvec> crds;
	for( unsigned int i=0; i < atoms.size(); ++i )
	{
	    crds.push_back( atoms[i].get_v(POSITION) );

            //std::cout << "name,pos: " << atoms[i].get_s(NAME) << " " << crds[i][0] << " " << crds[i][1] << " " << crds[i][2] << std::endl;
	}

        assert( crds.size()>0 );

        if( crds.size()==1 )
	{
	    return makevec( crds[0][0] + LENGTH*cos(THETA), crds[0][1]+LENGTH*sin(THETA), crds[0][2] );
        }

	numvec v10 = crds[1] - crds[0];
        numvec v20 = crds.size()==2 ? makevec(1.0, 0.0, 0.0) : (crds[2]-crds[0]);
	numvec axs = normalcpy( cross(v20, v10) );
	numvec v30 = rotcpy( v10, axs, DEGRAD );
	return v30 + crds[0];
    }



    void place( morf_t& prev, morf_t& curt, morf_t& tail, morf_t& head )
    {
        assert( prev.cmpid()==RESD && curt.cmpid()==RESD );

        atomvec_t atms;
        for( int i=curt.relid(); i < curt.getmol().nresd(); ++i )
        {
            atms.push_atom( curt.getmol().resds()[i] );
        }

        numvec hpos = get_merge_pos( get_merge_pts(prev,tail) );


        numvec offset = hpos - head.get_v(POSITION);
       // translate( curt, offset );
        for(unsigned int i=0; i < atms.size(); ++i )
        {
            numvec v = atms[i].get_v(POSITION);
            atms[i].set_v(POSITION, v+offset);
        }


        numvec tref = get_nbr( tail ).get_v(POSITION);
        numvec tpos = tail.get_v(POSITION);

        numvec ref1 = get_ang_ref( head );
        impose_angl( tpos, hpos, ref1, 0.0, atms );

        numvec ref2 = get_tor_ref( head );
        impose_tors( tref, tpos, hpos, ref2, 180.0, atms );
    }

    
    merge_command::merge_command( const string& action )
        : command_i( action )
    {
    }
    
    merge_command::merge_command( const string& action, const string& name, const string& list )
        : m_action( action ), m_name( name ), m_list( list )
    {
    }

    merge_command::~merge_command( )
    {
    }
    
    bool merge_command::exec( )
    {
        vector<string> args;

        split( args, m_list, is_any_of( " {}" ), token_compress_on );
        
        if( args.size() == 0 )
        {
            throw std::runtime_error( "Error : can not understand list: " + m_list );
        }

        molecule_ptr pmol( new molecule_t() );

        morf_t prev( *pmol, RESD, -1);

        int head=0, tail=0;
        for( unsigned int i=1; i < args.size() - 1; i++ )
        {
            std::cout << "Sequence: " << args[i] << std::endl;
            molecule_ptr pseg = content().get_mol( args[i] );


            morf_t resd = merge(*pmol, *pseg, -1);
            resd.get_i(HEAD, head);

            if( m_action == "sequence" && i > 1)
            {
                 if( tail>0 && head>0 )
                 {
                    atom_t tail_atom( *pmol, tail-1);
                    atom_t head_atom( *pmol, head-1);
                    std::cout << "    tail: " << tail_atom.get_i(ID) << " " << tail_atom.get_s(NAME) << std::endl;
                    std::cout << "    head: " << head_atom.get_i(ID) << " " << head_atom.get_s(NAME) << std::endl;
 
                    place( prev, resd, tail_atom, head_atom );
                    bond_t bond = bond_t::create( tail_atom, head_atom );
	            bond.set_i(ORDER, 1);
                }
                else if( tail==0 )
                {
                    std::cout << "    tail atom not set." << std::endl;
                }
                else
                {
                    std::cout << "    head atom not set." << std::endl;
                }
            }

            int tailresd=0;
            if( pseg->nresd()==1 )
            {
                prev = resd;
            }
            else if( pseg->get_i("tailresd", tailresd) && tailresd !=0 )
            {
                prev = resd_t( *pmol, resd.absid()+tailresd-1 );
            }
            else
            {
                prev = resd_t( *pmol, resd.absid()+pseg->nresd()-1 );
            }

            tail = 0;
            prev.get_i( TAIL, tail );
        }
     
        pmol->set_s( NAME, m_name );
        content().set( m_name, pmol );
	return true;
    }
    

    void merge_command::undo( )
    {
        throw std::runtime_error( "sorry, not implemented yet." );
    }
    
    const char* merge_command::info( ) const
    {
        return " usage: variable = merge list ";
    }
    
    shared_ptr< command_i > merge_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 3 )
        {
            throw std::runtime_error( "Error: wrong number of arguments" );
        }
        
        return shared_ptr< command_i >( new merge_command( args[0], args[1], args[2] ) );
    }
    
} // namespace amber

amber::merge_command g_combine_command( "combine" );
amber::merge_command g_sequence_command( "sequence" );

