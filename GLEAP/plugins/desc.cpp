#include <sstream>
#include <object.hpp>
#include "desc.hpp"
#include <sstream>

namespace amber
{
    using namespace std;
    
    string atom_info( const morf_t& atom )
    {
        ostringstream os;
        
        os << ".R< " << atom.get_s(RESNAME);
        os << " " << atom.get_i(RESID) << ">";
        os << ".A< " << atom.get_s(NAME);

        if( atom.getmol().natom()==0 ) 
        {
            os << " " << atom.get_i(ID) << ">";
        }
        else
        {
            int relate = atom.absid() - atom.resd_begin()->atom_begin()->absid() + 1;
            os << " " << relate << ">";
        }
            
        return os.str();
    }

    void desc_unit( const molecule_t& mol )
    {
        leaplog_t::putline( "UNIT name: " + mol.get_s(NAME) );
        
	int head;
        if( mol.get_i(HEAD,head) && head != 0 )
        {
            morf_t head_atom(mol, ATOM, head-1);
            leaplog_t::putline( "Head atom: " + atom_info(head_atom) );
        }
        
	int tail;
        if( mol.get_i(TAIL,tail) && tail != 0 )
        {
            morf_t tail_atom(mol, ATOM, tail-1);            
            leaplog_t::putline( "Tail atom: " + atom_info(tail_atom) );
        }
        
        leaplog_t::putline( "Contents:" );
        
        if( mol.nresd() == 0 )
        {
            leaplog_t::putline( "R<" + mol.get_s(NAME) + " 1>" );
        }
        else
        {
            resditer_t resd = mol.resd_begin();
            for( ; resd != mol.resd_end(); ++resd )
            {
                ostringstream os;
                os << "R<";
                os << resd->get_s(NAME) << " ";
                os << resd->get_i(ID) << ">";
                leaplog_t::putline( os.str() );
            }
        }
    }

    void desc_resd( const morf_t& resd )
    {
        leaplog_t::putline( "RESIDUE name: " + resd.get_s(NAME) );
        
	int head;
        if( resd.get_i(HEAD,head) && head != 0 )
        {
            morf_t head_atom(resd.getmol(), ATOM, head-1);
            leaplog_t::putline( "Connect atom 0: " + atom_info(head_atom) );            
        }
        
	int tail;
        if( resd.get_i(TAIL,tail) && tail != 0 )
        {
            morf_t tail_atom(resd.getmol(), ATOM, tail-1);
            leaplog_t::putline( "Connect atom 1: " + atom_info(tail_atom) );
        }

        int conn;
        if( resd.get_i(CONN1,conn) && conn != -1 )
        {
            morf_t conn_atom(resd.getmol(), ATOM, conn);
            leaplog_t::putline( "Connect atom 2: " + atom_info(conn_atom) );
        }

        leaplog_t::putline( "Contents:" );
        
        atomiter_t atom = resd.atom_begin();
        for( ; atom != resd.atom_end(); ++atom )
        {
            ostringstream os;
            os << "A<";
            os << atom->get_s(NAME) << " ";
            os << atom->get_i(ID) << ">";
            leaplog_t::putline( os.str() );
        }
    }

    void desc_atom( const morf_t& atom )
    {
        leaplog_t::putline( "ATOM" );

        leaplog_t::putline( "Name: " + atom.get_s(NAME) );

        leaplog_t::putline( "Type: " + atom.get_s(TYPE) );
        
        ostringstream os_chrg;
        os_chrg << "Charge: " << atom.get_d(PCHG);
        leaplog_t::putline( os_chrg.str() );
        
        leaplog_t::putline( "Element: " + atom.get_s(SYMBOL) );
        
        ostringstream os_pose;
	numvec pos = atom.get_v(POSITION);
        os_pose << "Atom position: " << pos[0] << " " << pos[1] << " " << pos[2];
        leaplog_t::putline( os_pose.str() );
        
        atomiter_t nbr = atom.atom_begin();
        for( ; nbr != atom.atom_end(); ++nbr )
        {
            leaplog_t::putline( "Bonded to " + atom_info( *nbr ) + " by a single bond " );
        }
        
    }
    
    desc_command::desc_command( )
        : command_i( "desc" )
    {
    }
    
    desc_command::desc_command( const string& object )
        : m_object( object )
    {
    }

    desc_command::~desc_command( )
    {
    }
    
    bool desc_command::exec( )
    {
        int ndot = count( m_object.begin(), m_object.end(), '.' );
        
        if( ndot > 2 )
        {
            throw logic_error( "Error: can't understand " + m_object );
        }

        if( ndot == 0 )
        {
            molecule_ptr pmol = content().get_mol( m_object );
            assert( pmol != NULL );
            desc_unit( *pmol );
        }
        else if( ndot == 1 )
        {
            desc_resd( content().get_resd(m_object) );
        }
        else 
        {
            assert( ndot == 2 );
            desc_atom( content().get_atom(m_object) );
        }
        return true;
    }
    
    void desc_command::undo( )
    {
    }
    
    const char* desc_command::info( ) const
    {
        return "  usage: desc object";
    }
    
    shared_ptr< command_i > desc_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 2 )
        {
            throw logic_error( "Error: wrong number of arguments" );
        }
        
        return shared_ptr< command_i >( new desc_command( args[1] ) );
    }
    
};

amber::desc_command g_desc_command;
