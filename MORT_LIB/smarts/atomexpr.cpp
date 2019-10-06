#include <common.hpp>
#include "ring.hpp"
#include "atomexpr.hpp"

namespace mort
{
    namespace atomexpr
    {
        
        node_i::node_i( int level )
        {
            m_level = level;
            m_father= NULL;
        }

	node_i::~node_i()
	{
	}
        
        int node_i::level( )
        {
            return m_level;
        }
        
        node_i* node_i::father( )
        {
            return m_father;
        }
        
        void node_i::set_father( node_i* father )
        {
            m_father = father;
        }

        // parm_node implementation
	iparm_node::iparm_node(const hashid_t& parmid, int ivalue )
	    : node_i(SINGLE)
	{
	    m_parmid = parmid;
	    m_ivalue = ivalue;
	}

	bool iparm_node::match( const morf_t& a ) const
	{
	    return a.get_i(m_parmid) == m_ivalue;
	}


	// env_node implementation
        
        env_node::env_node(const hashid_t& env, int value )
            :node_i( SINGLE )
        {
            m_env = env;
            m_value = value;
        }

        env_node::~env_node()
        {
        }
        
        bool env_node::match( const morf_t& atom ) const 
        {
            int result = -1;
            
            switch( m_env )
            {
            case HYBRID:
                result = get_atom_hybrid( atom );
                break;
                
            case DEGREE:
                result = get_atom_degree( atom );
                break;
                
            case VALENCE:
                result = get_atom_valence( atom );
                break;
                
            case HYDROGEN:
                result = get_atom_hydrogen( atom );
                break;
                
            case NEBR:
                result = atom.natom();
                break;
                
            default:
                assert( false );                
            }

            // std::cout << "matching env " << m_env << " " << m_value << " " << result << "\n";

            return result == m_value;            
        }

        logic_node::logic_node( int level, const hashid_t& logic )
            :node_i( level )
        {
            m_logic = logic;
        }

        logic_node::~logic_node()
        {
        }

        void logic_node::push( const shared_ptr< node_i >& node )
        {
            assert( m_logic != NOT || m_children.empty() );
    
            m_children.push_back( node );

            node->set_father( this );
        }
        
        void logic_node::pop()
        {
            assert( !m_children.empty() );
            m_children.pop_back();
        }
         
        shared_ptr< node_i > logic_node::top() const
        {
            return m_children.back();
        }

        int logic_node::logic() const
        {
            return m_logic;
        }
        
        bool logic_node::match( const morf_t& atom ) const
        {
            if( m_logic == NOT )
            {
                return !m_children[0]->match( atom );
            }
            
            assert( m_logic == AND || m_logic == OR );
            
            for( int i=0; i < (int)m_children.size(); i++ )
            {
                bool result = m_children[ i ]->match( atom );
                
                if( ( m_logic == AND ) && !result )
                {
                    return false;
                }
                
                if( ( m_logic == OR ) && result )
                {
                    return true;
                }
            }

            return m_logic == AND;
        }

        ring_node::ring_node( int arom, int size )
            :node_i( SINGLE )
        {
            m_arom = arom;
            m_size = size;
        }

        ring_node::~ring_node()
        {
        }

        bool ring_node::match( const morf_t& atom ) const
        {
            int ringed;

            if( ! atom.get_i(RING, ringed) )
            {
                ringed = has_ring(atom);
            }

            if( !ringed )
            {
                return false;
            }

            if( m_size!=0 && !has_ring(atom, m_size) )
            {
                return false;
            }
            
            int arom = is_arom( atom );

            return arom == m_arom;
        }

    } // namespace atomexpr

    using namespace atomexpr;

    atomexpr_t::atomexpr_t()
    {
        m_curt = NULL;
    }
        
    atomexpr_t::atomexpr_t( const atomexpr_t& )
    {
    }
        
    atomexpr_t::~atomexpr_t()
    {
    }
        
    void atomexpr_t::change_to( int level, const hashid_t& logic )
    {
        while( m_curt && m_curt->level() < level )
        {
            m_curt = m_curt->father();
        }

        if( m_curt == NULL )
        {
            shared_ptr< logic_node > node( new logic_node( level, logic ) );

            if( m_root )
            {                
                node->push( m_root );
            }

            m_root = node;
            m_curt = node.get();
        }
        else if( m_curt->level() > level )
        {
            shared_ptr< logic_node > node( new logic_node( level, logic ) );

            logic_node* curt = dynamic_cast< logic_node* >( m_curt );
            
            assert( curt != NULL );

            if( logic != NOT )
            {
                node->push( curt->top() ); 
                curt->pop();
            }

            curt->push( node );
            m_curt = node.get();
        }

    }        

    void atomexpr_t::insert( const shared_ptr< node_i >& node )
    {
        assert( node->level() == SINGLE );
            
        if( m_curt == NULL )
        {
            assert( m_root == NULL );
            m_root = node;
        }
        else 
        {
            if( m_curt->level() == SINGLE )
            {
                change_to( atomexpr::CONN, AND );
            }

            logic_node* curt = dynamic_cast< logic_node* >( m_curt );

            assert( curt != NULL );

            curt->push( node );
        }

        m_curt = node.get();

    }

    bool atomexpr_t::match( const morf_t& atom ) const
    {
        return m_root->match( atom );
    }
        
} // namespace mort

 
