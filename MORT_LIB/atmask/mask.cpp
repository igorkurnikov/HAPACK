#include <stdexcept>
#include <format.hpp>
#include <smarts.hpp>
#include "atom.hpp"
#include "resd.hpp"
#include "dist.hpp"
#include "logic.hpp"
#include "paren.hpp"

std::vector< const mort::atmask::node_i* >* insts()
{
    static std::vector<const mort::atmask::node_i*> holder;

    if( holder.size()==0 )
    {
         static const mort::atmask::atom_node g_atom_node;
	 static const mort::atmask::resd_node g_resd_node;
	 static const mort::atmask::dist_node g_dist_node;
	 static const mort::atmask::logic_node g_logic_node;
	 static const mort::atmask::paren_node g_paren_node;
	 holder.push_back(&g_atom_node);
	 holder.push_back(&g_resd_node);
	 holder.push_back(&g_dist_node);
	 holder.push_back(&g_logic_node);
	 holder.push_back(&g_paren_node);
    }

    return &holder;
}

namespace mort
{
    using std::logic_error;

    namespace atmask
    {

        node_i::node_i( bool reg )
        {
            // disabled if statement
            if( false && reg )
            {
                assert( insts != NULL );
                insts()->push_back( this );
            }
        }
        
        node_i::~node_i()
        {
        }

        const char* read_node(const char* ptr, shared_ptr<node_i>& curt )
        {
            assert(insts()->size() > 0); 

            for( int i=0; i < (int)insts()->size(); i++ )
            {
                if( insts()->at(i)->readable(ptr) )
                {
                    return insts()->at(i)->parse( ptr, curt );
                }
            }

            throw logic_error( "Error: can't understand mask " + string(ptr) );
        }

        const char* read_until(const char* ptr, char end, shared_ptr<node_i>& curt )
        {
            while( *ptr != end )
            {
                ptr = read_node( ptr, curt );
            }
            
            return ptr;
        }

        shared_ptr< node_i > read_mask( const string& mask )
        {
            shared_ptr< node_i > curt;
            
            read_until( mask.c_str(), '\0', curt );

            return curt;
        }

    } // namespace atmask

    atomvec_t mask_atom( const molecule_t& mol, const string& mask )
    {        
        shared_ptr< atmask::node_i > root = atmask::read_mask( mask );
        return ( root != NULL ) ? root->match(mol) : atomvec_t();
    }

    void mask_atom( const molecule_t& mol, const string& mask, vector<int>& maskvec )
    {
        maskvec.resize( mol.natom() );
	
	std::fill( maskvec.begin(), maskvec.end(), 0 );

        atomvec_t atoms = mask_atom( mol, mask );

	for( int i=0; i < (int)atoms.size(); ++i )
        {
            maskvec.at( atoms[i].absid() ) = 1;
        }
    }

    atomvec_t smarts_mask_atom(const molecule_t& m, const string& expr)
    {
        molecule_t sub;
	read_smarts(expr.c_str(), sub);

        atomvec_t result;
	atomiter_t ai = m.atom_begin();
	for( ; ai != m.atom_end(); ++ai )
	{
	    if( has_subst(*ai, sub) )
	    {
	        result.push_back(*ai);
	    }
        }

	return result;
    }   

    void smarts_mask_atom(const molecule_t& m, const string& expr, vector<int>& maskvec)
    {
        molecule_t sub;
	read_smarts(expr.c_str(), sub);

	maskvec.clear();
	atomiter_t ai = m.atom_begin();
	for( ; ai != m.atom_end(); ++ai)
	{
	    int masked = has_subst(*ai, sub) ? 1 : 0;
	    maskvec.push_back(masked);
	}
    }

} // namespace mort

