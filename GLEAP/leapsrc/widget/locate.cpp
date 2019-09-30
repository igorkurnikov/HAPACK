#include <common.hpp>
#include "parmfun.hpp"
#include "atomvec.hpp"
#include "atomfun.hpp"
#include "morange.hpp"
#include "molecule.hpp"
#include "database.hpp"

namespace mort
{

    bool inside( const numvec& rect, double x, double y, double z )
    {
        int dim = rect.size()/2;

        if( dim==1 )
        {
            double xmin = std::min( rect[0], rect[1] );
            double xmax = std::max( rect[0], rect[1] );
            return x > xmin && x < xmax;
        }
        

        if( dim==2 )
        {
            double xmin = std::min( rect[0], rect[2] );
            double xmax = std::max( rect[0], rect[2] );
            double ymin = std::min( rect[1], rect[3] );
            double ymax = std::max( rect[1], rect[3] );
            return x > xmin && x < xmax && y > ymin && y < ymax;
        }
        
        assert( dim==3 );
        
        double xmin = std::min( rect[0], rect[3] );
        double xmax = std::max( rect[0], rect[3] );
        double ymin = std::min( rect[1], rect[4] );
        double ymax = std::max( rect[1], rect[4] );
        double zmin = std::min( rect[2], rect[5] );
        double zmax = std::max( rect[2], rect[5] );
        return x > xmin && x < xmax && y > ymin && y < ymax && z > zmin && z < zmax;
    }
    
    bool locate_atom(const molecule_t& mol, const numvec& rect, atomvec_t& atms, int policy )
    {
        assert( atms.size()==0 );
        
        const vector<double> pos = get_vvec( mol, ATOM, POSITION );
        
        for( int i=0; i < mol.natom(); ++i )
        {
            double x = pos[3*i];
            double y = pos[3*i+1];
            double z = pos[3*i+2];

            if( inside(rect, x, y, z) )
            {
                atms.push_back( mol.atoms()[i] );
                
                if(policy==FIND_ONE)
                {
                    assert( atms.size()==1 );
                    return true;
                }
            }
        }

        return atms.size()>0;
    }

    bool locate_bond(const molecule_t& mol, const numvec& rect, bondvec_t& bnds, int policy )
    {
        assert( rect.size()==4 );

        const vector<double>& pos = get_vvec( mol, ATOM, POSITION );
        for( int i=0; i < mol.nbond(); ++i )
        {
            moref_t b = mol.bonds()[i];
            moref_t atm1 = atom_1st( b );
            moref_t atm2 = atom_2nd( b );
            
            int id1 = atm1.absid();
            int id2 = atm2.absid();

            double x = ( pos[3*id1  ] + pos[3*id2  ] ) * 0.5;
            double y = ( pos[3*id1+1] + pos[3*id2+1] ) * 0.5;
            double z = ( pos[3*id1+2] + pos[3*id2+2] ) * 0.5;
            
            if( inside(rect, x, y, z) )
            {
                bnds.push_back(b);
                if( policy==FIND_ONE )
                {
                    assert( bnds.size()==1 );
                    return true;
                }
            }    
        }

        return bnds.size()>0;
    }



                
    bool locate_atom( database_t& mdb, const numvec& rect, atomvec_t& atms, int policy )
    {
        database_t::iterator i = mdb.begin();
        for( ; i != mdb.end(); ++i )
        {
            molecule_ptr pm = dynamic_pointer_cast< molecule_t >( i->second );
            if( i->first[0] !='_' && pm != NULL )
            {
                locate_atom( *pm, rect, atms, policy );

                if( policy==FIND_ONE && atms.size()>0 )
                {
                    assert( atms.size()==1 );
                    return true;
                }        
            }
        }

        return atms.size()>0;
    }

    bool locate_bond( database_t& mdb, const numvec& rect, bondvec_t& bnds, int policy )
    {
        database_t::iterator i = mdb.begin();
        for( ; i != mdb.end(); ++i )
        {
            molecule_ptr pm = dynamic_pointer_cast< molecule_t >( i->second );
            if( i->first[0]!='_' && pm != NULL )
            {
                locate_bond( *pm, rect, bnds, policy );

                if( policy==FIND_ONE && bnds.size()>0 )
                {
                    assert( bnds.size()==1 );
                    return true;
                }        
            }
        }

        return bnds.size()>0;
    }

} // namespace mort

