#include <fstream>
#include <sstream>
#include <common.hpp>
#include <object.hpp>
#include "atomfun.hpp"
#include "parmfun.hpp"
#include "bondfun.hpp"

namespace mort
{

    ostream& alog()
    {
        static std::ofstream inst( "addhs.log" );
        return inst;
    }
    

    bool addable( int elem )
    {
        static const int ELEMS[]={CARBON, NITROGEN, OXYGEN, PHOSPHORUS, SULFUR};
        int id = std::find( ELEMS, ELEMS+5, elem )-ELEMS;
        return id < 5;
    }

    int maxbo( atom_t& a )
    {
        int r = 0;

        bonditer_t bi = a.bond_begin();
        for( ; bi != a.bond_end(); ++bi )
        {
            r = std::max( r, bi->order() );
        }
        
        return r;
    }

    int valence( int elem )
    {
        static const int ELEMS[]={CARBON, NITROGEN, OXYGEN, PHOSPHORUS, SULFUR};
        static const int VALES[]={     4,        3,      2,          3,      2};
        
        int id = std::find( ELEMS, ELEMS+5, elem )-ELEMS;
        if( id == 5 )
        {
            throw std::runtime_error( "Error: do not know valence of this element" );
        }
   
        return VALES[id];
    }
        

    void addH( atom_t& a, const numvec& po )
    {
        atom_t h = atom_t::create( a.getmol(), "" );
        bond_t b = bond_t::create( a, h );        
        h.set_v( POSITION, po );
        h.set_i( ELEMENT,  HYDROGEN );
        h.set_s( TYPE,     "H" );        
        h.set_s( NAME,     uniq_name(h) );
        b.set_i( ORDER,    1 );

	/*
        if( a.natom() > 0 )
        {
            resd_t r = a.resd();
            h.connect( r );
            r.connect( h );
        }
	*/
    }

    bool too_close( const vector<numvec>& crds, const numvec& po, double mindis )
    {
        double mdis2 = mindis * mindis;
        
        for( int i=0; i < (int)crds.size(); ++i )
        {
            if( dis2(crds[i], po) < mdis2 )
                return true;
        }
        
        return false;
    }
    
    void addHs_impl( atom_t& a, int addval, double intang, double rotang, vector<numvec>& hpos )
    {
        assert( hpos.size()==0 );
        
        static const double BOND_LENGTH = 1.3;
        int elem = a.element();
        int nhyb = valence(elem) - a.natom() - addval;
        if( nhyb<=0 ) return;

        vector< numvec > crds( 1, a.get_v(POSITION) );
        transform( a.atom_begin(), a.atom_end(), back_inserter(crds), vparm_getter(POSITION) );        
        assert( crds.size() > 0 );
        
        if( crds.size()==1 )
        {
            alog() << "        no nbr, add one to x axis" << std::endl;
            numvec po = crds[0] + makevec(BOND_LENGTH, 0.0, 0.0);
            hpos.push_back( po );
            crds.push_back( po );
            nhyb--;
            if( nhyb==0 ) return;
        }
        
        assert( nhyb>0 && crds.size() >= 2 );
        numvec rot = crds[1] - crds[0];
        normalize(rot);
        rot *= BOND_LENGTH;

        if( crds.size()==2 )
        {
            alog() << "        one nbr, add another one by rotate " << intang << std::endl;

            numvec ax = makevec( -rot[1], rot[0], 0.0);
            normalize( ax );

            numvec po = crds[0] + rotcpy( rot, ax, intang );
            hpos.push_back( po );
            crds.push_back( po );
            nhyb--;
            if( nhyb==0 ) return;
        }


        alog() << "        more than one nbr, add rest by rotate " << rotang << std::endl;        
        assert( nhyb>0 && crds.size() >= 3 );
        numvec axs = crds[2]-crds[0];
        normalize(axs);

        assert( nhyb > 0 );

        int i=1;
        while( nhyb > 0  && i < 4)
        {
            numvec po = crds[0] + rotcpy( rot, axs, rotang*i );
            i++;
            
            if( !too_close(crds, po, 1.0) )
            {
                crds.push_back( po );
                hpos.push_back( po );
                nhyb--;
            }
        }

        assert( i<4);
    }

    void addHs_sp3( atom_t& a, vector<numvec>& hpos )
    {
        addHs_impl( a, 0, 109.5, 120.0, hpos );
    }
    
    void addHs_sp2( atom_t& a, vector<numvec>& hpos )
    {
        addHs_impl( a, 1, 120.0, 180.0, hpos );
    }
    
    void addHs_sp1( atom_t& a, vector<numvec>& hpos )
    {
        addHs_impl( a, 2, 180.0, 180.0, hpos );
    }    

    void addHs( atom_t& a, vector<numvec>& hpos )
    {
        if( !addable(a.get_i(ELEMENT)) )
            return;

        alog() << "add hydrogen to atom: " << a.absid() << std::endl;
   
        int mbo = maxbo( a );
        
        if( mbo==0 || mbo==1 )
        {
            alog() << "    type sp3 " << std::endl;
            addHs_sp3( a, hpos );
            alog() << "    finished " << std::endl;
        }
        else if( mbo==2 )
        {
            alog() << "    type sp2 " << std::endl;
            addHs_sp2( a, hpos );
            alog() << "    finished " << std::endl;
        }
        else
        {
            alog() << "    type sp1 " << std::endl;
            assert( mbo==3 );
            addHs_sp1( a, hpos );
            alog() << "    finished " << std::endl;
        }
    }


    void addHs( atom_t& a )
    {
        assert( a.cmpid()==ATOM );
        
        vector<numvec> hpos;
        addHs( a, hpos );
        
        for( int i=0; i < (int)hpos.size(); ++i )
        {
            addH( a, hpos[i] );
        }
    }
    
    void addHs( molecule_t& m )
    {
        vector< vector<numvec> > hpos( m.natom() );
        for( int i=0; i < m.natom(); ++i )
        {
            atom_t a = m.atoms()[i];
            addHs( a, hpos[i] );
        }
        
        for( unsigned int i=0; i < hpos.size(); ++i )
        {
            for( unsigned int j=0; j < hpos[i].size(); ++j )
            {
                atom_t a = m.atoms()[i];
                addH( a, hpos[i][j] );
            }
        }
        
    }
    
    void addHs( atomvec_t& avec )
    {
        vector< vector<numvec> > hpos( avec.size() );
        for( unsigned int i=0; i < avec.size(); ++i )
        {
            atom_t a = avec[i];
            addHs( a, hpos[i] );
        }
        
        for( unsigned int i=0; i < avec.size(); ++i )
        {
            for( unsigned int j=0; j < hpos[i].size(); ++j )
            {
                atom_t a = avec[i];
                addH( a, hpos[i][j] );
            }
        }

    }

    void addHs( database_t& mdb )
    {
        database_t::iterator i = mdb.begin();
        for( ; i != mdb.end(); ++i )
        {
            molecule_ptr pm = dynamic_pointer_cast< molecule_t >( i->second );
            if( pm != NULL )
            {
                addHs( *pm );
            }
        }
    }
 

    void delHs( molecule_t& m )
    {
        atomvec_t hs;

	atomiter_t ai = m.atom_begin();
	for( ; ai != m.atom_end(); ++ai )
	{
	    if( ai->get_i(ELEMENT) == HYDROGEN )
	    {
	        hs.push_back( *ai );
	    }
	}

	for( unsigned int i=0; i < hs.size(); ++i )
	{
	    m.remove_atom( hs[i] );
	}
    }

  
} // namespace mort
