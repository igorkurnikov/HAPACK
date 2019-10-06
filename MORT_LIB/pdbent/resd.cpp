#include <common.hpp>
#include <object.hpp>

#include "atom.hpp"

namespace mort
{

    namespace pdbent
    {
        double bondcut()
        {
            string d;
            return mortenv().get_s("pdbbondcut",d) ? atof(d.c_str()):1.875;
        }



        struct resd_s
        {
            resd_s(const string& line)
            {
//				printf("resd_s::resd_s()  line = %s \n",line.c_str());   // IGOR_TMP
                if( line.length() >= 10) name = line.substr(RESD_TYPE, 10);
                if( line.length() >= 3) type = line.substr(RESD_TYPE, 3);
                if( line.length() >= RESD_CHAIN) chain = line[RESD_CHAIN];
                rseq = atoi( line.c_str() + RESD_SEQ );
                if( line.length() >= RESD_ICODE) icode = line[RESD_ICODE];
                name.erase( std::remove( name.begin(), name.end(), ' ' ), name.end() );
                type.erase( std::remove( type.begin(), type.end(), ' ' ), type.end() );
            }

            bool operator==( const resd_s& rhs )
            {
                return name == rhs.name;
            }

            bool operator!=( const resd_s& rhs )
            {
                return name != rhs.name;
            }

            string name;
            string type;
            int    rseq;
            char   chain;
            char   icode;
        };


        resd_s peek_resd( istream& is )
        {
            //assert( is != NULL );
            
            char line[MAX_LINE_WIDTH];

            std::ios::pos_type pos = is.tellg();
        
            is.getline(line, MAX_LINE_WIDTH);
            
            is.seekg(pos);
        
            return resd_s( line );
        }

        resd_t read_resd( istream& is, molecule_t& m )
        {    
            resd_t resd = resd_t::create( m );

            string keyw = peek_keyw( is );
            resd_s rinf = peek_resd( is );

            resd.set_s(NAME, rinf.name );
            resd.set_s(TYPE, rinf.type );
            resd.set_i(CHAIN, rinf.chain );
            
            while( (keyw=="ATOM"||keyw == "HETA") && peek_resd(is)==rinf )
            {
                read_atom( is, resd );
                keyw = peek_keyw( is );
            }

            //bond_bydis( resd, bondcut() );
            return resd;
        }

        // make bond between two neighbouring residues
        bool bond_resd( resd_t& ri, resd_t& rj )
        {
            assert( ri.cmpid()==RESD && rj.cmpid()==RESD );

            double cut = bondcut();

            atomiter_t ai = ri.atom_begin();
            atomiter_t ei = ri.atom_end();
        
            atmvec conns;
            for( ; ai != ei; ++ai )
            {
                atomiter_t aj = rj.atom_begin();
                atomiter_t ej = rj.atom_end();
                if( ai->get_i(ELEMENT)==HYDROGEN ) 
                    continue;            

                for( ; aj != ej; ++aj )
                {
                    if( aj->get_i(ELEMENT)==HYDROGEN )
                        continue;
                    
                    if( bond_bydis(*ai, *aj, cut) )
                    {
                        if( conns.size()==0 )
                        {
                            ri.set_i( TAIL, ai->absid()+1 );
                            rj.set_i( HEAD, aj->absid()+1 );
                        }
                        else
                        {
                            ri.set_i( CONN1+conns.size()/2-1, ai->absid()+1 );
                            rj.set_i( CONN1+conns.size()/2-1, aj->absid()+1 );
                        }

                        conns.push_back( *ai );
                        conns.push_back( *aj );
        
                    }
                }
            }

            if( conns.size()>2 )
            {
                std::cout << "Warning: there is more than one bond made between residue ";
                std::cout << ri.name() << " and " << rj.name() << ":" << std::endl;
                for( unsigned int i=0; i < conns.size()/2; ++i )
                {
                    std::cout << "Warning:        bond " << i+1 << ": " << conns[2*i].name() << "-" << conns[2*i+1].name() << std::endl;
                }
                std::cout << "Warning: Maybe current pdbbondcut (" << cut << ") is too large?" << std::endl;
            }

            return conns.size() > 0;
        }
 

        void write_resd( ostream& os, const resd_t& resd )
        {
            funstack_t::push( "write_resd" );

            atomiter_t ai = resd.atom_begin();
            atomiter_t ae = resd.atom_end();
            for( ; ai != ae; ++ai )
            {
                write_atom( os, *ai );
            }

            funstack_t::pop();
        }        
 
        void setbgn( resd_t& r )
        {
            r.set_i( AAPOS, NTERM );
            r.set_i( HEAD,  0 );
        }

        void setend( resd_t& r )
        {
            int pos = r.get_i(AAPOS);
            if( pos==NTERM )
            {
                r.set_i( AAPOS, ALONE );
            }
            else
            {
                r.set_i( AAPOS, CTERM );
            }

            r.set_i( TAIL, 0 );
        }

        void setpos( resd_t& rprev, resd_t& rcurt )
        {
            if( rprev.absid()==-1 )
            {
                setbgn( rcurt );
                return;
            }

            int prevpos = rprev.get_i(AAPOS);
            if( prevpos == CTERM || prevpos == ALONE )
            {
                setbgn( rcurt );
                return;
            }

            if( rprev.chain()!=rcurt.chain() )
            {
                setend( rprev );
                setbgn( rcurt );
                return;
            }

            rcurt.set_i( AAPOS, OTHER );
        }
 
    } // namespace pdbent

} // namespace mort

