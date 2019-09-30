#include <object.hpp>

#include "les.hpp"

namespace mort
{
    
    namespace prmtop
    {

        void atom2angl( const molecule_t* pff, const atom_t& ac, bool fixed_water_model )
        {
            if( ac.type()=="HW" )
            {
                return;
            }

            if( fixed_water_model && ac.get_s(TYPE)=="OW" )
            {
                return;
            }

//			std::cout << "atom2angl() " << "atoms bonded to " << ac.get_s(RESNAME) << ":" << ac.get_s(NAME) << std::endl; 
            
            atomiter_t ai = ac.atom_begin();
            atomiter_t ae = ac.atom_end();

//			for( ; ai != ae; ++ai )
//            {
//				std::cout << "          " << (*ai).get_s(RESNAME) << ":" << (*ai).get_s(NAME) << std::endl;
//			}
//			ai = ac.atom_begin();

            for( ; ai != ae; ++ai )
            {
                atomiter_t aj = ai+1;
                for( ; aj != ae; ++aj )
                {
                    if( les_forbid(*ai, *aj) )
                    {
                        continue;
                    }

//					std::cout << "atom2angle() : " << (*ai).get_s(RESNAME) << ":" << (*ai).get_s(NAME) << "-";
//					std::cout << ac.get_s(RESNAME) << ":" << ac.get_s(NAME) << "-";
//					std::cout << (*aj).get_s(RESNAME) << ":" << (*aj).get_s(NAME) << std::endl;

                    atmvec as(*ai, ac, *aj);
                    parm_angl(pff, as);
                }
            }
        }

        // bc: central bond
        void bond2dihe( const molecule_t* pff, const bond_t& bc )
        {
            atom_t a1 = bc.atoms()[0];
            atom_t a2 = bc.atoms()[1];
            assert( a1 != a2 );

            atomiter_t ai = a1.atom_begin();
            for( ; ai != a1.atom_end(); ++ai )
            {
                if( *ai==a2 )
                {
                    continue;
                }

                if( les_forbid(*ai, a2) )
                {
                    continue;
                }       

                atomiter_t aj = a2.atom_begin();
                for( ; aj != a2.atom_end(); ++aj )
                {
                    if( les_forbid(*ai, *aj) || les_forbid(a1, *aj) )
                    {
                        continue;
                    }

                    if( (*aj!=a1) && (*aj!=*ai) )
                    {
                        atmvec as(*ai, a1, a2, *aj);
                        parm_dihe( pff, as );
                    }
                }
            }
        }
        
        void atom2impr( const molecule_t* pff, const atom_t& ac )
        {
            atomiter_t ab = ac.atom_begin();
            atomiter_t ae = ac.atom_end();

            int max_score=0;
            
            for(atomiter_t ai=ab; ai!=ae; ++ai)
            {
                for(atomiter_t aj=ab; aj!=ae; ++aj)
                {
                    if( aj==ai )
                    {
                        continue;
                    }

		    if( les_forbid( *aj, *ai ) )
                    {
                        continue;
                    }                 

                    for( atomiter_t ak=ab; ak!=ae; ++ak )
                    {
                        if( ak==ai || ak==aj )
                        {
                            continue;
                        }
                     
                        if( les_forbid(*ak,*ai) || les_forbid(*ak,*aj) )
                        {
                            continue;
                        }
  
                        atmvec as(*ai, *aj, ac, *ak);
                        if( pff->noops() > 0 )
                        {
                            parm_sander_impr(pff, as, max_score);
                        }
                        else
                        {
                            parm_amoeba_impr(as);
                        }
                    }
                }
            }
        }

        void angl2tor2( const molecule_t* pff, const morf_t& angl )
        {
            morf_t atm1 = atom_1st(angl);
            morf_t atm2 = atom_2nd(angl);
            morf_t atm3 = atom_3rd(angl);
            
            assert( atm1!=atm2 && atm2!=atm3 );

            atomiter_t atm0 = atm1.atom_begin();
            for( ; atm0 != atm1.atom_end(); ++atm0 )
            {
                if( *atm0==atm2 || *atm0==atm3 )
                {
                    continue;
                }

                atomiter_t atm4 = atm3.atom_begin();
                for( ; atm4 != atm3.atom_end(); ++atm4 )
                {
                    if( *atm4==atm1 || *atm4==atm2 || *atm4==*atm0 )
                    {
                        continue;
                    }

                    atomvec_t atoms(*atm0, atm1, atm2, atm3, *atm4);
                    parm_tor2(pff, atoms);
                }
            }
        }
        
        void bond2ptor( const morf_t& bond )
        {
            assert( bond.has_d(PITORS) );

            if( bond.get_d(PITORS) == 0.0 )
            {
                return;
            }

            morf_t atm1 = atom_1st( bond );
            morf_t atm2 = atom_2nd( bond );

            if( atm1.natom() != 3 )
            {
                return;
            }
            
            if( atm2.natom() != 3 )
            {
                return;
            }

            atomvec_t nbr1;
            atomiter_t atm0 = atm1.atom_begin();
            for( ; atm0 != atm1.atom_end(); ++atm0 )
            {
                if( *atm0 != atm2 )
                {
                    nbr1.push_back( *atm0 );
                }
            }

            assert( nbr1.size() == 2 );
            
            atomvec_t nbr2;
            atomiter_t atm3 = atm2.atom_begin();
            for( ; atm3 != atm2.atom_end(); ++atm3 )
            {
                if( *atm3 != atm1 )
                {
                    nbr2.push_back( *atm3 );
                }
            }
            
            assert( nbr2.size() == 2 );

            atomvec_t atoms( nbr1[0], nbr1[1], atm1, atm2, nbr2[0], nbr2[1] );

            morf_t ptor = create_ptor( atoms, 2 );           
            ptor.set_d(FORCE, bond.get_d(PITORS) );
            ptor.set_d(EQUIL, M_PI);
        }
        
    } // namespace prmtop
    
} // namespace mort

