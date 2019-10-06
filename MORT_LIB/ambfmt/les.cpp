

#include <map>
#include <iomanip>

#include <object.hpp>
#include <common/fortran.hpp>
#include "les.hpp"
#include "parmset.hpp"

namespace mort
{
    using std::map;
    using std::for_each;

    namespace prmtop
    {
        void set_les(morf_t& obj)
        {
            int lestype;
            if( ! obj.getmol().get_i(LESTYPE, lestype) || lestype == 0 )
            {
                return;
            }

            obj.set_i(LESTYPE, 0);

            atomiter_t atom = obj.atom_begin();
            for( ; atom != obj.atom_end(); ++atom )
            {
                if( atom->get_i(LESTYPE, lestype) && lestype > 0 )
                {
                    obj.set_i( LESTYPE, lestype);
                }
            }
        }

        bool les_forbid( const morf_t& atm1, const morf_t& atm2 )
        {
            
            if( ! atm1.has_i(LESTYPE) )
            {
                return false;
            }
            
            if( atm1.get_i(LESTYPE) != atm2.get_i(LESTYPE) )
            {
                return false;
            }
            
            return atm1.get_i(LESCOPY) != atm2.get_i(LESCOPY);
        }


        void les_exclude( const morf_t& atom, vector<int>& list, vector<int>& dist )
        {
            // for les atom, put atoms from other copy into its exclusion list 
	    int lestype;
            if( atom.get_i(LESTYPE,lestype) && lestype > 0 )
            {
                atomiter_t excl = atom.getmol().atom_begin();
                for( ; excl != atom.getmol().atom_end(); ++excl )
                {
                    if( ! les_forbid( *excl, atom ) )
                    {
                        continue;
                    }

                    if( excl->get_i(ID) <= atom.get_i(ID) )
                    {
                        continue;
                    }

                    if( std::find( list.begin(), list.end(), excl->get_i(ID) ) == list.end() )
                    {
                        list.push_back( excl->get_i(ID) );
                        dist.push_back( 9 );
                    }
                }
            }    
        }

        struct les_scale_each
        {
             les_scale_each(int ncopy)
             {
                 m_ncopy = ncopy;
             }

             void operator()(morf_t& obj) const
             {
	          int lestype;
		  double force;
                  assert( obj.get_i(LESTYPE,lestype) && obj.get_d(FORCE,force) );

                  if( lestype > 0 )
                  {
                      obj.set_d(FORCE, force/m_ncopy);
                  }
             }

             int m_ncopy;
        };

        void les_scale( molecule_t& mol, parmset_t& )
        {
	    int mollestype;
            assert( mol.get_i(LESTYPE, mollestype) && mollestype == 1 );
            
            int ncopy = mol.get_i(LESCOPY);

            for_each( mol.bond_begin(), mol.bond_end(), les_scale_each( ncopy ) );

            for_each( mol.angl_begin(), mol.angl_end(), les_scale_each( ncopy ) );
            
            for_each( mol.dihe_begin(), mol.dihe_end(), les_scale_each( ncopy ) );
            
            for_each( mol.impr_begin(), mol.impr_end(), les_scale_each( ncopy ) );
            
	    atomiter_t atom = mol.atom_begin();
	    for( ; atom != mol.atom_end(); ++atom )
	    {
	        if( atom->get_i(LESTYPE) > 0 )
		{
                    atom->set_d(DEPTH, atom->get_d(DEPTH)/(ncopy*ncopy) );
		    atom->set_d(PCHG,  atom->get_d(PCHG) / ncopy );
                }
            }
        }

        void write_les( ostream& os,  const molecule_t& mol )
        {
	    os << "%FLAG LES_NTYP" << std::endl;
	    os << "%FORMAT(10I8)" << std::endl;
            os << std::setw( 8 ) << 2 << std::endl;

            os << "%FLAG LES_TYPE" << std::endl;
            os << "%FORMAT(10I8)"  << std::endl;
            fortran_t i_format( os, "10I8" );
	    atomiter_t atom = mol.atom_begin();
            for( ; atom != mol.atom_end(); ++atom )
            {
	        i_format( atom->get_i(LESTYPE) + 1 );
            } 
            i_format.end();

            os << "%FLAG LES_FAC" << std::endl;
            os << "%FORMAT(5E16.8)" << std::endl;
            fortran_t d_format( os, "5E16.8" );
            d_format( 1.0 );
            d_format( 1.0 );
            d_format( 1.0 );
            d_format( mol.get_i(LESCOPY) );
            d_format.end();

            os << "%FLAG LES_CNUM" << std::endl;
	    os << "%FORMAT(10I8)" << std::endl;
	    i_format.begin();
            atomiter_t ai = mol.atom_begin();
            atomiter_t ae = mol.atom_end();
            for( ; ai != ae; ++ai )
            {
                fortran_write_iparm( &i_format, *ai, LESCOPY );
            }
	    i_format.end();

            os << "%FLAG LES_ID" << std::endl;
	    os << "%FORMAT(10I8)" << std::endl;
	    i_format.begin();
            ai = mol.atom_begin();
            ae = mol.atom_end();
            for( ; ai != ae; ++ai )
            {
                fortran_write_iparm( &i_format, *ai, LESTYPE );
            }
	    i_format.end();
        }

    } // namespace prmtop

} // namespace mort


