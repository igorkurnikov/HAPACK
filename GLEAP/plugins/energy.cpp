#include <sstream>
#include <object.hpp>
#include <ambfmt.hpp>
#include <atmask.hpp>
#include <enefrc.hpp>

#include "energy.hpp"


namespace amber
{
    energy_command::energy_command( )
        : command_i( "energy" )
    {
    }
    
    energy_command::energy_command(const string& name, const string& parm)
        : m_name(name), m_parm(parm)
    {
    }
    
    energy_command::energy_command(const string& name, const string& parm, const string& vec1, const string& vec2)
        : m_name(name), m_parm(parm), m_vec1(vec1), m_vec2(vec2)
    {
    }

    energy_command::~energy_command( )
    {
    }
    
    bool energy_command::exec( )
    {
        //leaplog_t::putline( "Energy of " + m_name );
        
        //std::ostringstream os;
        ostream& os(std::cout);

        parmset_t params;

        molecule_ptr pmol = content().get_mol(m_name);
        molecule_ptr parm = content().get_mol("_amber-atom");

        if( pmol== NULL )
        {
            throw std::runtime_error( "Error: molecule " + m_name + " does not exist" );
        }

        if( parm== NULL )
        {
            throw std::runtime_error( "Error: amber force field parm not found" );
        }

        assert( pmol != NULL && parm != NULL );
        parametrize( *pmol, *parm, params );
       
        if( m_vec1.empty() )
	{
            os << "Bond:   " << format( "%20.4f" ) % eval_bond( *pmol ) << " kcal/mol" << std::endl;
            os << "Angle:  " << format( "%20.4f" ) % eval_angl( *pmol ) << " kcal/mol" << std::endl;
            os << "Torsion:" << format( "%20.4f" ) % ( eval_tors(*pmol)+eval_oops(*pmol) )<< " kcal/mol" << std::endl;
        }

        ctrlparm_t p(m_parm);

        if(p.igb==0)
	{
	    if( m_vec1.length()> 0 )
	    {
	        throw std::runtime_error( "Error: can't decompose PBC nonbond energy" );
            }

            numvec enb = get_dir(*pmol, p);
            os << "Direct: " << format( "%20.4f" ) % enb[0] << " kcal/mol" << std::endl;
            os << "Adjust: " << format( "%20.4f" ) % enb[1] << " kcal/mol" << std::endl;
            //double e_rec = eval_regewald( *pmol, 8.0 );
            //os << "Recip:  " << format( "%20.4f" ) % e_rec << " kcal/mol" << std::endl;
        }
	else
	{
	    numvec enb(3);
	    if( m_vec1.length()>0 )
	    {
	        atomvec_t vec1 = mask_atom(*pmol, m_vec1);
		atomvec_t vec2 = mask_atom(*pmol, m_vec2);
                enb = nonbond_egb(vec1, vec2, p);
	    }
            else
	    {
                enb = nonbond_egb(*pmol, p);
	    }

            os << "EELT:   " << format( "%20.4f" ) % enb[0] << " kcal/mol" << std::endl;
            os << "EVDW:   " << format( "%20.4f" ) % enb[1] << " kcal/mol" << std::endl;
	    os << "EGB :   " << format( "%20.4f" ) % enb[2] << " Kcal/mol" << std::endl;
        }

        //leaplog_t::putline( os.str() );

        return true;
    }

    void energy_command::undo( )
    {
    }
    
    const char* energy_command::info( ) const
    {
        return "  usage: energy unit ctrlparm";
    }

    shared_ptr< command_i > energy_command::clone( const vector<string>& args ) const
    {
        if(args.size()==3)
        {
            return command_ptr(new energy_command(args[1], args[2]));
        }
        else if(args.size()==5)
	{
	    return command_ptr(new energy_command(args[1], args[2], args[3], args[4]) );
	}
	else
	{
            throw std::runtime_error( "Error: wrong number of arguments" );
        }
    }

} // namespace amber

amber::energy_command g_energy_command;




