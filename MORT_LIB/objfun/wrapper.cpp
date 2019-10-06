#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>

#include <object.hpp>
#include <ambfmt.hpp>
#include <format.hpp>
#include "parmfun.hpp"

#if defined(_MSC_VER)
#include <direct.h>

	static int mkdir(const char* dir_name, int mode_t)
	{
		return _mkdir(dir_name);
	}

#endif

namespace mort
{
    using std::runtime_error;

    void setpchg( molecule_t& m )
    {
        char* amhome = getenv( "AMBERHOME" );
        if( amhome==NULL )
        {
            throw runtime_error( "Error: environment AMBERHOME not set" );
        }

        // prepare input for antechamber
        mkdir( ".antechamber", 0755 );
        save_mol( ".antechamber/ac_tmp_in.mol2", m );

        // execute antechamber
        chdir( ".antechamber" );
        string exe(amhome);
        exe += "/bin/antechamber";
        string cmd = "\"" + exe + "\" -i ac_tmp_in.mol2 -fi mol2 -o ac_tmp_ot.mol2 -fo mol2 -c bcc";
        if( system(cmd.c_str()) !=0 )
        {
            throw runtime_error( "Error: encountered problem while running antechamber" );
        }
        chdir( ".." );
            
        // read the output.
        molecule_t o;
        load_mol( ".antechamber/ac_tmp_ot.mol2", o );
        if( o.natom() != m.natom() )
        {
            throw runtime_error( "Error: the output mol has different size from the input mol." );
        }

	m.atoms()[0].set_d(PCHG, 0.0);

        // copy the output.
        vector<string>& typsrc = get_svec( o, ATOM, TYPE );
        vector<double>& chgsrc = get_dvec( o, ATOM, PCHG );
        vector<string>& typdst = get_svec( m, ATOM, TYPE );
        vector<double>& chgdst = get_dvec( m, ATOM, PCHG );
        std::copy( typsrc.begin(), typsrc.end(), typdst.begin() );
        std::copy( chgsrc.begin(), chgsrc.end(), chgdst.begin() );
    }
    
    void parmchk( molecule_t& m, molecule_t& ff )
    {
        // ACHOME does not exist in Amber10 or 11; SRB July 2008.
        // ACHOME is put here temporarily for chimera, which has a special 
        // version of antechamber with it.
        char* achome = getenv( "AC_HOME" );
        char* amhome = (achome!=NULL) ? achome : getenv( "AMBERHOME" );
        if( amhome == NULL )
        {
            throw std::runtime_error( "Error: neither AC_HOME nor AMBERHOME has been set, cannot find parmchk executable" );
        }

        // prepare input file for parmchk
        mkdir( ".parmchk", 0755 );
        save_mol( ".parmchk/ac_tmp_in.mol2", m );

        string exe(amhome);
        exe += "/bin/parmchk";
        string cmd = "\"" + exe + "\" -i .parmchk/ac_tmp_in.mol2 -f mol2 -o .parmchk/frcmod.tmp";
        if( system( cmd.c_str() )!=0) 
        {
            throw runtime_error( "Error: encoutered problem while running parmchk" );
        }
        
        // read the output
        load_frc( ".parmchk/frcmod.tmp", ff );
    }

} // namespace mort

