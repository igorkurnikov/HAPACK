#include <sstream>
#include <fstream>
#include <stdexcept>
#include <object.hpp>
#include <format.hpp>
#include <ambfmt.hpp>
#include <sstream>

#include "saveprm.hpp"

namespace amber
{
    using namespace mort;
    saveprm_command::saveprm_command( const string& name )
        :command_i( name )
    {
    }
    
    saveprm_command::saveprm_command( const string& action, const string& name, const string& top, const string& xyz )
        : m_name( name ), m_action( action ), m_topfile( top ), m_xyzfile( xyz )
    {
    }
    
    saveprm_command::~saveprm_command()
    {
    }


    string itoa( int i )
    {
        std::ostringstream os;
        os << i;
        return os.str();
    }

    void switch_amoeba_type( molecule_t& mol, const molecule_t& poleff )
    {
        int switched = 0;
        if( mol.get_i("amoeba-switched", switched) && switched>0 )
        {
            return;
        }

        atomiter_t ai = mol.atom_begin();
        for( ; ai != mol.atom_end(); ++ai )
        {
            string pole_type = ai->get_s( TYPE );
            int pole_type_id = atoi( pole_type.c_str() );

            atom_t pole( poleff, pole_type_id-1 );

            int atom_type_id = pole.get_i(TYPEID);
            string atom_type = itoa( atom_type_id );
           
            ai->set_s( TYPE, atom_type );
            ai->set_s( POLTYPE, pole_type );
	    ai->set_i( POLTYPEID, pole_type_id );
        }

        mol.set_i("amoeba-switched", 1);
    }

    bool saveprm_command::exec()
    {
        if( !content().has(m_name) )
        {
            throw std::runtime_error( "Error: molecule " + m_name + " doesn't exist." );
        }

        if( m_action == "savetinkerparm" )
        {
            molecule_ptr pmol = content().get_mol( m_name );
            assert( pmol != NULL );
            save_mol( m_xyzfile, *pmol, TINKERXYZ );
            return true;
        }

        std::ofstream os( m_topfile.c_str() );
        if( !os )
        {
            throw std::runtime_error( "Error: can't open file " + m_topfile + " for writing" );
        }


        assert( m_action == "saveamberparm" || m_action == "saveamberparmpol" || m_action == "saveamoebaparm" );

        string frcfld = (m_action=="saveamberparm"||m_action=="saveamberparmpol") ? "_amberffp" : "_amoeba-atom";
            
        parmset_t params;
             
        molecule_ptr pmol = content().get_mol( m_name );
        molecule_ptr parm = content().get_mol( frcfld );
        assert( pmol != NULL );
        if( parm == NULL )
        {
            throw std::runtime_error( "Error: no force field parameter loaden" );
        }

        if( m_action=="saveamoebaparm" )
        {
            molecule_ptr poleff = content().get_mol( "_amoeba-pole" );
            assert( poleff != NULL );
            switch_amoeba_type( *pmol, *poleff );
            write_amoeba_prmtop( os, *pmol, *parm, *poleff );
        }
        else
        {

            if( m_action=="saveamberparmpol" )
            {
                pmol->set_i( POLAR, 1 );
            }
            write_amber_prmtop( os, *pmol, *parm );
        }

                
        save_mol( m_xyzfile, *pmol, AMBERXYZ );
        return true;
    }
    
    void saveprm_command::undo()
    {
    }
    
    command_ptr saveprm_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 4 )
        {
            throw std::runtime_error( "Error: wrong number of arguments" );   
        }

        return command_ptr( new saveprm_command( args[0], args[1], args[2], args[3] ) );
    }

} // namespace amber 

amber::saveprm_command g_saveamberparm_command( "saveamberparm" );
amber::saveprm_command g_saveamoebaprm_command( "saveamoebaparm" );
amber::saveprm_command g_saveamberparmpol_command( "saveamberparmpol" );
amber::saveprm_command g_savetinkerprm_command( "savetinkerparm" );

            
