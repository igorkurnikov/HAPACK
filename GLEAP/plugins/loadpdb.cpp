#include <fstream>
#include <object.hpp>
#include <format.hpp>
#include <pdbent.hpp>

#include "loadpdb.hpp"

namespace amber
{
    loadpdb_command::loadpdb_command(const string& cmdname)
        :command_i(cmdname)
    { 
    }

    loadpdb_command::loadpdb_command(const string& action, const string& name, const string& file, const string& seq)
        :m_action(action), m_name(name), m_file(file), m_seq(seq)
    {
    }

    loadpdb_command::~loadpdb_command()
    {
    }


    bool loadpdb_command::exec( )
    {
        std::ifstream is( m_file.c_str() );

        if( ! is )
        {
            throw std::runtime_error( "Error: file " + m_file + " doesn't exist" );
        }

        molecule_ptr pmol( new molecule_t() );

        read_pdb( is, *pmol );

        string disulf = "auto";
	mortenv().get_s("disulfide", disulf);
        if( disulf!="off" )
        {
            console_t* pcon = (disulf=="auto") ? NULL : console().get();

            double disulfcut = 2.1;
            string tmp;
            if( mortenv().get_s("disulfcut", tmp) )
            {
	        disulfcut = atof( tmp.c_str() );
            }

            disulfide( *pmol, disulfcut, pcon );
        }

        namemap_ptr nmap = content().get_nmap( "_namemap" );
        if( nmap != NULL )
        {
            if( m_action=="loadpdb" )
            {
                mdlize_mdb( *pmol, content() );
            }
            else 
            {
                assert( m_action =="loadpdbusingseq" );

                database_ptr seq = content().get_mdb( m_seq );
                if( seq == NULL )
                {
                    throw std::runtime_error("Error: undefined sequence " + m_seq );
                }
            
                mdlize_seq( *pmol, *seq, *nmap );
            }
        }

        pmol->set_s(NAME, m_name);
    
        content().set( m_name, pmol );

        return true;        
    }

    void loadpdb_command::undo( )
    {
        content().remove( m_name.c_str() );    
    }

    shared_ptr< command_i > loadpdb_command::clone( const vector< string >& args ) const
    {
        if( args.size()!=3 && args.size()!=4 )
        {
            throw std::runtime_error( "Error: wrong arguments " + args[0] + " " + args[1] + " " + args[2] );
        }

        string seq = (args.size()==3)? "" : args[3];

        return command_ptr( new loadpdb_command(args[0], args[1], args[2], seq) );
    }

}

static amber::loadpdb_command g_loadpdb_command( "loadpdb" );
static amber::loadpdb_command g_loadpdbusingseq_command( "loadpdbusingseq" );


