#include <capbox.hpp>
#include <object.hpp>
#include "solvate.hpp"

namespace amber
{
    
    solvate_command::solvate_command( const string& type )
        :command_i( "solvate" + type )
    {
    }
    
    solvate_command::solvate_command( const string& action, const string& solute, const string& solvent, 
                                      const string& center, const string& extent, double closeness )
        : m_action( action ), m_solute( solute ), m_center( center ), m_solvent( solvent ),
          m_extent( extent ), m_closeness( closeness )
    {
    }
    
    solvate_command::~solvate_command()
    {
    }
    
    bool solvate_command::exec( )
    {
        molecule_ptr psolute = content().get_mol( m_solute );
            
        molecule_ptr psolvent = content().get_mol( m_solvent );    

        assert( psolute != NULL && psolvent != NULL );

        if( m_action == "solvatebox" )
        {
            double length = atof( m_extent.c_str() );
            
            solvatebox( *psolute, *psolvent, length, m_closeness );

            return true;
        }
        
	if( m_action == "solvateoct" )
	{
	    double length = atof( m_extent.c_str() );

	    solvateoct( *psolute, *psolvent, length, m_closeness );

	    return true;
	}


        if( m_action == "solvatecap" )
        {
            numvec cnt(3);

            int ndot = std::count( m_center.begin(), m_center.end(), '.' );
            
            if( ndot == 0 )
            {
                cnt = center( *psolute );
            }
            else if( ndot == 1 )
            {
                cnt = center( content().get_resd( m_center ) );
            }
            else if( ndot == 2 )
            {
                cnt = content().get_atom( m_center ).get_v(POSITION);
            }
            else
            {
                throw std::runtime_error( "Error: can not understand mask " + m_center );
            }
            
            double radius = atof( m_extent.c_str() );
            
            solvatecap( *psolute, *psolvent, cnt, radius, m_closeness );

            return true;
        }
        
        if( m_action == "solvateshell" )
        {        
            double extent = atof( m_extent.c_str() );
            solvateshl( *psolute, *psolvent, extent, m_closeness );
            return true;
        }


        throw std::runtime_error( "Error: unknown solvate command: " + m_action );
    }
    
    void solvate_command::undo( )
    {
        throw std::runtime_error( "Sorry: not implemented yet" );
    }
    
    const char* solvate_command::info() const
    {
        return "  usage: solvateBox solute solvent buffer\n"
            "         solvateCap solute solvent center radius\n"
            "         solvateShell solute solvent buffer\n";
    }

    shared_ptr< command_i > solvate_command::clone( const vector< string >& args ) const
    {
        if( args[0]=="solvatebox" || args[0]=="solvateshell" || args[0]=="solvateoct" )
        {
            if( args.size()!=4 && args.size()!=5 )
            {
                throw std::runtime_error( "Error: wrong number of arguments" );
            }

	    double closeness = (args.size()==4) ? 1.0 : atof( args[4].c_str() ); 
            return shared_ptr< command_i >( new solvate_command( args[0], args[1], args[2], "", args[3], closeness ) );
        }
        
        assert( args[0] == "solvatecap" );
        
        if( args.size()!=5 && args.size()!=6 )
        {
            throw std::runtime_error( "Error: wrong number of arguments" );
        }
        
        double closeness = (args.size()==5) ? 1.0 : atof( args[5].c_str() ); 
        return shared_ptr< command_i >( new solvate_command( args[0], args[1], args[2], args[3], args[4], closeness ) );
    }
    
} // namespace amber

amber::solvate_command g_solvatebox_command( "box" );
amber::solvate_command g_solvateoct_command( "oct" );
amber::solvate_command g_solvatecap_command( "cap" );
amber::solvate_command g_solvateshl_command( "shell" );


            
