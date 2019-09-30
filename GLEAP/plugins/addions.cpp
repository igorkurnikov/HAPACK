#include <object.hpp>
#include <capbox.hpp>
#include "addions.hpp"

namespace amber
{
    
    addions_command::addions_command( )
        :command_i( "addions" )
    {
    }
    
    addions_command::addions_command( const string& dest, const string& ion, int number )
        :m_dest( dest ), m_ion( ion ), m_number( number )
    {
    }
    
    addions_command::~addions_command()
    {
    }
    
    bool addions_command::exec()
    {
        molecule_ptr pdst = content().get_mol( m_dest );
        molecule_ptr pion = content().get_mol( m_ion );

        if( m_number==0 )
        {
            double dstchg = charge( *pdst );
            double ionchg = charge( *pion );
            if( dstchg * ionchg >= 0 )
            {
                throw std::runtime_error( "Error: molecule and ion are both positive (or negative), can not make neutral" );
            }
            
            m_number = int( -dstchg/ionchg+0.5);
        }
        
        
        assert( pdst != NULL && pion != NULL );

        double shell_extent = 4.0;
        double resolution   = 1.0;

        addions( *pdst, *pion, m_number, shell_extent, resolution);
	return true;
    }
    
    void addions_command::undo()
    {
        throw std::runtime_error( "Sorry: not implemented yet" );
    }

    const char* addions_command::info( ) const
    {
        return "  usage: addion mol ion1 numIon1 ";
    }
    
    shared_ptr< command_i > addions_command::clone( const vector< string >& args ) const
    {
        if( args.size()!=3 && args.size()!=4 )
        {
            throw std::runtime_error( "Error: wrong number of arguments " );
        }

        int numb = args.size()==3 ? 0 : atoi(args[3].c_str());
    
        return shared_ptr< command_i >( new addions_command(args[1], args[2], numb) );
    }

    addions_command g_addions_command;
    
} // namespace amber



