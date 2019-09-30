#include <sstream>
#include <object.hpp>
#include "measure.hpp"
#include <sstream>

namespace amber
{
    using namespace std;

    measure_command::measure_command( )
        :command_i( "measuregeom" )
    {
    }
    
    measure_command::measure_command( const vector< string >& atoms )
        : m_atoms( atoms )
    {
    }
    
    measure_command::~measure_command( )
    {
    }
    
    bool measure_command::exec( )
    {
	for( int i=0; i < (int)m_atoms.size(); i++ )
        {
            int ndot = count( m_atoms[i].begin(), m_atoms[i].end(), '.' );

            if( ndot != 2 )
            {
                throw logic_error( m_atoms[i] + " is not a valid atom" );
            }
        }

        if( m_atoms.size() == 2 )
        {
            morf_t atom0 = content().get_atom( m_atoms[0] );
            morf_t atom1 = content().get_atom( m_atoms[1] );
            double dist  = mort::dist( atom0, atom1 );
            ostringstream os;
            os << "Distance: " << dist << " angstroms.";
            leaplog_t::putline( os.str() );
            return true;
        }
        
        if( m_atoms.size() == 3 )
        {
            morf_t atom0 = content().get_atom( m_atoms[0] );
            morf_t atom1 = content().get_atom( m_atoms[1] );
            morf_t atom2 = content().get_atom( m_atoms[2] );

            double value = angl( atom0, atom1, atom2 );
            ostringstream os;
            os << "Angle: " << value << " degrees.";
            leaplog_t::putline( os.str() );
            return true;
        }

        if( m_atoms.size() == 4 )
        {
            morf_t atom0 = content().get_atom( m_atoms[0] );
            morf_t atom1 = content().get_atom( m_atoms[1] );
            morf_t atom2 = content().get_atom( m_atoms[2] );
            morf_t atom3 = content().get_atom( m_atoms[3] );

            double value = tors( atom0, atom1, atom2, atom3 );
            ostringstream os;
            os << "Torsion angle: " << value << " degrees.";
            leaplog_t::putline( os.str() );
            return true;
        }

        return true;
    }
    
    void measure_command::undo( )
    {
        
    }

    const char* measure_command::info( ) const
    {
        return "  usage: measureGemo atom1 atom2 [ atom3 [ atom4 ] ] ";
    }
    
    shared_ptr< command_i > measure_command::clone( const vector< string >& args ) const
    {
        if( args.size() > 5 || args.size() < 3 )
        {
            throw logic_error( "Error: wrong number of arguments" );
        }
        
        return shared_ptr< command_i >( new measure_command( vector< string >( args.begin() + 1 , args.end() ) ) );
    }
    
} // namespace amber

amber::measure_command g_measure_command;

