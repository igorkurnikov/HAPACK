#include <object.hpp>
#include "bond.hpp"

namespace amber
{
    void connect( root_t& atom1, root_t& atom2, int order )
    {
        if( atom1.has_s(NBRLIST) )
        {
            assert( atom1.has_s(NBRLIST) );
            string nbrlist = atom1.get_s(NBRLIST);
            string nbrtype = atom1.get_s(NBRTYPE);
                
            nbrlist += atom2.get_s(NAME) + " ";
            nbrtype.append( 1, char('0'+order) );
            nbrtype.append( 1, char(' ') );
   
            atom1.set_s(NBRLIST, nbrlist);
            atom1.set_s(NBRTYPE, nbrtype);
        }
        else
        {
            atom1.set_s(NBRLIST, atom2.get_s(NAME) + " " );
            string nbrtype;
            nbrtype.append( 1, char('0'+order) );
            nbrtype.append( 1, char(' ') );
            atom1.set_s(NBRTYPE, nbrtype );
        }
    }

    bond_command::bond_command( )
        : command_i( "bond" )
    {
    }
    
    bond_command::bond_command( const string& atom1, const string& atom2, int order )
        :m_atom1( atom1 ), m_atom2( atom2 )
    {
        m_order = order;
    }
    
    bond_command::~bond_command( )
    {
    }

    bool bond_command::exec( )
    {
        if( m_atom1.find( '.' ) != string::npos )
        {
            atom_t atom1 = content().get_atom(m_atom1);
            atom_t atom2 = content().get_atom(m_atom2);
        
            if( &atom1.getmol() != &atom2.getmol() )
            {
                throw std::runtime_error( "Error: make bond between atoms from different molecule: " + m_atom1 + " " + m_atom2 );
            }
        
            bond_t::create( atom1, atom2 ).set_i(ORDER, m_order);
        }
        else
        {
            entity_ptr atom1 = content().get( m_atom1 );
            entity_ptr atom2 = content().get( m_atom2 );
            connect( *atom1, *atom2, m_order );
            connect( *atom2, *atom1, m_order );
        }

	return true;
    }
    
    void bond_command::undo( )
    {
        throw std::runtime_error( "sorry not implemented yet" );
    }
    
    const char* bond_command::info( ) const
    {
        return "  usage: bond atom1 atom2 [order]";
    }
    
    shared_ptr< command_i > bond_command::clone( const vector< string >& args ) const
    {
        int order;
        
        if( args.size() == 3 )
        {
            order = 1;
        }
        else if( args.size() == 4 )
        {
            if( isdigit( args[3][0] ) ) 
            {
                order = atoi( args[2].c_str() );
            }
            else if( args[3] == "S" )
            {
                order = 1;
            }
            else if( args[3] == "D" )
            {
                order = 2;
            }
            else if( args[3] == "T" )
            {
                order = 3;
            }
            else if( args[3] == "A" )
            {
                order = 4;
            }
            else
            {
                throw std::runtime_error( "Error: unkown tag for order: " + args[3] );
            }    
        }
        else
        {
            throw std::runtime_error( string( "Error: wrong number of arguments\n " ) + info() );
        }
        
        return shared_ptr< command_i >( new bond_command( args[1], args[2], order ) );
    }

    bondbydis_command::bondbydis_command( )
        : command_i( "bondbydistance" )
    {
    }
    
    bondbydis_command::bondbydis_command( const string& object, double cutoff )
        : m_object( object ), m_cutoff( cutoff )
    {
    }
    
    bondbydis_command::~bondbydis_command( )
    {
    }
    
    bool bondbydis_command::exec( )
    {
        int dot = m_object.find( '.' );
        
        if( dot == (int)string::npos )
        {
            molecule_ptr pmol = content().get_mol(m_object);
            assert( pmol != NULL );
            bond_bydis( *pmol, m_cutoff );
        }
        else
        {
            resd_t resd = content().get_resd( m_object );
            bond_bydis( resd, m_cutoff );
        }

	return true;
    }
    
    void bondbydis_command::undo( )
    {
        throw std::runtime_error( "Sorry, not implemented yet" );
    }

    const char* bondbydis_command::info( ) const
    {
        return " usage: bondbydistance container [cutoff]";
    }
    
    shared_ptr< command_i > bondbydis_command::clone( const vector< string >& args ) const
    {
        if( args.size() < 2 || args.size() > 3 )
        {
            throw std::runtime_error( string( "Error: wrong number of argument " ) + info() );
        }

        double cutoff = ( args.size() == 2 ) ? 1.6 : atof( args[2].c_str() );

        return shared_ptr< command_i >( new bondbydis_command( args[1], cutoff ) );
    }
    
} // namespace amber

amber::bond_command g_bond_command;
amber::bondbydis_command g_bondbydis_command;

