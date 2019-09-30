#include <sstream>
#include <object.hpp>
#include "set.hpp"

namespace amber
{
    numvec atov(const string& a, int len)
    {
        std::istringstream is(a);

        char c;
	is >> c;

	if( c != '{' )
	{
	    is.putback(c);
	}

	numvec v(len);
	for( int i=0; i<len; ++i )
	{
	    is >> v[i];
	}

	return v;
    }

    numvec atov3d(const string& a)
    {
        return atov(a, 3);
    }

    numvec atov4d(const string& a)
    {
        return atov(a, 4);
    }

    template< typename T >
    bool set_parm( T& obj, const string& parm, const string& value )
    {       
        if( parm == "name" )
        {
            obj.set_s(NAME, value);
            return true;
        }
        else if( parm == "type" )
        {
            obj.set_s(TYPE, value);
            return true;
        }
        else if( parm == "charge" )
        {
            obj.set_d(PCHG, atof( value.c_str() ) );
            return true;
        }
        else if( parm == "element" )
        {
            obj.set_i(ELEMENT, pertab_t::get_element( value.c_str() ) );
            return true;
        }
        else if( parm == "position" )
        {
	    numvec pos = atov3d(value);
            obj.set_v(POSITION, atov3d(value) );
            return true;
        }
        else if( parm == "connect0" )
        {
	    atom_t a = content().get_atom(value);
            obj.set_i(HEAD, a.get_i(ID) );
            return true;
        }
        else if( parm == "connect1" )
        {
	    atom_t a = content().get_atom(value);
            obj.set_i(TAIL, a.get_i(ID) );
            return true;
        }
        else if( parm == "connect2" )
        {
	    atom_t a = content().get_atom(value);
            obj.set_i(CONN1, a.get_i(ID) );
            return true;
        }
        else if( parm == "head" )
        {
	    atom_t a = content().get_atom( value );
            obj.set_i(HEAD, a.get_i(ID) );
            return true;
        }
        else if( parm == "tail" )
        {
	    atom_t a = content().get_atom( value );
            resd_t r = a.resd();
            obj.set_i("tailresd", r.absid()+1 );
            r.set_i(TAIL, a.absid()+1 );
            return true;
        }
        else if( parm == "box" )
        {
            obj.set_i(SOLUTE, BOX);
            obj.set_v(BOX, atov4d(value) );
            return true;
        }
        else if( parm == "cap" )
        {
	    obj.set_i(SOLUTE, CAP);
            obj.set_v(CAP, atov4d(value) );
            return true;
        }

        return false;
    }

    set_command::set_command( )
        : command_i( "set" )
    {
    }
    
    set_command::set_command( const string& object, const string& parm, const string& value )
        : m_object( object ), m_parm( parm ), m_value( value )
    {
    }
    
    set_command::~set_command( )
    {
    }
    
    bool set_command::exec( )
    {
        if( m_object == "default" )
        {
            mortenv().set_s( m_parm.c_str(), m_value.c_str() );
            return true;
        }

        int ndot = std::count( m_object.begin(), m_object.end(), '.' );
        
        if( ndot == 0 )
        {
            entity_ptr obj = content().get( m_object );

            if( ! set_parm( *obj, m_parm, m_value ) )
            {
                obj->set_s( m_parm, m_value );
            }   
        }
        else if( ndot == 1 )
        {
            resd_t resd = content().get_resd( m_object );
            
            if( ! set_parm( resd, m_parm, m_value ) )
            {
                throw std::runtime_error( "Error: can't set paramerter " + m_parm + " for a residue" );
            }
        }
        else
        {
            atom_t atom = content().get_atom( m_object );
            
            if( ! set_parm( atom, m_parm, m_value ) )
            {
                throw std::runtime_error( "Error: can't set parameter " + m_parm + " for an atom" );
            }
        }

	return true;
    }

    void set_command::undo( )
    {
        throw std::runtime_error( "sorry, undo not (yet) implemented"  );
    }

    const char* set_command::info( ) const
    {
        return "  usage: set default variable value \n"
            "                        or \n "
            "            set variable parameter value\n";    
    }

    shared_ptr< command_i > set_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 4 )
        {
            throw std::runtime_error( "Error: wrong number of arguments" );
        }
        
        return shared_ptr< command_i >( new set_command( args[1], args[2], args[3] ) );
    }       
        

} // namespace amber

amber::set_command g_set_command;


