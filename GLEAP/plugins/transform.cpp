#include <boost/algorithm/string.hpp>
#include <common.hpp>
#include <object.hpp>
#include "transform.hpp"

namespace amber
{
    using std::runtime_error;

    atomvec_t get_atomlist( const string& obj )
    {
        int ndot = std::count( obj.begin(), obj.end(), '.' );

        atomvec_t r;
        if( ndot==0 )
        {
            molecule_ptr pmol = content().get_mol( obj );
            r.insert( r.end(), pmol->atom_begin(), pmol->atom_end() );
        }
        else if( ndot==1 )
        {
            r.push_atom( content().get_resd(obj) );
        }
        else if( ndot==2 )
        {
            r.push_back( content().get_atom(obj) );
        }
        return r;
    }

    void transform( const matrix& m, numvec& v )
    {
        assert( m.size1()==4 && m.size2()==4 );
        numvec r(3);
        r[0] = m(0,0)*v[0] + m(0,1)*v[1] + m(0,2)*v[2];
        r[1] = m(1,0)*v[0] + m(1,1)*v[1] + m(1,2)*v[2];
        r[2] = m(2,0)*v[0] + m(2,1)*v[1] + m(2,2)*v[2];
        v[0] = r[0] * m(3,3);
        v[1] = r[1] * m(3,3);
        v[2] = r[2] * m(3,3);
    }

    
    transform_command::transform_command( const string& action )
        : command_i( action )
    {
    }
    
    transform_command::transform_command( const string& action, const string& object, const string& offset )
        : m_action( action ), m_object( object ), m_offset( offset )
    {
    }
    
    transform_command::~transform_command( )
    {
    }
    
    bool transform_command::exec( )
    {
        atomvec_t atoms = get_atomlist( m_object );

        vector< string > values;
        
        split( values, m_offset, is_any_of( " {}" ), token_compress_on );
        
        if( m_action == "translate" )
        {
            if( values.size() != 5 )
            {
                throw runtime_error( "Error: the second argument of translate must be a 3d vector" );
            }
            
            numvec offset = makevec( atof( values[1].c_str() ),
                             atof( values[2].c_str() ), 
                             atof( values[3].c_str() ) );

            for( unsigned int i=0; i < atoms.size(); ++i )
            {
                numvec pos = atoms[i].get_v(POSITION);
                pos += offset;
                atoms[i].set_v(POSITION, pos );
            }

            return true;
        }
        
        if( m_action == "rotate" )
        {
            assert( values.size() == 6 );
            
            numvec axis = makevec( atof( values[1].c_str() ),
                            atof( values[2].c_str() ),
                            atof( values[3].c_str() ) );
            double degree = atof( values[4].c_str() );
            
            for( unsigned int i=0; i < atoms.size(); ++i )
            {
                numvec pos = atoms[i].get_v(POSITION);
                rotate( pos, axis, degree );
                atoms[i].set_v(POSITION, pos );
            }

            return true;
        }
        
        if( m_action == "transform" )
        {
            if( values.size()!=11 && values.size()!=18 )
            {
                throw runtime_error( "Error: transform matrix should have either 9 or 16 elements");
            }
            
            matrix mat(4,4);
            
            int n = (values.size() == 11) ? 3 : 4;
            int k = 1;
            for( int i=0; i < n; ++i )
            {
                for( int j=0; j < n; ++j )
                {
                    mat(i,j)  = atof( values[k].c_str() );
                    k++;
                }
            }
            
            if( n == 3 ) mat(3,3) = 1.0;

            for( unsigned int i=0; i < atoms.size(); ++i )
            {
                numvec pos = atoms[i].get_v(POSITION);
                transform( mat, pos );
                atoms[i].set_v(POSITION, pos );
            }

            return true;
        }
    
        throw runtime_error( "Error: unknown action for command transform: " + m_action );
    }

    void transform_command::undo( )
    {
    }

    const char* transform_command::info( ) const
    {
        return "  usage: transform atoms matrix\n    or\n rotate atoms vector degree\n or\n translate atoms vector";
    }

    shared_ptr< command_i > transform_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 3 )
        {
            throw runtime_error( "Error: wrong number of argument" );
        }
        
        return shared_ptr< command_i >( new transform_command( args[0], args[1], args[2] ) );
    }
    
} // namespace amber


