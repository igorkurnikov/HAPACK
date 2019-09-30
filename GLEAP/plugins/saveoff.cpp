#include <fstream>
#include <object.hpp>
#include <format.hpp>
#include "saveoff.hpp"

namespace amber
{
    
    saveoff_command::saveoff_command( )
        : command_i("saveoff")
    {
    }
    
    saveoff_command::saveoff_command( const string& unit, const string& file )
        : m_unit(unit), m_file(file)
    {
    }
    
    saveoff_command::~saveoff_command( )
    {
    }
    
    bool saveoff_command::exec( )
    {
        std::ofstream os( m_file.c_str() );
        if(os == NULL)
        {
            throw std::runtime_error( "Error: can't open file " + m_file + " for write" );
        }

        entity_ptr pe = content().get(m_unit);
        if(pe == NULL)
        {
            throw std::runtime_error("Error: no such unit " + m_unit);
        }
        
        molecule_ptr pm = dynamic_pointer_cast< molecule_t >(pe);
        if( pm != NULL)
        {
            os << "!!index array str" << std::endl;
            os << " \"" <<  m_unit << "\"" << std::endl;
            write_off(os, *pm);
            return true;
        }

        database_ptr pdb = dynamic_pointer_cast< database_t >(pe);
        if( pdb != NULL )
        {
            os << "!!index array str" << std::endl;

            database_t::iterator mi, me;
            mi = pdb->begin();
            me = pdb->end();
            for( ; mi != me; ++mi )
            {
                os << " \"" <<  mi->first << "\"" << std::endl;
            }
            
            mi = pdb->begin();
            for( ; mi != me; ++mi )
            {
                molecule_ptr pm = dynamic_pointer_cast<molecule_t>(mi->second);
                if( pm != NULL )
                {
                    write_off(os, *pm);
                }
            }

            return true;
        }

        throw std::runtime_error("Error: unit " + m_unit + " is neither a molecule nor a database" );
    }
    
    void saveoff_command::undo( )
    {
    }
    
    const char* saveoff_command::info( ) const
    {
        return "  usage: saveoff unit file";
    }
    
    shared_ptr< command_i > saveoff_command::clone( const vector< string >& args ) const
    {
        if( args.size() != 3 )
        {
            for(int i=0; i < (int)args.size(); ++i)
                std::cout << args[i] << std::endl;
            throw std::runtime_error( "Error: wrong number of arguments." );
        }

        return shared_ptr< command_i >( new saveoff_command( args[1], args[2] ) );
    }

    
} // namespace amber


amber::saveoff_command g_saveoff_command;

