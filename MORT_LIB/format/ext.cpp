#include <common.hpp>

namespace mort
{

    bool is_dbfmt( const hashid_t& format )
    {
        return format==SDF || format==OFF || format==MOL2;
    }

    string get_ext( const string& file )
    {
        int last_dot = file.find_last_of( '.' );
        if( last_dot == (int)string::npos )
        {
            return "";
        }

        return file.substr( last_dot + 1, file.length() - last_dot - 1 );
    }

    hashid_t get_fmt( const string& file )
    {
        string ext = get_ext( file );
        
        if( ext == "sd" || ext == "sdf" || ext == "mol" )
        {
            return SDF;
        }

        if( ext == "pdb" || ext == "ent" )
        {
            return PDB;
        }

        if( ext == "off" || ext == "lib" )
        {
            return OFF;
        }
        
        if( ext == "mol2" )
        {
            return MOL2;
        }        
        
        throw std::runtime_error( "Error: unknown file name extension: " + ext );
    }

} // namespace mort
 
