#include <sys/types.h>
#include <sys/stat.h>
#if !defined(_MSC_VER)
#include <unistd.h>
#endif

#include <fstream>
#include <object.hpp>
#include <format.hpp>
#include "xyz.hpp"

namespace mort
{
    void load_mol( istream& is, molecule_t& mol, const hashid_t& format )
    {
        if( format==SDF )
        {
            read_sdf(is, mol);
        }
        else if( format==PDB )
        {
            read_pdb(is, mol);
        }
        else if( format==OFF )
        {
            read_off(is, mol);
        }
        else if( format==MOL2 )
        {
            read_mol2(is, mol);
        }
        else
        {
            throw std::runtime_error( "Error: unknown file type " );
        }
    }

    void save_mol( ostream& os, const molecule_t& mol, const hashid_t& format )
    {
        if( format==SDF )
        {
            write_sdf( os, mol );
        }
        else if( format==PDB )
        {
            write_pdb( os, mol );
        }
        else if( format==OFF )
        {
            write_off( os, mol );
        }
        else if( format==MOL2 )
        {
            write_mol2( os, mol );
        }
        else if( format==AMBERXYZ )
        {
            write_amber_xyz( os, mol );
        }
        else if( format==TINKERXYZ )
        {
            write_tinkerxyz( os, mol );
        }
        else
        {
            throw std::runtime_error( "Error: unknown file type " );
        }
    }

    void load_mol( const string& file, molecule_t& mol, const hashid_t& format )
    {
        hashid_t fmtcpy = (format==UNKNOWN) ? get_fmt(file) : format;

        std::ifstream is( file.c_str() );
        
        if( ! is )
        {
            throw std::runtime_error( "Error: can't open file " + file + " for read." );
        }
        
        load_mol( is, mol, fmtcpy );
    }

    void save_mol( const string& file, const molecule_t& mol, const hashid_t& format )
    {
        hashid_t fmtcpy = (format==UNKNOWN) ? get_fmt(file) : format;

        std::ofstream os( file.c_str() );
        
        if( ! os )
        {
            throw std::runtime_error( "Error: can't open file " + file + " for written" );
        }
        
        save_mol( os, mol, fmtcpy );
    }

} // namespace mort



