#include <sstream>
#include <fstream>
#include <format.hpp>

#include "parse.hpp"
#include "subst.hpp"
#include "reader.hpp"

namespace mort
{
    bool read_smiles( const char* exp, molecule_t& mol )
    {
        static const daylight::reader_set readers( SMILES );

        molecule_t temp;
        
        daylight::parse( exp, temp, readers, SMILES );

        mol.swap( temp );

        return true;
    }
    
    bool read_smarts( const char* exp, molecule_t& mol )
    {        
        static const daylight::reader_set readers( SMARTS );

        molecule_t temp;
        
        daylight::parse( exp, temp, readers, SMARTS );

        mol.swap( temp );

	return true;
    }
    
} // namespace mort
