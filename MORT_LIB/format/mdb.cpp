#include <fstream>
#include <common.hpp>
#include <object.hpp>
#include "mol.hpp"
#include "ext.hpp"

namespace mort
{


    void load_mdb( istream& is, database_t& db, const hashid_t& format )
    {
        while( is )
        {
            shared_ptr< molecule_t > pmol( new molecule_t() );
            load_mol( is, *pmol, format );
            if( pmol->natom() == 0 )
            {
                continue;
            }
            
            string name( "noname" );
            pmol->get_s(NAME,name);
            db.add( name, pmol );
        } 

    }    


    void save_mdb( ostream& os, const database_t& db, const hashid_t& format )
    {
        database_t::const_iterator i = db.begin();
        
        for( ; i != db.end(); ++i )
        {           
            shared_ptr< molecule_t > pmol = dynamic_pointer_cast< molecule_t >( i->second );

            if( pmol )
            {
                save_mol( os, *pmol, format );
            }
        }
    }
    
    void load_mdb( const string& file, database_t& db, const hashid_t& format )
    {
        hashid_t fmtcpy = (format==UNKNOWN) ? get_fmt(file) : format;

        std::ifstream is( file.c_str() );
        
        if( ! is )
        {
            throw std::runtime_error( "Error: can not open file " + file + " for read" );
        }

        load_mdb( is, db, fmtcpy );
    }
    
    void save_mdb( const string& file, const database_t& db, const hashid_t& format )
    {
        hashid_t fmtcpy = (format==UNKNOWN) ? get_fmt(file) : format;
       
        std::ofstream os( file.c_str() );
        
        if( !os )
        {
            throw std::runtime_error( "can't open file " + file + " for written" );
        }

        save_mdb( os, db, fmtcpy );
    }

} // namespace mort
 
