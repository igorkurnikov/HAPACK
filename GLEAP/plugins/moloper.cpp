#include <object.hpp>
#include "moloper.hpp"

namespace amber
{
    
    moloper_command::moloper_command( const string& action )
        : command_i( action )
    {
    }
    
    moloper_command::moloper_command( const string& action, const string& molname )
        : m_action( action )
    {
        m_pmol = content().get_mol( molname );
        
        if( m_pmol==NULL )
        {
            throw std::runtime_error( "Error: no such a molecule: " + molname );
        }
    }
    
    moloper_command::~moloper_command()
    {
    }

    bool moloper_command::exec()
    {
        assert( m_pmol !=NULL );

        namemap_ptr nmap = content().get_nmap( "_namemap" );
        if( nmap==NULL )
        {
            throw std::runtime_error( "Error: name map not loaden" );
        }

        int nresd = m_pmol->nresd();
        resditer_t ri = m_pmol->resd_begin();
        for( ; ri != m_pmol->resd_end(); ++ri )
        {
            string type = ri->get_s(TYPE);

            try
            {
                molecule_ptr mdl = content().get_mol( type );
                continue;
            }
            catch( std::exception& )
            {
            }

            std::cout << m_action << " residue " << type << std::endl;

            molecule_t tmp;
            if( m_pmol->nresd() > 1 )
            {
                merge( tmp, *ri ); 
            }

            if( m_action=="fixbond" )
            {
                if( nresd==1 ) 
                {
                     
                    fixbond( *m_pmol );
                }
                else
                {
                    fixbond( tmp );
                }
            }
            else if( m_action=="addhydr" )
            {
                if( nresd==1 )
                {
                    addHs( *m_pmol );
                }
                else
                {
                    addHs( tmp );
                }
            }
            else if( m_action=="setpchg" )
            {
                if( nresd==1 )
                {
                    setpchg( *m_pmol );
                }
                else
                {
                    setpchg( tmp );
                }
            }
            else
            {
                assert( m_action=="parmchk" );

                molecule_ptr pff = content().get_mol( "_amberffp" );
        
                if( pff==NULL )
                {
                    throw std::runtime_error( "Error: can not find pre-loaden amber force field parameter" );
                }

                if( nresd==1 )
                {
                    parmchk( *m_pmol, *pff );
                }
                else
                {
                    parmchk( tmp, *pff );
                }
            }
        }
            
        return true;
    }

    
    void moloper_command::undo()
    {
    }
    
    command_ptr moloper_command::clone( const vector<string>& args ) const
    {
        if( args.size()!=2 )
        {
            throw std::runtime_error( "Error: wrong number of argument" );
        }
        
        return command_ptr( new moloper_command( args[0], args[1] ) );
    }
    

} // namespace amber
