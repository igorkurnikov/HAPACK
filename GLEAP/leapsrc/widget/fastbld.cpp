#include <common.hpp>
#include <object.hpp>
#include <pdbent.hpp>
#include "fastbld.hpp"


namespace mortgtk
{
    class fastbld_command : public command_i
    {
    public:

        fastbld_command()
        {
        }
        
        fastbld_command( const string& mname )
            : m_mname( mname )
        {
        }
        
        virtual bool exec()
        {
            molecule_ptr pmol = content().get_mol( m_mname );
            molecule_ptr pfrc = content().get_mol( "_amber-atom" );
            namemap_ptr pnmap = content().get_nmap( "_namemap" );
            assert( pmol != NULL );
            assert( pnmap != NULL );

            mortenv().set_s( "fixbo", "on" );
            
            build_bymdl( *pmol, content(), *pnmap, pfrc.get() );
        }

        virtual void undo()
        {
        }
        
        virtual command_ptr clone(const vector<string>& args) const
        {
            return command_ptr( new fastbld_command(args[1]) );
        }
        
        
    private:

        string m_mname;
        
    };
    
        
    fastbld_t::fastbld_t()
        : basicdlg_t( "Build up a biopolymer in fastest way" )
    {
        dlgfactory_t::add( "fastbld", this );
    }
    

    fastbld_t::~fastbld_t()
    {
    }
    
    command_ptr fastbld_t::get_command()
    {
        return command_ptr( new fastbld_command( molname() ) );
    }
}
