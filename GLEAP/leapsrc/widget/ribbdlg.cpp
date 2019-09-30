#include <mortgl.hpp>
#include "mainwin.hpp"
#include "ribbdlg.hpp"

namespace mortgtk
{

    class ribbon_command : public command_i
    {
    public:

        ribbon_command()
        {
        }

        ribbon_command( const string& molname )
            : m_molname( molname )
        {
        }

        virtual bool exec()
        {
	    molecule_ptr pmol = content().get_mol( m_molname );

            drawing()->add( "ribbon", *pmol, "");
        }

        virtual void undo()
        {
        }

        virtual command_ptr clone( const vector<string>& args ) const
        {
            assert( args.size()==2 );
            return command_ptr( new ribbon_command(args[1]) );
        }
            
    private:

        string m_molname;
    };
  
    ribbdlg_t::ribbdlg_t()
        : basicdlg_t( "Close data" )
    {
        dlgfactory_t::add( "ribbdlg", this );
    }

    ribbdlg_t::~ribbdlg_t()
    {
    }

    command_ptr ribbdlg_t::get_command()
    {
        return command_ptr( new ribbon_command( molname() ) );
    }

    basicdlg_ptr ribbdlg_t::clone() const
    {
        return basicdlg_ptr( new ribbdlg_t() );
    }

}

