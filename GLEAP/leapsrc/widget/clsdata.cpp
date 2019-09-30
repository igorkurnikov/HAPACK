#include <mortgl.hpp>
#include "mainwin.hpp"
#include "clsdata.hpp"

namespace mortgtk
{

    class clsdata_command : public command_i
    {
    public:

        clsdata_command()
        {
        }

        clsdata_command( const string& molname )
            : m_molname( molname )
        {
        }

        virtual bool exec()
        {
            content().remove( m_molname );

            drawing()->remove( m_molname );
        }

        virtual void undo()
        {
        }

        virtual command_ptr clone( const vector<string>& args ) const
        {
            assert( args.size()==2 );
            return command_ptr( new clsdata_command(args[1]) );
        }
            
    private:

        string m_molname;
    };
  
    clsdata_t::clsdata_t()
        : basicdlg_t( "Close data" )
    {
        dlgfactory_t::add( "clsdata", this );
    }

    clsdata_t::~clsdata_t()
    {
    }

    command_ptr clsdata_t::get_command()
    {
        return command_ptr( new clsdata_command( molname() ) );
    }

    basicdlg_ptr clsdata_t::clone() const
    {
        return basicdlg_ptr( new clsdata_t() );
    }

}

