#include <object.hpp>
#include <pdbent.hpp>
#include <guilib.hpp>
#include "addmap.hpp"

namespace amber
{

    addmap_command::addmap_command(const string& type)
        : command_i(type)
    {
    }

    addmap_command::addmap_command(const string& type, const string& snmap)
        : m_type(type), m_nmap( interpret2d(snmap) )
    {
    } 

    addmap_command::~addmap_command()
    {
    }

    bool addmap_command::exec()
    {
        namemap_ptr nmap = content().get_nmap("_namemap");
	if( nmap == NULL )
	{
	    nmap = namemap_ptr( new namemap_t() );
            content().set("_namemap", nmap);
        }
 
        // std::cout << "m_type:" << m_type << std::endl;

        if( m_type=="addpdbresmap")
	{
            for(unsigned int i=0; i < m_nmap.size(); ++i )
                nmap->add_resd_map( m_nmap[i] );
	}
	else
	{
	    assert( m_type == "addpdbatommap" );
            for(unsigned int i=0; i < m_nmap.size(); ++i )
                nmap->add_atom_map( m_nmap[i] );
        }

	return true;
    }

    void addmap_command::undo()
    {
    }

    const char* addmap_command::info() const
    {
        return "addPdbResMap { {[0|1] nameA nameB} ...} \n addPdbAtomMap { {nameA nameB} ... }";
    }

    shared_ptr<command_i> addmap_command::clone(const vector<string>& args) const
    {
        assert( args.size()==2 );
        return shared_ptr<command_i>( new addmap_command(args[0], args[1]) );
    }

} // namespace amber

