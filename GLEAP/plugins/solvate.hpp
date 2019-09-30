#ifndef GLEAP_PLUGINS_SOLVATE_HPP
#define GLEAP_PLUGINS_SOLVATE_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class solvate_command : public command_i
    {
    public:

        solvate_command( const string& type );
        
        solvate_command( const string& action, const string& solute, const string& solvent, const string& center, const string& extent, double closeness );
        
        virtual ~solvate_command( );
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
        
    private:

        string m_action;

        string m_solute;
        
        string m_center;
        
        string m_solvent;
        
        string m_extent;

	double m_closeness;
    };

} // namespace amber

#endif
