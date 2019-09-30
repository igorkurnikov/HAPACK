#include <guilib.hpp>

namespace amber
{
    using namespace mort;
   
    class bond_command : public command_i
    {
    public:

        bond_command( );

        bond_command( const string& atom1, const string& atom2, int order );
        
        virtual ~bond_command( );
        
        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;
        
        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;
        
    private:

        string m_atom1;
        
        string m_atom2;
        
        int m_order;
    
    };

    class bondbydis_command : public command_i
    {
    public:

        bondbydis_command( );
        
        bondbydis_command( const string& object, double cutoff );
        
        virtual ~bondbydis_command( );

        virtual bool exec( );
        
        virtual void undo( );
        
        virtual const char* info( ) const;

        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;

    private:

        string m_object;
        
        double m_cutoff;
    };
    

} // namespace amber
