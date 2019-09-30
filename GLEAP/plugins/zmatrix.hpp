#ifndef GLEAP_PLUGINS_ZMATRIX_HPP
#define GLEAP_PLUGINS_ZMATRIX_HPP

#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class zmatrix_command : public command_i
    {
    public:

        zmatrix_command( );

        zmatrix_command( const string& object, const vector< vector< string > >& inters );

        virtual ~zmatrix_command( );

        virtual bool exec( );

        virtual void undo( );

        virtual const char* info( ) const;

        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;

    private:

        string m_object;

        vector< vector< string > > m_inters;
    };

} // namespace amber

#endif
