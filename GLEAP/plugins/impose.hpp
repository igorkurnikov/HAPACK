#ifndef GLEAP_PLUGINS_IMPOSE_HPP
#define GLEAP_PLUGINS_IMPOSE_HPP
#include <common.hpp>
#include <object.hpp>
#include <guilib.hpp>

namespace amber
{
    using namespace mort;

    class impose_command : public command_i
    {
    public:

        impose_command( );

        impose_command( const string& unit, const vector< string >& units, const vector< vector< string > >& inters );

        virtual ~impose_command( );

        virtual bool exec( );

        virtual void undo( );

        virtual const char* info( ) const;

        virtual shared_ptr< command_i > clone( const vector< string >& args ) const;

    private:

        string m_unit;

        vector< string > m_resds;

        vector< vector< string > > m_inters;

    };

    void impose_dist( const atomvec_t& atoms, double dist );

    void impose_angl( const atomvec_t& atoms, double angl );

    void impose_angl( const numvec& v0, const numvec& v1, const numvec& v2, double angl, atomvec_t& set );
 
    void impose_angl( const morf_t& a0, const morf_t& a1, const morf_t& a2, double angl, atomvec_t& set );
 
    void impose_tors( const atomvec_t& atoms, double tors, bool middle_start );

    void impose_tors( const numvec& v0, const numvec& v1, const numvec& v2, const numvec& v3, double tors, atomvec_t& set );

    void impose_tors( const morf_t& a0, const morf_t& a1, const morf_t& a2, const morf_t& a3, double tors, atomvec_t& set );

    void impose_dist( atomvec_t& atoms, const vector< string >& inter );

    void impose_angl( atomvec_t& atoms, const vector< string >& inter );

    void impose_tors( atomvec_t& atoms, const vector< string >& inter );

} // namespace amber

#endif
