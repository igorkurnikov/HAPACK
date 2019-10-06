#ifndef DDUTILS_MORT_DAYLIGHT_READER_H
#define DDUTILS_MORT_DAYLIGHT_READER_H
#include <boost/function.hpp>
#include <common.hpp>

namespace mort
{
    class molecule_t;

    namespace daylight
    {
        struct status_t;
        
        typedef boost::function< const char* ( const char*, status_t&, molecule_t& ) > reader_t;

        class reader_set
        {
        public:

            reader_set( int lan );

            virtual ~reader_set();

            void add( char c, int brack, reader_t reader );

            void add( const char* ptr, int brack, reader_t reader );
        
            reader_t get( char c, int bracket ) const;
        
        private:

            reader_t m_readers[256];
        };    

        void init_smiles( reader_set& readers );
        
        void init_smarts( reader_set& readers );

    } // namespace smarts
    
} // namespace mort

#endif
