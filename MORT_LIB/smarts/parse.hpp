#ifndef DDUTILS_MORT_DAYLIGHT_PARSE_H
#define DDUTILS_MORT_DAYLIGHT_PARSE_H

namespace mort
{
    class molecule_t;

    namespace daylight
    {
        
        class status_t;

        class reader_set;

        void parse( const char* ptr, molecule_t& mol, const reader_set& readers, int lan );

        const char* parse_env( const char* ptr, status_t& status, molecule_t& mol );

        const char* parse_bond( const char* ptr, status_t& status, molecule_t& mol );

        const char* parse_ring( const char* ptr, status_t& status, molecule_t& mol );

        const char* parse_alpha( const char* ptr, status_t& status, molecule_t& mol );

        const char* parse_brack( const char* ptr, status_t& status, molecule_t& mol );

        const char* parse_sharp( const char* ptr, status_t& status, molecule_t& mol );

        const char* parse_paren( const char* ptr, status_t& status, molecule_t& mol );

        const char* parse_logic( const char* ptr, status_t& status, molecule_t& mol );

        const char* parse_rinfo( const char* ptr, status_t& status, molecule_t& mol );
        
        const char* parse_charge( const char* ptr, status_t& status, molecule_t& mol );

        const char* parse_weight( const char* ptr, status_t& status, molecule_t& mol );

        const char* parse_asterisk( const char* ptr, status_t& status, molecule_t& mol );

        const char* parse_aromatic( const char* ptr, status_t& status, molecule_t& mol );
        
    } // namespace daylight
    
} // namespace mort

#endif
