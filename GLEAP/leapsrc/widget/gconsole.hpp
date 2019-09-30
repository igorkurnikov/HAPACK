#ifndef GLEAP_LEAPSRC_WIDGET_GCONSOLE_HPP
#define GLEAP_LEAPSRC_WIDGET_GCONSOLE_HPP

#include <common.hpp>
#include <guilib.hpp>


namespace mortgtk
{
    using namespace mort;
    
    class gtk_console : public console_t
    {
    public:

        gtk_console();
        
        virtual ~gtk_console();
        
        virtual void print( const std::string& str );
        
        virtual char getchar();
        
        virtual bool getline( string& line, const char* prompt=NULL );
    };
    
} // namespace mortgtk



#endif
