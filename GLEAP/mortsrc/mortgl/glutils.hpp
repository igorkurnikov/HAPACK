#ifndef MORT_GLUTIL_HPP
#define MORT_GLUTIL_HPP

#include <object.hpp>
#include <guilib.hpp>

namespace mort
{
    class graphic_t;

    /// molecular display style base class
    /// 
    namespace GL
    {
        void light();
        
        void anti_alias( bool flag );
        
        void depth_cue( bool flag );

        void new_list( const string& name );
    
        bool has_list( const string& name );

        void end_list( );
    
        void call_list( const string& name );

        void set_material( int matid );
    
    } // namespace GL
 

} // namespace mort

#endif

