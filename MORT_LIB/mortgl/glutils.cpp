#include <iomanip>
#include <GL/gl.h>
#include <GL/glu.h>
#include <object.hpp>
#include <common.hpp>
#include "glutils.hpp"
#include "graphic.hpp"

namespace mort
{
    namespace GL
    {        
        void light( )
        {        
            GLfloat position[] = { 0.0, 0.0, 500.0, 1.0 };
            GLfloat ambient[] = { 0.5, 0.5, 0.5, 1.0 };
            GLfloat diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
            GLfloat local_view = 0.0;

  
            glLightfv(GL_LIGHT0, GL_POSITION, position);
            glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse );
            glLightfv(GL_LIGHT0, GL_AMBIENT, ambient );

            glLightModelfv(GL_LIGHT_MODEL_LOCAL_VIEWER, &local_view);
    

            glFrontFace(GL_CW);
            glEnable(GL_LIGHTING);
            glEnable(GL_LIGHT0);
        
            glEnable(GL_AUTO_NORMAL);
            glEnable(GL_NORMALIZE);

        }    

        void depth_cue( bool flag )
        {
            if( flag )
            {
                GLfloat color[4] = { 0.1, 0.1, 0.1, 0.5 };

                glEnable( GL_FOG );            
                glEnable( GL_DEPTH_TEST );
                glHint( GL_FOG_HINT, GL_NICEST );
                glFogf( GL_FOG_MODE, GL_EXP );
                glFogf(GL_FOG_DENSITY, 0.0025f);
                glFogfv( GL_FOG_COLOR, color);
            }
            else
            {
                glDisable( GL_FOG );
            }
        }

        void anti_alias( bool flag )
        {
            if( flag )
            {
                glEnable(GL_BLEND);
                glEnable(GL_LINE_SMOOTH);
                glEnable(GL_POLYGON_SMOOTH);
                glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
                glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
            }
            else
            {
                glDisable( GL_BLEND );
                glDisable( GL_LINE_SMOOTH );
                glDisable( GL_POLYGON_SMOOTH );
            }
        }

        void new_list( const string& name )
        {
            glNewList( hash( name.c_str() ), GL_COMPILE );
        }
    
        void call_list( const string& name )
        {
            glCallList( hash( name.c_str() ) );
        }

        bool has_list( const string& name )
        {
            return glIsList( hash( name.c_str() ) );
        }
    
        void end_list( )
        {
            glEndList();
        }


        void set_material( int element )
        {
            const double ambient_ratio = 0.4, diffuse_ratio=0.8, specular_ratio = 0.2;
            double red   = (element==int(SELECTED)) ? 1.0 : pertab_t::get_red( element );
            double blue  = (element==int(SELECTED)) ? 0.0 : pertab_t::get_blue( element );
            double green = (element==int(SELECTED)) ? 1.0 : pertab_t::get_green( element );
            float shineness=50.0;
        
            float ambient[4], diffuse[4], specular[4];
        
            ambient[0] = ambient_ratio*red;
            ambient[1] = ambient_ratio*green;
            ambient[2] = ambient_ratio*blue; 
            ambient[3]=1.0;

            diffuse[0] = diffuse_ratio*red; 
            diffuse[1] = diffuse_ratio*green; 
            diffuse[2] = diffuse_ratio*blue; 
            diffuse[3]=1.0;
            
            specular[0]=specular_ratio*red; 
            specular[1]=specular_ratio*green; 
            specular[2]=specular_ratio*blue; 
            specular[3]=1.0;


            glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT, ambient);
            glMaterialfv(GL_FRONT_AND_BACK,GL_DIFFUSE, diffuse);
            glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR, specular);
            glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&shineness);
        }

    } // namespace GL



} // namespace mort

