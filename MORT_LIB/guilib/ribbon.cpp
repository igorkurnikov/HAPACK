#include <GL/gl.h>
#include <GL/glu.h>
#include "object.hpp"
#include "glutils.hpp"

namespace mort
{
    using std::find_if;

    static const int NSEGMENT = 6;

    void ctrl_points( const molecule_t& mol, vector< numvec >& ctrls_a, vector< numvec >& ctrls_b )
    {
        numvec prev_o;
        numvec curt_o;
        numvec curt_ca;
        numvec prev_ca;
        numvec prev_off;
        numvec curt_off;

        iterator_t resd = mol.resd_begin();
        for( ; resd != mol.resd_end(); ++resd )
        {
            iterator_t ca = find_if( resd->related_atom_begin(), resd->related_atom_end(), 
                                     sparm_comper1(NAME, "CA") );

            iterator_t o  = find_if( resd->related_atom_begin(), resd->related_atom_end(), 
                                     sparm_comper1(NAME, "O") );
            
            if( ca == resd->related_atom_end() || o == resd->related_atom_end() )
            {
                continue;
            } 

            curt_ca = ca->get_v(POSITION);
            curt_o = o->get_v(POSITION);
          
	    int head;
            if( resd->get_i(HEAD,head) && head == 0 )
            {
                ctrls_a.push_back( curt_ca );
                ctrls_b.push_back( curt_ca );          
                prev_o  = curt_o;
                prev_ca = curt_ca;
                prev_off = numvec( 0.0, 0.0, 0.0 );
            }
            else
            {
                numvec a = curt_ca - prev_ca;
                numvec b = prev_o  - prev_ca;
                numvec c = cross( a, b );
                numvec d = cross( c, a );

                if( inner_prod( d, prev_off ) < 0 )
                {
                    curt_off = normalize_copy( prev_off - d );
                }
                else
                {
                    curt_off = normalize_copy( prev_off + d );
                }

                double width = 0.1;

                iterator_t prev = resd - 1;
                iterator_t next = resd + 1;

                int second = resd->get_i(SECOND);
                if( resd == mol.resd_begin( ) || next == mol.resd_end() )
                {
                    width = 0.1;
                }
                else if( second!=LOOP && prev->get_i(SECOND)==second && next->get_i(SECOND)==second )
                {
                    width = 1.0;
                    std::cout << "resd " << resd->get_s(NAME) << " helix|sheet" << std::endl;
                }
                else
                {                    
                    width = 0.1;
                }

                ctrls_a.push_back( curt_ca + width*curt_off );
                ctrls_b.push_back( curt_ca - width*curt_off );

                prev_o = curt_o;
                prev_ca = curt_ca;                    
                prev_off = curt_off;
                
            }
        }
    }       

    void side_points( const vector< numvec >& ctrls, component_t* , vector< int >::iterator si )
    {
        //static const double inc =  1.0 / NSEGMENT;    
        
        vector< numvec >::const_iterator ci = ctrls.begin();

        for( ; ci != ctrls.end() - 3; ++ci )
        {
            for( int i=0; i < NSEGMENT; ++i, si += 2 )
            {
	    /*
                double t1 = i * inc;
                double t2 = t1 * t1;
                double t3 = t1 * t2;
                
                double s0 = ( 1 - t1 ) * ( 1 - t1 ) * ( 1 - t1 ) / 6.0;
                double s1 = ( 3 * t3 - 6 * t2 + 4 ) / 6.0;
                double s2 = (-3 * t3 + 3 * t2 + 3 * t1 + 1 ) / 6.0;
                double s3 = t3 / 6.0;
              */  
                // ribbon->at< position_t >( *si ) = *ci * s0 + *( ci + 1 ) * s1 + *( ci + 2 ) * s2 + *( ci + 3 ) * s3;
            }
        }    
    }
    numvec get_normal( const numvec& a, const numvec& b, const numvec& c )
    {
        numvec ab = a - b;
        numvec ac = c - a;   
        return normalize_copy( cross( ab, ac ) );
    }


/*
    void side_normals( component_t* ribbon, vector< int >::iterator si, vector< int >::iterator se )
    {
        numvec curt_1;
        numvec curt_2;
        numvec next_1 = get_normal( ribbon->at< position_t >( *( si + 2 ) ),
                                      ribbon->at< position_t >( *( si + 1 ) ),
                                      ribbon->at< position_t >( *( si ) ) );

        numvec next_2 = get_normal( ribbon->at< position_t >( *( si + 1 ) ),
                                      ribbon->at< position_t >( *( si + 2 ) ),
                                      ribbon->at< position_t >( *( si + 3 ) ) );
        
        ribbon->at< normal_t >( *( si + 0 ) ) = next_1;
        ribbon->at< normal_t >( *( si + 1 ) ) = normalize_copy( next_1 + next_2 );

        for( ; si < se - 4; si += 2 )
        {
            curt_1 = next_1;
            curt_2 = next_2;
            
            next_1 = get_normal( ribbon->at< position_t >( *( si + 4 ) ),
                                 ribbon->at< position_t >( *( si + 3 ) ),
                                 ribbon->at< position_t >( *( si + 2 ) ) );

            next_2 = get_normal( ribbon->at< position_t >( *( si + 3 ) ),
                                 ribbon->at< position_t >( *( si + 4 ) ),
                                 ribbon->at< position_t >( *( si + 5 ) ) );

            ribbon->at< normal_t >( *( si + 2 ) ) = normalize_copy( curt_1 + curt_2 + next_1 );
            ribbon->at< normal_t >( *( si + 3 ) ) = normalize_copy( curt_2 + next_1 + next_2 );
        }

//        *( ni + 4 ) = normalize( next_1 + next_2 );        
//        *( ni + 5 ) = next_2;
    }


    
    bool create_ribbon( molecule_t& mol )
    {
        vector< numvec > ctrls_a, ctrls_b;

        if( mol.has_component( RIBBON ) )
        {
            return false;
        }

        ctrl_points( mol, ctrls_a, ctrls_b );

        mol.add_component( RIBBON, 12 * ctrls_a.size() - 24 );
        
        component_t* ribbon = mol.get_component( RIBBON );
        
        side_points( ctrls_a, ribbon, ribbon->begin() );
        
        side_points( ctrls_b, ribbon, ribbon->begin() + 1 );       

        side_normals( ribbon, ribbon->begin(), ribbon->end() );
        
        return true;
    }
*/



    void glpoint( const numvec& posi, const numvec& norm )
    {
        glNormal3d( norm[ 0 ], norm[ 1 ], norm[ 2 ] );
        glVertex3d( posi[ 0 ], posi[ 1 ], posi[ 2 ] );
    }


    /* 
    void render_ribbon( const molecule_t& mol )
    {
        float ambient[] = { 0.0, 0.2, 0.0, 1.0 };
        float diffuse[] = { 0.0, 0.6, 0.0, 1.0 };
        float shiness[] = {30.0};
        
        glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );
        glEnable( GL_POLYGON_SMOOTH );
        
        
        glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT, ambient );
        glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE, diffuse );
        glMaterialfv( GL_FRONT_AND_BACK, GL_SHININESS, shiness );

        BOOST_ASSERT( mol.has_component( RIBBON ) );

        const component_t* ribbon = mol.get_component( RIBBON );

        vector< int >::iterator pi = ribbon->begin();
        vector< int >::iterator pe = ribbon->end();

        for( ; pi < pe - 2; pi += 2 )
        {
            glBegin( GL_TRIANGLES );
            {
                glpoint( ribbon->at< position_t >( *( pi + 2 ) ), ribbon->at< normal_t >( *( pi + 2 ) ) );
                glpoint( ribbon->at< position_t >( *( pi + 1 ) ), ribbon->at< normal_t >( *( pi + 1 ) ) );
                glpoint( ribbon->at< position_t >( *( pi ) ),     ribbon->at< normal_t >( *( pi ) ) );
                glpoint( ribbon->at< position_t >( *( pi + 1 ) ), ribbon->at< normal_t >( *( pi + 1 ) ) );
                glpoint( ribbon->at< position_t >( *( pi + 2 ) ), ribbon->at< normal_t >( *( pi + 2 ) ) );
                glpoint( ribbon->at< position_t >( *( pi + 3 ) ), ribbon->at< normal_t >( *( pi + 3 ) ) );            
            } 
            glEnd();
        }
    }

    

    
    void render_ribbon( const molecule_t& mol )
    {
        GLfloat ctrls[][3] = { { -4.0, -8.0, 0.0 }, { -4.0,  8.0, 0.0 }, {  4.0,  8.0, 0.0 }, {  4.0, -8.0, 0.0 }, 
                               { -4.0, -8.0, 0.0 }, { -4.0,  8.0, 0.0 }, {  4.0,  8.0, 0.0 }, {  4.0, -8.0, 0.0 }, 
                               { -6.0, -12.0, -2.0 }, { -6.0,  12.0, -2.0 }, {  6.0,  12.0, -2.0 }, {  6.0, -12.0, -2.0 }, 
                               { -6.0, -12.0, -2.0 }, { -6.0,  12.0, -2.0 }, {  6.0,  12.0, -2.0 }, {  6.0, -12.0, -2.0 }, 
                               { -6.0, -12.0, -4.0 }, { -6.0,  12.0, -4.0 }, {  6.0,  12.0, -4.0 }, {  6.0, -12.0, -4.0 }, 
                               { -6.0, -12.0, -4.0 }, { -6.0,  12.0, -4.0 }, {  6.0,  12.0, -4.0 }, {  6.0, -12.0, -4.0 }, 
                               { -4.0, -8.0, -6.0 }, { -4.0,  8.0, -6.0 }, {  4.0,  8.0, -6.0 }, {  4.0, -8.0, -6.0 }, 
                               { -4.0, -8.0, -6.0 }, { -4.0,  8.0, -6.0 }, {  4.0,  8.0, -6.0 }, {  4.0, -8.0, -6.0 } };
        


        GLfloat knots_1[] = { 0.0, 1.0/7, 2.0/7, 3.0/7, 4.0/7, 5.0/7, 6.0/7, 1.0 };

        GLUnurbsObj* nurb = gluNewNurbsRenderer();

        for( int i=0; i < 4; ++i )
        {
            gluBeginSurface( nurb );

                gluNurbsSurface( 
                    nurb, 
                    8, knots_1,
                    8, knots_1,
                    3, 8*3,
                    &ctrls[i][0], 
                    4, 4,
                    GL_MAP2_VERTEX_3);

            gluEndSurface( nurb );
        }
                
    }                              

    */

    void render_ribbon( vector< GLfloat >& ctrls, vector< int >& mats )
    {
        int npoint = ctrls.size() / 3;

        if( npoint < 32 )
        {
            return;
        }

        GLfloat knots_1[] = { 0.0, 1.0/7, 2.0/7, 3.0/7, 4.0/7, 5.0/7, 6.0/7, 1.0 };
        //GLfloat knots_2[] = { 0.0,   0.0,   0.0,   0.0,   1.0,   1.0,   1.0, 1.0 };


        GLUnurbsObj* nurb = gluNewNurbsRenderer();

 	for( int i=0; i < npoint - 28; i += 8 )
        {
            set_material( mats[i/8] );

            for( int j=0; j < 4; j++ )
            {
                gluBeginSurface( nurb );

                    gluNurbsSurface( 
                        nurb, 
                        8, knots_1,
                        8, knots_1,
                        3,8*3,
                        &ctrls[3*i+3*j], 
                        4,4,
                        GL_MAP2_VERTEX_3);

                gluEndSurface( nurb );
            }
        }
        
        gluDeleteNurbsRenderer( nurb );
    }

   
    void render_ribbon( const molecule_t& mol )
    {
        vector< numvec >  ctrls_a, ctrls_b;
        ctrl_points( mol, ctrls_a, ctrls_b );

        vector< int > maters;
        vector< GLfloat > glctrls;

        
        for( int i=0; i < (int)ctrls_a.size(); ++i )
        {
            pointer_t resd( mol, ATOM, i );

            numvec normal_a1, normal_a2, normal_a3;
            numvec normal_b1, normal_b2, normal_b3;

            int head;
            if( resd.get_i(HEAD, head) && head == 0 )
            {
                normal_a1 = numvec( 0.0, 0.0, 0.0 );
                normal_b1 = numvec( 0.0, 0.0, 0.0 );
                normal_b2 = numvec( 0.0, 0.0, 0.0 );
            }
            else
            {
                normal_a1 = get_normal( ctrls_a[i], ctrls_b[i], ctrls_a[i-1] );
                normal_b1 = get_normal( ctrls_b[i], ctrls_b[i-1], ctrls_a[i-1] );
                normal_b2 = get_normal( ctrls_b[i], ctrls_a[i-1], ctrls_a[i] );
            }

            int tail;
            if( resd.get_i(TAIL, tail) && head == 0 )
            {
                normal_a2 = numvec( 0.0, 0.0, 0.0 );
                normal_a3 = numvec( 0.0, 0.0, 0.0 );
                normal_b3 = numvec( 0.0, 0.0, 0.0 );
            }
            else
            {
                normal_a2 = get_normal( ctrls_a[i], ctrls_b[i+1], ctrls_b[i] );
                normal_a3 = get_normal( ctrls_a[i], ctrls_a[i+1], ctrls_b[i+1] );
		normal_b3 = get_normal( ctrls_b[i], ctrls_a[i], ctrls_b[i+1] );
           }

            numvec normal_a = normalize_copy( normal_a1 + normal_a2 + normal_a3 );
            numvec normal_b = normalize_copy( normal_b1 + normal_b2 + normal_b3 );

            numvec a1 = ctrls_a[i] + 0.1*normal_a;
            numvec b1 = ctrls_b[i] + 0.1*normal_b;
            numvec a2 = ctrls_a[i] - 0.1*normal_a;
            numvec b2 = ctrls_b[i] - 0.1*normal_b;


            glctrls.push_back( GLfloat( a2[0] ) );
            glctrls.push_back( GLfloat( a2[1] ) );
            glctrls.push_back( GLfloat( a2[2] ) );

            glctrls.push_back( GLfloat( a1[0] ) );
            glctrls.push_back( GLfloat( a1[1] ) );
            glctrls.push_back( GLfloat( a1[2] ) );

            glctrls.push_back( GLfloat( b1[0] ) );
            glctrls.push_back( GLfloat( b1[1] ) );
            glctrls.push_back( GLfloat( b1[2] ) );

            glctrls.push_back( GLfloat( b2[0] ) );
            glctrls.push_back( GLfloat( b2[1] ) );
            glctrls.push_back( GLfloat( b2[2] ) );

            glctrls.push_back( GLfloat( a2[0] ) );
            glctrls.push_back( GLfloat( a2[1] ) );
            glctrls.push_back( GLfloat( a2[2] ) );

            glctrls.push_back( GLfloat( a1[0] ) );
            glctrls.push_back( GLfloat( a1[1] ) );
            glctrls.push_back( GLfloat( a1[2] ) );

            glctrls.push_back( GLfloat( b1[0] ) );
            glctrls.push_back( GLfloat( b1[1] ) );
            glctrls.push_back( GLfloat( b1[2] ) );


            glctrls.push_back( GLfloat( b2[0] ) );
            glctrls.push_back( GLfloat( b2[1] ) );
            glctrls.push_back( GLfloat( b2[2] ) );

            if( i+1 == (int)ctrls_a.size() )
            {
                maters.push_back( 6 );
            }
            else
            {
                pointer_t next(mol, ATOM, i+1);

                int second = resd.get_i(SECOND);
                
                if( second != LOOP && next.get_i(SECOND) == second )
                {
                    assert( second == SHEET || second == HELIX );

                    maters.push_back( (second == HELIX? 8 : 16 ) );
                }
                else 
                {
                    maters.push_back( 6 );
                }
            }
            
            if( resd.get_i(TAIL,tail) && tail == 0 )
            {
                render_ribbon( glctrls, maters );
                glctrls.clear();
                maters.clear();
            }
        }

	if( glctrls.size() != 0 )
        {
            render_ribbon( glctrls, maters );
        }

    } 



} // namespace mort
