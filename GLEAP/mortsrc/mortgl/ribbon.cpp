#include <GL/gl.h>
#include <GL/glu.h>
#include "object.hpp"
#include "glutils.hpp"

namespace mort
{
    using std::find_if;

    static const int NSEGMENT = 6;

    void ctrl_points( const molecule_t& mol, vector<numvec>& ctrls_a, vector<numvec>& ctrls_b )
    {
        numvec prev_o(3);
        numvec curt_o(3);
        numvec curt_ca(3);
        numvec prev_ca(3);
        numvec prev_off(3);
        numvec curt_off(3);

        moiter_t resd = mol.resd_begin();
        for( ; resd != mol.resd_end(); ++resd )
        {
            moiter_t ca = find_if( resd->related_atom_begin(), resd->related_atom_end(), 
                                     sparm_comper1(NAME, "CA") );

            moiter_t o  = find_if( resd->related_atom_begin(), resd->related_atom_end(), 
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

                moiter_t prev = resd - 1;
                moiter_t next = resd + 1;

                std::cout << "resd " << resd->get_s(NAME) << " " << "head, tail: " << resd->get_i(HEAD) << " " << resd->get_i(TAIL) << std::endl;

                int second = resd->get_i(SECOND);
                if( resd == mol.resd_begin( ) || next == mol.resd_end() )
                {
                    width = 0.1;
                }
                else if( second!=LOOP )
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

    numvec get_normal( const numvec& a, const numvec& b, const numvec& c )
    {
        numvec ab = a - b;
        numvec ac = c - a;   
        return normalize_copy( cross(ab, ac) );
    }

    void glpoint( const numvec& posi, const numvec& norm )
    {
        glNormal3d( norm[0], norm[1], norm[2] );
        glVertex3d( posi[0], posi[1], posi[2] );
    }

    void render_ribbon( vector< GLfloat >& ctrls, vector< int >& mats )
    {
        int npoint = ctrls.size() / 3;
        if( npoint < 32 )
        {
            return;
        }

        //GLfloat knots_1[] = { 0.0, 1.0/15,  2.0/15,  3.0/15,  4.0/15,  5.0/15,  6.0/15.0, 7.0/15.0,
	//                 8.0/15.0, 9.0/15, 10.0/15, 11.0/15, 12.0/15, 13.0/15, 14.0/15.0,   1.0};
        //GLfloat knots_2[] = { 0.0,   0.0,   0.0,   0.0,   1.0,   1.0,   1.0, 1.0 };
	GLfloat knots_1[] = { 0.0, 1.0/3, 2.0/3, 1.0 };
        GLfloat knots_2[] = { 0.0, 1.0/7, 2.0/7, 3.0/7, 4.0/7, 5.0/7, 6.0/7, 1.0 };


        GLUnurbsObj* nurb = gluNewNurbsRenderer();

 	for( int i=0; i < npoint - 28; i += 8 )
        {
            GL::set_material( mats[i/8] );

            for( int j=0; j < 4; j++ )
            {
                gluBeginSurface( nurb );

                    gluNurbsSurface( 
                        nurb, 
                        8, knots_2,
                        8, knots_2,
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
	std::cout << "end of ctrl_points" << std::endl;

        vector< int > maters;
        vector< GLfloat > glctrls;

        
        for( int i=0; i < (int)ctrls_a.size(); ++i )
        {
	    std::cout << "ctrl " << i << std::endl;

            moref_t resd( mol, RESD, i );

            numvec normal_a1(3), normal_a2(3), normal_a3(3);
            numvec normal_b1(3), normal_b2(3), normal_b3(3);

            int head;
            if( resd.get_i(HEAD, head) && head == 0 )
            {
	        std::cout << "head init" << std::endl;
                normal_a1 = numvec( 0.0, 0.0, 0.0 );
                normal_b1 = numvec( 0.0, 0.0, 0.0 );
                normal_b2 = numvec( 0.0, 0.0, 0.0 );
            }
            else
            {
	        assert( i > 0 );
	        std::cout << "head " << head << " get normal" << std::endl;
                normal_a1 = get_normal( ctrls_a[i], ctrls_b[i], ctrls_a[i-1] );
                normal_b1 = get_normal( ctrls_b[i], ctrls_b[i-1], ctrls_a[i-1] );
                normal_b2 = get_normal( ctrls_b[i], ctrls_a[i-1], ctrls_a[i] );
            }

	    std::cout << "after head" << std::endl;

            int tail;
            if( resd.get_i(TAIL, tail) && tail==0 )
            {
                normal_a2 = numvec( 0.0, 0.0, 0.0 );
                normal_a3 = numvec( 0.0, 0.0, 0.0 );
                normal_b3 = numvec( 0.0, 0.0, 0.0 );
            }
            else
            {
                assert( i+1 < ctrls_b.size() );
                normal_a2 = get_normal( ctrls_a[i], ctrls_b[i+1], ctrls_b[i] );
                normal_a3 = get_normal( ctrls_a[i], ctrls_a[i+1], ctrls_b[i+1] );
		normal_b3 = get_normal( ctrls_b[i], ctrls_a[i], ctrls_b[i+1] );
            }

	    std::cout << "after tail" << std::endl;

            numvec normal_a = normalize_copy( normal_a1 + normal_a2 + normal_a3 );
            numvec normal_b = normalize_copy( normal_b1 + normal_b2 + normal_b3 );

            numvec a1 = ctrls_a[i] + 0.1*normal_a;
            numvec b1 = ctrls_b[i] + 0.1*normal_b;
            numvec a2 = ctrls_a[i] - 0.1*normal_a;
            numvec b2 = ctrls_b[i] - 0.1*normal_b;

            std::cout << "after norm" << std::endl;

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

            if( i==0 || i+1 == (int)ctrls_a.size() )
            {
                maters.push_back( 6 );
            }
            else
            {
	        moref_t prev(mol, RESD, i-1);
                moref_t next(mol, RESD, i+1);

                int prev_second = prev.get_i(SECOND);
                int curt_second = resd.get_i(SECOND);
                
                if( curt_second!=LOOP )
                {
                    maters.push_back( (curt_second==HELIX? 8 : 16) );
                }
		else if( prev_second!=LOOP )
		{
		    maters.push_back( (prev_second==HELIX? 8 : 16) );
		}
                else 
                {
                    maters.push_back( 6 );
                }
            }
            
            if( resd.get_i(TAIL,tail) && tail == 0 )
            {

	        std::cout << "render ribbon" << std::endl << std::endl;
                render_ribbon( glctrls, maters );
                glctrls.clear();
                maters.clear();
            }
        }

	std::cout << "end glctrls" << std::endl;

	if( glctrls.size() != 0 )
        {
            render_ribbon( glctrls, maters );
        }

    } 



} // namespace mort


