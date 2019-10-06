#include <stdio.h>
#include <iostream>
#include <boost/format.hpp>
#include "geometry.hpp"

#if !defined(M_PI)
#define M_PI 3.141592653589793238462643
#endif 

namespace mort
{
    using std::cout;

    double dis2( const numvec& v1, const numvec& v2 )
    {
        numvec v12 = v1 - v2;
        return dotpd( v12, v12 );
    }
 
    double dist( const numvec& v1, const numvec& v2 )
    {
        return sqrt( dis2(v1, v2) );
    }

    double angl( const numvec& v1, const numvec& v2, const numvec& v3 )
    {
        numvec v12 = normalcpy( v1 - v2 );
        numvec v32 = normalcpy( v3 - v2 );
        return acos(dotpd(v12, v32))*180.0/M_PI;
    }

    double tors( const numvec& v1, const numvec& v2, const numvec& v3, const numvec& v4 )
    {
        numvec v12 = v1 - v2;
        numvec v32 = v3 - v2;
        numvec v42 = v4 - v2;
        numvec ax1 = normalcpy( cross(v12, v32) );
        numvec ax2 = normalcpy( cross(v42, v32) );
        numvec ax3 = cross(ax1, ax2);

        double cot = dotpd(ax1, ax2);

        if( cot >= 1.0 )
        {
            return 0;
        }
        else if( cot <= -1.0 )
        {
            return 180.0;
        }
        
        double theta = acos(cot) * 180.0/ M_PI;

        return dotpd(ax3, v32) > 0 ? theta : 360.0 -theta;
    }

    numvec& rotate( numvec& src, const numvec& axs, double ang )
    {
        numvec dst = rotcpy(src, axs, ang);
        src = dst;
        return src;
    }

    double* rotate( double* vec, const matrix& mat )
    {
        assert( mat.size1()==3 && mat.size2()==3 );
        double x = vec[0];
        double y = vec[1];
        double z = vec[2];
        vec[0] = mat(0,0)*x + mat(1,0)*y + mat(2,0)*z;
        vec[1] = mat(0,1)*x + mat(1,1)*y + mat(2,1)*z;
        vec[2] = mat(0,2)*x + mat(1,2)*y + mat(2,2)*z;
        return vec;
    }

    numvec rotcpy( const numvec& src, numvec axs, double ang )
    {
        double rad = ang * M_PI / 180.0;
        normalize( axs );
        numvec dst = src*cos(rad) + axs * dotpd(axs, src) * (1 - cos(rad)) - cross(src, axs) * sin(rad);
        return dst;
    }

    void printvec( const std::vector<numvec>& vecs )
    {
		char buf[80];
        for(unsigned int i=0; i < vecs.size(); ++i )
        {
			sprintf(buf,"vec %2d,%d: %10.3f %10.3f %10.3f",i,vecs[i].size(),vecs[i][0],vecs[i][1],vecs[i][2]); 
            std::cout << buf << std::endl;
        }

        for( unsigned int i=0; i < vecs.size(); ++i )
        {
            for( unsigned int j=i+1; j < vecs.size(); ++j)
            {
                std::cout << "dist ";
				sprintf(buf,"%2d %2d: %10.3f ",i,j,dist( vecs[i], vecs[j] ) );
                std::cout << buf << std::endl;
            }
        }
    }
 
    numvec impose( std::vector< numvec > srcpos, const std::vector< numvec >& dstpos, const numvec& input )
    {
        assert( srcpos.size() == dstpos.size() );
        assert( srcpos.size() > 0 );
        assert( srcpos.size() <= 3 );

        srcpos.push_back( input );

        //
        //  std::cout << "before impose:" << std::endl;
        //  printvecs( srcpos );

        assert( srcpos.size() >= 2 );
        
        numvec dist_0 = dstpos[0] - srcpos[0];
        
        for( int i=0; i < (int)srcpos.size(); i++ )
        {
            srcpos[i] += dist_0;
        }

        if( srcpos.size() > 2 )
        {
            numvec src = srcpos[1] - srcpos[0];
            numvec dst = dstpos[1] - dstpos[0];
            numvec axs = cross( src, dst );
            double n = norm( axs );

            if( n > 1e-5 )
            {
                axs /= n;
                double ang = acos( dotpd( normalcpy(src), normalcpy(dst) ) ) * 180.0  / M_PI;
                for( int i=1; i < (int)srcpos.size(); i++ )
                {
                    srcpos[i] = rotcpy( srcpos[i]-srcpos[0], axs, ang) + dstpos[0];
                }
            }
        }
            

        if( srcpos.size() > 3 )
        {
            numvec face_0 = cross( srcpos[2] - srcpos[0], srcpos[1] - srcpos[0] );
            numvec face_1 = cross( dstpos[2] - dstpos[0], dstpos[1] - dstpos[0] );
            normalize(face_0);
	    normalize(face_1);
            numvec axs = cross(face_0, face_1); 

            double n = norm( axs );
            if( n > 1e-5 )
            { 
                axs /= n;
                double ang = acos( dotpd(face_0, face_1) ) * 180.0  / M_PI;

                for( int i=2; i < (int)srcpos.size(); i++ )
                {
                    srcpos[i] = rotcpy( srcpos[i] - srcpos[0], axs, ang ) + dstpos[0];
                }
            }
        }


       
        //
        // std::cout << "after impose:" << std::endl;
        // printvecs( srcpos );
        //
	
        return srcpos.back();
    }

    numvec plane( const numvec& p0, const numvec& p1, const numvec& p2 )
    {
        numvec a = p1 - p0;
        numvec b = p2 - p0;
        numvec n = cross( a, b );
        double d = dotpd( normalize(n), p0 );
        return makevec( n[0], n[1], n[2], d );
    }

    bool same_line( const std::vector<numvec>& pts )
    {
        assert( pts.size() >= 3 );

        numvec line1 = pts[1] - pts[0];
        double norm1 = norm( line1 );

        for( int i=2; i < (int)pts.size(); ++i )
        {
            numvec line2 = pts[i] - pts[0];
            double norm2 = norm( line2 );

            double d = dotpd(line1, line2) / norm1 / norm2;

            if( std::abs(d) < 0.95 || std::abs(d) > 1.05 )
            {
                return false;
            }
        }

        return true;
    }

    bool same_plane( const numvec& c, const numvec& p )
    {
        assert( c.size()==3 && p.size()==4 );
        double d = p[0]*c[0] + p[1]*c[1] + p[2]*c[2];
        double r = std::abs(d-p[3]);
        return r < 0.2;
    }

    bool same_plane( const std::vector<numvec>& crds )
    {
        numvec dummy(4);
        return same_plane( crds, dummy );
    }

    bool same_plane( const std::vector<numvec>& crds, numvec& rp )
    {
        assert( crds.size() >= 4 );

        numvec p = plane( crds[0], crds[1], crds[2] );
        for( int i=3; i < (int)crds.size(); ++i )
        {
            if( ! same_plane(crds[i], p) )
            {
                return false;
            }
        }

        rp = p;
        return true;
    }    
 
    numvec& rotate_axis(numvec& src, int cidx, int sidx, double ang)
    {
        double arc = ang*M_PI/180.0;
        double c = src[cidx];
        double s = src[sidx];
        double newc = c*cos(arc) - s*sin(arc);
        double news = s*cos(arc) + c*cos(arc);
        src[cidx] = newc;
        src[sidx] = news;
        return src;
    }

    numvec& rotate_x(numvec& src, double ang)
    {
        return rotate_axis(src, 1, 2, ang);
    }
    
    numvec& rotate_y(numvec& src, double ang)
    {
        return rotate_axis(src, 2, 0, ang);
    }
    
    numvec& rotate_z(numvec& src, double ang)
    {
        return rotate_axis(src, 0, 1, ang);
    }

    numvec zmatrix(const numvec& p1, double dist)
    {
        return makevec( p1[0]+dist, p1[1], p1[2] );
    }
    
    numvec zmatrix(const numvec& p1, double dist, const numvec& p2, double angl)
    {
        numvec p21 = p2 - p1;
        double yang = atan2(p21[2], p21[0])*180.0/M_PI;
        rotate_y(p21, yang);
        // std::cout << "yang,p21:" << yang << " " << p21[0] << " " << p21[1] << " " << p21[2] << std::endl;
        
        double zang = atan2(p21[1], p21[0])*180.0/M_PI;
        double rad = angl*M_PI/180.0;
        numvec p31 = makevec( dist*cos(rad+zang), dist*sin(rad+zang), 0.0);

        rotate_y(p31, -yang);        
        return p31+p1;
    }

    numvec zmatrix(const numvec& p1, double dist, const numvec& p2, double angl, const numvec& p3, double tors)
    {
        numvec p32 = p3 - p2;
        numvec p12 = p1 - p2;
        rotate(p32, p12, tors);
        double curt_angl = acos( dotpd(p32, p12)/norm(p32)/norm(p12) )/M_PI*180.0;
        
        numvec axs = cross(p32, p12);
        rotate(p32, axs, angl+curt_angl-180.0);
  
        numvec p41 = p32/norm(p32)*dist;
        return p41 + p1;
    }

    void transform( double* v1, const double* mat, bool invert )
    {
        if( invert )
        {
            double x = mat[15]*v1[0] - mat[12];
            double y = mat[15]*v1[1] - mat[13];
            double z = mat[15]*v1[2] - mat[14];
            v1[0] = x*mat[0] + y*mat[1] + z*mat[2];
            v1[1] = x*mat[4] + y*mat[5] + z*mat[6];
            v1[2] = x*mat[8] + y*mat[9] + z*mat[10];
        }
        else
        {
            double x = v1[0];
            double y = v1[1];
            double z = v1[2];
            v1[0] = (x*mat[0] + y*mat[4] + z*mat[8] + mat[12])/mat[15];
            v1[1] = (x*mat[1] + y*mat[5] + z*mat[9] + mat[13])/mat[15];
            v1[2] = (x*mat[2] + y*mat[6] + z*mat[10]+ mat[14])/mat[15];
        }
    }

    void translate( double* v, const numvec& sft )
    {
        v[0] += sft[0];
        v[1] += sft[1];
        v[2] += sft[2];
    }

} // namespace mort

