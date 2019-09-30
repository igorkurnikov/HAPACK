#include <iostream>
#include <cassert>
#include "numvec.hpp"
#include "matrix.hpp"
#include "crdalgo.hpp"
#include "constant.hpp"
#include "funstack.hpp"
#include "geometry.hpp"

namespace mort
{
    using std::slice;

    numvec center( const vector<double>& crd )
    {
        numvec r(zero_numvec(3));

        int natom = crd.size()/3;
        for( int i=0; i < natom; ++i )
        {
            r[0] += crd[3*i  ];
            r[1] += crd[3*i+1];
            r[2] += crd[3*i+2];
        }

        r /= double(natom);
        return r;
    }

    numvec center( const vector<double>& crd, const vector<int>& ids )
    {
        numvec r(zero_numvec(6));

        for( int i=0; i < (int)ids.size(); ++i )
        {
            int id = ids[i];
            r[0] += crd[3*id  ];
            r[1] += crd[3*id+1];
            r[2] += crd[3*id+2];
        }

        r /= double(ids.size());
        return r;
    }


    // r[0:3] store min
    // r[3:6] store max
    numvec minmax( const vector<double>& crd )
    {
        assert( crd.size() > 0 );

        numvec r(zero_numvec(6));

        r[0] = r[3] = crd[0];
        r[1] = r[4] = crd[1];
        r[2] = r[5] = crd[2];

        int natom = crd.size()/3;
        for( int i=1; i < natom; ++i )
        {
            r[0] = std::min( r[0], crd[3*i  ] );
            r[1] = std::min( r[1], crd[3*i+1] );
            r[2] = std::min( r[2], crd[3*i+2] );
            r[3] = std::max( r[3], crd[3*i  ] );
            r[4] = std::max( r[4], crd[3*i+1] );
            r[5] = std::max( r[5], crd[3*i+2] );
        }

        return r;
    }
    
    numvec minmax( const vector<double>& crd, const vector<double>& vdwr )
    {
        assert( crd.size() > 0 );

        numvec r(zero_numvec(6));

        r[0] = DOUBLEMAX; r[3] = -DOUBLEMAX;
        r[1] = DOUBLEMAX; r[4] = -DOUBLEMAX;
        r[2] = DOUBLEMAX; r[5] = -DOUBLEMAX;

        int natom = crd.size()/3;
        for( int i=0; i < natom; ++i )
        {
            r[0] = std::min( r[0], crd[3*i  ]-vdwr[i] );
            r[1] = std::min( r[1], crd[3*i+1]-vdwr[i] );
            r[2] = std::min( r[2], crd[3*i+2]-vdwr[i] );
            r[3] = std::max( r[3], crd[3*i  ]+vdwr[i] );
            r[4] = std::max( r[4], crd[3*i+1]+vdwr[i] );
            r[5] = std::max( r[5], crd[3*i+2]+vdwr[i] );
        }

        return r;
    }

    numvec region( const vector<double>& crd )
    {
        numvec mmx = minmax( crd );
        numvec rgn(zero_numvec(6));

        rgn[0] = (mmx[0] + mmx[3])*0.5;
        rgn[1] = (mmx[1] + mmx[4])*0.5;
        rgn[2] = (mmx[2] + mmx[5])*0.5;
        rgn[3] = (mmx[3] - mmx[0]);
        rgn[4] = (mmx[4] - mmx[1]);
        rgn[5] = (mmx[5] - mmx[2]);

        assert( mmx.size()==6u );
        return rgn;
    }

    numvec region( const vector<double>& crd, const vector<double>& vdwr )
    {
        numvec mmx = minmax( crd, vdwr );
        numvec rgn(zero_numvec(6));
        rgn[0] = (mmx[0] + mmx[3])*0.5;
        rgn[1] = (mmx[1] + mmx[4])*0.5;
        rgn[2] = (mmx[2] + mmx[5])*0.5;
        rgn[3] = (mmx[3] - mmx[0]);
        rgn[4] = (mmx[4] - mmx[1]);
        rgn[5] = (mmx[5] - mmx[2]);

        assert( mmx.size()==6u && rgn.size()==6u );
        return rgn;
    }

    matrix moment( const vector<double>& crd )
    {
        int natom = crd.size()/3;

        double sxx = 0.0;
	double syy = 0.0;
        double szz = 0.0;
        double sxy = 0.0;
	double syz = 0.0;
	double szx = 0.0;

	for( int i=0; i < natom; ++i )
	{
	    double x = crd[3*i  ];
	    double y = crd[3*i+1];
	    double z = crd[3*i+2];

	    sxx += y*y + z*z;
	    syy += z*z + x*x;
	    szz += x*x + y*y;
	    sxy += x*y;
	    syz += y*z;
	    szx += z*x;
	}

        matrix mat(3,3);
	mat(0,0) =  sxx; mat(0,1) = -sxy; mat(0,2) = -szx;
	mat(1,0) = -sxy; mat(1,1) =  syy; mat(1,2) = -syz;
	mat(2,0) = -szx; mat(2,1) = -syz; mat(2,2) =  szz;
	return mat;
    }


    double octbufsize( const vector<double>& crd, double uinput )
    {
        int natom = crd.size()/3;

        numvec mmx = minmax( crd );

        double xhalf = 0.5*( mmx[3] - mmx[0] ) + uinput;
	double yhalf = 0.5*( mmx[4] - mmx[1] ) + uinput;
        double zhalf = 0.5*( mmx[5] - mmx[2] ) + uinput;

        double x = yhalf * zhalf;
	double y = zhalf * xhalf;
	double z = xhalf * yhalf;

	double tmp = 1.0 / sqrt(x*x + y*y + z*z);

        double xunit = x * tmp;
	double yunit = y * tmp;
	double zunit = z * tmp;

        double max = 0.0;
	for(int i=0; i < natom; ++i )
	{
	    tmp = std::abs(crd[3*i])*xunit  + std::abs(crd[3*i+1]) * yunit + std::abs(crd[3*i+2]) * zunit;
            if( max < tmp )
	        max = tmp;
	}

	double bmax = 0.5 * sqrt( xhalf*xhalf + yhalf*yhalf + zhalf*zhalf );
	tmp = max + uinput;

	if( tmp < bmax )
	{
	    return uinput;
	}

        return uinput * tmp / bmax;
    }

    void rotate( vector<double>& crd, const matrix& mat )
    {
        int natom = crd.size()/3;
        for( int i=0; i < natom; ++i )
        {
            rotate( &crd[3*i], mat );
        }
    }

    void rotate( vector<double>& crd, const matrix& mat, const vector<int>& ids )
    {
        for( int i=0; i < (int)ids.size(); ++i )
        {
            int id = ids.size();
            rotate( &crd[3*id], mat );
        }
    }

    void transform( vector<double>& crd, const double* mat )
    {
        int natom = crd.size()/3;
        for( int i=0; i < natom; ++i )
        {
            transform( &crd[3*i], mat );
        }
    }

    void transform( vector<double>& crd, const double* mat, const vector<int>& ids )
    {
        for( int i=0; i < (int)ids.size(); ++i )
        {
            int id = ids[i];
            transform( &crd[3*id], mat );
        }
    }

    void translate( vector<double>& crd, const numvec& sft )
    {
        int natom = crd.size()/3;
        for( int i=0; i < natom; ++i )
        {
            crd[3*i  ] += sft[0];
            crd[3*i+1] += sft[1];
            crd[3*i+2] += sft[2];
        }
    }

    void translate( vector<double>& crd, const numvec& sft, const vector<int>& ids )
    {
        for( int i=0; i < (int)ids.size(); ++i )
        {
            int id = ids[i];
            crd[3*id  ] += sft[0];
            crd[3*id+1] += sft[1];
            crd[3*id+2] += sft[2];
        }
    }

    void centralize( vector<double>& crd )
    {
        numvec cnt = center( crd );
        translate( crd, -cnt );
    }

    numvec regionlize( vector<double>& crd )
    {
        funstack_t::push( "regionlize" );
        numvec rgn = region( crd );
        translate( crd, -subvec(rgn, 0, 3));
        funstack_t::pop();
        return subvec(rgn, 3, 3);
    }

    numvec regionlize( vector<double>& crd, const vector<double>& vdwr )
    {
        numvec rgn = region( crd, vdwr );
        assert( rgn.size()==6u );
        translate( crd, -subvec(rgn, 0, 3));
        return subvec(rgn, 3, 3);
    }

    void rotatewald( vector<double>& crd )
    {
        double pi=3.1415927;
        double phi=pi/4.;
        double cos1=cos(phi);
	double sin1=sin(phi);
	double cos2=sqrt(2.)/sqrt(3.);
	double sin2=1./sqrt(3.);

        matrix mat(3,3);

        mat(0,0) = cos2*cos1; mat(1,0) = -cos2*sin1; mat(2,0) = -sin2;
	mat(0,1) =-sin2*cos1; mat(1,1) =  sin2*sin1; mat(2,1) = -cos2;
	mat(0,2) = sin1;      mat(1,2) =  cos1;      mat(2,2) = 0.0;

        rotate( crd, mat );
    }

    // rotate to the best orientation, that longest distance is along x,y,z axis
    void rotatelong( vector<double>& crd )
    {
        matrix mom = moment(crd );
	numvec eigval(3);  // eigen values
        matrix rotmat(3, 3); // eigen vectors
	jacobi( mom, eigval, rotmat );

        rotate( crd, rotmat );
    }

} // namespace mort
