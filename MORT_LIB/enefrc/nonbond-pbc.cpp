#include <iostream>
#include <common.hpp>
#include <object.hpp>
#include <ambfmt.hpp>
#include "nonbond.hpp"
#include "ctrlparm.hpp"

#if defined(_MSC_VER)
//#include <boost/math/special_functions/erf.hpp>
//using boost::math::erfc;
//static double erfc(double x) { return boost::math::erfc(x); }
#else
#ifndef __GNUC__
extern "C" double erfc(double);
#endif
#endif

namespace mort
{
    double img_dis2(const double* v1, const double* v2, const numvec& box)
    {
         double v12[3];

         v12[0] = v1[0] - v2[0];
         v12[1] = v1[1] - v2[1];
         v12[2] = v1[2] - v2[2];
 
         assert( v12[0] >= -box[0] && v12[0] < box[0] );
         assert( v12[1] >= -box[1] && v12[1] < box[1] );
         assert( v12[2] >= -box[2] && v12[2] < box[2] );

         v12[0] = std::abs( v12[0] );
         v12[1] = std::abs( v12[1] );
         v12[2] = std::abs( v12[2] );
    
         if( v12[0] > 0.5*box[0] ) v12[0] -= box[0];
         if( v12[1] > 0.5*box[1] ) v12[1] -= box[1];
         if( v12[2] > 0.5*box[2] ) v12[2] -= box[2];

         return v12[0]*v12[0] + v12[1]*v12[1] + v12[2]*v12[2];
    }

    numvec img_frac(const double* frac, const numvec& box)
    {
        numvec imgcrd(3);
        imgcrd[0] = (frac[0]+0.5)*box[0];
        imgcrd[1] = (frac[1]+0.5)*box[1];
        imgcrd[2] = (frac[2]+0.5)*box[2];
        return imgcrd;
    }

    void get_frac(int natm3, const double* pos, const numvec& box, double* frac)
    {
        for( int i=0; i < natm3; ++i )
        {
            int j = (i%3);
            frac[i] = pos[i]/box[j];
            frac[i] = frac[i] - floor(frac[i]);
            assert( frac[i] >= 0.0 && frac[i] < 1.0 );
            if( frac[i] >= 0.5 ) frac[i] -= 1.0;
            assert( frac[i] >= -0.5 && frac[i] < 0.5 );
        }
    }

    double get_ewcoef( double cut, double dtol )
    {
        double lo = 0.0;
        double hi = 0.5;

        while( erfc(hi)/cut > dtol )
        {
            hi *= 2.0;
        }
        
        double mi = (lo+hi) * 0.5;
        double vm = erfc(mi) / cut;
        while( std::abs( vm-dtol )>1e-10 )
        {
            if( vm > dtol )
            {
                lo = mi;
            }
            else
            {
                hi = mi;
            }
            
            mi = (lo+hi) * 0.5;
            vm = erfc(mi)/ cut;
        }

        return mi/cut;
    }

    numvec get_dir(const molecule_t& m, const ctrlparm_t& p)
    {
        excl_t excl;
        exclude(m, excl, 3);

        nonbond_t nb(m);
        nb.list_for_pbc(p.cut);

        const vector<double>& pchg = get_dvec(m, ATOM, PCHG);
	
        double eadj=0.0;
	double edir=0.0;
        double dsum_tol = 1.e-5;
        double ewcoef = get_ewcoef(p.cut, dsum_tol);
        for(int i=0; i<(int)nb.list.size(); ++i)
	{
	    double chgi = pchg[i];
	    for(int k=0; k<(int)nb.list[i].size();++k)
	    {
	        int j = nb.list[i][k];
                double chgj = pchg[j];
	        double dist = nb.dist[i][k];
		if( excl.found(i, j) )
		{
	            eadj += INVCHG2*chgi*chgj*(erfc(dist*ewcoef)-1 )/dist;
		}
		else
		{
		    edir += INVCHG2*chgi*chgj*erfc(dist*ewcoef)/dist;
		}
	    }
	}

	return makevec(edir, eadj);
    }

    namespace regewald
    {

        int max_of_three(const vector<int>& v)
        {
            assert(v.size()==3);
            int res = std::max( v[0], v[1] );
            return std::max( res, v[2] );
        }

        double eval(const numvec& box, int natm3, const double* frac, const double* chg,
                    const vector<int>& mlimit, double maxexp2, double ewaldcoef)
        {
            int natom = natm3/3;
            int mmax  = max_of_three(mlimit) + 1;
       
            double fac = PI2/(ewaldcoef*ewaldcoef);
            double volume = box[0]*box[1]*box[2];

            vector< vector<double> > cosf1( natom, vector<double>(mmax) );
            vector< vector<double> > cosf2( natom, vector<double>(mmax) );
            vector< vector<double> > cosf3( natom, vector<double>(mmax) );
            vector< vector<double> > sinf1( natom, vector<double>(mmax) );
            vector< vector<double> > sinf2( natom, vector<double>(mmax) );
            vector< vector<double> > sinf3( natom, vector<double>(mmax) );
        
            for( int i=0; i < natom; ++i )
            {
                cosf1[i][0] = 1.0;
                cosf2[i][0] = 1.0;
                cosf3[i][0] = 1.0;
            
                cosf1[i][1] = cos( TWOPI*frac[3*i] );
                cosf2[i][1] = cos( TWOPI*frac[3*i+1] );
                cosf3[i][1] = cos( TWOPI*frac[3*i+2] );
            
                sinf1[i][0] = 0.0;
                sinf2[i][0] = 0.0;
                sinf3[i][0] = 0.0;
            
                sinf1[i][1] = sin( TWOPI*frac[3*i] );
                sinf2[i][1] = sin( TWOPI*frac[3*i+1] );
                sinf3[i][1] = sin( TWOPI*frac[3*i+2] );

                for( int m=2; m < mmax; ++m )
                {
                    cosf1[i][m] = cosf1[i][m-1]*cosf1[i][1] - sinf1[i][m-1]*sinf1[i][1];
                    cosf2[i][m] = cosf2[i][m-1]*cosf2[i][1] - sinf2[i][m-1]*sinf2[i][1];
                    cosf3[i][m] = cosf3[i][m-1]*cosf3[i][1] - sinf3[i][m-1]*sinf3[i][1];
                
                    sinf1[i][m] = sinf1[i][m-1]*cosf1[i][1] + cosf1[i][m-1]*sinf1[i][1];
                    sinf2[i][m] = sinf2[i][m-1]*cosf2[i][1] + cosf2[i][m-1]*sinf2[i][1];
                    sinf3[i][m] = sinf3[i][m-1]*cosf3[i][1] + cosf3[i][m-1]*sinf3[i][1];
                }    
            }

            double e_recip = 0.0;
            for( int m1=0; m1 <= mlimit[0]; m1++)
            {
                double mult = (m1==0)? 1.0 :2.0;
            
                for( int m2=-mlimit[1]; m2 <= mlimit[1]; ++m2 )
                {
                    vector< double > c12( natom );
                    vector< double > s12( natom );

                    if( m2 <0 )
                    {
                        for( int i=0; i < natom; ++i )
                        {
                            c12[i] = cosf1[i][m1]*cosf2[i][-m2] + sinf1[i][m1]*sinf2[i][-m2];
                            s12[i] = sinf1[i][m1]*cosf2[i][-m2] - cosf1[i][m1]*sinf2[i][-m2];
                        }
                    }
                    else
                    {
                        for( int i=0; i < natom; ++i )
                        {
                            c12[i] = cosf1[i][m1]*cosf2[i][m2] - sinf1[i][m1]*sinf2[i][m2];
                            s12[i] = sinf1[i][m1]*cosf2[i][m2] + cosf1[i][m1]*sinf2[i][m2];
                        }
                    }

                    for( int m3=-mlimit[2]; m3 <= mlimit[2]; ++m3 )
                    {
                        double mhat1 = m1/box[0];
                        double mhat2 = m2/box[1];
                        double mhat3 = m3/box[2];
                        double msq = mhat1*mhat1 + mhat2*mhat2 + mhat3*mhat3;
                        double denom = M_PI*volume*msq;
                        double eterm = (m1*m1+m2*m2+m3*m3>0) ? mult*exp(-fac*msq)/denom : 0.0;

                        if( msq <maxexp2 )
                        {
                            vector< double > c(natom);
                            vector< double > s(natom);
                        
                            if( m3 < 0 )
                            {
                                for(int i=0; i < natom; ++i )
                                {
                                    c[i] = c12[i]*cosf3[i][-m3] + s12[i]*sinf3[i][-m3];
                                    s[i] = s12[i]*cosf3[i][-m3] - c12[i]*sinf3[i][-m3];
                                }
                            }
                            else
                            {
                                for(int i=0; i < natom; ++i )
                                {
                                    c[i] = c12[i]*cosf3[i][m3] - s12[i]*sinf3[i][m3];
                                    s[i] = s12[i]*cosf3[i][m3] + c12[i]*sinf3[i][m3];
                                }
                            }
                        
                            double cstruct = 0.0;
                            double sstruct = 0.0;
                            for( int i=0; i < natom; ++i )
                            {
                                cstruct += 18.2223*chg[i]*c[i];
                                sstruct += 18.2223*chg[i]*s[i];
                            }
                            double struc2 = cstruct*cstruct + sstruct*sstruct;

                            e_recip += eterm*struc2;
                        }
                    }
                }
            }
        
            return e_recip*0.5;
        }
    
        double rene( double x, double ewaldcoef )
        {
            double y = M_PI*x/ewaldcoef;
            return 2.0*ewaldcoef*erfc(y)*INVSQRTPI;
        }

        double find_maxexp( double ewaldcoef, double rtol )
        {
            double xlo = 0.0;
            double xhi = 1.0;
            while( rene( xhi, ewaldcoef ) > rtol )
            {
                xhi *= 2.0;
            }
            
            double xmi = (xlo+xhi)*0.5;
            double emi = rene(xmi, ewaldcoef);
            
            while( std::abs(emi-rtol) > 1e-8 && std::abs(xhi-xlo) > 1e-8 )
            {
                if(emi > rtol)
                {
                    xlo = xmi;
                }
                else
                {
                    xhi = xmi;
                }
                
                xmi = (xlo+xhi)*0.5;
                emi = rene(xmi, ewaldcoef);
            }
            
            return xmi;
        }
        
        vector<int> find_mlimit(const numvec& box, double maxexp, double eigmin )
        {
            int mtop1 = int( box[0]*maxexp/sqrt(eigmin) );
            int mtop2 = int( box[1]*maxexp/sqrt(eigmin) );
            int mtop3 = int( box[2]*maxexp/sqrt(eigmin) );
            
            vector<int> mlimit(3, 0);
            
            int nvecs=0;
            for( int m1=-mtop1; m1 <= mtop1; ++m1 )
            {
                for( int m2=-mtop2; m2 <= mtop2; ++m2 )
                {
                    for( int m3=-mtop3; m3 <= mtop3; ++m3 )
                    {
                        double mhat1 = m1/box[0];
                        double mhat2 = m2/box[1];
                        double mhat3 = m3/box[2];                    
                        double msq = mhat1*mhat1 + mhat2*mhat2 + mhat3*mhat3;

                        if( msq <= maxexp*maxexp ) 
                        {
                            nvecs++;
                            if( std::abs(m1) > mlimit[0] ) mlimit[0] = std::abs(m1);
                            if( std::abs(m2) > mlimit[1] ) mlimit[1] = std::abs(m2);
                            if( std::abs(m3) > mlimit[2] ) mlimit[2] = std::abs(m3);
                        }
                    }
                }
            }

            std::cout << "nvecs: " << nvecs << std::endl;
            return mlimit;
        }    

    } // namespace regewald

    using namespace regewald;
   

    double eval_regewald( const molecule_t& mol, double cut )
    {
        numvec box = mol.get_v(BOX);
        
        const vector<double>& pos = get_vvec(mol, ATOM, POSITION);
        const vector<double>& chg = get_dvec(mol, ATOM, PCHG);

        vector<double> frac( 3*mol.natom() );
        get_frac(frac.size(), &pos[0], box, &frac[0]);
        std::cout << "frac got" << std::endl;

        double dsum_tol = 1.0e-5;
        double rsum_tol = 5.0e-5;
        double eigmin = 0.5;

        double ewaldcoef = get_ewcoef( cut, dsum_tol );
        std::cout << "ewaldcoef: " << ewaldcoef << std::endl;

        double maxexp = find_maxexp( ewaldcoef, rsum_tol );
        std::cout << "maxexp: " << maxexp << std::endl;
       
        vector<int> mlimit = find_mlimit( box, maxexp, eigmin );
        std::cout << "mlimit:" << mlimit[0] << " "  << mlimit[1] << " " << mlimit[2] << std::endl;
       
        double e_recip = regewald::eval( box, frac.size(), &frac[0], &chg[0], mlimit, maxexp*maxexp, ewaldcoef );

        std::cout << "Recip: " << format( "%20.4f" ) % e_recip << std::endl;

	return e_recip;
    }


} // namespace mort
