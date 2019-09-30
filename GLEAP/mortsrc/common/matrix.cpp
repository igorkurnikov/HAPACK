#include <stdexcept>
#include <cassert>
#include "matrix.hpp"

namespace mort
{
    using std::out_of_range;

    using std::runtime_error;

    matrix::matrix(int n1, int n2) 
        : m_size1(n1), 
          m_size2(n2), 
          m_data(0.0, n1*n2) 
    {
    }

    matrix::matrix(const matrix& rhs) 
        : m_size1(rhs.m_size1), 
          m_size2(rhs.m_size2), 
          m_data( rhs.m_data )
    {
    }

    double& matrix::operator()( int i, int j ) 
    {
        if( i<0 || i >= m_size1 )
            throw out_of_range( "Error: matrix row id out of range" );

        if( j<0 || j >= m_size2 )
            throw out_of_range( "Error: matrix col id out of range" );

        return m_data[i*m_size2 + j];
    }
                
    double matrix::operator()(int i, int j ) const
    {
        if( i<0 || i >= m_size1 )
            throw out_of_range( "Error: matrix row id out of range" );

        if( j<0 || j >= m_size2 )
            throw out_of_range( "Error: matrix col id out of range" );

        return m_data[i*m_size2 + j];
    }

    bool is_identity( const double* mat, int size )
    {
        for( int i=0; i <size; ++i )
        {
            for( int j=0; j < size; ++j )
            {
                int id= size*i + j;
                if( i==j && std::abs(mat[id]-1.0) > 0.01 )
                    return false;
                
                if( i!=j && std::abs(mat[id]) > 0.01 )
                    return false;
            }
        }
        return true;
    }

    inline void ROTATE(matrix& a, int i,int j,int k,int l, double s, double tau)
    {
        double g=a(i,j);
        double h=a(k,l);
        a(i,j)=g-s*(h+g*tau);
        a(k,l)=h+s*(g-h*tau);
    }

    void jacobi( matrix& mA, numvec& vD, matrix& mV)
    {
        assert( mA.size1() == mA.size2()&& mV.size1() == mV.size2() );
        assert( mA.size1() == (int)vD.size() && mA.size1() == mV.size1() );
        
        int VSIZE = mA.size1();
        numvec vB(VSIZE), vZ(VSIZE);
        int	iNrot,	j, iq, ip, i;
        double	thresh, theta, tau, t, sm, s, h, g, c;

        for ( ip=0; ip<VSIZE; ip++ ) {
            for ( iq=0; iq<VSIZE; iq++ ) mV(ip,iq) = 0.0;
            mV(ip,ip) = 1.0;
        }

        for ( ip=0; ip<VSIZE; ip++ ) {
            vB[ip] = vD[ip] = mA(ip,ip);
            vZ[ip] = 0.0;
        }

        iNrot = 0;
        for ( i=1; i<=50; i++ ) {
            sm = 0.0;
            for ( ip=0; ip< VSIZE-1; ip++ ) {
                for ( iq=ip+1; iq<VSIZE; iq++ )
                    sm += fabs(mA(ip,iq));
            }
            if ( sm == 0.0 ) return;
            if ( i<4 )	thresh = 0.2*sm/(VSIZE*VSIZE);
            else		thresh = 0.0;
            for ( ip=0; ip<VSIZE-1; ip++ ) {
                for ( iq=ip+1; iq<VSIZE; iq++ ) {
                    g = 100.0*fabs(mA(ip,iq));
                    if ( i>4 && fabs(vD[ip])+g == fabs(vD[ip])
                         && fabs(vD[iq])+g  == fabs(vD[iq]))
                        mA(ip,iq) = 0.0;
                    else if ( fabs(mA(ip,iq)) > thresh ) {
                        h = vD[iq]-vD[ip];
                        if ( fabs(h)+g == fabs(h))
                            t = (mA(ip,iq))/h;
                        else {
                            theta = 0.5*h/(mA(ip,iq));
                            t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
                            if ( theta < 0.0 ) t = -t;
                        }
                        c = 1.0/sqrt(1+t*t);
                        s = t*c;
                        tau = s/(1.0+c);
                        h = t*mA(ip,iq);
                        vZ[ip] -= h;
                        vZ[iq] += h;
                        vD[ip] -= h;
                        vD[iq] += h;
                        mA(ip,iq) = 0.0;
                        for ( j=0; j<=ip-1; j++ ) {
                            ROTATE( mA, j, ip, j, iq, s, tau );
                        }
                        for ( j=ip+1; j<=iq-1; j++ ) {
                            ROTATE( mA, ip, j, j, iq, s, tau );
                        }
                        for ( j=iq+1; j<VSIZE; j++ ) {
                            ROTATE( mA, ip, j, iq, j, s, tau );
                        }
                        for ( j=0; j<VSIZE; j++ ) {
                            ROTATE( mV, j, ip, j, iq, s, tau );
                        }
                        ++iNrot;
                    }
                }
            }
            for ( ip=0; ip<VSIZE; ip++ ) {
                vB[ip] += vZ[ip];
                vD[ip] = vB[ip];
                vZ[ip] = 0.0;
            }
        }

        throw runtime_error( "Error: jacobi, too many iterations" );
    }

} // namespace mort

