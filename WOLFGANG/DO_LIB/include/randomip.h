#ifndef RANDOMIP_H
#define RANDOMIP_H

#include <math.h>

class Random
{
    int inext;
    int inextp;
    long ma[56];
    long mj;

public:

    Random(long idum)
    {
	const long   MBIG  = 1000000000;
	const long   MSEED = 161803398;
	const long   MZ    = 0;
	const double FAC   = 1.0/MBIG;

	int i,ii,k;
	long mk;
	
	mj=MSEED-(idum < 0 ? -idum : idum);
	mj %= MBIG;
	ma[55]=mj;
	mk=1;
	for (i=1;i<=54;i++) 
	{
	    ii=(21*i) % 55;
	    ma[ii]=mk;
	    mk=mj-mk;
	    if (mk < MZ) mk += MBIG;
	    mj=ma[ii];
	}
	for(k=1;k<=4;k++)
	    for (i=1;i<=55;i++) 
	    {
		ma[i] -= ma[1+(i+30) % 55];
		if (ma[i] < MZ) ma[i] += MBIG;
	    }
	inext=0;
	inextp=31;
    }

double operator()()
    {
	const long   MBIG  = 1000000000;
	const double FAC   = 1.0/MBIG;
	if (++inext  == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < 0L) mj += MBIG;
	ma[inext] =mj;
	return mj*FAC;
    }

double GetRndNum()
   {
	return (*this)();
   }
};



class RandomGauss : public Random
{
    int iset;
    double gset;
    
public:

    RandomGauss(long idum) : Random(idum)
    {
	iset = 0;
    }
    
double operator()()
    {
	double fac,rsq,v1,v2;

	if  (iset == 0) 
	{
	    do 
	    {
		v1=2.0*Random::operator()()-1.0;
		v2=2.0*Random::operator()()-1.0;
		rsq=v1*v1+v2*v2;
	    } while (rsq >= 1.0 || rsq == 0.0);
	    fac= sqrt(-2.0*log(1.0*rsq)/rsq);
	    gset=v1*fac;
	    iset=1;
	    return v2*fac;
	} else 
	{
	    iset=0;
	    return gset;
	}
    }
};


#endif
