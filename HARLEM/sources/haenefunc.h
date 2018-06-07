/*!  \file haenefunc.h

    Classes to model energy functional in HARLEM  
 
    \author Igor Kurnikov 
    \date 1999-2009
*/

#ifndef HAENEFUNC_H
#define HAENEFUNC_H

class harlem::Coord;

class HaEnergyFunc
{
public:
	virtual double ComputeEnergy(harlem::Coord* pcrd) = 0;
};


#endif // !defined(HAENEFUNC_H)


