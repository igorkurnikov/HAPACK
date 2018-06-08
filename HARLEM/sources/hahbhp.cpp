#include "hahbhp.h"
#include <string>
#include "command.h"


HBondAvg::HBondAvg(HaAtom* donor_atom, HaAtom* acceptor_atom, HaAtom* h_atom): 
   HaHBond(donor_atom,acceptor_atom,h_atom)
{
	min_distance = 1.0;
	max_distance = 5.0;
	avg_distance = 0.0;

	dutyCycle = 0.0;
}

HBondAvg::~HBondAvg()
{

}

void HBondAvg::SetAvgDistance(double dist)
{
	if(dist < 0.0)
		ErrorInMod("HBondAvg::SetAvgDistance()", "Average distances cannot be negative");
	this->avg_distance = dist;
}

void HBondAvg::SetDutyCycle(double duty_cycle_new)
{
	this->dutyCycle = duty_cycle_new;
}

void HBondAvg::SetMaxDistance(double maxDis)
{
	if(maxDis < 0.0)
		ErrorInMod("HBondAvg::SetMaxDistance()", "Maximum distance cannot be negative");
	this->max_distance = maxDis;
}

void HBondAvg::SetMinDistance(double minDis)
{
	if(minDis < 0.0)
		ErrorInMod("HBondAvg::SetMinDistance()", "Minimum distance cannot be negative");
	this->min_distance = minDis;
}

double HBondAvg::GetAvgDistance()
{
	return this->avg_distance;
}

double HBondAvg::GetDutyCycle()
{
	return this->dutyCycle;
}

double HBondAvg::GetMaxDistance()
{
	return this->max_distance;
}

double HBondAvg::GetMinDistance()
{
	return this->min_distance;
}

bool HBondAvg::operator == (const HBondAvg& hbond)
{
	if(this->src == hbond.src && this->dst == hbond.dst )
		return true;
	return false;
}

/* *****************************************************************************
 *	///////////////////////////////////////////////////////////////////////////
 *	// Implementation of HaHydrophobicTether
 *	///////////////////////////////////////////////////////////////////////////
 * ******************************************************************************
 */

HaHydrophobicTether::HaHydrophobicTether(HaAtom* fst_atom, HaAtom* sec_atom)
{
	firstAtom = fst_atom;
	secondAtom = sec_atom;

	min_distance = 1.0;
	max_distance = 5.0;
	avg_distance = 0.0;
	dutyCycle = 0.0;
}

HaHydrophobicTether::~HaHydrophobicTether()
{
}


double HaHydrophobicTether::GetAvgDistance()
{
	return this->avg_distance;
}


double HaHydrophobicTether::GetDevDistance()
{
	return this->dev_distance;
}

double HaHydrophobicTether::GetDutyCycle()
{
	return this->dutyCycle;
}

HaAtom* HaHydrophobicTether::GetFirstAtom()
{
	return this->firstAtom;
}

HaAtom* HaHydrophobicTether::GetSecondAtom()
{
	return this->secondAtom;
}

double HaHydrophobicTether::GetMaxDistance()
{
	return this->max_distance;
}

double HaHydrophobicTether::GetMinDistance()
{
	return this->min_distance;
}



void HaHydrophobicTether::SetAvgDistance(double dist)
{
	if(dist < 0.0)
		ErrorInMod("HaHydrophobicTether::SetAvgDistance()", "Average distances cannot be negative");
	this->avg_distance = dist;
}

void HaHydrophobicTether::SetDutyCycle(double duty_cycle_new)
{
	this->dutyCycle = duty_cycle_new;
}

void HaHydrophobicTether::SetMaxDistance(double maxDis)
{
	if(maxDis < 0.0)
		ErrorInMod("HaHydrophobicTether::SetMaxDistance()", "Maximum distance cannot be negative");
	this->max_distance = maxDis;
}

void HaHydrophobicTether::SetMinDistance(double minDis)
{
	if(minDis < 0.0)
		ErrorInMod("HaHydrophobicTether::SetMinDistance()", "Minimum distance cannot be negative");
	this->min_distance = minDis;
}

