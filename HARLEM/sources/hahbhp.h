#ifndef HAHBPH_H
#define HAHBPH_H

#include <iostream>
#include <string>
using namespace std;

#include "haatom.h"
#include "habond.h"


class HBondAvg: public HaHBond
//! Class for H-Bond averaged over MD trajectory
{
protected:

	double	min_distance;		// Smallest distance during MD
	double	max_distance;		// Largest distance during MD
	double	avg_distance;		// Averge distance during MD
	double	dev_distance;		// Metric for distance deviation

	double	dutyCycle;			// Duty Cycle

public:
	HBondAvg(HaAtom* donor_atom = NULL, HaAtom* acceptor_atom = NULL, HaAtom* h_atom = NULL );
	~HBondAvg();

	double GetMinDistance();
	double GetMaxDistance();
	double GetAvgDistance();
	double GetDevDistance();
	double GetDutyCycle();

	void SetMinDistance(double);
	void SetMaxDistance(double);
	void SetAvgDistance(double);
	void SetDutyCycle(double duty_cycle_new);

	// operators
	bool operator == (const HBondAvg& hbond);
};

class HaHydrophobicTether
{
protected:
	HaAtom*	firstAtom;		// The first atom in the tether
	HaAtom*	secondAtom;		// The second Atom

	double	min_distance;		// Minimum distance
	double	max_distance;		// Maximum distance
	double	avg_distance;		// Average distance
	double	dev_distance;		// Distance Deviation
	
	double	dutyCycle;			// Duty Cycle of the PH

public:
	HaHydrophobicTether(HaAtom* atom1=NULL, HaAtom* atom2=NULL);
	~HaHydrophobicTether();

	HaAtom* GetFirstAtom();
	HaAtom* GetSecondAtom();

	double GetMinDistance();
	double GetMaxDistance();
	double GetAvgDistance();
	double GetDevDistance();
	double GetDutyCycle();

	void SetMinDistance(double);
	void SetMaxDistance(double);
	void SetAvgDistance(double);
	void SetDutyCycle(double duty_cycle_new);

};



#endif


