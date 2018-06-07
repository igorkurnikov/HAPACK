// HARLEM
// Igor Kurnikov
// University of Pittsburgh
//
// hapseudopot.cpp: 
// implementation of classes to manage core pseudopotentials
// 
// Created: December 1998
// Revised: December 1998
//////////////////////////////////////////////////////////////////////

#define HAPSEUDOPOT_CPP

#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

#include "haio.h"
#include "hapseudopot.h"


PseudoTerm::PseudoTerm()
{
	npower = 0;
	expon  = 0.0;
	coef   = 0.0; 
}



PseudoTerm::PseudoTerm(int new_npower, double new_expon, double new_coef)
{
	npower = new_npower;
	expon  = new_expon;
	coef   = new_coef; 
}

PseudoTerm::~PseudoTerm()
{
	
}

PseudoBlock::PseudoBlock()
{
	
}

PseudoBlock::~PseudoBlock()
{
	
}




HaPseudoPot::HaPseudoPot()
{

}

HaPseudoPot::~HaPseudoPot()
{

}

bool HaPseudoPot::SaveGaussInp(std::ostream& os) const
{
	int i;	
	if(at_symbol.empty()) return false;
	char buf[120];

	os << at_symbol.c_str() << " 0 " << std::endl;
	sprintf(buf,"%s %d %d ",pot_name.c_str(),max_angm,ncore_el);
	os << buf << std::endl;

	if(blocks.size() == 0) 
	{
		PrintLog(" Error in  HaPseudoPot::SaveGaussInp \n");
		PrintLog(" Empty  Description Blocks in the Psedopotential \n"); 
		PrintLog("Name %s for atom %s \n",pot_name.c_str(), at_symbol.c_str());
		return false;
	}

	vector<PseudoBlock>::const_iterator bitr;
	for(bitr=blocks.begin(); bitr != blocks.end(); bitr++)
	{
		os << (*bitr).description.c_str() << std::endl;

		int nt= (*bitr).terms.size();
		sprintf(buf,"%3d",nt);
		os << buf << std::endl;

		for( i=0; i < nt; i++)
		{
			const PseudoTerm* ptrm = &((*bitr).terms[i]);
			sprintf(buf,"%3d %16.9f %16.9f ",ptrm->npower,ptrm->expon,ptrm->coef);
			os << buf << std::endl;
		}
	}
	return true;
}

HaPseudoPotRef::HaPseudoPotRef()
{

}
 
HaPseudoPotRef::HaPseudoPotRef(const std::string & new_pot_name, const std::string & new_at_label)
{
	pot_name=new_pot_name; 
	boost::to_upper(pot_name);
	at_label=new_at_label; 
	boost::to_upper(at_label);
}
	
HaPseudoPotRef::~HaPseudoPotRef()
{

}

bool 
HaPseudoPotRef::operator == (const HaPseudoPotRef &  rhs) const
{
	if(pot_name != rhs.pot_name) 
		return false;
	return(at_label == rhs.at_label);
}

bool 
HaPseudoPotRef::operator < (const HaPseudoPotRef &  rhs) const 
{
	if(pot_name != rhs.pot_name) 
		return (pot_name < rhs.pot_name);
	return(at_label < rhs.at_label);
}




HaPseudoPotDB::HaPseudoPotDB()
{

}

HaPseudoPotDB::~HaPseudoPotDB()
{
	map<HaPseudoPotRef,HaPseudoPot*, less<HaPseudoPotRef> >::iterator pitr;
	
	for(pitr = dat.begin(); pitr != dat.end(); pitr++ )
	{
		delete( (*pitr).second);
	}
}


bool
HaPseudoPotDB::Init()
{
	InitHayWadt_1();
	return true;
}


HaPseudoPot* HaPseudoPotDB::Extract(const std::string & pot_name, const std::string & at_label) 
{
	HaPseudoPotRef pref(pot_name,at_label);

	map<HaPseudoPotRef,HaPseudoPot*, less<HaPseudoPotRef> >::iterator pitr;
	pitr=dat.find(pref);

	if(pitr == dat.end())
		return NULL;
	
	return ( (*pitr).second );

}

bool
HaPseudoPotDB::InitHayWadt_1()
{
	PseudoBlock block;
	HaPseudoPot* at_pot;
	
	int debug_level = 1;
	if(debug_level > 5)
		PrintLog(" HaPseudoPotDB::InitHayWadt_1() Pt 1 \n");

	at_pot = new HaPseudoPot;
	at_pot->pot_name="HAY_1";

// Be
// Pseudo potential for artificial donor/acceptor with E=-0.14 a.u.

	at_pot->at_symbol="Be";
	at_pot->max_angm = 1;
	at_pot->ncore_el = 2;

	at_pot->blocks.clear();

	block.terms.clear();
	block.description=         "be P potential fit";
 
	if(debug_level > 5)
		PrintLog(" HaPseudoPotDB::InitHayWadt_1() Pt 2 \n");

	block.terms.push_back(PseudoTerm( 2,        1.88230000,   -0.25828800));
	block.terms.push_back(PseudoTerm( 2,        6.05660000,   -2.09550900));
	block.terms.push_back(PseudoTerm( 1,       17.07800000,   -1.39610200));

	if(debug_level > 5)
		PrintLog(" HaPseudoPotDB::InitHayWadt_1() Pt 3 \n");

	at_pot->blocks.push_back(block);

	if(debug_level > 5)
		PrintLog(" HaPseudoPotDB::InitHayWadt_1() Pt 4 \n");

	block.terms.clear();
	block.description=         "be SP potential fit";
// coef of the first two terms were modified from 
// correspondingly (-28.574705, 35.735227) to push S-level up 0.1 a.u.  
    block.terms.push_back(PseudoTerm( 2,        1.51470000,   -25.0        ));
    block.terms.push_back(PseudoTerm( 2,        1.58470000,    40.0        ));
    block.terms.push_back(PseudoTerm( 1,        2.30910000,    -1.31942100 ));
    block.terms.push_back(PseudoTerm( 0,        4.20010000,     2.96972500 ));
	at_pot->blocks.push_back(block);
	
	if(debug_level > 5)
		PrintLog(" HaPseudoPotDB::InitHayWadt_1() Pt 5 \n");

	dat[HaPseudoPotRef(at_pot->pot_name,at_pot->at_symbol)]=at_pot;

	if(debug_level > 5)
		PrintLog(" HaPseudoPotDB::InitHayWadt_1() Pt 6 \n");

// Fe
	at_pot = new HaPseudoPot;
	at_pot->pot_name="HAY_1";

	at_pot->at_symbol="Fe";
	at_pot->max_angm = 3;
	at_pot->ncore_el = 18;

	at_pot->blocks.clear();

	block.terms.clear();
	block.description=         "fe f potential fit";
	block.terms.push_back(PseudoTerm( 1,      305.9908267,    -18.0000000));
	block.terms.push_back(PseudoTerm( 2,       57.5401809,   -111.5658307));
	block.terms.push_back(PseudoTerm( 2,       12.8787208,    -28.4011758));
	block.terms.push_back(PseudoTerm( 2,        3.3671745,    -10.1669409));
	block.terms.push_back(PseudoTerm( 2,        1.0778632,     -1.0631887));
	at_pot->blocks.push_back(block);

	block.terms.clear();
	block.description=         "fe s-f potential fit";
	block.terms.push_back(PseudoTerm( 0,       16.8558207,      3.0000000));
	block.terms.push_back(PseudoTerm( 1,        4.4141716,     17.4313169));
	block.terms.push_back(PseudoTerm( 2,        0.9390046,     29.4117858));
	block.terms.push_back(PseudoTerm( 2,        0.8294794,    -16.8563282));
	at_pot->blocks.push_back(block);

	block.terms.clear();
	block.description=         "fe p-f potential fit";
	block.terms.push_back(PseudoTerm( 0,       28.7494405,      5.0000000));
	block.terms.push_back(PseudoTerm( 1,       13.5008313,     18.1533743));
	block.terms.push_back(PseudoTerm( 2,       14.0936860,     54.0296157));
	block.terms.push_back(PseudoTerm( 2,        4.1439477,     48.7957799));
	block.terms.push_back(PseudoTerm( 2,        0.7948771,      6.2143968));
	at_pot->blocks.push_back(block);

 	block.terms.clear();
	block.description=         "fe d-f potential fit";
	block.terms.push_back(PseudoTerm( 2,        3.4156612,      0.2738868));
	block.terms.push_back(PseudoTerm( 2,        0.5623151,     -0.3947896));
	at_pot->blocks.push_back(block);

	dat[HaPseudoPotRef(at_pot->pot_name,at_pot->at_symbol)]=at_pot;

// Ru
	at_pot = new HaPseudoPot;
	at_pot->pot_name="HAY_1";

	at_pot->at_symbol="Ru";
	at_pot->max_angm = 3;
	at_pot->ncore_el = 36;

	at_pot->blocks.clear();

	block.terms.clear();
	block.description=         "Ru f potential fit";
	block.terms.push_back(PseudoTerm(  0,      239.2060503,     -0.0515270));
	block.terms.push_back(PseudoTerm(  1,       77.7081665,    -24.2111980));
	block.terms.push_back(PseudoTerm(  2,       24.1792246,    -90.1324868));
	block.terms.push_back(PseudoTerm(  2,       12.9011817,    -18.7412056));
	block.terms.push_back(PseudoTerm(  2,        5.5165833,    -21.3351883));
	block.terms.push_back(PseudoTerm(  2,        1.6355769,     -5.9996563));
	block.terms.push_back(PseudoTerm(  2,        0.4519993,     -0.2994102));
	at_pot->blocks.push_back(block);
	
	block.terms.clear();
	block.description=         "Ru s-f potential fit";
	block.terms.push_back(PseudoTerm(  0,       40.1276685,      2.9206085));
	block.terms.push_back(PseudoTerm(  1,       12.1750285,     34.7473053));
	block.terms.push_back(PseudoTerm(  2,        4.0290193,     41.6989147));
	block.terms.push_back(PseudoTerm(  2,        0.7120459,     29.8034158));
	block.terms.push_back(PseudoTerm(  2,        0.6306855,    -16.9340057));
	at_pot->blocks.push_back(block);

	block.terms.clear();
	block.description=         "Ru p-f potential fit";
	block.terms.push_back(PseudoTerm(  0,       32.6793541,      1.9652461));
	block.terms.push_back(PseudoTerm(  1,        9.7452864,     30.4935522));
	block.terms.push_back(PseudoTerm(  2,        2.7434318,     28.0446947));
	block.terms.push_back(PseudoTerm(  2,        0.6652098,      9.0079429));
	at_pot->blocks.push_back(block);

	block.terms.clear();
	block.description=         "Ru d-f potential fit";
	block.terms.push_back(PseudoTerm(  0,       57.6689621,      2.9793035));
	block.terms.push_back(PseudoTerm(  1,       62.3738213,     20.6159067));
	block.terms.push_back(PseudoTerm(  2,       21.9924916,    175.8361030));
	block.terms.push_back(PseudoTerm(  2,        2.8136093,     64.1834035));
	block.terms.push_back(PseudoTerm(  2,        2.4287338,    -40.1987590));
	block.terms.push_back(PseudoTerm(  2,        0.3765219,     -0.1880375));
	at_pot->blocks.push_back(block);

	dat[HaPseudoPotRef(at_pot->pot_name,at_pot->at_symbol)]=at_pot;

	if(debug_level > 5)
		PrintLog(" HaPseudoPotDB::InitHayWadt_1() Pt end \n");

	return true;

}
