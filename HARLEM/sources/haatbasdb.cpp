/*! \file haatbasdb.cpp
 
    Classes to manage Atomic basis set DataBase.

    \author Igor Kurnikov 
	\date 1998-2003
 
*/

#include <mpi.h>

#include <math.h>

#include "hastl.h"
#include "hastring.h"

#include <boost/algorithm/string.hpp>

#include "harlemapp.h"
#include "haatbasdb.h"
#include "haqchem.h"
#include "hapseudopot.h"
#include "haatombasis.h"


HaAtBasDB::HaAtBasDB(void)
{

}


HaAtBasDB::~HaAtBasDB(void)
{
	StrPtrMap::iterator itr;
	for( itr = bas_map.begin(); itr != bas_map.end(); itr++ )
	{
		GauAtomBasis* ptr = (GauAtomBasis*) (*itr).second;
		delete ptr;
	}
}

bool HaAtBasDB::Init(void)
{
	int debug_level = 1;
	Init321G();

	Init321PPG();
	InitHay1_dz();
	InitHay1_dzPP();	

	return true;
}


GauAtomBasis* HaAtBasDB::Extract(const std::string & BName, const std::string & Atomlbl)
{
	int debug_level = 5;
	if(!qc_db_init_flag)
	{
		qc_db_init_flag=true;
		pseudo_db.Init();
		bas_db.Init(); 
	}
	
	std::string tmp1 = BName;
	std::string tmp2 = Atomlbl;

	boost::trim(tmp1);
	boost::trim(tmp2);
	
	std::string ref = tmp1 + ":" + tmp2;
	
	boost::to_upper(ref);
	
	StrPtrMap::iterator bitr;
	
	bitr=bas_map.find(ref);	
	if(bitr != bas_map.end())
	{
		return( (GauAtomBasis*) (*bitr).second);
	}
	
	bool result=SetFromDaltonFile(BName,Atomlbl); // If basis is not in the database try
	if(result)                                    // to load form the Dalton Basis directory
	{
		bitr=bas_map.find(ref);	
		if(bitr != bas_map.end())
		{
			return((GauAtomBasis*)(*bitr).second);
		}	
	}
	return NULL;
}
  
bool HaAtBasDB::AddAtomBasis(const GauAtomBasis& atbas)
{
	std::string ref = atbas.GetBasName();
	ref += ":";
	ref += atbas.GetAtomType();
	boost::to_upper(ref);

	GauAtomBasis* p_atbas_new = new GauAtomBasis(atbas);

    if(bas_map.count(ref) != NULL ) 
	{
		GauAtomBasis* ptr = (GauAtomBasis*) bas_map[ref];
		delete ptr;
	}

	bas_map[ref]= p_atbas_new;
	return true;
}


bool
HaAtBasDB::Init321G()
{

// 3-21G in DALTON style First S then P the D

  const char* at_str[12]= {

"-H"
" S 2 1.0 "
" 5.4471780   0.156285"
" 0.82454724  0.904691"
" S 1 1.0 "
" 0.18319158  1.0"
"**** ",

"-He"
" S 2 1.0 "
" 13.62670  0.175230 "
"  1.99935  0.893483 "
" S 1 1.00 "
" 0.382993  1.0 "
"****",

"-Li "
" S 3 1.00 "
" 36.83820  0.0696686 "
"  5.48172  0.3813460 "
"  1.11327  0.6817020 "
" S 2 1.00 "
" 0.540205 -0.2631270 "
" 0.102255  1.1433900 "
" S 1 1.00 "
" 0.0285645 1.0 "
" P 2 1.00 "
" 0.540205  0.161546  "
" 0.102255  0.915663  "
" P 1 1.00 "
" 0.0285645 1.0  "
"**** ",

"-Be "
" S 3 1.0 "
" 71.88760 0.0644263 "
" 10.72890 0.3660960 "
"  2.22205 0.6959340 "
" S 2 1.0 "
"  1.29548 -0.421064 "
"  0.268881 1.224070 "
" S 1 1.0 "
" 0.0773501 1.0"
" P 2 1.0 "
" 1.29548  0.205132 "
" 0.268881 0.882528 "
" P 1 1.0"
" 0.0773501 1.0"
"****",

"-B"
" S 3 1.0 " 
" 116.434   0.0629605"
"  17.4314  0.3633040"
"   3.68016 0.6972550"
" S 2 1.0 "
" 2.281870 -0.368662 "
" 0.465248  1.199440 "
" S 1 1.0 "
" 0.124328  1.0 "
" P 2 1.0 "
" 2.281870  0.231152 "
" 0.465248  0.866764 "
" P 1 1.0 "
" 0.124328 1.0"
"****",

"-C"
" S 3 1.0 "
" 172.25600  0.0617669"
"  25.91090  0.358794 "
"   5.53335  0.700713 "
" S 2 1.0 "
"   3.664980 -0.395897"
"   0.770545  1.215840"
" S 1 1.0 "
" 0.195857  1.0"
" P 2 1.0 "
" 3.664980  0.236460"
" 0.770545  0.860619"
" P 1 1.0 "
" 0.195857  1.0"
"****",

"-N"
" S 3 1.0"
" 242.766 0.0598657"
" 36.4851 0.3529550"
" 7.81449 0.7065130"
" S 2 1.0"
" 5.42522 -0.413301"
" 1.14915  1.224420"
" S 1 1.0"
" 0.283205 1.0"
" P 2 1.0"
" 5.42522 0.237972"
" 1.14915 0.858953"
" P 1 1.0"
" 0.283205 1.0"
"****",

"-O"
" S 3 1.0"
" 322.037 0.0592394"
" 48.4308 0.3515000"
" 10.4206 0.7076580"
" S 2 1.0"
" 7.40294 -0.404453"
" 1.57620  1.221560"
" S 1 1.0"
" 0.373684 1.0"
" P 2 1.0"
" 7.40294 0.244586"
" 1.57620 0.853955"
" P 1 1.0"
" 0.373684 1.0"
"****",

"-F"
" S 3 1.0"
" 413.801 0.0585483"
" 62.2446 0.3493080"
" 13.4340 0.7096320"
" S 2 1.0"
" 9.77759 -0.407327"
" 2.08617  1.223140"
" S 1 1.0"
" 0.482383 1.0"
" P 2 1.0"
" 9.77759 0.246680"
" 2.08617 0.852321"
" P 1 1.0"
" 0.482383 1.0"
"****",

"-Ne"
" S 3 1.0"
" 515.724  0.058143"
" 77.6538  0.347951"
" 16.8136  0.710714"
" S 2 1.0"
" 12.48300 -0.409922"
"  2.66451  1.224310"
" S 1 1.0"
" 0.60625 1.0"
" P 2 1.0"
" 12.4830 0.247460"
" 2.66451 0.851743"
" P 1 1.0"
" 0.60625 1.0 "
"****",

"-P"
" S 3 1.0"
" 1054.9  0.0655407"
" 159.195 0.384036"
" 34.5304 0.674541"
" S 3 1.0"
" 44.2866 -0.10213"
" 10.1019 0.0815922"
" 2.73997 0.969788"
" S 2 1.0"
" 1.21865 -0.371495"
" 0.395546 1.27099"
" S 1 1.0"
" 0.122811 1.0"
" P 3 1.0"
" 44.2866 0.110851"
" 10.1019 0.456495"
" 2.73997 0.606936"
" P 2 1.0"
" 1.21865  0.0915823"
" 0.395546 0.934924"
" P 1 1.00"
" 0.122811 1.0"
"****",

"-S"
" S 3 1.0"
" 1210.62 0.0650071"
" 182.747 0.38204"
" 39.6673 0.676545"
" S 3 1.0"
" 52.2236 -0.10031"
" 11.9629 0.0650877"
" 3.28911 0.981455"
" S 2 1.0"
" 1.22384 -0.286089"
" 0.457303 1.22806"
" S 1 1.0"
" 0.142269 1.0"
" P 3 1.0"
" 52.2236 0.109646"
" 11.9629 0.457649"
" 3.28911 0.604261"
" P 2 1.0"
" 1.22384  0.164777"
" 0.457303 0.870855"
" P 1 1.0"
" 0.142269 1.0"
"****" 
};

  
for(int i=0; i < 12; i++)
{
	istrstream is(at_str[i]);
	GauAtomBasis  tmpbas;
	tmpbas.SetBasName("3-21G");
    tmpbas.SetFromGaussianInp(is);
	this->AddAtomBasis(tmpbas);
}

	return true;
}

bool
HaAtBasDB::Init321PPG()
{

// 3-21++G in DALTON style First S then P the D

  const char* at_str[12]= {

"-H"
" S 2 1.0 "
" 5.4471780   0.156285"
" 0.82454724  0.904691"
" S 1 1.0 "
" 0.18319158  1.0"
" S 1 1.0 "
" 0.03600000  1.0"

"**** ",

"-He"
" S 2 1.0 "
" 13.62670  0.175230 "
"  1.99935  0.893483 "
" S 1 1.00 "
" 0.382993  1.0 "
"****",

"-Li "
" S 3 1.00 "
" 36.83820  0.0696686 "
"  5.48172  0.3813460 "
"  1.11327  0.6817020 "
" S 2 1.00 "
" 0.540205 -0.2631270 "
" 0.102255  1.1433900 "
" S 1 1.00 "
" 0.0285645 1.0 "
" S 1 1.00 "
" 0.0074000 1.0 "
" P 2 1.00 "
" 0.540205  0.161546  "
" 0.102255  0.915663  "
" P 1 1.00 "
" 0.0285645 1.0  "
" P 1 1.00 "
" 0.0074000 1.0  "
"**** ",

"-Be "
" S 3 1.0 "
" 71.88760 0.0644263 "
" 10.72890 0.3660960 "
"  2.22205 0.6959340 "
" S 2 1.0 "
"  1.29548 -0.421064 "
"  0.268881 1.224070 "
" S 1 1.0 "
" 0.0773501 1.0"
" S 1 1.0 "
" 0.0207000 1.0"
" P 2 1.0 "
" 1.29548  0.205132 "
" 0.268881 0.882528 "
" P 1 1.0"
" 0.0773501 1.0"
" P 1 1.0"
" 0.0207000 1.0"
"****",

"-B"
" S 3 1.0 " 
" 116.434   0.0629605"
"  17.4314  0.3633040"
"   3.68016 0.6972550"
" S 2 1.0 "
" 2.281870 -0.368662 "
" 0.465248  1.199440 "
" S 1 1.0 "
" 0.124328  1.0 "
" S 1 1.0 "
" 0.031500  1.0 "
" P 2 1.0 "
" 2.281870  0.231152 "
" 0.465248  0.866764 "
" P 1 1.0 "
" 0.124328 1.0"
" P 1 1.0 "
" 0.031500 1.0 "
"****",

"-C"
" S 3 1.0 "
" 172.25600  0.0617669"
"  25.91090  0.358794 "
"   5.53335  0.700713 "
" S 2 1.0 "
"   3.664980 -0.395897"
"   0.770545  1.215840"
" S 1 1.0 "
" 0.195857  1.0"
" S 1 1.0 "
" 0.043800  1.0"
" P 2 1.0 "
" 3.664980  0.236460"
" 0.770545  0.860619"
" P 1 1.0 "
" 0.195857  1.0"
" P 1 1.0 "
" 0.043800  1.0"
"****",

"-N"
" S 3 1.0"
" 242.766 0.0598657"
" 36.4851 0.3529550"
" 7.81449 0.7065130"
" S 2 1.0"
" 5.42522 -0.413301"
" 1.14915  1.224420"
" S 1 1.0"
" 0.283205 1.0"
" S 1 1.0"
" 0.063900 1.0"
" P 2 1.0"
" 5.42522 0.237972"
" 1.14915 0.858953"
" P 1 1.0"
" 0.283205 1.0"
" P 1 1.0"
" 0.063900 1.0"
"****",

"-O"
" S 3 1.0"
" 322.037 0.0592394"
" 48.4308 0.3515000"
" 10.4206 0.7076580"
" S 2 1.0"
" 7.40294 -0.404453"
" 1.57620  1.221560"
" S 1 1.0"
" 0.373684 1.0"
" S 1 1.0"
" 0.084500 1.0"
" P 2 1.0"
" 7.40294 0.244586"
" 1.57620 0.853955"
" P 1 1.0"
" 0.373684 1.0"
" P 1 1.0"
" 0.084500 1.0"
"****",

"-F"
" S 3 1.0"
" 413.801 0.0585483"
" 62.2446 0.3493080"
" 13.4340 0.7096320"
" S 2 1.0"
" 9.77759 -0.407327"
" 2.08617  1.223140"
" S 1 1.0"
" 0.482383 1.0"
" S 1 1.0"
" 0.107600 1.0"
" P 2 1.0"
" 9.77759 0.246680"
" 2.08617 0.852321"
" P 1 1.0"
" 0.482383 1.0"
" P 1 1.0"
" 0.107600 1.0"
"****",

"-Ne"
" S 3 1.0"
" 515.724  0.058143"
" 77.6538  0.347951"
" 16.8136  0.710714"
" S 2 1.0"
" 12.48300 -0.409922"
"  2.66451  1.224310"
" S 1 1.0"
" 0.60625 1.0"
" S 1 1.0"
" 0.13000 1.0"
" P 2 1.0"
" 12.4830 0.247460"
" 2.66451 0.851743"
" P 1 1.0"
" 0.60625 1.0 "
" P 1 1.0"
" 0.13000 1.0"
"****",

"-P"
" S 3 1.0"
" 1054.9  0.0655407"
" 159.195 0.384036"
" 34.5304 0.674541"
" S 3 1.0"
" 44.2866 -0.10213"
" 10.1019 0.0815922"
" 2.73997 0.969788"
" S 2 1.0"
" 1.21865 -0.371495"
" 0.395546 1.27099"
" S 1 1.0"
" 0.122811 1.0"
" S 1 1.0"
" 0.034800 1.0"
" P 3 1.0"
" 44.2866 0.110851"
" 10.1019 0.456495"
" 2.73997 0.606936"
" P 2 1.0"
" 1.21865  0.0915823"
" 0.395546 0.934924"
" P 1 1.00"
" 0.122811 1.0"
" P 1 1.0"
" 0.034800 1.0"
"****",

"-S"
" S 3 1.0"
" 1210.62 0.0650071"
" 182.747 0.38204"
" 39.6673 0.676545"
" S 3 1.0"
" 52.2236 -0.10031"
" 11.9629 0.0650877"
" 3.28911 0.981455"
" S 2 1.0"
" 1.22384 -0.286089"
" 0.457303 1.22806"
" S 1 1.0"
" 0.142269 1.0"
" S 1 1.0"
" 0.040500 1.0"
" P 3 1.0"
" 52.2236 0.109646"
" 11.9629 0.457649"
" 3.28911 0.604261"
" P 2 1.0"
" 1.22384  0.164777"
" 0.457303 0.870855"
" P 1 1.0"
" 0.142269 1.0"
" P 1 1.0"
" 0.040500 1.0"
"****" 
};

  
for(int i=0; i < 12; i++)
{
	istrstream is(at_str[i]);
	GauAtomBasis  tmpbas;
	tmpbas.SetBasName("3-21++G");
    tmpbas.SetFromGaussianInp(is);
	this->AddAtomBasis(tmpbas);
}

	return true;
}


bool HaAtBasDB::InitHay1_dz()
{
	GauShell tmpshell;

// 3-21G for Light Atoms and Hay-Wadt double-zeta basis set 
// for Hays And Wadt pseudo-potenntial

	const std::string BNAME("HAY1_DZ");
	std::string atid;
	const int NLAT=12;
	int i;
	std::string lat_ids[NLAT] = { "H","He",
		                       "Li","Be","B","C","N","O","F","Ne",
							                     "P","S" };
	const GauAtomBasis* patb;

	const std::string bname_ref="3-21G";
// Set Basis for Light Atoms as in 3-21G basis

	for(i=0; i < NLAT; i++)
	{
		patb= Extract(bname_ref,lat_ids[i]);
		if( patb == NULL)
		{
			cerr << " HaAtBasDB::InitHay1_dz() " << endl; 
			cerr << " Can't Find in DB basis " << bname_ref << " for atom " 
				 << lat_ids[i] << endl;
			continue;
	
		}
		GauAtomBasis tbas(*patb);
		if(i==3) tbas.SetPseudoPotName("HAY_1"); // set Pseudopotential for Be 
		tbas.SetBasName(BNAME);
		this->AddAtomBasis(tbas);
	}

// Set Fe - Atom
	HaMat_double cf;

	GauAtomBasis  tmpbas;
	tmpbas.SetBasName(BNAME);
	tmpbas.SetPseudoPotName("HAY_1");


const char* Fe_str=
"-fe                 "
" S 2 1.0             "
" 0.5736   -0.4585154 "                                                
" 0.1021    1.2091119 "
" S 1 1.0             "
" 0.0363    1.0000000 "
" P 1 1.0             "
" 0.0740    1.0000000 "
" P 1 1.0             "
" 0.0220    1.0000000 " 
" D 4 1.0             "
" 37.0800   0.0316328 "
" 10.1000   0.1764010 "
" 3.2200    0.4527502 "
" 0.9628    0.5852844 "
" D 1 1.0             "
" 0.2262    1.0000000 "
"****                 "; 
        
if(1)
{ 
	istrstream is(Fe_str);
    tmpbas.SetFromGaussianInp(is);
}
	this->AddAtomBasis(tmpbas);
	tmpbas.ClearCoef();

// Set Ru Atom

const char* Ru_str=
"-Ru                  "
" S 2 1.0             "
" 0.3816   -1.1960626 "
" 0.1362    1.7273666 "
" S 1 1.0             "
" 0.0417    1.0000000 "
" P 2 1.0             "
" 0.5725   -0.0880864 "
" 0.0830    1.0283970 "
" P 1 1.0             "
" 0.0250    1.0000000 "
" D 3 1.0             "
" 4.1950    0.0583381 "
" 1.3770    0.4960916 "
" 0.4828    0.5824427 "        
" D 1 1.0             "
" 0.1501    1.0000000 "
" ****                ";

if(1)
{ 
	istrstream is(Ru_str);
 	tmpbas.SetFromGaussianInp(is);
}
	this->AddAtomBasis(tmpbas);
	tmpbas.ClearCoef();


	return true;
}

bool HaAtBasDB::InitHay1_dzPP()
{
	GauShell tmpshell;

// 3-21++G for Light Atoms and Hay-Wadt double-zeta basis set 
// for Hays And Wadt pseudo-potenntial

	const std::string BNAME("HAY1_DZPP");
	std::string atid;
	const int NLAT=12;
	int i;
	std::string lat_ids[NLAT] = { "H","He",
		                       "Li","Be","B","C","N","O","F","Ne",
							                     "P","S" };
	const GauAtomBasis* patb;

	const std::string bname_ref="3-21++G";
// Set Basis for Light Atoms as in 3-21++G basis

	for(i=0; i < NLAT; i++)
	{
		patb= Extract(bname_ref,lat_ids[i]);
		if( patb == NULL)
		{
			cerr << " HaAtBasDB::InitHay1_dzPP() " << endl; 
			cerr << " Can't Find in DB basis " << bname_ref << " for atom " 
				 << lat_ids[i] << endl;
			continue;
	
		}
		GauAtomBasis tbas(*patb);
		if(i==3) tbas.SetPseudoPotName("HAY_1"); // set Pseudopotential for Be to shift its IP energy
		tbas.SetBasName(BNAME);
		this->AddAtomBasis(tbas);
	}

// Set Basis for Heavy Atoms as in Double zeta Hay-Wadt pseudopotential basis 
	const int NMET=2;
	std::string met_ids[NMET] = { "Fe","Ru" };

	for(i=0; i < NMET; i++)
	{
		patb= Extract("HAY1_DZ", met_ids[i]);
		if( patb == NULL)
		{
			cerr << " HaAtBasDB::InitHay1_dzPP() " << endl; 
			cerr << " Can't Find in DB basis " << bname_ref << " for atom " 
				 << met_ids[i] << endl;
			continue;
	
		}
		GauAtomBasis tbas(*patb);
		tbas.SetPseudoPotName("HAY_1"); // set Pseudopotential for Be to shift its energy 
		tbas.SetBasName(BNAME);
		this->AddAtomBasis(tbas);
	}

	return true;
}


static char bufl[120];

static char* 
get_Dalton_Bas_file_line(ifstream & is)
{
	for(;;)
	{
		is.getline(bufl,120);
		if(is.eof() || is.fail() ) 
			return NULL;
		bool empty=true;
		for(int i=0; i < strlen(bufl); i++)
		{
			if(!isspace((int)bufl[i]))
			{
				empty=false;
				break;
			}
		}
		if(empty) continue;
		if(bufl[0] != '$' ) 
			return bufl;
	}
	return NULL;
}

bool HaAtBasDB::SetFromDaltonFile(const std::string& BName, const std::string & Atomlbl)
{
	std::string BasDir;

    PrintLog(" Harlem Home Dir = %s\n",pApp->harlem_home_dir.c_str());

	if(!pApp->harlem_home_dir.empty())
	{
		BasDir= pApp->harlem_home_dir;
		BasDir+="basis/";
	}
	
	std::string BName_full(BasDir);
	BName_full+=BName;
        
        PrintLog(" Full Name of the basis file: %s \n",BName_full.c_str());
        
//        BName_full = "/usr/local/lib/harlem/basis/3-21G";	
	ifstream is(BName_full.c_str());
	if(is.fail())
	{
		PrintLog("Error in HaAtBasDB::SetFromDaltonFile \n");
		PrintLog("File %s is not found \n",BName_full.c_str());
		PrintLog("File QUQU  is not found \n");
		return false;
	}
	
	char* buf;
	for(;;)
	{
		buf=get_Dalton_Bas_file_line(is);
		if(buf == NULL) 
			return false;
		// Find the beginning of first Atom Basis description
		if(*buf == 'A' || *buf == 'a')
		{
			istrstream str(++buf);
			int iat;
			str >> iat;	
			std::string Atomlbl_2=HaAtom::GetStdSymbolElem(iat);
			if(Atomlbl_2 == Atomlbl)
			{
				break;
			}
		}
			
	}
	
	// Start to read atom basis description 
	GauAtomBasis tmpbas;
	
	tmpbas.SetBasName(BName);
	tmpbas.SetAtomType(Atomlbl);
	
	// Cycle on Shell types 
	for( int l=0; l < 5; l++)
	{
		buf=get_Dalton_Bas_file_line(is);
		if( buf == NULL) 
			break;
		if(*buf == 'A' || *buf == 'a' )
			break; // start of the next atom basis description
		
		istrstream str2(buf);
		int ng; int maxexp; int nsh;
		
		str2 >> maxexp >> nsh;
		if(str2.fail())
		{
			break;
		}
		assert(maxexp > 0 && nsh > 0);
		HaVec_double exp(maxexp); 
		HaMat_double cf(maxexp,nsh);
		int iexp;
		for(iexp=1; iexp <= maxexp; iexp++)
		{
			buf=get_Dalton_Bas_file_line(is);
			istrstream istr3(buf);
			istr3 >> exp(iexp);
			
			// read subshell gaussian expansion coefficients:
			
			for(int  ms=1; ms <= nsh ; ms++)
			{
				istr3 >> cf(iexp,ms);
				if(istr3.fail())
				{
					cerr << " Error in HaAtBasDB::SetFromDaltonFile() " << endl;
					cerr << " Setting Basis " << BName << " for Atom " << Atomlbl << endl;
					return false;
				}
			}
		}
		
		
		// cut zero coefficients and fill shell coefficients:
		
		double dt[100];
		
		for( int ms =1; ms <= nsh ; ms++)
		{
			ng=0; 
			
			for( iexp=1; iexp <= maxexp; iexp++)
			{
				if( fabs(cf(iexp,ms)) < 1.0e-10)
					continue;
				ng++;
				dt[ng*2-2]=exp(iexp);
				dt[ng*2-1]=cf(iexp,ms);
			}
			GauShell tmpshell(l,ng);
			tmpshell.SetCoef(dt);
			tmpbas.AddShell(tmpshell);
		}
		
	}
	
	this->AddAtomBasis(tmpbas);
	
	return true;

}
