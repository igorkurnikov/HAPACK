/*!  \file mm_driver_gromacs.cpp

    Molecular Mechanics simulations using GROMACS package  

    \author Igor Kurnikov 
    \date 2010-

*/
#include <mpi.h>

#include <tinyxml.h>

#include <boost/algorithm/string.hpp>

#include "harlemapp.h"
#include "hamolmech.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "mm_driver_gromacs.h"

MMDriverGromacs::MMDriverGromacs(HaMolMechMod* p_mm_mod_new)
{
	p_mm_mod = p_mm_mod_new;
	p_mm_model    = p_mm_mod->p_mm_model;
	pmset = p_mm_mod->GetMolSet();

	to_save_input_files = TRUE;
}

MMDriverGromacs::~MMDriverGromacs()
{

}

int MMDriverGromacs::SaveAllInpFiles()
{	
	SaveTopFile();
	to_save_input_files = FALSE;
	return TRUE;
}

int MMDriverGromacs::SaveTopFile()
{
	std::string gmx_top_fname = "system_gmx.itp";
	if(p_mm_model->to_init_mm_model) p_mm_mod->InitMolMechModel();
	std::ofstream os(gmx_top_fname.c_str());
	if(os.fail()) 
	{
		PrintLog("Error in MMDriverAmber::SaveTopFile() \n");
		PrintLog("Can't create file %s \n", gmx_top_fname.c_str());
		return FALSE;
	}
	return SaveTopToStream(os);
}

int MMDriverGromacs::SaveTopToStream(std::ostream& os)
{
	int i;
	char buf[256];
	
	try
	{
		AtomIntMap& at_idx_map = p_mm_model->GetAtIdxMap(TRUE);
		os << "; system.itp created by HARLEM" << std::endl;
		os << std::endl;

		SaveAtomTypesToStream(os);

		os << "[ moleculetype ]" << std::endl;
		os << ";name            nrexcl" << std::endl;
		os << " REF            3 " << std::endl;
		os << "  " << std::endl;

		SaveAtomsToStream(os);
		SaveBondsToStream(os);
		SavePairsToStream(os);
		SaveAnglesToStream(os);
		SaveDihedralsToStream(os);
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in MMDriverGromacs::SaveTopToStream() \n");
		PrintLog("%s\n",ex.what());
		return FALSE;
	}
	return FALSE;
}

int MMDriverGromacs::SaveAtomTypesToStream(std::ostream& os)
{
		char buf[256];
		std::set<std::string> at_ff_names;

		os << "[ atomtypes ]" << std::endl;
		os << ";name   bond_type     mass     charge   ptype   sigma         epsilon       Amb" << std::endl;
//		os << ";   c        c           0.00000  0.00000   A     3.39967e-01   3.59824e-01 ; 1.91  0.0860" << std::endl;

		int i;
		for( i = 0; i < p_mm_model->Atoms.size(); i++ )
		{
			std::string ff_s = p_mm_model->Atoms[i]->GetFFSymbol();
			boost::trim(ff_s);
			if( ff_s.empty() ) continue;
			if( at_ff_names.count(ff_s) > 0 ) continue;
			at_ff_names.insert(ff_s);

			double vdw_rad = p_mm_model->Atoms[i]->vdw_rad;
			double ew = p_mm_model->Atoms[i]->ew;

			double sigma   = 2*vdw_rad*0.1/pow(2.0,1.0/6.0);
			double epsilon = 4.184*ew;
			
			sprintf( buf, "%4s %4s    0.00000    0.00000  A  %14.5e %14.5e  ;  %9.5f  %9.6f ",
				           ff_s.c_str(),ff_s.c_str(), sigma, epsilon, vdw_rad, ew );
			os << buf << std::endl;
		}
		os << "  " << std::endl;
		return TRUE;
}

int MMDriverGromacs::SaveAtomsToStream(std::ostream& os)
{
	char buf[256];
		
	os << "[ atoms ]" << std::endl;
	os << ";   nr  type   resi  res  atom    cgnr     charge      mass       ; qtot   bond_type" << std::endl;
	double qtot = 0.0;
	for( int i = 0; i < p_mm_model->Atoms.size(); i++)
	{
		HaAtom* aptr = p_mm_model->Atoms[i];
		std::string atn = aptr->GetName();
		double ch = aptr->GetCharge();
		double mass = aptr->GetMass();
		HaResidue* pres = aptr->GetHostRes();
		std::string ff_s = aptr->GetFFSymbol();
		std::string res_name = pres->GetName();
		int resi = pres->GetSerNo();
		qtot += ch;
		sprintf(buf,"%6d %5s%6d %5s %5s %6d %12.6f %12.6f ; qtot %7.3f ",i+1,ff_s.c_str(),resi,res_name.c_str(),atn.c_str(),i+1,ch,mass,qtot);
		os << buf << std::endl;
	}
	os << "  " << std::endl;

	return TRUE;
}

int MMDriverGromacs::SaveBondsToStream(std::ostream& os)
{
	char buf[256];
	AtomIntMap& at_idx_map = p_mm_model->GetAtIdxMap(0);
	os << "[ bonds ]" << std::endl;
	os << ";   ai     aj funct   r             k   " << std::endl;

	set<MMBond, less<MMBond> >::iterator mbitr = p_mm_model->MBonds.begin();
	for(; mbitr != p_mm_model->MBonds.end(); mbitr++)
	{
		MMBond& bnd = (MMBond&) *mbitr;		
		HaAtom* aptr1 = (HaAtom*) bnd.pt1;
		HaAtom* aptr2 = (HaAtom*) bnd.pt2;

		if( at_idx_map.count(aptr1) == 0 ) continue;
		if( at_idx_map.count(aptr2) == 0 ) continue;

		int idx1 = at_idx_map[aptr1]+1;  // convert to 1-based index
		int idx2 = at_idx_map[aptr2]+1;  // convert to 1-based index

		double r0 = bnd.r0 * 0.1; // Ang -> nm
		double fc = bnd.fc*4.184*4.184*100; // kcal/mol/Ang^2  ->kJ/nm^2

		std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);

		sprintf(buf,"%6d  %6d   1 %14.4e%14.4e ; %s - %s ",idx1,idx2, r0,fc, at_lbl_1.c_str(), at_lbl_2.c_str());
		os << buf << std::endl;
	}
	os << "  " << std::endl;
	return TRUE;
}

int MMDriverGromacs::SavePairsToStream(std::ostream& os)
{
	char buf[256];
	AtomIntMap& at_idx_map = p_mm_model->GetAtIdxMap(0);
	os << "[ pairs ]" << std::endl;
	os << ";   ai     aj    funct" << std::endl;

	vector<MMDihedral>::iterator ditr;
	for( ditr = p_mm_model->Dihedrals.begin(); ditr != p_mm_model->Dihedrals.end(); ditr++)
	{
		HaAtom* aptr1 = (HaAtom*) (*ditr).pt1;
		HaAtom* aptr2 = (HaAtom*) (*ditr).pt4;

		if( at_idx_map.count(aptr1) == 0 ) continue;
		if( at_idx_map.count(aptr2) == 0 ) continue;

		int idx1 = at_idx_map[aptr1]+1;  // convert to 1-based index
		int idx2 = at_idx_map[aptr2]+1;  // convert to 1-based index

		std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);

		sprintf(buf,"%6d  %6d  1 ; %s - %s ",idx1,idx2, at_lbl_1.c_str(), at_lbl_2.c_str());
		os << buf << std::endl;
	}
	os << "  " << std::endl;
	return TRUE;
}

int MMDriverGromacs::SaveAnglesToStream(std::ostream& os)
{
	char buf[256];
	AtomIntMap& at_idx_map = p_mm_model->GetAtIdxMap(0);
	os << "[ angles ]" << std::endl;
	os << ";   ai     aj     ak    funct   theta         cth " << std::endl;

	set<MMValAngle, less<MMValAngle> >::iterator vaitr = p_mm_model->ValAngles.begin();
	for(; vaitr != p_mm_model->ValAngles.end(); vaitr++)
	{
		MMValAngle& ang = (MMValAngle&)(*vaitr);
		
		double a0 = ang.a0;
		double fc = ang.fc * 4.184;

		HaAtom* aptr1 = (HaAtom*) ang.pt1;
		HaAtom* aptr2 = (HaAtom*) ang.pt2;
		HaAtom* aptr3 = (HaAtom*) ang.pt3;

		if( at_idx_map.count(aptr1) == 0 ) continue;
		if( at_idx_map.count(aptr2) == 0 ) continue;
		if( at_idx_map.count(aptr3) == 0 ) continue;

		int idx1 = at_idx_map[aptr1]+1;  // convert to 1-based index
		int idx2 = at_idx_map[aptr2]+1;  // convert to 1-based index
		int idx3 = at_idx_map[aptr3]+1;  // convert to 1-based index

		std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_3 = aptr3->GetRef(HaAtom::ATOMREF_STD);

		sprintf(buf,"%6d %6d %6d  1  %14.4e%14.4e; %s - %s - %s ",idx1,idx2,idx3, a0,fc,at_lbl_1.c_str(), at_lbl_2.c_str(),at_lbl_3.c_str());
		os << buf << std::endl;
	}
	os << "  " << std::endl;
	return TRUE;
}

int MMDriverGromacs::SaveDihedralsToStream(std::ostream& os)
{
	char buf[256];
	AtomIntMap& at_idx_map = p_mm_model->GetAtIdxMap(0);
	
	os << "[ dihedrals ] ; proper " << std::endl;
	os << ";    i      j      k      l   func   phase     kd      pn" << std::endl;

	int i;
	std::vector<MMDihedral>::iterator ditr;		
	for(ditr = p_mm_model->Dihedrals.begin(); ditr != p_mm_model->Dihedrals.end(); ditr++)
	{
		MMDihedral& dih = (MMDihedral&)(*ditr);
		
		std::vector<double> pk = dih.pk;
		std::vector<double> pn = dih.pn;
		std::vector<double> phase = dih.phase;

		HaAtom* aptr1 = (HaAtom*) dih.pt1;
		HaAtom* aptr2 = (HaAtom*) dih.pt2;
		HaAtom* aptr3 = (HaAtom*) dih.pt3;
		HaAtom* aptr4 = (HaAtom*) dih.pt4;

		if( at_idx_map.count(aptr1) == 0 ) continue;
		if( at_idx_map.count(aptr2) == 0 ) continue;
		if( at_idx_map.count(aptr3) == 0 ) continue;
		if( at_idx_map.count(aptr4) == 0 ) continue;

		int idx1 = at_idx_map[aptr1]+1;  // convert to 1-based index
		int idx2 = at_idx_map[aptr2]+1;  // convert to 1-based index
		int idx3 = at_idx_map[aptr3]+1;  // convert to 1-based index
		int idx4 = at_idx_map[aptr4]+1;  // convert to 1-based index

		std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_3 = aptr3->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_4 = aptr4->GetRef(HaAtom::ATOMREF_STD);

		for( i = 0; i < pk.size(); i++)
		{
			pk[i] *= 4.184;
			sprintf(buf,"%6d %6d %6d %6d   1  %14.4e%14.4e %2.0f  ; %s - %s - %s - %s ",idx1,idx2,idx3,idx4, phase[i], pk[i], pn[i], 
				at_lbl_1.c_str(), at_lbl_2.c_str(),at_lbl_3.c_str(),at_lbl_4.c_str());
			os << buf << std::endl;
		}
	}
	os << "  " << std::endl;

	os << "[ dihedrals ] ; improper " << std::endl;		
	os << ";    i      j      k      l   func   phase     kd      pn" << std::endl;
	
	for(ditr = p_mm_model->ImprDihedrals.begin(); ditr != p_mm_model->ImprDihedrals.end(); ditr++)
	{
		MMDihedral& dih = (MMDihedral&)(*ditr);
		
		std::vector<double> pk = dih.pk;
		std::vector<double> pn = dih.pn;
		std::vector<double> phase = dih.phase;

		HaAtom* aptr1 = (HaAtom*) dih.pt1;
		HaAtom* aptr2 = (HaAtom*) dih.pt2;
		HaAtom* aptr3 = (HaAtom*) dih.pt3;
		HaAtom* aptr4 = (HaAtom*) dih.pt4;

		if( at_idx_map.count(aptr1) == 0 ) continue;
		if( at_idx_map.count(aptr2) == 0 ) continue;
		if( at_idx_map.count(aptr3) == 0 ) continue;
		if( at_idx_map.count(aptr4) == 0 ) continue;

		int idx1 = at_idx_map[aptr1]+1;  // convert to 1-based index
		int idx2 = at_idx_map[aptr2]+1;  // convert to 1-based index
		int idx3 = at_idx_map[aptr3]+1;  // convert to 1-based index
		int idx4 = at_idx_map[aptr4]+1;  // convert to 1-based index

		std::string at_lbl_1 = aptr1->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_2 = aptr2->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_3 = aptr3->GetRef(HaAtom::ATOMREF_STD);
		std::string at_lbl_4 = aptr4->GetRef(HaAtom::ATOMREF_STD);

		for( i = 0; i < pk.size(); i++)
		{
			pk[i] *= 4.184;
			sprintf(buf,"%6d %6d %6d %6d   1  %14.4e%14.4e %2.0f  ; %s - %s - %s - %s ",idx1,idx2,idx3,idx4, phase[i], pk[i], pn[i], 
				at_lbl_1.c_str(), at_lbl_2.c_str(),at_lbl_3.c_str(),at_lbl_4.c_str());
			os << buf << std::endl;
		}	
	}
			
	os << "  " << std::endl;
	return TRUE;
}
