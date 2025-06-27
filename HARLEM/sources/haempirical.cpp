#include "hampi.h"
#include <stdio.h>
#include <memory>

#include "haatom.h"
#include "hamolset.h"
#include "hamolecule.h"
#include "haresdb.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "haempirical.h"
#include <time.h>
#include "canvas3d.h"
#include "math.h"


//temporal fix on Linux
#if !defined(_MSC_VER)
	#define SEEK_SET	0	/* Seek from beginning of file. from stdio */
#endif

HaEmpiricalMod::HaEmpiricalMod(MolSet* new_phost_mset):
HaCompMod(COMP_MOD_EMPIRICAL,new_phost_mset)
{
	SetStdParams();
}


HaEmpiricalMod::~HaEmpiricalMod()
{

}

int HaEmpiricalMod::SetStdParams()
{
	npt_dist = 0;
	sigma_constr = 0.5;
	weight_constraints = 10.0 ;
	pack_dist_com = 12.8;
	sigma_pack_dist_com = 5.3;
	weight_pack_distance = 0.6 ;
	pack_dist_axis = 10.7;
	sigma_pack_dist_axis= 5.2;
	pack_angle = 30.9;
	sigma_pack_angle = 16.4;
	weight_pack_angle = 10.0;
	dist_contct_com = 18.6 ;
	sigma_dist_contct_com = 7.32;
	num_contact_com = 2;
	weight_num_contact = 0.6;
	weight_vdw = 1.0;
	weight_bured = 1.0;
//	pack_dens = 14.0;
	pack_dens = 24.0;
	sigma_pack_dens = 1.0;
    weight_pack_dens = 0.6;
	sigma_sym = 0.5 ;
	weight_sym = 0.6 ;
	face_up_bound = 85.0 ;
	weight_sa = 0.6 ;
	dist_neighborhood = 15.0;
	lenght_factor = 1.0;
	com_flag = TRUE ;
	module_to_init_flag = TRUE;
	curr_state = 0;
	return TRUE;
}

int HaEmpiricalMod::Initialize()
{
	EstablishChains() ;
	LoadEmpParam();
	ResidueTypeList();
	CbettaSetUp(); //due to dummy atoms
	//LoadEmpConstrains();
	//LoadSolventAccessibleAtoms();
	//LineSegments();
	//FindCentralAxis(); //
	//Neighborhood();
	//InitCylinders(); // init cylinders toy
	//GetMaxDimension();

	module_to_init_flag = FALSE;
	return TRUE;
}
double
HaEmpiricalMod::ScoreEnergy()
{
	int output_yes = 1;
	if(module_to_init_flag)
	{
		Initialize();
	}

	char buf[256]; 
	if (com_flag)
	{
		CenterOfMass();   // commented using cylinders
		//FindAxes();
		//GetGeomCenterToy();
		com_flag = FALSE;
	}
	double score_ene = 0.0;
	double p_acs = 0.0 ;
	double p_bured= 0.0;
	double p_ene = 0.0;
double p_constr =0;
double p_packdist=0;
double p_packangle =0;
double p_contct =0;
double p_vdw =0;
double p_den =0;
double	p_sym = 0;
double p_pair =0;
	double p_rep = 0;
	double p_btc = 0;
	double pack_dens_calculated =0;
	double center_attract = 0.0;
	double p_helicepack = 0.0;  //added by Jose
	double p_toy = 0.0;
//clock_t tt1 = clock();
//	time_t t3=time(NULL);
 //  center_attract = PenaltyCentralAttract(); //commented by Jose
 //  p_bured= BuriedEnergyCyl();
//	dist_terminal = PenaltyLoop();
//   p_helicepack = PenaltyHelicePack();//PenaltyHelicePack();
//clock_t tt2 = clock();
//printf("%2.2f clock ticks between  tt2-tt1 CentralAttract \n", (double)(tt2-tt1)/CLOCKS_PER_SEC);

//printf(" CentralAttract %3.4f\n",  CentralAttract());

//	time_t t1=time(NULL);
//clock_t tt1 = clock();
//printf("clock_t tt1\n");
//10.3.3.p_den =  PenaltyDensity();
//clock_t tt2 = clock();
//printf("%2.2f clock ticks between  tt2-tt1 PenaltyDensity \n", (double)(tt2-tt1)/CLOCKS_PER_SEC);

//p_constr = PenaltyConstraints();
//	time_t t2=time(NULL);
//	printf("%d seconds elapsed between  t2-t1 score\n", t2-t1);
	
//p_packdist=  PenaltyPackDistance();
//	time_t t3=time(NULL);
//	printf("%d seconds elapsed between  t3-t2 score\n", t3-t2);


//  p_packangle = PenaltyPackAngle() ; //commented by Jose
//    p_rep = PenaltyRepulsionCyl();
//	p_toy = ToyEnergy();
//	p_toy = MinEnergy();
//    p_toy = LJState_ene();
//	p_toy = LJ_eneCylinder();	
// clock_t tt3 = clock();
//PrintLog("%f PackAngle + Repulsion + LJ_ene\n", (double)(tt3-tt1));///CLOCKS_PER_SEC);

//printf(" PenaltyPackAngle %3.4f\n",  PenaltyPackAngle());
	//time_t t4=time(NULL);
//	printf("%f seconds elapsed between  t4-t3 score\n", t4-t3);

//p_contct = PenaltyContact(); 
//	time_t t5=time(NULL);
//	printf("%d seconds elapsed between  t5-t4 score\n", t5-t4);


//p_vdw =  PenaltyVDW();

//		double p_vdw =  PenaltyVDW_Bured();
//	time_t t6=time(NULL);
//	printf("%d seconds elapsed between  t6-t5 score\n", t6-t5);

//	time_t t7=time(NULL);

//p_sym = PenaltySymmetry(); 
//	time_t t8=time(NULL);
//	printf("%d seconds elapsed between  t8-t1 score\n", t8-t1);

 
//	time_t t9=time(NULL);
//	clock_t tt2 = clock();
//	printf("%2.2f clock ticks between  tt2-tt1 \n", (double)(tt2-tt1)/CLOCKS_PER_SEC);
   //	printf("%d seconds elapsed between  t7-t6 score\n", t7-t6);

//p_acs = PenaltySolventAccessible();
//	p_btc = PenaltyBringToCenter();

//	time_t t10, t11;

//	if (pack_dens_calculated <= pack_dens+sigma_pack_dens)
//	{

// p_bured= PenaltyBured();

  p_pair = PenaltyPairwise();
//clock_t tt4 = clock();
//printf("%2.2f PenaltyPairwise\n", (double)(tt4-tt3)/CLOCKS_PER_SEC);

//     t10=time(NULL);

//	t11=time(NULL);
//		p_ene = p_den + p_acs + p_bured + p_sym + p_btc; 

//	}
//	else
//	{
//		double p_btc = PenaltyBringToCenter();
//		p_ene = p_btc + p_acs + p_sym;
//		p_ene = p_btc ;
//		PrintLog ("PPPp_ene %2.1f  \n", p_ene);
//	}
//	sprintf(buf,"Time 1= %d , 2= %d , 3= %d , 4= %d , 5= %d , 6= %d , 7= %d , 8= %d , 9= %d , 10= %d ",
//		 t2-t1, t3-t2,t4-t3,t5-t4,t6-t5,t7-t6,t8-t7,t9-t8,t10-t9,t11-t10);
//	component_file << buf << endl;

	score_ene = center_attract + p_packdist+ p_packangle + p_contct + p_vdw + p_bured +p_pair + p_sym + p_constr + p_helicepack + p_toy + p_rep;
 //   score_ene = p_vdw ;
	if(output_yes){
	std::fstream component_file;
	component_file.open("energy_components.dat", std::ios::out | std::ios::app);
//	sprintf(buf,"Constr %3.2f PDist %3.2f PAngl %3.2f Conct %3.2f Bured %3.2f VDW %3.4f Pair %3.2f PSym %3.2f PSAcc %3.2f ",
//	p_constr, p_packdist, p_packangle, p_contct, p_bured, p_vdw, p_pair, p_sym, p_acs);
//	sprintf(buf,"PAngl %14.9f Pair %14.9.2f  Phel %14.9f ",
//	p_packangle, p_pair, center_attract);
	sprintf(buf,"Ppack %14.9f Prep %14.9f PLJ %14.9f Pburied % 14.9f", center_attract, p_rep, p_toy, p_bured);
	component_file << buf << std::endl;
	component_file.close();
	}
//	sprintf(buf,"PACK DENSITY %4.3f BringToCenter %3.2f \n", pack_dens_calculated, p_btc); 
//	component_file << buf << endl; // Remove

//PrintLog("score_ene %4.3f \n", score_ene);
	com_flag = TRUE ;
//clock_t tt2 = clock();
//printf("%2.2f clock ticks between  tt2-tt1, p_pair %12.1f \n", (double)(tt2-tt1)/CLOCKS_PER_SEC, p_pair);
//printf("%2.2f TOTAL\n", (double)(tt4-tt1)/CLOCKS_PER_SEC);

	return score_ene; 
}

double
HaEmpiricalMod::GeometryScoreEnergy()
{
	if(module_to_init_flag)
	{
		Initialize();
	}
	double score_ene = 0.0;
	score_ene  =   PenaltyCentralAttract() + PenaltyRepulsionCyl();//+ PenaltyPackAngle();
//	score_ene = PenaltyHelicePack() + PenaltyRepulsion();//+PenaltyPairwise();
	return score_ene;
}
double HaEmpiricalMod::LennardJonesEnergy()
{
	if(module_to_init_flag)
	{
		Initialize();
	}
	double score_ene = 0.0;
	score_ene  =   PenaltyCentralAttract() + PenaltyRepulsion() + LJ_ene();//+ PenaltyPackAngle();
//	score_ene = PenaltyHelicePack() + PenaltyRepulsion();//+PenaltyPairwise();
	return score_ene;
}

/*double 
HaEmpiricalMod::PenaltyLoop()   // Function added by Jose
{
	double weight_loop = 0.5;
	HaAtom* aptr1;
	HaAtom* aptr2;
	int imol = 0;
	int nmol = pmset->GetNMol();
	int i = 0;
	int ires;
	
	for (imol=0; imol < nmol; imol++)
	{
		pMol = pmset->GetMolByIdx(imol);
		num_res = pMol -> GetNRes();
		ires = 0;
		HaMolecule::ResidueIterator ritr(pMol);
		for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
		{
			ires = res_ptr-> GetSerNo();
			if (ires == )
			aptr = (HaAtom*) atm_sc_array[i];
			aptrca1 = (HaAtom*) atm_ca_array[i];
			i++;
		}
		
	}

}
*/
double
HaEmpiricalMod::PenaltyConstraints()
{
	HaAtom* aptr1;
	HaAtom* aptr2;
	double constrains_ene = 0.0;
	double cur_constrains_ene = 0.0;
	double cur_low_bound = 0.0;
	double cur_up_bound = 0.0;

	MolSet* pmset = GetMolSet();	
	int i;
 	for(i =0; i< npt_dist; i++)
	{
		if(emp_dist[i] != 0.0 )
		{
			aptr1 = (HaAtom*) atm_dist[i];
			aptr2 = (HaAtom*) atm2_dist[i];
			if (aptr1 != NULL && aptr2 != NULL )
			{
				double dist = Vec3D::CalcDistance(aptr1,aptr2,ANGSTROM_U);
			//	PrintLog ("emp_dist[i] %4.2f , dist  %4.2f \n", emp_dist[i], dist);
			//	cur_low_bound=  emp_dist[i]- sigma_constr;
			//	cur_up_bound =  emp_dist[i] + sigma_constr;
				cur_constrains_ene= SqrPotential(dist,emp_dist[i],sigma_constr, weight_constraints);
			}
				constrains_ene += cur_constrains_ene;
		}
	}
//	PrintLog ("constrains_ene %4.2f \n", constrains_ene );
	return constrains_ene;
}


double
HaEmpiricalMod::PenaltyPackDistance()
{
//	fstream component_file;
//	component_file.open("energy_components.dat",ios::out | ios::app);
//  char buf[256];
	MolSet* pmset = GetMolSet();
	double packdistance_ene = 0.0;
	int ichain = 0;
	int ichain_old =0;
	double pack_dist_ene = 0.0;
	double cur_pack_dist_ene = 0.0;
   	int nmol = pmset->GetNMol();
	if (com_flag)
	{
		CenterOfMass();
	}
	for(  int imol =0; imol < nmol-1 ; imol++)
	{	
		ichain = chain_arr[imol];
		ichain_old= chain_arr[imol+1];
		if (ichain == ichain_old)
		{
			Vec3D* vec1 = &center_arr[imol];
			Vec3D* vec2 = &center_arr[imol+1];
			double dist = Vec3D::CalcDistance(vec1,vec2,ANGSTROM_U);
			cur_pack_dist_ene = SqrPotential(dist,pack_dist_com, sigma_pack_dist_com,weight_pack_distance);
			pack_dist_ene += cur_pack_dist_ene;	
//sprintf(buf,"PACK DISTANCE %3.2f imol1 %d imol2 %d \n", dist, imol, imol+1);
//	component_file << buf << endl;
		}
	}
//	PrintLog ("pack_dist_ene %4.2f  \n", pack_dist_ene);
//  component_file.close();
	return pack_dist_ene;
}

double
HaEmpiricalMod::PenaltyPackAngle()
{
//	fstream component_file;
//	component_file.open("energy_components.dat",ios::out | ios::ate);
//	char buf[256];
	MolSet* pmset = GetMolSet();
	int ichain = 0;
	int ichain_old =0;
	int nmol = pmset->GetNMol();
	double pack_angle_ene =0.0;
    double ux,uy,uz;
    double vx,vy,vz;
    double temp;
 	double angle = 0.0; 
	

	for(int imol =0; imol < nmol-1 ; imol++)
	{
		ichain = chain_arr[imol];
		ichain_old= chain_arr[imol+1];
		if (ichain != ichain_old)  // ichain != ichain_old but this is working
		{
			Vec3D vec1 = axis_arr[imol];
			Vec3D vec2 = axis_arr[imol+1];
			ux =  vec1[0];
			uy = vec1[1];
			uz = vec1[2];
			
			vx =  vec2[0];
			vy = vec2[1];
			vz = vec2[2];

			temp = (ux*vx + uy*vy + uz*vz) ;
			angle = (acos(temp))*RAD_TO_DEG; 
// easy to cut
//			PrintLog("#angle %2.3f\n", angle);
			if (topology_arr[imol] == 1 || topology_arr[imol+1] == 1) angle = 180.0 - angle ;
// easy to cut
			if (angle > 90.0) angle = 180.0 - angle;
//			PrintLog("#angle %2.3f\n", angle);
			angle = fabs(angle);
//			pack_angle_ene += (angle - pack_angle)*(angle - pack_angle)*weight_pack_angle;
			pack_angle_ene += SoftSqrWellPotential(angle, pack_angle, sigma_pack_angle, weight_pack_angle);
//			PrintLog("#angle %2.3f, pack_angle %2.3f, sigma_pack_angle %2.3f pack_angle_ene %2.3f\n", angle, pack_angle, sigma_pack_angle, pack_angle_ene);
		}
//		else
//			PrintLog("ichain != ichain_old\n");
	}

//	component_file.close();
	return pack_angle_ene;
}



int
HaEmpiricalMod::CalcPackAngleForceTorque(Vec3DValArray& torque_array)
{
	
//	fstream component_file;
//	component_file.open("energy_components.dat",ios::out | ios::ate);
//	char buf[256];
	MolSet* pmset = GetMolSet();
	int ichain = 0;
	int ichain_old =0;
	int nmol = pmset->GetNMol();
	double deriv_pack_angle;
    double ux,uy,uz;
    double vx,vy,vz;
    double temp;
 	double angle = 0.0; 
	Quaternion quat1;
	Quaternion quat2;
	Vec3D trans_v;
	HaMat_double rot_mat;
	rot_mat.newsize(3,3);
	Vec3DValArray temp_torque_array;
	temp_torque_array.resize(nmol);
	Vec3D temp_vec1;
	Vec3D temp_vec2;
	HaMolecule* pMol;
	int imol;
	NumVector<double> torque(3);
	temp_vec1[0] =0;
	temp_vec1[1] =0;
	temp_vec1[2] =0;
	static HaVec_double vec(3);
	double qang, qx, qy, qz;

	for(imol =0; imol < nmol ; imol++)
	{
		 temp_torque_array[imol] = temp_vec1;
	}

	for(imol =0; imol < nmol-1 ; imol++)
	{
		ichain = chain_arr[imol];
		ichain_old= chain_arr[imol+1];

		if (ichain == ichain_old)
		{
			temp_vec1 = temp_torque_array[imol];
			temp_vec2 = temp_torque_array[imol+1];
			Vec3D vec1 = axis_arr[imol];
			Vec3D vec2 = axis_arr[imol+1];
			ux =  vec1[0];
			uy = vec1[1];
			uz = vec1[2];
			
			vx =  vec2[0];
			vy = vec2[1];
			vz = vec2[2];

			temp = (ux*vx + uy*vy + uz*vz) ;
 			angle = (acos(temp))*RAD_TO_DEG; 
			double cos_val = temp;
			double sin_val = 1-temp;

//			rot_mat.SetVal(1,1,cos_val); rot_mat.SetVal(1,2,-sin_val);rot_mat.SetVal(1,3,0);
//			rot_mat.SetVal(2,1,sin_val); rot_mat.SetVal(2,2,cos_val);rot_mat.SetVal(2,3,0);
//			rot_mat.SetVal(3,1,0 ); rot_mat.SetVal(1,3,0); rot_mat.SetVal(1,3,1);

			if (topology_arr[imol] == 1 || topology_arr[imol+1] == 1) angle = 180.0 - angle ;
			deriv_pack_angle = 2*(angle - pack_angle)*weight_pack_angle * 100.0;
			quat1.QuaternionFromAxis(1, vec1[0], vec1[1], vec1[2]);
			quat2.QuaternionFromAxis(1, vec2[0], vec2[1], vec2[2]);
			quat1.operator *=(quat2);
			quat1.GetQuaternion(qang, qx, qy, qz);
			vec1[0]= qx *deriv_pack_angle;
			vec1[1]= qx *deriv_pack_angle;
			vec1[2]= qz *deriv_pack_angle;
			PrintLog("vec1[0], vec1[1], vec1[2] %2.1f, %2.1f, %2.1f \n",vec1[0], vec1[1], vec1[2]);
			temp_vec1[0] += vec1[0];
			temp_vec1[1] += vec1[1];
			temp_vec1[2] += vec1[2];
			temp_vec2[0] -= vec1[0];;
			temp_vec2[1] -= vec1[1];
			temp_vec2[2] -= vec1[2];

			temp_torque_array[imol].SetCoordFrom(temp_vec1);
			temp_torque_array[imol+1].SetCoordFrom(temp_vec2);
			//pack_angle_ene += SoftSqrWellPotential(angle, pack_angle, sigma_pack_angle, weight_pack_angle);
//		PrintLog("#angle %2.3f, pack_angle %2.3f, sigma_pack_angle %2.3f pack_angle_ene %2.3f\n", angle, pack_angle, sigma_pack_angle, pack_angle_ene);
		}
	}



	for(imol =0; imol < nmol ; imol++)
	{
		pMol = pmset->GetMolByIdx(imol);
		pMol->GetQuaternionTrans(quat1,trans_v);

		Quaternion::QuaternionToRotMat(quat1, rot_mat);
		temp_vec1 = temp_torque_array[imol];
		vec[0] = temp_vec1[0];
		vec[1] = temp_vec1[1];
		vec[2] = temp_vec1[2];
		torque = matmult(rot_mat, vec);
		temp_vec2[0] = torque(1);
		temp_vec2[1] = torque(2);
		temp_vec2[2] = torque(3);
		torque_array[imol].SetCoordFrom(temp_vec2);
	}

//	component_file.close();
	return True;
}


double
HaEmpiricalMod::PenaltyContact()
{
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	int n_contact ;
	double cur_contact_ene =0.0;
	double contact_ene =0.0;
	if (com_flag)
	{
		CenterOfMass();
	}

	double up_com = dist_contct_com - 1.5*sigma_dist_contct_com;
	for(  int imol =0; imol < nmol; imol++)
	{	
		Vec3D* vec1 = &center_arr[imol];
		n_contact = 0;
		for(  int jmol =0; jmol < nmol; jmol++)
		{
			if (imol != jmol)
			{
				Vec3D* vec2 = &center_arr[jmol];
				double dist = Vec3D::CalcDistance(vec1,vec2,ANGSTROM_U);
				if (dist < up_com)
				{
					n_contact ++;
				}
			}
		}
		if (n_contact >= num_contact_com)
		{
			cur_contact_ene = 0.0;
		}
		else
		{
			cur_contact_ene = weight_num_contact*(num_contact_com - n_contact);
		}
		contact_ene += cur_contact_ene; 
	}
	return contact_ene;
}


double
HaEmpiricalMod::PenaltyBured()
{
	int output_yes = 1;
	std::fstream component_file;
	component_file.open("energy_components.dat", std::ios::out | std::ios::ate);
	char buf[256];


	MolSet* pmset = GetMolSet();
	wxString name, name1 ;
	HaAtom* aptr;
	HaAtom* aptr1;
	HaAtom* aptr2;
	HaResidue* res_ptr;
	HaResidue* res_ptr1;
	HaMolecule* pMol;
	HaMolecule* pMol1;
	int exposed_flag ;
	double cur_bured_ene = 0.0;
	double cur_la_value = 0.0;
	int ires, ires1, num_res, num_res1;
	int found_flag = FALSE; // residue is found/not found in solventaccesible.dat file
	double sigma1, sigma2;
	double cur_dist_bured = 0.0;
	double dist_bured = 0.0;
	int raw,raw_attract;
	int i = 0;
	int j = 0;
	int nmol = pmset->GetNMol();
	double pairwise =0;
	double vdw_ene = 0.0;
	int	ncount_bured	;
	int imol;

	if (residue_arr[0] == " ") // break if the LA scale file was not read
	{
		return cur_bured_ene;
	}
	else
	{
		for( imol =0; imol < nmol ; imol++)
		{
			pMol = pmset->GetMolByIdx(imol);
			num_res = pMol -> GetNRes();
			std::string mol_name = pMol -> GetRef();
			ires = 0;
			HaMolecule::ResidueIterator ritr(pMol);
			for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
			{	
				found_flag = FALSE;
				ires = res_ptr-> GetSerNo();
                if (marker_res_sa[i])	found_flag = TRUE ;
				aptr = (HaAtom*) atm_sc_array[i];
				name = aptr ->GetName();
				raw = mol_res_correspond_la[i];
				raw_attract = mol_res_correspond[i];
				i++;
				sigma1 = aptr -> radius ;
				exposed_flag = TRUE;
				
				if(found_flag) continue;
				
				ncount_bured = 0;
				int jmol;
				for(  jmol = 0; jmol < nmol ; jmol++)
				{
					j = 0;
					if (imol != jmol && CheckNeighbor(imol, jmol) == 1.0 )
					{
						int number_one = first_res_mol[jmol] ;
						pMol1 = pmset->GetMolByIdx(jmol);
						std::string mol_name1 = pMol1 -> GetRef();
						HaMolecule::ResidueIterator ritr1(pMol1);
						for(res_ptr1 = ritr1.GetFirstRes(); res_ptr1; res_ptr1 = ritr1.GetNextRes())
						{	
							ires1 = res_ptr1 -> GetSerNo();
							aptr1 = (HaAtom*) atm_sc_array[j+ number_one ];
							cur_la_value = 0.0;
							int	col = mol_res_correspond[number_one + j];
							
							j++;
							dist_bured = Vec3D::CalcDistance(aptr,aptr1,ANGSTROM_U);
							sigma2 = aptr1 -> radius ;
							cur_dist_bured = (sigma1 + sigma2)* 1.2;
							
							
							//if (dist_bured <= cur_dist_bured )
							if (dist_bured <= 7.0 )
							{
								exposed_flag = FALSE;
								ncount_bured ++;								
								//PrintLog("Bured@@: $%s$%s%d-$%s$%s%d, dist %2.1f,dist_bured %2.1f Ene %3.3f la_value[%d] %3.3f\n", mol_name.c_str() , res_ptr->GetName(), res_ptr->GetSerNo(),mol_name1.c_str() , res_ptr1->GetName(), res_ptr1->GetSerNo(),dist_bured, cur_dist_bured,cur_bured_ene,raw,la_value[raw]);
								//									cur_bured_ene += cur_la_value ;
								
								
								/*								int bin;
								if(dist_bured <2.0) 
								{
								bin = 0;
								}
								if(dist_bured >= 2.0 && dist_bured <8.0)	bin = (int) ((dist_bured-2.0)/0.5) +1;
								if(dist_bured >= 8.0 && dist_bured <15.0)	
								{
								bin = (int) (dist_bured-8.0)+13;
								}
								//int bin = (int) ( ((dist-2.0)/delta_bin) + 0.5) ;
								double pairwise_ene;
								pairwise_ene += pairwise_energy_arr.GetValue(raw_attract, col,bin);
								*/								
								//								pairwise -= pairwise_energy_arr.GetVal_idx0(raw_attract,col)/(dist_bured+0.000001);
								//								char buf[256];
								//PrintLog("Pair@@: $%s$%s%d-$%s$%s%d, pair %3.3f\n", mol_name.c_str() , res_ptr->GetName(), res_ptr->GetSerNo(),mol_name1.c_str() , res_ptr1->GetName(), res_ptr1->GetSerNo(),pairwise_ene );
								//PrintLog("Bured@@: $%s$%s%d-$%s$%s%d, dist %2.1f,dist_bured %2.1f Ene %3.3f la_value[%d] %3.3f\n", mol_name.c_str() , res_ptr->GetName(), res_ptr->GetSerNo(),mol_name1.c_str() , res_ptr1->GetName(), res_ptr1->GetSerNo(),dist_bured, cur_dist_bured,cur_bured_ene,raw,la_value[raw]);
								//	break;
							}
							////////	This part for taking pairwise interaction on pore lining residues
							////////									else if (found_flag && dist_bured>= sigma && pairwise_energy_arr.GetVal_idx0(0,0) != -1000000.0)
							//									{
							//										res_ptr = aptr-> GetHostRes();
							//										res_ptr1 = aptr1-> GetHostRes();
							//										double pairwise = pairwise_energy_arr.GetVal_idx0(raw_attract,col);
							//										atract_vdw_ene += pairwise*exp(cur_dist_bured-dist_bured);
							//	PrintLog("@@POREresidue@@: residue %s%d,residue1 %s%d pairwise %4.3f, raw %d,col %d \n",  res_ptr->GetName(), res_ptr->GetSerNo(),res_ptr1->GetName(), res_ptr1->GetSerNo(), pairwise,raw_attract,col);
							////////									}
							
							
						}
					}
				}				
				//				if(exposed_flag)
				//				{	
				//					cur_la_value = - la_value[raw] * la_weight_value[raw];
				//					cur_bured_ene += cur_la_value  ;
				//					//PrintLog("More Dist: residue %s%d, la_value %3.3f \n",  res_ptr->GetName(), res_ptr->GetSerNo(), la_value[raw]);
				//					//PrintLog("@@BuredExposed@@: residue %s%d, dist %3.3f, cur_bured_ene %3.3f \n",  res_ptr->GetName(), res_ptr->GetSerNo(),dist_bured, cur_bured_ene);
				//				}
				// PrintLog("@@@@: residue %s%d,  \n",  res_ptr->GetName(), res_ptr->GetSerNo());
if (output_yes){
				sprintf(buf,"residue %s%d Ncontacs %d",res_ptr->GetName(), res_ptr->GetSerNo(),ncount_bured); // Remove
				component_file << buf<< std::endl ; // Remove
			}

				cur_la_value = - la_value[raw] * la_weight_value[raw];
				if (ncount_bured < 1) 
					cur_bured_ene += cur_la_value ;
				
				if (ncount_bured >=1 && ncount_bured < 2) 
					cur_bured_ene += cur_la_value * 0.7   ;
				
				if (ncount_bured >= 2 && ncount_bured < 4) 
					cur_bured_ene += cur_la_value * 0.5 ;
				
				if (ncount_bured >= 4 && ncount_bured < 6) 
					cur_bured_ene += cur_la_value * 0.3 ;
				
				if (ncount_bured >= 6) 
					cur_bured_ene += cur_la_value * 0.0 ;
			}			
		}
		
	}
	//	PrintLog ("&&&cur_bured_ene %5.3f \n", cur_bured_ene);
	//	PrintLog("@@VDW@@ atract_vdw_ene: %3.3f \n", atract_vdw_ene);
if (output_yes){
	component_file.close();
}
	return cur_bured_ene ;
}

/*
double
HaEmpiricalMod::PenaltyPairwise()
{
	MolSet* pmset = GetMolSet();
	wxString name, name1 ;
	HaAtom* aptr;
	HaAtom* aptr1;
	HaAtom* aptr2;
	HaAtom* aptrca;
	HaResidue* res_ptr;
	HaResidue* res_ptr1;
	HaMolecule* pMol;
	HaMolecule* pMol1;
	double cur_bured_ene = 0.0;
	double cur_la_value = 0.0;
	int ires, ires1, num_res, num_res1;
	int found_flag = FALSE; // residue is found/not found in solventaccesible.dat file
	int found_flag1 = FALSE; // residue is found/not found in solventaccesible.dat file
	double sigma1, sigma2, sigma_ca;
	double cur_dist = 0.0;
	double dist = 0.0;
	double distca;
	int raw,col, raw_la;
	int i = 0;
	int j = 0;
	int nmol = pmset->GetNMol();
	double pairwise_ene =0;
	double vdw_ene = 0.0;
	double delta_bin;
	int bin;
	double overcount_coef =1.0;
double depth = 20/(double) nmol;
//	double depth = 1.0;

	int ncount_bured;
	double ene;

//double depth = 20.0;
	if (residue_arr[0] == " ") // break if the Pairwise scale file was not read
	{
		return pairwise_ene= 0;
	}
	else
	{
		for( int imol =0; imol < nmol ; imol++)
		{
			pMol = pmset->GetMolByIdx(imol);
			num_res = pMol -> GetNRes();
			std::string mol_name = pMol -> GetRef();
			ires = 0;
			HaMolecule::ResidueIterator ritr(pMol);
			for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
			{	
				found_flag = FALSE;
				ires = res_ptr-> GetSerNo();
				aptr = (HaAtom*) atm_sc_array[i];
				raw = mol_res_correspond[i];
				raw_la = mol_res_correspond_la[i];
                if (marker_res_sa[i])
				{
					found_flag = TRUE ;
				}
				sigma1 = aptr -> radius ;
				name = res_ptr ->GetName();

				i++;				
//				if(found_flag) continue;
				ncount_bured = 0;

				for( int jmol = 0; jmol < nmol ; jmol++)
				{
					j = 0;
					//if (imol != jmol)
					if (imol != jmol && CheckNeighbor(imol, jmol) == 1.0 )
					{

						int number_one = first_res_mol[jmol] ;
						pMol1 = pmset -> GetMolByIdx(jmol);
						std::string mol_name1 = pMol1 -> GetRef();
						HaMolecule::ResidueIterator ritr1(pMol1);
						for(res_ptr1 = ritr1.GetFirstRes(); res_ptr1; res_ptr1 = ritr1.GetNextRes())
						{	
							found_flag1 = FALSE;
							ires1 = res_ptr1 -> GetSerNo();
							aptr1 = (HaAtom*) atm_sc_array[j+ number_one ];
							aptrca= (HaAtom*) atm_ca_array[j+ number_one ];
							col = mol_res_correspond[number_one + j];
							if (marker_res_sa[number_one + j])
							{
								found_flag1 = TRUE ;
							}
							j++;
							dist = Vec3D::CalcDistance(aptr,aptr1,ANGSTROM_U);
							sigma2 = aptr1 -> radius;
							sigma_ca = aptrca -> radius;

							if (dist <= 7.0 )
							{
								ncount_bured ++;								
							}

							double sigma = (sigma1 + sigma2)* 0.65;

							distca = Vec3D::CalcDistance(aptr,aptrca,ANGSTROM_U);
							double dsigma_ca = (sigma_ca+sigma1)*0.8;
							if( distca < dsigma_ca )
							{
								pairwise_ene += (dsigma_ca - distca)*(dsigma_ca - distca);
								//									PrintLog("sigma_ca - distca %2.1f\n", (dsigma_ca - distca)*(dsigma_ca - distca));
							}

							name1 = res_ptr1 ->GetName();
							
							if(name == name1) overcount_coef = 0.5;
							else overcount_coef = 1.0;

							if ( dist <= 15.0)
							{
								if(dist < sigma ) pairwise_ene += (sigma-dist)*(sigma-dist);
//PrintLog("sigma - dist %2.1f\n", (sigma - dist)*(sigma - dist));
								
							//	if(dist < 3.2 ) pairwise_ene += dist*dist;
								if(dist <2.0) 
									bin = 0;
								
								if(dist >= 2.0 && dist <8.0)
									bin = (int) ((dist-2.0)/0.5) +1;

								if(dist >= 8.0 && dist <15.0)	
									bin = (int) (dist-8.0)+13;

								if (!found_flag && !found_flag1)
								{
									ene = pairwise_energy_arr.GetValue(raw, col,bin)*overcount_coef* depth;
								}
								if (found_flag || found_flag1)
								{
//PrintLog("Found %d %d pairwise_ene=%2.6f %s%d-%s%d dist %2.3f bin %d\n", found_flag, found_flag1, pairwise_energy_arr.GetValue(raw, col,bin), res_ptr->GetName(), res_ptr->GetSerNo(), res_ptr1->GetName(), res_ptr1->GetSerNo(), dist, bin);
//									ene = 0.0;
									ene = pairwise_energy_arr_sa.GetValue(raw, col,bin)*overcount_coef* depth;
								}

								pairwise_ene += ene;

//								name  = res_ptr -> GetName();
//								name1 = res_ptr1->GetName();
//								if (name == "GLY" && name1 == "VAL")
//PrintLog("%d %d pairwise_ene=%2.6f %s%d-%s%d dist %2.3f bin %d\n", imol, jmol, pairwise_energy_arr.GetValue(raw, col,bin), res_ptr->GetName(), res_ptr->GetSerNo(), res_ptr1->GetName(), res_ptr1->GetSerNo(), dist, bin);
							}
						}
					}
				}
//				if (!found_flag)
//				{
//	PrintLog("%d %d  %s%d \n", imol, jmol,  res_ptr->GetName(), res_ptr->GetSerNo());
			
//					cur_la_value = - la_value[raw_la] * la_weight_value[raw_la];
//					if (ncount_bured < 1) 
//						cur_bured_ene += cur_la_value ;
//					
//					if (ncount_bured >=1 && ncount_bured < 2) 
//						cur_bured_ene += cur_la_value * 0.7   ;
//					
//					if (ncount_bured >= 2 && ncount_bured < 4) 
//						cur_bured_ene += cur_la_value * 0.5 ;
//					
//					if (ncount_bured >= 4 && ncount_bured < 6) 
//						cur_bured_ene += cur_la_value * 0.3 ;
//					
//					if (ncount_bured >= 6) 
//						cur_bured_ene += cur_la_value * 0.0 ;
//				}
			}
		}
		
	}
	
//PrintLog("pairwise_ene=%2.6f cur_bured_ene %2.6f\n", pairwise_ene , cur_bured_ene);
	return pairwise_ene+ cur_bured_ene;
}
*/


double
HaEmpiricalMod::PenaltyPairwise()
{
	MolSet* pmset = GetMolSet();  // assign the Van der Waals radius of each CA and SC residue on HaMolMech
	wxString name, name1 ;
	HaAtom* aptr;
	HaAtom* aptr1;
	HaAtom* aptr2;
	HaAtom* aptrca1;
	HaAtom* aptrca2;
	HaResidue* res_ptr;
	HaResidue* res_ptr1;
	HaMolecule* pMol;
	HaMolecule* pMol1;
	std::fstream pairwi_file;
	pairwi_file.open("pairwi_file.dat", std::ios::out | std::ios::app);
	char buf[256];

	double cur_bured_ene = 0.0;
	double cur_la_value = 0.0;
	int ires, ires1, num_res, num_res1;
	int sol_acces_flag = FALSE; // residue is found/not found in solventaccesible.dat file
	int sol_acces_flag1 = FALSE; // residue is found/not found in solventaccesible.dat file
	double sigma1, sigma2, sigma_ca1,sigma_ca2;
	double cur_dist = 0.0;
	double dist = 0.0;
	double dist_ca_sc,dist_ca_ca;
	int raw,col, raw_la;
	int i = 0;
	int j = 0;
	int nmol = pmset->GetNMol();
	double pairwise_ene =0.0;
	double vdw_ene = 0.0;
	double delta_bin;
	int bin;
	double overcount_coef =1.0;
//double depth = 20/(double) nmol;
	double depth = 63.69427;
 
	double ene;
	double surf_tension = 0.005;
	double solv_rad = 1.4;

	int natom = 2*atm_ca_array.size();
	HaMat_double coord;
	coord.newsize(3,natom);
	static HaVec_double radia;
	radia.newsize(natom);
	static HaVec_double weight_hydrophil;
	weight_hydrophil.newsize(natom);
	
	int imol;
	int jmol;
	int icount;
	for( imol =0; imol < nmol ; imol++)
	{
		pMol = pmset->GetMolByIdx(imol);
		num_res = pMol -> GetNRes();
		std::string mol_name = pMol -> GetRef();
		ires = 0;
		HaMolecule::ResidueIterator ritr(pMol);
		for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
		{	
			sol_acces_flag = FALSE; //control of pairwise potential jose
			ires = res_ptr-> GetSerNo();
			aptr = (HaAtom*) atm_sc_array[i];
			aptrca1 = (HaAtom*) atm_ca_array[i];
			raw = mol_res_correspond[i];
			raw_la = mol_res_correspond_la[i];
			if (marker_res_sa[i])
			{
				sol_acces_flag = TRUE ;
			}
			sigma1 = aptr -> radius ;
			sigma_ca1 = aptrca1 -> radius ;
			
			coord.SetVal_idx0( 0,2*i, aptr-> GetX());
			coord.SetVal_idx0(1,2*i, aptr-> GetY());
			coord.SetVal_idx0(2,2*i, aptr-> GetZ());
			radia(2*i+1)   = sigma1;
			//if (la_value[raw_la] < 0) weight_hydrophil(2*i+1) = 1.0;
			//else weight_hydrophil(2*i+1) = 0.0;
			weight_hydrophil(2*i+1) = la_value[raw_la];
			if (sol_acces_flag) weight_hydrophil(2*i+1) = 0.0;

			if (nmol <=12)
			{
				coord.SetVal_idx0(0,2*i+1, aptrca1-> GetX());
				coord.SetVal_idx0(1,2*i+1, aptrca1-> GetY());
				coord.SetVal_idx0(2,2* i+1, aptrca1-> GetZ()); 
				radia(2*i+2)   = sigma_ca1;
				//if (la_value[raw_la] < 0) weight_hydrophil(2*i+2) = 1.0;
				//else weight_hydrophil(2*i+2) = 0.0;
				weight_hydrophil(2*i+2) = 0.0;
				if (sol_acces_flag) weight_hydrophil(2*i+2) = 0.0;
			}

			name = res_ptr ->GetName();
			
			i++;				
			icount =i;			

			for( jmol = imol; jmol < nmol ; jmol++)
			{
				j = 0;
				if (imol != jmol && CheckNeighbor(imol, jmol) == 1.0 )
				{
					
					int number_shift = first_res_mol[jmol] ;
					pMol1 = pmset -> GetMolByIdx(jmol);
					std::string mol_name1 = pMol1 -> GetRef();
					HaMolecule::ResidueIterator ritr1(pMol1);
					for(res_ptr1 = ritr1.GetFirstRes(); res_ptr1; res_ptr1 = ritr1.GetNextRes())
					{	
						sol_acces_flag1 = FALSE; //control of pairwise potential jose
						ires1 = res_ptr1 -> GetSerNo();
						aptr1 = (HaAtom*) atm_sc_array[j+ number_shift ];
						aptrca2= (HaAtom*) atm_ca_array[j+ number_shift ];
						col = mol_res_correspond[number_shift + j];
						if (marker_res_sa[number_shift + j])
						{
							sol_acces_flag1 = TRUE ;
						}
						j++;
						dist = Vec3D::CalcDistance(aptr,aptr1,ANGSTROM_U);
						sigma2 = aptr1 -> radius;
						sigma_ca2 = aptrca2 -> radius;
						
						double sigma = (sigma1 + sigma2)* 0.5;
						
						dist_ca_sc = Vec3D::CalcDistance(aptr,aptrca2,ANGSTROM_U);
						double sigma_ca_sc = (sigma_ca2+sigma1)*0.9;
						
						name1 = res_ptr1 ->GetName();
						
						if( dist_ca_sc < sigma_ca_sc )
						{
//						PrintLog("%s-%s distca %2.2f dsigma_ca %2.2f vdw %2.3f \n",name, name1, distca, dsigma_ca,(dsigma_ca - distca)*(dsigma_ca - distca) );
							vdw_ene += (sigma_ca_sc - dist_ca_sc)*(sigma_ca_sc - dist_ca_sc)*50.0;
						}

						dist_ca_ca = Vec3D::CalcDistance(aptrca1,aptrca2,ANGSTROM_U);
						double sigma_ca_ca = (sigma_ca1+sigma_ca2);
						if( dist_ca_ca < sigma_ca_ca )
						{
//						PrintLog("%s-%s distca %2.2f dsigma_ca %2.2f vdw %2.3f \n",name, name1, distca, dsigma_ca,(dsigma_ca - distca)*(dsigma_ca - distca) );
							vdw_ene += (sigma_ca_ca - dist_ca_ca)*(sigma_ca_ca - dist_ca_ca)*50.0;
						}
						
						if ( dist < 15.0) //modified by jose <=
						{
							if(dist < sigma )
							{
								vdw_ene += (sigma-dist)*(sigma-dist)*1.0; //*1.0
							}

							
							//PrintLog("sigma - dist %2.1f\n", (sigma - dist)*(sigma - dist));
							
							//	if(dist < 3.2 ) pairwise_ene += dist*dist;
							if(dist <2.0) 
								bin = 0;
							
							if(dist >= 2.0 && dist <8.0)
								bin = (int) ((dist-2.0)/0.5) +1;
							
							if(dist >= 8.0 && dist <15.0)	
								bin = (int) (dist-8.0)+13;
							
							if (!sol_acces_flag && !sol_acces_flag1)
							{
								ene = pairwise_energy_arr.GetValue(raw, col,bin)*overcount_coef* depth;
							}
							if (sol_acces_flag || sol_acces_flag1)
							{
							PrintLog("Found %d %d pairwise_ene=%2.6f %s%d-%s%d dist %2.3f bin %d\n", sol_acces_flag, sol_acces_flag1, pairwise_energy_arr_sa.GetValue(raw, col,bin), res_ptr->GetName(), res_ptr->GetSerNo(), res_ptr1->GetName(), res_ptr1->GetSerNo(), dist, bin);
								//									ene = 0.0;
								ene = pairwise_energy_arr_sa.GetValue(raw, col,bin)*overcount_coef* depth;
							}
							
							pairwise_ene += ene;
							
							
							//								name  = res_ptr -> GetName();
							//								name1 = res_ptr1->GetName();
							//								if (name == "GLY" && name1 == "VAL")
				//			PrintLog("%d %d pairwise_ene=%2.6f %s%d-%s%d dist %2.3f bin %d\n", imol, jmol, pairwise_energy_arr.GetValue(raw, col,bin), res_ptr->GetName(), res_ptr->GetSerNo(), res_ptr1->GetName(), res_ptr1->GetSerNo(), dist, bin);
						}
					}
				}
			}

		}
	}
		
//	if (pairwise_ene != 0.0 )
//	{
		if (nmol >12)
		{
			CalculateCoarseGrainedBackbone();
			natom	= rad_cgbb_array.size();
			for( i = 0; i < natom; i++)
			{
				Vec3D* ptr = &atm_cgbb_array[i];
//PrintLog("SSSS %d %d %2.3f %2.3f %2.3f\n",i, natom, ptr-> GetX(), ptr-> GetY(),ptr-> GetZ());
				coord.SetVal_idx0(0, i+icount, ptr-> GetX());
				coord.SetVal_idx0(1, i+icount, ptr-> GetY());
				coord.SetVal_idx0(2, i+icount, ptr-> GetZ()); 
				radia(i+icount +1)   = rad_cgbb_array(i+1);

				weight_hydrophil(i+icount +1) = 0.0;
			}
		}
//PrintLog("SSSS %d %2.3f %2.3f %2.3f\n", natom, ptr-> GetX(), ptr-> GetY(),ptr-> GetZ());
		
//PrintLog("coordPE %12.3f %12.3f %12.3f \n", coord.GetVal(1,1),coord.GetVal(2,1),coord.GetVal(3,1) );
//clock_t tt1 = clock();
/*

		int switchvalue = 0;
		HaSurface surf;
		surf.CalcMolSurfAlpha(switchvalue, solv_rad, coord, radia);
//clock_t tt2 = clock();
//printf("%2.5f clock ticks between  tt2-tt1 PAIRWISE \n", (double)(tt2-tt1)/CLOCKS_PER_SEC);
//PrintLog("surface_alpha_total =%2.6f \n", surf.surface_alpha_total);
//		double ca_area =0;
//		double cb_area =0;
		for(i=0; i < natom; i++)
		{
//			if (i%2 == 0) cb_area += surf.surface_alpha(i+1);
//			else ca_area += surf.surface_alpha(i+1);
			cur_bured_ene += surf.surface_alpha(i+1) * weight_hydrophil(i+1);
//			PrintLog("weight_hydrophil(%d)=%2.6f\n", i+1, weight_hydrophil(i+1));	
		}  
		cur_bured_ene *= surf_tension;
//		PrintLog(" ca_area =%2.6f cb_area =%2.6f \n", ca_area, cb_area);	
	PrintLog("cur_bured_ene =%2.6f pairwise_ene=%2.6f, vdw_ene =%2.6f\n",cur_bured_ene, pairwise_ene, vdw_ene);	

*/     
//PrintLog("pairwise_ene=%2.6f cur_bured_ene %2.6f vdw_ene %2.6f\n", pairwise_ene , cur_bured_ene, vdw_ene);

//	}
// end comments
	sprintf(buf, "%3.2f %3.2f %3.2f", pairwise_ene , cur_bured_ene, vdw_ene);
	pairwi_file << buf << std::endl;
//	pairwi_file.close();
	return pairwise_ene;// + cur_bured_ene + vdw_ene;
}






double
HaEmpiricalMod::PenaltyVDW_Bured()
{
	MolSet* pmset = GetMolSet();
	wxString name, name1 ;
	HaAtom* aptr;
	HaAtom* aptr1;
	HaAtom* aptr2;
	HaMolecule* pMol;
	HaMolecule* pMol1;
	int exposed_flag;
	int found_flag;
	double vdw_ene= 0.0;
	double cur_bured_ene = 0.0;
	double cur_la_value;
	double repulse_vdw_ene =0.0;
	double atract_vdw_ene =0.0;
	int nmol = pmset -> GetNMol();
	double cur_dist_vdw;
	double sigma, sigma1, sigma2;
	int i, j, raw, col, raw1;
	double dist;
	i=0;
	int abort_bured = FALSE;
	if (residue_arr[0] == " ") // break if the LA scale file was not read
	{
		abort_bured = TRUE;
	}

	for( int imol =0; imol < nmol ; imol++)
	{
		pMol = pmset->GetMolByIdx(imol);
		AtomIteratorMolecule aitr(pMol);
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			name = aptr -> GetName();
			raw = mol_res_correspond[i];
			raw1 = mol_res_correspond_la[i];
			i++;
			if (name == "CB" || name == "CA")
			{
				sigma1 = aptr -> radius ;
				exposed_flag = TRUE;
				found_flag = FALSE;
				for(j=0 ; j < n_sa_atoms ; j++)
				{
					aptr2 = (HaAtom*) atm_solacces[j];
					if (aptr2 != NULL)
					{
						if( aptr2 == aptr)
						{
							found_flag= TRUE;
							break;
						}
					}
				}
				for( int jmol = 0; jmol < nmol ; jmol++)
				{	
					j=0;
					if (imol != jmol && CheckNeighbor(imol, jmol) == 1.0 )
					{
						int number_one = first_res_mol[jmol] ;
						pMol1 = pmset->GetMolByIdx(jmol);
						AtomIteratorMolecule aitr1(pMol1);
						for(aptr1 = aitr1.GetFirstAtom(); aptr1; aptr1 = aitr1.GetNextAtom())
						{
							dist = Vec3D::CalcDistance(aptr,aptr1,ANGSTROM_U);
							name1 = aptr1 ->GetName();
							sigma2 = aptr1 -> radius ;
							sigma = sigma1 + sigma2;
							if (jmol > imol)
							{
								col = mol_res_correspond[number_one + j];
								j++;
								if (name1 == "CB" || name1 == "CA")
								{
									cur_dist_vdw = sigma*0.56125;
									if (dist < cur_dist_vdw)
									{
										repulse_vdw_ene +=	(cur_dist_vdw*cur_dist_vdw- dist*dist)*(cur_dist_vdw*cur_dist_vdw- dist*dist);
									}
//									else if (dist < sigma && name == "CB" && name1 == "CB" && pairwise_energy_arr.GetVal_idx0(0,0) != -1000000.0)
//									{
//										double pairwise = pairwise_energy_arr.GetVal_idx0(raw,col);
//										if (pairwise != 0.0) atract_vdw_ene += pairwise/dist;
//									}
								}
							}
							//Bured part
							if (dist <= sigma && !found_flag && name == "CB" && name1 == "CB" && !abort_bured)
							{
								cur_la_value = la_value[raw1] * la_weight_value[raw1];
								cur_bured_ene += cur_la_value * weight_bured;
								exposed_flag = FALSE;
							}
						}
					}
				}
				if(exposed_flag && !found_flag && name == "CB" && name1 == "CB" && !abort_bured)
				{	
					cur_la_value = - la_value[raw1] * la_weight_value[raw1];
					cur_bured_ene += cur_la_value * weight_bured ;
				}
			}
		}
	}
//	PrintLog ("^^^^^cur_vdw_ene %5.3f \n", (repulse_vdw_ene + atract_vdw_ene) * weight_vdw);
//	PrintLog ("^^^^^cur_Bured_ene %5.3f \n", cur_bured_ene);
	return vdw_ene = (repulse_vdw_ene + atract_vdw_ene) * weight_vdw + cur_bured_ene;
}

double HaEmpiricalMod::PenaltyDensity()
{
	double pack_dens_ene = 0.0;
	double svolume= 0.0;
	pack_dens_calculated = 0.0;
	int surf_type = 1;
	double solv_rad = 1.4;
	MolSet* pmset = GetMolSet();
    int natom = pmset -> GetNAtoms();

	HaMat_float coord1(3,natom);
	HaVec_float radii1(natom);
	HaMat_double coord;
	coord.newsize(3,natom);
	HaVec_double radii;
	radii.newsize(natom);

	HaAtom* aptr;
	int iat = 0;
//	double low_pack = pack_dens - sigma_pack_dens;
//	double up_pack = pack_dens + sigma_pack_dens;
	AtomIteratorMolSet aitr(pmset);
	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	{
		
		coord.SetVal_idx0( 0,iat, aptr-> GetX());
		coord.SetVal_idx0(1,iat, aptr-> GetY());
		coord.SetVal_idx0(2,iat, aptr-> GetZ());
		radii(iat+1)   = aptr->radius;
		coord1(1,iat+1) = aptr->GetX();
		coord1(2,iat+1) = aptr->GetY();
		coord1(3,iat+1) = aptr->GetZ();
		radii1(iat+1)   = aptr->radius;
		iat++;
//		PrintLog ("radius  %s %2.1f \n" , aptr->GetName(), radii(iat) );
	} 
clock_t tt1 = clock();
	int sw =0 ;
	HaSurface surf;
	surf.CalcMolSurfAlpha(sw, solv_rad, coord, radii);
	double summm= 0;
	int i;
	for(i=0; i < natom; i++)
	{
//    summm = summm + surf.surface_alpha[i] ;
		}  
	PrintLog("summm in PenaltyDensity %12.1f \n", summm);
clock_t tt2 = clock();
printf("%2.2f clock ticks between  tt2-tt1 CalcMolSurfAlpha \n", (double)(tt2-tt1)/CLOCKS_PER_SEC);

//	HaSurface surf1;
//	surf1.CalcMolSurf(surf_type, solv_rad, coord1, radii1);
//    svolume= surf.GetSurfVolume();
//	pack_dens_calculated = svolume/atomic_volume;
//	pack_dens_ene = SqrPotential(pack_dens_calculated, pack_dens, sigma_pack_dens, weight_pack_dens);
clock_t tt3 = clock();
printf("%2.2f clock ticks between  tt2-tt1 CalcMolSurf \n", (double)(tt3-tt2)/CLOCKS_PER_SEC);

	return pack_dens_ene;
}

double
HaEmpiricalMod::PenaltySymmetry()
{
	//		fstream component_file;// Remove
	//	component_file.open("energy_components.dat",ios::out | ios::app);// Remove
	//	char buf[256];// Remove
	
	if (com_flag)
	{
		CenterOfMass();
		FindAxes();
	}
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	Vec3D mol_center ;
	Vec3D* pmolcenter ;
	double numer_x = 0.0;
	double numer_y = 0.0;
	double numer_z = 0.0;
	double com_x = 0.0;
	double com_y = 0.0;
	double com_z = 0.0;
	double denumer = 0.0;
	double sym_ene =0.0 ;
	double cur_sym_ene =0.0 ;
	HaVec_int nmol_chain_arr;
	int ichain_old =1 ;
	int ichain =1 ;
	double delta_dist, dist1, dist2;
	double low_bound = 0.0 ;
	double up_bound = 0.0 + sigma_sym ; 
	int j = 0;
	nmol_chain_arr.newsize(nchain);
	
	// Calculate COM of molecule
	int i;
	for( i =0; i < nmol ; i++ )
	{
		Vec3D* ptr = &center_arr[i];
		numer_x += (ptr -> GetX());
		numer_y += (ptr -> GetY());
		numer_z += (ptr -> GetZ());
		
	}
	double denom = 1/nmol;
	com_x = (numer_x*denom) ;
	com_y = (numer_y*denom) ;
	com_z = (numer_z*denom) ;
	
	mol_center.SetX(com_x);
	mol_center.SetY(com_y);
	mol_center.SetZ(com_z);
	pmolcenter = &mol_center;
	
	int n =0;
	for( i =0; i < nmol ; i++ )
	{
	
		ichain = chain_arr[i];
		if (ichain == ichain_old)
		{
			n++;
		}
		else
		{
			nmol_chain_arr[j] = n; // define how many molecules in each chain
			ichain_old = ichain ;
			n = 1;
			j++ ;
		}
		if(i+1 == nmol)
		{
			nmol_chain_arr[j] = n;
		}
	}
	
	int n_count = nmol/nchain;
	for(j =0; j < n_count ; j++ )
	{
		for( int k =j; k < nmol ; k += n_count )
		{
			Vec3D* ptr1 = &center_arr[k];
			dist1 = Vec3D::CalcDistance(ptr1,pmolcenter,ANGSTROM_U);
			for( int m = k; m < nmol ; m += n_count )
			{
				if(k != m)
				{
					Vec3D* ptr2 = &center_arr[m];
					dist2 = Vec3D::CalcDistance(ptr2,pmolcenter,ANGSTROM_U);
					delta_dist = sqrt((dist1-dist2)*(dist1-dist2));
					cur_sym_ene += SqrPotential(delta_dist, low_bound, sigma_sym, weight_sym);
//			PrintLog ("k %d m %d \n", k, m);
//			PrintLog ("dist1 %5.2f, dist2 %3.2f  \n", dist1, dist2);
//			PrintLog ("C-H delta_dist %5.2f, cur_sym_ene comToCenter %3.2f  \n", delta_dist, cur_sym_ene);
				}
			}
		}
	}
	
//	double len_mol_center = sqrt(com_x*com_x + com_y*com_y + com_z*com_z);
	//	PrintLog ("Length of mol_center %2.3f\n", len_mol_center);

	Vec3D vec1;
	Vec3D vec2;
	Vec3D* pnt_vec1;
	Vec3D* pnt_vec2;
	HaVec_double side_arr, hel_central_arr;
	double temp1, temp2;
	double angle1, delta_angle, angle2;
	side_arr.newsize(nchain);
	hel_central_arr.newsize(nchain);

	for(j =0; j < n_count ; j++ )
	{
		i = 0;
		for( n =j; n < nmol ; n +=n_count )
		{
			//PrintLog ("ichain1 %i, chain2 %i  \n", ichain, n);
			vec1 = axis_arr[n];
			temp1 = vec1[0]*central_axis[0] + vec1[1]*central_axis[1] + vec1[2]*central_axis[2] ;
			angle1 = (acos(temp1))*RAD_TO_DEG; 
			
			//			Vec3D center_helix1 =  center_arr[n] ;
			
			//			cur_x = center_helix1[0] ;
			//			cur_y = center_helix1[1] ;
			//			cur_z = center_helix1[2] ;
			//			len1 = sqrt(cur_x*cur_x + cur_y*cur_y + cur_z*cur_z); 
			//			double ux = vec1[0]* len1;
			//			double uy = vec1[1]* len1;
			//			double uz = vec1[2]* len1;
//			vec1.SetX(vec1[0]);
//			vec1.SetY(vec1[1]);
//			vec1.SetZ(vec1[2]);
			pnt_vec1= &vec1;
			
			//			Vec3D center_helix2 ;		
			if (n + nmol/nchain < nmol )
			{
				vec2 = axis_arr[n+n_count];
				//				center_helix2 =  center_arr[n+nmol/nchain] ;
			}
			else
			{
				vec2 = axis_arr[j];
				//				center_helix2 =  center_arr[j] ;	
			}
			//			cur_x = center_helix2[0] ;
			//			cur_y = center_helix2[1] ;
			//			cur_z = center_helix2[2] ;
			//			len2 = sqrt(cur_x*cur_x + cur_y*cur_y + cur_z*cur_z);
			//			PrintLog ("LEN1 %2.1f LEN2 %2.1f Mol_N %2d Mol_N+ %2d\n", len1, len2, n, n + nmol/nchain);
			
			//			double vx = vec2[0]* len2;
			//			double vy = vec2[1]* len2;
			//			double vz = vec2[2]* len2;
//			vec2.SetX(vec2[0]);
//			vec2.SetY(vec2[1]);
//			vec2.SetZ(vec2[2]);
			pnt_vec2= &vec2;
			//PrintLog ("vec1 x %2.3f, y %2.3f, z %2.3f \n", vec1[0], vec1[1], vec1[2]);
			//PrintLog ("vec2 x %2.3f, y %2.3f, z %2.3f \n", vec2[0], vec2[1], vec2[2]);
//			double dist = Vec3D::CalcDistance(pnt_vec1,pnt_vec2,ANGSTROM_U);
			temp2 = vec2[0]*central_axis[0] + vec2[1]*central_axis[1] + vec2[2]*central_axis[2] ;

			angle2 = (acos(temp2))*RAD_TO_DEG; 
			delta_angle = sqrt( (angle2 - angle1)*(angle2 - angle1));
			hel_central_arr[i] = delta_angle;
//		PrintLog ("delta_angle %2.1f\n", delta_angle);
			temp1 = (vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]) ;
			angle1 = (acos(temp1))*RAD_TO_DEG; 

			side_arr[i] = angle1;
			i++;
		}
		//	sprintf(buf,"Molecule %d", j+1); // Remove
		//	component_file << buf<< endl ; // Remove
//		PrintLog ("SYMMETRY before E= %2.3f  \n", cur_sym_ene);
		for(n =0; n < nchain -1 ; n++ )
		{
//	PrintLog ("side_arr[n] %2.1f\n", side_arr[n]);
//		PrintLog ("hel_central_arr[n] %2.1f\n", hel_central_arr[n]);

			delta_angle = sqrt( (side_arr[n]-side_arr[n+1])*(side_arr[n]-side_arr[n+1]) );
			if (pack_dens_calculated <= pack_dens+sigma_pack_dens)
			{
				cur_sym_ene += SqrPotential(delta_angle, low_bound, sigma_sym, weight_sym);
//PrintLog ("H-H delta_angle %5.2f, cur_sym_ene comToCenter %3.2f  \n", delta_angle, cur_sym_ene);
				cur_sym_ene += SqrPotential(hel_central_arr[n], low_bound, sigma_sym, weight_sym);
//			PrintLog ("n %d \n", n);
//			PrintLog ("C-H hel_central_arr[n] %5.2f, cur_sym_ene comToCenter %3.2f  \n", hel_central_arr[n], cur_sym_ene);

			}
			else
			{
				cur_sym_ene += SqrPotential(delta_angle, low_bound, sigma_sym, weight_sym);
				cur_sym_ene += SqrPotential(hel_central_arr[n], low_bound, sigma_sym, weight_sym);
			}
			//			PrintLog ("SYMMETRY1 Angle axis E= %2.3f  \n", cur_sym_ene);

//		PrintLog ("SYMMETRY2  Angle axis to center E= %2.3f  \n", cur_sym_ene);
//			PrintLog ("delta_angle1 = %2.3f  \n", delta_angle);
//			PrintLog ("delta_angle2 = %2.3f  \n", hel_central_arr[n]);

		//	sprintf(buf,"Side lenght %4.3f ", side_arr[n]); // Remove
			//	component_file << buf ; // Remove
		}
//			PrintLog ("side_arr[nchain-1] %2.1f\n", side_arr[nchain-1]);
		delta_angle = sqrt( (side_arr[nchain -1]-side_arr[0])*(side_arr[nchain -1]-side_arr[0]) );
		if (pack_dens_calculated <= pack_dens+sigma_pack_dens)
		{
			cur_sym_ene += SqrPotential(delta_angle, low_bound, sigma_sym, weight_sym);
			cur_sym_ene += SqrPotential(hel_central_arr[nchain-1], low_bound, sigma_sym, weight_sym);
		}
		else
		{
			cur_sym_ene += SqrPotential(delta_angle, low_bound, sigma_sym, weight_sym);
			cur_sym_ene += SqrPotential(hel_central_arr[nchain-1], low_bound, sigma_sym, weight_sym);
		}

		//	sprintf(buf,"Side lenght %4.3f ", side_arr[nchain -1]); // Remove
		//	component_file << buf << endl; // Remove
	}
	for(j =0; j < n_count ; j++ )
	{
		i = 0;
		for( n =j; n < nmol ; n +=n_count )
		{
			vec1 = center_arr[n];
			pnt_vec1 = &vec1;
			if (n + n_count < nmol )
			{
				vec2 = center_arr[n+n_count];
			}
			else
			{
				vec2 = center_arr[j];
			}
			double dist = Vec3D::CalcDistance(pnt_vec1,pnt_vec2,ANGSTROM_U);
			if (i%2 == 0)
			{
				side_arr[i] = dist;
			//	PrintLog ("side1 %2.3f , dist %2.3f\n", side_arr[i], dist);
			}
			else
			{
				side_arr[i] = dist* lenght_factor;
			//	PrintLog ("side2 %2.3f , dist %2.3f \n", side_arr[i], dist);
			}
			i++;
		}
		for(n =0; n < nchain -1 ; n++ )
		{
			delta_dist = sqrt( (side_arr[n]-side_arr[n+1])*(side_arr[n]-side_arr[n+1]) );
			cur_sym_ene += SqrPotential(delta_dist, low_bound, sigma_sym, weight_sym);
//			PrintLog ("SYMMETRY3  com to com n %d delta_dist %2.3f E= %2.3f  \n", n,delta_dist, cur_sym_ene);
		}
		delta_dist = sqrt( (side_arr[nchain -1]-side_arr[0]) *(side_arr[nchain -1]-side_arr[0]) );
		cur_sym_ene += SqrPotential(delta_dist, low_bound, sigma_sym, weight_sym);
//		PrintLog ("SYMMETRY3  com to com n %d delta_dist %2.3f E= %2.3f  \n", n,delta_dist, cur_sym_ene);

	} 
	return sym_ene = cur_sym_ene ;
	//component_file.close(); // Remove
}


int
HaEmpiricalMod::LoadEmpConstrains()
{
	char buf[256];
	int ires = 0;
	npt_dist = 0;
	HaAtom* aptr1;
	HaAtom* aptr2;
	
	MolSet* pmset = GetMolSet();	
	int nmol = pmset -> GetNMol();
	std::string fname = "restrains.dat";
	FILE* finfo = fopen(fname.c_str(),"r");
	topology_arr.newsize(nmol);
	int imol;
	for( imol =0; imol < nmol ; imol++)
	{
		topology_arr[imol]= 0;
	}
	
	if(finfo == NULL)
	{
		PrintLog("No distance constrains loaded.\n");
		
		finfo = fopen(fname.c_str(),"w");
		
		fprintf(finfo,"%s\n", "DISTANCE_CONSTRAINS");
		fprintf(finfo,"%s\n", "# Format: \n# $MolID$ResIDResNumber:ChainID.AtomID $MolID$ResIDResNumber:ChainID.AtomID Distance");
		fprintf(finfo,"%s\n", "DONE");
		
		fprintf(finfo,"%s\n", "TOPOLOGY");
		fprintf(finfo,"%s\n", "# Format: \n# $MolID$ Inside");
		fprintf(finfo,"%s\n", "DONE");
		fclose(finfo);
		
		return FALSE;
	}
	
	char* cres = NULL;
	
	
	for(;;)
	{
		cres = fgets(buf,255,finfo);
		if( strncmp(buf,"DISTANCE_CONSTRAINS",19) == 0)
		{
			for(;;)
			{
				cres = fgets(buf,255,finfo); 
				npt_dist++;
				//				PrintLog("%s\n", cres);
				
				if(buf[0] == '#') continue;
				if( strncmp(buf,"DONE",4) == 0) break;
			}
		}
		if(buf[0] == '#') continue;
		if( strncmp(buf,"DONE",4) == 0) break;
	}
	
	emp_dist.newsize(npt_dist);
	atm_dist.resize(npt_dist);
	atm2_dist.resize(npt_dist);
	wxString str;
	int i;
	//	PrintLog("npt_dist %d\n", npt_dist);


	fseek(finfo,0,SEEK_SET);
	
	for(;;)
	{
		cres = fgets(buf,255,finfo); 
		if( strncmp(buf,"DISTANCE_CONSTRAINS",19) == 0)
		{
			PrintLog ("Read Distance Constraints \n");
			for(;;)
			{
				cres = fgets(buf,255,finfo); 
				if(buf[0] == '#') continue;
				if( strncmp(buf,"DONE",4) == 0) break;
				
				str = buf;  
				wxArrayString sub_str;
				wxStringTokenizer tkz(str," ");
				while ( tkz.HasMoreTokens())
				{
					wxString token = tkz.GetNextToken();
					sub_str.Add(token);
				}
				int nsub_str = sub_str.Count();
				for(i =0; i < nsub_str; i++)
				{
					aptr1 = pmset->GetAtomByRef( sub_str[0].ToStdString() );
					aptr2 = pmset->GetAtomByRef( sub_str[1].ToStdString() );
					if (aptr1 != NULL && aptr2 != NULL )
					{
						atm_dist[ires]=  aptr1 ;
						atm2_dist[ires]=  aptr2 ;
						bool res = sub_str[2].ToDouble(&emp_dist[ires]) ;
					}
					else
					{
						atm_dist[ires] = NULL;
						atm2_dist[ires] = NULL;
						emp_dist[ires] = 0.0 ;
					}
				}
				HaResidue* pres1 = aptr1->GetHostRes(); 
				HaResidue* pres2 = aptr2->GetHostRes(); 
				
				PrintLog ("Constraint %i: Atom1 %s.%s , Atom2 %s.%s, Distance %4.2f Angstrom\n", ires+1, 
					pres1->GetName(), aptr1 ->GetName(), pres1->GetName(), 
					aptr2 ->GetName(), emp_dist[ires]);  // Printout of Constraints
				ires++;
			}
		if (ires ==0) PrintLog("No distance constrains loaded.\n\n");
		}
		
		wxString  mol_name;
		wxString  mol_name1;
		HaMolecule* pmol;
		
		if( strncmp(buf,"TOPOLOGY",8) == 0)
		{
			PrintLog ("Read Topology \n");
			for(;;)
			{
				cres = fgets(buf,255,finfo); 
				if(buf[0] == '#') continue;
				if( strncmp(buf,"DONE",4) == 0) break;
				
				str = buf;  
				wxArrayString sub_str;
				wxStringTokenizer tkz(str," ");
				while ( tkz.HasMoreTokens())
				{
					wxString token = tkz.GetNextToken();
					sub_str.Add(token);
				}
				mol_name = sub_str[0];
				
				for( imol =0; imol < nmol ; imol++)
				{
					pmol = pmset -> GetMolByIdx(imol);
					mol_name1 = (pmol ->GetRef()).c_str();
					mol_name1.Append("$") ;
					mol_name1.Prepend("$") ;
					if (mol_name == mol_name1)
					{
						PrintLog ("Molecule# %d %s   %s \n", imol, mol_name.ToStdString().c_str(), (sub_str[1]).ToStdString().c_str());
						if ( strncmp(sub_str[1].mb_str(),"Inside",6) == 0 )	topology_arr[imol] = 1;
						if ( strncmp(sub_str[1].mb_str(),"inside",6) == 0 )	topology_arr[imol] = 1;
						if ( strncmp(sub_str[1].mb_str(),"INSIDE",6) == 0 )	topology_arr[imol] = 1;
						
					}
				}
			}
			if (mol_name.Len() == 0) PrintLog("No topology loaded.\n\n"); 
			
		}
		if(cres == NULL) break ;
	}
	fclose(finfo);
	return TRUE;
}

int 
HaEmpiricalMod::LoadSolventAccessibleAtoms()
{
	char buf[256];
	int ires = 0;
	n_sa_atoms = 0;
	HaAtom* aptr;
	HaAtom* aptr1;
	int imol;
	int i;
	int j;
	HaMolecule* pMol; 
	HaResidue* res_ptr;
	MolSet* pmset = GetMolSet();	
	int nmol = pmset -> GetNMol();
	int nres = pmset -> GetNRes();
	int natom = nres*nmol;
	marker_res_sa.newsize(natom);


	if (atm_sc_array[0] == NULL) CbettaSetUp();

	std::string fname = "solventaccesible.dat";
	FILE* finfo = fopen(fname.c_str(),"r");
	i = 0;
	if(finfo == NULL)
	{
		PrintLog("No solvent accesible residues loaded\n");
		finfo = fopen(fname.c_str(),"w");
		fprintf(finfo,"%s\n", "#Format: \n# $MOL_ID$Res_IDRes_Number:Chain.X Weight");
		fprintf(finfo,"%s\n", "#Must be X type atoms!\n");
		fclose(finfo);

		for( imol =0; imol < nmol ; imol++)
		{
			pMol = pmset->GetMolByIdx(imol);
			HaMolecule::ResidueIterator ritr(pMol);
			for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
			{
				marker_res_sa[i] = 0;
				i++;	
			}

		}
		return FALSE;
	}
	
	char* cres = NULL;


	for(;;)
	{
		cres = fgets(buf,255,finfo);
		if(buf[0] == '#') continue;
		if(cres == NULL) break;
		if(buf[0] == '\n') break;
		n_sa_atoms ++;
	}

	atm_solacces.resize(n_sa_atoms );

	wxString str;
	i = 0;
	fseek(finfo,0,SEEK_SET);
	PrintLog ("Load Solvent Accessible Residues \n");
	for(;;)
	{
		cres = fgets(buf,255,finfo); 
		
		if(cres == NULL) break;
		if(buf[0] == '#') continue;
		str = buf;  
		wxArrayString sub_str;
		wxStringTokenizer tkz(str," ");
		while ( tkz.HasMoreTokens())
		{
			wxString token = tkz.GetNextToken();
			sub_str.Add(token);
		}
		int nsub_str = sub_str.Count();
		if (sub_str[0] == "\n") break; 
		for(int i =0; i < nsub_str; i++)
		{
			aptr1 = pmset->GetAtomByRef(sub_str[0].ToStdString() );
			if (aptr1 != NULL )
			{
				atm_solacces[ires]=  aptr1 ;
			}
			else
			{
				atm_solacces[ires] = NULL;
			}
		}
		HaResidue* pres1 = aptr1->GetHostRes(); 
		
		PrintLog ("Solvent Accessible Residue %i: %s%d \n", ires+1, 
			pres1->GetName(), pres1-> GetSerNo() );  // Printout of solvent accessible residues
		ires++;
	}
	if (n_sa_atoms  == 0) PrintLog("No solvent accesible residues loaded\n\n");
	i = 0;	
	for( i =0; i < nres ; i++)
	{
		marker_res_sa[i] = 0;
		aptr = (HaAtom*) atm_sc_array[i];
		for(j=0 ; j < n_sa_atoms ; j++)
		{
			aptr1 = (HaAtom*) atm_solacces[j];
			if (aptr1 != NULL)
			{
				if( aptr1 == aptr) marker_res_sa[i] = 1;
			}
		}
	}
	fclose(finfo);
	return TRUE;
}



int 
HaEmpiricalMod::EstablishChains()
{
	// This function assigns a number for each different chain
	int ichain = 0;
	MolSet* pmset = GetMolSet();
	HaChain*    phost_chain = NULL ;
	HaMolecule* pMol;
	HaAtom* aptr;
	char chain_letter = 0;
	char chain_letter_old = 0;
	int nmol = pmset->GetNMol();

	chain_arr.newsize(nmol);
	for( int imol =0; imol < nmol ; imol++)
	{
		pMol = pmset->GetMolByIdx(imol);
		AtomIteratorMolecule aitr(pMol);
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			phost_chain= aptr->GetHostChain();
			chain_letter = phost_chain -> ident;
			if( chain_letter != chain_letter_old )
			{
				chain_letter_old = chain_letter;
				ichain ++ ;
				//PrintLog (" letter %c \n", chain_letter ); 
			}	
		}
		chain_arr[imol] = ichain;
		nchain = ichain;
		PrintLog("chain %2d\n",ichain); //added by Jose
	}
	return TRUE;
}

double
HaEmpiricalMod::SoftSqrWellPotential(double& current_value, double& average, double& stdev, double& weight)
{
	double potential =0.0 ;
	double denom = stdev*stdev;
	double up_bound = average + stdev ;
	double low_bound = average - stdev ;
	potential = pow((current_value - average),2)/(denom*denom);

	if (current_value < low_bound)
	{	
		potential += pow((current_value - stdev),2)/denom;
	}
	if (current_value > up_bound )
	{
		potential +=  pow((up_bound - current_value),2)/denom;
	}
	return potential*weight;
}

double
HaEmpiricalMod::SqrPotential(double& current_value, double& average, double& stdev, double& weight)
{
	double potential =0.0 ;
	double denom = stdev*stdev;
	potential = (current_value - average)*(current_value - average)/denom;
	return potential*weight;
}

/*
Vec3DValArray
HaEmpiricalMod::FindAxes()
{ 
	if (com_flag)
	{
		CenterOfMass();
	}

	HaAtom* aptr;
	wxString name ;
	Vec3D vec_com1;
	Vec3D vec_com2;
	HaMolecule* pMol;
	Vec3D center_helix;
	HaMat_double fmat;
	HaMat_double cc(3,3);
	HaVec_double eig_val;
	//fmat.newsize(nr,nc);
	fmat.newsize(3,3);
	eig_val.newsize(3);
	axis_vec.newsize(3);

	
	Vec3D point;
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	
	axis_arr.resize(nmol);
	int imol;
	int ind1 = 0;
	int ind2 = 0;
	int j;

	for( imol =0; imol < nmol ; imol++)
	{
		pMol = pmset->GetMolByIdx(imol);
		double c11 =0.0;
		double c12 =0.0;
		double c13 =0.0;
		double c21 =0.0;
		double c22 =0.0;
		double c23 =0.0;
		double c31 =0.0;
		double c32 =0.0;
		double c33 =0.0;
		double cur_x = 0.0;
		double cur_y = 0.0;
		double cur_z = 0.0;
		double  eig_val_minim = 9999999999999.0 ;
		int indx = 0;
		double delta_x = 0.0;
		double delta_y = 0.0;
		double delta_z = 0.0;
		
		//		buf1 = pMol-> GetRef();
		//		sprintf(buf,"Molname %s, imol %i \n", buf1.c_str(), imol); // Remove
		//		component_file << buf << endl; // Remove
		//PrintLog(" Molname  %s, imol %i \n", buf1.c_str(), imol); // Remove
		
		//		FindAxis(pMol);
		//		point.SetX(axis_vec[0]);
		//		point.SetY(axis_vec[1]);
		//		point.SetZ(axis_vec[2]);
		
		
		
		center_helix =  center_arr[imol] ;
		AtomIteratorMolecule aitr(pMol);
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			name = aptr->GetName();
			if(name == "C" || name == "CA" || name == "N")
			{
				cur_x = aptr->GetX() - center_helix[0];
				cur_y=  aptr->GetY() - center_helix[1];
				cur_z = aptr->GetZ() - center_helix[2];
				
				c11+= (cur_y*cur_y+ cur_z*cur_z)*(aptr->GetMass());
				c22+= (cur_x*cur_x+ cur_z*cur_z)*(aptr->GetMass());
				c33+= (cur_x*cur_x+ cur_y*cur_y)*(aptr->GetMass());
				c12 = c12 - cur_x*cur_y*(aptr->GetMass());
				c13 = c13 - cur_x*cur_z*(aptr->GetMass());
				c23 = c23 - cur_y*cur_z*(aptr->GetMass());
			}
		}
		
		fmat.SetVal_idx0(0,0,c11);
		fmat.SetVal_idx0(0,1,c12);
		fmat.SetVal_idx0(0,2,c13);
		fmat.SetVal_idx0(1,0,c12);
		fmat.SetVal_idx0(1,1,c22);
		fmat.SetVal_idx0(1,2,c23);
		fmat.SetVal_idx0(2,0,c13);
		fmat.SetVal_idx0(2,1,c23);
		fmat.SetVal_idx0(2,2,c33);
		HaMat_double::mat_sdiag(fmat, cc, eig_val);
		
		for(j = 1; j <= 3; j++)
		{
			if (eig_val(j) < eig_val_minim)
				//if (eig_val(j)	> eig_val_max)
			{
				indx = j;
				eig_val_minim = eig_val(j);
				//eig_val_max = eig_val(j);
			}
			//PrintLog(" E_minim(indx=%3d)= %12.6f , eig_val(j) %12.6f\n",indx,eig_val_minim,eig_val(j) );
		}
		
//		for( i = 1; i <= 3; i++)
//		{
			
			for( j = 1; j <= 3; j++)
			{
				axis_vec[j-1] = cc(j, indx); 
				//PrintLog(" cc( %3d %3d)= %4.6f \n",i, j, cc(i,j) );
				// PrintLog(" axis_vec( %3d)= %4.6f \n",j,axis_vec[j-1] );
			}
//		}
		point.SetX(axis_vec[0]);
		point.SetY(axis_vec[1]);
		point.SetZ(axis_vec[2]);
		axis_arr[imol].SetCoordFrom(point);
		
	}
	for( imol =0; imol < nmol-1 ; imol++)
	{
		
		Vec3D vec1 = axis_arr[imol];
		Vec3D vec2 = axis_arr[imol+1];
		double ux =  vec1[0];
		double uy = vec1[1];
		double uz = vec1[2];

		double vx =  vec2[0];
		double vy = vec2[1];
		double vz = vec2[2];
				
		double temp = ux*vx + uy*vy + uz*vz;
		double angle = (acos(temp))*RAD_TO_DEG; 
		if (angle > 90.0)
		{
			point.SetX(-vec2[0]);
			point.SetY(-vec2[1]);
			point.SetZ(-vec2[2]);
			axis_arr[imol+1].SetCoordFrom(point);
		}
	}
	return axis_arr ;
}
*/

Vec3DValArray
HaEmpiricalMod::FindAxes()
{ 
	HaAtom* aptr;
	Vec3D cur_com_atom1;
	Vec3D cur_com_atom2;

	double numer_x,numer_y, numer_z;
	double com_x, com_y, com_z;
	Vec3D point;
	int imol, i;
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	axis_arr.resize(nmol);

	for( imol =0; imol < nmol ; imol++)
	{
		int istart = imol * 8;
		numer_x = 0.0;
		numer_y = 0.0;
		numer_z = 0.0;
		com_x = 0.0;
		com_y = 0.0;
		com_z = 0.0;
		for(i = istart; i < istart+4 ; i++)
		{
			aptr = (HaAtom*) segment_vec[i];
			numer_x += aptr -> GetX();
			numer_y += aptr -> GetY();
			numer_z += aptr -> GetZ();
		}
		com_x = (numer_x*0.25) ;
		com_y = (numer_y*0.25) ;
		com_z = (numer_z*0.25) ;
		
		cur_com_atom1.SetX(com_x);
		cur_com_atom1.SetY(com_y);
		cur_com_atom1.SetZ(com_z);
		
		numer_x = 0.0;
		numer_y = 0.0;
		numer_z = 0.0;
		com_x = 0.0;
		com_y = 0.0;
		com_z = 0.0;
		for(i = istart+4; i <istart+8 ; i++)
		{
			aptr = (HaAtom*) segment_vec[i];
			numer_x += aptr -> GetX();
			numer_y += aptr -> GetY();
			numer_z += aptr -> GetZ();
		}
		com_x = (numer_x*0.25) ;
		com_y = (numer_y*0.25) ;
		com_z = (numer_z*0.25) ;
		
		cur_com_atom2.SetX(com_x);
		cur_com_atom2.SetY(com_y);
		cur_com_atom2.SetZ(com_z);
	

		double length ;
		com_x = cur_com_atom1.GetX() - cur_com_atom2.GetX();
		com_y = cur_com_atom1.GetY() - cur_com_atom2.GetY();
		com_z = cur_com_atom1.GetZ() - cur_com_atom2.GetZ();
		length = 1/(sqrt(com_x*com_x + com_y*com_y +com_z*com_z) );

		point.SetX(com_x*length);
		point.SetY(com_y*length);
		point.SetZ(com_z*length);
		axis_arr[imol].SetCoordFrom(point);
	}

	return axis_arr ;
}





Vec3DValArray
HaEmpiricalMod::CenterOfMass()
{ 
	MolSet* pmset = GetMolSet();
	HaMolecule* pMol;
	HaAtom* aptr;
	int nmol = pmset->GetNMol();
	center_arr.resize(nmol);
	Vec3D cur_com_atom;
	double numer_x;
	double numer_y;
	double numer_z;
	double com_x;
	double com_y;
	double com_z;
	double com_x1;
	double com_y1;
	double com_z1;

	double denumer = 0.0;
	int i, imol;

	for( imol =0; imol < nmol ; imol++)
	{
		int istart = imol * 8;
		numer_x = 0.0;
		numer_y = 0.0;
		numer_z = 0.0;
		com_x = 0.0;
		com_y = 0.0;
		com_z = 0.0;
		for(i = istart; i < istart+4 ; i++)
		{
			aptr = (HaAtom*) segment_vec[i];
			numer_x += aptr -> GetX();
			numer_y += aptr -> GetY();
			numer_z += aptr -> GetZ();
		}
		com_x = (numer_x*0.25) ;
		com_y = (numer_y*0.25) ;
		com_z = (numer_z*0.25) ;

	
		numer_x = 0.0;
		numer_y = 0.0;
		numer_z = 0.0;
		com_x1 = 0.0;
		com_y1 = 0.0;
		com_z1 = 0.0;
		for(i = istart+4; i <istart+8 ; i++)
		{
			aptr = (HaAtom*) segment_vec[i];
			numer_x += aptr -> GetX();
			numer_y += aptr -> GetY();
			numer_z += aptr -> GetZ();
		}
		com_x1 = (numer_x*0.25) ;
		com_y1 = (numer_y*0.25) ;
		com_z1 = (numer_z*0.25) ;
		
		com_x = (com_x + com_x1) * 0.5;
		com_y = (com_y + com_y1) * 0.5;
		com_z = (com_z + com_z1) * 0.5;
	
		cur_com_atom.SetX(com_x);
		cur_com_atom.SetY(com_y);
		cur_com_atom.SetZ(com_z);
   	   	center_arr[imol].SetCoordFrom(cur_com_atom);
	}
	return center_arr;
}



int HaEmpiricalMod::LoadEmpParam()
{
	MolSet* pmset = GetMolSet();	
	std::string fname = "empir_param.dat";
	char buf[256];

	char* cres = NULL;
	int i = 0;
	int j = 0;
	int k = 0;
	long col; // column index
	long raw; // raw index
	long bin; // raw index
	double ene; // energy

	la_value.newsize(32);
	la_weight_value.newsize(32);
	sc_vdwradius.newsize(32);
	pairwise_energy_arr.SetDimensions(32,32,32);
	pairwise_energy_arr_sa.SetDimensions(32,32,32);
	wxString str;

	FILE* finfo = fopen(fname.c_str(),"r");
	if(finfo == NULL)
	{
		ErrorInMod("HaEmpiricalMod::LoadEmpParam",
			"Can not open EMPIR_PARAM.DAT file");

		for (i=0; i<32; i++)
		{
			la_value[i] = 0.0;
			la_weight_value[i] =  0.0;
			sc_vdwradius[i] =  0.0;
		}
		pairwise_energy_arr.FillZeros();
		pairwise_energy_arr_sa.FillZeros();

		return FALSE;
	}

	int ncount_pair ;

	for(;;)
	{
		cres = fgets(buf,255,finfo); 
		if( strncmp(buf,"LIPID_ACCESS_SCALE",18) == 0)
		{
			PrintLog ("Read Lipid Propensity Scale \n");
			i = 0 ;
			for(;;)
			{

				cres = fgets(buf,255,finfo); 
				if(buf[0] == '#') continue;
				if( strncmp(buf,"DONE",4) == 0) break;

				str = buf;  
				wxArrayString sub_str;
				wxStringTokenizer tkz(str," ");
				while ( tkz.HasMoreTokens())
				{
					wxString token = tkz.GetNextToken();
					sub_str.Add(token);
				}

				pairwise_name_arr.Add(sub_str[0]);
				residue_arr.Add(sub_str[0]);
				bool res = sub_str[1].ToDouble(&la_value[i]) ;
				bool res1 = sub_str[2].ToDouble(&la_weight_value[i]) ;
				PrintLog ("%s  Value %5.2f Weigth %5.2f\n", residue_arr[i].ToStdString().c_str(), la_value[i], la_weight_value[i]);  // Printout of Read Values
				i++;
			}
		}
		if( strncmp(buf,"UNRES_PARAMS",12) == 0)
		{
			PrintLog ("Read vdW radia \n");
			i = 0 ;
			for(;;)
			{
				cres = fgets(buf,255,finfo); 
				if(buf[0] == '#') continue;
				if( strncmp(buf,"DONE",4) == 0) break;

				str = buf;  
				wxArrayString sub_str;
				wxStringTokenizer tkz(str," ");
				while ( tkz.HasMoreTokens())
				{
					wxString token = tkz.GetNextToken();
					sub_str.Add(token);
				}
	
				residue_unres_arr.Add(sub_str[0]);
				bool res3 = sub_str[1].ToDouble(&sc_vdwradius[i]) ;
				
				PrintLog ("%s  Side Chain vdW radius %5.2f\n", residue_unres_arr[i].ToStdString().c_str(), sc_vdwradius[i]);  // Printout of Values
				i++;
			}

			if (residue_unres_arr.GetCount() != residue_arr.GetCount()) 
			{
				PrintLog("Only %d radii is recorded! Add radius to UNRES_PARAMS.\n", residue_unres_arr.GetCount() );
			}
		}
		if( strncmp(buf,"PAIRWISE_PARAMS",15) == 0)
		{
			PrintLog ("Read pairwise energies \n");
			int j;
			double tmp;
			str = buf;  
			wxArrayString sub_str;
			wxStringTokenizer tkz(str," ");
			pairwise_energy_arr.FillZeros();
			i=0;
			ncount_pair = pairwise_name_arr.GetCount();

			if (ncount_pair ==0)
			{
				PrintLog ("# LIPID ACCESSIBLE SCALE MUST BE THE FIRST \n# PART IN EMPIR_PARAM.DAT FILE\n");
			}

			for(;;)
			{
				cres = fgets(buf,255,finfo); 
				col = 100;
				raw = 100;
				if(buf[0] == '#') continue;
				if( strncmp(buf,"DONE",4) == 0) break;
				str = buf;  
				wxArrayString sub_str;
				wxStringTokenizer tkz(str," ");
				while ( tkz.HasMoreTokens())
				{
					wxString token = tkz.GetNextToken();
					sub_str.Add(token);
				}
				for(j=0; j< ncount_pair; j++)
				{
					if (pairwise_name_arr[j] == sub_str[0]) col =j ;
				}
				for(j=0; j< ncount_pair; j++)
				{
					if (pairwise_name_arr[j] == sub_str[1]) raw =j ;
				}

				if (col == 100) 
				{
					PrintLog("Found unknown residue %s. Add residue to lipid accessible scale.\n", sub_str[0].ToStdString().c_str());
					break;
				}
				if (raw == 100) 
				{
					PrintLog("Found unknown residue %s. Add residue to lipid accessible scale.\n", sub_str[1].ToStdString().c_str());
					break;
				}

				sub_str[3].ToLong(&bin);
				bin -= 1; 
				sub_str[2].ToDouble(&ene);

				pairwise_energy_arr.SetValue( col, raw, bin, ene);
//				PrintLog("pairwise_energy_arr( %d, %d,%d)= %2.3f \n", col, raw, bin, pairwise_energy_arr.GetValue( col, raw, bin) );
				i++;
			}
						if (i <4199) PrintLog("ERROR. Not complete PAIRWISE potential. Only %d parameters\n", i);
						else PrintLog("Successfully loaded core pairwise potential\n"); 
		}

		if( strncmp(buf,"SA_PAIRWISE_PARAMS",18) == 0)
		{
			PrintLog ("Read pairwise energies for solvent-accessible residues \n");
			double tmp;
			str = buf;  
			wxArrayString sub_str;
			wxStringTokenizer tkz(str," ");
			pairwise_energy_arr_sa.FillZeros();
			i=0;
			j=0;
			ncount_pair = pairwise_name_arr.GetCount();

			if (ncount_pair ==0)
			{
				PrintLog ("# LIPID ACCESSIBLE SCALE MUST BE THE FIRST \n # PART IN EMPIR_PARAM.DAT FILE\n");
			}


			for(;;)
			{
				cres = fgets(buf,255,finfo); 
				if(buf[0] == '#') continue;
				if( strncmp(buf,"DONE",4) == 0) break;
				str = buf; 
				col = 100;
				raw =100;
				wxArrayString sub_str;
				wxStringTokenizer tkz(str," ");
				while ( tkz.HasMoreTokens())
				{
					wxString token = tkz.GetNextToken();
					sub_str.Add(token);
				}
				for(j=0; j< ncount_pair; j++)
				{
					if (pairwise_name_arr[j] == sub_str[0]) col =j ;
				}
				for(j=0; j< ncount_pair; j++)
				{
					if (pairwise_name_arr[j] == sub_str[1]) raw =j ;
				}
			
				
				sub_str[3].ToLong(&bin);
				bin -= 1; 
				sub_str[2].ToDouble(&ene);

				if (col == 100) 
				{
					PrintLog("Found unknown residue %s. Add residue to lipid accessible scale.\n", sub_str[0].ToStdString().c_str());
					break;
				}
				if (raw == 100) 
				{
					PrintLog("Found unknown residue %s. Add residue to lipid accessible scale.\n", sub_str[1].ToStdString().c_str());
					break;
				}

				pairwise_energy_arr_sa.SetValue( col, raw, bin, ene);
//				PrintLog("pairwise_energy_arr_sa( %d, %d,%d)= %2.3f \n", col, raw, bin, pairwise_energy_arr_sa.GetValue( col, raw, bin) );
				i++;
			}
						if (i <4199) PrintLog("ERROR. Not complete PAIRWISE_SA potential. Only %d parameters\n", i);
						else PrintLog("Successfully loaded solvent-accessible pairwise potential\n\n"); 
		}



		if(cres == NULL) break ;
	}
	
	ncount_pair = pairwise_name_arr.GetCount();
	for(i=0; i< ncount_pair; i++)
	{
		for(j=0; j< ncount_pair; j++)
		{
			for(k=0; k< ncount_pair; k++)
			{
				ene = pairwise_energy_arr.GetValue( i, j, k);
				pairwise_energy_arr.SetValue( j, i, k, ene);
			}
		}
	}

	for(i=0; i< ncount_pair; i++)
	{
		for(j=0; j< ncount_pair; j++)
		{
			for(k=0; k< ncount_pair; k++)
			{
				ene = pairwise_energy_arr_sa.GetValue( i, j, k);
				pairwise_energy_arr_sa.SetValue( j, i, k, ene);
			}
		}
	}
				


	fclose(finfo);
	return TRUE;
}

int HaEmpiricalMod::CbettaSetUp()
{
	MolSet* pmset = GetMolSet();
	char buf[256];
	std::string mol_name;
	std::string res_name;
	std::string atm_name;
	int res_num; 
	wxString atm_ref; 
	HaAtom* aptr ;
	HaAtom* aptr1 ;
	HaResidue* res_ptr ;
	HaMolecule* pMol ;
	HaChain* phost_chain ;
	double dist_unres ;
	char chain_letter ; 
	double vdwradius_unres ;
	char atn[5];
	double atm_mass;
	double cogx, cogy, cogz;
	int natoms;
	int i;
	int sec_count =0;
	int sec_count1 =0;
	HaAtom* aptr_new;
//	fstream out_file;
//	out_file.open("atom_radia.txt",ios::out );
	int nres = pmset -> GetNRes();
	atm_sc_array.resize(nres);
	atm_ca_array.resize(nres);
	int ncont_res = residue_unres_arr.GetCount();

	
	MolSet::ResidueIterator ritr(pmset);
	for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
	{
		
		cogx =0;
		cogy =0;
		cogz =0;
		natoms = 0;
		res_name =res_ptr ->GetName();
		bool new_atm_flag = true;
		if (res_name != "GLY")
		{
			AtomIteratorResidue aitr(res_ptr);
			for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
			{
				atm_name = aptr-> GetName();
				if (!(aptr->IsHydrogen()) )
				{
					if (atm_name != "CA" && atm_name != "C" && atm_name != "N" && atm_name != "O" && atm_name != "OXT") //no peptide atoms
					{
						cogx += aptr->GetX();
						cogy += aptr->GetY();
						cogz += aptr->GetZ();
						natoms++;
					}
				}
				if (atm_name == "X")
				{
					new_atm_flag = false;
					HaColor green_c(2, 200, 200);
					aptr -> col= green_c.cidx;
					atm_sc_array[sec_count]=  aptr;
					for ( i=0; i< ncont_res; i++)
					{
						if (residue_unres_arr[i] == res_ptr ->GetName())
						{
							vdwradius_unres = sc_vdwradius[i] ;
						}
					}
					aptr -> radius= vdwradius_unres;
					
					sec_count++;
				}
				if (atm_name == "CA")
				{
					atm_ca_array[sec_count1]=  aptr;
					sec_count1++;
				}
			}
		}
		else
		{
			AtomIteratorResidue aitr(res_ptr);
			for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
			{
				atm_name = aptr-> GetName();
				if (atm_name == "CA")
				{
					cogx = aptr->GetX();
					cogy = aptr->GetY();
					cogz = aptr->GetZ();
					natoms =1;
					atm_ca_array[sec_count1]=  aptr;
					sec_count1++;
				}			
				if (atm_name == "X")
				{
					new_atm_flag = false;
					HaColor green_c(2, 200, 200);
					aptr -> col= green_c.cidx;
					atm_sc_array[sec_count]=  aptr;
					
					for ( i=0; i< ncont_res; i++)
					{
						if (residue_unres_arr[i] == res_ptr ->GetName())
						{
							vdwradius_unres = sc_vdwradius[i] ;
						}
					}
					aptr -> radius= vdwradius_unres;
					sec_count++;
				}
			}
			
		}	
		if(new_atm_flag)
		{				
			double denom = (double) 1/natoms;
			cogx *= denom;
			cogy *= denom;
			cogz *= denom;
			
			aptr_new = res_ptr ->AddNewAtom();
			aptr_new -> SetName("X");
			aptr_new->SetDummy();
			aptr_new -> SetX(cogx);
			aptr_new -> SetY(cogy);
			aptr_new -> SetZ(cogz);

			for ( i=0; i< ncont_res; i++)
			{
				if (residue_unres_arr[i] == res_ptr ->GetName())
				{
					vdwradius_unres = sc_vdwradius[i];
				}
			}
			aptr_new -> radius= vdwradius_unres;
			HaColor green_c(2, 200, 200);
			aptr_new -> col= green_c.cidx;
			atm_sc_array[sec_count]=  aptr_new;
			sec_count++;
		}
	}
	pmset->AnnounceGeomChange();
	return TRUE;
//	out_file.close;
//	aptr_new->SetElemNo(DUMMY_ELEM);
}


double HaEmpiricalMod::PenaltySolventAccessible()
{
	if (com_flag)
	{
		FindAxes();
		CenterOfMass();
	}
	std::fstream space_file;// Remove
	space_file.open("sampling_coord_REGxy_2.dat", std::ios::out | std::ios::app);// Remove
	char buf[256];// Remove

	HaAtom* aptr ;
//	HaAtom* aptr1 ;
	MolSet* pmset = GetMolSet();
	HaMolecule* pMol ;
	HaMolecule* pMol1 ;
	HaResidue* res_ptr;
	std::string mol_name;
	std::string mol_name1;
	std::string res_name;
	std::string res_name1;
	int res_num;
	double coord_x, coord_y, coord_z;
	Vec3D prj;
	Vec3D* prj_pnt;
	Vec3D prj1;
	Vec3D* prj_pnt1;
	Vec3D vec ;
	Vec3D vec1 ;
	Vec3D center_helix;
	Vec3D center_mol;
	Vec3D* pnew_vec;
	Vec3D* pnew_vec1;
	Vec3D new_vec;
	Vec3D new_vec1;
	Vec3D zero_pnt;
	
	double angle;
	double dot_prd; 
	int nmol = pmset->GetNMol();
	double cb_x, cb_y, cb_z;
	double cur_sa_ene = 0.0;
	double sa_ene = 0.0 ;
	double face_low_bound = 0.0 ;
	center_mol[0] = 0.0;
	center_mol[1] = 0.0;
	center_mol[2] = 0.0;
	int n = 0;
	

	for(int imol=0 ; imol < nmol; imol++)
	{
		center_helix =  center_arr[imol] ;
		center_mol[0] += center_helix[0] ;
		center_mol[1] += center_helix[1] ;
		center_mol[2] += center_helix[2] ;
	}
	double denum =1/nmol;
	center_mol[0] = center_mol[0]* denum ;
	center_mol[1] = center_mol[1]* denum ;
	center_mol[2] = center_mol[2]* denum ;

	for(int i=0 ; i < n_sa_atoms ; i++)
	{
		aptr = (HaAtom*) atm_solacces[i];
		res_ptr = aptr ->GetHostRes();
		res_num = res_ptr -> GetSerNo();
		res_name = res_ptr -> GetName();
		pMol = aptr ->GetHostMol();
		mol_name = pMol -> GetRef();
		
		for(int imol=0 ; imol < nmol; imol++)
		{
			pMol1 = pmset->GetMolByIdx(imol);
			mol_name1 = pMol1 -> GetRef();
			if(mol_name == mol_name1)
			{
				vec = axis_arr[imol];
				center_helix =  center_arr[imol] ;
			}
		}

		cb_x = (aptr->GetX()) - center_helix[0];
		cb_y = (aptr->GetY()) - center_helix[1];
		cb_z = (aptr->GetZ()) - center_helix[2];
		dot_prd = cb_x*vec[0] + cb_y*vec[1] + cb_z*vec[2];

		sprintf(buf,"%d %2.1f %2.1f %2.1f ", i, center_helix[0],center_helix[1], center_helix[2]); //REMOVE
		space_file<< buf; //REMOVE
//PrintLog("Atom             %3.2f %3.2f %3.2f \n", cb_x,cb_y,cb_z);
		new_vec.SetX(cb_x);
		new_vec.SetY(cb_y);
		new_vec.SetZ(cb_z);
		pnew_vec = &new_vec; // CB atom 

		coord_x = dot_prd * vec[0];
		coord_y = dot_prd * vec[1];
		coord_z = dot_prd * vec[2];
		
		prj.SetX(coord_x*ANG_TO_BOHR);
		prj.SetY(coord_y*ANG_TO_BOHR);
		prj.SetZ(coord_z*ANG_TO_BOHR);
		prj_pnt = &prj;    // projection of atom on helix axis
//PrintLog("Atom on helix    %3.2f %3.2f %3.2f \n", prj[0],prj[1],prj[2]);
		
		cb_x = center_mol[0]- center_helix[0];
		cb_y = center_mol[1]- center_helix[1];
		cb_z = center_mol[2]- center_helix[2];
		dot_prd = cb_x*vec[0] + cb_y*vec[1] + cb_z*vec[2];

		coord_x = dot_prd * vec[0];
		coord_y = dot_prd * vec[1];
		coord_z = dot_prd * vec[2];
		
		prj1.SetX(coord_x*ANG_TO_BOHR);
		prj1.SetY(coord_y*ANG_TO_BOHR);
		prj1.SetZ(coord_z*ANG_TO_BOHR);
		prj_pnt1 = &prj1; // projection of center of mol on helix axis
//PrintLog("Cmol on helix axis %3.2f %3.2f %3.2f \n", prj1[0],prj1[1],prj1[2]);

		new_vec1.SetX(cb_x*ANG_TO_BOHR);
		new_vec1.SetY(cb_y*ANG_TO_BOHR);
		new_vec1.SetZ(cb_z*ANG_TO_BOHR);
//PrintLog("Center molecule    %3.2f %3.2f %3.2f \n", new_vec1[0],new_vec1[1],new_vec1[2]);
		pnew_vec1 = &new_vec1; // center of molset
		double average_face = face_up_bound/2;
		double average_face_minus = -average_face;
		angle = Vec3D::CalcTorsion(pnew_vec, prj_pnt, prj_pnt1, pnew_vec1 )*RAD_TO_DEG;
		if (angle > 0)
		{
			cur_sa_ene += SqrPotential(angle, average_face, average_face, weight_sa);
//			PrintLog("angle %3.3f pMol %s cur_sa_ene %3.3f\n", angle, (pMol -> GetRef()).c_str(),cur_sa_ene );
		}	
		else
		{
			//double low_bound_new = 0 - face_up_bound ;
			//double up_bound_new = 0 -  face_low_bound; 
			cur_sa_ene += SqrPotential(angle, average_face_minus, average_face, weight_sa);
//			PrintLog("angle %3.3f pMol %s cur_sa_ene %3.3f\n", angle, (pMol -> GetRef()).c_str(),cur_sa_ene );
		}
		sprintf(buf,"%2.1f ", angle); //REMOVE
		space_file<< buf; //REMOVE
	} 
	space_file<< buf << std::endl; //REMOVE
	sa_ene = cur_sa_ene;
	space_file.close(); //REMOVE
	return sa_ene;
}

Vec3D
HaEmpiricalMod::FindCentralAxis()
{
	if (com_flag)
	{
		FindAxes();
	}
	
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	Vec3D vec1 ;

	central_axis[0] = 0;
	central_axis[1] = 0;
	central_axis[2] = 0;

	for(int imol=0 ; imol < nmol; imol++)
	{
		vec1 = axis_arr[imol];
		central_axis[0] += vec1[0];
		central_axis[1] += vec1[1];
		central_axis[2] += vec1[2];
	}
	central_axis[0] = central_axis[0]/nmol;
	central_axis[1] = central_axis[1]/nmol;
	central_axis[2] = central_axis[2]/nmol;
	return central_axis ;
}



int
HaEmpiricalMod::Neighborhood()
{
	
	MolSet* pmset = GetMolSet();
	wxString name, name1 ;
	HaAtom* aptr;
	HaAtom* aptr1;
	HaAtom* aptr_ca;
	HaMolecule* pMol;
	HaMolecule* pMol1;
	HaResidue* res_ptr;
	HaResidue* res_ptr1;
	int nmol = pmset -> GetNMol();
	neighbor_mat.newsize(nmol,nmol);
	int i, j;

	for( i =0; i < nmol ; i++)
	{
		for( int j =0; j < nmol ; j++)
		{
			neighbor_mat.SetVal_idx0(i,j,0.0);
		}
	}

	for( int imol =0; imol < nmol ; imol++)
	{
		i=0;
		pMol = pmset->GetMolByIdx(imol);
		HaMolecule::ResidueIterator ritr(pMol);
		for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
		{	
			aptr = (HaAtom*) atm_sc_array[i];
			//PrintLog ("i %d\n", i);
			i++;

			for( int jmol = 0; jmol < nmol ; jmol++)
			{
				if (imol != jmol)
				{
					j=0;
					pMol1 = pmset->GetMolByIdx(jmol);
					HaMolecule::ResidueIterator ritr1(pMol1);
					for(res_ptr1 = ritr1.GetFirstRes(); res_ptr1; res_ptr1 = ritr1.GetNextRes())
					{	
						int number_one = first_res_mol[jmol] ;
						aptr1 = (HaAtom*) atm_sc_array[number_one + j];
						aptr_ca = (HaAtom*) atm_ca_array[number_one + j];

						j++;
						double dist = Vec3D::CalcDistance(aptr,aptr1,ANGSTROM_U);
						double dist_ca = Vec3D::CalcDistance(aptr,aptr_ca,ANGSTROM_U);
						if (dist < dist_neighborhood || dist_ca < dist_neighborhood)
						{
							neighbor_mat.SetVal_idx0(imol,jmol,1.0);
							break;
						}
						
					}
				}
			}
			
		}
	}
	
	return TRUE;
}


double
HaEmpiricalMod::CheckNeighbor(int i, int j)
{
	return neighbor_mat.GetVal_idx0(i,j);
}

void
HaEmpiricalMod::ResidueTypeList()
{
	//	vector <wxString> res_type;
	HaMolecule* pMol;
	HaResidue* res_ptr;
	HaAtom* aptr;
	wxString res_name;
	wxString name;
	int res_num, res_num_old;
	std::string fname = "empir_param.dat";
	FILE* finfo = fopen(fname.c_str(),"r");
	if(finfo == NULL)
	{
		ErrorInMod("HaEmpiricalMod::ResidueTypeList",
		"Can not open EMPIR_PARAM.DAT file");

		pairwise_name_arr.Add("ALA") ;
		pairwise_name_arr.Add("ARG") ;
		pairwise_name_arr.Add("ASN") ;
		pairwise_name_arr.Add("ASP") ;
		pairwise_name_arr.Add("CYS") ;
		pairwise_name_arr.Add("GLN");
		pairwise_name_arr.Add("GLU") ;
		pairwise_name_arr.Add("GLY") ;
		pairwise_name_arr.Add("HIS") ;
		pairwise_name_arr.Add("ILE") ;
		pairwise_name_arr.Add("LEU") ;
		pairwise_name_arr.Add("LYS") ;
		pairwise_name_arr.Add("MET") ;
		pairwise_name_arr.Add("PHE") ;
		pairwise_name_arr.Add("PRO") ;
		pairwise_name_arr.Add("SER") ;
		pairwise_name_arr.Add("THR") ;
		pairwise_name_arr.Add("TRP") ;
		pairwise_name_arr.Add("TYR") ;
		pairwise_name_arr.Add("VAL") ;
		residue_arr.Add("ALA") ;
		residue_arr.Add("ARG") ;
		residue_arr.Add("ASN") ;
		residue_arr.Add("ASP") ;
		residue_arr.Add("CYS") ;
		residue_arr.Add("GLN") ;
		residue_arr.Add("GLU") ;
		residue_arr.Add("GLY") ;
		residue_arr.Add("HIS") ;
		residue_arr.Add("ILE") ;
		residue_arr.Add("LEU");
		residue_arr.Add("LYS");
		residue_arr.Add("MET");
		residue_arr.Add("PHE");
		residue_arr.Add("PRO");
		residue_arr.Add("SER");
		residue_arr.Add("THR");
		residue_arr.Add("TRP");
		residue_arr.Add("TYR");
		residue_arr.Add("VAL");
		residue_unres_arr.Add("ALA") ;
		residue_unres_arr.Add("ARG") ;
		residue_unres_arr.Add("ASN") ;
		residue_unres_arr.Add("ASP") ;
		residue_unres_arr.Add("CYS") ;
		residue_unres_arr.Add("GLN") ;
		residue_unres_arr.Add("GLU") ;
		residue_unres_arr.Add("GLY") ;
		residue_unres_arr.Add("HIS") ;
		residue_unres_arr.Add("ILE") ;
		residue_unres_arr.Add("LEU");
		residue_unres_arr.Add("LYS");
		residue_unres_arr.Add("MET");
		residue_unres_arr.Add("PHE");
		residue_unres_arr.Add("PRO");
		residue_unres_arr.Add("SER");
		residue_unres_arr.Add("THR");
		residue_unres_arr.Add("TRP");
		residue_unres_arr.Add("TYR");
		residue_unres_arr.Add("VAL");
	}

	MolSet* pmset = GetMolSet();
	int nmol = pmset -> GetNMol();
	int imol;
	int i;
	int type_number =0;
	int type_number_la =0;
	int j= 0;
	int nres = pmset-> GetNRes();
	int ncount_pair = pairwise_name_arr.GetCount();
	int ncount_res = residue_arr.GetCount();

	first_res_mol.newsize(nmol);
	mol_res_correspond.newsize(nres*nmol);
	mol_res_correspond_la.newsize(nres*nmol);

	for( imol =0; imol < nmol ; imol++)
	{
		first_res_mol[imol] = j; 
		pMol = pmset->GetMolByIdx(imol);
		HaMolecule::ResidueIterator ritr(pMol);
		for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
		{	
			res_num = res_ptr -> GetSerNo();
			for(i=0; i<ncount_pair; i++)
			{
				if(pairwise_name_arr[i] == res_ptr->GetName() ) type_number = i;
			}
			for(i=0; i<ncount_res; i++)
			{
				if(residue_arr[i] == res_ptr->GetName() ) type_number_la = i;
			}
			mol_res_correspond[j] = type_number;
			mol_res_correspond_la[j] = type_number_la ;
			j++;
		}
	}

if(finfo != NULL) fclose(finfo);
}



int
HaEmpiricalMod::LineSegments()
{
	MolSet* pmset = GetMolSet();
	wxString name ;
	HaAtom* aptr;
	HaMolecule* pMol;
	HaResidue* res_ptr;
	int nmol = pmset -> GetNMol();
	segment_vec.resize(8*nmol);
	int i,j,k;
	VecPtr temp_vec;
	temp_vec.resize(4);
	k= 0;
	int ires_count =0;
	
	for( int imol = 0; imol < nmol ; imol++)
	{
		pMol = pmset->GetMolByIdx(imol);
		AtomIteratorMolecule aitr(pMol);
		int if_flag = TRUE;
		i=0;
		j = 0;
		HaMolecule::ResidueIterator ritr(pMol);
		for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
		{	
			aptr = (HaAtom*) atm_ca_array[ires_count];
			if (i == 4) i = 0;
			temp_vec[i] = aptr;
			if(j > 3) if_flag = FALSE;
			if (if_flag)
			{
				segment_vec[k]= aptr;
				k++;
			}
			i++;
			j++;
			ires_count++;
		}
		int m;
		for(m= 0 ; m< 4;m++)
		{
			segment_vec[k] = temp_vec[m]; 
			k++;
		}
	}
	return TRUE;
}

double
HaEmpiricalMod::PenaltyVDW()
{ 
//	fstream dist_file;
//	dist_file.open("dist_vdw.dat",ios::out | ios::app);
	char buf[256];
	MolSet* pmset = GetMolSet();
	HaAtom* aptr;
	int nmol = pmset->GetNMol();
	Vec3DValArray end1_arr;
	Vec3DValArray end2_arr;
	end1_arr.resize(nmol);
	end2_arr.resize(nmol);
	Vec3D cur_com_atom;
	double numer_x = 0.0;
	double numer_y = 0.0;
	double numer_z = 0.0;
	double com_x = 0.0;
	double com_y = 0.0;
	double com_z = 0.0;
	double centr_com_x1 = 0.0;
	double centr_com_y1 = 0.0;
	double centr_com_z1 = 0.0;
	double centr_com_x2 = 0.0;
	double centr_com_y2 = 0.0;
	double centr_com_z2 = 0.0;
	double denumer = 0.0;
	int i;
	double dist, cur_dist_vdw, repulse_vdw_ene;
	double atract_vdw_ene = 0.0;
	double vdw_ene = 0.0;
	int imol,jmol;
	cur_dist_vdw = 7.5;
	double coef1 = 0.3/sqrt((double) nmol);
	double coef2 = 1.0;
 	double coef3 = 1.0/ (double) nmol;

	repulse_vdw_ene = 0.0;
//	HaAtom* aptr1;
//	double dist1 = 0 ;
	wxString name ;
	HaResidue* res;
	wxString resname ;

	for( imol =0; imol < nmol ; imol++)
	{
		int istart = imol * 8;
		numer_x = 0.0;
		numer_y = 0.0;
		numer_z = 0.0;
		com_x = 0.0;
		com_y = 0.0;
		com_z = 0.0;
//		HaAtom* aptr1;
		for(i = istart; i < istart+4 ; i++)
		{
			aptr = (HaAtom*) segment_vec[i];
			numer_x += aptr -> GetX();
			numer_y += aptr -> GetY();
			numer_z += aptr -> GetZ();
//			aptr1 = (HaAtom*) segment_vec[i+4];
//			dist1 = Vec3D::CalcDistance(aptr,aptr1,ANGSTROM_U);
//		PrintLog ("####DISTANCE0000 ### = %2.3f \n", dist1);
		}
		com_x = (numer_x*0.25) ;
//		PrintLog ("####GetX 000 ### = %2.3f imol= %d\n", com_x, imol);
		com_y = (numer_y*0.25) ;
		com_z = (numer_z*0.25) ;

//		centr_com_x1 += com_x;
//		centr_com_y1 += com_y;
//		centr_com_z1 += com_z;

		cur_com_atom.SetX(com_x);
		cur_com_atom.SetY(com_y);
		cur_com_atom.SetZ(com_z);
   	   	end1_arr[imol].SetCoordFrom(cur_com_atom);
		
		numer_x = 0.0;
		numer_y = 0.0;
		numer_z = 0.0;
		com_x = 0.0;
		com_y = 0.0;
		com_z = 0.0;
		for(i = istart+4; i <istart+8 ; i++)
		{
			aptr = (HaAtom*) segment_vec[i];
			numer_x += aptr -> GetX();
			numer_y += aptr -> GetY();
			numer_z += aptr -> GetZ();
		}
		com_x = (numer_x*0.25) ;
		com_y = (numer_y*0.25) ;
		com_z = (numer_z*0.25) ;
//		PrintLog ("####GetX 111 ### = %2.3f imol= %d\n", com_x, imol);
//		centr_com_x2 += com_x;
//		centr_com_y2 += com_y;
//		centr_com_z2 += com_z;

		cur_com_atom.SetX(com_x);
		cur_com_atom.SetY(com_y);
		cur_com_atom.SetZ(com_z);
   	   	end2_arr[imol].SetCoordFrom(cur_com_atom);
		
	}
	double p1[3];
	double p2[3];
	Vec3D temp1;
	Vec3D temp2;
	double p3[3];
	double p4[3];
	Vec3D temp3;
	Vec3D temp4;
	Vec3D temp_atr1;
	Vec3D temp_atr2;

	int imol_one = 0;
	int	imol_two = 0;

	Vec3D* vec1;
	Vec3D* vec2;
	Vec3D vec3;
	double center_attract = 0;

//	temp_atr1[0] = (centr_com_x1 + centr_com_x2)* 0.5 *coef3;
//	temp_atr1[1] = (centr_com_y1 + centr_com_y2)* 0.5 *coef3;
//	temp_atr1[2] = (centr_com_z1 + centr_com_z2)* 0.5 *coef3;
	temp_atr1[0] = 0.0;
	temp_atr1[1] = 0.0;
	temp_atr1[2] = 0.0;

	vec1 = &temp_atr1;

	for(imol =0; imol < nmol ; imol++)
	{
		temp1 = end1_arr[imol];
		temp2 = end2_arr[imol];
		for (i = 0;i<3; i++)
		{
			p1[i] = temp1[i];
			p2[i] = temp2[i];
		}
		
		for (i = 0; i<3; i++)
		{
			temp_atr2[i] = (p1[i]+p2[i])*0.5;
		}

		vec2 = &temp_atr2;
		dist = Vec3D::CalcDistance(vec1,vec2,ANGSTROM_U); // Distance from (0, 0, 0) point
		center_attract +=pow(dist*coef1 - coef2, 4)* weight_vdw;
//		if(imol == 1 ||imol == 3)
//		{
//	sprintf(buf,"p1[0] %2.4f p1[1] %2.4f p1[2] %2.4f p2[0] %2.4f p2[1] %2.4f p0[2] %2.4f", p1[0], p1[1],p1[2], p2[0],p2[1], p2[2]);
//	dist_file<< buf << endl;
//	sprintf(buf,"distAt %2.4f center_attract %2.4f imol %d" , dist,center_attract, imol);
//	dist_file<< buf << endl;
//		}
		for( jmol = imol; jmol < nmol ; jmol++)
		{
			temp3 = end1_arr[jmol];
			temp4 = end2_arr[jmol];
			for (i = 0; i<3; i++)
			{
				p3[i] = temp3[i];
				p4[i] = temp4[i];
			}
			for (i = 0; i<3; i++)
			{
				vec3[i] = (p3[i]+p4[i])*0.5;
			}

			if (jmol != imol) 
			{
				//dist = 1/ Segments_Dist_3D(p1, p2, p3, p4); // Distance between closest points of line segments
				dist = 1/ Vec3D::CalcDistance(&vec3,vec2,ANGSTROM_U); // Distance between centers
				if(dist>0)
				{
					repulse_vdw_ene = pow(cur_dist_vdw*dist,4);
//					atract_vdw_ene = pow(cur_dist_vdw/dist,3);
					// PrintLog("###VDW### dist %2.3f, vdw_ene %2.3f imol_one %d, imol_two %d\n", dist, repulse_vdw_ene-atract_vdw_ene, imol, jmol);
					//vdw_ene +=(repulse_vdw_ene - atract_vdw_ene)*weight_vdw;
					vdw_ene += repulse_vdw_ene*weight_vdw;
					//				sprintf(buf,"dist %2.4f repulse_vdw_ene %2.4f atract_vdw_ene %2.4f imol %d jmol %d ", dist,repulse_vdw_ene, atract_vdw_ene, imol, jmol);
					//				dist_file<< buf << endl;
				}
				else vdw_ene = 10e16;
			}
		}

  //		PrintLog("###VDW### vdw_ene %2.3f center_attract %2.3f imol_one %d\n", vdw_ene, center_attract, imol);
	}

	
	return vdw_ene + center_attract;
//	return center_attract;
}


double 
// HaEmpiricalMod::Segments_Dist_3D ( Vec3D p1, Vec3D p2, Vec3D p3, Vec3D p4)
HaEmpiricalMod::Segments_Dist_3D (double p1[3], double p2[3], double p3[3],   double p4[3] )

//********************************************************************
//
//  Purpose:
//
//    SEGMENTS_DIST_3D computes the distance between two line segments in 3D.
//
//  Discussion:
//
//    A line segment is the portion of an infinite line that lies between
//    two given points.  The behavior of the distance function is a bit
//    complicated.  
//
//  Modified:
//
//    03 November 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double P1[3], P2[3], the endpoints of the first segment.
//
//    Input, double P3[3], P4[3], the endpoints of the second segment.
// 
//    Output, double SEGMENTS_DIST_3D, the distance between the line segments.
//
{


  double d1;
  double d2;
  double dist;
  double dl;
  double dm;
  double dr;
//  Vec3D pn1;
//  Vec3D pn2;
//  Vec3D pt;
  double pn1[3]; 
  double pn2[3];
  double pt[3];
  bool result;
  double t1;
  double t2;
  double tl;
  double tm;
  double tmin;
  double tr;
//
//  Find the nearest points on line 2 to the endpoints of line 1.
//
  Segment_Point_Near_3D ( p3, p4, p1, pn1, &d1, &t1 );
  Segment_Point_Near_3D ( p3, p4, p2, pn2, &d2, &t2 );
  if ( t1 == t2 )
  {
//	PrintLog("t1 %2.2f t2 %2.2f d1 %2.2f d2 %2.2f\n", t1, t2, d1, d2);
//	PrintLog("PN %2.2f %2.2f %2.2f \n", pn1[0],pn1[1],pn1[2]);
    dist = Segment_Point_Dist_3D ( p1, p2, pn1 );
    return dist;
  }
//
//  On line 2, over the interval between the points nearest to line 1, 
//  the square of the distance of any point to line 1 is a quadratic function.  
//  Evaluate it at three points, and seek its local minimum.
//
  dl = Segment_Point_Dist_3D ( p1, p2, pn1 );

  pt[0] = 0.5 * ( pn1[0] + pn2[0] );
  pt[1] = 0.5 * ( pn1[1] + pn2[1] );
  pt[2] = 0.5 * ( pn1[2] + pn2[2] );

  dm = Segment_Point_Dist_3D ( p1, p2, pt );

  dr = Segment_Point_Dist_3D ( p1, p2, pn2 );

  tl = 0.0;
  tm = 0.5;
  tr = 1.0;

  dl = dl * dl;
  dm = dm * dm;
  dr = dr * dr;

  result = Minquad ( tl, dl, tm, dm, tr, dr, &tmin, &dist );

  if ( !result )
  {
		ErrorInMod("HaEmpiricalMod::Segments_Dist_3D",
			       "ERROR. Can not calculate distance");
		return dist= 0.0;
  }

  dist = sqrt ( dist );
  return dist;
}

void 
//HaEmpiricalMod::Segment_Point_Near_3D ( Vec3D p1, Vec3D p2, Vec3D p,
//  Vec3D pn, double *dist, double *t )
HaEmpiricalMod::Segment_Point_Near_3D ( double p1[3], double p2[3], double p[3],
  double pn[3], double *dist, double *t )

//********************************************************************
//
//  Purpose:
//
//    SEGMENT_POINT_NEAR_3D finds the point on a line segment nearest a point in 3D.
//
//  Modified:
//
//    04 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double P1[3], P2[3], the two endpoints of the line segment.
//
//    Input, double P[3], the point whose nearest neighbor
//    on the line segment is to be determined.
//
//    Output, double PN[3], the point on the line segment which is nearest to P.
// 
//    Output, double *DIST, the distance from the point to the nearest point
//    on the line segment.
//
//    Output, double *T, the relative position of the nearest point
//    PN to the defining points P1 and P2.
//
//      PN = (1-T)*P1 + T*P2.
//
//    T will always be between 0 and 1.
//
//
{

  double bot;

  if ( Dvec_Eq ( 3, p1, p2 ) )
  {
    *t = 0.0;
    Dvec_Copy ( 3, p1, pn );
  }
  else
  {

    bot = 
        ( p1[0] - p2[0] ) * ( p1[0] - p2[0] )
      + ( p1[1] - p2[1] ) * ( p1[1] - p2[1] )
      + ( p1[2] - p2[2] ) * ( p1[2] - p2[2] );

    *t = (
      + ( p1[0] - p[0] ) * ( p1[0] - p2[0] )
      + ( p1[1] - p[1] ) * ( p1[1] - p2[1] )
      + ( p1[2] - p[2] ) * ( p1[2] - p2[2] ) ) / bot;

    if ( *t < 0.0 )
    {
      *t = 0.0;
      Dvec_Copy ( 3, p1, pn );
    }
    else if ( 1.0 < *t )
    {
      *t = 1.0;
      Dvec_Copy ( 3, p2, pn );
    }
    else
    {
      pn[0] = p1[0] + *t * ( p2[0] - p1[0] );
      pn[1] = p1[1] + *t * ( p2[1] - p1[1] );
      pn[2] = p1[2] + *t * ( p2[2] - p1[2] );
    }
  }
  *dist = sqrt ( 
      ( pn[0] - p[0] ) * ( pn[0] - p[0] ) 
    + ( pn[1] - p[1] ) * ( pn[1] - p[1] ) 
    + ( pn[2] - p[2] ) * ( pn[2] - p[2] ) );
//PrintLog("$$$Dist %2.2f \n",&dist);
  return;
}


double 
//HaEmpiricalMod::Segment_Point_Dist_3D ( Vec3D p1, Vec3D p2, Vec3D p )
HaEmpiricalMod::Segment_Point_Dist_3D ( double p1[3], double  p2[3], double  p[3] )

//********************************************************************
//
//  Purpose:
//
//    SEGMENT_POINT_DIST_3D computes the distance from a point to a line segment in 3D.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double P1[3], P2[3], the endpoints of the line segment.
//
//    Input, double P[3], the point whose nearest neighbor on the line 
//    segment is to be determined.
//
//    Output, double SEGMENT_POINT_DIST_3D, the distance from the point to the line segment.
//
{


  double bot;
  double dist;
  double t;
//  Vec3D pn;
  double pn[3];
//
//  If the line segment is actually a point, then the answer is easy.
//
  if ( Dvec_Eq ( 3, p1, p2 ) )
  {
    Dvec_Copy ( 3, p1, pn );
  }
  else
  {
    bot = pow ( p2[0] - p1[0], 2 ) 
        + pow ( p2[1] - p1[1], 2 )
        + pow ( p2[2] - p1[2], 2 );

	t = (
		(p1[0] - p[0])*(p1[0] - p2[0])
		+(p1[1] - p[1])*(p1[1] - p2[1])
		+(p1[2] - p[2])*(p1[2] - p2[2]) )/bot;

    t = D_Max ( t, 0.0 );
    t = D_Min ( t, 1.0 );

    pn[0] = p1[0] + t * ( p2[0] - p1[0] );
    pn[1] = p1[1] + t * ( p2[1] - p1[1] );
    pn[2] = p1[2] + t * ( p2[2] - p1[2] );
//	PrintLog("PN0 %2.2f %2.2f %2.2f \n", pn[0],pn[1],pn[2]);
  }
 
  dist = sqrt (
      pow ( pn[0] - p[0], 2 )
    + pow ( pn[1] - p[1], 2 )
    + pow ( pn[2] - p[2], 2 ) );
//PrintLog("DIST in SEGMENT_POINT_DIST_3D %2.2f \n",dist);
//PrintLog("P1 %2.2f %2.2f %2.2f \n", p1[0],p1[1],p1[2]);
//PrintLog("P2 %2.2f %2.2f %2.2f \n", p2[0],p2[1],p2[2]);
//PrintLog("P %2.2f %2.2f %2.2f \n", p[0],p[1],p[2]);
//PrintLog("PN %2.2f %2.2f %2.2f \n", pn[0],pn[1],pn[2]);

  return dist;
}


bool 
//HaEmpiricalMod::Dvec_Eq ( int n, Vec3D  a1, Vec3D a2 )
HaEmpiricalMod::Dvec_Eq ( int n, double  a1[], double  a2[] )

//******************************************************************************
//
//  Purpose:
//
//    DVEC_EQ is true if every pair of entries in two vectors is equal.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], two vectors to Compare.
//
//    Output, bool DVEC_EQ.
//    DVEC_EQ is TRUE if every pair of elements A1(I) and A2(I) are equal,
//    and FALSE otherwise.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] != a2[i] )
    {
      return false;
    }
  }
  return true;
}

void 
//HaEmpiricalMod::Dvec_Copy ( int n, Vec3D  a1, Vec3D a2 )
HaEmpiricalMod::Dvec_Copy ( int n, double  a1[], double a2[] )
//******************************************************************************
//
//  Purpose:
//
//    DVEC_COPY copies a real vector.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Input, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}

double 
HaEmpiricalMod::D_Max ( double x, double y )

//*********************************************************************
//
//  Purpose:
//
//    D_MAX returns the maximum of two real values.
//
//  Modified:
//
//    10 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double D_MAX, the maximum of X and Y.
//
{
  if ( y < x )
  {
    return x;
  } 
  else
  {
    return y;
  }
}
//*********************************************************************

double 
HaEmpiricalMod::D_Min ( double x, double y )

//*********************************************************************
//
//  Purpose:
//
//    D_MIN returns the minimum of two real values.
//
//  Modified:
//
//    09 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double D_MIN, the minimum of X and Y.
//
{
  if ( y < x )
  {
    return y;
  } 
  else
  {
    return x;
  }
}

bool 
HaEmpiricalMod::Minquad ( double x1, double y1, double x2, double y2, double x3, double y3, 
  double *xmin, double *ymin )

//********************************************************************
//
//  Purpose:
//
//    MINQUAD finds a local minimum of F(X) = A * X**2 + B * X + C.
//
//  Discussion:
//
//    MINQUAD is primarily intended as a utility routine for use by
//    DISLSLS3.  The square of the distance function between a point
//    and a line segment has the form of F(X).  Hence, we can seek
//    the line on the second segment which minimizes the square of
//    the distance to the other line segment.
//
//  Modified:
//
//    02 November 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X1, Y1, X2, Y2, X3, Y3, are three sets of data
//    of the form ( X, F(X) ).  The three X values must be distinct.
//
//    Output, double *XMIN, *YMIN.  XMIN is a point within the interval
//    spanned by X1, X2 and X3, at which F takes its local minimum
//    value YMIN.
//
//    Output, bool MINQUAD, 
//    true if no error, 
//    false if error because X values are not distinct.
//
{
  int ierror;
  double x;
  double xleft;
  double xrite;
  double y;

  *xmin = 0.0;
  *ymin = 0.0;
//
//  Refuse to deal with coincident data.
//
  if ( x1 == x2 || x2 == x3 || x3 == x1 )
  {
    return false;
  }
//
//   the interval endpoints.
//
  xleft = x1;
  if ( x2 < xleft )
  {
    xleft = x2;
  }
  if ( x3 < xleft )
  {
    xleft = x3;
  }
  xrite = x1;
  if ( xrite < x2 )
  {
    xrite = x2;
  }
  if ( xrite < x3 )
  {
    xrite = x3;
  }
//
//  Find the minimizer and its function value over the three input points.
//
  if ( y1 <= y2 && y1 <= y3 )
  {
    *xmin = x1;
    *ymin = y1;
  }
  else if ( y2 <= y1 && y2 <= y3 )
  {
    *xmin = x2;
    *ymin = y2;
  }
  else if ( y3 <= y1 && y3 <= y2 )
  {
    *xmin = x3;
    *ymin = y3;
  }
//
//  Find the minimizer and its function value over the real line.
//
  ierror = Parabola_Ex ( x1, y1, x2, y2, x3, y3, &x, &y );

  if ( ierror != 2 && y < *ymin && xleft < x && x < xrite )
  {
    *xmin = x;
    *ymin = y;
  }

  return true;
}


int 
HaEmpiricalMod::Parabola_Ex ( double x1, double y1, double x2, double y2, double x3, 
  double y3, double *x, double *y )

//********************************************************************
//
//  Purpose:
//
//    PARABOLA_EX finds the extremal point of a parabola determined by three points.
//
//  Modified:
//
//    17 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X1, Y1, X2, Y2, X3, Y3, the coordinates of three points
//    on the parabola.  X1, X2 and X3 must be distinct.
//
//    Output, double *X, *Y, the X coordinate of the extremal point of the
//    parabola, and the value of the parabola at that point.
//
//    Output, int PARABOLA_EX, error flag.
//    0, no error.
//    1, two of the X values are equal.
//    2, the data lies on a straight line; there is no finite extremal
//    point.
//    3, the data lies on a horizontal line; every point is "extremal".
//
{
  double bot;

  *x = 0.0;
  *y = 0.0;

  if ( x1 == x2 || x2 == x3 || x3 == x1 )
  {
    return 1;
  }

  if ( y1 == y2 && y2 == y3 && y3 == y1 )
  {
    *x = x1;
    *y = y1;
    return 3;
  }

  bot = ( x2 - x3 ) * y1 - ( x1 - x3 ) * y2 + ( x1 - x2 ) * y3;

  if ( bot == 0.0 )
  {
    return 2;
  }

  *x = 0.5 * ( 
      x1 * x1 * ( y3 - y2 )
    + x2 * x2 * ( y1 - y3 )
    + x3 * x3 * ( y2 - y1 ) ) / bot;

  *y =  (
      ( *x - x2 ) * ( *x - x3 ) * ( x2 - x3 ) * y1
    - ( *x - x1 ) * ( *x - x3 ) * ( x1 - x3 ) * y2
    + ( *x - x1 ) * ( *x - x2 ) * ( x1 - x2 ) * y3 ) /
    ( ( x1 - x2 ) * ( x2 - x3 ) * ( x1 - x3 ) );

  return 0;
}
double HaEmpiricalMod::MinEnergy() // added by jose 06-04-08
{
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	HaMolecule* pMol;
	HaMolecule* pMol1;
	HaResidue* res_ptr;
	HaResidue* res_ptr1;
	//double Ma_dist[8][8];
	double Ma_dist[8][8]={
	11.9353,  8.14251,  8.62619,  11.3768,  9.73541,  7.40371,  10.8149,  12.6421,
	13.4311,  9.9104,  8.92935,  11.4249,  10.2228,  6.80975,  9.31407,  11.6569,
	14.5244,  11.1717,  9.84305,  13.0715,  12.736,  9.53072,  11.5909,  14.5456,
	11.2351,  8.08982,  7.00536,  10.6476,  10.7761,  8.37356,  10.634, 13.5925,
	9.96735,  7.19219,  4.75002,  7.81988,  8.28271,  5.45053,  6.96824,  10.2211,
	13.3565,  10.9006,  8.01736,  10.7489,  11.6604,  8.47903,  8.85643,  12.4801,
	13.1786,  11.0912,  8.47234,  11.7686,  13.2567,  10.7649,  11.3849,  15.0237,
	10.1269,  8.87282,  5.78664,  8.73555,  10.99,  9.34263, 9.46895,  12.9722,
	};
	HaAtom* aptr;
	HaAtom* aptr1;
	int imol, jmol;
	int i=0,j=0;
	double weigth = 0.005;
	//fstream constraints_file;
	//constraints_file.open("constraints.dat", ios::out | ios::app);
	std::fstream energy_file;
	energy_file.open("energy_file.dat", std::ios::out | std::ios::app);
	double ene = 0.0;
	double energy_total = 0.0, dist;

 	for (imol= 0; imol < nmol; imol++)
	{
		pMol = pmset -> GetMolByIdx(imol);
		HaMolecule::ResidueIterator ritr(pMol);
		for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
		{
			aptr = (HaAtom*) atm_sc_array[i];
			//aptrca = (HaAtom*) atm_ca_array[i];
			for (jmol= imol+1; jmol < nmol;jmol++)
			{
				int number_shift = first_res_mol[jmol];
				pMol1 = pmset -> GetMolByIdx(jmol);
				HaMolecule::ResidueIterator ritr1(pMol1);
				j = 0;
				for(res_ptr1 = ritr1.GetFirstRes(); res_ptr1; res_ptr1 = ritr1.GetNextRes())
				{
					aptr1 = (HaAtom*) atm_sc_array[j+ number_shift ];
					//aptrca1= (HaAtom*) atm_ca_array[j+ number_shift ];
					double dist_sc = Vec3D::CalcDistance(aptr, aptr1, ANGSTROM_U);
					//Ma_dist[i][j]=dist_sc;
					dist =Ma_dist[i][j];  
					ene +=(dist-dist_sc)*(dist-dist_sc)*weigth;
					//PrintLog("ene %f\n",ene);
					j++;
				}
			}
			//constraints_file<<Ma_dist[i][0]<<"  "<<Ma_dist[i][1]<<"  "<<Ma_dist[i][2]<<"  "<<Ma_dist[i][3]<<"  "<<Ma_dist[i][4]<<"  "<<Ma_dist[i][5]<<"  "<<Ma_dist[i][6]<<"  "<<Ma_dist[i][7]<<endl;
			i++;
		} 
	}
	energy_file<<ene<< std::endl;
	return ene;
}
double HaEmpiricalMod::ToyEnergy() // added by jose 05-21-08 it will be deleted soon
{
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	HaMolecule* pMol;
	HaMolecule* pMol1;
	HaAtom* aptr;
	HaAtom* aptr1;
	HaAtom* aptrca;
	HaAtom* aptrca1;
	HaResidue* res_ptr;
	HaResidue* res_ptr1;
	double energy_final = 0.0;
	int i, j;
	int imol=0, jmol =0;
	char name1[16][16] = {"$_3$GLY8:A.CA","$_2$GLY1:A.CA","$_3$GLY1:A.CA","$_2$GLY8:A.CA","$_2$GLY8:A.CA","$_3$GLY8:A.CA","$_2$GLY1:A.CA","$_3$GLY4:A.CA","$_2$GLY3:A.CA","$_3$GLY4:A.CA","$_2$GLY6:A.CA","$_3$GLY7:A.CA","$_2$GLY3:A.CA","$_3$GLY1:A.CA","$_2$GLY6:A.CA","$_3$GLY4:A.CA"};
	//{"$_3$GLY6:A.CA","$_2$GLY2:A.CA","$_3$GLY3:A.CA","$_2$GLY8:A.CA","$_3$GLY4:A.CA","$_2$GLY6:A.CA"};
	Vec3DValArray pos_arr1;
	pos_arr1.resize(16); //number of atoms interacting
	Vec3DValArray tor_arr1;
	tor_arr1.resize(8);
	Vec3D cur_atom;
	Vec3D vec1 ;
	Vec3D vec2 ;
	double distance[8];
	double sigma =0.0, sigma1= 0.0, sigma_ca = 0.0, sigma_ca1=0.0;
	double sigma2 = 0.0, sigma_ca2 = 0.0, sigma_sc_ca = 0.0;
	double dist_ca, dist_sc, dist_sc_ca;
	double temp_x = 0;
	double temp_y = 0;
	double temp_z = 0;
	double weight = 1;
	double dist;
	double ene = 0.0;
	double ene_tor= 0.0;
	double vdw_ene = 0.0;
	double ene_val = 0.0;
	std::fstream harmonic_file;
	harmonic_file.open("harmonic_file.dat", std::ios::out | std::ios::app);
	std::fstream distance_constraint;
	distance_constraint.open("distance_constraint.dat", std::ios::out | std::ios::app);
	char buf[256];
	for (i=0; i <16;i++)
	{
		aptr = pmset->GetAtomByRef(name1[i]);
		temp_x = aptr->GetX();
		temp_y = aptr->GetY();
		temp_z = aptr->GetZ();
		cur_atom.SetX(temp_x);
		cur_atom.SetY(temp_y);
		cur_atom.SetZ(temp_z);
		pos_arr1[i].SetCoordFrom(cur_atom);
//		PrintLog("x %2.1f, y %2.1f, z %2.1f\n", temp_x, temp_y, temp_z); 
	}
	int k=0;
	for (i=0; i<16; i=i+2)
	{
		vec1 = pos_arr1[i];
		vec2 = pos_arr1[i+1];
//		PrintLog("i=%d\n",i);
		double dist = Vec3D::CalcDistance(&vec1,&vec2,ANGSTROM_U);
		if (i==0){
			ene = (dist - 12.642)*(dist - 12.642)*weight;//(dist- 6.8097)*(dist- 6.8097)*weight;
			//PrintLog("i %d ene %f dist=%f\n", i, ene, dist);
			distance[k] = dist;
			tor_arr1[k].SetCoordFrom(vec1);
			k++;
		}
		if (i==2){
			ene = (dist - 10.127)*(dist - 10.127)*weight;//(dist- 5.7866)*(dist- 5.7866)*weight;
			//PrintLog("i %d ene %f dist=%f\n", i, ene, dist);
		    distance[k] = dist;
			tor_arr1[k].SetCoordFrom(vec1);
			k++;}
		if (i==4){
			ene = (dist-12.972)*(dist-12.972)*weight;
			//PrintLog("i %d ene %f dist=%f\n", i, ene, dist);
		    distance[k] = dist;
			tor_arr1[k].SetCoordFrom(vec1);
			k++;}
		if (i==6){
			ene = (dist-11.377)*(dist-11.377)*weight;
			distance[k] = dist;
			tor_arr1[k].SetCoordFrom(vec1);
			k++;
		}
		if (i==8){
			ene = (dist-13.701)*(dist-13.701)*weight;
			distance[k] = dist;
			//PrintLog("distance %f", distance[k]);
			tor_arr1[k].SetCoordFrom(vec1);
			k++;
		}
		if (i==10)
		{ ene = (dist-8.856)*(dist-8.856)*weight;
		  distance[k] = dist;
		  tor_arr1[k].SetCoordFrom(vec1);
		  k++;
		}
		if (i==12)
		{
			ene = (dist-14.524)*(dist-14.524)*weight;
			distance[k] = dist;
			tor_arr1[k].SetCoordFrom(vec1);
			k++;
		}
		if (i==14)
		{
			ene = (dist-10.749)*(dist-10.749)*weight;
			distance[k] = dist;
			tor_arr1[k].SetCoordFrom(vec1);
			k++;
		}
		energy_final += ene;
//		PrintLog("X=%2.1f, Y=%2.1f, Z=%2.1f\n", vec1->GetX(), vec1->GetY(), vec1->GetZ());
//		PrintLog("Dist %2.1f, Energy %2.1f\n", dist, energy_final);
	}
	//double tor_ang = Vec3D::CalcTorsion(&tor_arr1[0], &tor_arr1[1], &tor_arr1[2], &tor_arr1[3]);
	//double val_ang1 = Vec3D::CalcAngle(&tor_arr1[0], &tor_arr1[1], &tor_arr1[2]);
	//double val_ang2 = Vec3D::CalcAngle(&tor_arr1[1], &tor_arr1[2], &tor_arr1[3]);
	//ene_tor = (1.25181-tor_ang)*(1.25181-tor_ang)*50;
	//ene_val = (1.35400-val_ang1)*(1.35400-val_ang1)*50 + (1.22843-val_ang2)*(1.22843-val_ang2)*50;
	//PrintLog("ener_tor= %f", ene_tor);
	distance_constraint << distance[0] << "  "<< distance[1] <<"  "<<distance[2] <<"  "<< distance[3]<<"  "<< distance[4]<<"  "<<distance[5]<<"  "<< distance[6]<<"  "<<distance[7]<< std::endl;

// clash energy & probably in the future pairwise 
 	for (imol= 0; imol < nmol; imol++)
	{
		i=0;
		pMol = pmset -> GetMolByIdx(imol);
		HaMolecule::ResidueIterator ritr(pMol);
		for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
		{
			aptr = (HaAtom*) atm_sc_array[i];
			aptrca = (HaAtom*) atm_ca_array[i];
			sigma = aptr->radius;
			sigma_ca = aptrca->radius;
			//PrintLog("sigma %3.3f\n", sigma);
			i++;
			for (jmol= imol+1; jmol < nmol;jmol++)
			{
				int number_shift = first_res_mol[jmol];
				pMol1 = pmset -> GetMolByIdx(jmol);
				HaMolecule::ResidueIterator ritr1(pMol1);
				j = 0;
				for(res_ptr1 = ritr1.GetFirstRes(); res_ptr1; res_ptr1 = ritr1.GetNextRes())
				{
					aptr1 = (HaAtom*) atm_sc_array[j+ number_shift ];
					aptrca1= (HaAtom*) atm_ca_array[j+ number_shift ];
					sigma1 = aptr1->radius;
					sigma_ca1 = aptrca1->radius;
					//PrintLog("sigma_ca1 %3.3f\n",sigma_ca1);
					//PrintLog("sigma_ca %3.3f\n",sigma_ca);
					// according to pairwise function
					sigma2 = (sigma+sigma1); //soft interaction
					sigma_ca2 = (sigma_ca+sigma_ca1); // hard interaction
					//PrintLog("sigma_ca2 %3.3f\n",sigma_ca2);
					sigma_sc_ca = (sigma_ca1+sigma); 
					double dist_sc = Vec3D::CalcDistance(aptr, aptr1, ANGSTROM_U);
					double dist_ca = Vec3D::CalcDistance(aptrca, aptrca1, ANGSTROM_U);
					double dist_sc_ca = Vec3D::CalcDistance(aptr, aptrca1, ANGSTROM_U);
					if (dist_sc_ca < sigma_sc_ca)
					{
						vdw_ene +=(sigma_sc_ca-dist_sc_ca)*(sigma_sc_ca-dist_sc_ca)*50;
					}
					if (dist_ca < sigma_ca2)
					{
						vdw_ene +=(sigma_ca2-dist_ca)*(sigma_ca2-dist_ca)*50;
						//PrintLog("dist_ca %3.3f vdw_ene %3.3f\n",dist_ca, vdw_ene);
					}
					if (dist_sc < sigma2)
					{
						vdw_ene +=(sigma2-dist_sc)*(sigma2-dist_sc)*1;
					}
					j++;
				}
			}
		} 

	}
	//sprintf(buf,"Eharm %3.2f Evdw %3.2f", energy_final, vdw_ene);
	//harmonic_file << buf << endl;
	harmonic_file << energy_final << "  "<< vdw_ene <<"  "<<ene_tor<<"  "<<ene_val<< std::endl;

	return energy_final+vdw_ene;//+ene_tor+ene_val;
}

double
HaEmpiricalMod::HarmonicEnergy()
{
	MolSet* pmset = GetMolSet();
   	int nmol = pmset->GetNMol();
	HaMolecule* pMol;
	HaAtom* aptr ;
	HaVec_double ene(nmol);
	double dist;
	int imol = 0;
	int j=0;
	Vec3DValArray center_arr1;
	center_arr1.resize(nmol);
	Vec3D cur_com_atom;
	Vec3D* vec1 ;
	Vec3D* vec2 ;
	std::fstream harmonic_file;
	harmonic_file.open("harmonic_file.dat", std::ios::out | std::ios::app);
	char buf[256];
//	HaEmpiricalMod* emp_mod = pmset->GetEmpiricalMod(true);
//	center_arr1 = emp_mod->CenterOfMass();
	for( imol =0; imol < nmol; imol++)
	{	
		double x_aver =0 ;
		double y_aver = 0;
		double z_aver = 0;
		int count= 0;
		pMol = pmset->GetMolByIdx(imol);
		AtomIteratorMolecule aitr(pMol);
		for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
		{
			double x_coor = aptr->GetX();
			double y_coor = aptr->GetY();
			double z_coor = aptr->GetZ();
			x_aver += x_coor ;
			y_aver += y_coor ;
			z_aver += z_coor ;
			count++;
		}
		cur_com_atom.SetX(x_aver/count);
		cur_com_atom.SetY(y_aver/count);
		cur_com_atom.SetZ(z_aver/count);
   	   	center_arr1[imol].SetCoordFrom(cur_com_atom);
	}

	for( imol =0; imol < nmol-1; imol++)
	{
		vec1 = &center_arr1[imol];
		vec2 = &center_arr1[imol+1];
		double dist = Vec3D::CalcDistance(vec1,vec2,ANGSTROM_U);
		ene[imol] = (dist- 6.0)*(dist- 6.0);
//		PrintLog("DISt %f, vec1[0]=%2.1f,vec2[0]=%2.1f\n", dist, vec1[0], vec2[0]);
	} 
/*		vec1 = &center_arr1[0];
		vec2 = &center_arr1[nmol-1];
		dist = Vec3D::CalcDistance(vec1,vec2,ANGSTROM_U);
		ene[nmol-1] = (dist- 6.0)*(dist- 6.0); not valid for 2 molecules*/
		double ene_final = 0;
	for( imol =0; imol < nmol; imol++)
	{
		ene_final = ene_final+ene[imol];
	}
/*commented by Jose
	double ene_final =0.0;
	double ene_r = 0.0;
	for( imol = 0; imol < nmol-1; imol++)
	{
		vec1 = &center_arr1[imol];
		for( j = imol+1; j < nmol-1; j++)
			vec2 = &center_arr1[imol+1];
			double dist = Vec3D::CalcDistance(vec1,vec2,ANGSTROM_U);
			ene_r = (dist- 6.0)*(dist- 6.0);
			ene_final += ene_r;
			PrintLog("DISt %2.1f, vec1[0]=%2.1f, vec2[0]=%2.1f, ene=%2.1f\n", dist, vec1[0], vec2[0], ene_r);
	}*/
	sprintf(buf, "%2.1f", ene_final);
	harmonic_file << buf << std::endl;

	return ene_final;
}

double
HaEmpiricalMod::PenaltyCentralAttract()
{
	MolSet*  pmset =  GetMolSet();
	int nmol = pmset->GetNMol();
	int i;
	int imol;
	double weight_cntr_attract = 0.1;
	if (com_flag)
	{
//		CenterOfMass();
		GetGeomCenterToy();

	}
	Vec3D temp2;

	double center_attract = 0;
	double maximum_dimension_2 = maximum_dimension*maximum_dimension;

	for(imol =0; imol < nmol ; imol++)
	{
		//temp2 = center_arr[imol];
		temp2 = center_arr_t[imol];
		double	dist_2 = (temp2[0]*temp2[0] + temp2[1]*temp2[1] + temp2[2]*temp2[2]);
		if (dist_2 > maximum_dimension_2){
			center_attract += pow(maximum_dimension_2 - dist_2, 2)* weight_cntr_attract; //modified because this is not clear how is going to minimize the energy
//			center_attract += pow(maximum_dimension-sqrt(dist_2),2)*weight_cntr_attract;
		    //PrintLog("maximum_dim %2.2f center_attract %2.2f dist %2.2f\n",maximum_dimension, center_attract, sqrt(dist_2));
		}
	}
//	PrintLog("Center energy: %f\n",center_attract);
	return center_attract;
}

double
HaEmpiricalMod::PenaltyHelicePack() // added by jose
{
	// pack helices
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	int a = GetGeomCenter();
	double ene_pack = 0.0;
	double weigth = 0.1;
	double maximum_dimension_2 = maximum_dimension*maximum_dimension;
	double dist = 0.0, dist_2 = 0.0;
	//double ri = 1/(maximum_dimension_2*maximum_dimension_2*maximum_dimension_2);
	Vec3D* vec1;
	int imol = 0;
	for (imol=0; imol<nmol;imol++)
	{
		vec1 = &center_arr[imol];
		double dist = Vec3D::CalcDistance(&geom_center_t,vec1,ANGSTROM_U);
		//dist_2 = vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2];
		dist_2 = dist*dist;
		if (dist_2 >= maximum_dimension_2)
			//ene_pack += 1/(dist_2*dist_2*dist_2)-ri;
			ene_pack += (dist_2-maximum_dimension_2)*(dist_2-maximum_dimension_2)*weigth;
		//PrintLog("ene %2.1f\n", ene_pack);
	}
	return ene_pack;
}

int
HaEmpiricalMod::CalcForceCentralAttract(Vec3DValArray& force_cntl_array)
{
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	int i;
	int imol;
	double weight_cntr_attract = 0.1;
	if (com_flag)
	{
		CenterOfMass();
	}

	force_cntl_array.resize(nmol);

	Vec3D temp1;
	Vec3D temp2;
	Vec3D tempv;
	double dist;

	double deriv_center_attract = 0;
	double maximum_dimension_2 = maximum_dimension*maximum_dimension;

	for (i = 0; i<3; i++)
	{
		temp1[i] = 0.0;
		tempv[i] = 0.0;
	}

	for(imol =0; imol < nmol ; imol++)
	{
		temp2 = center_arr[imol];
		double	dist_2 = temp2[0]*temp2[0] + temp2[1]*temp2[1] + temp2[2]*temp2[2];
		if (dist_2 > maximum_dimension_2)
			dist  = Vec3D::CalcDistance(&center_arr[imol],&temp1,ANGSTROM_U);
			deriv_center_attract = 4* dist* (maximum_dimension_2 - dist_2)* weight_cntr_attract;
			tempv[0] = temp2[0] * deriv_center_attract;
			tempv[1] = temp2[1] * deriv_center_attract;
			tempv[2] = temp2[2] * deriv_center_attract;
		force_cntl_array[imol].SetCoordFrom( tempv );
//	PrintLog("center_attract %2.2f dist %2.2f\n", center_attract, dist);
	}
	return TRUE;
}

void 
HaEmpiricalMod::CalculateCoarseGrainedBackbone()
{
	MolSet* pmset = GetMolSet();
	int nmol = pmset -> GetNMol();
	HaAtom* aptr;
	HaMolecule* pMol; 
	Vec3D ptr_new;
	double  numer_x;
	double 	numer_y;
	double 	numer_z;
	double 	com_x;
	double 	com_y;
	double 	com_z;
	int nres_tmp ;
	int nres_delta;
	int nres;
	static int nres_new =0;
	int i, j, imol;
	int icount =0;
	int jcount =0;
	
	if (nres_new ==  0)
	{
		for( imol =0; imol < nmol ; imol++)
		{
			pMol = pmset->GetMolByIdx(imol);
			nres = pMol -> GetNRes();
			nres_tmp = (int) floor(nres* 0.25);
			nres_delta = nres - nres_tmp*4;
			nres_new += nres_tmp + nres_delta;
		}
		atm_cgbb_array.resize(nres_new);
		rad_cgbb_array.resize(nres_new);
	}
	
	for( imol =0; imol < nmol ; imol++)
	{
		pMol = pmset->GetMolByIdx(imol);
		nres = pMol -> GetNRes();
		nres_tmp = (int) floor(nres* 0.25);
		nres_delta = nres - nres_tmp*4;

		for(i =0; i < nres_tmp; i++)
		{
			numer_x = 0;
			numer_y = 0;
			numer_z = 0;
			com_x = 0;
			com_y = 0;
			com_z = 0;
			for(j = 0; j < 4; j++)
			{
				aptr = (HaAtom*) atm_ca_array[jcount];
				numer_x += aptr -> GetX();
				numer_y += aptr -> GetY();
				numer_z += aptr -> GetZ();
				jcount++;
			}
			com_x = (numer_x*0.25);
			com_y = (numer_y*0.25);
			com_z = (numer_z*0.25);
			
			ptr_new.SetX(com_x);
			ptr_new.SetY(com_y);
			ptr_new.SetZ(com_z);
			
			atm_cgbb_array[icount] = ptr_new;
			rad_cgbb_array(icount+1) = 1.6 * ANG_TO_BOHR;
			icount++;
		}
		
		for(i =0; i < nres_delta; i++)
		{
			aptr = (HaAtom*) atm_ca_array[jcount];
			com_x = aptr -> GetX();
			com_y = aptr -> GetY();
			com_z = aptr -> GetZ();
			
			ptr_new.SetX(com_x);
			ptr_new.SetY(com_y);
			ptr_new.SetZ(com_z);
			
			atm_cgbb_array[icount] = ptr_new;
			rad_cgbb_array(icount+1) = aptr -> radius;
			jcount++;
			icount++;
//			PrintLog("nres_delta %d\n", nres_delta);
		}
	}
//	nres = pmset -> GetNRes();
//	for(i =0; i < nres; i++)
//	{
//		aptr = (HaAtom*) atm_ca_array[i];
//		PrintLog("SSSS %d  %2.3f %2.3f %2.3f\n",i, aptr->GetX(), aptr->GetY(),aptr->GetZ());
//		
//	}
}

double
HaEmpiricalMod::PenaltyRepulsion()
{
	HaAtom* aptr;
	HaAtom* aptr1;
	HaAtom* aptr2;
	HaAtom* aptrca1;
	HaAtom* aptrca2;
	HaResidue* res_ptr;
	HaResidue* res_ptr1;
	HaMolecule* pMol;
	HaMolecule* pMol1;
	int i = 0;
	int j = 0;
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();

	
	int imol;
	int jmol;
	double ene = 0;

	for( imol =0; imol < nmol ; imol++)
	{
		pMol = pmset->GetMolByIdx(imol);
		HaMolecule::ResidueIterator ritr(pMol);
		for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
		{	
			aptr = (HaAtom*) atm_sc_array[i];
			aptrca1 = (HaAtom*) atm_ca_array[i];
			
			i++;				

			for( jmol = imol; jmol < nmol ; jmol++)
			{
				j = 0;
				if (imol != jmol )
				{
					int number_shift = first_res_mol[jmol] ;
					pMol1 = pmset -> GetMolByIdx(jmol);
					HaMolecule::ResidueIterator ritr1(pMol1);
					for(res_ptr1 = ritr1.GetFirstRes(); res_ptr1; res_ptr1 = ritr1.GetNextRes())
					{	
						aptr1 = (HaAtom*) atm_sc_array[j+ number_shift ];
						aptrca2= (HaAtom*) atm_ca_array[j+ number_shift ];
						j++;
						
						ene += GetRepulEnergy(aptr,aptr1) ;
						ene += GetRepulEnergy(aptr,aptrca2) ;
						ene += GetRepulEnergy(aptrca1,aptrca2) ;
					}
				} 
			} //End of for( jmol)

		}
	}//End of for( imol)
    //PrintLog("ENERGY: %4.2f \n",ene);
	return ene;
}


int
HaEmpiricalMod::CalcRepulForceTorque(Vec3DValArray& force_array, Vec3DValArray& torque_array)
{

	HaAtom* aptr;
	HaAtom* aptr1;
	HaAtom* aptr2;
	HaAtom* aptrca1;
	HaAtom* aptrca2;
	HaResidue* res_ptr;
	HaResidue* res_ptr1;
	HaMolecule* pMol;
	HaMolecule* pMol1;
	int i = 0;
	int j = 0;
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	force_array.resize(nmol);
	torque_array.resize(nmol);
	Vec3D force_vec1;
	Vec3D force_vec2;
	Vec3D torque_vec1;
	Vec3D torque_vec2;
	Vec3D force_v;
	Vec3D torque_v;
	
	int imol;
	int jmol;
	double force;
	
	Vec3D trans;
	Vec3D trans1;
	Vec3D vec_com1;
	Vec3D vec_com2;
	Vec3D vec_com_ca1;
	Vec3D vec_com_ca2;

	Quaternion quat;
	Quaternion quat1;
	HaMat_double rmat;
	rmat.newsize(3,3);
	HaMat_double tranf_mat;
	tranf_mat.newsize(3,3);
	HaMat_double rot_mat;
	rot_mat.newsize(3,3);
	NumVector<double> torque(3);
	double cx, cy, cz;
	double cx1, cy1, cz1;

	Vec3DValArray com_array;
	com_array.resize(nmol);
	force_vec1[0]= 0;
	force_vec1[1]= 0;
	force_vec1[2]= 0;
	for( imol =0; imol < nmol ; imol++)
	{
		force_array[imol].SetCoordFrom(force_vec1);
		torque_array[imol].SetCoordFrom(force_vec1);
		pMol = pmset -> GetMolByIdx(imol);
		pMol->GetAverageCoord(cx, cy, cz);
		trans[0] = cx;
		trans[1] = cy;
		trans[2] = cz;
		com_array[imol].SetCoordFrom(trans);
	}

	for( imol =0; imol < nmol ; imol++)
	{
		force_vec1= force_array[imol];
		torque_vec1 = torque_array[imol];
		pMol = pmset->GetMolByIdx(imol);
		trans = com_array[imol];
		cx = trans[0];
		cy = trans[1];
		cz = trans[2];

//		PrintLog("imol %d    com: x=%4.2f y=%4.2f z=%4.2f\n", imol, trans[0]/(trans.length()), trans[1]/(trans.length()), trans[2]/(trans.length()) );
		HaMolecule::ResidueIterator ritr(pMol);
		for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
		{	
			aptr = (HaAtom*) atm_sc_array[i];
			aptrca1 = (HaAtom*) atm_ca_array[i];
			vec_com1[0] = aptr->GetX() - cx;
			vec_com1[1] = aptr->GetY() - cy;
			vec_com1[2] = aptr->GetZ() - cz;
			vec_com_ca1[0] = aptrca1->GetX() - cx;
			vec_com_ca1[1] = aptrca1->GetY() - cy;
			vec_com_ca1[2] = aptrca1->GetZ() - cz;

			
			i++;				

			for( jmol = imol+1; jmol < nmol ; jmol++)
			{
				j = 0;
				force_vec2= force_array[jmol];
				torque_vec2 = torque_array[jmol];
				trans = com_array[jmol];
				cx1 = trans[0];
				cy1 = trans[1];
				cz1 = trans[2];

				int number_shift = first_res_mol[jmol] ;
				pMol1 = pmset -> GetMolByIdx(jmol);
				
				HaMolecule::ResidueIterator ritr1(pMol1);
				for(res_ptr1 = ritr1.GetFirstRes(); res_ptr1; res_ptr1 = ritr1.GetNextRes())
				{	
					aptr1 = (HaAtom*) atm_sc_array[j+ number_shift ];
					aptrca2= (HaAtom*) atm_ca_array[j+ number_shift ];
					vec_com2[0] = aptr1->GetX() - cx1;
					vec_com2[1] = aptr1->GetY() - cy1;
					vec_com2[2] = aptr1->GetZ() - cz1;
					vec_com_ca2[0] = aptrca2->GetX() - cx1;
					vec_com_ca2[1] = aptrca2->GetY() - cy1;
					vec_com_ca2[2] = aptrca2->GetZ() - cz1;
					
					j++;
					
					force = GetRepulDerivative(aptr, aptr1) ;
					force_v.SetX( force* (aptr1->GetX() - aptr->GetX() ) );
					force_v.SetY( force* (aptr1->GetY() - aptr->GetY() ) );
					force_v.SetZ( force* (aptr1->GetZ() - aptr->GetZ() ) );
					force_vec1.SetX( force_vec1.GetX() + force_v.GetX() );
					force_vec1.SetY( force_vec1.GetY() + force_v.GetY() );
					force_vec1.SetZ( force_vec1.GetZ() + force_v.GetZ() );
					force_vec2.SetX( force_vec2.GetX() - force_v.GetX() );
					force_vec2.SetY( force_vec2.GetY()  - force_v.GetY());
					force_vec2.SetZ( force_vec2.GetZ()  - force_v.GetZ());
					
					Vec3D::VecProduct(torque_v, vec_com1, force_v);
					torque_vec1[0] += torque_v.GetX();
					torque_vec1[1] += torque_v.GetY();
					torque_vec1[2] += torque_v.GetZ();
					Vec3D::VecProduct(torque_v, vec_com2, force_v);
					torque_vec2[0] -= torque_v.GetX();
					torque_vec2[1] -= torque_v.GetY();
					torque_vec2[2] -= torque_v.GetZ();
					
					
					force = GetRepulDerivative(aptr,aptrca2) ;
					force_v.SetX( force* (aptrca2->GetX() - aptr->GetX() ) );
					force_v.SetY( force* (aptrca2->GetY() - aptr->GetY() ) );
					force_v.SetZ( force* (aptrca2->GetZ() - aptr->GetZ() ) );
					force_vec1.SetX( force_vec1.GetX() + force_v.GetX() );
					force_vec1.SetY( force_vec1.GetY() + force_v.GetY() );
					force_vec1.SetZ( force_vec1.GetZ() + force_v.GetZ() );
					force_vec2.SetX( force_vec2.GetX() - force_v.GetX() );
					force_vec2.SetY( force_vec2.GetY()  - force_v.GetY());
					force_vec2.SetZ( force_vec2.GetZ()  - force_v.GetZ());
					
					Vec3D::VecProduct(torque_v, vec_com1, force_v);
					torque_vec1[0] += torque_v.GetX();
					torque_vec1[1] += torque_v.GetY();
					torque_vec1[2] += torque_v.GetZ();
					Vec3D::VecProduct(torque_v, vec_com_ca2, force_v);
					torque_vec2[0] -= torque_v.GetX();
					torque_vec2[1] -= torque_v.GetY();
					torque_vec2[2] -= torque_v.GetZ();
					
					
					force = GetRepulDerivative(aptrca1,aptrca2) ;
					force_v.SetX( force* (aptrca2->GetX() - aptrca1->GetX() ) );
					force_v.SetY( force* (aptrca2->GetY() - aptrca1->GetY() ) );
					force_v.SetZ( force* (aptrca2->GetZ() - aptrca1->GetZ() ) );
					force_vec1.SetX( force_vec1.GetX() + force_v.GetX() );
					force_vec1.SetY( force_vec1.GetY() + force_v.GetY() );
					force_vec1.SetZ( force_vec1.GetZ() + force_v.GetZ() );
					force_vec2.SetX( force_vec2.GetX() - force_v.GetX() );
					force_vec2.SetY( force_vec2.GetY()  - force_v.GetY());
					force_vec2.SetZ( force_vec2.GetZ()  - force_v.GetZ());
					
					Vec3D::VecProduct(torque_v, vec_com_ca1, force_v);
					torque_vec1[0] += torque_v.GetX();
					torque_vec1[1] += torque_v.GetY();
					torque_vec1[2] += torque_v.GetZ();
					Vec3D::VecProduct(torque_v, vec_com_ca2, force_v);
					torque_vec2[0] -= torque_v.GetX();
					torque_vec2[1] -= torque_v.GetY();
					torque_vec2[2] -= torque_v.GetZ();
					 
				} 
				force_array[imol].SetCoordFrom(force_vec1);
				force_array[jmol].SetCoordFrom(force_vec2);
				torque_array[imol].SetCoordFrom(torque_vec1);
				torque_array[jmol].SetCoordFrom(torque_vec2);
			} //End of for( jmol)
		}

		pMol->GetStdMomInertRotMat(rot_mat);
		
		static HaVec_double vec(3);
		static Vec3D vec1;
		vec1 = torque_array[imol];
		
		torque = matmult(rot_mat, vec );
		vec1[0] = torque(1);
		vec1[1] = torque(2);
		vec1[2] = torque(3);
		torque_array[imol].SetCoordFrom( vec1 );
		Vec3D tempv= force_array[imol];
		tempv[0]= tempv[0];
		tempv[1]= tempv[1];
		tempv[2]= tempv[2];
		force_array[imol].SetCoordFrom( tempv );
//	PrintLog("imol %d    Force: x=%4.2f y=%4.2f z=%4.2f\n", imol, tempv[0]/(tempv.length()), tempv[1]/(tempv.length()), tempv[2]/(tempv.length()) );
//	PrintLog("imol %d    Force: x=%4.2f y=%4.2f z=%4.2f\n", imol, tempv[0], tempv[1], tempv[2] );

//	PrintLog("imol %d    Torque: x=%4.2f y=%4.2f z=%4.2f\n", imol, vec1[0], vec1[1], vec1[2] );
			

	}//End of for( imol)
		
	return TRUE;
}

double
HaEmpiricalMod::GetRepulDerivative(HaAtom* aptr1, HaAtom* aptr2)
{
	double sigma1 = aptr1 -> radius;
	double sigma2 = aptr2 -> radius;
	double force;
	double dx = aptr1 ->GetX() - aptr2 ->GetX() ;
	double dx_2 = dx*dx;
	double dy = aptr1 ->GetY() - aptr2 ->GetY() ;
	double dy_2 = dy*dy;
	double dz = aptr1 ->GetZ() - aptr2 ->GetZ() ;
	double dz_2 = dz*dz;
	double dist_2 = dx_2+ dy_2+ dz_2;
	double dist_6 = dist_2*dist_2*dist_2;
	double sigma = sigma1 + sigma2;
	double sigma_2 = sigma*sigma;

	if (dist_2 < sigma_2)
		force = (4.0* sigma_2* sigma_2)/ dist_6;
	else 
		force = 0.0f;

	return force;
}	

double HaEmpiricalMod::GetRepulEnergy(HaAtom* aptr1, HaAtom* aptr2)
{
	double ene;
	double sigma1 = (aptr1 -> radius);
	double sigma2 = (aptr2 -> radius);
	double dx = aptr1 ->GetX() - aptr2 ->GetX() ;
	double dx_2 = dx*dx;
	double dy = aptr1 ->GetY() - aptr2 ->GetY() ;
	double dy_2 = dy*dy;
	double dz = aptr1 ->GetZ() - aptr2 ->GetZ() ;
	double dz_2 = dz*dz;
	double dist_2 = dx_2+ dy_2+ dz_2;
	double dist_4 = dist_2*dist_2;
	double sigma = sigma1 + sigma2;
	double sigma_2 = sigma*sigma;

	//PrintLog("dist_4  %2.3f sigma1 %2.3f\n", dist_4, sigma1);
	if (dist_2 <= sigma_2)
		// if distance b/n 2 points less than sigma
		// ene = ( (sigma_2 *sigma_2) / dist_4) - 1.0; old version does not have a first derivative
		ene = ((sigma_2*sigma_2)/dist_4)-1.0+4*(sqrt(dist_2)/sigma-1.0);
	else 
		ene = 0.0f;
	return ene;
}
int HaEmpiricalMod::GetGeomCenterToy()  // added by jose
{
	// calculates the center of mass of each cylinder 
	//then the center of mass of the molecule
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	//PrintLog("Number of molecules: %d\n", nmol);
	HaAtom* aptr = NULL;
	int imol=0;
	center_arr_t.resize(nmol);
	double x_temp = 0;
	double y_temp = 0;
	double z_temp = 0;

	MoleculesType::iterator mol_itr;
	for( mol_itr=pmset->HostMolecules.begin(); mol_itr != pmset->HostMolecules.end(); mol_itr++)
	{
		aptr= (*mol_itr)->GetAtomByRef("GEN6.X");
		x_temp = aptr->GetX();
		y_temp = aptr->GetY();
		z_temp = aptr->GetZ();
		center_arr_t[imol].SetX(x_temp);
		center_arr_t[imol].SetY(y_temp);
		center_arr_t[imol].SetZ(z_temp);
		imol++;

	}
	center_arr_t.GetAverageCoord(x_temp,y_temp,z_temp);
	geom_center_t.SetX(x_temp);
	geom_center_t.SetY(y_temp);
	geom_center_t.SetZ(z_temp);

	//MolSet::ResidueIterator ritr(pmset);
	//for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
	//{
	//	AtomIteratorResidue aitr(res_ptr);
	//	for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
	//	{

	//	}
	//}
	return TRUE;
}
int HaEmpiricalMod::GetGeomCenter()  //added by jose 05-22-08
{
	// Find geometric center of all helices
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	int imol;
	if (com_flag)
	{
		CenterOfMass();
	}

	Vec3D* vec1;
	double temp_x = 0.0, temp_y = 0.0, temp_z = 0.0;

	for (imol = 0; imol < nmol; imol++)
	{
		vec1 = &center_arr[imol];
		temp_x += vec1->GetX();
		temp_y += vec1->GetY();
		temp_z += vec1->GetZ();
	}
	temp_x = double(temp_x / nmol);
	temp_y = double(temp_y / nmol);
	temp_z = double(temp_z / nmol);
	geom_center.SetX(temp_x);
	geom_center.SetY(temp_y);
	geom_center.SetZ(temp_z);
	return TRUE;
}


int HaEmpiricalMod::GetMaxDimension()
{
	//  units are in Angstroms
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	int imol;
	HaMolecule* pMol;
	double xmin, ymin, zmin;
	double xmax, ymax, zmax;
	Vec3D v1;
	Vec3D v2;
	maximum_dimension = -10e9;
	state_old.newsize(nmol); 
	state_new.newsize(nmol);
	exch12.newsize(nmol);
	for (imol=0;imol<nmol;imol++)
	{
		state_old.SetVal(imol+1,imol);
		exch12.SetVal(imol+1,0);
		state_new.SetVal(imol+1,imol);
	}
	for( imol =0; imol < nmol ; imol++)
	{
		pMol = pmset->GetMolByIdx(imol);
		pMol->GetMinMaxCrd(xmin, ymin, zmin, xmax, ymax, zmax);
		v1[0]= xmin;
		v1[1]= ymin;
		v1[2]= zmin;
		v2[0]= xmax;
		v2[1]= ymax;
		v2[2]= zmax;
		double length  = Vec3D::CalcDistance(&v1,&v2,ANGSTROM_U);
		if (length > maximum_dimension) maximum_dimension = length;
	}
//  maximum_dimension = nmol* 4.0 +8.0; //plane
//maximum_dimension = maximum_dimension +6.0;//
//	maximum_dimension = 7.5; // 7.5 used for toy model 
	maximum_dimension =  35.0; // r0 = 27.0 bohr, minimum radio is 11A.
	// 6.35A; //9.66A;// 18.25 bohr; used for 4 helices
//	PrintLog("maximum_dimension %2.1f \n", maximum_dimension);
// 8.0 A is diameter of the helix
	return TRUE;
}
int HaEmpiricalMod::QuantSampling() // added by jose
{	// This function evaluates the rate of exchange of positions of rigid helices on the plane XY
	// the array state_new gives you the current positions of the rigid helices with respect to the first molecule
	// the angle most negative corresponds to the position 1
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	int imol=0;
	Vec3D vec1;
	Vec3D vec2;
	Vec3D temp;
	double xnew, ynew, znew, imol_angle, ang, temp3;
	int temp2;
	Vec3D point0;
	Vec3D point1;
	Vec3D point2;
	HaMat_double rot_mat;
	HaVec_double angle_vec, angle_vec1;
	HaVec_int ind_array; // Index array used by the function QuantSampling
	HaVec_double angle_array; // angle array used by the function QuantSampling
	ind_array.resize(nmol);
	angle_array.resize(nmol);
	angle_vec.newsize(nmol); angle_vec1.newsize(nmol); 
	rot_mat.resize(3,3);
	rot_mat.SetVal(1,3,0);
	rot_mat.SetVal(2,3,0);
	rot_mat.SetVal(3,3,1);
	rot_mat.SetVal(3,1,0);
	rot_mat.SetVal(3,2,0);
	int ind_save=0;
	Vec3DValArray center_arr_temp;
	char buf[256];
	std::fstream state_file;
	state_file.open("state_file.dat", std::ios::out | std::ios::app);
	center_arr_temp.resize(nmol);
	/*this->GetGeomCenter();
	*/
	/*if (com_flag)
	{
		CenterOfMassToy();
	}*/  // this for molecule with residues templates 
	if (com_flag)	this->GetGeomCenterToy(); // calculates center_arr_t and geom_center_t
	//translate to nmol=0

	vec1.SetX(center_arr_t[0].GetX());
	vec1.SetY(center_arr_t[0].GetY());
	vec1.SetZ(center_arr_t[0].GetZ());
	for (imol=0;imol<nmol;imol++)
	{
		temp.SetX(center_arr_t[imol].GetX()-vec1.GetX());
		temp.SetY(center_arr_t[imol].GetY()-vec1.GetY());
		temp.SetZ(center_arr_t[imol].GetZ()-vec1.GetZ());
		center_arr_temp[imol] = temp;
		//PrintLog("imol %d, center_arr x= %f y= %f z= %f\n",imol,center_arr[imol].GetX(),center_arr[imol].GetY(),center_arr[imol].GetZ());
	}
	//PrintLog("geom_center x= %g y= %g z= %g\n",geom_center.GetX(),geom_center.GetY(),geom_center.GetZ());
	vec2.SetX(geom_center_t.GetX()-vec1.GetX());
	vec2.SetY(geom_center_t.GetY()-vec1.GetY());
	vec2.SetZ(0.0);
	point1.SetX(1.0); point1.SetY(0.0); point1.SetZ(0.0);
	point0.SetX(0.0); point0.SetY(0.0); point0.SetZ(0.0);
	point2=vec2;
	ang = Vec3D::CalcAngle(&point1, &point0, &point2);
	if ((geom_center_t.GetY()-vec1.GetY())>0.0)
		ang=2*PI-ang;
	rot_mat.SetVal(1,1,cos(ang));
	rot_mat.SetVal(1,2,-sin(ang));
	rot_mat.SetVal(2,1,sin(ang));
	rot_mat.SetVal(2,2,cos(ang));
	for (imol=0;imol<nmol;imol++)
	{
		temp.SetX(center_arr_temp[imol].GetX());
		temp.SetY(center_arr_temp[imol].GetY());
		temp.SetZ(center_arr_temp[imol].GetZ());
		xnew= temp.GetX()*rot_mat.GetVal(1,1) + temp.GetY()*rot_mat.GetVal(1,2);
		ynew= temp.GetX()*rot_mat.GetVal(2,1) + temp.GetY()*rot_mat.GetVal(2,2);
		znew= temp.GetZ()*rot_mat.GetVal(3,1) + temp.GetY()*rot_mat.GetVal(3,2);
		point2.SetX(xnew); point2.SetY(ynew); point2.SetZ(znew);
		ang = Vec3D::CalcAngle(&point1,&point0,&point2);
		if (ynew>0.0)
			ang= -ang;
		//PrintLog("ang %f\n", ang);
		angle_vec.SetVal(imol+1,ang);
		angle_vec1.SetVal(imol+1,ang);
	}

	double max_val=1000000;

	for (int ind=0;ind<nmol;ind++)
	{
		double min_angle = 100000;
		for (imol=0;imol<nmol;imol++)
		{
			if (imol>0)
			{
				imol_angle = angle_vec1.GetVal(imol+1);
				if (imol_angle < min_angle)
				{
					min_angle = imol_angle;
					ind_save = imol+1;
					//PrintLog("ind_save %d\n",ind_save);
				}
			}
		}
		angle_vec1.SetVal(ind_save, max_val);
		ind_array.SetVal(ind+1,ind_save);
		angle_array.SetVal(ind+1,min_angle);
		//PrintLog("angle_array %g, ind_array %d\n", angle_array.GetVal(ind+1),ind_array.GetVal(ind+1));
	}
	
	for (imol=0;imol<nmol;imol++)
	{
		if (imol>0)  // avoid the first molecule 
		{
			temp3= angle_vec.GetVal(imol+1);
			for (int ind1=0;ind1<nmol-1;ind1++)
			{
				if (temp3 == angle_array.GetVal(ind1+1))
				{
					state_new.SetVal(imol+1,ind1+1);
					//PrintLog("imol %d, state_new %d\n",imol, state_new.GetVal(imol+1));
				}
			}
		}

	}
	for (imol=0;imol<nmol;imol++)
	{
		//PrintLog("imol %d, state_old %d\n",imol,state_old.GetVal(imol+1));
		if (state_old.GetVal(imol+1)!=state_new.GetVal(imol+1))
		{
			temp2=exch12.GetVal(imol+1)+1;
			exch12.SetVal(imol+1,temp2);
		}
	}
	for (imol=0;imol<nmol;imol++)
	{
		//print >>state_out,state_old.GetVal(imol+1);
		temp2= state_new.GetVal(imol+1);
		state_old.SetVal(imol+1,temp2);

	}
	for (imol=0;imol<nmol;imol++)
	{
	//	if (imol>0)
	//		PrintLog("imol %d exch %d \n", imol, exch12.GetVal(imol+1));
	//}
		sprintf(buf,"%2d  ", state_new.GetVal(imol+1));
		state_file<< buf;
	}
	state_file<< std::endl;
	return TRUE;
}

double HaEmpiricalMod::LJ_eneCylinder() //added by jose
{
	// Calculates Lennard Jones interactions between sticky points of cylinders
	// there are two sticky points located at the extremes of the cylinders
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	double e=4.0;   // well parameter
	double r0=sqrt(89.0);  // 3 (angstroms) distance at the minimum
	HaAtom* aptr = NULL;
	HaAtom* aptr1 = NULL;
	HaAtom* aptr2 = NULL;
	HaAtom* aptr3 = NULL;
	HaAtom* aptr4 = NULL;
	HaAtom* aptr5 = NULL;
	HaAtom* aptr6 = NULL;
	HaAtom* aptr7 = NULL;
	double dist_2 = 0.0;
	double dist = 0.0;
	double lambda = 0.8;
	double add = lambda*100;
	int i=0;
	double A=1.4;

	this->QuantSampling();
	//fstream energy_LJ_file;
	//energy_LJ_file.open("energy_LJstat_file.dat", ios::out | ios::app);
	double ene = 0.0;

	// we start a first iteration loop through all molecules in molset
	MoleculesType::iterator mol_itr;
	for( mol_itr=pmset->HostMolecules.begin(); mol_itr != pmset->HostMolecules.end(); mol_itr++)
	{
		i++;
		aptr = (*mol_itr)->GetAtomByRef("GEN1.XA1");
		aptr1 = (*mol_itr)->GetAtomByRef("GEN1.XA2");
		aptr2 = (*mol_itr)->GetAtomByRef("GEN11.XA1");
		aptr3 = (*mol_itr)->GetAtomByRef("GEN11.XA2");

		// we start a second iteration loop through mol_itr to all molecules in molset
		MoleculesType::iterator mol_itr1;
		for (mol_itr1 = pmset->HostMolecules.begin()+i; mol_itr1 != pmset->HostMolecules.end(); mol_itr1++)
		{
			aptr4 = (*mol_itr1)->GetAtomByRef("GEN1.XA1");
			aptr5 = (*mol_itr1)->GetAtomByRef("GEN1.XA2");
			aptr6 = (*mol_itr1)->GetAtomByRef("GEN11.XA1");
			aptr7 = (*mol_itr1)->GetAtomByRef("GEN11.XA2");

			//dist_2 = Vec3D::CalcDistanceSq(aptr, aptr4,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));
			//
			//dist_2 = Vec3D::CalcDistanceSq(aptr, aptr5,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));

			//dist_2 = Vec3D::CalcDistanceSq(aptr, aptr6,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));
			//
			//dist_2 = Vec3D::CalcDistanceSq(aptr, aptr7,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));
			//
			//dist_2 = Vec3D::CalcDistanceSq(aptr1, aptr4,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));
			//
			//dist_2 = Vec3D::CalcDistanceSq(aptr1, aptr5,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));

			//dist_2 = Vec3D::CalcDistanceSq(aptr1, aptr6,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));
			//
			//dist_2 = Vec3D::CalcDistanceSq(aptr1, aptr7,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));

			//dist_2 = Vec3D::CalcDistanceSq(aptr2, aptr4,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));
			//
			//dist_2 = Vec3D::CalcDistanceSq(aptr2, aptr5,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));

			//dist_2 = Vec3D::CalcDistanceSq(aptr2, aptr6,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));
			//
			//dist_2 = Vec3D::CalcDistanceSq(aptr2, aptr7,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));

			//dist_2 = Vec3D::CalcDistanceSq(aptr3, aptr4,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));
			//
			//dist_2 = Vec3D::CalcDistanceSq(aptr3, aptr5,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));

			//dist_2 = Vec3D::CalcDistanceSq(aptr3, aptr6,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));
			//
			//dist_2 = Vec3D::CalcDistanceSq(aptr3, aptr7,ANGSTROM_U)+add;
			//ene += e*(pow((r0*r0/dist_2),6)-2*pow((r0*r0/dist_2),3));

			// ++++++++++++++++++++++++++++++++++++++++++++++++++++++
			dist = Vec3D::CalcDistanceSq(aptr, aptr4,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));
			
			dist = Vec3D::CalcDistanceSq(aptr, aptr5,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));

			dist = Vec3D::CalcDistanceSq(aptr, aptr6,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));
			
			dist = Vec3D::CalcDistanceSq(aptr, aptr7,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));
			
			dist = Vec3D::CalcDistanceSq(aptr1, aptr4,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));
			
			dist = Vec3D::CalcDistanceSq(aptr1, aptr5,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));

			dist = Vec3D::CalcDistanceSq(aptr1, aptr6,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));
			
			dist = Vec3D::CalcDistanceSq(aptr1, aptr7,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));

			dist = Vec3D::CalcDistanceSq(aptr2, aptr4,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));
			
			dist = Vec3D::CalcDistanceSq(aptr2, aptr5,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));

			dist = Vec3D::CalcDistanceSq(aptr2, aptr6,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));
			
			dist = Vec3D::CalcDistanceSq(aptr2, aptr7,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));

			dist = Vec3D::CalcDistanceSq(aptr3, aptr4,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));
			
			dist = Vec3D::CalcDistanceSq(aptr3, aptr5,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));

			dist = Vec3D::CalcDistanceSq(aptr3, aptr6,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));
			
			dist = Vec3D::CalcDistanceSq(aptr3, aptr7,ANGSTROM_U);
			ene += 4*e*(pow((r0*r0/dist_2),6)-r0/sqrt(dist_2));
			//---------------------------------------------------

		}
	}
	//energy_LJ_file<< ene<<endl;
	//energy_LJ_file << ene*(1-lambda)<<endl;
	//PrintLog("LJ energy: %f\n",ene);
	return /*ene;*/ ene*(1-lambda);


}
double HaEmpiricalMod::LJ_ene() //added by jose
{ // This function evaluates the Lennard Jones potential energy of a bundle of N helical rigid bodies
  // with respect to the first molecule. There will at least (N-1)! basins of energy.
  // The interacting residues are located at the extremes of the helices
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	HaMolecule* pMol;
	HaMolecule* pMol1;
	HaResidue* res_ptr;
	HaResidue* res_ptr1;

	double e=0.50;   // well parameter

	double r0=27.0;  // distance at the minimum, 14.0 for two gly helices
	double dist_sc_2;
	HaAtom* aptr;
	HaAtom* aptr1;
	HaAtom* aptrca;
	HaAtom* aptrca1; 
	int nMol, nMol1;
	int imol, jmol;
	int i=0,j=0,h1,h2;
	int s, s1;
	this->QuantSampling();
	std::fstream energy_LJ_file;
	energy_LJ_file.open("energy_LJ_file.dat", std::ios::out | std::ios::app);
	double ene = 0.0;
	//double energy_total = 0.0, dist;
    for (imol= 0; imol < nmol; imol++)
	{
		pMol = pmset -> GetMolByIdx(imol);
		nMol = pMol ->GetNRes();
		HaMolecule::ResidueIterator ritr(pMol);
		s = state_new(imol+1);
		h1 = 0;
		for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
		{
			aptr = (HaAtom*) atm_sc_array[i]; //sc
			//PrintLog("rsc= %g\n",aptr->radius);
			//aptrca = (HaAtom*) atm_ca_array[i]; //ca
			//PrintLog("rca= %g\n",aptrca->radius);
			if (h1==0 || h1==(nMol-1))
			{
				for (jmol= imol+1; jmol < nmol;jmol++)
				{
					int number_shift = first_res_mol[jmol];
					pMol1 = pmset -> GetMolByIdx(jmol);
					nMol1 = pMol1->GetNRes();
					s1 = state_new(jmol+1);
					HaMolecule::ResidueIterator ritr1(pMol1);
					j = 0;
					h2 = 0;
					if (abs(s-s1)==1 || abs(s-s1)==nmol-1)
					{
						for(res_ptr1 = ritr1.GetFirstRes(); res_ptr1; res_ptr1 = ritr1.GetNextRes())
						{
							aptr1 = (HaAtom*) atm_sc_array[j+ number_shift ];
							//aptrca1= (HaAtom*) atm_ca_array[j+ number_shift ];
							//if (((h1==0 & h2==0 ) | (h1==(nMol-1) & h2==(nMol1-1))) &(abs(s-s1)==1 | abs(s-s1)==nmol-1))
							if ((h1==0 && h2==0 ) || (h1==(nMol-1) && h2==(nMol1-1)))
							//if (h2==0 | h2==(nMol1-1))
							{
								dist_sc_2 = Vec3D::CalcDistanceSq(aptr, aptr1);
								ene += e*(pow((r0*r0/dist_sc_2),6)-2*pow((r0*r0/dist_sc_2),3));
								//PrintLog("h1= %d h2= %d\n", h1,h2);
							}
							h2++;
							j++;
						}
					//double dist_ijmol = Vec3D::CalcDistanceSq(&center_arr[imol], &center_arr[jmol]);
					//energy_LJ_file<<imol<<" "<<jmol<<" "<<sqrt(dist_ijmol)<<"  "<< ene<<endl;
					}
				}
			}
			h1++;
			i++;
		} 
	}
	energy_LJ_file<<sqrt(dist_sc_2)<<"  "<< ene<< std::endl;
	return ene;
}

double HaEmpiricalMod::BuriedEnergyCyl() // added by jose for cylinders 
{
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	HaAtom* aptr;
	HaAtom* aptr1;
	double weight = 1.0;
	int i=0;
	double ene= 0.0;
	double zup= 15, zdown= -15, z1= 0, z2= 0;

	// we start a iteration through all molecules in molset
	MoleculesType::iterator mol_itr;
	for( mol_itr=pmset->HostMolecules.begin(); mol_itr != pmset->HostMolecules.end(); mol_itr++)
	{
		aptr = (*mol_itr)->GetAtomByRef("GEN2.X");
		aptr1 = (*mol_itr)->GetAtomByRef("GEN10.X");
		z1 = aptr->GetZ();
		z2 = aptr1->GetZ();
		if (z1<zup)		ene += weight*(z1-zup)*(z1-zup);
		if (z2>zdown)	ene += weight*(z2-zdown)*(z2-zdown);
	}
//	PrintLog("Buried energy: %f\n", ene);
	return ene;
}
double HaEmpiricalMod::BuriedEnergy()  //added by jose
{
	// It gives a harmonic energy potential if the resiudes located at the exterior of the membrane are buried
	// the width of the membrane is 30 A and the helices length is 50 A
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	HaAtom* aptr;
	HaAtom* aptr1;
	double weight = 1.0;
	int i=0;
	double ene= 0.0;
	double z = 0.0 , zup= 26.46, zdown = -24.57;
	char name1[8][19] = {"$HELIX_1$GLY9.X","$HELIX_2$GLY9.X","$HELIX_3$GLY9.X","$HELIX_4$GLY9.X","$HELIX_1$ALA32.X","$HELIX_2$ALA32.X","$HELIX_3$ALA32.X","$HELIX_4$ALA32.X"};
	//char name1[8][19] = {"GLY8.X","GLY46.X","GLY84.X","GLY122.X","ALA31.X","ALA69.X","ALA107.X","ALA145.X"};
	for (i=0; i <8;i++)
	{   
		aptr = pmset->GetAtomByRef(name1[i]);
		z= aptr->GetZ();
		//double x = aptr->GetX();
		//double y = aptr->GetY();
		//PrintLog("nmol= %d x= %f y= %f z= %f\n",nmol, x,y,z);
		if (i<nmol)
		{
			if (z < zup)
			{
				ene += (z-zup)*(z-zup)*weight;
				//PrintLog("z= %f zup= %f ene %f\n",z, zup, ene);
			}
		}
		else
		{
			if (z>zdown)
			{	ene += (z-zdown)*(z-zdown)*weight;
				//PrintLog("ene %f\n",ene);
			}
		}
	}
	return ene;



}
double HaEmpiricalMod::LJState_ene() //added by jose
{ // This function evaluates the Lennard Jones potential energy of a bundle of N helical rigid bodies
  // with respect to the first molecule. There will at least (N-1)! basins of energy.
  // The interacting residues are located at the extremes of the helices
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	HaMolecule* pMol;
	HaMolecule* pMol1;
	HaResidue* res_ptr;
	HaResidue* res_ptr1;
	double e=2.0;   // well parameter
	double r0=27.0;  // distance at the minimum, 14.0 for two gly helices
	double dist_sc_2;
	HaAtom* aptr;
	HaAtom* aptr1;
	HaAtom* aptrca;
	HaAtom* aptrca1; 
	int nMol, nMol1;
	int imol, jmol;
	int h1,h2;
	int s, s1;
	this->QuantSampling();
	std::fstream energy_LJ_file;
	energy_LJ_file.open("energy_LJstat_file.dat", std::ios::out | std::ios::app);
	double ene = 0.0;
	//int state[4]={0,1,2,3};
	int	state[6][4]={0,1,2,3,
					 0,1,3,2,
					 0,2,1,3,
					 0,3,1,2,
					 0,3,2,1,
					 0,2,3,1};
	double energy_total = 0.0, dist;
    for (imol= 0; imol < nmol-1; imol++)
	{
		pMol = pmset -> GetMolByIdx(state[curr_state-1][imol]);
		nMol = pMol ->GetNRes();
		HaMolecule::ResidueIterator ritr(pMol);
		h1 = 0;
		int number_shift1 = first_res_mol[state[curr_state-1][imol]];
		for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
		{
			aptr = (HaAtom*) atm_sc_array[h1+number_shift1]; //sc
			//aptrca = (HaAtom*) atm_ca_array[i]; //ca
			if (h1==0 || h1==(nMol-1))
			{
					jmol= state[curr_state-1][imol+1];
					int number_shift = first_res_mol[jmol];
					pMol1 = pmset -> GetMolByIdx(jmol);
					nMol1 = pMol1->GetNRes();
					HaMolecule::ResidueIterator ritr1(pMol1);
					h2 = 0;
					for(res_ptr1 = ritr1.GetFirstRes(); res_ptr1; res_ptr1 = ritr1.GetNextRes())
					{
						aptr1 = (HaAtom*) atm_sc_array[h2+ number_shift ];
						//aptrca1= (HaAtom*) atm_ca_array[j+ number_shift ];
						//if (((h1==0 & h2==0 ) | (h1==(nMol-1) & h2==(nMol1-1))) &(abs(s-s1)==1 | abs(s-s1)==nmol-1))
						if ((h1==0 && h2==0 ) || (h1==(nMol-1) && h2==(nMol1-1)))
						//if (h2==0 | h2==(nMol1-1))
						{
							dist_sc_2 = Vec3D::CalcDistanceSq(aptr, aptr1);
							ene += e*(pow((r0*r0/dist_sc_2),6)-2*pow((r0*r0/dist_sc_2),3));
						}
						h2++;
					}
					//double dist_ijmol = Vec3D::CalcDistanceSq(&center_arr[imol], &center_arr[jmol]);
					//energy_LJ_file<<imol<<" "<<jmol<<" "<<sqrt(dist_ijmol)<<"  "<< ene<<endl;
			}
			h1++;
		} 
	}
	pMol = pmset -> GetMolByIdx(state[curr_state-1][nmol-1]);
	nMol = pMol ->GetNRes();
	HaMolecule::ResidueIterator ritr(pMol);
	h1 = 0;
	int number_shift1 = first_res_mol[state[curr_state-1][nmol-1]];
	for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
	{
		aptr = (HaAtom*) atm_sc_array[h1+number_shift1]; //sc
		//aptrca = (HaAtom*) atm_ca_array[i]; //ca
		if (h1==0 || h1==(nMol-1))
		{
				jmol= 0;
				int number_shift = first_res_mol[jmol];
				pMol1 = pmset -> GetMolByIdx(jmol);
				nMol1 = pMol1->GetNRes();
				HaMolecule::ResidueIterator ritr1(pMol1);
				h2 = 0;
				for(res_ptr1 = ritr1.GetFirstRes(); res_ptr1; res_ptr1 = ritr1.GetNextRes())
				{
					aptr1 = (HaAtom*) atm_sc_array[h2+ number_shift ];
					//aptrca1= (HaAtom*) atm_ca_array[j+ number_shift ];
					//if (((h1==0 & h2==0 ) | (h1==(nMol-1) & h2==(nMol1-1))) &(abs(s-s1)==1 | abs(s-s1)==nmol-1))
					if ((h1==0 && h2==0 ) || (h1==(nMol-1) && h2==(nMol1-1)))
					//if (h2==0 | h2==(nMol1-1))
					{
						dist_sc_2 = Vec3D::CalcDistanceSq(aptr, aptr1);
						ene += e*(pow((r0*r0/dist_sc_2),6)-2*pow((r0*r0/dist_sc_2),3));
					}
					h2++;
				}
				//double dist_ijmol = Vec3D::CalcDistanceSq(&center_arr[imol], &center_arr[jmol]);
				//energy_LJ_file<<imol<<" "<<jmol<<" "<<sqrt(dist_ijmol)<<"  "<< ene<<endl;
		}
		h1++;
	} 
	energy_LJ_file<< ene<< std::endl;
	return ene;
}
int HaEmpiricalMod::InitCylinders()  //added by jose
{
	// function the initialize the radius values of dummy atoms on Cylinders
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	HaAtom* aptr=NULL;
	double radius_vw = 4.0; 

	MoleculesType::iterator mol_itr;
	for( mol_itr=pmset->HostMolecules.begin(); mol_itr != pmset->HostMolecules.end(); mol_itr++)
	{
		AtomIteratorMolecule aitr(*mol_itr);
		for(aptr= aitr.GetFirstAtom(); aptr;aptr= aitr.GetNextAtom())
		{
			if (!stricmp_loc(aptr->GetName(),"X")) aptr->radius = radius_vw;
			else
				aptr->radius = 0.5;
//			aptr->image_radius = radius_vw;
		}
	}
    SetMolHost(pmset);
	return TRUE;
}
double
HaEmpiricalMod::PenaltyRepulsionCyl() //added by jose
{
	MolSet* pmset = GetMolSet();
	HaAtom* aptr = NULL;
	HaAtom* aptr1 = NULL;
	double ene = 0.0;
	int i=0;

	// first loop iteration
	MoleculesType::iterator mol_itr;
	for( mol_itr=pmset->HostMolecules.begin(); mol_itr != pmset->HostMolecules.end(); mol_itr++)
	{
		i++;
		AtomIteratorMolecule aitr(*mol_itr);
		for (aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
		{
			if (stricmp_loc(aptr->GetName(),"X")) continue;
			// second loop iteration
			MoleculesType::iterator mol_itr1;
			for( mol_itr1=pmset->HostMolecules.begin()+i; mol_itr1 != pmset->HostMolecules.end(); mol_itr1++)\
			{
				AtomIteratorMolecule aitr1(*mol_itr1);
				for (aptr1 = aitr1.GetFirstAtom(); aptr1; aptr1 = aitr1.GetNextAtom())
				{
					if (stricmp_loc(aptr1->GetName(),"X")) continue;
					ene += GetRepulEnergy(aptr,aptr1);
				}
			}
		}
	}
//	PrintLog("Repulsion Energy: %f\n", ene);
	return ene;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// HaMolMembraneMod Class
// knowledge based potential, future coarse grained models structure for membrane protein interactions 
// ScoreEnergy() function used for Intermolecular interactions, jose November 2008
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////

HaMolMembraneMod::HaMolMembraneMod(MolSet* new_phost_mset):
HaCompMod(COMP_MOD_MEMBRANE,new_phost_mset)
{
	MolMechModule = new HaMolMechMod(new_phost_mset);
	SetStdParams();
}

//HaMolMembraneMod::HaMolMembraneMod(MolSet* new_phost_mset)//: MolMechModule(new HaMolMechMod(new_phost_mset))
//{
//	MolMechModule = new HaMolMechMod(new_phost_mset); // prefer composition over inheritance 
//	MolMechModule = new_phost_mset->GetMolMechMod(true);
//	//MolMechModule = HaMolMechMod(new_phost_mset);
//	SetStdParams();
//}

HaMolMembraneMod::~HaMolMembraneMod()
{
	delete MolMechModule;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Set Standard parameters of HaMolMembraneMod
//
///////////////////////////////////////////////////////////////////////////////////////////////

int HaMolMembraneMod::SetStdParams()
{
	module_to_init_flag = TRUE;
	module_to_init_HaMolMechMod = FALSE;
	build_nb_coarsegrained_contact_list = TRUE;
	pairwiseDfire_flag = true;

	pairwiseDfire_core = TRUE;
	pairwiseDfire_sa = FALSE;

	display_results_flag = FALSE;

	MolMechModule->SetStdParams();
	MolMechModule->p_mm_model->calc_electr_flag = FALSE; // Turn off the flag of Electrostatic calculations

	nonbond_DFIRE_cutoff_dist = 15.0;

	MAX_radius = 15.0; // angstroms

	vdw_HardSphere = false;
	vdw_x4_flag = true;
	vdw_x12_flag = false;
	vdw_x4_f_flag = false;

	angle_stat = 40.1;
	stdev = 5.4;

	return TRUE;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Initialization of HaMolMembraneMod (might be modified)
//		AtomsCentroids
//		Residues
//		CentreAtoms
//		LipidInterfaceAtoms
//
///////////////////////////////////////////////////////////////////////////////////////////////

int
HaMolMembraneMod::Initialize()
{
	ClearMembraneModel();

	if ( module_to_init_HaMolMechMod ) 
	{
		MolMechModule->Initialize();
		module_to_init_HaMolMechMod = FALSE;
	}

	HaAtom* aptr;

	AtomIteratorMolSet aitr(phost_mset); //phost_mset
	for(aptr = aitr.GetFirstAtom(); aptr ; aptr = aitr.GetNextAtom())
	{
		std::string at_name = aptr->GetName();

		if ((at_name == "X") || (at_name == "CA"))
				AtomsCentroids.push_back(aptr);
	}

	HaResidue* rptr;  // jose October 22, 2008

	ResidueIteratorMolSet res_ptr(phost_mset); // jose October 22, 2008
	for (rptr = res_ptr.GetFirstRes(); rptr ; rptr = res_ptr.GetNextRes()) // jose October 22, 2008
	{
		Residues.push_back(rptr); // jose October 22, 2008

		/*if (rptr->GetSerNo()==16)   //2zta
		{
			aptr = rptr->GetAtomByName("CA");
		    CentreAtoms.push_back(aptr);
		}*/

		if (rptr->GetSerNo()==84)   //glycophorin
		{
			aptr = rptr->GetAtomByName("CA");
		    CentreAtoms.push_back(aptr);
		}
		if (rptr->GetSerNo()==95)
		{
			aptr = rptr->GetAtomByName("CA");
		    LipidInterfaceAtoms.push_back(aptr);
		}
		if (rptr->GetSerNo()==73)
		{
			aptr = rptr->GetAtomByName("CA");
		    LipidInterfaceAtoms.push_back(aptr);
		}

		/*if (rptr->GetSerNo()==13)   //zeta-zeta
		{
			aptr = rptr->GetAtomByName("CA");
		    CentreAtoms.push_back(aptr);
		}
		if (rptr->GetSerNo()==23)
		{
			aptr = rptr->GetAtomByName("CA");
		    LipidInterfaceAtoms.push_back(aptr);
		}
		if (rptr->GetSerNo()==4)
		{
			aptr = rptr->GetAtomByName("CA");
		    LipidInterfaceAtoms.push_back(aptr);
		}*/

		/*if (rptr->GetSerNo()==65 || rptr->GetSerNo() == 165)   //Erb2 51-77
		{
			aptr = rptr->GetAtomByName("CA");
		    CentreAtoms.push_back(aptr);
		}
		if (rptr->GetSerNo()==175 || rptr->GetSerNo() == 75)
		{
			aptr = rptr->GetAtomByName("CA");
		    LipidInterfaceAtoms.push_back(aptr);
		}
		if (rptr->GetSerNo()==153 || rptr->GetSerNo() == 53)
		{
			aptr = rptr->GetAtomByName("CA");
		    LipidInterfaceAtoms.push_back(aptr);
		}*/

		/*if (rptr->GetSerNo()==65 || rptr->GetSerNo() == 165)   //Erb2 52-73
		{
			aptr = rptr->GetAtomByName("CA");
		    CentreAtoms.push_back(aptr);
		}
		if (rptr->GetSerNo()==173 || rptr->GetSerNo() == 73)
		{
			aptr = rptr->GetAtomByName("CA");
		    LipidInterfaceAtoms.push_back(aptr);
		}
		if (rptr->GetSerNo()==154 || rptr->GetSerNo() == 54)
		{
			aptr = rptr->GetAtomByName("CA");
		    LipidInterfaceAtoms.push_back(aptr);
		}*/

		/*if (rptr->GetSerNo()==555)   //EphA1 PH=6.3 544-569
		{
			aptr = rptr->GetAtomByName("CA");
		    CentreAtoms.push_back(aptr);
		}
		if (rptr->GetSerNo()==544 )
		{
			aptr = rptr->GetAtomByName("CA");
		    LipidInterfaceAtoms.push_back(aptr);
		}
		if (rptr->GetSerNo()==568 )
		{
			aptr = rptr->GetAtomByName("CA");
		    LipidInterfaceAtoms.push_back(aptr);
		}*/

		/*if (rptr->GetSerNo()==555)   //EphA1 PH=4.3 544-569
		{
			aptr = rptr->GetAtomByName("CA");
		    CentreAtoms.push_back(aptr);
		}
		if (rptr->GetSerNo()==545 )
		{
			aptr = rptr->GetAtomByName("CA");
		    LipidInterfaceAtoms.push_back(aptr);
		}
		if (rptr->GetSerNo()==569 )
		{
			aptr = rptr->GetAtomByName("CA");
		    LipidInterfaceAtoms.push_back(aptr);
		}*/

	}

	SetCoarseGrainedDFireCoreParams();
	module_to_init_flag = FALSE;
	this->CbettaSetUp();
	this->LineSegments();

	return TRUE;
}

int 
HaMolMembraneMod::ClearMembraneModel()
{
	AtomsCentroids.clear();
	Residues.clear();
	LipidInterfaceAtoms.clear();
	CentreAtoms.clear();
	
	//module_to_init_flag = TRUE;

	return TRUE;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Set DFIRE parameters
//  Force centroids atoms:
//     X  side chain force centroid
//	   CA alpha carbon force centroid
//
///////////////////////////////////////////////////////////////////////////////////////////////

int HaMolMembraneMod::SetCoarseGrainedDFireCoreParams() // jose October 22, 2008
{
	//if( module_to_init_flag ) 
	//	Initialize();

	LoadDFireCoreParams(); 	// Load DFIRE_SCM Core and Adamian 2006 empirical lipid potential parameters

	HaResidue* res_ptr;
	HaAtom* aptr;

	MolSet::ResidueIterator ritr(phost_mset); //phost_mset
	for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
	{
		std::string res_name = res_ptr->GetName();
		if (res_name != "GLY")	
			aptr = res_ptr->GetAtomByName("X"); 		// set Side Chain vdw radius parameter
			//HaAtom* aptr1 = res_ptr->GetAtomByName("CA");
			//aptr1->vdw_rad = 1.9;}
			//PrintLog("radius CA %s, %4.2f\n", res_name.c_str(), aptr1->vdw_rad);}
		else 
			aptr = res_ptr->GetAtomByName("CA");

		aptr->SetVdWRad( sc_vdwradius.GetVal(res_name.c_str()));
		//PrintLog("radius %s , %4.2f\n", res_name.c_str(), sc_vdwradius.GetVal(res_name.c_str()));
	}

	BuildNonBondSCContactList();
	BuildClashAtomList();
	return TRUE;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Build non-bonded coarse grained atom list for rigid body calculations 
//  atoms are from different rigid body
//
///////////////////////////////////////////////////////////////////////////////////////////////

bool HaMolMembraneMod::BuildNonBondSCContactList() //@ jose October 22, 2008 
{
	if (!nonbond_SC_contact_list.empty()) nonbond_SC_contact_list.clear();

	int i , j;
	int nn = Residues.size();

	nonbond_SC_contact_list.resize(nn);

	double cut2 = nonbond_DFIRE_cutoff_dist*nonbond_DFIRE_cutoff_dist;
	
	for ( i = 0; i < nn-1; i++)
	{
		HaResidue* rptr1 = Residues[i];
		HaMolecule* mol1 = rptr1->GetHostMol();
		std::string mol_name1 = mol1->GetRef();

		HaAtom* pt1;
		std::string res_name1 = rptr1->GetName();
		if (res_name1 != "GLY") 	
			pt1 = rptr1->GetAtomByName("X");
		else 
			pt1 = rptr1->GetAtomByName("CA");

		double x2 = pt1->GetX();
		double y2 = pt1->GetY();
		double z2 = pt1->GetZ();

		std::set<HaAtom*>& pt_nonb_SC_contact_list = nonbond_SC_contact_list[i];

		for (j = i+1; j < nn; j++)
		{
			HaResidue* rptr2 = Residues[j];
			HaMolecule* mol2 = rptr2->GetHostMol();
			std::string mol_name2 = mol2->GetRef();

			if (mol_name1 == mol_name2) continue; // check if they are from the same molecule

			HaAtom* pt2;
			std::string res_name2 = rptr2->GetName();
			if (res_name2 != "GLY")	
				pt2 = rptr2->GetAtomByName("X");
			else 
				pt2 = rptr2->GetAtomByName("CA");

			double xx2 = pt2->GetX() - x2;
			double r2 = xx2*xx2;
			if (r2 > cut2)
				continue;

			double yy2 = pt2->GetY() - y2;
			r2 += yy2*yy2;
			if (r2 > cut2)
				continue;

			double zz2 = pt2->GetZ() - z2;
			r2 += zz2*zz2;
			if (r2 > cut2)
				continue;

			pt_nonb_SC_contact_list.insert(pt2);
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Build clash atoms list for repulsion forces interactions (in the future this method might be deleted)
//
///////////////////////////////////////////////////////////////////////////////////////////////


bool HaMolMembraneMod::BuildClashAtomList()
{
	if (!nonbond_atom_clash_list.empty()) nonbond_atom_clash_list.clear();

	int i , j;
	int nn = AtomsCentroids.size();
	double cut2 = 0.0;

	nonbond_atom_clash_list.resize(nn);
	
	for ( i = 0; i < nn-1; i++)
	{
		HaAtom* pt1 = AtomsCentroids[i];
		HaMolecule* mol1 = pt1->GetHostMol();
		std::string mol_name1 = mol1->GetRef();

		double sigma1 = pt1->GetVdWRad();

		double x1 = pt1->GetX();
		double y1 = pt1->GetY();
		double z1 = pt1->GetZ();

		std::set<HaAtom*>& pt_nonb_clash_list = nonbond_atom_clash_list[i];

		for (j = i+1; j < nn; j++)
		{
			HaAtom* pt2 = AtomsCentroids[j];
			HaMolecule* mol2 = pt2->GetHostMol();
			std::string mol_name2 = mol2->GetRef();

			if (mol_name1 == mol_name2) continue; // check if they are from the same molecule

			double sigma2 = pt2->GetVdWRad();
			cut2 = sigma1 + sigma2;
			
			cut2 = cut2*cut2;

			double xx2 = pt2->GetX() - x1;
			double r2 = xx2*xx2;
			if (r2 >= cut2)
				continue;

			double yy2 = pt2->GetY() - y1;
			r2 += yy2*yy2;
			if (r2 >= cut2)
				continue;

			double zz2 = pt2->GetZ() - z1;
			r2 += zz2*zz2;
			if (r2 >= cut2)
				continue;

			pt_nonb_clash_list.insert(pt2);
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Load DFIRE_SCM CORE and DFIRE_SCM standard parameters 
//
///////////////////////////////////////////////////////////////////////////////////////////////

int 
HaMolMembraneMod::LoadDFireCoreParams()  //jose modified Nov 4th, rewritten October 22, 2008
{
	if (!la_value.empty())				 la_value.clear();
	if (!la_weight_value.empty())		 la_weight_value.clear(); 
	if (!sc_vdwradius.empty())			 sc_vdwradius.clear(); 
	if (!pairwise_energy_vec.empty())    pairwise_energy_vec.clear();
	if (!pairwise_energy_vec_sa.empty()) pairwise_energy_vec_sa.clear(); 

	std::string fname = "empir_param.dat"; // the file name can be obtained through an interface
	char buf[256], str[256];

	FILE* fp = fopen(fname.c_str(),"r");
	if(fp == NULL)
	{
		ErrorInMod("HaEmpiricalMod::LoadEmpParam",
			"Can not open EMPIR_PARAM.DAT file");
		return FALSE;
	}

	char* pch = fgets(str,255,fp); 
	if( pch == NULL )
	{
		sprintf(buf," Error reading line:  %s ",str);
		ErrorInMod("HaMolMembraneMod::LoadDFireCoreParams()",buf);
		fclose(fp);
		return FALSE;
	}

	enum READ_MODE {LIPID_ACCESS_SCALE, UNRES_PARAMS, PAIRWISE_PARAMS, SA_PAIRWISE_PARAMS} read_mode;

	read_mode = LIPID_ACCESS_SCALE;
	PrintLog ("Reading Lipid Propensity Scale \n");

	for(;;)
	{
		pch = fgets(str,255,fp); 

		std::istrstream is(str);

		if (read_mode == LIPID_ACCESS_SCALE)
		{
			if( strncmp(str,"UNRES_PARAMS",12) == 0)
			{
				read_mode = UNRES_PARAMS;
				PrintLog("Reading Van der Waals radii of Side Chains Force centroids \n");
				continue;
			}

			if( str[0] == '#') continue;
			if( strncmp(str,"DONE",4) == 0) continue;

			std::string rt1s;
			int ich = is.get();
			if (!isspace(ich) && !is.fail()) rt1s += (char) ich;
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt1s += (char) ich;
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt1s += (char) ich;

			double rt_weight;
			double rt_value;
			is >> rt_value;
			is >> rt_weight;

			if(is.fail()) continue;
			//{
			//	ErrorInMod( "HaMolMechMod::InitParamSet()",
			//	            "Error Reading Lipid Propensity Scale");
			//	fclose(fp);
			//	return FALSE;
			//}

			la_value.SetVal(rt1s.c_str(), rt_value);
			la_weight_value.SetVal(rt1s.c_str(), rt_weight);

			PrintLog("%s  Value %5.2f Weigth %5.2f\n",rt1s.c_str(), la_value.GetVal(rt1s.c_str()), la_weight_value.GetVal(rt1s.c_str()));
		}

		if (read_mode == UNRES_PARAMS)
		{
			if( strncmp(str,"PAIRWISE_PARAMS",15) == 0)
			{
				read_mode = PAIRWISE_PARAMS;
				PrintLog ("Reading pairwise energies \n");
				continue;
			} 

			if(str[0] == '#') continue;
			if( strncmp(str,"DONE",4) == 0) continue;

			std::string rt1s;
			int ich = is.get();
			if (!isspace(ich) && !is.fail()) rt1s += (char) ich;
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt1s += (char) ich;
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt1s += (char) ich;

			double rt1_vdw;
			is >> rt1_vdw;

			if(is.fail()) continue;
			//{
			//	ErrorInMod( "HaMolMechMod::InitParamSet()",
			//	            "Error Reading Van der Waals radii");
			//	fclose(fp);
			//	return FALSE;
			//}

			sc_vdwradius.SetVal(rt1s.c_str(), rt1_vdw); // map vdw radius of each residue to their three code letter

			PrintLog ("%s  Side Chain vdW radius %5.2f\n", rt1s.c_str(), sc_vdwradius.GetVal(rt1s.c_str()));		
		}

		if (read_mode == PAIRWISE_PARAMS)
		{
			if( strncmp(str,"SA_PAIRWISE_PARAMS",18) == 0)
			{
				read_mode = SA_PAIRWISE_PARAMS;
				PrintLog ("Read pairwise energies for solvent-accessible residues\n");
				continue;
			}
			
			if(str[0] == '#') continue;
			if( strncmp(str,"DONE",4) == 0) continue;
			
			int i=0;
			HaVec_double pair_vec(20, 0.0);
			std::string rt1s;
			std::string rt2s;
			
			int ich = is.get();
			if (!isspace(ich) && !is.fail()) rt1s += (char) ich;
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt1s += (char) ich;
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt1s += (char) ich;
			
			ich = is.get();
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt2s += (char) ich;
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt2s += (char) ich;
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt2s += (char) ich;
		
			if(rt2s < rt1s)
			{
				std::string tmp = rt2s;
				rt2s = rt1s;
				rt1s = tmp;
			}

			double pair_ene;
			int bin;

			is >> pair_ene;
			is >> bin;
			
			if(is.fail()) continue;
		
			bin -=1;
			pair_vec[bin] = pair_ene;

			i=1;
			while (strncmp(str,"DONE",4)!=0)
			{	
				pch = fgets(str,255,fp);
				std::istrstream is2(str);

				rt1s.clear();
				rt2s.clear();

				i++;
				int ich = is2.get();
				if (!isspace(ich) && !is2.fail()) rt1s += (char) ich;
				ich = is2.get();
				if (!isspace(ich) && !is2.fail()) rt1s += (char) ich;
				ich = is2.get();
				if (!isspace(ich) && !is2.fail()) rt1s += (char) ich;
				
				ich = is2.get();
				ich = is2.get();
				if (!isspace(ich) && !is2.fail()) rt2s += (char) ich;
				ich = is2.get();
				if (!isspace(ich) && !is2.fail()) rt2s += (char) ich;
				ich = is2.get();
				if (!isspace(ich) && !is2.fail()) rt2s += (char) ich;
				
				if(rt2s < rt1s)
				{
					std::string tmp = rt2s;
					rt2s = rt1s;
					rt1s = tmp;
				}

				is2 >> pair_ene;
				is2 >> bin;
				
				//if(is.fail()) continue;
		
				bin -=1;
				pair_vec[bin] = pair_ene;
				std::string res_res = rt1s + "-" + rt2s;
				if (i%20==0)
				{
					pairwise_energy_vec[res_res] = pair_vec;	// Map the residue-residue vector of interactions}
					pair_vec.set(0.0); //clear
				}
			} 
			PrintLog("Number of Pairwise interactions %d\n", pairwise_energy_vec.size());
			//if (i <4199) PrintLog("ERROR. Not complete PAIRWISE potential. Only %d parameters\n", i);
			//else PrintLog("Successfully loaded core pairwise potential\n"); 
		}

		if (read_mode == SA_PAIRWISE_PARAMS)
		{
			if( strncmp(str,"DONE",4) == 0)
			{
				break;
				//PrintLog ("Read pairwise energies for solvent-accessible residues \n");
			}

			if(str[0] == '#') continue;
			if( strncmp(str,"DONE",4) == 0) break;

			int i=0;
			HaVec_double pair_vec(20, 0.0);
			std::string rt1s;
			std::string rt2s;
			
			int ich = is.get();
			if (!isspace(ich) && !is.fail()) rt1s += (char) ich;
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt1s += (char) ich;
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt1s += (char) ich;
			
			ich = is.get();
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt2s += (char) ich;
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt2s += (char) ich;
			ich = is.get();
			if (!isspace(ich) && !is.fail()) rt2s += (char) ich;
		
			if(rt2s < rt1s)
			{
				std::string tmp = rt2s;
				rt2s = rt1s;
				rt1s = tmp;
			}

			double pair_ene;
			int bin;

			is >> pair_ene;
			is >> bin;
			
			if(is.fail()) continue;
		
			bin -=1;
			pair_vec[bin] = pair_ene;

			i=1;
			while (strncmp(str,"DONE",4)!=0)
			{	
				pch = fgets(str,255,fp);
				std::istrstream is2(str);

				rt1s.clear();
				rt2s.clear();

				i++;
				int ich = is2.get();
				if (!isspace(ich) && !is2.fail()) rt1s += (char) ich;
				ich = is2.get();
				if (!isspace(ich) && !is2.fail()) rt1s += (char) ich;
				ich = is2.get();
				if (!isspace(ich) && !is2.fail()) rt1s += (char) ich;
				
				ich = is2.get();
				ich = is2.get();
				if (!isspace(ich) && !is2.fail()) rt2s += (char) ich;
				ich = is2.get();
				if (!isspace(ich) && !is2.fail()) rt2s += (char) ich;
				ich = is2.get();
				if (!isspace(ich) && !is2.fail()) rt2s += (char) ich;
				
				if(rt2s < rt1s)
				{
					std::string tmp = rt2s;
					rt2s = rt1s;
					rt1s = tmp;
				}

				is2 >> pair_ene;
				is2 >> bin;
				
				//if(is.fail()) continue;

				bin -=1;
				pair_vec[bin] = pair_ene;
				std::string res_res = rt1s + "-" + rt2s;
				if (i%20==0)
				{
					pairwise_energy_vec_sa[res_res] = pair_vec;	// Map the residue-residue vector of interactions
					pair_vec.set(0.0); //clear
				}
			} 
//						if (i <4199) PrintLog("ERROR. Not complete PAIRWISE_SA potential. Only %d parameters\n", i);
//						else PrintLog("Successfully loaded solvent-accessible pairwise potential\n\n"); 
		}

		if(pch == NULL) break;
	}

	fclose(fp);
	return TRUE;
}
///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Calculates the DFIRE Score Energy
//
///////////////////////////////////////////////////////////////////////////////////////////////

int
HaMolMembraneMod::ScoreEnergy()
{
	if( module_to_init_flag ) 
		Initialize();
	//if (build_nb_coarsegrained_contact_list) 
	BuildNonBondSCContactList(); // frecuency update
	BuildClashAtomList();

	this->FindAxes();
	pairwise_ene_cg = 0.0;

	vdw_at_repul = 0.0;

	vdw_ene_repul = 0.0;

	lipid_polar_ene = 0.0;

	constraint_ene_mol = 0.0;
	constraint_ene = 0.0; 

	angle_pack = 0.0;

	tot_energy = 0.0;

	int i;
	
	HaAtom* pt1;
	HaAtom* pt2;

	HaResidue* rptr1; 
	HaResidue* rptr2; 

	int nn = AtomsCentroids.size();

	int nr = Residues.size(); //!< jose October 22, 2008. computation of DFIRE Pairwise potential energies 

	if( nr != nonbond_SC_contact_list.size() )
	{
		ErrorInMod("HaMolMechMod::CalculateEnergy()",
		" The number of residues is not equal to the size of nonbonded SC force centroids contact list");
		return FALSE;
	}

	if( nn != nonbond_atom_clash_list.size() )
	{
		ErrorInMod("HaMolMechMod::CalculateEnergy()",
		" The number of atoms CA and SC is not equal to the size of nonbonded atom force centroids clash list");
		return FALSE;
	}

	if (pairwiseDfire_flag)     //!< jose October 22, 2008. computation of DFIRE Pairwise potential energies 
	{
		for ( i =0; i <nr; i++)
		{
			rptr1 = Residues[i];
			std::string res_name1 = rptr1->GetName();
			if (res_name1 != "GLY")
				pt1 = rptr1->GetAtomByName("X");
			else
				pt1 = rptr1->GetAtomByName("CA");

			std::set<HaAtom*>::iterator scitr;
			for( scitr = nonbond_SC_contact_list[i].begin(); scitr != nonbond_SC_contact_list[i].end(); scitr++)
			{
				pt2 = *scitr;
				pairwise_ene_cg += PairwiseEnergy(pt1, pt2);
//				PrintLog("Pairwise energy = %12.6f Pair %s SerNo %d - %d dist= %12.6f \n", pairwise_energy_vec.find(res_res)->second[bin] * 63.69427, res_res.c_str(), rptr1->GetSerNo(), rptr2->GetSerNo(), dist);
			}
		}
	}

	for ( i = 0; i< nn; i++)
	{
		pt1 = AtomsCentroids[i];
		rptr1 = pt1->GetHostRes();
		std::string res_nam1 = rptr1->GetName();
		std::set<HaAtom*>::iterator mitr;
		for ( mitr =  nonbond_atom_clash_list[i].begin(); 
			mitr !=  nonbond_atom_clash_list[i].end(); mitr++)
		{
			pt2 = *(mitr);
			rptr2 = pt2->GetHostRes();
			std::string res_nam2 = rptr2->GetName();
			CalcVdwRep(pt1, pt2, vdw_at_repul);
			vdw_ene_repul += vdw_at_repul;

			//if (res_nam1 == "GLY" && res_nam2 == "GLY")
			//	PrintLog("vdw gly-gly: %4.2f\n", vdw_at_repul);
		}
	}

	double d2= MAX_radius*MAX_radius;
	std::vector<HaAtom* >::iterator aitr = CentreAtoms.begin();
	for (; aitr != CentreAtoms.end(); aitr ++)
	{
		pt1 = (*aitr);
		double temp  = pt1->GetX();
		double r2 = temp * temp;

		temp = pt1->GetY();
		r2 += temp*temp;

		temp = pt1->GetZ();
		r2 += temp*temp;

		if (r2 > d2)
		{
			constraint_ene_mol = (d2 - r2);
			constraint_ene_mol = constraint_ene_mol*constraint_ene_mol;
			constraint_ene += constraint_ene_mol;
		}
	}

	double weight_lipid = 10.0;
	std::vector<HaAtom* >::iterator aitr1 = LipidInterfaceAtoms.begin();
	for (; aitr1 != LipidInterfaceAtoms.end(); aitr1 ++)
	{
		pt1 = (*aitr1);
		HaResidue* rptr1 = pt1->GetHostRes();
		int res_num = rptr1->GetSerNo();
		double temp = pt1->GetZ(); // In this case the helices are parallel to the Y axis
		temp = fabs(temp);

		double width_core1 = 12.9; // angstroms
		double width_core2 = 17.0;
		
		if (res_num == 73)
		{
			if (temp < ( width_core1) )
			{
				//sigma6/(2*r6*r6) - 6.5/sigma6+6.0*sqrt(r2)/sigma7
				//lipid_polar_ene += pow(width_core,6)/(2*pow(temp,12))- 6.5/pow(width_core,6)+6.0*temp/pow(width_core,7);
				//lipid_polar_ene += 1000000000000000000000000.0;
				lipid_polar_ene += weight_lipid*( temp - width_core1 ) *  ( temp - width_core1 );
			}
		}
		else
		{
			if (temp < ( width_core2) )
			{
				//sigma6/(2*r6*r6) - 6.5/sigma6+6.0*sqrt(r2)/sigma7
				//lipid_polar_ene += pow(width_core,6)/(2*pow(temp,12))- 6.5/pow(width_core,6)+6.0*temp/pow(width_core,7);
				//lipid_polar_ene += 1000000000000000000000000.0;
				lipid_polar_ene += weight_lipid*( temp - width_core2 ) *  ( temp - width_core2 );
			}
		}


	}
	//-------------------------------------------------------------------------//
	double weight = 10.0;
	double ux,uy,uz,vx,vy,vz;
	Vec3D vec1 = axis_arr[0];
	Vec3D vec2 = axis_arr[1];
	ux =  vec1[0];
	uy = vec1[1];
	uz = vec1[2];
	
	vx =  vec2[0];
	vy = vec2[1];
	vz = vec2[2];
	double temp = (ux*vx + uy*vy + uz*vz) ;
	anglevar = (acos(temp));
	anglevar *= RAD_TO_DEG; 

	angle_pack = SoftSqrWellPotentiala(anglevar, angle_stat, stdev, weight); 

	//---------------------------------------------------------------------------//

	tot_energy = pairwise_ene_cg + vdw_ene_repul + constraint_ene + lipid_polar_ene + angle_pack;

	if (display_results_flag)
	{
		if (pairwiseDfire_flag)     //!< jose October 22, 2008
		{
			PrintLog(" Pairwise energy = %12.6f kcal/mol \n", pairwise_ene_cg);
		}
		PrintLog(" vdw side chains and alpha carbons centroids = %12.6f Kcal/mol \n", vdw_ene_repul);

		PrintLog(" constraint energy = %12.6f \n", constraint_ene);

		PrintLog(" lipid polar energy = %12.6f \n", lipid_polar_ene);

		PrintLog(" Energy = %12.6f \n", tot_energy);
	}

	return TRUE;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Calculates van der Waals repulsion forces 
//
///////////////////////////////////////////////////////////////////////////////////////////////

bool HaMolMembraneMod::CalcVdwRep( HaAtom* pt1,  HaAtom* pt2, double &vdw_at_ene)
{
	double r2 = 0.0;

	double tmp = pt1->GetX() - pt2->GetX();

	int vdw_4modf =1;

	r2 += tmp*tmp;
	tmp = pt1->GetY() - pt2->GetY();
	r2 += tmp*tmp;
	tmp = pt1->GetZ() - pt2->GetZ();
	r2 += tmp*tmp;

	double sigma2 = pt1->GetVdWRad() + pt2->GetVdWRad();
	double sigma_hardsphere = sigma2* 0.7509; // hard sphere contact from Kussell JMB 2001, 311, 183-193
  	sigma2 = sigma2*sigma2;

	double hardsphere_at_ene = 0.0;

	if (vdw_HardSphere)
	{
		if (r2 < sigma_hardsphere * sigma_hardsphere)
			hardsphere_at_ene = 100000000000000000.0;
	}

	if (vdw_x4_flag) 
	{// calculate repulsion forces between two atoms 
		if (vdw_4modf)
		{
			//sigma2 = sigma2/4.0;
			vdw_at_ene = 0.0157*(sigma2*sigma2)/(r2*r2);//+4*sqrt(r2/sigma2) - 5.0;
		}
		else 
		{
			sigma2 = sigma2/4.0;
			vdw_at_ene = (sigma2*sigma2)/(r2*r2);//+4*sqrt(r2/sigma2) - 5.0;
		}
	}
	if (vdw_x12_flag)
	{
		//weight*(sigma^6/(2*r^12)-6.5/sigma^6 + 6*r/sigma^7)
		double sigma6 = sigma2*sigma2*sigma2;
		double sigma7 = sigma6*sqrt(sigma2);
		double r6 = r2*r2*r2;
		//vdw_at_ene = sigma6/(2*r6*r6) - 6.5/sigma6+6.0*sqrt(r2)/sigma7;
		vdw_at_ene = 0.0157*(sigma6/r6); //0.0157*((sigma)^6/(r^6)-1)
	}
//	PrintLog(" Rep d = %4.3f sigma = %4.3f\n", sqrt(r2), sqrt(sigma2));
//	PrintLog("%s%d %s%d\n", pt1->GetHostRes()->GetName(), pt1->GetHostRes()->GetSerNo(), pt2->GetHostRes()->GetName() ,pt2->GetHostRes()->GetSerNo());

	if (vdw_x4_f_flag)
	{
		double sigma3 = sigma2*sqrt(sigma2);
		vdw_at_ene = sigma2/(2*r2*r2)-2.5/sigma2+2*sqrt(r2)/sigma3;
	}
	// might be added another lennard jones functional form not implemented in HaMolMechMod

	vdw_at_ene += hardsphere_at_ene;
	//PrintLog("%s\n", pt1->GetHostRes()->GetName());
	//HaResidue* rptr1 = pt1->GetHostRes();
	//std::string res_nam1 = rptr1->GetName();
	//HaResidue* rptr2 = pt2->GetHostRes();
	//std::string res_nam2 = rptr2->GetName();
	//if (res_nam1== "GLY" && res_nam2 == "GLY")
	//	PrintLog("Gly-Gly distance = %4.2f\n",sqrt(r2));

	return true;
}
///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Get pairwise interaction 
//
///////////////////////////////////////////////////////////////////////////////////////////////

double
HaMolMembraneMod::PairwiseEnergy(HaAtom* pt1, HaAtom* pt2)
{
	double dist = Vec3D::CalcDistance(pt1, pt2, ANGSTROM_U); // Determination of the bin index 0,1,3...19 span 0 to 15 Ang
	HaResidue* rptr1 = pt1->GetHostRes();
	std::string res_name1 = rptr1->GetName();
	HaResidue* rptr2 = pt2->GetHostRes();
	std::string res_name2 = rptr2->GetName();

	if(res_name2 < res_name1)
	{ 
		std::string tmp = res_name2;
		res_name2 = res_name1;
		res_name1 = tmp;
	}
	int bin=100;
	if(dist <2.0) 
		bin = 0;
	
	if(dist >= 2.0 && dist <8.0)
		bin = (int) (floor((dist-2.0)/0.5)) +1;
	
	if(dist >= 8.0 && dist < nonbond_DFIRE_cutoff_dist)	
		bin = (int)(floor( dist-8.0))+13;
	
	std::string res_res = res_name1 + "-" + res_name2;
	
	if (pairwiseDfire_core)
		return pairwise_energy_vec.find(res_res)->second[bin] *63.69427/3.5;
	else if (pairwiseDfire_sa)
		return pairwise_energy_vec_sa.find(res_res)->second[bin] * 63.69427/3.5;
	else
	{
		ErrorInMod("HaMolMechMod::PairwiseEnergy()",
		" Select one type of pairwise interactions");
		return 0.0;
	}
		
}
///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Loop Limit Sampling (this might be different for each type of membrane protein)
//
///////////////////////////////////////////////////////////////////////////////////////////////

int
HaMolMembraneMod::LoopClosure()
{
	// Loop constraint
	return TRUE;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Calculates SideChain Entropy (possible rotamers library, Dunbrack and Karplus 1993, 1995
//
///////////////////////////////////////////////////////////////////////////////////////////////

int
HaMolMembraneMod::EntropySCM()
{
	// Side chain rotamer entropies
	return TRUE;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Build non-bonded Calpha coarse grained atom list for rigid body calculations 
//  atoms are from different rigid body
//
///////////////////////////////////////////////////////////////////////////////////////////////

bool
HaMolMembraneMod::BuildNonBondCAContactList() 
{
	if (!nonbond_CA_contact_list.empty()) nonbond_CA_contact_list.clear();

	int i , j;
	int nn = Residues.size();

	nonbond_CA_contact_list.resize(nn);

	//double cut2 = nonbond_DFIRE_cutoff_dist*nonbond_DFIRE_cutoff_dist;
	
	for ( i = 0; i < nn-1; i++)
	{
		HaResidue* rptr1 = Residues[i];
		HaMolecule* mol1 = rptr1->GetHostMol();
		std::string mol_name1 = mol1->GetRef();

		HaAtom* pt1;
		std::string res_name1 = rptr1->GetName();
		if (res_name1 != "GLY") 	
			pt1 = rptr1->GetAtomByName("CA");

		/*double x2 = pt1->GetX();
		double y2 = pt1->GetY();
		double z2 = pt1->GetZ();*/

		std::set<HaAtom*>& pt_nonb_CA_contact_list = nonbond_CA_contact_list[i];

		for (j = i+1; j < nn; j++)
		{
			HaResidue* rptr2 = Residues[j];
			HaMolecule* mol2 = rptr2->GetHostMol();
			std::string mol_name2 = mol2->GetRef();

			if (mol_name1 == mol_name2) continue; // check if they are from the same molecule

			HaAtom* pt2;
			std::string res_name2 = rptr2->GetName();
			if (res_name2 != "GLY")	
				pt2 = rptr2->GetAtomByName("CA");

/*			double xx2 = pt2->GetX() - x2;
			double r2 = xx2*xx2;
			if (r2 > cut2)
				continue;

			double yy2 = pt2->GetY() - y2;
			r2 += yy2*yy2;
			if (r2 > cut2)
				continue;

			double zz2 = pt2->GetZ() - z2;
			r2 += zz2*zz2;
			if (r2 > cut2)
				continue;*/

			pt_nonb_CA_contact_list.insert(pt2);
		}
	}
	return true;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Set Coarse-grained OPEP parameters
//
///////////////////////////////////////////////////////////////////////////////////////////////

int HaMolMembraneMod::SetCoarseGrainedOPEPParams()  // jose 11/04/2008 under construction
{
	if ( module_to_init_HaMolMechMod ) 
	{
		MolMechModule->Initialize();
		module_to_init_HaMolMechMod = FALSE;
	}
	//if( module_to_init_flag ) Initialize();

	std::vector<HaAtom*>::iterator aitr;
	for (aitr = MolMechModule->p_mm_model->Atoms.begin(); aitr != MolMechModule->p_mm_model->Atoms.end(); aitr++)
	{
		HaAtom* aptr = (*aitr);
		HaResidue* rptr = aptr->GetHostRes();
		std::string res_name = rptr->GetName();
		std::string res_mod = rptr->GetNameModifier();

		if (res_mod != "OPEP" ) continue; // under decision
		
		std::string at_name = aptr->GetName();
		
		if (at_name =="X" || (res_name == "GLY" && at_name == "CA"))
		{
			if (res_name == "ALA")
			{
				aptr->SetVdWRad(2.287); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);   //< future modification from OPEP radius for each residue
			}
			if (res_name == "ARG")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "ASN")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "ASP")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "GLY")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "GLU")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "HIS")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "LEU")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "LYS")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "ILE")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "MET")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "SER")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "TYR")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "THR")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "PRO")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "CYS")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "GLN")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "VAL")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "TRP")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
			if (res_name == "PHE")
			{
				aptr->SetVdWRad(1.0); //< future modification from OPEP radius for each residue
				aptr->SetVdWEne(0.0);		 //< future modification from OPEP radius for each residue
			}
		}
	}
	
	PrintLog("Set OPEP Bonds params \n");

	std::set<MMBond>::iterator bitr;
	for( bitr = MolMechModule->p_mm_model->MBonds.begin(); bitr != MolMechModule->p_mm_model->MBonds.end(); bitr++)
	{
        MMBond* bptr = (MMBond*)&(*bitr);
		HaAtom* aptr1 = bptr->pt1; 
		HaAtom* aptr2 = bptr->pt2;

		HaResidue* pres1 = aptr1->GetHostRes();
		HaResidue* pres2 = aptr2->GetHostRes();

		//std::string res_name_1 = pres1->GetFullName();
		//std::string res_name_2 = pres2->GetFullName();

		std::string res_mod1 = pres1->GetNameModifier();
		std::string res_mod2 = pres2->GetNameModifier();

		if (res_mod1 != "OPEP" || res_mod2 != "OPEP")	continue;

		std::string at_name1 = aptr1->GetName();
		std::string at_name2 = aptr2->GetName();

		if (( at_name1 == "X" && at_name2 == "CA") || ( at_name1 == "CA" && at_name2 == "X"))
		{
			bptr->set_type = MolMechModel::SET_FF_FIELD;
			bptr->fc = 400.0;	//< Kcal/mol*A^2 Derreumaux JCP, 1999 
			bptr->r0 = 0.0;		//< future modification from OPEP parameter for each residue
		}
	}

	PrintLog("Set OPEP Valence params \n");

	std::set<MMValAngle>::iterator vaitr;
	for( vaitr = MolMechModule->p_mm_model->ValAngles.begin(); vaitr != MolMechModule->p_mm_model->ValAngles.end(); vaitr++)
	{
        MMValAngle* pang = (MMValAngle*) &(*vaitr);
		HaAtom* aptr1 = pang->pt1; 
		HaAtom* aptr2 = pang->pt2;
		HaAtom* aptr3 = pang->pt3;

		HaResidue* pres1 = aptr1->GetHostRes();
		HaResidue* pres2 = aptr2->GetHostRes();
		HaResidue* pres3 = aptr3->GetHostRes();

		//std::string res_name_1 = pres1->GetFullName();
		//std::string res_name_2 = pres2->GetFullName();
		//std::string res_name_3 = pres3->GetFullName();

		std::string res_mod1 = pres1->GetNameModifier();
		std::string res_mod2 = pres2->GetNameModifier();
		std::string res_mod3 = pres3->GetNameModifier();

		if (res_mod1 != "OPEP" || res_mod2 != "OPEP" || res_mod3 != "OPEP")	continue;

		std::string at_name1 = aptr1->GetName();
		std::string at_name2 = aptr2->GetName();
		std::string at_name3 = aptr3->GetName();

		if ((at_name1 == "X" && at_name2 == "CA" && at_name3 == "N") || (at_name1 == "N" && at_name2 == "CA" && at_name3 == "X"))
		{
			pang->set_type = MolMechModel::SET_FF_FIELD;
			pang->fc = 12.0;	//< Kcal/mol*A^2 Derreumaux JCP, 1999
			pang->a0 = 0.0;		//< Degrees, OPEP parameter for each residue
		}
	}

	PrintLog("Set OPEP Dihedral Angles params \n");

	for(std::shared_ptr<MMDihedral> daitr : MolMechModule->p_mm_model->Dihedrals)
	{
		HaAtom* aptr1 = (*daitr).pt1; 
		HaAtom* aptr2 = (*daitr).pt2;
		HaAtom* aptr3 = (*daitr).pt3;
		HaAtom* aptr4 = (*daitr).pt4;

		HaResidue* pres1 = aptr1->GetHostRes();
		HaResidue* pres2 = aptr2->GetHostRes();
		HaResidue* pres3 = aptr3->GetHostRes();
		HaResidue* pres4 = aptr4->GetHostRes();

		//std::string res_name_1 = pres1->GetFullName();
		//std::string res_name_2 = pres2->GetFullName();
		//std::string res_name_3 = pres3->GetFullName();
		//std::string res_name_4 = pres4->GetFullName();

		std::string res_mod1 = pres1->GetNameModifier();
		std::string res_mod2 = pres2->GetNameModifier();
		std::string res_mod3 = pres3->GetNameModifier();
		std::string res_mod4 = pres4->GetNameModifier();

		if (res_mod1 != "OPEP" || res_mod2 != "OPEP" || res_mod3 != "OPEP" || res_mod4 != "OPEP")	continue;

		std::string at_name1 = aptr1->GetName();
		std::string at_name2 = aptr2->GetName();
		std::string at_name3 = aptr3->GetName();
		std::string at_name4 = aptr4->GetName();

		if ((at_name1 == "X" && at_name2 == "N" && at_name3 == "CA" && at_name4 == "C") ||
			(at_name1 == "C" && at_name2 == "CA" && at_name3 == "N" && at_name4 == "X"))
		{
			(*daitr).ClearParams();
			(*daitr).set_type = MolMechModel::SET_FF_FIELD;
			(*daitr).pn.push_back(1.0);
			(*daitr).pk.push_back(3);		// Kcal/mol*rad^2 OPEP parameter Derreumaux JCP, 1999
			(*daitr).idivf.push_back(1.0);
			(*daitr).phase.push_back(0.0);	// OPEP geometric parameter parameter for each residue
		}
	}

	//PrintLog("Set OPEP nonbonded params \n");

	//int i, j;

	//HaAtom* pt1;
	//HaAtom* pt2;

	//HaResidue* rptr1; 
	//HaResidue* rptr2; 

	//int nr = Residues.size();

	//for (i=0; i< nr; i++)
	//{
	//	rptr1 = Residues[i];
	//	std::string res_name1 = rptr1->GetName();

	//	pt1 = rptr1->GetAtomByName("CA");
	//	if res_name1 != "GLY"
	//		pt2 = rptr1->GetAtomByName("X");
	//	

	//}




	//PrintLog("Set OPEP hydrogen-bonded params \n");

	//// future coding

	return TRUE;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Evaluates OPEP potential
//
///////////////////////////////////////////////////////////////////////////////////////////////
int
HaMolMembraneMod::OPEP()
{
	// Evaluates Intermolecular Interactions

	return TRUE;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Any function might be added below 
//
///////////////////////////////////////////////////////////////////////////////////////////////

double
HaMolMembraneMod::SoftSqrWellPotentiala(double& current_value, double& average, double& stdev, double& weight)
{
	double potential =0.0 ;
	double denom = stdev*stdev;
	double up_bound = average + stdev ;
	double low_bound = average - stdev ;
	potential = pow((current_value - average),2)/(denom*denom);

	if (current_value < low_bound)
	{	
		potential += pow((current_value - low_bound),2)/denom;
	}
	if (current_value > up_bound )
	{
		potential +=  pow((up_bound - current_value),2)/denom;
	}
	return potential*weight;
}
///////////////////////////////////////////////////////////////////////////////////////////////
//
//	Any function might be added below 
//
///////////////////////////////////////////////////////////////////////////////////////////////

Vec3DValArray
HaMolMembraneMod::FindAxes()
{ 
	HaAtom* aptr;
	Vec3D cur_com_atom1;
	Vec3D cur_com_atom2;

	double numer_x,numer_y, numer_z;
	double com_x, com_y, com_z;
	Vec3D point;
	int imol, i;
	MolSet* pmset = GetMolSet();
	int nmol = pmset->GetNMol();
	axis_arr.resize(nmol);

	for( imol =0; imol < nmol ; imol++)
	{
		int istart = imol * 8;
		numer_x = 0.0;
		numer_y = 0.0;
		numer_z = 0.0;
		com_x = 0.0;
		com_y = 0.0;
		com_z = 0.0;
		for(i = istart; i < istart+4 ; i++)
		{
			aptr = (HaAtom*) segment_vec[i];
			numer_x += aptr -> GetX();
			numer_y += aptr -> GetY();
			numer_z += aptr -> GetZ();
		}
		com_x = (numer_x*0.25) ;
		com_y = (numer_y*0.25) ;
		com_z = (numer_z*0.25) ;
		
		cur_com_atom1.SetX(com_x);
		cur_com_atom1.SetY(com_y);
		cur_com_atom1.SetZ(com_z);
		
		numer_x = 0.0;
		numer_y = 0.0;
		numer_z = 0.0;
		com_x = 0.0;
		com_y = 0.0;
		com_z = 0.0;
		for(i = istart+4; i <istart+8 ; i++)
		{
			aptr = (HaAtom*) segment_vec[i];
			numer_x += aptr -> GetX();
			numer_y += aptr -> GetY();
			numer_z += aptr -> GetZ();
		}
		com_x = (numer_x*0.25) ;
		com_y = (numer_y*0.25) ;
		com_z = (numer_z*0.25) ;
		
		cur_com_atom2.SetX(com_x);
		cur_com_atom2.SetY(com_y);
		cur_com_atom2.SetZ(com_z);
	

		double length ;
		com_x = cur_com_atom1.GetX() - cur_com_atom2.GetX();
		com_y = cur_com_atom1.GetY() - cur_com_atom2.GetY();
		com_z = cur_com_atom1.GetZ() - cur_com_atom2.GetZ();
		length = 1/(sqrt(com_x*com_x + com_y*com_y +com_z*com_z) );

		point.SetX(com_x*length);
		point.SetY(com_y*length);
		point.SetZ(com_z*length);
		axis_arr[imol].SetCoordFrom(point);
	}

	return axis_arr ;
}

int
HaMolMembraneMod::LineSegments()
{
	MolSet* pmset = GetMolSet();
	wxString name ;
	HaAtom* aptr;
	HaMolecule* pMol;
	HaResidue* res_ptr;
	int nmol = pmset -> GetNMol();
	segment_vec.resize(8*nmol);
	int i,j,k;
	VecPtr temp_vec;
	temp_vec.resize(4);
	k= 0;
	int ires_count =0;
	
	for( int imol = 0; imol < nmol ; imol++)
	{
		pMol = pmset->GetMolByIdx(imol);
		AtomIteratorMolecule aitr(pMol);
		int if_flag = TRUE;
		i=0;
		j = 0;
		HaMolecule::ResidueIterator ritr(pMol);
		for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
		{	
			aptr = (HaAtom*) atm_ca_array[ires_count];
			if (i == 4) i = 0;
			temp_vec[i] = aptr;
			if(j > 3) if_flag = FALSE;
			if (if_flag)
			{
				segment_vec[k]= aptr;
				k++;
			}
			i++;
			j++;
			ires_count++;
		}
		int m;
		for(m= 0 ; m< 4;m++)
		{
			segment_vec[k] = temp_vec[m]; 
			k++;
		}
	}
	return TRUE;
}

int 
HaMolMembraneMod::CbettaSetUp()
{
	MolSet* pmset = GetMolSet();
	char buf[256];
	std::string mol_name;
	std::string res_name;
	std::string atm_name;
	int res_num; 
	wxString atm_ref; 
	HaAtom* aptr ;
	HaAtom* aptr1 ;
	HaResidue* res_ptr ;
	HaMolecule* pMol ;
	HaChain* phost_chain ;
	double dist_unres ;
	char chain_letter ; 
	double vdwradius_unres ;
	char atn[5];
	double atm_mass;
	double cogx, cogy, cogz;
	int natoms;
	int i;
	int sec_count =0;
	int sec_count1 =0;
	HaAtom* aptr_new;
//	fstream out_file;
//	out_file.open("atom_radia.txt",ios::out );
	int nres = pmset -> GetNRes();
	atm_sc_array.resize(nres);
	atm_ca_array.resize(nres);
	int ncont_res = residue_unres_arr.GetCount();

	
	MolSet::ResidueIterator ritr(pmset);
	for(res_ptr = ritr.GetFirstRes(); res_ptr; res_ptr = ritr.GetNextRes())
	{
		
		cogx =0;
		cogy =0;
		cogz =0;
		natoms = 0;
		res_name =res_ptr ->GetName();
		bool new_atm_flag = true;
		if (res_name != "GLY")
		{
			AtomIteratorResidue aitr(res_ptr);
			for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
			{
				atm_name = aptr-> GetName();
				/*if (!(aptr->IsHydrogen()) )
				{
					if (atm_name != "CA" && atm_name != "C" && atm_name != "N" && atm_name != "O" && atm_name != "OXT") //no peptide atoms
					{
						cogx += aptr->GetX();
						cogy += aptr->GetY();
						cogz += aptr->GetZ();
						natoms++;
					}
				}
				if (atm_name == "X")
				{
					new_atm_flag = false;
					HaColor green_c(2, 200, 200);
					aptr -> col= green_c.cidx;
					atm_sc_array[sec_count]=  aptr;
					for ( i=0; i< ncont_res; i++)
					{
						if (residue_unres_arr[i] == res_ptr ->GetName())
						{
							vdwradius_unres = sc_vdwradius[i] ;
						}
					}
					aptr -> radius= vdwradius_unres;
					
					sec_count++;
				}*/
				if (atm_name == "CA")
				{
					atm_ca_array[sec_count1]=  aptr;
					sec_count1++;
				}
			}
		}
		else
		{
			AtomIteratorResidue aitr(res_ptr);
			for(aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom())
			{
				atm_name = aptr-> GetName();
				if (atm_name == "CA")
				{
					cogx = aptr->GetX();
					cogy = aptr->GetY();
					cogz = aptr->GetZ();
					natoms =1;
					atm_ca_array[sec_count1]=  aptr;
					sec_count1++;
				}			
				/*if (atm_name == "X")
				{
					new_atm_flag = false;
					HaColor green_c(2, 200, 200);
					aptr -> col= green_c.cidx;
					atm_sc_array[sec_count]=  aptr;
					
					for ( i=0; i< ncont_res; i++)
					{
						if (residue_unres_arr[i] == res_ptr ->GetName())
						{
							vdwradius_unres = sc_vdwradius[i] ;
						}
					}
					aptr -> radius= vdwradius_unres;
					sec_count++;
				}*/
			}
			
		}	
	}
//	pmset->AnnounceGeomChange();
	return TRUE;
}


