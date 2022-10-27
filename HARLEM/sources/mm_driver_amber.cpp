/*!  \file mm_driver_amber.cpp

    Molecular Mechanics simulations using AMBER package  

    \author Igor Kurnikov 
    \date 2008-

*/

#include <ctime>
#include <exception>

#define HARLEM_MPI 1
#include <mpi.h>

#include <math.h>

#include "haio.h"
#include "hastring.h"

#include "hampi.h"

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
//#include <boost/process.hpp>
//#include <boost/process/child.hpp>

#include <wx/event.h>
#include <wx/filename.h>
#include <wx/process.h>

#include "hawx_add.h"

#include "rapidxml.hpp"

#include "FCMangle.h"

#include "harlemapp.h"
#include "haintcrd.h"
#include "hamolset.h"
#include "hamolmech.h"
#include "mm_elements.h"
#include "mm_model.h"
#include "mm_traj_anal.h"
#include "mm_driver_amber.h"

int MMDriverAmber::stack_limit = 0;

extern "C"
{
extern char FC_FUNC_MODULE(file_io_dat_mod,mdout_name)[80];
//extern char FC_FUNC_MODULE(file_io_dat_mod,mdinfo_name)[80];
extern char FC_FUNC_MODULE(file_io_dat_mod,mdinfo_name)[80];
extern char FC_FUNC_MODULE(file_io_dat_mod,refc_name)[80];
extern char FC_FUNC_MODULE(file_io_dat_mod,mdcrd_name)[80];
extern char FC_FUNC_MODULE(file_io_dat_mod,mdvel_name)[80];
extern char FC_FUNC_MODULE(file_io_dat_mod,mden_name)[80];
extern char FC_FUNC_MODULE(file_io_dat_mod,restrt_name)[80];
extern char FC_FUNC_MODULE(file_io_dat_mod,prmtop_name)[80];

extern char FC_FUNC_MODULE(file_io_dat_mod,error_msg)[300];
extern char FC_FUNC_MODULE(file_io_dat_mod,result_str)[80];

extern int FC_FUNC_MODULE(parallel_dat_mod,lib_mpi_comm);
extern int FC_FUNC_MODULE(parallel_dat_mod,lib_mpi_group);
extern int FC_FUNC_MODULE(parallel_dat_mod,master);
extern int FC_FUNC_MODULE(parallel_dat_mod,numtasks);
extern int FC_FUNC_MODULE(parallel_dat_mod,mytaskid);

extern void FC_FUNC_MODULE(parallel_mod,mpi_gathervec)(int* atm_cnt, double* vec, int* atm_owner_map, int* my_atm_lst, int* my_atm_cnt);
extern void FC_FUNC_MODULE(parallel_mod,mpi_allgathervec)(int* atm_cnt, double* vec, int* atm_owner_map, int* my_atm_lst, int* my_atm_cnt);
extern void FC_FUNC_MODULE(parallel_mod,gb_mpi_gathervec)(int* atm_cnt, double* vec);
extern void FC_FUNC_MODULE(parallel_mod,gb_mpi_allgathervec)(int* atm_cnt, double* vec);
extern void FC_FUNC_MODULE(parallel_mod,distribute_crds_proxy)(int* atm_cnt, int* my_atm_cnt, double* crd);

extern int FC_FUNC_MODULE(prmtop_dat_mod,natom);

extern void FC_FUNC_MODULE(prmtop_dat_mod,alloc_prmtop_mem)();
// extern void FC_FUNC_MODULE(prmtop_dat_mod,alloc_atm_nsp)(int*);
extern void FC_FUNC_MODULE(prmtop_dat_mod,alloc_atm_nsp)();
extern void FC_FUNC_MODULE(prmtop_dat_mod,calc_dihedral_parms)(int*, double*, double*, double*, 
												   double*, double*, int*, double*);

extern void FC_FUNC_MODULE(mdin_amoeba_dat_mod,init_mdin_amoeba_dat)();

extern void open_iunit_debug_(int* irank);

extern void FC_FUNC_MODULE(loadbal_mod,check_new_list_limit)(int* new_list);
extern void FC_FUNC_MODULE(loadbal_mod,do_load_balancing)(int* new_list, int* atm_cnt, double* crd, double* frc,
											  int* atm_owner_map, int* my_atm_lst, int* my_atm_cnt, double* mass);
extern int FC_FUNC_MODULE(loadbal_mod,get_atm_redist_needed)();
extern int FC_FUNC_MODULE(loadbal_mod,get_fft_slab_redist_needed)();

extern void write_mdout_(const char* str,int len);
extern void format_int_fort_(int* i, const char* str_fmt, int len);
extern void format_double_fort_(double* x, const char* str_fmt, int len);

extern void FC_FUNC_MODULE(dynamics_dat_mod,init_dynamics_dat)(int* natom, int* nres, int* nspm, int* res_atms, 
											      int* atm_nsp, double* mass, double* crd);

extern void FC_FUNC_MODULE(dist_constr_mod,dist_constr_setup)( int* use_atm_map, int* atm_owner_map);
extern void FC_FUNC_MODULE(dist_constr_mod,dist_constr_alloc)( int* num_constr );


// prmtop_dat structures

extern void alloc_atm_nsp_(int* n);

extern void set_atm_qterm_(double* vals,int* n);
extern void set_atm_iac_(int* vals,int* n);
extern void set_typ_ico_(int* vals,int* n);
extern void set_gbl_cn1_(double* vals,int* n);
extern void set_gbl_cn2_(double* vals,int* n);
extern void set_gbl_rk_(double* vals,int* n);
extern void set_gbl_req_(double* vals,int* n);
extern void set_gbl_tk_(double* vals,int* n);
extern void set_gbl_teq_(double* vals,int* n);
extern void set_gbl_asol_(double* vals,int* n);
extern void set_gbl_bsol_(double* vals,int* n);
extern void set_gbl_pk_(double* vals,int* n);
extern void set_gbl_pn_(double* vals,int* n);
extern void set_gbl_phase_(double* vals,int* n);
extern void set_gbl_gamc_(double* vals,int* n);
extern void set_gbl_gams_(double* vals,int* n);
extern void set_gbl_fmn_(double* vals,int* n);
extern void set_gbl_ipn_(int* vals,int* n);
extern void set1_atm_igraph_(int* i, const char* vals,int len);
extern void set_gbl_res_atms_(int* vals,int* n);

extern void set_atm_gb_fs_(double* vals,int* n);
extern void set_atm_gb_radii_(double* vals,int* n);

extern void set_atm_numex_(int* vals,int* n);
extern void set_gbl_natex_(int* vals,int* n);
extern void set_prmtop_ititl_(const char* vals, int n);

extern void set_prmtop_int_(int* vals);
extern void set_mdin_ctrl_int_(int* vals);
extern void set_mdin_ctrl_dbl_(double* vals);
extern void FC_FUNC_MODULE(mdin_ewald_dat_mod,set_add_pmemd_int_pars)(int* vals);
extern void FC_FUNC_MODULE(mdin_ewald_dat_mod,set_pme_int_pars_to_fort)(int* vals);
extern void FC_FUNC_MODULE(mdin_ewald_dat_mod,set_pme_dbl_pars_to_fort)(double* dvals);
extern void FC_FUNC_MODULE(mdin_ewald_dat_mod,compute_even_nfft)(double* box, double* fft_grids_per_ang, int* nfft);
extern void FC_FUNC_MODULE(mdin_ewald_dat_mod,compute_nfft)(double* box, double* fft_grids_per_ang, int* nfft);
extern void FC_FUNC_MODULE(mdin_ewald_dat_mod,compute_ew_coeff)(double* ew_coeff,double* es_cutoff, double* dsum_tol);
extern void FC_FUNC_MODULE(mdin_ewald_dat_mod,set_netfrc)(int* netfrc);
extern void derfcfun_(double* y, double* erfc_val);

extern void set1_gbl_bond_(int* i,int* atm_i, int* atm_j, int* parm_idx);
extern void set1_gbl_angle_(int* i,int* atm_i, int* atm_j, int* atm_k,int* parm_idx);
extern void set1_gbl_dihed_(int* i,int* atm_i, int* atm_j, int* atm_k, int* atm_l, int* parm_idx);
extern void set1_gbl_dist_constr_(int* i,int* atm_i,int* atm_j, double* acoef,double* bcoef,int* nb_type);

extern void get_pbc_data_(int* is_orthog_n, double* pbc_alpha_n, double* pbc_beta_n, double* pbc_gamma_n, 
                          double* recip_n, double* ucell_n, double* cut_factor_n, double* reclng_n,
                          double* uc_volume_n, double* uc_sphere_n);

extern void set_pbc_data_(int* is_orthog_n, double* pbc_alpha_n, double* pbc_beta_n, double* pbc_gamma_n, 
                          double* recip_n, double* ucell_n, double* cut_factor_n, double* reclng_n,
                          double* uc_volume_n, double* uc_sphere_n);

extern void FC_FUNC_MODULE(amoeba_bonds_mod,set_bond_list)( int* list_new, int* numbnd );
extern void FC_FUNC_MODULE(amoeba_bonds_mod,set_bond_params)( double* fc, double* eq, int* num_params );
extern void FC_FUNC_MODULE(amoeba_bonds_mod,set_ftable)(double* ftab_coef, int* ftab_degree);
extern void FC_FUNC_MODULE(amoeba_bonds_mod,set_valid_bit)(int* ival);

extern void FC_FUNC_MODULE(amoeba_ureyb_mod,set_ureyb_list)( int* list_new, int* numbnd );
extern void FC_FUNC_MODULE(amoeba_ureyb_mod,set_ureyb_params)( double* fc, double* eq, int* num_params );
extern void FC_FUNC_MODULE(amoeba_ureyb_mod,set_ftable)(double* ftab_coef, int* ftab_degree);
extern void FC_FUNC_MODULE(amoeba_ureyb_mod,set_valid_bit)(int* ival);

extern void FC_FUNC_MODULE(amoeba_reg_angles_mod,set_reg_angles_list)( int* list_new, int* num_angle );
extern void FC_FUNC_MODULE(amoeba_reg_angles_mod,set_reg_angles_params)( double* fc, double* eq, int* num_params );
extern void FC_FUNC_MODULE(amoeba_reg_angles_mod,set_ftable)(double* ftab_coef, int* ftab_degree);
extern void FC_FUNC_MODULE(amoeba_reg_angles_mod,set_valid_bit)(int* ival);

extern void FC_FUNC_MODULE(amoeba_trig_angles_mod,set_trig_angles_list)( int* list_new, int* num_angle );
extern void FC_FUNC_MODULE(amoeba_trig_angles_mod,set_trig_angles_params)( double* fc, double* eq, int* num_params );
extern void FC_FUNC_MODULE(amoeba_trig_angles_mod,set_ftable)(double* ftab_coef, int* ftab_degree);
extern void FC_FUNC_MODULE(amoeba_trig_angles_mod,set_valid_bit)(int* ival);

extern void FC_FUNC_MODULE(amoeba_opbend_angles_mod,set_opbend_angles_list)( int* list_new, int* num_angle );
extern void FC_FUNC_MODULE(amoeba_opbend_angles_mod,set_opbend_angles_params)( double* fc, int* num_params );
extern void FC_FUNC_MODULE(amoeba_opbend_angles_mod,set_ftable)(double* ftab_coef, int* ftab_degree);
extern void FC_FUNC_MODULE(amoeba_opbend_angles_mod,set_valid_bit)(int* ival);

extern void FC_FUNC_MODULE(amoeba_torsions_mod,set_torsions_list)( int* list_new, int* num_torsions );
extern void FC_FUNC_MODULE(amoeba_torsions_mod,set_torsions_params)( double* fc, double* per, double* phase, int* num_params );
extern void FC_FUNC_MODULE(amoeba_torsions_mod,set_valid_bit)(int* ival);

extern void FC_FUNC_MODULE(amoeba_pitorsions_mod,set_pitorsions_list)( int* list_new, int* num_pitorsions );
extern void FC_FUNC_MODULE(amoeba_pitorsions_mod,set_pitorsions_params)( double* fc, double* per, double* phase, int* num_params );
extern void FC_FUNC_MODULE(amoeba_pitorsions_mod,set_valid_bit)(int* ival);

extern void FC_FUNC_MODULE(amoeba_stretch_bend_mod,set_stretch_bend_list)( int* list_new, int* num_strbend );
extern void FC_FUNC_MODULE(amoeba_stretch_bend_mod,set_stretch_bend_params)( double* fc, double* aeq, double* beq1, double* beq2, int* num_params );
extern void FC_FUNC_MODULE(amoeba_stretch_bend_mod,set_valid_bit)(int* ival);

extern void FC_FUNC_MODULE(amoeba_torsion_torsion_mod,set_torsion_torsion_list)( int* list_new, int* num_strbend );
extern void FC_FUNC_MODULE(amoeba_torsion_torsion_mod,set_torsion_torsion_num_params)(int* n);
extern void FC_FUNC_MODULE(amoeba_torsion_torsion_mod,set_torsion_torsion_params)( int* i, int* dim1, int* dim2, double* angle1, double* angle2,
						 double* func, double* dfunc_dangle1, double* dfunc_dangle2, double* dfunc_dangle1_dangle2 );
extern void FC_FUNC_MODULE(amoeba_torsion_torsion_mod,set_valid_bit)(int* ival);

extern void FC_FUNC_MODULE(amoeba_multipoles_mod,set_local_multipoles_list)( double* list_new, int* n);
extern void FC_FUNC_MODULE(amoeba_multipoles_mod,set_chiral_frame_list)( int* list_new, int* n);
extern void FC_FUNC_MODULE(amoeba_multipoles_mod,set_reg_frame_list)( int* list_new, int* n);
extern void FC_FUNC_MODULE(amoeba_multipoles_mod,set_valid_bit)(int* ival);

extern void FC_FUNC_MODULE(amoeba_adjust_mod,set_amoeba_adjust_list)( int* list_new, int* n );
extern void FC_FUNC_MODULE(amoeba_adjust_mod,set_adjust_weights)( double* vdw_weight_new, double* mpole_weight_new, 
													  double* direct_weight_new, double* polar_weight_new, double* mutual_weight_new );
extern void FC_FUNC_MODULE(amoeba_adjust_mod,set_valid_bit)(int* ival);

extern void FC_FUNC_MODULE(amoeba_vdw_mod,set_vdw_types_list)( int* vdw_atom_type_new, int* vdw_atom_parent_new, double* vdw_atom_parent_crd_wt_new, int* n);
extern void FC_FUNC_MODULE(amoeba_vdw_mod,set_vdw_params)( double* vdw_epsilon_new, double* vdw_radius_new, 
											   double* vdw_buf_delta_new, double* vdw_buf_gamma_new, int* n);
extern void FC_FUNC_MODULE(amoeba_vdw_mod,set_valid_bit)(int* ival);

extern void FC_FUNC_MODULE(amoeba_recip_mod,amoeba_recip_final_setup)();

extern void FC_FUNC_MODULE(amoeba_induced_mod,set_atm_polar)(double* polarizability, int* is_polarizable, double* atm_screen_polar, 
	                                             double* damp_polar_coef1, double* damp_polar_coef2, double* damp_polar_rad, double* hpolar, int* n);
extern void FC_FUNC_MODULE(amoeba_induced_mod,set_valid_bit)(int* ival);

extern void FC_FUNC_MODULE(pbc_mod,pressure_scale_pbc_data)(double* box, double* factor);
extern void pressure_scale_crds_proxy_(double* crd, double* mass);
extern void pressure_scale_restraint_crds_proxy_(double* atm_xc, double* mass);

extern void FC_FUNC_MODULE(img_mod,alloc_img_mem)(int* natom);
extern void FC_FUNC_MODULE(nb_pairlist_mod,alloc_nb_pairlist_mem)(int* natom, double* cutlist);

extern void FC_FUNC_MODULE(cit_mod,set_cit_tbl_dims)(double* pbc_box, double* list_cutoff);

extern void FC_FUNC_MODULE(runfiles_mod,flush_mdout)();
extern void FC_FUNC_MODULE(runfiles_mod,close_mdout)();

extern void open_mdcrd_form_();
extern void open_mdvel_form_();
extern void open_mdcrd_bin4_(int* ntb, int* ntwprt);
extern void open_mdvel_bin4_(int* ntwprt);
extern void open_mden_();
extern void open_restart_form_();
extern void open_restart_bin_();
extern void close_mdcrd_();
extern void close_mdvel_();
extern void close_mden_();
extern void close_restart_();

extern void FC_FUNC_MODULE(timers_mod,zero_time)();
extern void FC_FUNC_MODULE(timers_mod,update_time)(int*);

extern void FC_FUNC_MODULE(shake_mod,shake_setup)(int* atm_owner_map );
extern void FC_FUNC_MODULE(shake_mod,shake)(double* frc, double* crd, int* ntb, int* igroup, double* box, double* mass_inv,double* tol);

extern void set_irespa_(int* irespa_new);

extern void FC_FUNC_MODULE(dynamics_mod,vrand_set_velocities)(int* atm_cnt, double* vel, double* mass_inv, 
												  double* temp, int* atm_owner_map );
extern void FC_FUNC_MODULE(dynamics_mod,langevin_setvel)(int* atm_cnt, double* vel, double* frc, double* mass, double* mass_inv, double* dt, 
											 double* temp0, double* gamma_ln, int* atm_owner_map);



extern void FC_FUNC_MODULE(nb_pairlist_mod,check_all_atom_movement)(int* atm_cnt, double* crd, double* saved_crd, double* box, double* saved_box,
														double* skinnb, int* ntp, int* new_list);
extern void FC_FUNC_MODULE(nb_pairlist_mod,check_my_atom_movement)(double* crd, double* saved_crd, double* box, double* saved_box, int* my_atm_lst, int* my_atm_cnt, 
                                                       double* skinnb, int* ntp, int* new_list);

extern void mm_run_();

extern void FC_FUNC_MODULE(master_setup_mod,master_setup)();
extern void FC_FUNC_MODULE(alltasks_setup_mod,alltasks_setup)(double* crd, double* box, double* vel,int* igroup, int* imin, int* ntb, int* atm_owner_map, 
												  int* my_atm_lst, int* my_atm_cnt,double* mass_inv);
extern void FC_FUNC_MODULE(dynamics_mod,all_atom_setvel)(int* natom, double* vel, double* mass_inv, double* tempi);
extern void FC_FUNC_MODULE(random_mod,amrset)(int* ig);

extern void start_mdout_log_(double* box);

extern logical FC_FUNC_MODULE(pme_force_mod,setup_not_done);
extern int FC_FUNC_MODULE(pme_force_mod,irespa);

extern void FC_FUNC_MODULE(pme_setup_mod,final_pme_setup)(int* igroup, double* tranvec );
extern void FC_FUNC_MODULE(pme_force_mod,alloc_pme_force_mem)(int* natom, int* nb_list_length, int* ntypes);

extern void FC_FUNC_MODULE(pme_direct_mod,init_pme_direct_dat)();

extern void FC_FUNC_MODULE(ene_frc_splines_mod,init_ene_frc_splines_dat)();

extern void FC_FUNC_MODULE(pbc_mod,init_pbc)(double* box1, double* box2, double* box3,
			double* box_alpha, double* box_beta, double* box_gamma, double* cut_vdw_with_skin);

extern void FC_FUNC_MODULE(gb_ene_mod,final_gb_setup)( int* natom );
extern void FC_FUNC_MODULE(nmr_calls_mod,bcast_nmr_dat)();

extern double FC_FUNC_MODULE(timers_mod,run_end_cputime);    
extern double FC_FUNC_MODULE(timers_mod,run_setup_end_cputime);
extern void   FC_FUNC_MODULE(timers_mod,profile_cpu)(int* imin, int* igb);
extern void   FC_FUNC_MODULE(timers_mod,init_test_timers)();
extern void   FC_FUNC_MODULE(timers_mod,print_test_timers_mpi)();
extern void   FC_FUNC_MODULE(timers_mod,print_test_timers)();

extern void unlimit_stack_(int* new_limit);

extern void gb_force_proxy_(int* atm_cnt, double* crd, double* frc, double* mass, double* sys_info, int* ncalls,
							int* atm_jrc, double* atm_xc, double* atm_weight, int* igroup, int* belly_atm_cnt, int* natc, 
                            int* atm_owner_map, int* my_atm_lst, int* my_atm_cnt );
extern void pme_force_proxy_(int* atm_cnt, double* crd, double* saved_crd, double* box, double* saved_box, 
							 double* vel, double* frc, double* mass, int* new_list, 
							 int* atm_jrc, double* atm_xc, double* atm_weight, int* igroup, int* natc, double* sys_info,
							 double* virial, double* ekcmt, double* pme_err_est, int* atm_owner_map, int* my_atm_lst, int* my_atm_cnt, 
							 double* tranvec, int* imin );

extern void get_wall_time_(int* current_sec, int* current_usec);

extern void FC_FUNC_MODULE(fft1d_mod,fft1d_create_test) (int* fft_size);
extern void FC_FUNC_MODULE(fft1d_mod,fft1d_destroy_test)();
extern void FC_FUNC_MODULE(fft1d_mod,fft1d_forward_test)(double* fft_array);
extern void FC_FUNC_MODULE(fft1d_mod,fft1d_back_test)   (double* fft_array);

void setup_alloc_error_() 
{ 
	throw std::runtime_error("Setup Allocation Error"); 
}

void mol_mech_error_() 
{
	std::string msg;
	msg.resize(300);
	int i;
	for(i=0; i < 299; i++)
		msg[i] = FC_FUNC_MODULE(file_io_dat_mod,error_msg)[i];
	throw std::runtime_error(msg); 
}

}; // extern "C"

int MMDriverAmber::OpenAmberMDTrajFortran(int iunit, const std::string& fname,  bool write_flag, bool formatted  )
{
	int ires;

	char traj_name_buf[120];
	strncpy(traj_name_buf,fname.c_str(),119);

	char fstat[8];
	strcpy(fstat,"UNKNOWN");
	fstat[7] = 0;

	char file_format[12];
	if( formatted )
	{
		strcpy(file_format,"FORMATTED");
		file_format[9] = 0;
	}
	else
	{
		strcpy(file_format,"UNFORMATTED");
		file_format[11] = 0;
	}

	ires = openforf_(&iunit,traj_name_buf, strlen(traj_name_buf),
		file_format, strlen(file_format),
		fstat, strlen(fstat));

	char title_str[120]= "THIS IS MODIFED MD TRAJECTORY FILE                                          ";
	wrtfstr_(&iunit,title_str,80);

	return TRUE;
}


int MMDriverAmber::WriteCrdToAmberMDTrajFortran( int iunit, PointContainer* pt_cont, bool save_box, bool formatted )
{
	char form_1[20]="(10F8.3)";
	char form_2[20]="(3F8.3)";

	MolSet* pmset = dynamic_cast<MolSet*>(pt_cont);

	int npt = pt_cont->GetNumPt();
	int ncrd = 3*npt;

	HaVec_double crd_grp(3*npt);
	HaVec_double box(3);

	pt_cont->SaveCrdToArray(crd_grp);

	int ires = wrtfarrd_(&iunit,crd_grp.begin(),&ncrd,form_1,strlen(form_1));
	if(!ires) return FALSE;

	if( pmset != NULL && save_box &&  pmset->per_bc->IsSet() )
	{
		box[0] = pmset->per_bc->GetA();
		box[1] = pmset->per_bc->GetB();
		box[2] = pmset->per_bc->GetC();

		int nbx = 3;
		ires = wrtfarrd_(&iunit,box.begin(),&nbx,form_2,strlen(form_2));
		if(!ires) return FALSE;
	}
	return TRUE;
}

int MMDriverAmber::CloseAmberMDTrajFortran(int iunit)
{
	clsforf_(&iunit);
	return TRUE;
}


MMDriverAmber::MMDriverAmber(HaMolMechMod* p_mm_mod_new) 
{
	p_mm_mod = p_mm_mod_new; 
	p_mm_model    = p_mm_mod->p_mm_model;
	p_amber_model = p_mm_model->p_amber_model; 
	pmset = p_mm_mod->GetMolSet();
	p_tm  = new TimerAmber(this);

	to_save_input_files = TRUE;

//		sander_exe_fname = "sander";
    sander_exe_fname = "pmemd";
	amber_version = 10;

	mdcrd = 12;
    mdvel = 13;
	mdin  = 35;
	mdout = 36;

	dirfrc_efs = TRUE;  // corresponds to DIRFRC_EFS in PMEMD_LIB - need to remove DIRFRC_EFS
	emulate_ext_amber = TRUE;

	master   = TRUE;
	numtasks = 1;
	mytaskid = 0;
	my_atm_cnt = 0;

	driver_mpi_comm   = MPI_COMM_NULL;
	driver_mpi_group  = MPI_GROUP_NULL;

	numtasks = 1;
	mytaskid = 0;
	master = 1;

	amber_log_file = "amber_log_file.log";

	sys_info.resize(SI_CNT);
	sys_info = 0.0;

	collect_crds    = TRUE;      
	all_crds_valid  = FALSE;    
	all_vels_valid  = FALSE;

	pbc_box.resize(3);
	pbc_box = 0.0;
	gbl_saved_box.resize(3);
	gbl_saved_box = 0.0;
	is_orthog = TRUE;
	is_octahedral = FALSE;
	pbc_alpha = 0.0;
	pbc_beta  = 0.0;
	pbc_gamma = 0.0;
	recip.resize(9); recip = 0.0;
	ucell.resize(9); ucell = 0.0;
	cut_factor.resize(3); cut_factor = 0.0;
	reclng.resize(3); reclng = 0.0;
	uc_volume = 0.0; 
	uc_sphere = 0.0; 

	tranvec.resize(3*18); 
	tranvec = 0.0;

	trj_coord_fp = NULL;
	trj_vel_fp = NULL;
	trj_ene_fp = NULL;
}

MMDriverAmber::~MMDriverAmber()
{
	delete p_tm;
}

void MMDriverAmber::BcastCrd(MPI_Comm& comm)
{
	int ires;
	int rank;
	int nat = p_amber_model->natom;

	ires = MPI_Comm_rank(comm,&rank);

	if(rank != 0) atm_crd.resize(3*nat);
	ires = MPI_Bcast(atm_crd.v(),3*nat,MPI_DOUBLE,0,comm);
}

void MMDriverAmber::BcastVel(MPI_Comm& comm)
{
	int ires;
	int rank;
	int nat = p_amber_model->natom;

	ires = MPI_Comm_rank(comm,&rank);
	if(rank != 0) atm_vel.resize(3*nat);
	ires = MPI_Bcast(atm_vel.v(),3*nat,MPI_DOUBLE,0,comm);
}

void MMDriverAmber::BcastFrc(MPI_Comm& comm)
{
	int ires;
	int rank;
	int nat = p_amber_model->natom;

	ires = MPI_Comm_rank(comm,&rank);

	if(rank != 0) atm_frc.resize(3*nat);
	ires = MPI_Bcast(atm_frc.v(),3*nat,MPI_DOUBLE,0,comm);
}

void MMDriverAmber::BcastPBox(MPI_Comm& comm)
{
	int ires;
	int rank;

#if defined(HARLEM_MPI)
	ires = MPI_Comm_rank(comm,&rank);

	if(rank == 0) GetPBoxDataFromFortran();

	pmset = p_mm_mod->GetMolSet();

	ires = MPI_Bcast(&pmset->per_bc->orthogonal_flag, 1, MPI_INT,0,comm);
	ires = MPI_Bcast(&pmset->per_bc->octahedron_flag, 1, MPI_INT,0,comm);

	ires = MPI_Bcast(pbc_box.v(),3,MPI_DOUBLE,0,comm);	
	ires = MPI_Bcast(&pbc_alpha,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&pbc_beta,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&pbc_gamma,1,MPI_DOUBLE,0,comm);

	ires = MPI_Bcast(&is_orthog,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&is_octahedral,1,MPI_INT,0,comm);

	ires = MPI_Bcast(recip.v(),9,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(ucell.v(),9,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(cut_factor.v(),3,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(reclng.v(),3,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&uc_volume,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&uc_sphere,1,MPI_DOUBLE,0,comm);
	
	if(rank != 0) SetPBoxDataToFortran();
#endif
}

void MMDriverAmber::SetupMasterNode()
{	
	if (p_amber_model->using_pme_potential) InitPMEParams();
	
	int inerr = FALSE;
     
    if (p_amber_model->using_pme_potential) 
	{
// Set up coordinate index table dimensions; 
// for partitioning of atoms into regions probably ??
// these dimensions will not change on small changes of box dimenstions
// so no need to recompute for small changes of periodical box 
		double tmp = p_amber_model->vdw_cutoff + p_mm_model->skin_nb;
		FC_FUNC_MODULE(cit_mod,set_cit_tbl_dims)(pbc_box.v(),&tmp);
	}

// Init dynamics data; 

   FC_FUNC_MODULE(dynamics_dat_mod,init_dynamics_dat)(&p_amber_model->natom, &p_amber_model->nres, &p_amber_model->nspm, 
	                                      p_amber_model->gbl_res_atms.v(), p_amber_model->atm_nsp.v(), p_amber_model->atm_mass.v(), atm_crd.v());

  if (p_amber_model->using_pme_potential) 
  {
    // Set up image dynamic memory:
    FC_FUNC_MODULE(img_mod,alloc_img_mem)(&p_amber_model->natom);

    // Set up pairlist memory:

	double tmp = p_amber_model->vdw_cutoff + p_mm_model->skin_nb;
	FC_FUNC_MODULE(nb_pairlist_mod,alloc_nb_pairlist_mem)(&p_amber_model->natom, &tmp);
    
    // Set up ewald variables and memory:

    FC_FUNC_MODULE(pme_force_mod,alloc_pme_force_mem)(&p_amber_model->natom, &p_amber_model->next, &p_amber_model->ntypes);
    FC_FUNC_MODULE(pme_direct_mod,init_pme_direct_dat)();
	if(dirfrc_efs) FC_FUNC_MODULE(ene_frc_splines_mod,init_ene_frc_splines_dat)();
  }

// Consistency checking:

  if ( p_mm_mod->period_bcond != p_mm_mod->period_bcond.NO_PERIODICITY && !pmset->per_bc->IsSet() ) 
  {
      PrintLog(" Periodical Box is not set while periodic boundary conditions is not NO_PERIODICITY \n");
      inerr = TRUE;
  }

  if (p_amber_model->using_pme_potential) 
  {
    if (p_amber_model->vdw_cutoff >= pbc_box[0] * 0.5 || p_amber_model->vdw_cutoff >= pbc_box[1] * 0.5 || 
        p_amber_model->vdw_cutoff >= pbc_box[2] * 0.5)
	{
		PrintLog("Error in MMDriverAmber::SetupMasterNode() \n");
        PrintLog(" max cut must be < half smallest box dimension! \n");
		PrintLog(" max cut= %12.6f\n", p_amber_model->vdw_cutoff);
		PrintLog(" pbc_box[0] = %12.6f \n",pbc_box[0]);
		PrintLog(" pbc_box[1] = %12.6f \n",pbc_box[1]);
		PrintLog(" pbc_box[2] = %12.6f \n",pbc_box[2]);
        inerr = TRUE;
	}
  }

// Warnings:

  if (p_amber_model->using_pme_potential && p_amber_model->ibelly > 0) 
  {
	  PrintLog(" Warning: Although EWALD will work with belly \n");
      PrintLog(" (for equilibration), it is not strictly correct' \n");
  }
  
  if (inerr) 
  {
	  throw std::runtime_error(" Error in MM parameters - Stop MM Setup \n");
  }

// All the bond, angle, and dihedral parameters may be changed here as the
// bond, angle, and dihedral arrays are repacked! Note in particular that
// diheda_idx may also be changed.  We also count atoms in the "belly" here,
// which is probably redundant (also done in rgroup() above).

   if (p_amber_model->ibelly > 0) 
   {  
	   if (!p_amber_model->using_pme_potential) 
	   {
		   // The only allowable belly here has just the first belly_atm_cnt atoms
		   // in the moving part.  Confirm this.

		   int i;
		   for(i = p_amber_model->belly_atm_cnt; i < p_amber_model->natom; i++)
		   {
			   if (p_amber_model->atm_igroup[i] != 0) 
			   {
				   std::string error_msg = " When ibelly != 0 and igb != 0, the moving part must \n";
				   error_msg += " be at the start of the molecule, which seems to not be the case! \n";
				   throw std::runtime_error(error_msg);
			   }
		   }
	   }
   }
}

int MMDriverAmber::SetMPICommAllProcs()
{
	int ires, itmp;
//	PrintLog(" MMDriverAmber::SetMPICommAllProcs() pt 1 \n");
	if(driver_mpi_comm  != MPI_COMM_NULL) MPI_Comm_free(&driver_mpi_comm);
	if(driver_mpi_group != MPI_GROUP_NULL) MPI_Group_free(&driver_mpi_group);
	if( pApp->mpi_driver->nprocs == 1)
	{
//		PrintLog(" MMDriverAmber::SetMPICommAllProcs() pt 2 \n");
		driver_mpi_comm  = MPI_COMM_NULL;
		driver_mpi_group = MPI_GROUP_NULL;
		numtasks = 1;
		mytaskid = 0;
		return TRUE;
	}
//	PrintLog(" MMDriverAmber::SetMPICommAllProcs() pt 3 \n");
	
	MPI_Comm_dup(MPI_COMM_WORLD,&driver_mpi_comm);
	ires = MPI_Comm_group(driver_mpi_comm,&driver_mpi_group);
	ires = MPI_Comm_size(driver_mpi_comm,&itmp); numtasks = itmp;
	ires = MPI_Comm_rank(driver_mpi_comm,&itmp); mytaskid = itmp;

	if(p_mm_mod->single_job_comm != MPI_COMM_NULL) MPI_Comm_free(&p_mm_mod->single_job_comm);
	ires = MPI_Comm_dup(MPI_COMM_WORLD,&p_mm_mod->single_job_comm);
	ires = MPI_Comm_rank(p_mm_mod->single_job_comm,&p_mm_mod->single_job_rank);

	return TRUE;
}


int MMDriverAmber::InitCtrlParams()
{
	ntave = 0;

	nsnb = p_mm_mod->nb_list_update_freq;

	t =      p_mm_mod->start_time;
	nrespa = 1;

	if( p_mm_mod->period_bcond == p_mm_mod->period_bcond.CONST_PRES 
		&& p_mm_mod->pressure_reg_method == p_mm_mod->pressure_reg_method.NO_CRD_SCALING)
	{
		p_mm_mod->pressure_reg_method = p_mm_mod->pressure_reg_method.ISOTROP_CRD_SCALING;
	}

	if( p_mm_mod->period_bcond != p_mm_mod->period_bcond.CONST_PRES )
	{
		p_mm_mod->pressure_reg_method = p_mm_mod->pressure_reg_method.NO_CRD_SCALING ;
	}

	ntp   = p_mm_mod->pressure_reg_method;
	ntc  = p_mm_mod->shake_constr;
	jfastw = p_mm_mod->fast_water_method;

	// PMEMD-specific options:

	mdinfo_flush_interval = 60;    // Flush mdinfo every 60 seconds by default
	mdout_flush_interval = 300;    // Flush mdout every 300 seconds by default

	// Retired options that we will complain about, but not fail on.

	dbg_atom_redistribution = 0;   

	loadbal_verbose = 0;     
	next_mdout_flush_sec = 0;

	if(p_mm_mod->press_relax_time <= 0.0)  p_mm_mod->press_relax_time = 0.2;
	if(p_mm_mod->temp_relax_time_solute <= 0.0 ) p_mm_mod->temp_relax_time_solute = 0.2;

	return true;
}

void MMDriverAmber::ResizeCrdVelFrcArrays(int natom)
{
	if(atm_crd.size() != 3*natom) 
	{
		atm_crd.resize(3*natom);
		atm_crd = 0.0;
	}
	if(atm_vel.size() != 3*natom) 
	{
		atm_vel.resize(3*natom);
		atm_vel = 0.0;
	}
	if(atm_last_vel.size() != 3*natom) 
	{
		atm_last_vel.resize(3*natom);
		atm_last_vel = 0.0;
	}
	if(atm_frc.size() != 3*natom) 
	{
		atm_frc.resize(3*natom);
		atm_frc = 0.0;
	}
	if(gbl_atm_saved_crd.size() != 3*natom) 
	{
		gbl_atm_saved_crd.resize(3*natom);
		gbl_atm_saved_crd = 0.0;
	}
	if(gbl_my_atm_lst.size() != natom)
	{
		gbl_my_atm_lst.resize(natom);
		int i;
		for(i = 0; i < natom; i++)
		{
			gbl_my_atm_lst[i] = i+1;
		}
	}
	if(gbl_atm_owner_map.size() != natom)
	{
		gbl_atm_owner_map.resize(natom);
		int i;
		for(i = 0; i < natom; i++)
		{
			gbl_atm_owner_map[i] = 0;
		}
	}
}

void MMDriverAmber::SetupShakePars()
{
	FC_FUNC_MODULE(shake_mod,shake_setup)( gbl_atm_owner_map.v() );
}

void MMDriverAmber::SetupPME()
{
	FC_FUNC_MODULE(pme_setup_mod,final_pme_setup)(p_amber_model->atm_igroup.v(), tranvec.v() );
	if( p_mm_model->IsAmoebaFF() ) FC_FUNC_MODULE(amoeba_recip_mod,amoeba_recip_final_setup)();
}

void MMDriverAmber::SetupGB()
{
	FC_FUNC_MODULE(gb_ene_mod,final_gb_setup)( &p_amber_model->natom );
}

void MMDriverAmber::ValidateCntrParams()
{
	int inerr = FALSE;
	std::string error_msg;
	  
	if  ((p_mm_mod->init_read_coord == p_mm_mod->init_read_coord.READ_X_BIN || 
		  p_mm_mod->init_read_coord == p_mm_mod->init_read_coord.READ_X_FORM) && p_mm_mod->restart_flag != 0 )  
	{
		inerr = TRUE;
		error_msg += "init_read_coord and restart_flag are inconsistent!\n";
	}


	if (p_mm_mod->wrt_vel_freq < -2) 
	{
		inerr = TRUE;
		error_msg += "wrt_vel_freq must be >= -1!\n";
	}

	if (p_mm_mod->wrt_vel_freq == -1 && p_mm_mod->wrt_coord_freq == 0) 
	{
		inerr = TRUE;
		error_msg += "wrt_vel_freq may be -1 only if wrt_coord_freq > 0!\n";
	}

	if (p_mm_mod->limit_wrt_atoms < 0) 
	{
		inerr = TRUE;
		error_msg += "limit_wrt_atoms must be >= 0!\n";
	}

	if( p_amber_model->scee == 0.0)
	{
		inerr = TRUE;
		error_msg += "scee must != 0.0!\n";
	}

	if (nsnb != 25) 
	{
		inerr = TRUE;
		error_msg += "The nsnb ctrl option does not affect nonbonded list update frequency.\n";
		error_msg += "It does affect steepest descent minimization freq if p_mm_mod->min_type == CONJ_GRAD \n";
	}
 
	if (p_amber_model->igb != 0 && p_amber_model->igb != 1 && p_amber_model->igb != 2 && 
		p_amber_model->igb != 5 && p_amber_model->igb != 7 && p_amber_model->igb != 20 ) 
	{
		inerr = TRUE;
		if (p_amber_model->igb == 10) 
		{
			error_msg += "HARLEM does not support MD with Poisson-Boltzmann simulations (igb == 10)! \n";
		}
		else
		{
			error_msg += "igb must be 0,1,2,5,7 or 20!\n";
		}
	}

	if (p_amber_model->alpb != 0 && p_amber_model->alpb != 1) 
	{
		inerr = TRUE;
		error_msg += "alpb must be 0 or 1!\n";
	}

	if (p_amber_model->alpb == 1 && !p_amber_model->using_gb_potential) 
	{
		inerr = TRUE;
		error_msg += "igb must be 1,2,5 or 7 if alpb == 1!\n";
	}

	if (p_amber_model->saltcon < 0.0) 
	{
		inerr = TRUE;
		error_msg += "saltcon must be a positive value!\n";
	}

	if (p_amber_model->rgbmax < 5.0 * p_amber_model->gb_fs_max) 
	{
		inerr = TRUE;
		error_msg += "rgbmax must be at least ";
		error_msg += (wxString::Format("%f\n",5.0 * p_amber_model->gb_fs_max)).c_str();
	}

	if (p_amber_model->rbornstat != 0 && p_amber_model->rbornstat != 1) 
	{
		inerr = TRUE;
		error_msg += "rbornstat must be 0 or 1!\n";
	}

	if (p_amber_model->gbsa != 0) 
	{
		inerr = TRUE;
		error_msg += "HARLEM does not support gbsa > 0 (Generalized Born/surface area simulations).\n";
	}

	if (p_amber_model->ibelly != 0 && p_amber_model->ibelly != 1) 
	{
		inerr = TRUE;
		error_msg += "ibelly must be set to 0 or 1!\n";
	}
	
	if (p_mm_mod->min_type < 0 || p_mm_mod->min_type > 2) 
	{
		inerr = TRUE;
		error_msg += "min_type must be set to 0, 1 or 2!\n";
	}

	if (nrespa < 1) 
	{
		inerr = TRUE;
		error_msg += "nrespa must be >= 1!\n";
	}

	if (nrespai < 1) 
	{
		inerr = TRUE;
		error_msg += "nrespai must be >= 1!\n";
	}

	if (p_mm_mod->temp_control_method < 0 || p_mm_mod->temp_control_method > 3) 
	{
		inerr = TRUE;
		error_msg += "ntt must be set to 0, 1, 2, or 3!\n";
	}

	if (ntp < 0 || ntp > 6) 
	{
		inerr = TRUE;
		error_msg += "ntp must be set to 0 - 6!\n";
	}

	if (ntc < 1 || ntc > 3) 
	{
		inerr = TRUE;
		error_msg += "ntc must be set to 1, 2 or 3!\n";
	}

	if (jfastw != 0 && jfastw != 4) 
	{
		inerr = TRUE;
		error_msg += "jfastw must be set to 0 or 4!\n";
	}

	if (mdinfo_flush_interval < 0) 
	{
		inerr = TRUE;
		error_msg += "mdinfo_flush_interval must be >= 0!\n";
	}
	else if (mdinfo_flush_interval > 60 * 60) 
	{
		inerr = TRUE;
		error_msg += "Excessive mdinfo_flush_interval reset to 1 hour!\n";
	}

	if (dbg_atom_redistribution != 0 && dbg_atom_redistribution != 1) 
	{
		inerr = TRUE;
		error_msg += "dbg_atom_redistribution must be set to 0 or 1!\n";
	}
// Consistency checks:

	if ((nrespa > 1 || nrespai > 1) && p_mm_mod->run_type != p_mm_mod->run_type.MD_RUN ) 
	{
		inerr = TRUE;
		error_msg += "For minimization, nrespa and nrespai must be 1 (default)!\n";
	}

	if (nrespa > 1 && ntp != 0) 
	{
		inerr = TRUE;
		error_msg += "For constant pressure MD, nrespa must be 1!\n";
	}

	if (p_mm_mod->temp_relax_time_solute < p_mm_mod->md_time_step && p_mm_mod->temp_control_method == 1) 
	{
		inerr = TRUE;
		error_msg += "temp_relax_time_solute must be >= md_time_step !\n";
	}

	if (ntp == 0 && p_mm_mod->period_bcond == p_mm_mod->period_bcond.CONST_PRES ) 
	{
		inerr = TRUE;
		error_msg += "ntp must be greater than 0 for CONST_PRESSURE simulations !\n";
	}

	if (ntp != 0 && p_mm_mod->period_bcond != p_mm_mod->period_bcond.CONST_PRES ) 
	{
		inerr = TRUE;
		error_msg += "ntp must be 0 for no CONST_PRESSURE simulations \n";
	} 

	if (p_mm_mod->press_relax_time < p_mm_mod->md_time_step && ntp != 0) 
	{
		inerr = TRUE;
		error_msg += "press_relax_time must be >= dt (step size)!\n";
	}

	if (p_mm_mod->period_bcond != p_mm_mod->period_bcond.NO_PERIODICITY && p_amber_model->igb > 0) 
	{
		inerr = TRUE;
		error_msg += "igb > 0 is only compatible NO_PERIODICITY boundary conditions \n";
	}

	if (p_amber_model->igb == 0)
	{
		// User must either specify no cutoffs, or cut only, or both es_cutoff and
		// vdw_cutoff.  If the user specified no cutoffs, then all values were set
		// to defaults in the input subroutine.  If es_cutoff and vdw_cutoff are
		// specified, es_cutoff must not be greater than vdw_cutoff.  We actually do
		// allow specification of cut, es_cutoff and vdw_cutoff if they are all
		// equal.
		if (p_mm_model->nb_cut_dist == 0.0) 
		{
			if (p_amber_model->es_cutoff == 0.0 || p_amber_model->vdw_cutoff == 0.0)
			{
				inerr = TRUE;
				error_msg += "Both es_cutoff and vdw_cutoff must be specified!\n";
			}
			else if (p_amber_model->es_cutoff > p_amber_model->vdw_cutoff) 
			{
				inerr = TRUE;
				error_msg += "vdw_cutoff must be greater than es_cutoff!\n";
			}
		}
		else
		{
			if ( fabs( p_amber_model->es_cutoff  - p_mm_model->nb_cut_dist) > 0.001 || 
				 fabs( p_amber_model->vdw_cutoff - p_mm_model->nb_cut_dist) > 0.001  )
			{
				inerr = TRUE;
				error_msg += "If cut is used then es_cutoff and vdw_cutoff must be consistent!\n";
			}
		}
		if (p_amber_model->es_cutoff <= 0.0 || p_amber_model->vdw_cutoff <= 0.0) 
		{
			inerr = TRUE;
			error_msg += "Cutoffs (cut, es_cutoff, vdw_cutoff) must be positive!\n";
		}
	}
	else if (p_amber_model->igb == 1 || p_amber_model->igb == 2 || p_amber_model->igb == 5 
		     || p_amber_model->igb == 7 || p_amber_model->igb == 20) 
	{
		if (p_amber_model->gb_cutoff < 8.05) 
		{
			inerr = TRUE;
			error_msg += "cut for Generalized Born simulation too small!\n";
		}
	}

	if( p_mm_model->IsAmoebaFF() ) 
	{
		if (p_amber_model->ibelly != 0) 
		{
			inerr = TRUE;
			error_msg += "Use of belly dynamics (ibelly != 0) with the \n" ;
			error_msg += "Amoeba forcefield is not yet supported. \n";
		}

		if (p_amber_model->beeman_integrator != 0) 
		{
			if (p_mm_mod->wrt_ener_freq != 0) 
			{
				inerr = TRUE;
				error_msg += "Use of mden file output (wrt_ener_freq != 0) with the \n";
				error_msg += "Beeman integrator is not supported. \n";
			}

			if ( p_mm_mod->wrap_coord != 0) 
			{
				inerr = TRUE;
				error_msg += "Use of coordinate wrapping (wrap_coord != 0) with the \n";
				error_msg += "Beeman integrator is not yet supported. \n";
			}

			if (p_mm_mod->limit_wrt_atoms != 0) 
			{
				inerr = TRUE;
				error_msg += "Use of the crd/vel archive limit flag (limit_wrt_atoms != 0) with the \n";
				error_msg += "Beeman integrator is not yet supported. \n";
			}

			if (nrespa != 1) 
			{
				inerr = TRUE;
				error_msg += "Use of respa (nrespa != 1) with the \n";
				error_msg += "Beeman integrator is not yet supported. \n";
			}

			if (ntp > 1) 
			{
				inerr = TRUE;
				error_msg += "Use of nonisotropic pressure scaling with the \n";
				error_msg += "Beeman integrator is not yet supported. \n";
				error_msg += "Please use ntp < 2 with the Beeman integrator. \n";
			}

			if (p_mm_mod->temp_control_method != 0 && p_mm_mod->temp_control_method != 1) 
			{
				inerr = TRUE;
				error_msg += "Use of ntt > 1 (Andersen thermostat or Langevin Dynamics) with the \n";
				error_msg += "Beeman integrator is not yet supported. \n";
			}
		}
	}

	if (inerr == TRUE) 
	{
		PrintLog(" MMDriverAmber::ValidateCntrParams() error_msg %s \n",error_msg.c_str());
		throw std::runtime_error(error_msg);
	}
}

void MMDriverAmber::PrintMDCntrData()
{
  if ( pmset->per_bc->IsSet() && pmset->per_bc->IsOrthogonal() ) PrintLogMDOUT("BOX TYPE: RECTILINEAR,\n"); 
  if ( pmset->per_bc->IsSet() && pmset->per_bc->IsOctahedron() ) PrintLogMDOUT("BOX TYPE: TRUNCATED OCTAHEDRON,\n");  
  if ( pmset->per_bc->IsSet() && !pmset->per_bc->IsOctahedron() && !pmset->per_bc->IsOrthogonal() ) PrintLogMDOUT("BOX TYPE: GENERAL,\n"); 
  
 // Print data characterizing the md-run. 

 // Print control data header:

  PrintLogMDOUT("-----------------------------------------------------------------------------\n");
  PrintLogMDOUT("   2.  CONTROL  DATA  FOR  THE  RUN\n");
  PrintLogMDOUT("-----------------------------------------------------------------------------\n\n");
  PrintLogMDOUT("%s\n",title.c_str());

	PrintLogMDOUT("General flags:\n");
	PrintLogMDOUT("     imin  = %8d  \n", p_mm_mod->run_type.value() );
	PrintLogMDOUT("Nature and format of input:\n");
	PrintLogMDOUT("     ntx   = %8d   irest  = %8d \n",p_mm_mod->init_read_coord.value(),p_mm_mod->restart_flag );
	PrintLogMDOUT("Nature and format of output:\n");
	PrintLogMDOUT(" Rst file fmt = %s   ntpr   = %8d  ntwr  = %8d \n",p_mm_mod->write_coord_format.label(),p_mm_mod->wrt_log_freq,p_mm_mod->wrt_rstrt_freq );
	PrintLogMDOUT("     iwrap = %8d   ntwx   = %8d  ntwv  = %8d ntwe  = %8d\n",p_mm_mod->wrap_coord,p_mm_mod->wrt_coord_freq,p_mm_mod->wrt_vel_freq,p_mm_mod->wrt_ener_freq);
	PrintLogMDOUT(" Trj file fmt = %s   ntwprt   = %7d   rbornstat  = %6d \n",p_mm_mod->traj_wrt_format.label(), p_mm_mod->limit_wrt_atoms, p_amber_model->rbornstat);
	PrintLogMDOUT("Potential function:\n");
	PrintLogMDOUT("     ntf   = %8d   ntb    = %8d  igb   = %8d nsnb  = %8d \n",p_amber_model->ntf,p_mm_mod->period_bcond.value(),p_amber_model->igb,nsnb);
	PrintLogMDOUT("     gbsa    = %8d\n",p_amber_model->gbsa);

	if (p_amber_model->igb == 0) 
	{
		// If there is only one cutoff, we use the old format for printout, just
		// to create fewer output deltas; otherwise a new format shows the different
		// values for es_cutoff and vdw_cutoff:

		if (p_amber_model->es_cutoff == p_amber_model->vdw_cutoff) 
		{
			PrintLogMDOUT("     dielc = %10.5f  vdw_cutoff = %10.5f  intdiel = %10.5f \n",
				p_amber_model->dielc,p_amber_model->vdw_cutoff,p_amber_model->intdiel);
		}
		else
		{
			PrintLogMDOUT("     es_cutoff = %10.5f  vdw_cutoff = %10.5f \n",p_amber_model->es_cutoff,p_amber_model->vdw_cutoff);
			PrintLogMDOUT("     dielc = %10.5f  intdiel = %10.5f \n",p_amber_model->dielc,p_amber_model->intdiel);
		}
	}
	else if (p_amber_model->igb == 1 || p_amber_model->igb == 2 || 
		     p_amber_model->igb == 5 || p_amber_model->igb == 7 || p_amber_model->igb == 20) 
	{
		PrintLogMDOUT("     dielc = %10.5f  gb_cutoff = %10.5f  intdiel = %10.5f \n",
			p_amber_model->dielc,p_amber_model->gb_cutoff,p_amber_model->intdiel);
		PrintLogMDOUT("     saltcon = %10.5f  offset = %10.5f  gb_alpha = %10.5f \n",
			p_amber_model->saltcon,p_amber_model->offset,p_amber_model->gb_alpha);
		PrintLogMDOUT("     gb_beta = %10.5f  gb_gamma = %10.5f  surften = %10.5f \n",
			p_amber_model->gb_beta,p_amber_model->gb_gamma,p_amber_model->surften);
		PrintLogMDOUT("     rgbmax = %10.5f  \n",p_amber_model->rgbmax);
		PrintLogMDOUT("     alpb  = %8d \n",p_amber_model->alpb);

		if (p_amber_model->alpb != 0) PrintLogMDOUT("     arad = %10.5f  \n",p_amber_model->arad);

		PrintLogMDOUT("     bbox_xmin = %10.5f  bbox_ymin = %10.5f  bbox_zmin = %10.5f \n",
			          p_amber_model->bbox_xmin,p_amber_model->bbox_ymin,p_amber_model->bbox_zmin);
		PrintLogMDOUT("     bbox_xmax = %10.5f  bbox_ymax = %10.5f  bbox_zmax = %10.5f \n",
			          p_amber_model->bbox_xmax,p_amber_model->bbox_ymax,p_amber_model->bbox_zmax);
	}

	PrintLogMDOUT("     scnb = %10.5f  scee = %10.5f \n",p_amber_model->scnb,p_amber_model->scee);

	PrintLogMDOUT("Frozen or restrained atoms:\n");
	int ntr_loc = 0; if( p_amber_model->natc > 0) ntr_loc = 1;
	PrintLogMDOUT("     ibelly  = %8d   ntr  = %8d\n",p_amber_model->ibelly, ntr_loc );

	if ( p_mm_mod->run_type == p_mm_mod->run_type.MIN_RUN ) 
	{
		PrintLogMDOUT("Energy minimization:\n");
    // print inputable variables applicable to all minimization methods.
		PrintLogMDOUT("     maxcyc  = %8d ncyc    = %8d   ntmin  = %8d \n",
			          p_mm_mod->max_num_minim_steps,p_mm_mod->num_steep_descent_steps,p_mm_mod->min_type.value() );
		PrintLogMDOUT("     dx0 = %10.5f  drms = %10.5f \n",p_mm_mod->init_min_step,p_mm_mod->grad_cnvrg_val);
	}
	else
	{
		PrintLogMDOUT("Molecular dynamics:\n");
		PrintLogMDOUT("     nstlim = %8d nscm    = %8d   nrespa  = %8d \n",p_mm_mod->num_md_steps,p_mm_mod->remove_rb_motion_freq,nrespa);
		PrintLogMDOUT("     dt     = %10.5f  vlimit = %10.5f  \n",p_mm_mod->md_time_step,p_mm_mod->vel_limit);
   
		if (p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_BERENDSEN ) 
		{
			PrintLogMDOUT("Berendsen (weak-coupling) temperature regulation:\n");
			PrintLogMDOUT("     temp0 = %10.5f  tempi = %10.5f  tautp = %10.5f \n",p_mm_mod->ref_temp,p_mm_mod->init_temp,p_mm_mod->temp_relax_time_solute);
		}
		else if (p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_RANDOMIZED)
		{
			PrintLogMDOUT("Anderson (strong collision) temperature regulation:\n");
			PrintLogMDOUT("     ig = %8d     vrand    = %8d \n",p_mm_mod->random_seed,p_mm_mod->rand_vel_freq);
			PrintLogMDOUT("     temp0 = %10.5f  tempi = %10.5f  tautp = %10.5f \n",p_mm_mod->ref_temp,p_mm_mod->init_temp,p_mm_mod->temp_relax_time_solute);
		}
		else if (p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN) 
		{
			PrintLogMDOUT("Langevin dynamics temperature regulation:\n");
			PrintLogMDOUT("     ig = %8d \n",p_mm_mod->random_seed);
			PrintLogMDOUT("     temp0 = %10.5f  tempi = %10.5f  gamma_ln = %10.5f \n",p_mm_mod->ref_temp,p_mm_mod->init_temp,p_mm_mod->langevin_dump_const);
		}

		if (ntp != 0) 
		{
			PrintLogMDOUT("Pressure regulation:\n");
			PrintLogMDOUT("     ntp = %8d \n",ntp);
			PrintLogMDOUT("     pres0 = %10.5f  comp = %10.5f  taup = %10.5f \n",
				p_mm_mod->ref_pressure,p_mm_mod->compressibility,p_mm_mod->press_relax_time);
		}
	}

	if (ntc != 1) 
	{
		PrintLogMDOUT("SHAKE:\n");
		PrintLogMDOUT("     ntc = %8d     jfastw    = %8d \n",ntc,jfastw);
		PrintLogMDOUT("     tol = %10.5f \n",p_mm_mod->shake_tol);
	}

	if(p_amber_model->using_pme_potential)
	{
		PrintLogMDOUT("Ewald parameters:\n");
		PrintLogMDOUT("     pme_verbose = %8d \n",p_mm_model->pme_verbose );
		PrintLogMDOUT("     vdw_correction_flag = %8d   netfrc  = %8d   \n",p_mm_model->vdw_correction_flag,p_mm_model->subtract_avg_force_flag);
		PrintLogMDOUT("     Box X = %10.5f  Box Y = %10.5f  Box Z = %10.5f \n",pbc_box[0],pbc_box[1],pbc_box[2]);
		PrintLogMDOUT("     Alpha = %10.5f  Beta  = %10.5f  Gamma = %10.5f \n",pbc_alpha,pbc_beta,pbc_gamma);
		PrintLogMDOUT("     pme_grid_nx = %8d   pme_grid_ny = %8d  pme_grid_nz = %8d   \n",p_mm_model->pme_grid_nx, p_mm_model->pme_grid_ny, p_mm_model->pme_grid_nz);
		PrintLogMDOUT("     Cutoff = %10.5f  Tol   = %10.5f \n",p_amber_model->es_cutoff, p_mm_model->pme_dsum_tol );
		PrintLogMDOUT("     Ewald Coefficient = %10.5f \n",p_mm_model->pme_ew_coeff );
		PrintLogMDOUT("     Interpolation order = %5d \n",p_mm_model->pme_spline_order );
	}
}

void MMDriverAmber::InitPMEParams()
{
	const double tollo = 1.0e-12; // low limit on pme_dsum_tol
	const double tolhi = 1.0e-2;  // high limit on pme_dsum_tol

//	const double ew_coefflo = 0.1; // low limit on ew_coeff
//  const double ew_coeffhi = 0.7; // high limit on ew_coeff

	const double ew_coefflo = 0.1e-6; // low limit on ew_coeff
    const double ew_coeffhi = 1.0; // high limit on ew_coeff

	const int gridlo = 6;    //  low limit on grid size 
	const int gridhi = 256;  //  high limit on grid size

	const int orderlo = 3;   // low limit on spline interpolation order
	const int orderhi = 25;  // high limit on spline interpolation order

	const double skinlo = 0.0; // low limit on skin_nb parameter
	const double skinhi = 5.0; // high limit on skin_nb parameter
  
	const double denslo = 100.0;
	const double denshi = 25000.0;

	wxString error_msg = "Error in MMDriverAmber::InitPMEParams() \n";
	int inerr = FALSE;

	if( pbc_box[0] < 1.0 || pbc_box[1] < 1.0 || pbc_box[2] < 1.0 ||
        pbc_alpha < 1.0 || pbc_beta < 1.0 || pbc_gamma < 1.0) 
	{
		error_msg += wxString::Format("pbc_box = %12.6f %12.6f %12.6f \n",
			                          pbc_box[0],pbc_box[1],pbc_box[2]);
		error_msg += wxString::Format("pbc_angles = %12.6f %12.6f %12.6f \n",
			                          pbc_alpha,pbc_beta,pbc_gamma);
			             
		throw std::runtime_error(error_msg.ToStdString().c_str());
	}

	if (p_mm_model->pme_grid_nx == 0) FC_FUNC_MODULE(mdin_ewald_dat_mod,compute_even_nfft)(&pbc_box[0], &p_mm_model->fft_grids_per_ang, &p_mm_model->pme_grid_nx);
	if (p_mm_model->pme_grid_ny == 0) FC_FUNC_MODULE(mdin_ewald_dat_mod,compute_nfft)(&pbc_box[1], &p_mm_model->fft_grids_per_ang, &p_mm_model->pme_grid_ny);
	if (p_mm_model->pme_grid_nz == 0) FC_FUNC_MODULE(mdin_ewald_dat_mod,compute_nfft)(&pbc_box[2], &p_mm_model->fft_grids_per_ang, &p_mm_model->pme_grid_nz);

  // Assign ewald coefficient.  The input values are checked here to insure
  // we don't have a fp error...

  if (p_mm_model->pme_ew_coeff < 1.0e-6) 
  {
	  if(p_mm_model->pme_dsum_tol < tollo || p_mm_model->pme_dsum_tol > tolhi)
	  {
		 error_msg += wxString::Format("Invalid value of pme_dsum_tol %10.5f\n",p_mm_model->pme_dsum_tol );
		 throw std::runtime_error(error_msg.ToStdString().c_str());
	  }
	  FC_FUNC_MODULE(mdin_ewald_dat_mod,compute_ew_coeff)(&p_mm_model->pme_ew_coeff,&p_amber_model->es_cutoff,&p_mm_model->pme_dsum_tol );
  }
  else
  {
	  if( p_mm_model->pme_ew_coeff < ew_coefflo || p_mm_model->pme_ew_coeff > ew_coeffhi)
	  {
		 error_msg += wxString::Format("Invalid value of pme_dsum_tol %10.5f\n",p_mm_model->pme_dsum_tol);
		 throw std::runtime_error(error_msg.ToStdString().c_str());
	  }
	  double erfc_val;
	  double y = p_mm_model->pme_ew_coeff * p_amber_model->es_cutoff;
	  derfcfun_(&y, &erfc_val);
 	  p_mm_model->pme_dsum_tol = erfc_val / p_amber_model->es_cutoff;
  }

  if( p_mm_model->pme_grid_nx < gridlo || p_mm_model->pme_grid_nx > gridhi )
  {
	inerr = TRUE;
	error_msg +=  wxString::Format("Invalid value of pme_grid_nx %d \n",p_mm_model->pme_grid_nx );
  }
  if( p_mm_model->pme_grid_ny < gridlo || p_mm_model->pme_grid_ny > gridhi )
  {
	inerr = TRUE;
	error_msg +=  wxString::Format("Invalid value of pme_grid_ny %d \n",p_mm_model->pme_grid_ny );
  }
  if( p_mm_model->pme_grid_nz < gridlo || p_mm_model->pme_grid_nz > gridhi )
  {
	inerr = TRUE;
	error_msg +=  wxString::Format("Invalid value of pme_grid_nz %d \n",p_mm_model->pme_grid_nz );
  }

  if( p_mm_model->pme_spline_order < orderlo || p_mm_model->pme_spline_order > orderhi)
  {
	inerr = TRUE;
	error_msg +=  wxString::Format("Invalid value of bspl_order %d \n",p_mm_model->pme_spline_order );
  }

  if( p_mm_model->skin_nb < skinlo || p_mm_model->skin_nb > skinhi)
  {
	inerr = TRUE;
	error_msg +=  wxString::Format("Invalid value of skinnb %12.6f \n",p_mm_model->skin_nb);
  }

  if( p_mm_model->vdw_correction_flag < 0 || p_mm_model->vdw_correction_flag > 1)
  {
	inerr = TRUE;
	error_msg +=  wxString::Format("Invalid value of vdw_correction_flag %d \n",p_mm_model->vdw_correction_flag);
  }

  if( p_mm_model->pme_eedtbdns < denslo || p_mm_model->pme_eedtbdns > denshi)
  {
	inerr = TRUE;
	error_msg +=  wxString::Format("Invalid value of eedtbdns %12.6e \n",p_mm_model->pme_eedtbdns);
  }

  if ( (p_amber_model->vdw_cutoff + p_mm_model->skin_nb) < 6.0)
  {
	inerr = TRUE;
	error_msg +=  wxString::Format("Pairlist cutoff is less than 6.0 Ang  = %12.6f \n",(p_amber_model->vdw_cutoff + p_mm_model->skin_nb));
  }

  if (inerr) 
  {
    throw std::runtime_error(error_msg.ToStdString().c_str());
  }

  SetPMEParsFortran();
}

void MMDriverAmber::InitAddMDCtrlParams()
{
	int global_master = FALSE;
	if( p_mm_mod->run_ti == TRUE )
	{
		if( (mytaskid == 0) && (p_mm_mod->inter_model_rank == 0)) global_master = TRUE;
		if( global_master )
		{
			if(p_mm_mod->p_ti_mod->cur_idx_lmb < 0 || p_mm_mod->p_ti_mod->cur_idx_lmb >= p_mm_mod->p_ti_mod->num_lmb )
			{
				p_mm_mod->p_ti_mod->cur_idx_lmb = 0;
			}
		}
	}
	else
	{
		if( mytaskid == 0) global_master = TRUE;
	}
	if( global_master )
	{
		if( p_mm_mod->run_type == p_mm_mod->run_type.MIN_RUN ) 
		{
			p_mm_model->subtract_avg_force_flag = FALSE;
		}
		else
		{
			p_mm_model->subtract_avg_force_flag = TRUE;
		}

		if(p_mm_mod->remove_rb_motion_freq > 0 && p_mm_mod->period_bcond == p_mm_mod->period_bcond.NO_PERIODICITY )
		{
			ndfmin = 6;   // both translation and rotation com motion removed
			if (p_amber_model->natom == 1) ndfmin = 3;
			if (p_amber_model->natom == 2) ndfmin = 5;
		}
		else if ( p_mm_mod->remove_rb_motion_freq > 0)
		{
			ndfmin = 3;    // just translation com will be removed
		}
		else
		{
			ndfmin = 0;
		}
		if(p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN) CalcNumDegFreedom();
		if (p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN ) ndfmin = 0; // No COM motion removal for LD simulation
		if( ndfmin > p_amber_model->num_deg ) ndfmin = p_amber_model->num_deg - 0.001;
		if( ndfmin > p_amber_model->num_deg_solute ) ndfmin = p_amber_model->num_deg_solute - 0.001;
	}
// This is needed for sander-consistent results
	FC_FUNC_MODULE(random_mod,amrset)(&p_mm_mod->random_seed);

	if(p_mm_mod->run_ti == TRUE || numtasks > 1)
	{
		BCastAddMDCtrlParams(p_mm_mod->single_job_comm);
	}
	SetAddMDCtrParamsFortran();
}

void MMDriverAmber::BCastAddMDCtrlParams(MPI_Comm& comm)
{
	int ires;
	ires = MPI_Bcast(&p_mm_mod->lambda_ti,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_mm_mod->run_ti,1,MPI_INT,0,comm);

	ires = MPI_Bcast(&p_mm_mod->p_ti_mod->ti_sync_freq,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&p_mm_mod->p_ti_mod->klambda_ti,1,MPI_DOUBLE,0,comm);

	ires = MPI_Bcast(&p_mm_mod->p_ti_mod->cur_idx_lmb,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&p_mm_mod->p_ti_mod->max_idx,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&p_mm_mod->p_ti_mod->num_lmb,1,MPI_INT,0,comm);

	ires = MPI_Bcast(&p_mm_model->subtract_avg_force_flag,1,MPI_INT,0,comm);

	ires = MPI_Bcast(&ndfmin,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_amber_model->num_deg,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_amber_model->num_deg_solvent,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_amber_model->num_deg_solute,1,MPI_DOUBLE,0,comm);
	
	ires = MPI_Bcast(&mdinfo_flush_interval,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&mdout_flush_interval,1,MPI_INT,0,comm);

	ires = MPI_Bcast(&p_mm_mod->p_mm_info->nstep,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&p_mm_mod->num_md_steps,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&p_mm_mod->wrt_log_freq,1,MPI_INT,0,comm);
}

void MMDriverAmber::SetAddMDCtrParamsFortran()
{
	FC_FUNC_MODULE(mdin_ewald_dat_mod,set_netfrc)( &p_mm_model->subtract_avg_force_flag );
}

void MMDriverAmber::AllGatherVec(HaVec_double& vec)
{
	FC_FUNC_MODULE(parallel_mod,mpi_allgathervec)(&p_amber_model->natom, vec.v(), gbl_atm_owner_map.v(), 
		                               gbl_my_atm_lst.v(), &my_atm_cnt);
}


void MMDriverAmber::ResetStackLimits()
{
	if( stack_limit == 0)
	{
		unlimit_stack_(&stack_limit);
		if ( master && stack_limit > 0)
		{
//			PrintLog( " Stack usage limited by a hard resource limit of %d bytes!",stack_limit);
//			PrintLog( " If segment violations occur, get your sysadmin to increase the limit. ");
		}
	}
}

static int write_string_array_chuncks(std::ostream& os, StrVec& str_vec, 
									  int field_len, int num_of_fields)
{
	int loc_idx  = 0;
	int glob_idx = 0;

	std::string buf_str; 
	int out_len = field_len * num_of_fields;
	buf_str.resize(out_len);

	int i,j;
	int ntot = str_vec.size();
	for(i= 0 ; i <= ntot; i++)                   
	{                                            
		if( i >= ntot)
		{
			if( loc_idx != 0 )
			{
				os << buf_str.c_str() << std::endl;
			}
			break;
		}
		if(loc_idx == 0)
		{
			for( j = 0; j < out_len; j++)
			{
				buf_str[j] = ' ';
			}
		}

		std::string& cur_str = str_vec[i];
		int copy_len = cur_str.size();
		if(field_len < copy_len) copy_len = field_len;
		
		for( j =0; j < copy_len; j++)
			buf_str[loc_idx + j] = cur_str[j];
		
		loc_idx += field_len;
		if( loc_idx >= out_len )
		{
			os << buf_str.c_str() << endl;
			loc_idx = 0;
		}
	}
	return TRUE;
}


void MMDriverAmber::SetAtomCrdToInternalArrays()
{
	if(p_amber_model->natom != p_mm_model->Atoms.size() )
	{
		ostrstream sstr;
		sstr << "Error in MMDriverAmber::SetAtomCrdToInternalArrays()" << endl;
		sstr << " natom = " << p_amber_model->natom << "not equal Atoms.size() = "; 
		sstr <<	p_mm_model->Atoms.size() << endl;
		throw sstr.str();
	}
	int i;

	ResizeCrdVelFrcArrays(p_amber_model->natom);

	for(i=0; i < p_amber_model->natom; i++)
	{
		atm_crd[3*i]   = p_mm_model->Atoms[i]->GetX_Ang();
		atm_crd[3*i+1] = p_mm_model->Atoms[i]->GetY_Ang();
		atm_crd[3*i+2] = p_mm_model->Atoms[i]->GetZ_Ang();
		if( emulate_ext_amber )
		{
			ModifyFormatVal(atm_crd[3*i]  ,FLOAT_F15_7);
			ModifyFormatVal(atm_crd[3*i+1],FLOAT_F15_7);
			ModifyFormatVal(atm_crd[3*i+2],FLOAT_F15_7);
		}
	}

	if( pmset->per_bc->IsSet())
	{
		pbc_box.resize(3);
		pbc_box[0] = pmset->per_bc->GetA();
		pbc_box[1] = pmset->per_bc->GetB();
		pbc_box[2] = pmset->per_bc->GetC();
		pbc_alpha = pmset->per_bc->GetAlpha()*RAD_TO_DEG;
		pbc_beta = pmset->per_bc->GetBeta()*RAD_TO_DEG;
		pbc_gamma = pmset->per_bc->GetGamma()*RAD_TO_DEG;

		if( emulate_ext_amber )
		{
			ModifyFormatVal( pbc_box[0] ,FLOAT_F15_7);
			ModifyFormatVal( pbc_box[1] ,FLOAT_F15_7);
			ModifyFormatVal( pbc_box[2] ,FLOAT_F15_7);
			ModifyFormatVal( pbc_alpha  ,FLOAT_F15_7);
			ModifyFormatVal( pbc_beta   ,FLOAT_F15_7);
			ModifyFormatVal( pbc_gamma  ,FLOAT_F15_7);
		}

		InitPBC();
	}
}

void MMDriverAmber::GetAtomCrdFromInternalArrays()
{
	if(p_amber_model->natom != p_mm_model->Atoms.size() )
	{
		ostrstream sstr;
		sstr << "Error in MMDriverAmber::GetAtomCrdFromInternalArrays()" << endl;
		sstr << " natom = " << p_amber_model->natom << "not equal Atoms.size() = "; 
		sstr <<	p_mm_model->Atoms.size() << endl;
		throw std::runtime_error(sstr.str());
	}
	int i;

	for(i=0; i < p_amber_model->natom; i++)
	{
		p_mm_model->Atoms[i]->SetX_Ang( atm_crd[3*i]   );
		p_mm_model->Atoms[i]->SetY_Ang( atm_crd[3*i+1] );
		p_mm_model->Atoms[i]->SetZ_Ang( atm_crd[3*i+2] );
	}

	if( pmset->per_bc->IsSet())
	{
		pmset->per_bc->SetBox(pbc_box[0],pbc_box[1],pbc_box[2]);
	}	
}

void MMDriverAmber::InitParallelDatMod()
{
	int ires,itmp;
	if( driver_mpi_comm  == MPI_COMM_NULL )
	{
		numtasks = 1;
		mytaskid = 0;
		master = 1;
	}
	else
	{
		ires = MPI_Comm_size(driver_mpi_comm,&itmp); numtasks = itmp;
		ires = MPI_Comm_rank(driver_mpi_comm,&itmp); mytaskid = itmp;
		master = (mytaskid == 0) ? TRUE: FALSE;
	}
#	if defined(WITH_LIB_PMEMD)
	FC_FUNC_MODULE(parallel_dat_mod,lib_mpi_comm)  = MPI_Comm_c2f(driver_mpi_comm);  
	FC_FUNC_MODULE(parallel_dat_mod,lib_mpi_group) = MPI_Group_c2f(driver_mpi_group);
	FC_FUNC_MODULE(parallel_dat_mod,numtasks) = numtasks;
	FC_FUNC_MODULE(parallel_dat_mod,mytaskid) = mytaskid;
	FC_FUNC_MODULE(parallel_dat_mod,master) = master;
#	endif
}

void AmberMMModel::FindResMolPartition()
{
	HaResidue* pres;
	max_res_size = 0;
	
	gbl_res_atms.clear();  // Array of indexes of first atoms of residues
	amber_residues.clear();
	nat_amber_residues.clear();
	res_labels.clear();

// Partition atoms into AMBER residues

	gbl_res_atms.reserve( pmset->GetNRes() + 1);

	ResidueIteratorMolSet ritr(pmset);
	HaResidue* prev_res = NULL;
	HaMolecule* p_cur_mol = NULL;

	int len,idiff,j;

	int fat = 1;
	
	HaAtom* aptr_add = NULL;

	for( pres = ritr.GetFirstRes(); pres; pres = ritr.GetNextRes() )
	{
		int na_res = pres->size();
		if( na_res == 0) continue;
		
		if( na_res == 1 && aptr_add == NULL )
		{
			HaAtom* aptr = (*pres)[0];
			if( aptr->GetElemNo() == 1 ) // PMEMD do not support residues consisting of one hydrogen atom
			{
				aptr_add = aptr;
				if( prev_res == NULL ) // If previous residue does not exist add hydrogen atoms to the next residue
				{
					continue;
				}
				if( prev_res->IsBonded(pres) ) // If hydrogen atom is bonded to atoms in the previous residue add it to the previous residue
				{
					fat++;
					int last_idx = amber_residues.size() - 1;
					nat_amber_residues[last_idx]++;
					aptr_add = NULL;
					continue;
				}
			}
		}
		
		res_labels.push_back(pres->GetName());
		std::string& last_lbl = res_labels.back();
		if( last_lbl == "HOH") last_lbl = "WAT";
		len = last_lbl.size();
		idiff = 4 - len;
		if( idiff > 0) 
		{
			for(j = 0; j < idiff; j++)
				last_lbl += ' '; 
		}

		gbl_res_atms.push_back(fat);   // 
		amber_residues.push_back(pres);

		int nat_curr_amber_res = pres->size();
		if(aptr_add)  nat_curr_amber_res++;
		if( max_res_size < nat_curr_amber_res ) max_res_size = nat_curr_amber_res;

		fat += nat_curr_amber_res;
		nat_amber_residues.push_back(nat_curr_amber_res);
		aptr_add = NULL;
		prev_res = pres;
	}
	nres = gbl_res_atms.size();
	gbl_res_atms.push_back(natom+1);

// Partition atoms into AMBER molecules

	int i;

	nspm = 0;
	atm_nsp.clear();

	int na = p_mm_model->Atoms.size();
	if(na == 0) return;

	vector<AtomGroup> res_conn(na); // array of atoms in the residue the current atom belongs to 

	int ir;
	for(ir = 0; ir < gbl_res_atms.size() - 1; ir++)
	{
		for(i = gbl_res_atms[ir] - 1; i < gbl_res_atms[ir+1] - 1; i++)
		{
			HaAtom* aptr1 = p_mm_model->Atoms[i];
			for(j = gbl_res_atms[ir] - 1; j < gbl_res_atms[ir+1] - 1; j++)
			{
				if( i == j) continue;
				HaAtom* aptr2 = p_mm_model->Atoms[j];
				res_conn[i].push_back(aptr2);
				res_conn[j].push_back(aptr1);
			}
		}
	}

	HaAtom* aptr;
	HaAtom* aptr2;

	std::vector<PtrSet> part;
	std::vector<HaAtom*> atoms_mol;
	HaVec_int at_mol_idx(na); // indexes of molecules of atoms
	at_mol_idx = -1;

	PtrIntMap atm_idx_map; 
	for(i = 0; i < na; i++)
	{
		aptr = p_mm_model->Atoms[i];
		atm_idx_map[aptr] = i;
	}

	PtrSet at_set_empty;
	part.push_back(at_set_empty);
	PtrSet at_set_2;

	int mol_idx;

	for(i = 0; i < na; i++)
	{
		if( at_mol_idx[i] > -1) continue; // Atom is already assigned to a molecule
		PtrSet& at_set_curr = part.back();  
		mol_idx = part.size() - 1;        
		
		aptr = p_mm_model->Atoms[i];
			
		for(;;)
		{
			if(at_set_curr.size() == 0 )
			{
				at_set_curr.insert(aptr);
				at_set_2 = at_set_curr;
				continue;
			}
			PtrSet::iterator aitr;
			for( aitr = at_set_curr.begin(); aitr != at_set_curr.end(); aitr++)
			{
				aptr = (HaAtom*) (*aitr);
				AtomGroup bonded_atoms;
				aptr->GetBondedAtoms(bonded_atoms);
				for(j = 0; j < bonded_atoms.size(); j++) // Adding to the current molecule all atoms bonded to atoms of the molecule
				{
					aptr2 = bonded_atoms[j];
					at_set_2.insert(aptr2);  
				}
				if( atm_idx_map.count(aptr) > 0) 
				{
					int at_idx = atm_idx_map[aptr];  // Prevent splitting a residue into different molecules 
					for(j = 0; j < res_conn[at_idx].size(); j++) // Adding to the current molecule all atoms belonging to the residues for atoms already in the molecule
					{	
						aptr2 = res_conn[at_idx][j];
						at_set_2.insert(aptr2);
					}
				}
			}
			int n1 = at_set_curr.size();
			int n2 = at_set_2.size();

			// PrintLog("at %d   n1 = %d   n2 = %d \n",i,n1,n2); 
			if( at_set_2.size() == at_set_curr.size()) // if no additional atoms appear in the molecule after iteration set to atoms of the molecule index of this molecule and continue with partitioning
			{
				for(aitr = at_set_curr.begin(); aitr != at_set_curr.end(); aitr++)
				{
					aptr = (HaAtom*) (*aitr);
					if( atm_idx_map.count(aptr) > 0) 
					{
						int at_idx = atm_idx_map[aptr];
						at_mol_idx[at_idx] = mol_idx;
					}
				}
				part.push_back(at_set_empty);
				break;
			}
			else
			{
				at_set_curr = at_set_2;
			}
		} // for(;;) 
	} // for( i = 0; i < na; i++)
		

	mol_idx = 0;
	int at_cnt = 0;
	int non_seq_mols = FALSE;
	for(i = 0; i < na; i++)
	{
		if( at_mol_idx[i] < mol_idx ) non_seq_mols = TRUE;
		if( at_mol_idx[i] != mol_idx || i == (na-1) )
		{
			if( i == (na - 1) ) at_cnt++;
			atm_nsp.push_back(at_cnt);	
			mol_idx++;
			at_cnt = 1;
		}
		else
		{
			at_cnt++;
		}
	}

	if( non_seq_mols )
	{
		PrintLog("MMDriverAmber::FindResMolPartion() WARNING: \n");
		PrintLog("Atoms are not in sequential order with molecule order \n");
	}

	nspm = atm_nsp.size();

	int nsolv_mol = 0; // Number of solvent molecules
    this->n_solute_mol = nspm;
	for( ir = nres - 1; ir >= 0; ir--)
	{	
		HaResidue* pres = amber_residues[ir];
		if( pres->IsSolvent() )
		{
			nsolv_mol++;
			this->n_solute_mol--;
		}
		else
		{
			break;
		}
	}
	n_solute_res = nres - nsolv_mol; // number of solute residues

/*
	if( pmset->per_bc->IsSet() ) // Additional output for periodic boundary conditions and solvent info
	{           
		
		n_solute_mol = 0;
		HaMolecule* prev_mol = NULL;
		int mol_size = 0;
		atm_nsp.clear();
		for( ir = 0; ir <  n_solute_res; ir++) // Determine Solute Molecules
		{
			HaResidue* pres = amber_residues[ir];
			HaMolecule* pmol = pres->GetHostMol();
			mol_size += nat_amber_residues[ir];
			if( prev_mol != pmol || ir == (n_solute_res - 1) )
			{
				n_solute_mol++;
				atm_nsp.push_back(mol_size);
				mol_size = 0;
				prev_mol = pmol;
			}
		}

		for( ir = n_solute_res; ir < nres; ir++) // Determine Solvent Molecules
		{
			mol_size = nat_amber_residues[ir];
			atm_nsp.push_back(mol_size);
		}
	}
	else
	{
		atm_nsp.resize(1);
		atm_nsp[0] = natom;
	}
*/
//	nspm = atm_nsp.size();
}

void AmberMMModel::CalcAddDihParams()
{
	int nptra = gbl_pk.size();

	gbl_gamc.resize(nptra);
	gbl_gams.resize(nptra);
	gbl_ipn.resize(nptra);
	gbl_fmn.resize(nptra);

	FC_FUNC_MODULE(prmtop_dat_mod,calc_dihedral_parms)(&nptra, gbl_pk.begin(), gbl_pn.begin(), gbl_phase.begin(), 
			                               gbl_gamc.begin(), gbl_gams.begin(), gbl_ipn.begin(), 
										   gbl_fmn.begin());
}

void MMDriverAmber::SaveModelToFortran()
{
//	PrintLog(" MMDriverAmber::SaveModelToFortran() pt 1 \n");
	int i;	
	SetMDinCtrlIntFortran();
	SetMDinCtrlDblFortran();
	SetAddIntParsFortran();
	SetPMEParsFortran();
	SetPrmTopIntFortran();

	if( p_mm_model->IsAmoebaFF() ) FC_FUNC_MODULE(mdin_amoeba_dat_mod,init_mdin_amoeba_dat)();

	FC_FUNC_MODULE(prmtop_dat_mod,alloc_prmtop_mem)();
	
	set_prmtop_ititl_(title.c_str(), title.size());
	for( i = 0; i < p_amber_model->natom; i++)
	{
		int ip1 = i+1;
		set1_atm_igraph_(&ip1,p_amber_model->atm_igraph[i].c_str(),4);
	}

	if( !p_mm_model->IsAmoebaFF() ) 
	{
		set_atm_qterm_(p_amber_model->atm_charge.begin(),&p_amber_model->natom);
	}

	set_atm_iac_(p_amber_model->atm_iac.begin(),&p_amber_model->natom);
	set_atm_numex_(p_amber_model->atm_numex.begin(),&p_amber_model->natom);
	
	int nres_p1 = p_amber_model->nres + 1;
	set_gbl_res_atms_(p_amber_model->gbl_res_atms.begin(),&nres_p1);

	set_gbl_rk_(p_amber_model->gbl_rk.begin(),&p_amber_model->numbnd);
	set_gbl_req_(p_amber_model->gbl_req.begin(),&p_amber_model->numbnd);

	set_gbl_tk_(p_amber_model->gbl_tk.begin(),&p_amber_model->numang);
	set_gbl_teq_(p_amber_model->gbl_teq.begin(),&p_amber_model->numang);

	set_gbl_pk_(p_amber_model->gbl_pk.begin(),&p_amber_model->nptra);
	set_gbl_pn_(p_amber_model->gbl_pn.begin(),&p_amber_model->nptra);
	set_gbl_phase_(p_amber_model->gbl_phase.begin(),&p_amber_model->nptra);

	set_gbl_cn1_(p_amber_model->gbl_cn1.begin(),&p_amber_model->nttyp);
	set_gbl_cn2_(p_amber_model->gbl_cn2.begin(),&p_amber_model->nttyp);

// Fix up the typ-ico array for more efficient processing of vdw and 10_12
// interactions.  Basically, if the ico array entry is 0, then vdw processing
// is unnecessary.

	int nt2 = p_amber_model->ntypes*p_amber_model->ntypes;
	HaVec_int iarr(nt2);
	for( i = 0; i < nt2; i++)
	{
		int ic = p_amber_model->typ_ico[i];
		iarr[i] = ic;
		if( ic > 0)
		{
			if( fabs(p_amber_model->gbl_cn2[ic-1]) < 1.0e-8 && fabs(p_amber_model->gbl_cn1[ic-1]) < 1.0e-8)
			{
				iarr[i] = 0;
			}
		}
		else if( ic < 0 )
		{
			ic = -ic;
			if( fabs(p_amber_model->gbl_asol[ic-1]) < 1.0e-8 && fabs(p_amber_model->gbl_bsol[ic-1]) < 1.0e-8)
			{
				iarr[i] = 0;
			}
		}
	}	
	set_typ_ico_(iarr.begin(),&nt2);
		
	for( i = 0; i < p_amber_model->gbl_bond_allocsize; i++)
	{
		int ip1 = i + 1;
		set1_gbl_bond_(&ip1,&p_amber_model->gbl_bond[3*i],&p_amber_model->gbl_bond[3*i+1],
			           &p_amber_model->gbl_bond[3*i+2]);
	}	

	for( i = 0; i < p_amber_model->gbl_angle_allocsize; i++)
	{
		int ip1 = i + 1;
		set1_gbl_angle_(&ip1,&p_amber_model->gbl_angle[4*i],&p_amber_model->gbl_angle[4*i+1],
				        &p_amber_model->gbl_angle[4*i+2], &p_amber_model->gbl_angle[4*i+3] );
	}

	for( i = 0; i < p_amber_model->gbl_dihed_allocsize; i++)
	{
		int ip1 = i + 1;
		set1_gbl_dihed_(&ip1,&p_amber_model->gbl_dihed[5*i],&p_amber_model->gbl_dihed[5*i+1],
		        		&p_amber_model->gbl_dihed[5*i+2], &p_amber_model->gbl_dihed[5*i+3], 
						&p_amber_model->gbl_dihed[5*i+4]);
	}

	set_gbl_natex_(p_amber_model->gbl_natex.begin(),&p_amber_model->next);
	set_gbl_asol_(p_amber_model->gbl_asol.begin(),&p_amber_model->nphb);
	set_gbl_bsol_(p_amber_model->gbl_bsol.begin(),&p_amber_model->nphb);

	HaVec_double farr;

    if( p_amber_model->num_dist_constr > 0 )               
	{                                                          
		for( i = 0; i < p_amber_model->num_dist_constr; i++)
		{
			int ip1 = i+1;
			set1_gbl_dist_constr_(&ip1,&p_amber_model->dist_constr_idx[3*i],&p_amber_model->dist_constr_idx[3*i+1],
				&p_amber_model->dist_constr_params[2*i],&p_amber_model->dist_constr_params[2*i+1],&p_amber_model->dist_constr_idx[3*i+2]);
		}			
	}
	
	if (p_amber_model->using_gb_potential) 
	{
	   set_atm_gb_radii_(p_amber_model->atm_gb_radii.begin(),&p_amber_model->natom);
	   set_atm_gb_fs_(p_amber_model->atm_gb_fs.begin(),&p_amber_model->natom);
	}

	set_gbl_gamc_(p_amber_model->gbl_gamc.begin(),&p_amber_model->nptra);
	set_gbl_gams_(p_amber_model->gbl_gams.begin(),&p_amber_model->nptra);
	set_gbl_ipn_(p_amber_model->gbl_ipn.begin(),&p_amber_model->nptra);
	set_gbl_fmn_(p_amber_model->gbl_fmn.begin(),&p_amber_model->nptra);

	if( p_mm_model->IsAmoebaFF() )
	{
		int iflag = 1;
        
		// Regular bonds
		FC_FUNC_MODULE(amoeba_bonds_mod,set_bond_list)( p_amber_model->gbl_bond_amoeba.v(),&p_amber_model->n_bond_amoeba );
		FC_FUNC_MODULE(amoeba_bonds_mod,set_bond_params)(p_amber_model->bond_amoeba_params[0].v(),p_amber_model->bond_amoeba_params[1].v(),
			                                 &p_amber_model->n_bond_amoeba_params);
		FC_FUNC_MODULE(amoeba_bonds_mod,set_ftable)(p_amber_model->bond_amoeba_ftab_coef.v(),&p_amber_model->bond_amoeba_ftab_degree);
		
		iflag = 1;
		if( p_amber_model->n_bond_amoeba == 0 || p_amber_model->n_bond_amoeba_params == 0 ) iflag = 0;
		FC_FUNC_MODULE(amoeba_bonds_mod,set_valid_bit)(&iflag);

		// Urey-Bradley bonds
		FC_FUNC_MODULE(amoeba_ureyb_mod,set_ureyb_list)( p_amber_model->gbl_bond_urey.v(),&p_amber_model->n_urey_bond );
		FC_FUNC_MODULE(amoeba_ureyb_mod,set_ureyb_params)(p_amber_model->bond_urey_params[0].v(),p_amber_model->bond_urey_params[1].v(),
			                                 &p_amber_model->n_urey_bond_params);
		FC_FUNC_MODULE(amoeba_ureyb_mod,set_ftable)(p_amber_model->bond_urey_ftab_coef.v(),&p_amber_model->bond_urey_ftab_degree);
		
		iflag = 1;
		if( p_amber_model->n_urey_bond == 0 || p_amber_model->n_urey_bond_params == 0 ) iflag = 0;
		FC_FUNC_MODULE(amoeba_ureyb_mod,set_valid_bit)(&iflag);

		// Regular Valence Angles
		FC_FUNC_MODULE(amoeba_reg_angles_mod,set_reg_angles_list)( p_amber_model->gbl_angle_amoeba_reg.v(), &p_amber_model->n_angle_amoeba );
		FC_FUNC_MODULE(amoeba_reg_angles_mod,set_reg_angles_params)(p_amber_model->angle_amoeba_params[0].v(),p_amber_model->angle_amoeba_params[1].v(),
			                                            &p_amber_model->n_angle_amoeba_params);
		FC_FUNC_MODULE(amoeba_reg_angles_mod,set_ftable)(p_amber_model->angle_amoeba_ftab_coef.v(),&p_amber_model->angle_amoeba_ftab_degree);
		
		iflag = 1;
		if( p_amber_model->n_angle_amoeba == 0 || p_amber_model->n_angle_amoeba_params == 0 ) iflag = 0;
		FC_FUNC_MODULE(amoeba_reg_angles_mod,set_valid_bit)(&iflag);

		// Trigonal Angles
		FC_FUNC_MODULE(amoeba_trig_angles_mod,set_trig_angles_list)( p_amber_model->gbl_angle_amoeba_trig.v(),&p_amber_model->n_trig_angles );
		FC_FUNC_MODULE(amoeba_trig_angles_mod,set_trig_angles_params)(p_amber_model->angle_amoeba_params[0].v(),p_amber_model->angle_amoeba_params[1].v(),
			                                              &p_amber_model->n_angle_amoeba_params);
		FC_FUNC_MODULE(amoeba_trig_angles_mod,set_ftable)(p_amber_model->angle_amoeba_ftab_coef.v(),&p_amber_model->angle_amoeba_ftab_degree);
		
		iflag = 1;
		if( p_amber_model->n_trig_angles == 0 || p_amber_model->n_angle_amoeba_params == 0 ) iflag = 0;
		FC_FUNC_MODULE(amoeba_trig_angles_mod,set_valid_bit)(&iflag);

		// Opbend Angles
		FC_FUNC_MODULE(amoeba_opbend_angles_mod,set_opbend_angles_list)( p_amber_model->gbl_opbend_angle.v(),&p_amber_model->n_opbend_angles );
		FC_FUNC_MODULE(amoeba_opbend_angles_mod,set_opbend_angles_params)(p_amber_model->opbend_angle_params.v(), &p_amber_model->n_opbend_angles_params );
		FC_FUNC_MODULE(amoeba_opbend_angles_mod,set_ftable)(p_amber_model->angle_amoeba_ftab_coef.v(),&p_amber_model->angle_amoeba_ftab_degree);
		
		iflag = 1;
		if( p_amber_model->n_opbend_angles == 0 || p_amber_model->n_opbend_angles_params == 0 ) iflag = 0;
		FC_FUNC_MODULE(amoeba_opbend_angles_mod,set_valid_bit)(&iflag);

		// Torsions
		FC_FUNC_MODULE(amoeba_torsions_mod,set_torsions_list)( p_amber_model->gbl_amoeba_tors_angle.v(),&p_amber_model->n_tors_amoeba );
		FC_FUNC_MODULE(amoeba_torsions_mod,set_torsions_params)(p_amber_model->tors_amoeba_params[0].v(),p_amber_model->tors_amoeba_params[1].v(),
			                                        p_amber_model->tors_amoeba_params[2].v(),&p_amber_model->n_tors_amoeba_params );
		
		iflag = 1;
		if( p_amber_model->n_tors_amoeba == 0 || p_amber_model->n_tors_amoeba_params == 0 ) iflag = 0;
		FC_FUNC_MODULE(amoeba_torsions_mod,set_valid_bit)(&iflag);

		// PI Torsions
		FC_FUNC_MODULE(amoeba_pitorsions_mod,set_pitorsions_list)( p_amber_model->gbl_pi_tors_angle.v(),&p_amber_model->n_pi_torsions );
		FC_FUNC_MODULE(amoeba_pitorsions_mod,set_pitorsions_params)(p_amber_model->pi_tors_params[0].v(),p_amber_model->pi_tors_params[1].v(),
			                                 p_amber_model->pi_tors_params[2].v(),&p_amber_model->n_pi_torsions_params);
		
		iflag = 1;
		if( p_amber_model->n_pi_torsions == 0 || p_amber_model->n_pi_torsions_params == 0 ) iflag = 0;
		FC_FUNC_MODULE(amoeba_pitorsions_mod,set_valid_bit)(&iflag);

		// Stretch-bend terms
		FC_FUNC_MODULE(amoeba_stretch_bend_mod,set_stretch_bend_list)( p_amber_model->gbl_str_bend_angle.v(),&p_amber_model->n_stretch_bend );
		FC_FUNC_MODULE(amoeba_stretch_bend_mod,set_stretch_bend_params)(p_amber_model->str_bend_params[0].v(),p_amber_model->str_bend_params[1].v(),
			                                 p_amber_model->str_bend_params[2].v(),p_amber_model->str_bend_params[3].v(),&p_amber_model->n_stretch_bend_params);
		
		iflag = 1;
		if( p_amber_model->n_stretch_bend == 0 || p_amber_model->n_stretch_bend_params == 0 ) iflag = 0;
		FC_FUNC_MODULE(amoeba_stretch_bend_mod,set_valid_bit)(&iflag);

		// Torsion Torsion terms
		FC_FUNC_MODULE(amoeba_torsion_torsion_mod,set_torsion_torsion_list)( p_amber_model->gbl_tors_tors.v(),&p_amber_model->n_tors_tors );
		FC_FUNC_MODULE(amoeba_torsion_torsion_mod,set_torsion_torsion_num_params)(&p_amber_model->n_tors_tors_params);
		
		for(i = 0; i < p_amber_model->n_tors_tors_params; i++)
		{
			int ip1 = i+1;
			int dim1 = p_amber_model->tors_tors_params[i][0].size();
			int dim2 = p_amber_model->tors_tors_params[i][1].size();
			FC_FUNC_MODULE(amoeba_torsion_torsion_mod,set_torsion_torsion_params)( &ip1, &dim1, &dim2, 
				      p_amber_model->tors_tors_params[i][0].v(), p_amber_model->tors_tors_params[i][1].v(),
				      p_amber_model->tors_tors_params[i][2].v(), p_amber_model->tors_tors_params[i][3].v(),
					  p_amber_model->tors_tors_params[i][4].v(), p_amber_model->tors_tors_params[i][5].v());
		}
		
		iflag = 1;
		if( p_amber_model->n_tors_tors == 0 || p_amber_model->n_tors_tors_params == 0 ) iflag = 0;
		FC_FUNC_MODULE(amoeba_torsion_torsion_mod,set_valid_bit)(&iflag);

		FC_FUNC_MODULE(amoeba_multipoles_mod,set_local_multipoles_list)( p_amber_model->atm_multipoles.v(), &p_amber_model->num_local_multipoles );
		FC_FUNC_MODULE(amoeba_multipoles_mod,set_chiral_frame_list)( p_amber_model->atm_chiral_frames.v(), &p_amber_model->num_chiral_frames );
		FC_FUNC_MODULE(amoeba_multipoles_mod,set_reg_frame_list)( p_amber_model->atm_reg_frames.v(), &p_amber_model->num_reg_frames );

		iflag = 1;
		if( p_amber_model->num_local_multipoles == 0 ) iflag = 0;
		FC_FUNC_MODULE(amoeba_multipoles_mod,set_valid_bit)(&iflag);

		if( p_amber_model->num_adjust_list != p_amber_model->next )
		{
			PrintLog(" Error in MMDriverAmber::SaveModelToFortran() \n");
			PrintLog(" num_adjust_list = %d   !=  excluded atoms list length = %d \n",
				       p_amber_model->num_adjust_list, p_amber_model->next);
		}
		FC_FUNC_MODULE(amoeba_adjust_mod,set_amoeba_adjust_list)( p_amber_model->atm_adjust_list.v(), &p_amber_model->num_adjust_list );
		FC_FUNC_MODULE(amoeba_adjust_mod,set_adjust_weights)( p_amber_model->adjust_vdw_weights.v(), p_amber_model->adjust_mpole_weights.v(), 
											      p_amber_model->adjust_direct_weights.v(), p_amber_model->adjust_polar_weights.v(), 
												  p_amber_model->adjust_mutual_weights.v() );

		iflag = 1;
		if( p_amber_model->num_adjust_list == 0 ) iflag = 0;
		FC_FUNC_MODULE(amoeba_adjust_mod,set_valid_bit)(&iflag);

		int vdw_atom_cnt = p_amber_model->atm_amoeba_vdw_type.size();
		if( p_amber_model->natom != vdw_atom_cnt )
		{
			PrintLog(" Error in MMDriverAmber::SaveModelToFortran() \n");
			PrintLog(" Size of array atm_amoeba_vdw_type = %d  is not equal to natom = %d \n",
				       vdw_atom_cnt,p_amber_model->natom);
	
		}
		FC_FUNC_MODULE(amoeba_vdw_mod,set_vdw_types_list)( p_amber_model->atm_amoeba_vdw_type.v(), p_amber_model->atm_parent_id.v(), 
			                                   p_amber_model->atm_parent_weight.v(), &vdw_atom_cnt);

		FC_FUNC_MODULE(amoeba_vdw_mod,set_vdw_params)( p_amber_model->amoeba_vdw_depths.v(), p_amber_model->amoeba_vdw_rstars.v(), 
			                               &p_amber_model->vdw_buffer_delta, &p_amber_model->vdw_buffer_gamma, &p_amber_model->n_vdw_params);

		iflag = 1;
		if( vdw_atom_cnt == 0 ) iflag = 0;
		FC_FUNC_MODULE(amoeba_vdw_mod,set_valid_bit)(&iflag);

		FC_FUNC_MODULE(amoeba_induced_mod,set_atm_polar)( p_amber_model->atm_polar.v(), p_amber_model->is_polarizable.v(), p_amber_model->atm_screen_polar.v(), 
			            p_amber_model->damp_polar_strength.v(), p_amber_model->damp_polar_sensitivity.v(), p_amber_model->damp_polar_rad.v(), 
						p_amber_model->atm_hpolar.v(), &p_amber_model->natom );
		set_atm_qterm_(p_amber_model->atm_qterm.v(),&p_amber_model->natom);

		iflag = 1;
		if( p_amber_model->natom == 0 ) iflag = 0;
		FC_FUNC_MODULE(amoeba_induced_mod,set_valid_bit)(&iflag);
	}
//	PrintLog(" MMDriverAmber::SaveModelToFortran() pt end \n");
}

void MMDriverAmber::SetPBoxDataToFortran()
{
#if defined(WITH_LIB_PMEMD)
	set_pbc_data_(&is_orthog, &pbc_alpha, &pbc_beta, &pbc_gamma, 
                   recip.v(), ucell.v(), cut_factor.v(), reclng.v(),
                   &uc_volume, &uc_sphere);

#endif
}

void MMDriverAmber::GetPBoxDataFromFortran()
{
#if defined(WITH_LIB_PMEMD)
	get_pbc_data_(&is_orthog, &pbc_alpha, &pbc_beta, &pbc_gamma, 
                   recip.v(), ucell.v(), cut_factor.v(), reclng.v(),
                   &uc_volume, &uc_sphere);

#endif
}

int MMDriverAmber::SaveAmberRunFile()
{	
	char buf[256];
	FILE* fp = fopen(amber_run_file.c_str(),"w");
	if(fp == NULL)
	{	
		sprintf(buf,"Can't create file %s ",amber_run_file.c_str());
		ErrorInMod("MMDriverAmber::SaveAmberRunFile()", buf);
	}

	fprintf(fp,"#!/bin/sh \n");

	fprintf(fp,"%s -O -i %s -o %s -p %s -c %s -r %s -ref %s -x %s -v %s -e %s \n",
		sander_exe_fname.c_str(),
		amber_inp_file.c_str(), amber_out_file.c_str(), 
		amber_top_file.c_str(), amber_crd_file.c_str(), 
		amber_rst_file.c_str(), amber_constr_crd_file.c_str(),
		amber_trj_coord_file.c_str(),
		amber_trj_vel_file.c_str(), amber_trj_ene_file.c_str() );

	fclose(fp);
	return TRUE;
}


int MMDriverAmber::SaveAmberInpFile()
{	
	InitCtrlParams();
	char buf[256];
	FILE* fp = fopen(amber_inp_file.c_str(),"w");
	if(fp == NULL)
	{	
		sprintf(buf,"Can't create file %s ",amber_inp_file.c_str());
		ErrorInMod("MMDriverAmber::SaveAmberInpFile()", buf);
		return FALSE;
	} 
	
	std::string title;
	title = " Amber Input file for the molecule set ";
	title += pmset->GetName();
	fprintf(fp,"%s\n",title.c_str());   // Write title string	
	fprintf(fp, "\n");                  // Empty String
	fprintf(fp, " &cntrl   \n");
    fprintf(fp, " imin=%d, ",  p_mm_mod->run_type.value() );
	if( p_mm_mod->run_type == p_mm_mod->run_type.MIN_RUN )
	{
		fprintf(fp, " ntmin=%d, ", p_mm_mod->min_type.value() );
		fprintf(fp, " maxcyc=%d,",  p_mm_mod->max_num_minim_steps); 
		fprintf(fp, " ncyc=%d, ",   p_mm_mod->num_steep_descent_steps );
		fprintf(fp, " dx0=%9.3f, ",  p_mm_mod->init_min_step );
		fprintf(fp, " drms=%12.6f, ",  p_mm_mod->grad_cnvrg_val );
	}
	else if( p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN )
	{
        fprintf(fp, " nstlim=%d, ", p_mm_mod->num_md_steps );
		if(amber_version < 7)
		{
		    fprintf(fp, " ntcm=%d, ",   p_mm_mod->remove_init_rb_motion_flag );
		}
        fprintf(fp, " nscm=%d, ", p_mm_mod->remove_rb_motion_freq);
		fprintf(fp, " \n");
        fprintf(fp, " t=%12.3f, ", p_mm_mod->start_time );
        fprintf(fp, " dt=%12.3f, ", p_mm_mod->md_time_step );
	}
	fprintf(fp," \n");
// Print out control:
	fprintf(fp," ntpr=%d, ",   p_mm_mod->wrt_log_freq);
	fprintf(fp," ntwr=%d, ",   p_mm_mod->wrt_rstrt_freq);
	fprintf(fp," ntwx=%d, ",   p_mm_mod->wrt_coord_freq);
	fprintf(fp," ntwv=%d, ",   p_mm_mod->wrt_vel_freq);
	fprintf(fp," ntwe=%d, ",   p_mm_mod->wrt_ener_freq);
	fprintf(fp," ioutfm=%d, ", p_mm_mod->traj_wrt_format.value());
	fprintf(fp," ntwprt=%d, ", p_mm_mod->limit_wrt_atoms);

	fprintf(fp," \n"); 

    fprintf(fp," ntf=%d, ",   p_mm_model->omit_interactions.value());
	fprintf(fp," ntb=%d, ",   p_mm_mod->period_bcond.value());

	if( p_mm_model->electr_method == p_mm_model->electr_method.GEN_BORN 
		|| p_mm_model->electr_method == p_mm_model->electr_method.SCREENED_COULOMB )
	{
		fprintf(fp," saltcon=%6.3f, ", p_mm_model->ion_strength );
	}

	fprintf(fp," dielc=%6.3f, ", p_mm_model->diel_const );
	if( p_mm_model->IsAmoebaFF() ) fprintf(fp," iamoeba=%d, ", 1 );

	fprintf(fp," \n");
	fprintf(fp," cut=%6.3f, ",  p_mm_model->nb_cut_dist);
	if( p_mm_mod->ext_mm_prog != MMExternalProg::PMEMD_12 && p_mm_mod->ext_mm_prog != MMExternalProg::SANDER_12 && p_mm_mod->ext_mm_prog != MMExternalProg::PMEMD_18 )
	{
		fprintf(fp," scnb=%6.3f, ", p_mm_model->scale_14_vdw);
		fprintf(fp," scee=%6.3f, ", p_mm_model->scale_14_electr);
	}

	fprintf(fp," iwrap=%d, ",  p_mm_mod->wrap_coord);	

	fprintf(fp," \n");
	fprintf(fp," irest= %d, ", p_mm_mod->restart_flag);
	fprintf(fp," ntx= %d, ",   p_mm_mod->init_read_coord.value() );
	if (p_mm_mod->ext_mm_prog == MMExternalProg::PMEMD_18)
	{
		fprintf(fp, " ntxo= %d, mdinfo_flush_interval = %d, ", 1, 2);
	}
	fprintf(fp,"\n"); 

	fprintf(fp,"igb=%d, \n", p_amber_model->igb);

// Frozen and restrained atoms:

	fprintf(fp,"ibelly=%d, ",p_amber_model->ibelly);
	int ntr_loc = 0; if( p_amber_model->natc > 0) ntr_loc = 1;
	fprintf(fp,"ntr=%d, ", ntr_loc );
	fprintf(fp,"\n"); 

// temperature regulation:

	fprintf(fp," temp0=%6.3f, ", p_mm_mod->ref_temp);
	fprintf(fp," tempi=%6.3f, ", p_mm_mod->init_temp);
	fprintf(fp," ntt=%d, ",      p_mm_mod->temp_control_method.value() );
	if( p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN )
	{
		fprintf(fp," gamma_ln=%6.3f, ", p_mm_mod->langevin_dump_const );
	}
	fprintf(fp,"\n"); 

//	fprintf(fp," isolvp=%d, ", last_solute_atom );
	if( p_mm_mod->ext_mm_prog != MMExternalProg::PMEMD_12 && p_mm_mod->ext_mm_prog != MMExternalProg::SANDER_12 && p_mm_mod->ext_mm_prog != MMExternalProg::PMEMD_18)
	{
		fprintf(fp," dtemp=%9.4f, ", p_mm_mod->temp_deviation);
	}
	fprintf(fp," tautp=%9.4f, ", p_mm_mod->temp_relax_time_solute);
	fprintf(fp,"\n");

// pressure regulation

    fprintf(fp," ntp=%d, ", p_mm_mod->pressure_reg_method.value() );
    fprintf(fp," pres0=%6.3f, ", p_mm_mod->ref_pressure);
	fprintf(fp," comp=%6.3f, ", p_mm_mod->compressibility);
	fprintf(fp," taup=%6.3f, ",p_mm_mod->press_relax_time);
	fprintf(fp,"\n");
	
// SHAKE bond length constraints:

	fprintf(fp," ntc=%d,",      p_mm_mod->shake_constr.value() );
	fprintf(fp," tol=%12.7f, ", p_mm_mod->shake_tol);
	fprintf(fp,"\n");

// Special Water treatment:

	fprintf(fp," jfastw=%d, ", p_mm_mod->fast_water_method);
	fprintf(fp," \n ");

// Water Cap treatment:

//	fprintf(fp," ivcap=%d, ", water_cap_flag);
//	fprintf(fp," matcap=%d, ", cap_atom_num);
//	fprintf(fp," fcap=%6.3f, ", cap_fconst);
//	fprintf(fp," \n ");
	
	fprintf(fp, " &end   \n");

	if(p_mm_mod->period_bcond > 0 && p_mm_model->IsAmoebaFF() )
	{
	   fprintf(fp, " &ewald \n");
	   fprintf(fp, " skinnb=%12.6f, nbtell=0, order=%d, ew_coeff=%12.6f \n",p_mm_model->skin_nb, 
		            p_mm_model->pme_spline_order, p_mm_model->pme_ew_coeff);
//	   fprintf(fp, " a=%12.6f, b=%12.6f, c= %12.6f, \n alpha=%12.6f, beta=%12.6f, gamma=%12.6f \n",
//	   		    pmset->per_bc->GetA(), pmset->per_bc->GetB(),pmset->per_bc->GetC(), 
//			    pmset->per_bc->GetAlpha()*RAD_TO_DEG, pmset->per_bc->GetBeta()*RAD_TO_DEG, pmset->per_bc->GetGamma()*RAD_TO_DEG);
	   fprintf(fp, " &end \n");
	}

	AtomIntMap& pt_idx_map = p_mm_model->GetAtIdxMap(TRUE);
		
	if( p_mm_model->DistConstraints.size() > 0)
	{	
		vector<AtomContact>::iterator ptr_itr; 
		for(ptr_itr = p_mm_model->DistConstraints.begin();ptr_itr != p_mm_model->DistConstraints.end();ptr_itr++)
		{
		    fprintf(fp, " &rst ");

			int indx1 = pt_idx_map[ptr_itr->pt1] + 1;
			int indx2 = pt_idx_map[ptr_itr->pt2] + 1;
			fprintf(fp,"iat = %d, %d, ", indx1,indx2);

			double dist   = ptr_itr->GetRMin();
			double frcpar = ptr_itr->GetHarmForceConst();

			fprintf(fp,"r2 = %9.5f, r3 = %9.5f, ",dist, dist);
			fprintf(fp,"r1 = 0.0, r4 = 200.0, ");//mostly parabolic potential
			fprintf(fp,"rk2 = %9.5f, rk3 = %9.5f ",frcpar, frcpar);
		    fprintf(fp, " &end   \n");
        }  
	}
			
	AtomGroup* mv_atlist =  p_mm_model->GetMovingAtoms();
	if(mv_atlist)
	{
		fprintf(fp,"%.80s\n",p_mm_model->moving_atoms.c_str());
		AtomIteratorAtomGroup aitr(mv_atlist);
		HaAtom* pt;
		int ibeg = -1; 
		int ifin = -1;
		for(pt = aitr.GetFirstAtom(); pt; pt = aitr.GetNextAtom())
		{
			int idx = pt_idx_map[pt] + 1;
			if( ibeg < 0) 
			{
				ibeg = idx;
				ifin = idx;
				continue;
			}
			if( idx == (ifin + 1))
			{
				ifin++;
				continue;
			}
			else
			{
				fprintf(fp,"ATOM  %6d %6d \n",ibeg,ifin);	
				ibeg = idx;
				ifin = idx;
			}
		}
		if( ibeg > 0) 
		{
			fprintf(fp,"ATOM  %6d %6d \n",ibeg,ifin);	
		}
		fprintf(fp,"END         \n");
		fprintf(fp,"END         \n");
	}

	AtomGroup* restr_atlist =  p_mm_model->GetRestrAtoms();
	if(restr_atlist)
	{
		std::string pos_restr_str = "Positional restraints group: " + p_mm_model->restrained_atoms;
		fprintf(fp,"%.80s\n",pos_restr_str.c_str());
		fprintf(fp,"%12.6f \n",p_mm_model->atom_restr_const);
		AtomIteratorAtomGroup aitr(restr_atlist);
		HaAtom* pt;
		int ibeg = -1; 
		int ifin = -1;
		for(pt = aitr.GetFirstAtom(); pt; pt = aitr.GetNextAtom())
		{
			int idx = pt_idx_map[pt] + 1;
			if( ibeg < 0) 
			{
				ibeg = idx;
				ifin = idx;
				continue;
			}
			if( idx == (ifin + 1))
			{
				ifin++;
				continue;
			}
			else
			{
				fprintf(fp,"ATOM  %6d %6d \n",ibeg,ifin);	
				ibeg = idx;
				ifin = idx;
			}
		}
		if( ibeg > 0) 
		{
			fprintf(fp,"ATOM  %6d %6d \n",ibeg,ifin);	
		}
		fprintf(fp,"END         \n");
		fprintf(fp,"END         \n");
	}

	fclose(fp);

	return TRUE;
}

int MMDriverAmber::SaveAmberCrdFile()
{
	int ires = SaveAmberRstFile( amber_crd_file.c_str());
	return ires;
}

int MMDriverAmber::SaveAmberTopFile()
{
	if(p_mm_model->to_init_mm_model) p_mm_mod->InitMolMechModel();
	InitCtrlParams();
	ofstream os(amber_top_file.c_str());
	if(os.fail()) 
	{
		PrintLog("Error in MMDriverAmber::SaveAmberTopFile() \n");
		PrintLog("Can't create file %s \n", amber_top_file.c_str());
		return FALSE;
	}
	return SaveAmberTopToStream(os);
}


int MMDriverAmber::SaveAmberTopToStream(ostream& os)
{
	if(os.fail()) return FALSE;

	if( p_mm_model->ff_type == ForceFieldType::AMOEBA && 
	/*    p_mm_mod->ext_mm_prog != MMExternalProg::PMEMD_10 &&  p_mm_mod->ext_mm_prog != MMExternalProg::SANDER_10 && */
	    p_mm_mod->ext_mm_prog != MMExternalProg::PMEMD_12 &&  p_mm_mod->ext_mm_prog != MMExternalProg::SANDER_12 && p_mm_mod->ext_mm_prog != MMExternalProg::PMEMD_18 )
	{
		PrintLog(" HARLEM doesn't support setup of AMOEBA force field for external program %s\n",p_mm_model->ff_type.label());
		return FALSE;
	}
	char buf[200];

	int nhparm = 0;                     // not used
	int nparm  = 0;                     // not used        
	int mbona  = p_amber_model->nbona;    // NBONA + number of constraint bonds     (not used) 
	int mtheta = p_amber_model->ntheta;   // MTHETA + number of constraint angles   (not used) 
	int mphia  = p_amber_model->nphia;    // MPHIA + number of constraint dihedrals (not used) 
    int natyp  = p_amber_model->ntypes;   // number of atom types in parameter file, see SOLTY below         
    int nfpert = 0;                 // set to 1 if perturbation info is to be read in
    int nbper  = 0;                 // number of bonds to be perturbed         (not used)
	int ngper  = 0;                 // number of angles to be perturbed        (not used)
	int ndper  = 0;                 // number of dihedrals to be perturbed     (not used)
	int mbper  = 0;                 // number of bonds with atoms completely in perturbed group      (not used)
	int mgper  = 0;                 // number of angles with atoms completely in perturbed group     (not used)
	int mdper  = 0;                 // number of dihedrals with atoms completely in perturbed groups (not used)
	int nmxrs  = p_amber_model->max_res_size;    // number of atoms in the largest residue  (not used ?)

// Write title string

	os << "%VERSION VERSION_STAMP = V0001.000 ";
    os << "DATE = 05/22/06  12:10:21" << std::endl;
        
    os << "%FLAG TITLE" << std::endl;
    os << "%FORMAT(20a4)" << std::endl;
    os << title.c_str() << std::endl;

    os << "%FLAG POINTERS" << std::endl;
    os << "%FORMAT(10I8)"  << std::endl;

	sprintf(buf,"%8d",p_amber_model->natom); 
	os << buf;  
	if( p_mm_model->ff_type == ForceFieldType::AMOEBA /* && 
		(p_mm_mod->ext_mm_prog == MMExternalProg::PMEMD_10 || p_mm_mod->ext_mm_prog != MMExternalProg::SANDER_10) */)
	{
		sprintf(buf,"%8d%8d%8d%8d%8d%8d%8d%8d%8d",1,1,1,1,1,1,1,0,0); 
		os << buf;
	}
	else
	{	
		sprintf(buf,"%8d%8d%8d%8d%8d%8d%8d%8d%8d",p_amber_model->ntypes,p_amber_model->nbonh,mbona,
				p_amber_model->ntheth,mtheta,p_amber_model->nphih,mphia,nhparm,nparm);  
		os << buf;

	}
    os << std::endl; 

	sprintf(buf,"%8d%8d",p_amber_model->next,p_amber_model->nres);
	os << buf;
	if( p_mm_model->ff_type == ForceFieldType::AMOEBA /* &&  
		(p_mm_mod->ext_mm_prog == MMExternalProg::PMEMD_10 || p_mm_mod->ext_mm_prog != MMExternalProg::SANDER_10) */ )
	{
		sprintf(buf,"%8d%8d%8d%8d%8d%8d%8d%8d",1,1,1,1,1,1,1,1);
		os << buf;
	}
	else
	{	
		sprintf(buf,"%8d%8d%8d%8d%8d%8d%8d%8d",p_amber_model->nbona,p_amber_model->ntheta,p_amber_model->nphia,
		       p_amber_model->numbnd,p_amber_model->numang,p_amber_model->nptra,natyp,p_amber_model->nphb);
		os << buf;
	}
	os << std::endl; 
	
	if( p_mm_model->ff_type == ForceFieldType::AMOEBA /* &&  
	(p_mm_mod->ext_mm_prog == MMExternalProg::PMEMD_10 || p_mm_mod->ext_mm_prog != MMExternalProg::SANDER_10) */ )
	{
		sprintf(buf,"%8d%8d%8d%8d%8d%8d%8d",0,0,0,0,0,0,0);
		os << buf;
	}
	else
	{
		sprintf(buf,"%8d%8d%8d%8d%8d%8d%8d",nfpert,nbper,ngper,ndper,mbper,mgper,mdper);
		os << buf;
	}
	int ifcap_dummy = 0;
	int ifbox_dummy = 0;
	if( pmset->per_bc->IsSet() )
	{
		ifbox_dummy = 3;
		if( pmset->per_bc->IsOrthogonal() ) ifbox_dummy = 1;
		if( pmset->per_bc->IsOctahedron() ) ifbox_dummy = 2;
	}
	sprintf(buf,"%8d%8d%8d",ifbox_dummy,nmxrs,ifcap_dummy);
    os << buf << std::endl;

	int numextra_dummy = 0;
	int ncopy_dummy = 0;
	sprintf(buf,"%8d%8d%8d", numextra_dummy, ncopy_dummy, p_amber_model->num_dist_constr);
	os << buf << std::endl;

// Write Atom Names Array 

	os << "%FLAG ATOM_NAME                                                                 " << std::endl;
	os << "%FORMAT(20a4)                                                                   " << std::endl;	
	write_string_array_chuncks(os, p_amber_model->atm_igraph, 4, 20);
		
	int i;

// Save Atom Charges:

	double dlc = 1.0; 
	if( p_amber_model->using_pme_potential && p_mm_model->diel_const > 1.0 )  dlc = p_mm_model->diel_const;
	HaVec_double farr(p_amber_model->natom);
	for( i = 0; i < p_amber_model->natom; i++)
	{
		farr[i] = p_amber_model->atm_charge[i]*dlc; //!< scale atom charges back by diel_const
	}
	os << "%FLAG CHARGE                                                                    " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,farr,5,FLOAT_E16_8); 
	
// Save Atom Masses:

	os << "%FLAG MASS                                                                      " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->atm_mass,5,FLOAT_E16_8); 

// IAC: indexes of atoms in the array of atom types (Fortran based 1)
	
	os << "%FLAG ATOM_TYPE_INDEX                                                           " << std::endl;
	os << "%FORMAT(10I8)                                                                   " << std::endl;
	write_int_array_chuncks(os,p_amber_model->atm_iac,10,"%8d");                     

// Excluded Atoms Arrays

	os << "%FLAG NUMBER_EXCLUDED_ATOMS                                                     " << std::endl;
	os << "%FORMAT(10I8)                                                                   " << std::endl;
	write_int_array_chuncks(os,p_amber_model->atm_numex,10,"%8d"); // NUMEX : total number of excluded atoms for atom "i"
	                                                   // NATEX below

// Non-bonded iteractions parameters index 

// Correction for H-bonds:
//   int idx_hbond = 0;
//   if( iprm.nphb > 0 )
//   {
//		for( i = idx_fst_hb_atom_type; i < (ntypes-1); i++)
//		{
//			idx = ntypes*i + i+1;
//			idx_hbond++;
//			iarr[idx] = -idx_hbond;
//			idx = ntypes*(i+1) + i;
//          iarr[idx] = -idx_hbond;
//			i++;
//		}
//   }

	os << "%FLAG NONBONDED_PARM_INDEX                                                      " << std::endl;
	os << "%FORMAT(10I8)                                                                   " << std::endl;
	write_int_array_chuncks(os,p_amber_model->typ_ico,10,"%8d"); // ICO: indexes of pairs of atom types in the arrays
	                                                  // of nonbonded parameters	
// Residue Labels:

	os << "%FLAG RESIDUE_LABEL                                                             " << std::endl;
	os << "%FORMAT(20a4)                                                                   " << std::endl;
	write_string_array_chuncks(os,p_amber_model->res_labels, 4, 20); // FORMAT(20A4) (LABRES(i), i=1,NRES)
                                                                     // LABRES : the residue labels
// Indexes of first atoms of residues (RESIDUE_POINTER)

	os << "%FLAG RESIDUE_POINTER                                                           " << std::endl;
	os << "%FORMAT(10I8)     " << std::endl;
	HaVec_int iaxx;
	if (p_amber_model->gbl_res_atms.size() > 0)
	{
		iaxx.resize(p_amber_model->gbl_res_atms.size() - 1);
		for (int i = 0; i < p_amber_model->gbl_res_atms.size() - 1; i++)
			iaxx[i] = p_amber_model->gbl_res_atms[i];
	}
	write_int_array_chuncks(os,iaxx,10,"%8d");   // FORMAT(10I8) (IPRES(i), i=1,NRES) 
                                                             // IPRES : the atom number of the first atom in residue
// Bonds constants:
	
	os << "%FLAG BOND_FORCE_CONSTANT                                                       " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->gbl_rk,5,FLOAT_E16_8);  // FORMAT(5E16.8) (RK(i), i=1,NUMBND)
	                                              // RK : force constant for the bonds of each type, kcal/mol
// Bond Equilibrium distances:

	os << "%FLAG BOND_EQUIL_VALUE                                                          " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->gbl_req,5,FLOAT_E16_8); // FORMAT(5E16.8) (REQ(i), i=1,NUMBND)
                                                       // REQ : equilibrium bond length for the bonds of each type, angstroms
// Valence Angle Force Constants

	os << "%FLAG ANGLE_FORCE_CONSTANT                                                      " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->gbl_tk,5,FLOAT_E16_8); // FORMAT(5E16.8) (TK(i), i=1,NUMANG)
			                                                   // TK : force constant for the angles of each type, kcal/mol A**2

// Valence Angle Equilibrium Values
	
	os << "%FLAG ANGLE_EQUIL_VALUE                                                         " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->gbl_teq,5,FLOAT_E16_8); // FORMAT(5E16.8) (TEQ(i), i=1,NUMANG)
											        	  // TEQ : the equilibrium angle for the angles of each type, degrees
// Dihedral Angle Force Constants
	os << "%FLAG DIHEDRAL_FORCE_CONSTANT                                                   " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->gbl_pk,5,FLOAT_E16_8); // FORMAT(5E16.8) (PK(i), i=1,NPTRA)
	        		                                     // PK : force constant for the dihedrals of each type, kcal/mol

// Dihedral Periodicity
	
	os << "%FLAG DIHEDRAL_PERIODICITY                                                      " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->gbl_pn,5,FLOAT_E16_8); // FORMAT(5E16.8) (PN(i), i=1,NPTRA)
			                                         // PN : periodicity of the dihedral of a given type
// Dihedral phase

	os << "%FLAG DIHEDRAL_PHASE                                                            " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->gbl_phase,5,FLOAT_E16_8); // FORMAT(5E16.8) (PHASE(i), i=1,NPTRA)
			                                                 // PHASE : phase of the dihedral of a given type
	
	if( p_mm_mod->ext_mm_prog == MMExternalProg::PMEMD_12 || p_mm_mod->ext_mm_prog == MMExternalProg::SANDER_12 )
	{
// 1-4 Elecstrostic interactions scale factors (SCEE)
	os << "%FLAG SCEE_SCALE_FACTOR                                                         " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	farr.resize( p_amber_model->gbl_pk.size() );
	farr = p_mm_model->scale_14_electr;
	write_double_array_chuncks(os,farr,5,FLOAT_E16_8); // FORMAT(5E16.8) (gbl_one_scee(i), i=1,NPTRA)
			                                           // gbl_one_scee : Electrostatic 1-4 scale factors for dihedral angles
	
// 1-4 non bonded interactions scale dactors (SCNB) 
	os << "%FLAG SCNB_SCALE_FACTOR                                                         " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	farr.resize( p_amber_model->gbl_pk.size() );
	farr = p_mm_model->scale_14_vdw;
	write_double_array_chuncks(os,farr,5,FLOAT_E16_8); // FORMAT(5E16.8) (gbl_one_scnb(i), i=1,NPTRA)
			                                           // gbl_one_scnb : Nonbonded 1-4 interactions scale factors for dihedral angles
	}

// SOLTY ARRAY - Currently Not used 

	farr.resize(natyp);
	for( i = 0; i < natyp; i++)
		farr[i] = 0.0;

	os << "%FLAG SOLTY                                                                     " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,farr,5,FLOAT_E16_8); // FORMAT(5E16.8) (SOLTY(i), i=1,NATYP)
                                                       // SOLTY : currently unused (reserved for future use)
// Atom Van der Waals Alpha (CN1) parameters:

	os << "%FLAG LENNARD_JONES_ACOEF                                                       " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->gbl_cn1,5,FLOAT_E16_8); // FORMAT(5E16.8) (CN1(i), i=1,NTYPES*(NTYPES+1)/2)
			                                         // CN1 : Lennard Jones r**12 terms for all possible atom type
				                                     // interactions, indexed by ICO and IAC; for atom i and j
					                                 // where i < j, the index into this array is as follows
						                             // (assuming the value of ICO(INDEX) is positive):
							                         // CN1(ICO(NTYPES*(IAC(i)-1)+IAC(j))).
	
// Atom Van der Waals Beta (CN2) parameters:

	os << "%FLAG LENNARD_JONES_BCOEF                                                       " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->gbl_cn2,5,FLOAT_E16_8); // FORMAT(5E16.8) (CN2(i), i=1,NTYPES*(NTYPES+1)/2)
			                                         // CN2 : Lennard Jones r**6 terms for all possible atom type
				                                     // interactions. Indexed like CN1 above.
	HaVec_int iarr;

// NOTE: the atom numbers in the arrays which follow that describe bonds,
// angles, and dihedrals are obfuscated by the following formula (for
// runtime speed in indexing arrays). The true atom number equals the
// absolute value of the number divided by three, plus one. In the case
// of the dihedrals, if the third atom is negative, this implies an
// improper torsion and if the fourth atom is negative, this implies that
// end group interactions are to be ignored. End group interactions are
// ignored, for example, in dihedrals of various ring systems (to prevent
// double counting) and in multiterm dihedrals.

	int i1, i2, i3, i4;
    
	HaVec_double bpar(2);
	HaVec_double vang_par(2);

// Bonds including hydrogen 

	iarr.resize(3*p_amber_model->nbonh);
	for( i=0; i < p_amber_model->nbonh; i++)
	{
		iarr[3*i]    =  (p_amber_model->gbl_bond[3*i]   - 1) * 3;
		iarr[3*i + 1] = (p_amber_model->gbl_bond[3*i+1] - 1) * 3;
		iarr[3*i + 2] =  p_amber_model->gbl_bond[3*i+2];
	}
	os << "%FLAG BONDS_INC_HYDROGEN                                                        " << std::endl;
	os << "%FORMAT(10I8)                                                                   " << std::endl;
	write_int_array_chuncks(os,iarr,10,"%8d");   // FORMAT(10I8) (IBH(i),JBH(i),ICBH(i), i=1,NBONH)
			                                             // IBH : atom involved in bond "i", bond contains hydrogen
				                                         // JBH : atom involved in bond "i", bond contains hydrogen
					                                     // ICBH : index into parameter arrays RK and REQ

// Bonds Not Including Hydrogen

	iarr.resize(3*p_amber_model->nbona);
	for( i=0; i < p_amber_model->nbona; i++)
	{
		int idx = p_amber_model->nbonh + i;
		iarr[3*i]    =  (p_amber_model->gbl_bond[3*idx]   - 1) * 3;
		iarr[3*i + 1] = (p_amber_model->gbl_bond[3*idx+1] - 1) * 3;
		iarr[3*i + 2] =  p_amber_model->gbl_bond[3*idx+2];
	}
	os << "%FLAG BONDS_WITHOUT_HYDROGEN                                                    " << std::endl;
	os << "%FORMAT(10I8)                                                                   " << std::endl;
	write_int_array_chuncks(os,iarr,10,"%8d"); // FORMAT(10I8) (IB(i),JB(i),ICB(i), i=1,NBONA)
			                                       // IB : atom involved in bond "i", bond does not contain hydrogen
				                                   // JB : atom involved in bond "i", bond does not contain hydrogen
					                               // ICB : index into parameter arrays RK and REQ
			
// Valence Angles contaning hydrogen

	iarr.resize(4*p_amber_model->ntheth);
	for( i=0; i < p_amber_model->ntheth; i++)
	{
		int idx = i;
		iarr[4*i]    =  (p_amber_model->gbl_angle[4*idx]   - 1) * 3;
		iarr[4*i + 1] = (p_amber_model->gbl_angle[4*idx+1] - 1) * 3;
		iarr[4*i + 2] = (p_amber_model->gbl_angle[4*idx+2] - 1) * 3;
		iarr[4*i + 3] =  p_amber_model->gbl_angle[4*idx+3];
	}
	os << "%FLAG ANGLES_INC_HYDROGEN                                                       " << std::endl;
	os << "%FORMAT(10I8)                                                                   " << std::endl;
	write_int_array_chuncks(os,iarr,10,"%8d"); // FORMAT(10I8) (ITH(i),JTH(i),KTH(i),ICTH(i), i=1,NTHETH)
			                                       // ITH : atom involved in angle "i", angle contains hydrogen
				                                   // JTH : atom involved in angle "i", angle contains hydrogen
					                               // KTH : atom involved in angle "i", angle contains hydrogen
						                           // ICTH : index into parameter arrays TK and TEQ for angle
							                       //

// Valence Angles not containing hydrogen

	iarr.resize(4*p_amber_model->ntheta);
	for( i=0; i < p_amber_model->ntheta; i++)
	{
		int idx = p_amber_model->ntheth + i;
		iarr[4*i]    =  (p_amber_model->gbl_angle[4*idx]   - 1) * 3;
		iarr[4*i + 1] = (p_amber_model->gbl_angle[4*idx+1] - 1) * 3;
		iarr[4*i + 2] = (p_amber_model->gbl_angle[4*idx+2] - 1) * 3;
		iarr[4*i + 3] =  p_amber_model->gbl_angle[4*idx+3];
	}
	os << "%FLAG ANGLES_WITHOUT_HYDROGEN                                                   " << std::endl;
	os << "%FORMAT(10I8)                                                                   " << std::endl;
	write_int_array_chuncks(os,iarr,10,"%8d"); // FORMAT(10I8) (IT(i),JT(i),KT(i),ICT(i), i=1,NTHETA)

			                                        // IT : atom involved in angle "i", angle does not contain hydrogen
				                                    // JT : atom involved in angle "i", angle does not contain hydrogen
					                                // KT : atom involved in angle "i", angle does not contain hydrogen
						                            // ICT : index into parameter arrays TK and TEQ for angle
							                        // IT(i)-JT(i)-KT(i)

// Dihedrals including hydrogens

	iarr.resize(5*p_amber_model->nphih);
	for( i=0; i < p_amber_model->nphih; i++)
	{
		int idx = i;
		iarr[5*i]    =  (p_amber_model->gbl_dihed[5*idx  ] - 1) * 3;
		iarr[5*i + 1] = (p_amber_model->gbl_dihed[5*idx+1] - 1) * 3;
		iarr[5*i + 2] = (p_amber_model->gbl_dihed[5*idx+2] - 1) * 3;
		if( p_amber_model->gbl_dihed[5*idx+2] < 0) iarr[5*i + 2] = (p_amber_model->gbl_dihed[5*idx+2] + 1) * 3; // not calculate 1-4 interactions 
		iarr[5*i + 3] = (p_amber_model->gbl_dihed[5*idx+3] - 1) * 3;
		if( p_amber_model->gbl_dihed[5*idx+3] < 0) iarr[5*i + 3] = (p_amber_model->gbl_dihed[5*idx+3] + 1) * 3; // improper dihedrals
		iarr[5*i + 4] =  p_amber_model->gbl_dihed[5*idx+4];
	}
	os << "%FLAG DIHEDRALS_INC_HYDROGEN                                                    " << std::endl;
	os << "%FORMAT(10I8)                                                                   " << std::endl;
	write_int_array_chuncks(os,iarr,10,"%8d"); // FORMAT(10I8) (IPH(i),JPH(i),KPH(i),LPH(i),ICPH(i), i=1,NPHIH)
			
// Dihedrals without hydrogens

	iarr.resize(5*p_amber_model->nphia);
	for( i=0; i < p_amber_model->nphia; i++)
	{
		int idx = p_amber_model->nphih + i;
		iarr[5*i]    =  (p_amber_model->gbl_dihed[5*idx  ] - 1) * 3;
		iarr[5*i + 1] = (p_amber_model->gbl_dihed[5*idx+1] - 1) * 3;
		iarr[5*i + 2] = (p_amber_model->gbl_dihed[5*idx+2] - 1) * 3;
		if( p_amber_model->gbl_dihed[5*idx+2] < 0) iarr[5*i + 2] = (p_amber_model->gbl_dihed[5*idx+2] + 1) * 3; // improper
		iarr[5*i + 3] = (p_amber_model->gbl_dihed[5*idx+3] - 1) * 3;
		if( p_amber_model->gbl_dihed[5*idx+3] < 0) iarr[5*i + 3] = (p_amber_model->gbl_dihed[5*idx+3] + 1) * 3; // not calculate 1-4 interactions ?? (or only calculate 1-4 interactions?)
		iarr[5*i + 4] =  p_amber_model->gbl_dihed[5*idx+4];
	}
	os << "%FLAG DIHEDRALS_WITHOUT_HYDROGEN                                                " << std::endl;
	os << "%FORMAT(10I8)                                                                   " << std::endl;
	write_int_array_chuncks(os,iarr,10,"%8d"); // FORMAT(10I8) (IP(i),JP(i),KP(i),LP(i),ICP(i), i=1,NPHIA)

// Nonbonded interations atom excluded list

	os << "%FLAG EXCLUDED_ATOMS_LIST                                                       " << std::endl;
	os << "%FORMAT(10I8)                                                                   " << std::endl;
	write_int_array_chuncks(os,p_amber_model->gbl_natex,10,"%8d"); // FORMAT(10I8) (NATEX(i), i=1,NEXT)
	
		                                           // NATEX : the excluded atom list. To get the excluded list for atom 
			                                       // "i" you need to traverse the NUMEX list, adding up all
				                                   // the previous NUMEX values, since NUMEX(i) holds the number
					                               // of excluded atoms for atom "i", not the index into the 
						                           // NATEX list. Let IEXCL = SUM(NUMEX(j), j=1,i-1), then
							                       // excluded atoms are NATEX(IEXCL) to NATEX(IEXCL+NUMEX(i)).

// H Bond ALPHA(CN1) parameters

	os << "%FLAG HBOND_ACOEF                                                               " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->gbl_asol,5,FLOAT_E16_8); // FORMAT(5E16.8) (ASOL(i), i=1,NPHB)
			                                         // ASOL : the value for the r**12 term for hydrogen bonds of all
				                                     // possible types. Index into these arrays is equivalent
					                                 // to the CN1 and CN2 arrays, however the index is negative.
						                             // For example, for atoms i and j, with i < j, the index is
							                         // -(NTYPES*(IAC(i)-1)+IAC(j)).
	
// H Bond BETA(CN2) parameters

	os << "%FLAG HBOND_BCOEF                                                               " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->gbl_bsol,5,FLOAT_E16_8); // FORMAT(5E16.8) (BSOL(i), i=1,NPHB)
                                                          // BSOL : the value for the r**10 term for hydrogen bonds of all
                                                          // possible types. Indexed like ASOL.
	
// H Bond Cutoff parameter (not used)

	farr.resize(p_amber_model->nphb);
	for(i = 0; i < p_amber_model->nphb; i++)
	{
		farr[i] = 0.0;
	}
		
	os << "%FLAG HBCUT                                                                     " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;	
	write_double_array_chuncks(os,farr,5,FLOAT_E16_8); // FORMAT(5E16.8) (HBCUT(i), i=1,NPHB)
			                                           // HBCUT : no longer in use

//  Non-bonded constrains info
   
    if( p_amber_model->num_dist_constr > 0 )               
	{                                              
		os << "%FLAG NON_BONDED_CONSTRAINTS                                                     " << std::endl;
		os << "%FORMAT(10I8)                                                                    " << std::endl;
		write_int_array_chuncks(os,p_amber_model->dist_constr_idx,10,"%8d");       
		
		os << "%FLAG NON_BONDED_CONSTRAINTS_COEF                                                " << std::endl;
		os << "%FORMAT(5E16.8)                                                                  " << std::endl;
		write_double_array_chuncks(os,p_amber_model->dist_constr_params,5,FLOAT_E16_8); 
	}

// AMBER ATOM FORCE FIELD SYMBOLS

	os << "%FLAG AMBER_ATOM_TYPE                                                           " << std::endl;
	os << "%FORMAT(20a4)                                                                   " << std::endl;
	write_string_array_chuncks(os, p_amber_model->atm_isymbl, 4, 20); // FORMAT(20A4) (ISYMBL(i), i=1,NATOM)
			                                            // ISYMBL : the AMBER atom types for each atom 
	
// ATOM TREE SYMBOLS
	
	os << "%FLAG TREE_CHAIN_CLASSIFICATION                                                 " << std::endl;
	os << "%FORMAT(20a4)                                                                   " << std::endl;
	write_string_array_chuncks(os, p_amber_model->atm_itree, 4, 20); // FORMAT(20A4) (ITREE(i), i=1,NATOM)
			                                          // ITREE : the list of tree joining information, classified into five
			                                          // types. M -- main chain, S -- side chain, B -- branch point, 
			  	                                      // 3 -- branch into three chains, E -- end of the chain
	
// TREE JOIN INFO ARRAYS  (Not USED)
		
	iarr.newsize(p_amber_model->natom);
	for( i=0; i < p_amber_model->natom; i++)
	{
		iarr[i] = 0;
	}

	os << "%FLAG JOIN_ARRAY                                                                " << std::endl;
	os << "%FORMAT(10I8)                                                                   " << std::endl;
	write_int_array_chuncks(os,iarr,10,"%8d"); // FORMAT(10I8) (JOIN(i), i=1,NATOM)
			                                       // JOIN : tree joining information, potentially used in ancient
				                                   // analysis programs. Currently unused in sander or gibbs

// IROTAT ARRAY (Not Used)

	os << "%FLAG IROTAT                                                                    " << std::endl;
	os << "%FORMAT(10I8)                                                                   " << std::endl;
	write_int_array_chuncks(os,iarr,10,"%8d"); // FORMAT(10I8) (IROTAT(i), i = 1, NATOM)
	                                           // IROTAT : apparently the last atom that would move if atom i was
                                               // rotated, however the meaning has been lost over time.
                                               // Currently unused in sander or gibbs.

// Periodic Boundary conditions and info

	if( pmset->per_bc->IsSet() ) // Additional output for periodic boundary conditions and solvent info
	{           
		iarr.newsize(3);
		iarr[0] = p_amber_model->n_solute_res;     // IPTRES
		iarr[1] = p_amber_model->nspm;             // NSPM
		iarr[2] = p_amber_model->n_solute_mol + 1; // NSPSOL - index of the first solvent molecule ??

		os << "%FLAG SOLVENT_POINTERS                                                          " << std::endl;
		os << "%FORMAT(10I8)                                                                   " << std::endl;
		write_int_array_chuncks(os,iarr,10,"%8d"); // FORMAT(10I8) IPTRES,NSPM,NSPSOL
		
		os << "%FLAG ATOMS_PER_MOLECULE                                                        " << std::endl;
		os << "%FORMAT(10I8)                                                                   " << std::endl;
		write_int_array_chuncks(os,p_amber_model->atm_nsp,10,"%8d"); // FORMAT(10I8) (IX(I+I70-1),I=1,NSPM)

		farr.newsize(4);
		farr[0] = pmset->per_bc->GetBeta() * RAD_TO_DEG;
		farr[1] = pmset->per_bc->GetA();
		farr[2] = pmset->per_bc->GetB();
		farr[3] = pmset->per_bc->GetC();

		os << "%FLAG BOX_DIMENSIONS                                                            " << std::endl;
		os << "%FORMAT(5E16.8)                                                                 " << std::endl;
		write_double_array_chuncks(os,farr,4,FLOAT_E16_8);
	}

// CAP information is not set yet

// Generalized Born parameters: 
		
// Atom Born radii 

	os << "%FLAG RADII                                                                     " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->atm_gb_radii,5,FLOAT_E16_8); // save atom radii

// Screening factors:
		
	os << "%FLAG SCREEN                                                                    " << std::endl;
	os << "%FORMAT(5E16.8)                                                                 " << std::endl;
	write_double_array_chuncks(os,p_amber_model->atm_gb_fs,5,FLOAT_E16_8); // save atom screening parameters

	if( p_mm_model->IsAmoebaFF() )
	{
		os << "%FLAG AMOEBA_FORCEFIELD" << std::endl;
        os << "%COMMENT This indicates that this parm file is specific to amoeba" << std::endl;
        os << "%COMMENT This must be present if do_amoeba(in mdin) is 1" << std::endl;
        os << "%COMMENT This must NOT be present if do_amoeba is 0" << std::endl;
        os << "%FORMAT(i5)" << std::endl;
        os << "    1" <<  std::endl;

		os << "%FLAG AMOEBA_ATOM_TYPE_INDEX" << std::endl;
        os << "%COMMENT   dimention = (" << p_amber_model->natom << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->atm_poltype,10,"%8d");

		os << "%FLAG AMOEBA_ATOMIC_NUMBER" << std::endl;  // ATOM ELEMENTS
        os << "%COMMENT   dimention = (" << p_amber_model->natom << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->atm_element,10,"%8d");

		os << "%FLAG AMOEBA_ATOM_CLASS_INDEX" << std::endl;
        os << "%COMMENT   dimention = (" << p_amber_model->natom << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		iarr.resize(p_amber_model->natom);
	
		for(i = 0; i < p_amber_model->natom; i++)
		{
			try { iarr[i] = boost::lexical_cast<int>(p_amber_model->atm_class_idx[i]); }
			catch(boost::bad_lexical_cast&)
			{
				PrintLog("Error in MMDriverAmber::SaveAmberTopToStream() \n");
				PrintLog("Can not convert atom class type to integer %s\n",p_amber_model->atm_class_idx[i].c_str());
				iarr[i] = 0;
			}
		}
		write_int_array_chuncks(os,iarr,10,"%8d");

		os << "%FLAG AMOEBA_REGULAR_BOND_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
		os << std::setw(8) << p_amber_model->n_bond_amoeba << std::endl;
		      
        os << "%FLAG AMOEBA_REGULAR_BOND_LIST" << std::endl;
        os << "%COMMENT dimension = (3," << p_amber_model->n_bond_amoeba << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->gbl_bond_amoeba,10,"%8d");

		os << "%FLAG AMOEBA_REGULAR_BOND_NUM_PARAMS" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
		os << std::setw(8) << p_amber_model->n_bond_amoeba_params << std::endl;

		os << "%FLAG " << "AMOEBA_REGULAR_BOND_FORCE_CONSTANT" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->bond_amoeba_params[0],5,FLOAT_E16_8);

		os << "%FLAG " << "AMOEBA_REGULAR_BOND_EQUIL_VALUE" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->bond_amoeba_params[1],5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_REGULAR_BOND_FTAB_DEGREE" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->bond_amoeba_ftab_degree << std::endl;

		os << "%FLAG AMOEBA_REGULAR_BOND_FTAB_COEFFS" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->bond_amoeba_ftab_coef,5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_UREY_BRADLEY_BOND_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_urey_bond << std::endl;

		os << "%FLAG AMOEBA_UREY_BRADLEY_BOND_LIST" << std::endl;
        os << "%COMMENT   dimension = (3," << p_amber_model->n_urey_bond << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->gbl_bond_urey,10,"%8d");

		os << "%FLAG AMOEBA_UREY_BRADLEY_BOND_NUM_PARAMS" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_urey_bond_params << std::endl;

		os << "%FLAG " << "AMOEBA_UREY_BRADLEY_BOND_FORCE_CONSTANT" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->bond_urey_params[0],5,FLOAT_E16_8);

		os << "%FLAG " << "AMOEBA_UREY_BRADLEY_BOND_EQUIL_VALUE" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->bond_urey_params[1],5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_UREY_BRADLEY_BOND_FTAB_DEGREE" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->bond_urey_ftab_degree << std::endl;
	        
        os << "%FLAG AMOEBA_UREY_BRADLEY_BOND_FTAB_COEFFS" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->bond_urey_ftab_coef,5,FLOAT_E16_8);

		p_amber_model->n_angle_amoeba =  p_amber_model->gbl_angle_amoeba_reg.size()/4;
		os << "%FLAG AMOEBA_REGULAR_ANGLE_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw( 8 ) << p_amber_model->n_angle_amoeba << std::endl;
		
		os << "%FLAG AMOEBA_REGULAR_ANGLE_LIST" << std::endl;
        os << "%COMMENT dimension = (4," << p_amber_model->n_angle_amoeba << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->gbl_angle_amoeba_reg,10,"%8d");

		os << "%FLAG AMOEBA_REGULAR_ANGLE_NUM_PARAMS" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_angle_amoeba_params << std::endl;

		os << "%FLAG " << "AMOEBA_REGULAR_ANGLE_FORCE_CONSTANT" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->angle_amoeba_params[0],5,FLOAT_E16_8);

		os << "%FLAG " << "AMOEBA_REGULAR_ANGLE_EQUIL_VALUE" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->angle_amoeba_params[1],5,FLOAT_E16_8);

        os << "%FLAG AMOEBA_REGULAR_ANGLE_FTAB_DEGREE" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->angle_amoeba_ftab_degree << std::endl;
        
        os << "%FLAG AMOEBA_REGULAR_ANGLE_FTAB_COEFFS" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->angle_amoeba_ftab_coef,5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_TRIGONAL_ANGLE_NUM_LIST" << std::endl;
		os << "%FORMAT(I8)" << std::endl;
		os << std::setw( 8 ) << p_amber_model->n_trig_angles << std::endl;

		os << "%FLAG AMOEBA_TRIGONAL_ANGLE_LIST" << std::endl;
        os << "%COMMENT dimension = (5," << p_amber_model->n_trig_angles << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->gbl_angle_amoeba_trig,10,"%8d");

		os << "%FLAG AMOEBA_TRIGONAL_ANGLE_NUM_PARAMS" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_angle_amoeba_params << std::endl;

		os << "%FLAG " << "AMOEBA_TRIGONAL_ANGLE_FORCE_CONSTANT" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->angle_amoeba_params[0],5,FLOAT_E16_8);

		os << "%FLAG " << "AMOEBA_TRIGONAL_ANGLE_EQUIL_VALUE" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->angle_amoeba_params[1],5,FLOAT_E16_8);

        os << "%FLAG AMOEBA_TRIGONAL_ANGLE_FTAB_DEGREE" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->angle_amoeba_ftab_degree << std::endl;
        
        os << "%FLAG AMOEBA_TRIGONAL_ANGLE_FTAB_COEFFS" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->angle_amoeba_ftab_coef,5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_OPBEND_ANGLE_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw( 8 ) << p_amber_model->n_opbend_angles << std::endl;

		os << "%FLAG AMOEBA_OPBEND_ANGLE_LIST" << std::endl;
        os << "%COMMENT dimension = (5," << p_amber_model->n_opbend_angles << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->gbl_opbend_angle,10,"%8d");

		os << "%FLAG AMOEBA_OPBEND_ANGLE_NUM_PARAMS" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_opbend_angles_params << std::endl;

		os << "%FLAG " << "AMOEBA_OPBEND_ANGLE_FORCE_CONSTANT" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->opbend_angle_params,5,FLOAT_E16_8);
            
        os << "%FLAG AMOEBA_OPBEND_ANGLE_FTAB_DEGREE" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->angle_amoeba_ftab_degree << std::endl;
        
        os << "%FLAG AMOEBA_OPBEND_ANGLE_FTAB_COEFFS" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->angle_amoeba_ftab_coef,5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_TORSION_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_tors_amoeba << std::endl;
        
        os << "%FLAG AMOEBA_TORSION_LIST" << std::endl;
        os << "%COMMENT   dimension = (5," << p_amber_model->n_tors_amoeba << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->gbl_amoeba_tors_angle,10,"%8d");

		os << "%FLAG AMOEBA_TORSION_NUM_PARAMS" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_tors_amoeba_params << std::endl;

		os << "%FLAG " << "AMOEBA_TORSION_FORCE_CONSTANT" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->tors_amoeba_params[0],5,FLOAT_E16_8);

		os << "%FLAG " << "AMOEBA_TORSION_PERIODICITY" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->tors_amoeba_params[1],5,FLOAT_E16_8);

		os << "%FLAG " << "AMOEBA_TORSION_PHASE" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->tors_amoeba_params[2],5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_PI_TORSION_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_pi_torsions << std::endl;
        
        os << "%FLAG AMOEBA_PI_TORSION_LIST" << std::endl;
        os << "%COMMENT   dimension = (7," << p_amber_model->n_pi_torsions << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->gbl_pi_tors_angle,10,"%8d");

		os << "%FLAG AMOEBA_PI_TORSION_NUM_PARAMS" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_pi_torsions_params << std::endl;

		os << "%FLAG " << "AMOEBA_PI_TORSION_FORCE_CONSTANT" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->pi_tors_params[0],5,FLOAT_E16_8);

		os << "%FLAG " << "AMOEBA_PI_TORSION_PERIODICITY" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->pi_tors_params[1],5,FLOAT_E16_8);

		os << "%FLAG " << "AMOEBA_PI_TORSION_PHASE" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->pi_tors_params[2],5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_STRETCH_BEND_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_stretch_bend << std::endl;
            
        os << "%FLAG AMOEBA_STRETCH_BEND_LIST" << std::endl;
        os << "%COMMENT   dimension = (4," << p_amber_model->n_stretch_bend << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->gbl_str_bend_angle,10,"%8d");

		os << "%FLAG AMOEBA_STRETCH_BEND_NUM_PARAMS" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_stretch_bend_params << std::endl;

		os << "%FLAG " << "AMOEBA_STRETCH_BEND_FORCE_CONSTANT" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->str_bend_params[0],5,FLOAT_E16_8);

		os << "%FLAG " << "AMOEBA_STRETCH_BEND_ANGLE_EQUIL_VALUE" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->str_bend_params[1],5,FLOAT_E16_8);

		os << "%FLAG " << "AMOEBA_STRETCH_BEND_BOND1_EQUIL_VALUE" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->str_bend_params[2],5,FLOAT_E16_8);

		os << "%FLAG " << "AMOEBA_STRETCH_BEND_BOND2_EQUIL_VALUE" << endl;
        os << "%FORMAT(5E16.8)" << endl;
		write_double_array_chuncks(os,p_amber_model->str_bend_params[3],5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_TORSION_TORSION_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_tors_tors << std::endl;

        os << "%FLAG AMOEBA_TORSION_TORSION_LIST" << std::endl;
        os << "%COMMENT   dimension = (6," << p_amber_model->n_tors_tors << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->gbl_tors_tors,10,"%8d");

		os << "%FLAG AMOEBA_TORSION_TORSION_NUM_PARAMS" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_tors_tors_params << std::endl;
            
		for( i = 0; i < p_amber_model->n_tors_tors_params; i++)
		{
			std::string prefix = "AMOEBA_TORSION_TORSION_TORTOR_TABLE_0";

			int id_par = p_amber_model->tors_tors_id_params[i];
			int ndim1 = p_amber_model->tors_tors_params[i][0].size();
			int ndim2 = p_amber_model->tors_tors_params[i][1].size();
			os << "%FLAG " << prefix << id_par << "_DIMS" << std::endl;
            os << "%COMMENT dimension = (2)" << std::endl;
            os << "%FORMAT(2I8)" << std::endl;
            os << std::setw(8) << ndim1;
            os << std::setw(8) << ndim2 << std::endl;

			os << "%FLAG " << prefix << id_par << "_ANGLE1" << std::endl;
            os << "%COMMENT   dimension = (25)" << std::endl;
            os << "%FORMAT(5E16.8)" << std::endl;
			write_double_array_chuncks(os,p_amber_model->tors_tors_params[i][0],5,FLOAT_E16_8);

			os << "%FLAG " << prefix << id_par << "_ANGLE2" << std::endl;
            os << "%COMMENT   dimension = (25)" << std::endl;
            os << "%FORMAT(5E16.8)" << std::endl;
			write_double_array_chuncks(os,p_amber_model->tors_tors_params[i][1],5,FLOAT_E16_8);

			os << "%FLAG " << prefix << id_par << "_FUNC" << std::endl;
            os << "%COMMENT   dimension = (25,25)" << std::endl;
            os << "%FORMAT(5E16.8)" << std::endl;
			write_double_array_chuncks(os,p_amber_model->tors_tors_params[i][2],5,FLOAT_E16_8);

			os << "%FLAG " << prefix << id_par << "_DFUNC_DANGLE1" << std::endl;
            os << "%COMMENT   dimension = (25,25)" << std::endl;
            os << "%FORMAT(5E16.8)" << std::endl;
			write_double_array_chuncks(os,p_amber_model->tors_tors_params[i][3],5,FLOAT_E16_8);

			os << "%FLAG " << prefix << id_par << "_DFUNC_DANGLE2" << std::endl;
            os << "%COMMENT   dimension = (25,25)" << std::endl;
            os << "%FORMAT(5E16.8)" << std::endl;
			write_double_array_chuncks(os,p_amber_model->tors_tors_params[i][4],5,FLOAT_E16_8);

			os << "%FLAG " << prefix << id_par << "_D2FUNC_DANGLE1_DANGLE2" << std::endl;
            os << "%COMMENT   dimension = (25,25)" << std::endl;
            os << "%FORMAT(5E16.8)" << std::endl;
			write_double_array_chuncks(os,p_amber_model->tors_tors_params[i][5],5,FLOAT_E16_8);
		}
		os << "%FLAG AMOEBA_VDW_ATOM_TYPES_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << p_amber_model->natom << std::endl;
            
        os << "%FLAG AMOEBA_VDW_ATOM_TYPES_LIST" << std::endl;
        os << "%COMMENT   dimension = (1," << p_amber_model->natom << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->atm_amoeba_vdw_type,10,"%8d");

		os << "%FLAG AMOEBA_VDW_ATOM_PARENT_LIST" << std::endl;
        os << "%COMMENT   dimension = (1," << p_amber_model->natom << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->atm_parent_id,10,"%8d");

		os << "%FLAG AMOEBA_VDW_PARENT_COORD_WEIGHT_LIST" << std::endl;
        os << "%COMMENT   dimension = (1," << p_amber_model->natom << ")" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->atm_parent_weight,5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_VDW_BUFFER_DELTA" << std::endl;
        os << "%FORMAT(E16.8)" << std::endl;
		farr.resize(1); farr[0] = p_amber_model->vdw_buffer_delta;
		write_double_array_chuncks(os,farr,5,FLOAT_E16_8);

        os << "%FLAG AMOEBA_VDW_BUFFER_GAMMA" << std::endl;
        os << "%FORMAT(E16.8)" << std::endl;
		farr.resize(1); farr[0] = p_amber_model->vdw_buffer_gamma;
		write_double_array_chuncks(os,farr,5,FLOAT_E16_8);
        
        os << "%FLAG AMOEBA_VDW_PARAMS_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->n_vdw_params << std::endl;

		os << "%FLAG AMOEBA_VDW_MIXED_RADII_LIST" << std::endl;
        os << "%COMMENT   dimension = (" << p_amber_model->n_vdw_params << "," << p_amber_model->n_vdw_params << ")" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->amoeba_vdw_rstars,5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_VDW_MIXED_EPSILONS_LIST" << std::endl;
        os << "%COMMENT   dimension = (" << p_amber_model->n_vdw_params << "," << p_amber_model->n_vdw_params << ")" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->amoeba_vdw_depths,5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_LOCAL_FRAME_MULTIPOLES_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->num_local_multipoles << std::endl;
        
        os << "%FLAG AMOEBA_LOCAL_FRAME_MULTIPOLES_LIST" << std::endl;
        os << "%COMMENT dimension = (10," << p_amber_model->num_local_multipoles << ")" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->atm_multipoles,5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_CHIRAL_FRAME_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->num_chiral_frames << std::endl;
                
        os << "%FLAG AMOEBA_CHIRAL_FRAME_LIST" << std::endl;
        os << "%COMMENT   dimension = (3," << p_amber_model->num_chiral_frames << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->atm_chiral_frames,10,"%8d");

		os << "%FLAG AMOEBA_FRAME_DEF_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->num_reg_frames << std::endl;
        
        os << "%FLAG AMOEBA_FRAME_DEF_LIST" << std::endl;
        os << "%COMMENT   dimension = (5," << p_amber_model->num_reg_frames << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->atm_reg_frames,10,"%8d");

		os << "%FLAG AMOEBA_ADJUST_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->num_adjust_list << std::endl;
            
        os << "%FLAG AMOEBA_ADJUST_LIST" << std::endl;
        os << "%COMMENT dimension = (3," <<  p_amber_model->num_adjust_list << ")" << std::endl;
        os << "%FORMAT(10I8)" << std::endl;
		write_int_array_chuncks(os,p_amber_model->atm_adjust_list,10,"%8d");

		os << "%FLAG AMOEBA_ADJUST_VDW_WEIGHTS_LIST" << std::endl;
        os << "%COMMENT   dimension = (9)" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->adjust_vdw_weights,5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_ADJUST_MPOLE_WEIGHTS_LIST" << std::endl;
        os << "%COMMENT   dimension = (9)" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->adjust_mpole_weights,5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_ADJUST_DIRECT_WEIGHTS_LIST" << std::endl;
        os << "%COMMENT   dimension = (9)" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->adjust_direct_weights,5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_ADJUST_POLAR_WEIGHTS_LIST" << std::endl;
        os << "%COMMENT   dimension = (9)" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->adjust_polar_weights,5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_ADJUST_MUTUAL_WEIGHTS_LIST" << std::endl;
        os << "%COMMENT   dimension = (9)" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->adjust_mutual_weights,5,FLOAT_E16_8);

		os << "%FLAG AMOEBA_POLARIZABILITY_NUM_LIST" << std::endl;
        os << "%FORMAT(I8)" << std::endl;
        os << std::setw(8) << p_amber_model->natom << std::endl;

        os << "%FLAG AMOEBA_POLARIZABILITY_LIST" << std::endl;
        os << "%COMMENT   dimension = (" << p_amber_model->natom << ")" << std::endl;
        os << "%FORMAT(5E16.8)" << std::endl;
		write_double_array_chuncks(os,p_amber_model->atm_polar,5,FLOAT_E16_8);

	} // !if(iamoeba)
		
    return TRUE;
}

int MMDriverAmber::CheckForStop()
{
	if( numtasks > 1 || p_mm_mod->run_ti ) MPI_Bcast(&p_mm_mod->to_stop_simulations, 1, MPI_INT , 0, p_mm_mod->single_job_comm);
	if( p_mm_mod->to_stop_simulations && mytaskid == 0) PrintLog (" Terminating MM Calculations: \n"); 
	return p_mm_mod->to_stop_simulations;
}

int MMDriverAmber::OpenOutputFiles()
{
	if(!master) return TRUE;
	if( p_mm_mod->run_ti ) p_mm_mod->p_ti_mod->SetTI_OutputFileNames();
	strcpy_to_fort(FC_FUNC_MODULE(file_io_dat_mod,mdout_name),amber_out_file.c_str(),80);
	strcpy_to_fort(FC_FUNC_MODULE(file_io_dat_mod,restrt_name),amber_rst_file.c_str(),80);
	strcpy_to_fort(FC_FUNC_MODULE(file_io_dat_mod,mdcrd_name),amber_trj_coord_file.c_str(),80);
	strcpy_to_fort(FC_FUNC_MODULE(file_io_dat_mod,mdvel_name),amber_trj_vel_file.c_str(),80);
	strcpy_to_fort(FC_FUNC_MODULE(file_io_dat_mod,mden_name),amber_trj_ene_file.c_str(),80);
	
	strcpy_to_fort(FC_FUNC_MODULE(file_io_dat_mod,prmtop_name),amber_top_file.c_str(),80);

	start_mdout_log_(pbc_box.v());
	PrintMDCntrData();

	if( p_mm_mod->run_ti )
	{
		PrintLogMDOUT("--------------------------------------------------------------------");
		PrintLogMDOUT(" Free energy  Thermodynamic Integration calculations \n");
		if( p_mm_mod->inter_model_rank == 1) 
		{
			PrintLogMDOUT(" \n");
			PrintLogMDOUT("     Axxiliary hamiltonian output \n");
			PrintLogMDOUT(" \n");
		}
		PrintLogMDOUT(" Current lambda index = %d   total number of lambdas = %d \n",
			            p_mm_mod->p_ti_mod->cur_idx_lmb,p_mm_mod->p_ti_mod->num_lmb);
		PrintLogMDOUT(" Current lambda value = %12.6f     Klambda = %8.2f \n",
			            p_mm_mod->p_ti_mod->GetCurLambda(),p_mm_mod->p_ti_mod->klambda_ti);
		PrintLogMDOUT(" \n");
		PrintLogMDOUT("--------------------------------------------------------------------");
		if( p_mm_mod->inter_model_rank != 0 ) return TRUE;
	}

	PrintLogMDOUT(" ");
	PrintLogMDOUT("File Assignments: \n");
	PrintLogMDOUT(" MDOUT = %s   RESTRT = %s\n", amber_out_file.c_str(),amber_rst_file.c_str());
	PrintLogMDOUT(" REFC =    refc.dat \n");
	PrintLogMDOUT(" MDVEL = %s   MDEN = %s\n",amber_trj_vel_file.c_str(),amber_trj_ene_file.c_str());
	PrintLogMDOUT(" MDCRD = %s   \n", amber_trj_coord_file.c_str());

	if ( p_mm_mod->traj_wrt_format == p_mm_mod->traj_wrt_format.FORMATTED )               // Formatted dumping:
	{
		if (p_mm_mod->wrt_coord_freq > 0) open_mdcrd_form_();
		if (p_mm_mod->wrt_vel_freq > 0) open_mdvel_form_();
	}
	else if (p_mm_mod->traj_wrt_format == p_mm_mod->traj_wrt_format.BINARY ) 
	{
		PrintLog(" No binary MD trajectory save in HARLEM \n");
		// FC_FUNC_MODULE(bintraj_mod,open_binary_files)();
	}
	else if ( p_mm_mod->traj_wrt_format.value() == 2)   // The new "bin4" efficiency format...
	{
		if (p_mm_mod->wrt_coord_freq > 0) open_mdcrd_bin4_( &p_mm_mod->period_bcond.value(), &p_mm_mod->limit_wrt_atoms );
		if (p_mm_mod->wrt_vel_freq > 0) open_mdvel_bin4_( &p_mm_mod->limit_wrt_atoms );
	}

// Open the energies file:

	if (p_mm_mod->wrt_ener_freq > 0) open_mden_();

	// Open the restart file:

	if ( p_mm_mod->write_coord_format == p_mm_mod->write_coord_format.BINARY ) 
	{		
		open_restart_bin_(); 
	}
	else
	{
		open_restart_form_();
	}

	if(p_mm_mod->run_ti) 
	{
		std::string dvdl_file_name = p_mm_mod->p_ti_mod->GetCurFilePrefix();
		dvdl_file_name += "_dvdl.dat";
		p_mm_mod->p_ti_mod->file_dvdl = fopen(dvdl_file_name.c_str(),"w");
	}
	if( numtasks > 1 ) PrintLogMDOUT("| Running MD/MIN simulations on %d  nodes \n\n", numtasks);

	PrintLogMDOUT("\n---------------\n  RESULTS \n-----------------");

	return TRUE;
} 

int MMDriverAmber::CloseOutputFiles()
{
	if(!master) return TRUE;
	PrintLog(" MMDriverAmber::CloseOutputFiles() pt 1 \n ");
	FC_FUNC_MODULE(runfiles_mod,close_mdout)(); 
	if( p_mm_mod->run_ti )
	{	
		if(p_mm_mod->single_job_rank == 0)
		{
			fclose(p_mm_mod->p_ti_mod->file_dvdl);
		}
		else
		{
			return TRUE;
		}
	}
	if ( p_mm_mod->traj_wrt_format == p_mm_mod->traj_wrt_format.BINARY ) 
	{
		PrintLog(" No binary MD trajectory save in HARLEM \n");
//		FC_FUNC_MODULE(bintraj_mod,close_binary_files)();
	}
	else
	{
		if (p_mm_mod->wrt_coord_freq > 0) close_mdcrd_();
		if (p_mm_mod->wrt_vel_freq > 0) close_mdvel_();
	}
	if (p_mm_mod->wrt_ener_freq > 0) close_mden_();
	close_restart_();
	return TRUE;
}

int MMDriverAmber::SaveAmberRstFile(const char* fname)
{
	char buf[256];
	ofstream os(fname);

	if(os.fail())
	{	
		PrintLog(" Error in MMDriverAmber::SaveAmberRstFile() \n");
		PrintLog(" Can not open file %s fro writing \n",fname);
		return FALSE;
	}
	std::string title;
	title = " Amber Coordinate file for the molecule set ";
	title += pmset->GetName();

	int extra_spaces = 80 - title.size();
	for(; extra_spaces > 0; extra_spaces-- )
		title += " ";
	
	os << title << std::endl;
	int npt = p_mm_model->Atoms.size();
	float time =0.0;
	sprintf(buf,"%5d%15.7E",npt, time);
	os << buf << std::endl;
	HaVec_double farr;
	farr.newsize(3*npt);
	int i;
	for( i = 0; i < npt; i++)
	{
		farr[3*i]   = p_mm_model->Atoms[i]->GetX_Ang();
		farr[3*i+1] = p_mm_model->Atoms[i]->GetY_Ang();
		farr[3*i+2] = p_mm_model->Atoms[i]->GetZ_Ang();
	}

	write_double_array_chuncks(os,farr,6,FLOAT_F12_7); // save coordinates

	if( pmset->per_bc->IsSet() )
	{
		farr.newsize(6);
		farr[0] = pmset->per_bc->GetA(); 
		farr[1] = pmset->per_bc->GetB();
		farr[2] = pmset->per_bc->GetC() ; 
		farr[3] = pmset->per_bc->GetAlpha() * RAD_TO_DEG; 
		farr[4] = pmset->per_bc->GetBeta() * RAD_TO_DEG; 
		farr[5] = pmset->per_bc->GetGamma() * RAD_TO_DEG;
		write_double_array_chuncks(os,farr,6,FLOAT_F12_7); // save periodical box description
	}
	return TRUE;
}

void MMDriverAmber::RunMinSlave()
{
//	PrintLog(" MMDriverAmber::RunMinSlave() pt 1 \n");
	int ncalls;
	int new_list;
	int not_first_loop;

	p_tm->ZeroTime();

	InitAddMDCtrlParams();

	sys_info = 0.0;

	not_first_loop = FALSE;
	ncalls = 0;

	for(;;)
	{
		if ( CheckForStop() ) 
		{
			return;
		}
		p_tm->UpdateTime(TimerAmber::FCVE_DIST_TIME);

     // For minimization, everyone gets a full set of crds from the master;
     // minimization is unfortunately not well parallelized.

		int ierr = MPI_Bcast(atm_crd.v(), 3 * p_amber_model->natom, MPI_DOUBLE, 0, driver_mpi_comm);
		          
		p_tm->UpdateTime(TimerAmber::FCVE_DIST_TIME);

		if( p_amber_model->using_pme_potential) 
		{
			if (not_first_loop) 
			{
        // Do a skin check to see if we will have to rebuild the pairlist.
				new_list = CheckAllAtomMovement(p_amber_model->natom, atm_crd);
				FC_FUNC_MODULE(loadbal_mod,check_new_list_limit)(&new_list);
				p_tm->UpdateTime(TimerAmber::NONBOND_TIME);
			}
			else
			{
				new_list = TRUE;
				not_first_loop = TRUE;
			}
		}
		ncalls = ncalls + 1;

		CalcForceAndEne(new_list,ncalls);

    // Potential energy info not used in slaves...
		if (p_amber_model->using_pme_potential)
		{
			FC_FUNC_MODULE(parallel_mod,mpi_gathervec)(&p_amber_model->natom, atm_frc.v(), gbl_atm_owner_map.v(), 
				                           gbl_my_atm_lst.v(), &my_atm_cnt);
		}
		else if (p_amber_model->using_gb_potential)
		{
			FC_FUNC_MODULE(parallel_mod,gb_mpi_gathervec)(&p_amber_model->natom, atm_frc.v());
		}
		p_tm->UpdateTime(TimerAmber::FCVE_DIST_TIME);
	}
}

void MMDriverAmber::RunMinMaster()
{	
  PrintLog(" Start Energy Minimization \n"); 
  OpenOutputFiles();

  double betax;
  double ddspln;
  double dfpr;
  double dxsth;
  double f;
  double fch;
  double fdmax;
  double finit;
  double gama;
  double gamden;
  double ginit;
  double gmin;
  double gnew;
  double gspln;
  double gsqrd;
  int i;
  int iatmax;
  int iretry;
  int iterrs;
  std::string labmax;
  int nfbeg;
  int nfopt;
  double  sbound;
  double  step;
  double  stepch;
  double  stmin;
  double  sum;
  double  swork;
  double  work;

  int atm_cnt = p_amber_model->natom;

  ZMatCrd* pzmat = NULL;
  if( p_mm_mod->IsZMatMin() )
  {
	 pzmat = pmset->GetZMat();
	 if( pzmat->IsEmpty() ) pzmat->InitStdZMat();
  }

  int ncrd = 3*p_amber_model->natom;
  if( p_mm_mod->IsZMatMin() ) 
  {
	  ncrd = pzmat->GetNCrdUnFrozen();
  }

// These are the work arrays, dynamically allocated

  HaVec_double w( ncrd );
  HaVec_double w_xopt( ncrd );
  HaVec_double w_gopt( ncrd );
  HaVec_double w_ginit( ncrd );
  HaVec_double w_rsdg( ncrd );
  HaVec_double w_rsdx( ncrd );

  HaVec_double crd_loc( ncrd );
  HaVec_double frc_loc( ncrd );

// Parameters of minimization, move to memberes of Minimization Simulator?

  int maxlin = 10;
  int mxfcon = 4;
  int kstcyc = 4;
  double dxstm  = 1.0e-05;
  double crits  = 1.0e-06;
  double dfpred = 1.e0;

  p_tm->ZeroTime();
  
  p_mm_mod->to_stop_simulations = FALSE;
  InitAddMDCtrlParams();

  sys_info = 0.0;

// Evaluate some constants:

  double fmin = 0.0e0;
  int nct = 0;
  if (ntc == 2 || ntc == 4) nct = p_amber_model->nbonh;  // Number of frozen bonds
  if (ntc == 3) nct = p_amber_model->nbonh + p_amber_model->nbona;      
  int ndfp = 3*p_amber_model->natom - nct;             // number of not frozen degrees of freedom
  if (p_amber_model->ibelly > 0) ndfp = 3 * p_amber_model->belly_atm_cnt - nct;
  int ntnb = 1;
  double fnq = sqrt((double)ndfp);
  double rms = 0.0e0;
  int steep  = FALSE;
  int newstr = FALSE;
  int nstcyc = 0;
  int mstcyc = kstcyc;
  if (p_mm_mod->min_type == p_mm_mod->min_type.STEEPEST_DESCENT) mstcyc = p_mm_mod->max_num_minim_steps;
  if (p_mm_mod->min_type == p_mm_mod->min_type.SD_AND_CG) mstcyc = p_mm_mod->num_steep_descent_steps;
  if (p_mm_mod->min_type != p_mm_mod->min_type.CONJ_GRAD) steep  = TRUE;
  double fold = 0.0e0;
  double dxst = p_mm_mod->init_min_step;
  int linmin = 0;

// Set some parameters to begin the calculation:

  int iterc = 0;
  int ncalls = 0;
  int iterfm = iterc;

// Let the initial search direction be minus the gradient vector. iterrs gives
// the iteration number of the most recent restart, but is set to zero when
// steepest descent direction is used:

// (Here is the beginning of a big loop:)

	int not_first_loop = FALSE;
	int new_list = TRUE;

	if( p_mm_mod->IsZMatMin() )
	{
		pzmat->GetElemCrdVal( crd_loc );
	}
	else
	{
		crd_loc = atm_crd;
	}

	for(;;) // Main LOOP: 
	{
		if( p_mm_mod->IsZMatMin() ) 
		{
			pzmat->SetFromElemCrdVal( crd_loc );
			pzmat->GetCrdSnapshot( atm_crd );
		}
		else
		{
			atm_crd = crd_loc;
		}

//		PrintLog(" MMDriverAmber::RunMinMaster() pt 3  ncrd = %d \n", ncrd );
//		return;

		if( CheckForStop() ) break;

		ncalls = ncalls + 1;
		if (ncalls % nsnb == 0) ntnb = 1;
		if (ntnb == 1 && ncalls > 1) steep = TRUE;

// Calculate the force and energy:

		p_tm->UpdateTime(TimerAmber::OTHER_TIME);

		int ierr;
  
  // For minimization, everyone gets a full set of crds from the master;
  // minimization is unfortunately not well parallelized.

		if( numtasks > 1 )
		{	
			ierr = MPI_Bcast(atm_crd.v(), 3 * p_amber_model->natom, MPI_DOUBLE, 0, driver_mpi_comm);
			p_tm->UpdateTime(TimerAmber::FCVE_DIST_TIME);
		}

		if (p_amber_model->using_pme_potential)
		{
			if (not_first_loop) 
			{
      // Now do a skin check to see if we will have to rebuild the pairlist.
				new_list = CheckAllAtomMovement(p_amber_model->natom, atm_crd);
				if(numtasks > 1 )
				{
					FC_FUNC_MODULE(loadbal_mod,check_new_list_limit)(&new_list);
				}
				p_tm->UpdateTime(TimerAmber::NONBOND_TIME);
			}
			else
			{
				new_list = TRUE;
				not_first_loop = TRUE; 
			} // (not_first_loop
		}

		CalcForceAndEne(new_list,ncalls);

		if(numtasks > 1)     // combine forces from different mpi nodes
		{
			if (p_amber_model->using_pme_potential)
			{
				FC_FUNC_MODULE(parallel_mod,mpi_gathervec)(&p_amber_model->natom, atm_frc.v(), gbl_atm_owner_map.v(), 
					                            gbl_my_atm_lst.v(), &my_atm_cnt);
			}
			else if (p_amber_model->using_gb_potential) 
			{
				FC_FUNC_MODULE(parallel_mod,gb_mpi_gathervec)(&p_amber_model->natom, atm_frc.v());
			}
			p_tm->UpdateTime(TimerAmber::FCVE_DIST_TIME);
		}
		
		GetAtomCrdFromInternalArrays();
		
		SetMMInfo((*p_mm_mod->p_mm_info),sys_info);
		p_mm_model->UpdateDataFromFort();
		p_mm_mod->p_mm_info->nstep = ncalls;

		f = sys_info[SI_POT_ENE];
		ntnb = 0;

		if( p_mm_mod->IsZMatMin() )
		{ 
			pzmat->GetElemCrdVal( crd_loc );
			pzmat->TransDerivToIntCrd( atm_frc, frc_loc );
		}
		else
		{
			crd_loc = atm_crd;
			frc_loc = atm_frc;
		}

		sum = dot_product(frc_loc,frc_loc);
		rms = sqrt(sum)/fnq;

// Print the intermediate results:

		if (  ncalls % p_mm_mod->wrt_log_freq == 0  || ncalls == 1 ) 
		{
			GrdMax(frc_loc,iatmax,fdmax);
			if( p_mm_mod->IsZMatMin() ) 
			{
				ElemCrd* pcrd = pzmat->GetCrdByIdx(iatmax);
				labmax = pcrd->GetTag();
			}
			else
			{
				iatmax = iatmax/3 + 1;
				labmax = p_amber_model->atm_igraph[iatmax-1];
			}
			PrintMinEneMDOUT(sys_info,ncalls,rms,fdmax,iatmax,labmax);
		}

// Do some steepest steps before entering the conjugate gradient method:

		if (steep)
		{
			nstcyc = nstcyc + 1;
			if (nstcyc <= mstcyc)
			{
				if (dxst <= crits) dxst = dxstm;
				dxst = dxst/2.0e0;
				if (f < fold) dxst = dxst*2.4e0;
				dxsth = dxst/sqrt(sum);
				if (nstcyc <= 1 || f <= fmin) 
				{
					fmin = f;
					nfopt = ncalls;
					w_xopt = crd_loc;
					w_gopt = frc_loc; w_gopt.scale(-1.0);
				}
// Check for convergence:
				if (rms <= p_mm_mod->grad_cnvrg_val) 
				{
					PrintLogMDOUT("  Steepest Descent Energy Minimization Converged. \n");
					p_mm_mod->to_stop_simulations = TRUE;	
				}

				if (ncalls >= p_mm_mod->max_num_minim_steps) 
				{
					PrintLogMDOUT("  Maximum number of minimization cycles reached. \n");
					p_mm_mod->to_stop_simulations = TRUE;
				}

				if( !p_mm_mod->to_stop_simulations )
				{
					fold = f;
					crd_loc += frc_loc*dxsth;
				}
				continue;
			}
			else
			{
// (arrive here when finished with this set of steepest descent cycles)
				steep = FALSE;
				newstr = TRUE;
				nstcyc = 0;
				mstcyc = kstcyc;
			} // if(nstcyc <= mstcyc)
		}   // if(steep)

// Start of conjugate gradient steps:
		frc_loc.scale(-1.0);

		if (!newstr && ncalls >= 2) goto LBL_82;

LBL_70:

		w = frc_loc; w.scale(-1.0);

		iterrs = 0;
		if (newstr) iterc = 0;

		if (iterc > 0) goto LBL_140;

LBL_82:

		gnew = dot_product(w,frc_loc);

		if (newstr || (ncalls == 1)) goto LBL_100;

		fch = f - fmin;

// Store the values of crd, f and g, if they are the best that have been
// calculated so far. Test for convergence:

		if(fch < 0.0 ) goto LBL_100;
		if(fch == 0.0 )goto LBL_90;
		if(fch > 0.0)  goto  LBL_130;
  
LBL_90:

		if ( gnew/gmin < -1.0) goto LBL_120;

LBL_100:

		fmin = f;
		gsqrd = sum;
		nfopt = ncalls;

		w_xopt = crd_loc;
		w_gopt = frc_loc;

LBL_120: 

		if (rms <= p_mm_mod->grad_cnvrg_val) 
		{
			PrintLogMDOUT(" Conjugate Gradient Energy minimization converged. \n");
			p_mm_mod->to_stop_simulations = TRUE;
			continue;
		}

// Test if the value of p_mm_mod->max_num_minim_steps allows another call of funct:

LBL_130:
		if (ncalls >= p_mm_mod->max_num_minim_steps) 
		{
			PrintLogMDOUT("  Maximum number of minimization cycles reached. \n");
			p_mm_mod->to_stop_simulations = TRUE;
			continue;
		}

		if (!newstr && ncalls > 1) goto LBL_180;

//  This section is executed at the beginning of a conjugate gradient set of
// minimization steps.

// Set dfpr to p_mm_mod->init_min_step*gsqrd. dfpr is the reduction in the function value. stmin is
// usually the step-length of the most recent line search that gives the least
// value of f:

// dac change, 10/91:  Return to original idea of trying to go downhill by the
//                     absolute amount, dfpred (which defaults to 1 kcal/mol,
//                     see data statement above).  This can eliminate very bad
//                    initial conjugate gradient steps.

		dfpr = dfpred;
		stmin = dfpred/gsqrd;

		newstr = FALSE;

// Begin the main conguate gradient iteration:

LBL_140:

		iterc = iterc + 1;

		finit = f;
		ginit = 0.0;

		w_ginit = frc_loc;

		ginit = dot_product(w, frc_loc);

		if (ginit >= 0.0) goto LBL_260;
		gmin = ginit;
		sbound = - 1.0;
		nfbeg = ncalls;
		iretry = - 1;

		stepch = MinFun(stmin, fabs(dfpr/ginit));
		stmin = dxstm;

LBL_160:
		step = stmin + stepch;
		dxst = step;
		crd_loc = w_xopt; crd_loc += w * stepch;

		swork = 0.0;
		for(i = 0; i < ncrd; i++)
		{
			swork = MaxFun(swork, fabs(crd_loc[i] - w_xopt[i]));
		}

		if (swork > 0.0) continue;  // Strange I think the condition is always true 

// "work = swork" may not be needed - wont hurt.  -gls

		work = swork;

// Terminate the line search if stepch is effectively zero:

		if ((ncalls > (nfbeg + 1)) || (fabs(gmin/ginit) > 0.2)) 
		{
			PrintLogMDOUT("    .... RESTARTED DUE TO LINMIN FAILURE ... \n");
			steep = TRUE;
			linmin = linmin + 1;
		}

		goto LBL_270;

LBL_180: 

		work = (fch + fch)/stepch - gnew - gmin;
		ddspln = (gnew - gmin)/stepch;
		if (ncalls > nfopt)
		{
			sbound = step;
		}
		else
		{
			if (gmin*gnew <= 0.0) sbound = stmin;
			stmin = step;
			gmin = gnew;
			stepch = - stepch;
		}

		if (fch != 0.0) ddspln = ddspln + (work + work)/stepch;

// Test for convergence of the line search, but force at least two steps to be
// taken in order not to lose quadratic termination:

		if (gmin == 0.0) goto LBL_270;
		if (ncalls <= (nfbeg + 1)) goto LBL_200;
		if (fabs(gmin/ginit) <= 0.2) goto LBL_270;

// Apply the test that depends on the parameter maxlin:

LBL_190:

		if (ncalls < (nfopt + maxlin) ) goto LBL_200;

// Possible non bonded update. make a restart:

		PrintLogMDOUT(" .... RESTARTED DUE TO LINMIN FAILURE ...\n");
		steep = TRUE;
		linmin = linmin + 1;

		goto LBL_270;

LBL_200:

		stepch = 0.5e0*(sbound - stmin);
		if (sbound < - 0.5e0) stepch = 9.0e0*stmin;
		gspln = gmin + stepch*ddspln;
		if (gmin*gspln < 0.0e0) stepch = stepch*gmin/(gmin - gspln);
   
		goto LBL_160;

// Calculate the value of betax in the new direction:

LBL_210:

		sum = dot_product(frc_loc, w_ginit);
		betax = (gsqrd-sum)/(gmin - ginit);

// Test that the new search direction can be made downhill.  If not then try to
// improve the accuracy of the line search:

		if (fabs(betax*gmin) < 0.2 * gsqrd) goto LBL_220;
		iretry = iretry + 1;
		if (iretry <= 0) goto LBL_190;

LBL_220:

		if (f < finit) iterfm = iterc;
		if (iterc >= (iterfm + mxfcon))
		{
			PrintLogMDOUT(" .... RESTARTED DUE TO LINMIN FAILURE ... \n");
			steep = TRUE;
			linmin = linmin + 1;
			goto LBL_270;
		}
		dfpr = stmin*ginit;

// Branch if a restart procedure is required due to the iteration number or due
// to the scalar product of consecutive gradients:

		if (iretry > 0) goto LBL_70;
		if (iterrs == 0) goto LBL_240;
		if ((iterc - iterrs) >= ncrd ) goto LBL_240;
		if (fabs(sum) >= 0.2e0*gsqrd) goto LBL_240;

// Calculate gama in the new search direction. gamden is set by the restart
// procedure:

		gama = dot_product(frc_loc, w_rsdg);
		sum  = dot_product(frc_loc, w_rsdx);
		gama = gama/gamden;

// Restart if the new search direction is not sufficiently downhill:

		if (fabs(betax*gmin + gama*sum) >= 0.2e0*gsqrd) goto LBL_240;

// Calculate the new search direction:

		w = frc_loc*(-1.0) + w*betax + w_rsdx*gama;

// Cycle back for more conjugate gradient steps:

		goto LBL_140;

// Apply the restart procedure:

LBL_240: 

		gamden = gmin - ginit;

		w_rsdx = w;
		w_rsdg = frc_loc - w_ginit;
		w = frc_loc*(-1.0) + w*betax;

		iterrs = iterc;

		goto LBL_140;

LBL_260: 

		steep = TRUE;
		PrintLogMDOUT(" .... RESTARTED DUE TO LINMIN FAILURE ...\n");
		linmin = linmin + 1;

// Ensure that f, crd and g are optimal:

LBL_270:

		if (ncalls != nfopt) 
		{
			f = fmin;

			crd_loc = w_xopt;
			frc_loc = w_gopt;
		}

		if (linmin > 4) 
		{
			PrintLogMDOUT("   ***** REPEATED LINMIN FAILURE ***** \n");
			p_mm_mod->to_stop_simulations = TRUE;
			continue;
		}

		if (steep) continue;
		goto LBL_210;
	} // Main LOOP

	PrintLogMDOUT("                  FINAL RESULTS        \n");
	GrdMax(frc_loc,iatmax,fdmax);
	if( p_mm_mod->IsZMatMin() ) 
	{
		ElemCrd* pcrd = pzmat->GetCrdByIdx(iatmax);
		labmax = pcrd->GetTag();
	}
	else
	{
		iatmax = iatmax/3 + 1;
		labmax = p_amber_model->atm_igraph[iatmax-1];
	}
	PrintMinEneMDOUT(sys_info,ncalls,rms,fdmax,iatmax,labmax);

	CloseOutputFiles();
	GetAtomCrdFromInternalArrays();
}

// MMDriverAmber operates in kcal/mol units for energy, amu for masses,
// and angstoms for distances.  To convert the input time parameters
// from picoseconds to internal units, multiply by 20.455
// (which is 10.0 * sqrt(4.184)).


void MMDriverAmber::RunMD()
{
	if (p_mm_mod->run_ti)
	{
		if (this->master)
		{

			PrintLog(" Start TI Simulations \n\n");
			PrintLog(" Number of lambdas = %d \n", p_mm_mod->p_ti_mod->num_lmb);
			PrintLog(" Current lambda idx = %d   lambda val = %12.6f \n",
				p_mm_mod->p_ti_mod->cur_idx_lmb, p_mm_mod->p_ti_mod->GetCurLambda());
		}
	}
	else
	{
		if (this->master)
		{
			PrintLog(" Start MD Simulations \n");
		}
	}
	
	double etot_save;
	sys_info = 0.0;
	Vec3D sys_cnt; // Center of the system (Reset in Langevin dynamics)

	int reset_velocities;

	p_mm_mod->to_stop_simulations = FALSE;

// Need to move init_dynamics_dat() to C++ and plug in here...

	p_tm->ZeroTime();
	
	int& nstep = p_mm_mod->p_mm_info->nstep;
	nstep = 0;
	InitAddMDCtrlParams();

    OpenOutputFiles();
	p_mm_mod->p_traj_io_agent->Init();

// Langevin dynamics setup:

	if (p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN && !pmset->per_bc->IsSet() ) 
	{
		GetPosition(atm_crd,sys_cnt);
	}
	
// The following flags keep track of coordinate/velocity update status when
// using mpi

	all_crds_valid = TRUE;    
	all_vels_valid = TRUE;

	if (p_amber_model->ibelly > 0) ZeroVelFrozenAtoms();

// Make a first dynamics step:
// We must build a new pairlist the first time we run force:
	
	int nvalid = 0;
	int new_list = TRUE;

	// PrintLog("MMDriverAmber::RunMD() pt 1 \n");
	if (p_mm_mod->restart_flag == 0) 
	{
		p_tm->UpdateTime(TimerAmber::RUNMD_TIME);
		CalcForceAndEne(new_list,nstep);

    // The coordinates will not be changed between here and the next
    // run of force, so we can just set new_list to .false.

		new_list = FALSE;

    // Reset quantities depending on ref_temp and temp_relax_time_solute (which may have been 
    // changed by modwt during force call):

		CalcKinEne(TRUE);
		sys_info[SI_TOT_ENE] = sys_info[SI_KIN_ENE] + sys_info[SI_POT_ENE];

    // Propagate velocities half step back to prepare leap frog integration 

		PropagateVelHalfStepBack(atm_vel, atm_frc);

		// Print the initial energies and temperatures 
		if (nstep <= 0 && master) PrintMDEneMDOUT(sys_info,nstep,t);

		if (p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_BERENDSEN ) 
		{
			double kin_ene_save = sys_info[SI_KIN_ENE];
			CalcKinEne(TRUE);
			sys_info[SI_KIN_ENE_MINUS_HALF_DT] = MaxFun(sys_info[SI_KIN_ENE], R_HALF_KCAL*(p_amber_model->num_deg - ndfmin)*10.0);
			sys_info[SI_KIN_ENE] = kin_ene_save;
		}

		if (p_mm_mod->num_md_steps == 0) return;
	}  //  end of   if (restart_flag == 0) ( setup )  
	else
	{
// restart_flag = 1:  Continuation of a previous trajectory:
//
// Note: if the last printed energy from the previous trajectory was
//       at time "t", then the restrt file has velocities at time
//       t + 0.5dt, and coordinates at time t + dt

		CalcKinEne(TRUE);
		sys_info[SI_KIN_ENE_MINUS_HALF_DT] = sys_info[SI_KIN_ENE];
	}
		
	atm_last_vel = atm_vel;

//  !=======================================================================
//  ! MAIN LOOP FOR PERFORMING THE DYNAMICS STEP:
//  ! At this point, the coordinates are a half-step "ahead" of the velocities;
//  ! sys_info[SI_KIN_ENE_MINUS_HALF_DT] holds the kinetic energy at these "-1/2" velocities,
//  ! which are stored in the array last_vel.
//  !=======================================================================

	for(;;)
	{
		if( CheckForStop() ) break;
		
		p_tm->UpdateTime(TimerAmber::RUNMD_TIME);
		double kin_ene_minus_half_dt = sys_info[SI_KIN_ENE_MINUS_HALF_DT];

		CalcForceAndEne(new_list,nstep);
		// PrintLog("Nstep = %d \n", nstep);
		
		IncrementVelAndCrd(nstep, atm_crd, atm_vel, atm_frc, reset_velocities);

  // Scale coordinates if constant pressure run:
		if( ntp > 0 ) 
		{
			CalcPressure(sys_info);
			ScaleCoordConstPress(atm_crd, pbc_box, sys_info);
		}

		new_list = CheckForNewNonBondList();

		if( p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_BERENDSEN || (nstep % nrespa) == 0 )  
		  	  CalcKinEne();

		// restore the value of the kinetic energy before calling of CalcForceAndEne() 
		// sys_info is set to zero there - need to fix
		sys_info[SI_KIN_ENE_MINUS_HALF_DT] = kin_ene_minus_half_dt;

       // Now distribute the coordinates:

		if( numtasks > 1 ) CollectCoords( new_list, nstep, collect_crds, all_crds_valid);
		if( all_crds_valid && p_mm_mod->single_job_rank == 0 ) GetAtomCrdFromInternalArrays();

		if (p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_BERENDSEN) 
		{
			ScaleVelConstTemp();
			if( (nstep % nrespa) == 0 ) sys_info[SI_KIN_ENE_MINUS_HALF_DT] = MaxFun(sys_info[SI_KIN_ENE_PLUS_HALF_DT], R_HALF_KCAL * (p_amber_model->num_deg - ndfmin) * 10.0);
		}
   
    // If velocities were reset, the KE is not accurate; fudge it here to keep
    // the same total energy as on the previous step.  Note that this only
    // affects printout and averages for Etot and KE -- it has no effect on the
    // trajectory, or on any averages of potential energy terms.
   
		if (reset_velocities) sys_info[SI_KIN_ENE] = etot_save - sys_info[SI_POT_ENE];
   
    // --- total energy is sum of KE + PE:
   
		sys_info[SI_TOT_ENE] = sys_info[SI_KIN_ENE] + sys_info[SI_POT_ENE];
		etot_save = sys_info[SI_TOT_ENE];

		SetMMInfo((*p_mm_mod->p_mm_info),sys_info);

// For periodic PME, Zero COM velocity if requested; used for preventing ewald
// "block of ice flying thru space" phenomenon. We make this correction for pme
// explicitly before collecting velocities to the master...

		if ( (nstep + 1) % p_mm_mod->remove_rb_motion_freq == 0 && pmset->per_bc->IsSet() ) RemoveCOMVelocity();
		
// It may be necessary for the master to have a complete copy of the velocities
// for archiving and writing the restart file. In any case, all processors must
// have a complete copy of velocities if atom redistribution is going to happen
// in the next call to force.
		if( numtasks > 1 ) CollectVelocities( new_list, nstep, collect_crds, all_vels_valid);

// Zero COM velocity if requested; here we are doing the adjustment for a
// nonperiodic system, not pme.

		if ((nstep + 1) % p_mm_mod->remove_rb_motion_freq == 0 && !pmset->per_bc->IsSet() ) RemoveCOMVelAndResetCenter(all_crds_valid, all_vels_valid, sys_cnt);

// Zero out any non-moving velocities if a belly is active:

		if (p_amber_model->ibelly > 0) ZeroVelFrozenAtoms();

// Save old, or last velocities...

		if( p_mm_mod->run_ti && nstep % p_mm_mod->p_ti_mod->ti_sync_freq == 0) p_mm_mod->p_ti_mod->SincCrdAndVelTI();
		
		SaveVelToLastVel(all_vels_valid);
        
// Update the step counter and the integration time:

		nstep = nstep + 1;
		t = t + p_mm_mod->md_time_step;

// Output for this step if required:
	   
		p_mm_mod->p_traj_io_agent->AnalyzePt(nstep,sys_info);

		if (nstep >= p_mm_mod->num_md_steps)
		{
//			p_mm_mod->p_traj_io_agent->Finalize(nstep);
//			CloseOutputFiles();
			if (p_mm_mod->single_job_rank == 0) GetAtomCrdFromInternalArrays();

			if (p_mm_mod->run_ti)
			{
				if (master)
				{
					p_mm_mod->p_ti_mod->cur_idx_lmb++;
					if (p_mm_mod->p_ti_mod->cur_idx_lmb >= p_mm_mod->p_ti_mod->num_lmb) p_mm_mod->to_stop_simulations = TRUE;
					if (p_mm_mod->p_ti_mod->max_idx >= 0 && p_mm_mod->p_ti_mod->cur_idx_lmb > p_mm_mod->p_ti_mod->max_idx) p_mm_mod->to_stop_simulations = TRUE;
					if (!p_mm_mod->to_stop_simulations) nstep = 0;
				}
			}
			else
			{
				p_mm_mod->to_stop_simulations = TRUE;
			}
		}
		if (p_mm_mod->run_ti)
		{
			if (numtasks > 1) BCastAddMDCtrlParams(p_mm_mod->single_job_comm);
			if (!p_mm_mod->to_stop_simulations)
			{
				OpenOutputFiles();
				PrintLog(" Current lambda idx = %d   lambda val = %12.6f \n",
					p_mm_mod->p_ti_mod->cur_idx_lmb, p_mm_mod->p_ti_mod->GetCurLambda());
			}
		}
	}  // Major cycle back to new step unless we have reached our limit:
	p_mm_mod->p_traj_io_agent->Finalize(nstep);
    CloseOutputFiles();

// Print averages:
	p_tm->UpdateTime(TimerAmber::RUNMD_TIME);
}

void MMDriverAmber::CalcCurrEne()
{
	if( master ) SetAtomCrdToInternalArrays();
	p_tm->ZeroTime();
	InitAddMDCtrlParams();
	if (ntc != 1)
	{
		atm_frc = atm_crd;
		FC_FUNC_MODULE(shake_mod,shake)(atm_frc.v(), atm_crd.v(),&p_mm_mod->period_bcond.value(),
			                p_amber_model->atm_igroup.v(),pbc_box.v(),p_amber_model->atm_mass_inv.v(),&p_mm_mod->shake_tol);
		p_tm->UpdateTime(TimerAmber::SHAKE_TIME);
	}
	int ierr;
	int not_first_loop = FALSE;
	int ncalls = 1;
	int new_list = TRUE;

	if(numtasks > 1)
	{
		ierr = MPI_Bcast(atm_crd.v(), 3 * p_amber_model->natom, MPI_DOUBLE, 0, driver_mpi_comm);
		p_tm->UpdateTime(TimerAmber::FCVE_DIST_TIME);
	}
	
	if (p_amber_model->using_pme_potential)
	{
		if (not_first_loop) 
		{
			// Now do a skin check to see if we will have to rebuild the pairlist.
			new_list = CheckAllAtomMovement(p_amber_model->natom, atm_crd);
			if(numtasks > 1 )
			{
				FC_FUNC_MODULE(loadbal_mod,check_new_list_limit)(&new_list);
			}
			p_tm->UpdateTime(TimerAmber::NONBOND_TIME);
		}
		else
		{
			new_list = TRUE;
			not_first_loop = TRUE; 
		} // (not_first_loop
	}

	sys_info = 0.0;
	CalcForceAndEne(new_list,ncalls);

	if(numtasks > 1)     // combine forces from different mpi nodes
	{
		if (p_amber_model->using_pme_potential)
		{
			FC_FUNC_MODULE(parallel_mod,mpi_gathervec)(&p_amber_model->natom, atm_frc.v(), gbl_atm_owner_map.v(), 
				                            gbl_my_atm_lst.v(), &my_atm_cnt);
		}
		else if (p_amber_model->using_gb_potential) 
		{
			FC_FUNC_MODULE(parallel_mod,gb_mpi_gathervec)(&p_amber_model->natom, atm_frc.v());
		}
		p_tm->UpdateTime(TimerAmber::FCVE_DIST_TIME);
	}

	if(master) SetMMInfo((*p_mm_mod->p_mm_info),sys_info);
	if(master) p_mm_model->UpdateDataFromFort();
}

void MMDriverAmber::SetMMInfo(MMSysInfo& info,HaVec_double& si_vec)
{          
	info.press     = si_vec[SI_TOT_PRESS];
	info.pres_x    = si_vec[SI_PRESS_0];
	info.pres_y    = si_vec[SI_PRESS_1];;
	info.pres_z    = si_vec[SI_PRESS_2];
//	info.press_scale_solute  = 0.0;
//	info.press_scale_solvent = 0.0;

//	info.tot_energy = si_vec[SI_TOT_ENE];

	info.kin_ene  = si_vec[SI_KIN_ENE];
	info.pot_ene  = si_vec[SI_POT_ENE]; 
	info.bond_ene = si_vec[SI_BOND_ENE];
	info.vang_ene = si_vec[SI_ANGLE_ENE];  
	info.dihed_ene = si_vec[SI_DIHEDRAL_ENE];
	
	info.vdw_ene    = si_vec[SI_VDW_ENE];
	info.vdw_ene_14 = si_vec[SI_VDW_14_ENE];    
	info.vdw_ene_nb = si_vec[SI_VDW_ENE] - si_vec[SI_VDW_14_ENE];  

	info.electr_ene = si_vec[SI_ELECT_ENE];  
	info.electr_ene_14 = si_vec[SI_ELECT_14_ENE];  
	info.electr_ene_nb = si_vec[SI_ELECT_ENE] - si_vec[SI_ELECT_14_ENE];   
	info.polar_ene      = si_vec[SI_POLAR];
	info.polar_dip_iter = si_vec[SI_DIPITER];
	info.polar_dip_rms  = si_vec[SI_DIPRMS];

	if(p_amber_model->using_pme_potential)
	{
		info.hbond_ene = si_vec[SI_HBOND_ENE];
	}
	else
	{
		info.gb_ene = si_vec[SI_HBOND_ENE];
	}

	info.constraints_ene = si_vec[SI_RESTRAINT_ENE]; 

	info.kin_ene_com   = si_vec[SI_TOT_EKCMT];
	info.kin_ene_com_x = si_vec[SI_EKCMT_0];
	info.kin_ene_com_y = si_vec[SI_EKCMT_1];
	info.kin_ene_com_z = si_vec[SI_EKCMT_2];

	info.kin_ene_plus_half_dt  = si_vec[SI_KIN_ENE_PLUS_HALF_DT];  
	info.kin_ene_minus_half_dt = si_vec[SI_KIN_ENE_MINUS_HALF_DT];
	info.kin_ene_pbs           = si_vec[SI_KIN_ENE_PBS];

	info.kin_ene_solute  = si_vec[SI_SOLUTE_KIN_ENE];
	info.kin_ene_solvent = si_vec[SI_SOLVENT_KIN_ENE];

	info.tot_energy = info.pot_ene + info.kin_ene;

	info.temp        = si_vec[SI_TEMP];
	info.temp_solute = si_vec[SI_TEMP_SOLUTE]; 
	info.temp_solv   = si_vec[SI_TEMP_SOLVENT];

	info.virial_tot = si_vec[SI_TOT_VIRIAL];
	info.virial_x   = si_vec[SI_VIR_0];
	info.virial_y   = si_vec[SI_VIR_1];
	info.virial_z   = si_vec[SI_VIR_2]; 

	info.volume  = si_vec[SI_VOLUME];
	info.density = si_vec[SI_DENSITY];  

	info.dv_dlambda = si_vec[SI_DV_DLAMBDA];

//	info.av_perm_moment = 0.0;
//	info.av_ind_moment  = 0.0;
//	info.av_tot_moment  = 0.0;
	
	info.pme_err_est = si_vec[SI_PME_ERR_EST];
//	info.rms_ene = 0.0;    
//	info.grad_ene_max = 0.0;   
}

int MMDriverAmber::InitSimulationsStep2()
{ 
	ResetStackLimits();
//	open_iunit_debug_(&pApp->mpi_driver->myrank); // FOR MPI_DEBUG

//	PrintLog("MMDriverAmber::InitSimulationsStep2() pt 1 \n");

	if ( master ) SetupMasterNode();

	if( numtasks > 1)
	{
		p_mm_mod->BcastCtrlParams(driver_mpi_comm);
		p_mm_model->Bcast(driver_mpi_comm);
		if(mytaskid > 0) SaveModelToFortran();
		ResizeCrdVelFrcArrays(p_amber_model->natom);
		BcastCrd(driver_mpi_comm);
		BcastVel(driver_mpi_comm);
		BcastPBox(driver_mpi_comm);
	}

	CalcNumDegFreedom();
	FC_FUNC_MODULE(alltasks_setup_mod,alltasks_setup)(atm_crd.v(),pbc_box.v(),atm_vel.v(),p_amber_model->atm_igroup.v(), 
		                                  &p_mm_mod->run_type.value(), &p_mm_mod->period_bcond.value(),
		                                  gbl_atm_owner_map.v(), gbl_my_atm_lst.v(),&my_atm_cnt,p_amber_model->atm_mass_inv.v());
	
	//  Initialize random number generator at same point in all processors. Then, if
	//  random initial velocities are needed, generate them in all processors. 
	//  In general, we must be careful to generate the same sequence of random
	//  numbers in all processors.

	FC_FUNC_MODULE(random_mod,amrset)(&p_mm_mod->random_seed);

	if (p_mm_mod->init_read_coord == p_mm_mod->init_read_coord.READ_X_BIN || 
		p_mm_mod->init_read_coord == p_mm_mod->init_read_coord.READ_X_FORM ) 
	{
		FC_FUNC_MODULE(dynamics_mod,all_atom_setvel)(&p_amber_model->natom, atm_vel.v(), p_amber_model->atm_mass_inv.v(), &p_mm_mod->init_temp);
	}

	if (p_amber_model->ibelly > 0) ZeroVelFrozenAtoms();
	
	if ( p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN || master ) SetupShakePars(); // Shake params setup

	if ( p_amber_model->using_pme_potential )
	{
		SetupPME();
	}
	else if ( p_amber_model->using_gb_potential )
	{
		SetupGB();
	}

// Validate the cntrl namelist in the mdin file:
	if( master) ValidateCntrParams();

// Set parameters to insure update of variables
// computed at first call of Force and energy Function

	FC_FUNC_MODULE(pme_force_mod,setup_not_done) = 1; 
	FC_FUNC_MODULE(pme_force_mod,irespa) = 0;

	p_tm->EndSetupTimers();

//	PrintLog("MMDriverAmber::InitSimulationsStep2() pt end \n");
	return TRUE;
}

void MMDriverAmber::CalcForceAndEne(int new_list, int ncalls)
{
	if (p_amber_model->using_pme_potential) 
	{
		if( !p_mm_model->IsAmoebaFF() )
		{ 
			FC_FUNC_MODULE(loadbal_mod,do_load_balancing)(&new_list, &p_amber_model->natom, atm_crd.v(), atm_frc.v(),
				                          gbl_atm_owner_map.v(), gbl_my_atm_lst.v(), &my_atm_cnt, p_amber_model->atm_mass.v() );
		}
//		PrintLog("MMDriverAmber::CalcForceAndEne() pt 1 \n");
		PMEForce(p_amber_model->natom, atm_crd, atm_vel, atm_frc, new_list, sys_info);
	}
	else if (p_amber_model->using_gb_potential)
	{
		GBForce(p_amber_model->natom,atm_crd,atm_frc,sys_info,ncalls);
	}

//	PrintLog(" MMDriverAmber::CalcForceAndEne() pt 1  sys_info[SI_BOND_ENE] = %12.6f \n",
//		         sys_info[SI_BOND_ENE] );

	if(p_mm_mod->run_ti) p_mm_mod->p_ti_mod->CollectForceAndEneTI(sys_info);

//	PrintLog(" MMDriverAmber::CalcForceAndEne() pt 2  sys_info[SI_BOND_ENE] = %12.6f \n",
//		         sys_info[SI_BOND_ENE] );

	if( p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN && ntp > 0 )
	{
		sys_info[SI_VOLUME] = uc_volume;
		sys_info[SI_DENSITY] = p_amber_model->tmass / (0.602204 * sys_info[SI_VOLUME]);
		sys_info[SI_TOT_EKCMT]  = sys_info[SI_EKCMT_0] + sys_info[SI_EKCMT_1] + sys_info[SI_EKCMT_2];
		sys_info[SI_TOT_VIRIAL] = sys_info[SI_VIR_0]   + sys_info[SI_VIR_1]   + sys_info[SI_VIR_2];
//		CalcPressure(sys_info);
	}

//	PrintLog(" MMDriverAmber::CalcForceAndEne() pot_ene = %12.6f  vdw_ene= %12.6f  elect_ene = %12.6f \n", 
//		      sys_info[SI_POT_ENE], sys_info[SI_VDW_ENE], sys_info[SI_ELECT_ENE]);
//	PrintLog(" MMDriverAmber::CalcForceAndEne() vir[0-2] = %12.6f  %12.6f %12.6f \n", sys_info[SI_VIR_0],sys_info[SI_VIR_1],sys_info[SI_VIR_2]);

// May be we should broadcast pressure for TI caclulations
// also may be we have to synchronize coordinate from time to time in TI runs

//		if( master)
//		{
//			PrintLogMDOUT(" nstep = %5d  \n", ncalls );
//			PrintLogMDOUT("             1st at crd = %16.9f %16.9f %16.9f \n",atm_crd[0],atm_crd[1],atm_crd[2]);
//			PrintLogMDOUT("                    box = %16.9f %16.9f %16.9f \n",pbc_box[0],pbc_box[1],pbc_box[2]);
//			PrintLogMDOUT("             1st at vel = %16.9f %16.9f %16.9f \n",atm_vel[0],atm_vel[1],atm_vel[2]);
//			PrintLogMDOUT("             1st at frc = %16.9f %16.9f %16.9f \n",atm_frc[0],atm_frc[1],atm_frc[2]);
//			PrintLogMDOUT("                 virial = %16.9f %16.9f %16.9f \n",sys_info[SI_VIR_0],sys_info[SI_VIR_1],sys_info[SI_VIR_2]);
//			PrintLogMDOUT("                  ekcmt = %16.9f %16.9f %16.9f \n",sys_info[SI_EKCMT_0],sys_info[SI_EKCMT_1],sys_info[SI_EKCMT_2]);
//			PrintLogMDOUT("  \n");
//		}
//		std::string str;
//		SetMMInfo((*p_mm_mod->p_mm_info),sys_info);
//		p_mm_mod->PrintEneStrAccurate((*p_mm_mod->p_mm_info),str);
//		PrintLogMDOUT("\n%s\n",str.c_str());
}

void MMDriverAmber::CalcKinEne(int only_cur_vel)
{
	// Get the kinetic energy, either for printing or for Berendsen.

	int j,atm_lst_idx;

	double half_dtx = p_mm_mod->md_time_step * 20.455 * 0.5;
	double gammai = p_mm_mod->langevin_dump_const / 20.455;
	double c_ave    = 1.0 + gammai * half_dtx;

    HaVec_double reduce_buf_in(3);
	HaVec_double reduce_buf_out(3);

	double eke = 0.0;
	double ekph = 0.0;
	double ekpbs = 0.0;  

	if( only_cur_vel )
	{
		if (p_mm_mod->temp_control_method != p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN)
		{
			for( atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++ )
			{
				j = gbl_my_atm_lst[atm_lst_idx] - 1;

				eke += p_amber_model->atm_mass[j] * ( atm_vel[3*j]  * atm_vel[3*j] + 
					                                  atm_vel[3*j+1]* atm_vel[3*j+1] +
									                  atm_vel[3*j+2]* atm_vel[3*j+2] );
			} 
		}
		else
		{
			for( atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++ )
			{
				j = gbl_my_atm_lst[atm_lst_idx] - 1;

				eke += p_amber_model->atm_mass[j] * 0.25 * c_ave * 
					                          ( atm_vel[3*j]  * atm_vel[3*j] + 
					                            atm_vel[3*j+1]* atm_vel[3*j+1] + 
					                            atm_vel[3*j+2]* atm_vel[3*j+2] );
			}
		} // (p_mm_mod->temp_control_method != HaMolMechMod::CONST_TEMP_LANGEVIN)
		ekpbs = eke;
		ekph  = eke;
	}
	else
	{
		if (p_mm_mod->temp_control_method != p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN)
		{
			for( atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++ )
			{
				j = gbl_my_atm_lst[atm_lst_idx] - 1;

				eke += p_amber_model->atm_mass[j] * ( (atm_vel[3*j] + atm_last_vel[3*j])*(atm_vel[3*j] + atm_last_vel[3*j]) + 
					                                  (atm_vel[3*j+1] + atm_last_vel[3*j+1])*(atm_vel[3*j+1] + atm_last_vel[3*j+1]) +
					                                  (atm_vel[3*j+2] + atm_last_vel[3*j+2])*(atm_vel[3*j+2] + atm_last_vel[3*j+2]) );

				ekpbs += p_amber_model->atm_mass[j] * ( atm_vel[3*j]  * atm_last_vel[3*j] + 
				                                    	atm_vel[3*j+1]* atm_last_vel[3*j+1] +
					                                    atm_vel[3*j+2]* atm_last_vel[3*j+2] );

				ekph  += p_amber_model->atm_mass[j] * ( atm_vel[3*j]  * atm_vel[3*j] + 
					                                    atm_vel[3*j+1]* atm_vel[3*j+1] +
					                                    atm_vel[3*j+2]* atm_vel[3*j+2] );
			}
			eke *= 0.25; 
		}
		else
		{
			for( atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++ )
			{
				j = gbl_my_atm_lst[atm_lst_idx] - 1;

				eke += p_amber_model->atm_mass[j] * 0.25 * c_ave * 
					                  ( (atm_vel[3*j] + atm_last_vel[3*j])*(atm_vel[3*j] + atm_last_vel[3*j]) + 
					                  (atm_vel[3*j+1] + atm_last_vel[3*j+1])*(atm_vel[3*j+1] + atm_last_vel[3*j+1]) +
					                  (atm_vel[3*j+2] + atm_last_vel[3*j+2])*(atm_vel[3*j+2] + atm_last_vel[3*j+2]) );
			}
		} // (p_mm_mod->temp_control_method != CONST_TEMP_LANGEVIN)
	}

	eke   *= 0.5;  // Kinetic energy at the current point 
	ekpbs *= 0.5;
	ekph  *= 0.5;  // Kinetic energy at the point t+ 1/2 * md_time_step

	if( numtasks > 1 )
	{
		// Sum up the partial kinetic energies and complete the skin check. 

		reduce_buf_in[0] = eke;
		reduce_buf_in[1] = ekph;
		reduce_buf_in[2] = ekpbs;

		p_tm->UpdateTime(TimerAmber::RUNMD_TIME);
		int ierr = MPI_Allreduce(reduce_buf_in.v(), reduce_buf_out.v(), 3, MPI_DOUBLE,
			MPI_SUM, driver_mpi_comm);
		p_tm->UpdateTime(TimerAmber::FCVE_DIST_TIME);

		eke = reduce_buf_out[0];
		ekph = reduce_buf_out[1];
		ekpbs = reduce_buf_out[2];
	}

	sys_info[SI_KIN_ENE] = eke;
	sys_info[SI_KIN_ENE_PLUS_HALF_DT] = ekph;
	sys_info[SI_KIN_ENE_PBS] = ekpbs;
	sys_info[SI_SOLVENT_KIN_ENE] = ekpbs + sys_info[SI_POT_ENE];
	sys_info[SI_SOLUTE_KIN_ENE] = eke;	

	sys_info[SI_TEMP]          = sys_info[SI_KIN_ENE]/(R_HALF_KCAL*(p_amber_model->num_deg - ndfmin));
	sys_info[SI_TEMP_SOLUTE]   = sys_info[SI_SOLUTE_KIN_ENE]/(R_HALF_KCAL*(p_amber_model->num_deg_solute - ndfmin));
	sys_info[SI_TEMP_SOLVENT]  = (p_amber_model->num_deg_solvent > 0.00001) ? sys_info[SI_SOLVENT_KIN_ENE]/(R_HALF_KCAL*p_amber_model->num_deg_solvent) : 0.0;
}

void MMDriverAmber::IncrementVelAndCrd(int nstep, HaVec_double& crd, HaVec_double& vel, HaVec_double& frc, int& reset_velocities)
{
	double dtx      = p_mm_mod->md_time_step * 20.455;
	double dtx_inv  = 1.0/dtx;
	double half_dtx = dtx * 0.5;

	int i,j,atm_lst_idx;
	 
	// Do randomization of velocities, if needed:
    // Assign new random velocities every Vrand steps, if temp_control_method == CONST_TEMP_RANDOMIZED:

	reset_velocities = FALSE;

	if (p_mm_mod->rand_vel_freq != 0 && p_mm_mod->temp_control_method == p_mm_mod->temp_control_method.CONST_TEMP_RANDOMIZED ) 
	{
		if ( ((nstep+1) % p_mm_mod->rand_vel_freq) ==  0) reset_velocities = TRUE;
	}

	if (reset_velocities) 
	{
		if (master)
		{
			PrintLogMDOUT ("Setting new random velocities at step %d \n", nstep + 1 );
		}
		double temp0_adjust = p_mm_mod->ref_temp*(p_amber_model->num_deg - ndfmin) / p_amber_model->num_deg;
		
		FC_FUNC_MODULE(dynamics_mod,vrand_set_velocities)(&p_amber_model->natom, atm_vel.v(), p_amber_model->atm_mass_inv.v(), 
			                                  &temp0_adjust, gbl_atm_owner_map.v() );
		if (p_amber_model->ibelly > 0) ZeroVelFrozenAtoms();

      // At this point in the code, the velocities lag the positions
      // by half a timestep.  If we intend for the velocities to be drawn
      // from a Maxwell distribution at the timepoint where the positions and
      // velocities are synchronized, we have to correct these newly
      // redrawn velocities by backing them up half a step using the
      // current force.
      // Note that this fix only works for Newtonian dynamics.

		if (p_mm_mod->temp_control_method != p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN)
		{
			for(atm_lst_idx = 0;  atm_lst_idx < my_atm_cnt; atm_lst_idx++)
			{
				j = gbl_my_atm_lst[atm_lst_idx]-1;

				double wfac = p_amber_model->atm_mass_inv[j] * half_dtx;
				vel[3*j]   -= frc[3*j]*wfac;
				vel[3*j+1] -= frc[3*j+1]*wfac;
				vel[3*j+2] -= frc[3*j+2]*wfac;
			}
		}
	}  // (reset_velocities)

// Do the velocity update:
   
    if (p_mm_mod->temp_control_method != p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN) 
	{
        // ---Newtonian dynamics:
		for(atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++)
		{
			j = gbl_my_atm_lst[atm_lst_idx] - 1;
			double wfac = p_amber_model->atm_mass_inv[j] * dtx;
			vel[3*j  ] += frc[3*j]  *wfac;
			vel[3*j+1] += frc[3*j+1]*wfac;
			vel[3*j+2] += frc[3*j+2]*wfac;
		}
	}
    else  //  doing langevin dynamics
	{
      // Simple model for Langevin dynamics, basically taken from
      // Loncharich, Brooks and Pastor, Biopolymers 32:523-535 (1992),
      // Eq. 11.  (Note that the first term on the rhs of Eq. 11b
      // should not be there.)
      
      FC_FUNC_MODULE(dynamics_mod,langevin_setvel)(&p_amber_model->natom, vel.v(), frc.v(), p_amber_model->atm_mass.v(), p_amber_model->atm_mass_inv.v(),
									   &p_mm_mod->md_time_step, &p_mm_mod->ref_temp, &p_mm_mod->langevin_dump_const, gbl_atm_owner_map.v() );

	}    // (p_mm_mod->temp_control_method != CONST_TEMP_LANGEVIN)

// Consider vlimit:

	if (p_mm_mod->vel_limit > 1.0e-7) 
	{
      // We here provide code that is most efficient if p_mm_mod->vel_limit is not exceeded.

		int vlimit_exceeded = FALSE;

		for(atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++)
		{
			j = gbl_my_atm_lst[atm_lst_idx] - 1;

			if (fabs(vel[3*j]) > p_mm_mod->vel_limit || fabs(vel[3*j+1]) > p_mm_mod->vel_limit ||
				fabs(vel[3*j+2])  > p_mm_mod->vel_limit) 
			{
				vlimit_exceeded = TRUE;
				break;
			}
		}

		if (vlimit_exceeded) 
		{
			double vmax = 0.0;
			for(atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++)
			{
				j = gbl_my_atm_lst[atm_lst_idx] - 1;

				for(i=0; i < 3; i++)
				{
					vmax = MaxFun(vmax, fabs(vel[3*j+i]));
					vel[3*j+i] = MinFun(fabs(vel[3*j+i]), p_mm_mod->vel_limit)* ( vel[3*j+i] > 0 ? 1.0 : -1.0 );
				}
			}
        // Only violations on the master node are actually reported
        // to avoid both MPI communication and non-master writes.
			if(master) PrintLogMDOUT ("vel_limit exceeded for step %d = ; vmax = %16.9f \n",nstep, vmax);
		}

	} // (p_mm_mod->vel_limit > 1.0e-7)

// Update the positions, putting the "old" positions into frc():

//	PrintLog(" MD Step = %d \n", nstep);
//	for( atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++ )
//	{
//
//	}

	for(atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++)
	{
		i = gbl_my_atm_lst[atm_lst_idx] - 1;

		frc[3*i]   = crd[3*i];
		frc[3*i+1] = crd[3*i+1];
		frc[3*i+2] = crd[3*i+2];
		
		crd[3*i]   += vel[3*i] * dtx;
		crd[3*i+1] += vel[3*i+1] * dtx;
		crd[3*i+2] += vel[3*i+2] * dtx;
	}

// If shake is being used, update new positions to fix bond lengths:

    if (ntc != 1) 
	{
		p_tm->UpdateTime(TimerAmber::RUNMD_TIME);
		FC_FUNC_MODULE(shake_mod,shake)(frc.v(), crd.v(),&p_mm_mod->period_bcond.value(),
			                p_amber_model->atm_igroup.v(),pbc_box.v(),p_amber_model->atm_mass_inv.v(),&p_mm_mod->shake_tol);
		p_tm->UpdateTime(TimerAmber::SHAKE_TIME);

      // Re-estimate velocities from differences in positions:

		for(atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++)
		{
			j = gbl_my_atm_lst[atm_lst_idx] - 1;
			vel[3*j]   = (crd[3*j]   - frc[3*j])  *dtx_inv;
			vel[3*j+1] = (crd[3*j+1] - frc[3*j+1])*dtx_inv;
			vel[3*j+2] = (crd[3*j+2] - frc[3*j+2])*dtx_inv;
		}
	}
	if(numtasks > 1 )
	{
		all_crds_valid = FALSE;
		all_vels_valid = FALSE;
	}
}

void MMDriverAmber::CalcPressure(HaVec_double& sys_info)
{
	const double pconv = 1.6604345e+4 * 4.184; // factor to convert pressure from kcal/mole to bar.

//    PrintLog(" kin_ene = %12.6f   vir_1 = %12.6f  vir_2 = %12.6f  vir_3 = %12.6f \n",
//		      sys_info[SI_KIN_ENE],sys_info[SI_VIR_0],sys_info[SI_VIR_1],sys_info[SI_VIR_2] );

// 	sys_info[SI_PRESS_0] = pconv*( sys_info[SI_KIN_ENE]*2.0/3.0 - sys_info[SI_VIR_0] )/sys_info[SI_VOLUME];
//  sys_info[SI_PRESS_1] = pconv*( sys_info[SI_KIN_ENE]*2.0/3.0 - sys_info[SI_VIR_1] )/sys_info[SI_VOLUME];
//  sys_info[SI_PRESS_2] = pconv*( sys_info[SI_KIN_ENE]*2.0/3.0 - sys_info[SI_VIR_2] )/sys_info[SI_VOLUME];

	if( p_mm_model->IsAmoebaFF() )
	{
	    sys_info[SI_PRESS_0] = pconv*( 2.0*sys_info[SI_EKCMT_0] - sys_info[SI_VIR_0] )/sys_info[SI_VOLUME];
		sys_info[SI_PRESS_1] = pconv*( 2.0*sys_info[SI_EKCMT_1] - sys_info[SI_VIR_1] )/sys_info[SI_VOLUME];
		sys_info[SI_PRESS_2] = pconv*( 2.0*sys_info[SI_EKCMT_2] - sys_info[SI_VIR_2] )/sys_info[SI_VOLUME];
	}
	else
	{
	    sys_info[SI_PRESS_0] = 2.0*pconv*( sys_info[SI_EKCMT_0] - sys_info[SI_VIR_0] )/sys_info[SI_VOLUME];
		sys_info[SI_PRESS_1] = 2.0*pconv*( sys_info[SI_EKCMT_1] - sys_info[SI_VIR_1] )/sys_info[SI_VOLUME];
		sys_info[SI_PRESS_2] = 2.0*pconv*( sys_info[SI_EKCMT_2] - sys_info[SI_VIR_2] )/sys_info[SI_VOLUME];
	}

	sys_info[SI_TOT_PRESS] = sys_info[SI_PRESS_0] + sys_info[SI_PRESS_1] + sys_info[SI_PRESS_2];
	sys_info[SI_TOT_PRESS] = sys_info[SI_TOT_PRESS] / 3.0;	
}

void MMDriverAmber::ScaleCoordConstPress(HaVec_double& crd, HaVec_double& box, HaVec_double& sys_info)
{
	HaVec_double rmu(3);

	double dtcp = p_mm_mod->compressibility * 1.0e-06 * p_mm_mod->md_time_step / p_mm_mod->press_relax_time;
	const double one_third = 1.0/3.0;

// Here is modification to allowed semiisotropic pressure coupling
// Use ntp 3 for semiisotropic
// Use ntp 4 for only z scaling

	if (p_mm_mod->pressure_reg_method == p_mm_mod->pressure_reg_method.ISOTROP_CRD_SCALING ) 
	{
      rmu[0] = pow((1.0 - dtcp * (p_mm_mod->ref_pressure - sys_info[SI_TOT_PRESS])),one_third);
      rmu[1] = rmu[0];
      rmu[2] = rmu[0];
	}
    else if (p_mm_mod->pressure_reg_method == p_mm_mod->pressure_reg_method.ANISOTROP_CRD_SCALING ) 
	{
		rmu[0] = pow((1.0 - dtcp * (p_mm_mod->ref_pressure - sys_info[SI_PRESS_0])),one_third);
		rmu[1] = pow((1.0 - dtcp * (p_mm_mod->ref_pressure - sys_info[SI_PRESS_1])),one_third);
		rmu[2] = pow((1.0 - dtcp * (p_mm_mod->ref_pressure - sys_info[SI_PRESS_2])),one_third);
	}
	else if (p_mm_mod->pressure_reg_method == p_mm_mod->pressure_reg_method.CRD_SCALING_XY_AND_Z ) 
	{
      rmu[0] = pow((1.0 - dtcp * (p_mm_mod->ref_pressure - 0.5*(sys_info[SI_PRESS_0]+sys_info[SI_PRESS_1]))),one_third);
      rmu[1] = rmu[0];
      rmu[2] = pow((1.0 - dtcp * (p_mm_mod->ref_pressure - sys_info[SI_PRESS_2])),one_third);
	}
	else if (p_mm_mod->pressure_reg_method == p_mm_mod->pressure_reg_method.CRD_SCALING_ONLY_Z ) 
	{
      rmu[0] = 1.0;
      rmu[1] = 1.0;
      rmu[2] = pow((1.0 - dtcp * (p_mm_mod->ref_pressure - sys_info[SI_PRESS_2])),one_third);
	}
	else if (p_mm_mod->pressure_reg_method == p_mm_mod->pressure_reg_method.CRD_SCALING_XZ_AND_Y ) 
	{
      rmu[0] = pow((1.0 - dtcp * (p_mm_mod->ref_pressure - 0.5*(sys_info[SI_PRESS_0]+sys_info[SI_PRESS_2]))),one_third);
      rmu[1] = pow((1.0 - dtcp * (p_mm_mod->ref_pressure - sys_info[SI_PRESS_1])),one_third); 
      rmu[2] = rmu[0];
	}
	else if (p_mm_mod->pressure_reg_method == p_mm_mod->pressure_reg_method.CRD_SCALING_YZ_AND_X ) 
	{
      rmu[0] = pow((1.0 - dtcp * (p_mm_mod->ref_pressure - sys_info[SI_PRESS_0])),one_third);
      rmu[1] = pow((1.0 - dtcp * (p_mm_mod->ref_pressure - 0.5*(sys_info[SI_PRESS_1]+sys_info[SI_PRESS_2]))),one_third); 
      rmu[2] = rmu[0];
	}
//    else if (ntp .gt. 1) then
//      rmu(1) = (1.d0 - dtcp * (p_mm_mod->ref_pressure - press(1)))**one_third
//      rmu(2) = (1.d0 - dtcp * (p_mm_mod->ref_pressure - press(2)))**one_third
//      rmu(3) = (1.d0 - dtcp * (p_mm_mod->ref_pressure - press(3)))**one_third
//    end if

//	if( master )
//	{
//		PrintLogMDOUT("                press = %16.9f %16.9f %16.9f \n",sys_info[SI_PRESS_0],sys_info[SI_PRESS_1],sys_info[SI_PRESS_2]);
//		PrintLogMDOUT("                  rmu = %16.9f %16.9f %16.9f \n",rmu[0],rmu[1],rmu[2]);
//	}

    if (ntp > 0) 
	{
	   pbc_box[0] *= rmu[0];
	   pbc_box[1] *= rmu[1];
	   pbc_box[2] *= rmu[2];

      // WARNING!!   This is not correct for non-orthogonal boxes if NTP > 1
      // (i.e. non-isotropic scaling).  Currently general cell updates which
      // allow cell angles to change are not implemented.  The virial tensor
      // computed for ewald is the general Nose Klein; however the cell response
      // needs a more general treatment.

		FC_FUNC_MODULE(pbc_mod,pressure_scale_pbc_data)(pbc_box.v(),rmu.v());
		GetPBoxDataFromFortran();
		
		double tmp = p_amber_model->vdw_cutoff + p_mm_model->skin_nb;
		FC_FUNC_MODULE(cit_mod,set_cit_tbl_dims)(pbc_box.v(),&tmp);

        pressure_scale_crds_proxy_(atm_crd.v(), p_amber_model->atm_mass.v());

        if (p_amber_model->natc > 0) 
		{
			pressure_scale_restraint_crds_proxy_(p_amber_model->atm_xc.v(), p_amber_model->atm_mass.v());
		}
	}    //  ntp > 0
}

void MMDriverAmber::PropagateVelHalfStepBack(HaVec_double& vel, HaVec_double& frc)
{
	int i,j;
	int atm_lst_idx;

	int use_vlimit = ( p_mm_mod->vel_limit > 1.0e-7 ) ? TRUE : FALSE;
	double half_dtx = p_mm_mod->md_time_step * 20.455 * 0.5;

//	PrintLogMDOUT(" Velocities 1/2 step back: \n");
	for(atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++)
	{
		j = gbl_my_atm_lst[atm_lst_idx] - 1;
		double winf = p_amber_model->atm_mass_inv[j] * half_dtx;
		for(i = 0; i < 3; i++)
		{
			vel[3*j+i] -= frc[3*j+i] * winf;
			if (use_vlimit) vel[3*j+i] =  MinFun(fabs(vel[3*j+i]),p_mm_mod->vel_limit)*( vel[3*j+i] < 0 ? -1.0 : 1.0 ); 
		}
//		PrintLogMDOUT(" %5i %16.11f %16.11f %16.11f \n",j+1,vel[3*j],vel[3*j+1],vel[3*j+2]);
	}
//	PrintLogMDOUT(" \n");
}

void MMDriverAmber::SaveVelToLastVel(int& all_vels_valid)
{
	int atm_lst_idx;
	int j;

	if( numtasks > 1) 
	{
		if (all_vels_valid) 
		{
			atm_last_vel = atm_vel;
		}
		else
		{
			for(atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++)
			{
				j = gbl_my_atm_lst[atm_lst_idx] - 1;
				atm_last_vel[3*j]   = atm_vel[3*j];
				atm_last_vel[3*j+1] = atm_vel[3*j+1];
				atm_last_vel[3*j+2] = atm_vel[3*j+2];
			}
		}
	}
	else
	{
		atm_last_vel = atm_vel;
	}	
}

void MMDriverAmber::InitPBC()
{
	double cut_vdw_with_skin = p_amber_model->vdw_cutoff + p_mm_model->skin_nb;

	FC_FUNC_MODULE(pbc_mod,init_pbc)(&pbc_box[0], &pbc_box[1], &pbc_box[2],
			   &pbc_alpha, &pbc_beta, &pbc_gamma, &cut_vdw_with_skin);

	GetPBoxDataFromFortran();
}

//! Get the degrees of freedom for the solute and solvent 
//! Then correct the solute value for the minimum degrees of freedom in the
//! system (due to removal of translational or rotational motion, etc.).
void MMDriverAmber::CalcNumDegFreedom()
{	
// ibelsl = number of moving atoms in solute, if ibelly > 0.
// rstssl = number of degrees of freedom in the solute lost to shake.
// rstssv = number of degrees of freedom in the solvent lost to shake.
//         (rstssl and rstssv omit degrees of freedom already lost in the
//          belly case because the bond corresponds to two non-moving atoms)

// count up degrees of freedom in the belly, if any.

  int ibelsl = 0;
  int ibelsv = 0;

  int nsolut = p_amber_model->natom; // change! assume all atoms is solute 

  int i;
  if (p_amber_model->ibelly > 0) 
  {
    for( i = 0; i < p_amber_model->natom; i++) 
	{
      if (p_amber_model->atm_igroup[i] > 0) 
	  {
        if (i < nsolut) 
		{
          ibelsl = ibelsl + 1;
		}
        else
		{
          ibelsv = ibelsv + 1;
		}
	  }
	}
  }

// Now, if shake is one, loop over the appropriate bond. Add up all bonds
// in the solute/solvent which cannot move because of shake. 
// For a belly run, do not count bonds which already cannot move because both
// atoms are part of the non-moving section.
// If a bond is 1/2 in solute, 1/2 in solvent, assign 1/2 to both counters.

  double rstssl = 0.0;
  double rstssv = 0.0;

  int ib,jb;

// bonds to H

  if (ntc >= 2) 
  {
	  for(i = 0; i < p_amber_model->nbonh; i++)
	  {
		  ib = p_amber_model->gbl_bond[3*i]   - 1;
		  jb = p_amber_model->gbl_bond[3*i+1] - 1;
		  if (p_amber_model->ibelly <= 0) 
		  {
			  if (ib < nsolut && jb < nsolut) 
			  {
				  rstssl = rstssl + 1.0;
			  }
			  else if (ib >= nsolut && jb >= nsolut) 
			  {
				  rstssv = rstssv + 1.0;
			  }
			  else
			  {
				  rstssl = rstssl + 0.5;
				  rstssv = rstssv + 0.5;
			  }
		  }
		  else if (p_amber_model->atm_igroup[ib] > 0 || p_amber_model->atm_igroup[jb] > 0) 
		  {
			  if (ib < nsolut && jb < nsolut) 
			  {
				  rstssl = rstssl + 1.0;
			  }
			  else if (ib >= nsolut && jb > nsolut) 
			  {
				  rstssv = rstssv + 1.0;
			  }
			  else
			  {
				  rstssl = rstssl + 0.5;
				  rstssv = rstssv + 0.5;
			  }
		  }
	  }
  }

  // bonds to heavy atoms

  if (ntc == 3) 
  {
	  for(i = 0; i < p_amber_model->nbona; i++)
	  {
		  ib = p_amber_model->gbl_bond[3*(i+p_amber_model->nbonh)] - 1;
		  jb = p_amber_model->gbl_bond[3*(i+p_amber_model->nbonh)+1] - 1;
		  if (p_amber_model->ibelly <= 0) 
		  {
			  if (ib < nsolut && jb < nsolut) 
			  {
				  rstssl = rstssl + 1.0;
			  }
			  else if (ib >= nsolut && jb >= nsolut) 
			  {
				  rstssv = rstssv + 1.0;
			  }
			  else
			  {
				  rstssl = rstssl + 0.5;
				  rstssv = rstssv + 0.5;
			  }
		  }
		  else if (p_amber_model->atm_igroup[ib] > 0 || p_amber_model->atm_igroup[jb] > 0) 
		  {
			  if (ib < nsolut && jb < nsolut) 
			  {
				  rstssl = rstssl + 1.0;
			  }
			  else if (ib  >= nsolut && jb >= nsolut) 
			  {
				  rstssv = rstssv + 1.0;
			  }
			  else
			  {
				  rstssl = rstssl + 0.5;
				  rstssv = rstssv + 0.5;
			  }
		  }
	  }
  }

  // Now determine num_deg_solute (net number of degrees of freedom for solute) and
  //               um_deg_solvent (net number of degrees of freedom for solvent).
  if (p_amber_model->ibelly <= 0) 
  {
	  p_amber_model->num_deg_solute   = 3*nsolut - rstssl;
	  p_amber_model->num_deg_solvent  = 3*(p_amber_model->natom - nsolut) - rstssv;
  }
  else
  {
	  p_amber_model->num_deg_solute  = 3*ibelsl - rstssl;
	  p_amber_model->num_deg_solvent = 3*ibelsv - rstssv;
  }
  p_amber_model->num_deg = p_amber_model->num_deg_solute + p_amber_model->num_deg_solvent;
}



void MMDriverAmber::CalcCenMassVel(HaVec_double& crd, HaVec_double& vel, 
								   HaVec_double& xcm, HaVec_double& vcm, HaVec_double& ocm)
{
	double tmassinv = 1.0/p_amber_model->tmass; // total mass inverse for all atoms.

	double aamass;
	double comvel;
	double det;
	int    i, j, m, n;
	int    i0;
	Vec3D acm;
	HaMat_double tcm(3, 3);
	double xx, xy, xz, yy, yz, zz;
	double x1, x2, x3;
	double ekcm, ekrot;

	const double crit = 1.0e-06;

	i0 = 3 * p_amber_model->natom;

	// calculate the center of mass coordinates:

	xcm.resize(3);
	xcm = 0.0;

	i = 0;
	for(j = 0; j < p_amber_model->natom; j++)
	{
		aamass = p_amber_model->atm_mass[j];
		xcm[0] += crd[i]   * aamass;
		xcm[1] += crd[i+1] * aamass;
		xcm[2] += crd[i+2] * aamass;
		i = i + 3;
	}

	xcm[0] *= tmassinv;
	xcm[1] *= tmassinv;
	xcm[2] *= tmassinv;

	// calculate velocity and translational kinetic energy of the center of mass:

	ekcm = 0.0;
	vcm.resize(3);
	vcm = 0.0;

	i = 0;
	for(j = 0; j < p_amber_model->natom; j++)
	{
		aamass = p_amber_model->atm_mass[j];
		vcm[0] += vel[i]   * aamass;
		vcm[1] += vel[i+1] * aamass;
		vcm[2] += vel[i+2] * aamass;
		i = i + 3;
	}

	for(i = 0; i < 3; i++)
	{
		vcm[i] *= tmassinv;
		ekcm += vcm[i] * vcm[i];
	}
	ekcm *= p_amber_model->tmass * 0.5;
	comvel = sqrt(vcm[0] * vcm[0] + vcm[1] * vcm[1] + vcm[2] * vcm[2]);

	// calculate the angular momentum about the center of mass:

	acm[0] = 0.0;
	acm[1] = 0.0;
	acm[2] = 0.0;

	i = 0;
	for(j = 0; j < p_amber_model->natom; j++)
	{
		aamass = p_amber_model->atm_mass[j];
		acm[0] += (crd[i+1] * vel[i+2] - crd[i+2] * vel[i+1]) * aamass;
		acm[1] += (crd[i+2] * vel[i]   - crd[i]   * vel[i+2]) * aamass;
		acm[2] += (crd[i]   * vel[i+1] - crd[i+1] * vel[i]  ) * aamass;
		i = i + 3;
	}

	acm[0] -= (xcm[1] * vcm[2] - xcm[2] * vcm[1]) * p_amber_model->tmass;
	acm[1] -= (xcm[2] * vcm[0] - xcm[0] * vcm[2]) * p_amber_model->tmass;
	acm[2] -= (xcm[0] * vcm[1] - xcm[1] * vcm[0]) * p_amber_model->tmass;

	// calculate the inertia tensor:

	xx = 0.0;
	xy = 0.0;
	xz = 0.0;
	yy = 0.0;
	yz = 0.0;
	zz = 0.0;

	i = 0;
	for(j = 0; j < p_amber_model->natom; j++)
	{
		x1 = crd[i]   - xcm[0];
		x2 = crd[i+1] - xcm[1];
		x3 = crd[i+2] - xcm[2];
		aamass = p_amber_model->atm_mass[j];
		xx += x1 * x1 * aamass;
		xy += x1 * x2 * aamass;
		xz += x1 * x3 * aamass;
		yy += x2 * x2 * aamass;
		yz += x2 * x3 * aamass;
		zz += x3 * x3 * aamass;
		i = i + 3;
	}

	tcm.r1(1,1) = yy + zz;
	tcm.r1(2,1) = -xy;
	tcm.r1(3,1) = -xz;
	tcm.r1(1,2) = -xy;
	tcm.r1(2,2) = xx + zz;
	tcm.r1(3,2) = -yz;
	tcm.r1(1,3) = -xz;
	tcm.r1(2,3) = -yz;
	tcm.r1(3,3) = xx + yy;

	// invert the inertia tensor:

	int ires = HaMat_double::mat_inverse(tcm);

	if (!ires) 
	{
		std::strstream os; 
		os << "Error in MMDriverAmber::CalcCenMassVel() \n";
		os << "Zero deteminant of inertia tensor matrix\n";
		throw std::runtime_error(os.str());
	}

	// calculate the angular velocity about the center of mass and the rotational
	// kinetic energy:

	ekrot = 0.0;
	for(m =0; m < 3; m++)
	{
		ocm[m] = 0.0;
		for(n = 0; n < 3; n++)
		{
			ocm[m] += tcm.r0(m,n) * acm[n];
		}
		ekrot = ekrot + ocm[m] * acm[m];
	}
	ekrot = ekrot * 0.5;

	if (master)
	{
		PrintLogMDOUT("  KE Trans = %11.4f  KE Rot = %11.4f   C.O.M. Vel = %12.6f \n",ekcm,ekrot,comvel);
	}
}


void MMDriverAmber::RemoveCOMVelocity()
{
	HaVec_double vcm(3);
	HaVec_double reduce_buf_in(3);
	HaVec_double reduce_buf_out(3);

	int j,k;
	int atm_lst_idx;

//	PrintLog(" MMDriverAmber::RemoveCOMVelocity() pt 1 \n");

	double tmassinv = 1.0/p_amber_model->tmass; // total mass inverse for all atoms.

	if (p_mm_mod->temp_control_method != p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN) 
	{
		vcm = 0.0;
		
		for(atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++)
		{
			j = gbl_my_atm_lst[atm_lst_idx] - 1;
			double aamass = p_amber_model->atm_mass[j];
			
			vcm[0] += aamass * atm_vel[3*j];
			vcm[1] += aamass * atm_vel[3*j+1];
			vcm[2] += aamass * atm_vel[3*j+2];
		}

		if( numtasks > 1) 
		{
			for(k = 0; k < 3; k++)
				reduce_buf_in[k] = vcm[k]; 

			p_tm->UpdateTime(TimerAmber::RUNMD_TIME);
			int ierr = MPI_Allreduce(reduce_buf_in.v(), reduce_buf_out.v(),3, MPI_DOUBLE,
				MPI_SUM, driver_mpi_comm);
			p_tm->UpdateTime(TimerAmber::FCVE_DIST_TIME);

			for(k = 0; k < 3; k++)
				vcm[k] = reduce_buf_out[k];
		} 
		for(k = 0; k < 3; k++)
			vcm[k] *= tmassinv;

		if (master) 
		{
			double velocity2 = vcm[0] * vcm[0] + vcm[1] * vcm[1] + vcm[2] * vcm[2];
			PrintLogMDOUT(" Check COM velocity, temp: %15.6f  %9.2f (Removed) \n",sqrt(velocity2),
				0.5*p_amber_model->tmass*velocity2/(R_HALF_KCAL*(p_amber_model->num_deg - ndfmin)));
		}

		for( atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++)
		{
			j = gbl_my_atm_lst[atm_lst_idx] - 1;
			for(k = 0; k < 3; k++)
				atm_vel[3*j+k] -= vcm[k];
		}
	}
}

void MMDriverAmber::RemoveCOMVelAndResetCenter(int& all_crds_valid, int& all_vels_valid, Vec3D& sys_cnt)
{
	double half_dtx = p_mm_mod->md_time_step * 20.455 * 0.5; 
	double tmassinv = 1.0/p_amber_model->tmass;  // total mass inverse for all atoms.

	HaVec_double ocm(3), xcm(3), vcm(3); // for COM velocities

	if(numtasks > 1) 
	{
		// WARNING - currently only GB code has !pmset->per_bc->IsSet(), and currently
		// all coordinates are always valid for GB.  We include the conditional
		// below more-or-less as maintenance insurance...
		if (!all_crds_valid)  // Currently always false...
		{
			FC_FUNC_MODULE(parallel_mod,gb_mpi_allgathervec)(&p_amber_model->natom, atm_crd.v());
			all_crds_valid  = TRUE;
		}
	}
	if (p_mm_mod->temp_control_method != p_mm_mod->temp_control_method.CONST_TEMP_LANGEVIN) 
	{
		if(numtasks > 1) 
		{
			// WARNING - currently GB code never updates all velocities unless
			//          forced by this scenario...
			if (!all_vels_valid) //  Currently always true...
			{
				FC_FUNC_MODULE(parallel_mod,gb_mpi_allgathervec)(&p_amber_model->natom, atm_vel.v());
				all_vels_valid = TRUE;
			}
		}
		// Nonperiodic simulation.  Remove both translation and rotation.
		// Back the coords up a half step so they correspond to the velocities,
		// temporarily storing them in frc(:,:).

		atm_frc = atm_crd - atm_vel*half_dtx;

		// Now compute COM motion and remove it; then recompute (sander
		// compatibility...).

		// NOTE - if mass can change, that has to be taken into account for
		//        tmass, tmassinv (say for TI).

		CalcCenMassVel(atm_frc,atm_vel,xcm,vcm,ocm);

		int j,m;
		int i = 0;
		for(j = 0; j < p_amber_model->natom; j++)
		{
			for(m = 0; m < 3; m++)
			{
				atm_vel[i] -= vcm[m];
				i = i + 1;
			}
		}  
	    // stop the rotation about the center of mass:
		i = 0;
		for(j = 0; j < p_amber_model->natom; j++)
		{
			double x1 = atm_frc[i]   - xcm[0];
			double x2 = atm_frc[i+1] - xcm[1];
			double x3 = atm_frc[i+2] - xcm[2];
			atm_vel[i]   -= (ocm[1] * x3 - ocm[2] * x2);
			atm_vel[i+1] -= (ocm[2] * x1 - ocm[0] * x3);
			atm_vel[i+2] -= (ocm[0] * x2 - ocm[1] * x1);
			i = i + 3;
		}
		if (master) PrintLogMDOUT("   Translational and rotational motion removed \n");

		CalcCenMassVel(atm_frc,atm_vel,xcm,vcm,ocm);
	}
	else  // (doing langevin dynamics )
	{
		// Recenter system to the original center:
		RePosition(atm_crd,sys_cnt);
	}  // (p_mm_mod->temp_control_method != CONST_TEMP_LANGEVIN)
}

void MMDriverAmber::GetPosition(HaVec_double& crd, Vec3D& cnt )
{
	cnt[0] = 0.0;
	cnt[1] = 0.0;
	cnt[2] = 0.0;

	Vec3D crd_min;
	Vec3D crd_max;

	if( p_amber_model->natom < 1 || crd.size() < 3 ) return;

	int i,j;
	for(i= 0; i < 3; i++)
	{
		crd_min[i] = crd[i];
		crd_max[i] = crd[i];
	}

	for(i = 0; i < p_amber_model->natom; i++)
	{
		for(j = 0; j < 3; j++)
		{
			if(crd_min[j] > crd[3*i + j]) crd_min[j] =  crd[3*i + j];
			if(crd_max[j] < crd[3*i + j]) crd_max[j] =  crd[3*i + j];
		}
	}

	cnt[0] = (crd_min[0] + crd_max[0]) * 0.5;
	cnt[1] = (crd_min[1] + crd_max[1]) * 0.5;
	cnt[2] = (crd_min[2] + crd_max[2]) * 0.5;
}

void MMDriverAmber::RePosition(HaVec_double& crd, Vec3D& cnt )
{
	Vec3D cnt_curr;
	GetPosition(atm_crd, cnt_curr);

	double xd = cnt[0] - cnt_curr[0];
	double yd = cnt[1] - cnt_curr[1];
	double zd = cnt[2] - cnt_curr[2];

	if (master) PrintLogMDOUT( "| RE_POSITION Moving by %10.6f %10.6f %10.6f \n", xd, yd, zd);

	int i;
	for(i = 0; i < p_amber_model->natom; i++)
	{
		crd[3*i]   += xd;
		crd[3*i+1] += yd;
		crd[3*i+2] += zd;
	}
}

void MMDriverAmber::GrdMax(HaVec_double& frc,int& iatmax, double& fdmax)
{
	fdmax = 0.0;
	iatmax = 0;
	int i;
	for(i = 0; i < frc.size(); i++)
	{
		double gi = fabs(frc[i]);
		if (gi > fdmax)
		{
			fdmax = gi;
			iatmax = i;
		}
	}
}

void MMDriverAmber::ZeroVelFrozenAtoms()
{
	int idx;
	for( idx = 0; idx < my_atm_cnt; idx++)
	{
		int i = gbl_my_atm_lst[idx];
		if (p_amber_model->atm_igroup[i] > 0) continue;
		atm_vel[3*i]   = 0.0;
		atm_vel[3*i+1] = 0.0;
		atm_vel[3*i+2] = 0.0;
	}
}

void MMDriverAmber::ZeroFrcFrozenAtoms()
{
	int idx;
	for( idx = 0; idx < my_atm_cnt; idx++)
	{
		int i = gbl_my_atm_lst[idx];
		if (p_amber_model->atm_igroup[i] > 0) continue;
		atm_frc[3*i]   = 0.0;
		atm_frc[3*i+1] = 0.0;
		atm_frc[3*i+2] = 0.0;
	}
}

void MMDriverAmber::CollectCoords(int& new_list, int& nstep, int& collect_crds, int& all_crds_valid)
{
	if( numtasks == 0 && !p_mm_mod->run_ti) return;

	p_tm->UpdateTime(TimerAmber::RUNMD_TIME);

	if( numtasks > 1) 
	{
		if (p_amber_model->using_gb_potential) 
		{
			// Generalized Born has it's own coordinate distribution scheme...
			FC_FUNC_MODULE(parallel_mod,gb_mpi_allgathervec)(&p_amber_model->natom, atm_crd.v());
			all_crds_valid = TRUE;
		}
		else if (new_list || p_mm_model->IsAmoebaFF() ) // IGOR TEMPORAL FIX
		{
			// This always gets done for a work redistribution cycle:
			AllGatherVec(atm_crd);

			// If this is a const pressure run and there are coordinate constraints,
			// we will need to send the adjusted coordinate constraints around when
			// we redistribute the work load.  Bummer.
			if (ntp > 0) 
			{
				if ( p_amber_model->natc > 0 ) 
				{
					if (FC_FUNC_MODULE(loadbal_mod,get_atm_redist_needed)() || FC_FUNC_MODULE(loadbal_mod,get_fft_slab_redist_needed)()) 
					{
						AllGatherVec(p_amber_model->atm_xc);
					}
				}
			}
			all_crds_valid = TRUE;
		}
		else // ( !using_gb_potential && !new_list ) 
		{
			FC_FUNC_MODULE(parallel_mod,distribute_crds_proxy)(&p_amber_model->natom, &my_atm_cnt, atm_crd.v());

			// It may be necessary for the master to have a complete copy of the
			// coordinates for archiving and writing the restart file.

			collect_crds = FALSE;

			if (p_mm_mod->wrt_coord_freq > 0) 
			{
				if ( (nstep + 1) % p_mm_mod->wrt_coord_freq == 0) collect_crds = TRUE;
			}
			if (p_mm_mod->wrt_constr_freq > 0) 
			{
				if ( (nstep + 1) % p_mm_mod->wrt_constr_freq == 0) collect_crds = TRUE;
			}
			if (p_mm_mod->wrt_rstrt_freq != 0) 
			{
				if ( (nstep + 1) % p_mm_mod->wrt_rstrt_freq == 0) collect_crds = TRUE;
			}
			if ( (nstep + 1) >= p_mm_mod->num_md_steps) collect_crds = TRUE;

			if (collect_crds) 
			{
				FC_FUNC_MODULE(parallel_mod,mpi_gathervec)(&p_amber_model->natom, atm_crd.v(), gbl_atm_owner_map.v(), 
					                           gbl_my_atm_lst.v(), &my_atm_cnt);
				all_crds_valid = TRUE;
			}
		} // ( using_gb_potential )
		
	} // ! (numtasks > 1)

	p_tm->UpdateTime(TimerAmber::FCVE_DIST_TIME);
}

void MMDriverAmber::CollectVelocities(int& new_list, int& nstep, int& collect_crds, int& all_vels_valid)
{
	if( numtasks > 1) 
	{
		// It may be necessary for the master to have a complete copy of the velocities
		// for archiving and writing the restart file. In any case, all processors must
		// have a complete copy of velocities if atom redistribution is going to happen
		// in the next call to force.

		p_tm->UpdateTime(TimerAmber::RUNMD_TIME);

		if (p_amber_model->using_pme_potential && new_list && 
			(FC_FUNC_MODULE(loadbal_mod,get_atm_redist_needed)() || FC_FUNC_MODULE(loadbal_mod,get_fft_slab_redist_needed)()) )
		{
			AllGatherVec(atm_vel);
			all_vels_valid = TRUE;
		}
		else
		{
			int collect_vels = FALSE;

			if (p_mm_mod->wrt_vel_freq > 0) 
			{
				if ( (nstep + 1) % p_mm_mod->wrt_vel_freq == 0) collect_vels = TRUE;
			}
			else if (collect_crds && p_mm_mod->wrt_vel_freq < 0) 
			{
				collect_vels = TRUE;
			}

			if (p_mm_mod->wrt_rstrt_freq != 0) 
			{
				if ( (nstep + 1) % p_mm_mod->wrt_rstrt_freq == 0) collect_vels = TRUE;
			}

			if ( (nstep + 1) >= p_mm_mod->num_md_steps) collect_vels = TRUE;

			if ( collect_vels ) 
			{
				if (p_amber_model->using_pme_potential)
				{
					FC_FUNC_MODULE(parallel_mod,mpi_gathervec)(&p_amber_model->natom, atm_vel.v(), gbl_atm_owner_map.v(), 
						                           gbl_my_atm_lst.v(), &my_atm_cnt);
				}
				else if (p_amber_model->using_gb_potential) 
				{
					FC_FUNC_MODULE(parallel_mod,gb_mpi_gathervec)(&p_amber_model->natom, atm_vel.v());
				}
			}
			all_vels_valid = FALSE;  // ?? should it be TRUE is collect_vels happened???
		}
		p_tm->UpdateTime(TimerAmber::FCVE_DIST_TIME);
	}  // (numtasks > 1)

//	if( p_mm_mod->run_ti )
//	{
//		if( mytaskid == 0) BcastVel(p_mm_mod->inter_model_comm);
//		MPI_Barrier(MPI_COMM_WORLD);
//		if( p_mm_mod->inter_model_comm == 1 )
//		{
//			BcastVel(driver_mpi_comm);
//		}
//	}
}


void MMDriverAmber::ScaleVelConstTemp()
{
//           --- following is from T.E. Cheatham, III and B.R. Brooks,
//              Theor. Chem. Acc. 99:279, 1998.
// Pastor, Brooks, Szabo conserved quantity for harmonic oscillator:
// Eqn. 4.7b of Mol. Phys. 65:1409-1419, 1988:

	// ekph - kinetic energy from atm_vel      - current  velocities ( t + md_time_step/2 ?)
	// ekmh - kinetic energy from atm_last_vel - previous velocities ( t - md_time_step/2 ?)

	double eke  = sys_info[SI_KIN_ENE];
	double ekmh = sys_info[SI_KIN_ENE_MINUS_HALF_DT];
	double ekph = sys_info[SI_KIN_ENE_PLUS_HALF_DT];

	int j,atm_lst_idx;
	
	double ekin0  = R_HALF_KCAL * (p_amber_model->num_deg - ndfmin) * p_mm_mod->ref_temp;  // kinetic energy corresponding to temperature to keep constant - temp0
	
	if( emulate_ext_amber )
	{
		ekin0  = (8.31441E-3 * 0.5 / 4.184 ) * (p_amber_model->num_deg - ndfmin) * p_mm_mod->ref_temp;
	}

	double scaltp = sqrt(1.0 + 2.0 * (p_mm_mod->md_time_step/p_mm_mod->temp_relax_time_solute) * (ekin0 - eke)/(ekmh + ekph));

 //   PrintLogMDOUT(" MMDriverAmber::ScaleVelConstTemp() \n");
	//PrintLogMDOUT("fac(1) = %16.11f   temp0 = %16.11f \n",R_HALF_KCAL * (p_amber_model->num_deg - ndfmin), p_mm_mod->ref_temp); 
	//PrintLogMDOUT("ekin0 = %16.11f   scaltp = %16.11f \n",ekin0, scaltp);
	//PrintLogMDOUT("eke  = %18.11f \n",eke);
	//PrintLogMDOUT("ekmh = %18.11f \n",ekmh);
	//PrintLogMDOUT("ekph = %18.11f \n",ekph);

	for( atm_lst_idx = 0; atm_lst_idx < my_atm_cnt; atm_lst_idx++)
	{
		j = gbl_my_atm_lst[atm_lst_idx] - 1;
		atm_vel[3*j]   *= scaltp;
		atm_vel[3*j+1] *= scaltp;
		atm_vel[3*j+2] *= scaltp;
	}
}

int MMDriverAmber::CheckForNewNonBondList()
{
	int new_list = FALSE;

	if (p_amber_model->using_pme_potential) 
	{
      // Now we can do a skin check to see if we will have to rebuild the
      // pairlist next time...

		if( numtasks > 1) 
		{
			FC_FUNC_MODULE(nb_pairlist_mod,check_my_atom_movement)(atm_crd.v(), gbl_atm_saved_crd.v(), pbc_box.v(), gbl_saved_box.v(), gbl_my_atm_lst.v(),&my_atm_cnt,
                                                       &p_mm_model->skin_nb, &ntp, &new_list);
		}
		else
		{
			new_list = CheckAllAtomMovement(p_amber_model->natom,atm_crd);
		}
	} // (use_pme_potential)
    // Nothing to do for Generalized Born...

	if( numtasks > 1) 
	{
		double new_list_cnt;                   
		double new_list_cnt_out;

		if (new_list)
		{
            new_list_cnt = 1.0;
		}
        else
		{
            new_list_cnt = 0.0;
		}
    
        p_tm->UpdateTime(TimerAmber::RUNMD_TIME);
		int ierr = MPI_Allreduce(&new_list_cnt, &new_list_cnt_out, 1, MPI_DOUBLE,
                             MPI_SUM, driver_mpi_comm);
        p_tm->UpdateTime(TimerAmber::FCVE_DIST_TIME);

    // Determine if any process saw an atom exceed the skin check.  We use the
    // comparison to 0.5 to handle any rounding error issues; we use
    // double precision for new_list_cnt in order to be able to piggyback the 
    // reduce.

		new_list = (new_list_cnt_out >= 0.5) ? TRUE : FALSE;

        if (p_amber_model->using_pme_potential) 
		{
            FC_FUNC_MODULE(loadbal_mod,check_new_list_limit)(&new_list);
		}
	}
	return new_list;
}

void MMDriverAmber::SetPrmTopIntFortran()
{
	HaVec_int vals(27);
	int iaxx = 0;

	vals(1) = p_amber_model->natom; 
	vals(2) = p_amber_model->ntypes; 
	vals(3) = p_amber_model->nbonh; 
	vals(4) = p_amber_model->ntheth; 
	vals(5) = p_amber_model->nphih; 
	vals(6) = p_amber_model->next; 
	vals(7) = p_amber_model->nres; 
    vals(8) = p_amber_model->nbona; 
    vals(9) = p_amber_model->ntheta; 
    vals(10) = p_amber_model->nphia; 
    vals(11) = p_amber_model->numbnd; 
    vals(12) = p_amber_model->numang; 
    vals(13) = p_amber_model->nptra;
    vals(14) = p_amber_model->nphb; 
    vals(15) = iaxx; 
    vals(16) = iaxx; 
    vals(17) = p_amber_model->nspm; 
    vals(18) = iaxx;
    vals(19) = iaxx; 
    vals(20) = p_amber_model->nttyp; 
    vals(21) = p_amber_model->bonda_idx; 
    vals(22) = p_amber_model->anglea_idx; 
    vals(23) = p_amber_model->diheda_idx; 
    vals(24) = p_amber_model->gbl_bond_allocsize; 
    vals(25) = p_amber_model->gbl_angle_allocsize; 
    vals(26) = p_amber_model->gbl_dihed_allocsize; 
    vals(27) = p_amber_model->num_dist_constr;
	
	set_prmtop_int_(vals.begin());
}


void MMDriverAmber::SetMDinCtrlIntFortran()
{
	HaVec_int vals(63);
	int ixx = 999999;

	vals(1) = ixx; 
	vals(2) = ixx;  
	vals(3) = ixx;   
	vals(4) = ixx;  
	vals(5) = ixx; 
	vals(6) = ixx;    
    vals(7) = ixx; 
    vals(8) = ixx; 
    vals(9)  = ixx;  
    vals(10) = ixx;
    vals(11) = ixx; 
    vals(12) = ixx; 
    vals(13) = ixx; 
    vals(14) = ixx; 
    vals(15) = ixx; 
    vals(16) = p_amber_model->ntf; 
    vals(17) = ixx; 
    vals(18) = nsnb; 
    vals(19) = ixx; 
    vals(20) = p_amber_model->ibelly; 
    vals(21) = ixx; 
    vals(22) = ixx; 
    vals(23) = ixx; 
    vals(24) = ixx; 
    vals(25) = ixx; 
    vals(26) = ixx; 
    vals(27) = nrespa; 
    vals(28) = ixx;
    vals(29) = ixx; 
    vals(30) = ixx; 
    vals(31) = ntp; 
    vals(32) = ntc; 
    vals(33) = jfastw; 
    vals(34) = ixx; 
    vals(35) = p_amber_model->igb; 
    vals(36) = p_amber_model->alpb; 
    vals(37) = p_amber_model->rbornstat; 
    vals(38) = p_amber_model->gbsa; 
    vals(39) = nrespai; 
    vals(40) = p_mm_model->IsAmoebaFF(); 
	vals(41) = p_amber_model->do_amoeba_valence; 
	vals(42) = p_amber_model->do_amoeba_nonbond; 
	vals(43) = p_amber_model->do_bond;
	vals(44) = p_amber_model->do_ureyb;
	vals(45) = p_amber_model->do_reg_angle; 
	vals(46) = p_amber_model->do_trig_angle; 
    vals(47) = p_amber_model->do_opbend; 
    vals(48) = p_amber_model->do_torsion; 
    vals(49) = p_amber_model->do_pi_torsion; 
	vals(50) = p_amber_model->do_strbend; 
	vals(51) = p_amber_model->do_torsion_torsion;
	vals(52) = p_amber_model->do_str_torsion;
	vals(53) = p_amber_model->do_recip;
	vals(54) = p_amber_model->do_adjust;
	vals(55) = p_amber_model->do_direct;
	vals(56) = p_amber_model->do_self;
	vals(57) = p_amber_model->do_vdw;
	vals(58) = p_amber_model->do_induced;
	vals(59) = p_amber_model->do_vdw_taper;
	vals(60) = p_amber_model->do_vdw_longrange;
	vals(61) = p_amber_model->beeman_integrator;
	vals(62) = p_mm_model->dipole_scf_iter_max;
	vals(63) = p_amber_model->amoeba_verbose;
 
	set_mdin_ctrl_int_(vals.begin());
}

void MMDriverAmber::SetMDinCtrlDblFortran()
{
	HaVec_double vals(48);

	double xx = 999999999.0;

	vals(1) = p_amber_model->dielc;
	vals(2) = p_amber_model->es_cutoff;  
	vals(3) = p_amber_model->vdw_cutoff; 
    vals(4) = p_amber_model->scnb;  
    vals(5) = p_amber_model->scee; 
    vals(6) = xx; 
    vals(7) = xx; 
    vals(8) = xx; 
    vals(9) = xx; 
    vals(10) = xx;  
    vals(11) = xx; 
    vals(12) = xx; 
    vals(13) = xx; 
    vals(14) = xx; 
    vals(15) = xx; 
    vals(16) = xx; 
    vals(17) = xx; 
    vals(18) = xx; 
    vals(19) = xx; 
    vals(20) = xx;
    vals(21) = p_amber_model->intdiel; 
    vals(22) = p_amber_model->extdiel; 
    vals(23) = p_amber_model->saltcon; 
    vals(24) = p_amber_model->rgbmax; 
    vals(25) = p_amber_model->offset;
    vals(26) = p_amber_model->surften; 
    vals(27) = p_amber_model->cut_inner; 
    vals(28) = p_amber_model->gb_cutoff; 
    vals(29) = p_amber_model->gb_alpha; 
    vals(30) = p_amber_model->gb_beta;  
    vals(31) = p_amber_model->gb_gamma; 
    vals(32) = p_amber_model->gb_fs_max; 
    vals(33) = p_amber_model->gb_kappa; 
    vals(34) = p_amber_model->gb_neckscale; 
    vals(35) = p_amber_model->arad; 
    vals(36) = p_amber_model->bbox_xmin; 
    vals(37) = p_amber_model->bbox_ymin; 
    vals(38) = p_amber_model->bbox_zmin;  
    vals(39) = p_amber_model->bbox_xmax; 
    vals(40) = p_amber_model->bbox_ymax; 
    vals(41) = p_amber_model->bbox_zmax;

	vals(42) = p_mm_mod->compressibility*1.0e-6;
	vals(43) = p_mm_model->dipole_scf_tol;
	vals(44) = p_mm_model->ee_dsum_cut;
	vals(45) = p_mm_model->ee_damped_cut;
	vals(46) = p_mm_model->sor_coefficient;
	vals(47) = p_mm_model->thole_expon_coeff; 
	vals(48) = p_mm_model->vdw_taper; 
	
	set_mdin_ctrl_dbl_(vals.begin());
}

void MMDriverAmber::SetAddIntParsFortran()
{
	HaVec_int vals(5);

	vals(1) = p_amber_model->using_pme_potential;
	vals(2) = p_amber_model->using_gb_potential;
	vals(3) = dbg_atom_redistribution;
	vals(4) = loadbal_verbose;
	vals(5) = master;

	FC_FUNC_MODULE(mdin_ewald_dat_mod,set_add_pmemd_int_pars)(vals.begin());
}

void MMDriverAmber::SetPMEParsFortran()
{
	HaVec_int vals(7);

// PME Int Params

	vals(1) = p_mm_model->pme_grid_nx;
	vals(2) = p_mm_model->pme_grid_ny;
	vals(3) = p_mm_model->pme_grid_nz;
	vals(4) = p_mm_model->pme_spline_order;
	vals(5) = p_mm_model->subtract_avg_force_flag;
	vals(6) = p_mm_model->vdw_correction_flag;
	vals(7) = p_mm_model->pme_verbose;

	FC_FUNC_MODULE(mdin_ewald_dat_mod,set_pme_int_pars_to_fort)(vals.v());

	HaVec_double dvals(5);

	dvals(1) = p_mm_model->skin_nb;
	dvals(2) = p_mm_model->pme_dsum_tol;
	dvals(3) = p_mm_model->pme_ew_coeff;
	dvals(4) = p_mm_model->pme_eedtbdns;
	dvals(5) = p_mm_model->fft_grids_per_ang;

	FC_FUNC_MODULE(mdin_ewald_dat_mod,set_pme_dbl_pars_to_fort)(dvals.v());
}

void MMDriverAmber::PrintLogMDOUT(const char* format, ... )
{
#	if defined(WITH_LIB_PMEMD)
	char buf[10000];
	char buf1[1000];
	int len;
	int ipos;
	int len1;

	va_list arg_list;
	va_start( arg_list, format );     /* Initialize variable arguments. */
	vsprintf( buf, format, arg_list);
	
	len = strlen(buf);
	if( len > 0 && buf[len-1] == '\n') len--;
	ipos = 0;
	len1 = 0;
	for(ipos = 0; ipos < len; ipos++ )
	{
		if(buf[ipos] == '\n' || ipos == (len - 1))
		{
			if( buf[ipos] != '\n')
			{
				buf1[len1] = buf[ipos];
				len1++;
			}
			if(len1 == 0)
			{
				buf1[0] = ' ';
				len1 = 1;
			}
			write_mdout_(buf1,len1);
			len1 = 0;
		}
		else
		{
			buf1[len1] = buf[ipos];
			len1++;
		}
	}
	va_end(arg_list);              /* Reset variable arguments.      */
	FlushMDOUT();
#	endif
}

void MMDriverAmber::FlushMDOUT()
{
#if defined(WITH_LIB_PMEMD)
  int current_sec;
  int current_usec;     // Dummy, not used.

  get_wall_time_(&current_sec, &current_usec);

  if (current_sec >= next_mdout_flush_sec)
  {
      next_mdout_flush_sec = current_sec + mdout_flush_interval;
	  FC_FUNC_MODULE(runfiles_mod,flush_mdout)();
  }
#endif
}

void MMDriverAmber::GBForce(int atm_cnt, HaVec_double& crd, HaVec_double& frc,
							HaVec_double& si_vec, int ncalls)
{
#	if defined(WITH_LIB_PMEMD)
	int n3 = atm_cnt*3;
	if(n3 != crd.size()) 
	{
		PrintLog("Error in MMDriverAmber::GBForce() \n");
		PrintLog(" crd.size() = %d  neq  3*natom = %d \n",crd.size(),n3 );
		return;
	}
	if( frc.size() != n3 ) frc.resize(n3);
	if( si_vec.size() != SI_CNT) si_vec.resize(SI_CNT);
	si_vec = 0.0;
	gb_force_proxy_(&atm_cnt,crd.begin(),frc.begin(), p_amber_model->atm_mass.v(), si_vec.v(),&ncalls,
		            p_amber_model->atm_jrc.v(), p_amber_model->atm_xc.v(), p_amber_model->atm_weight.v(), p_amber_model->atm_igroup.v(), &p_amber_model->belly_atm_cnt, &p_amber_model->natc,
                    gbl_atm_owner_map.v(), gbl_my_atm_lst.v(), &my_atm_cnt );

	if (p_amber_model->ibelly > 0) ZeroFrcFrozenAtoms();
#	endif
}

void MMDriverAmber::PMEForce(int atm_cnt, HaVec_double& crd, HaVec_double& vel, HaVec_double& frc,
					  		 int new_list, HaVec_double& sys_info )
{
	int n3 = atm_cnt*3;
	if(crd.size() != n3) 
	{
		ostrstream sstr; 
		sstr << "Error in MMDriverAmber::PMEForce() \n";
		sstr << " crd.size() = " << crd.size() << " neq  3*natom = \n" << n3;
		throw sstr.str();
	}
	if(frc.size() != n3) 
	{
		ostrstream sstr; 
		sstr << "Error in MMDriverAmber::PMEForce() \n";
		sstr << " frc.size() = " << frc.size() << " neq  3*natom = \n" << n3;
		throw sstr.str();
	}
	if(vel.size() != n3) 
	{
		ostrstream sstr; 
		sstr << "Error in MMDriverAmber::PMEForce() \n";
		sstr << " vel.size() = " << vel.size() << " neq  3*natom = \n" << n3;
		throw sstr.str();
	}
	if(sys_info.size() != SI_CNT) 
	{
		ostrstream sstr; 
		sstr << "Error in MMDriverAmber::PMEForce() \n";
		sstr << " sys_info.size() = " << sys_info.size() << " neq  SI_CNT = \n" << SI_CNT;
		throw sstr.str();
	}

	HaVec_double virial(3);
	HaVec_double ekcmt(3);
	double pme_err_est = 0.0;

	virial = 0.0;

	double kin_ene = sys_info[SI_KIN_ENE];

	ekcmt  = 0.0;
	sys_info = 0.0;

//	PrintLog(" PMEForce pt 1 \n");

	pme_force_proxy_(&atm_cnt,crd.v(),gbl_atm_saved_crd.v(), pbc_box.v(), gbl_saved_box.v(),
		             vel.v(),frc.v(),p_amber_model->atm_mass.v(),&new_list,
		             p_amber_model->atm_jrc.v(), p_amber_model->atm_xc.v(), p_amber_model->atm_weight.v(), 
					 p_amber_model->atm_igroup.v(), &p_amber_model->natc, sys_info.v(), 
					 virial.v(), ekcmt.v(), &pme_err_est, gbl_atm_owner_map.v(), gbl_my_atm_lst.v(), &my_atm_cnt,
					 tranvec.v(), &p_mm_mod->run_type.value() );

	if (p_amber_model->ibelly > 0) ZeroFrcFrozenAtoms();
	
	if( p_mm_model->IsAmoebaFF() )
	{
		ekcmt[0] = kin_ene*2.0/3.0;
		ekcmt[1] = kin_ene*2.0/3.0;
		ekcmt[2] = kin_ene*2.0/3.0;
	}

	sys_info[SI_PME_ERR_EST] = pme_err_est;
	sys_info[SI_VIR_0] = virial[0];
	sys_info[SI_VIR_1] = virial[1];
	sys_info[SI_VIR_2] = virial[2];
	sys_info[SI_EKCMT_0] = 0.5*ekcmt[0];
	sys_info[SI_EKCMT_1] = 0.5*ekcmt[1];
	sys_info[SI_EKCMT_2] = 0.5*ekcmt[2];
}


int MMDriverAmber::CheckAllAtomMovement(int atm_cnt, HaVec_double& crd)
{
	int new_list = 0;
#	if defined(WITH_LIB_PMEMD)
	if(crd.size() != 3*p_amber_model->natom) 
	{
		std::strstream os; 
		os << "Error in MMDriverAmber::CheckAllAtomMovement() \n";
		os << "Size of Coordinate array =" << crd.size() << "does not match 3*natom = %d \n" << (3*p_amber_model->natom);
		throw std::runtime_error(os.str());
	}
	FC_FUNC_MODULE(nb_pairlist_mod,check_all_atom_movement)(&atm_cnt, crd.v(), gbl_atm_saved_crd.v(), pbc_box.v(), gbl_saved_box.v(), 
		                                        &p_mm_model->skin_nb, &ntp, &new_list);
#	endif
	return new_list;
}

void MMDriverAmber::ModifyFormatVal(double& val, const std::string& format)
{
	char buf[256];
	sprintf(buf,format.c_str(),val);
	
	try
	{
		//format_double_fort_(&val,format.c_str(),format.size());

		//int found_digit = FALSE;
		//int i;
		//for(i=0; i < 80; i++)
		//{
		//	buf[i] = FC_FUNC_MODULE(file_io_dat_mod,result_str)[i];
		//	if( found_digit && isspace(buf[i]) ) break;
		//	if( !found_digit && isdigit(buf[i]) ) found_digit = TRUE;
		//}
		//if( i > 75 ) throw std::runtime_error(" can't find formatted number in result_str ");
		//buf[i] = 0;
		std::istringstream is(buf);
		is >> val;
		if( is.fail() ) throw std::runtime_error((std::string)" can't convert to double " + buf );  
	}
	catch( const std::exception& ex )
	{
		PrintLog(" Error in MMDriverAmber::ModifyFormatVal() converting %16.9e \n", val);
		PrintLog(" %s\n",ex.what());
	}
}

void MMDriverAmber::PrintMDEneMDOUT(HaVec_double& si_vec,int nstep, double t)
{
	MMSysInfo info(p_mm_mod);
	SetMMInfo(info,si_vec);
	p_mm_model->UpdateDataFromFort();
	info.nstep = nstep;
	info.time = t;
	std::string str;
	//p_mm_mod->PrintEneStr(info,str);
	p_mm_mod->PrintEneStrAccurate(info,str);
	PrintLogMDOUT("\n%s\n",str.c_str());
}

void MMDriverAmber::PrintMinEneMDOUT(HaVec_double& si_vec,int nstep,double rms, double fdmax, int iatmax, std::string labmax)
{
	MMSysInfo info(p_mm_mod);
	SetMMInfo(info,si_vec);
	p_mm_model->UpdateDataFromFort();

	PrintLogMDOUT(" ");
	PrintLogMDOUT("   NSTEP       ENERGY          RMS              GMAX          NAME     NUMBER\n");
	PrintLogMDOUT(" %6d    %13.4e  %13.4e  %13.4e     %s  %7d\n",
		          nstep, si_vec[SI_POT_ENE], rms, fdmax,labmax.c_str(), iatmax);
	PrintLogMDOUT(" ");
	PrintLogMDOUT(" BOND    = %13.4f  ANGLE   = %13.4f  DIHED      = %13.4f\n",
		          si_vec[SI_BOND_ENE], si_vec[SI_ANGLE_ENE], si_vec[SI_DIHEDRAL_ENE]);
	if (p_amber_model->using_pme_potential) 
	{
		PrintLogMDOUT(" VDWAALS = %13.4f  EEL     = %13.4f  HBOND      = %13.4f\n",
		          si_vec[SI_VDW_ENE], si_vec[SI_ELECT_ENE], si_vec[SI_HBOND_ENE]);
	}
	else
	{
		PrintLogMDOUT(" VDWAALS = %13.4f  EEL     = %13.4f  EGB        = %13.4f\n",
		          si_vec[SI_VDW_ENE], si_vec[SI_ELECT_ENE], si_vec[SI_HBOND_ENE]);
	}
	PrintLogMDOUT(" 1-4 VDW = %13.4f  1-4 EEL = %13.4f  RESTRAINT  = %13.4f\n",
		          si_vec[SI_VDW_14_ENE], si_vec[SI_ELECT_14_ENE], si_vec[SI_RESTRAINT_ENE]);
	
	if(si_vec[SI_RESTRAINT_ENE] != 0.0)
	{
		PrintLogMDOUT(" E_WITHOUT_RESTR  = %13.4f \n",si_vec[SI_POT_ENE] - si_vec[SI_RESTRAINT_ENE]);
	}
	if(info.polar_ene != 0.0)
	{
		PrintLogMDOUT(" E_POLAR = %14.4f \n ", info.polar_ene );
	}
	PrintLogMDOUT(" ");
}

int MMDriverAmber::LoadAmberRestartFile(std::string rst_file_name)
{	
	PrintLog("MMDriverAmber::LoadAmberRestartFile() \n " );
	char buf[256];

	if (rst_file_name.empty()) rst_file_name = this->amber_rst_file;

	if( !boost::filesystem::exists( (std::string) rst_file_name ) )
	{	
		PrintLog(" Fail to update coordinates from AMBER Restart File \n");
		PrintLog(" File %s does not exist \n", rst_file_name );
		return FALSE;
	}

//	PrintLog(" Update coordinates from AMBER Restart File %s \n",rst_file_name);

	ifstream is(rst_file_name);

	is.getline(buf,255);
	if( is.fail() )
	{
		ErrorInMod("HaMolMechMod::LoadAmberRestartFile()", 
			       "Error Reading Title Line of Amber Restart File");
		return FALSE;
	}
	is.getline(buf,255);
	if( is.fail() )
	{
		ErrorInMod("HaMolMechMod::LoadAmberRestartFile", 
			       "Error Reading Second Line of Restart File");
		return FALSE;
	}
       
	int npt, ires;
	double time_pt;
	ires = sscanf(buf,"%d %lf",&npt, &time_pt);
	if(ires == EOF)
	{
		PrintLog(" Error in MMDriverAmber::LoadAMberRestartFile() \n");
		PrintLog(" Can not read the number of atoms and time  from the Second Line of the Restart File");
		return FALSE;
	}
	// PrintLog(" npt = %d  time = %lf \n", npt, time_pt);
	
	if( npt != p_mm_model->Atoms.size() )
	{
		PrintLog(" Error in MMDriverAmber::LoadAMberRestartFile() \n");
		PrintLog(" The number of Atoms in the Molecular Mechanics module %d \n is not equal to number of Atoms in Restart file %d \n",
			       p_mm_model->Atoms.size(),npt);
		return FALSE;
	}
		
	HaVec_double darr(3*npt);

// Reading atom coordinate
	int i;
	for(i = 0 ; i < 3*npt; i++)
	{
		is >> darr[i];
		if(is.fail())
		{
			ErrorInMod("HaMolMechMod::LoadAMberRestartFile", 
				       "Error reading Amber Restart file");
			return FALSE;
		}
	}



	for(i = 0; i < npt; i++ )
	{
		p_mm_model->Atoms[i]->SetX( darr[3*i]   );
		p_mm_model->Atoms[i]->SetY( darr[3*i+1] );
		p_mm_model->Atoms[i]->SetZ( darr[3*i+2] );
	}

// Reading Atom Velocities:
	int vel_found = TRUE;
	for(i = 0 ; i < 3*npt; i++)
	{
		is >> darr[i];
		if(is.fail())
		{
			vel_found = FALSE;
			break;
		}
	}

	if( vel_found == FALSE  )
	{
		if( i == 3)
		{
			pmset->per_bc->SetBox(darr[0],darr[1],darr[2] );
			return TRUE;
		}
		else
		{
			return FALSE;
		}
	}
    
    darr.newsize(3);
	for( i = 0 ; i < 3; i++)
	{
		is >> darr[i];
		if(is.fail())
		{
			return FALSE;
		}
	}
	pmset->per_bc->SetBox(darr[0], darr[1], darr[2] );

	return TRUE;
}

int MMDriverAmber::LoadAmberMDInfoFile()
{
	char buf[256];
	std::string fname = "mdinfo";
	FILE* finfo = fopen(fname.c_str(),"r");

	if(finfo == NULL)
	{
		ErrorInMod("HaMolMechMod::LoadAmberMDInfoFile()",
			       "Unable to open mdinfo file");
		return FALSE;
	}
	
	char* cres = NULL;

	int irec = 0;
	for(;;)
	{
		cres = fgets(buf,255,finfo);
		irec++;
		if(cres == NULL) 
			break;
		
		int len = strlen(buf);
		int i;
		for( i = 0; i < len; i++)
		{
			if(buf[i] == '=') buf[i] = ' ';
			if(buf[i] == 10) buf[i] = ' ';
		}
		
		std::string str(buf);
		// PrintLog(" %s \n", str.c_str());

		char* spos;
		
		if( p_mm_mod->run_type == p_mm_mod->run_type.MD_RUN)
		{
			spos = strstr(buf,"NSTEP");
			if (spos != NULL)
			{
				sscanf(spos + 6, "%d", &p_mm_mod->p_mm_info->nstep);
				// PrintLog("nstep = %d \n", p_mm_mod->p_mm_info->nstep);
			}
			
			spos = strstr(buf,"TIME(PS)");
			if(spos != NULL) sscanf(spos+9,"%lf",&p_mm_mod->p_mm_info->time);
			
			spos = strstr(buf,"TEMP(K)");
			if(spos != NULL) sscanf(spos+8,"%lf",&p_mm_mod->p_mm_info->temp);
			
			spos = strstr(buf,"PRESS");
			if(spos != NULL) sscanf(spos+6,"%lf",&p_mm_mod->p_mm_info->press);

			spos = strstr(buf,"Etot");
			if(spos != NULL) sscanf(spos+5,"%lf",&p_mm_mod->p_mm_info->tot_energy);

			spos = strstr(buf,"EKtot");
			if(spos != NULL) sscanf(spos+6,"%lf",&p_mm_mod->p_mm_info->kin_ene);

			spos = strstr(buf,"EPtot");
			if(spos != NULL) sscanf(spos+6,"%lf",&p_mm_mod->p_mm_info->pot_ene);

			spos = strstr(buf,"EELEC");
			if(spos != NULL) sscanf(spos+6,"%lf",&p_mm_mod->p_mm_info->electr_ene);

			spos = strstr(buf,"1-4 NB");
			if(spos != NULL) sscanf(spos+7,"%lf",&p_mm_mod->p_mm_info->vdw_ene_14);

			spos = strstr(buf,"EHBOND");
			if(spos != NULL) sscanf(spos+7,"%lf",&p_mm_mod->p_mm_info->hbond_ene);
		}
		else if( p_mm_mod->run_type == p_mm_mod->run_type.MIN_RUN)
		{
			if( irec == 4)
			{
				sscanf(buf," %d %lf %lf %lf ", &p_mm_mod->p_mm_info->nstep, &p_mm_mod->p_mm_info->tot_energy,&p_mm_mod->p_mm_info->rms_ene,&p_mm_mod->p_mm_info->grad_ene_max);
			}

			spos = strstr(buf,"EEL");
			if(spos != NULL) sscanf(spos+4,"%lf",&p_mm_mod->p_mm_info->electr_ene);

			spos = strstr(buf,"HBOND");
			if(spos != NULL) sscanf(spos+6,"%lf",&p_mm_mod->p_mm_info->hbond_ene);

			spos = strstr(buf,"1-4 VDW");
			if(spos != NULL) sscanf(spos+8,"%lf",&p_mm_mod->p_mm_info->vdw_ene_14);

			p_mm_mod->p_mm_info->vdw_ene_nb = p_mm_mod->p_mm_info->vdw_ene - p_mm_mod->p_mm_info->vdw_ene_14;

		}

		spos = strstr(buf,"VDWAALS");
		if(spos != NULL) sscanf(spos+8,"%lf",&p_mm_mod->p_mm_info->vdw_ene);

		spos = strstr(buf," BOND ");
		if(spos != NULL) sscanf(spos+7,"%lf",&p_mm_mod->p_mm_info->bond_ene);
		
		spos = strstr(buf,"DIHED");
		if(spos != NULL) sscanf(spos+6,"%lf",&p_mm_mod->p_mm_info->dihed_ene);
		
		spos = strstr(buf,"ANGLE");
		if(spos != NULL) sscanf(spos+6,"%lf",&p_mm_mod->p_mm_info->vang_ene);

		spos = strstr(buf,"1-4 EEL");
		if(spos != NULL) sscanf(spos+8,"%lf",&p_mm_mod->p_mm_info->electr_ene_14);

		p_mm_mod->p_mm_info->electr_ene_nb = p_mm_mod->p_mm_info->electr_ene - p_mm_mod->p_mm_info->electr_ene_14;

		spos = strstr(buf,"CONSTRAINTS");
		if(spos != NULL) sscanf(spos+8,"%lf",&p_mm_mod->p_mm_info->constraints_ene);
	}

	fclose(finfo);
	return TRUE;
}

int MMDriverAmber::SaveAllInpFiles()
{
	int ires = p_mm_model->UpdateModel();
    if( p_mm_mod->to_init_simulations ) ires = p_mm_mod->InitMMSimulations();
	if(!ires) return FALSE;

	SaveAmberRunFile();
	SaveAmberInpFile();
	SaveAmberCrdFile();
    if( p_mm_model->GetRestrAtoms() )
	{
		SaveAmberRstFile(amber_constr_crd_file.c_str());
	}
	SaveAmberTopFile();
	to_save_input_files = FALSE;
	return TRUE;
}

int MMDriverAmber::DeleteOutputFiles()
{
	ha_delete_file( amber_out_file.c_str());
	ha_delete_file( amber_rst_file.c_str() );
	ha_delete_file( amber_trj_coord_file.c_str() );
	ha_delete_file( amber_trj_vel_file.c_str() );
	ha_delete_file( amber_trj_ene_file.c_str() );	
	return TRUE;
}

AmberMMModel::AmberMMModel(MolMechModel* p_mm_model_new)
{
	p_mm_model = p_mm_model_new;
	p_mm_mod   = p_mm_model->p_mm_mod;
	p_amber_driver = p_mm_mod->p_amber_driver;
	pmset = p_mm_mod->GetMolSet();
	
	SetAtomNum( pmset->GetNAtoms() );
	
	Clear();
}

AmberMMModel::~AmberMMModel()
{

}

int AmberMMModel::SetAtomNum(int natom_new)
{
	if( natom_new < 0 ) 
	{
		PrintLog("Error in AmberMMModel::SetAtomNum() \n");
		PrintLog("natom_new = %d < 0 \n",natom_new );
		return FALSE;
	}
	natom = natom_new;
	if( atm_xc.size() != natom ) atm_xc.resize( 3*natom );
	return TRUE;

}

int AmberMMModel::UpdateAmberData()
{
	int ires = TRUE;

	MolSet* pmset = p_mm_model->GetMolSet();

	p_amber_driver->title = " Amber topology and parameter file for the molecule set ";
	p_amber_driver->title += pmset->GetName();

	int extra_spaces = 80 - p_amber_driver->title.size();
	for(; extra_spaces > 0; extra_spaces-- )
	{
		p_amber_driver->title += " ";
	}
	 
	ntf = p_mm_model->omit_interactions;
	if( pmset->per_bc->IsSet() )
	{
		if( p_mm_model->GetNBCutDist() > 0.35*pmset->per_bc->GetA() ) p_mm_model->SetNBCutDist( 0.35*pmset->per_bc->GetA() );
		if( p_mm_model->GetNBCutDist() > 0.35*pmset->per_bc->GetB() ) p_mm_model->SetNBCutDist( 0.35*pmset->per_bc->GetB() );
		if( p_mm_model->GetNBCutDist() > 0.35*pmset->per_bc->GetC() ) p_mm_model->SetNBCutDist( 0.35*pmset->per_bc->GetC() );
	}

	alpb = 0;       // used for Analytical Linearized Poisson-Boltzmann
	arad = 15.0;    // used for Analytical Linearized Poisson-Boltzmann
	intdiel = 1.0;  // used for generalized Born
	extdiel = 78.5; // used for generalized Born
	saltcon = 0.0;  // used for generalized Born
	rgbmax = 25.0;  // used for generalized Born
	rbornstat = 0;  // used for generalized Born
	offset = 0.09;  // used for generalized Born
	gbsa = 0;        // used for generalized Born
	surften = 0.005; // used for generalized Born
	p_amber_driver->nrespai = 1;     // used for generalized Born
	cut_inner = 8.0; // used for generalized Born

	// Generalized Born variables, perhaps better factored elsewhere, but...

	gb_alpha = 0.0;
	gb_beta = 0.0;
	gb_gamma = 0.0;
	gb_fs_max = 0.0;
	gb_kappa = 0.0;
	gb_neckscale = 1.0/3.0;

	// Boundary box dimensions for non periodic caclulations restrained in a box:

	bbox_xmin = 0.0;
	bbox_ymin = 0.0;
	bbox_zmin = 0.0;
	bbox_xmax = 0.0;
	bbox_ymax = 0.0;
	bbox_zmax = 0.0;

	igb = 0;

	if(p_mm_model->electr_method == p_mm_model->electr_method.PME_METHOD)
	{
		if(p_mm_mod->period_bcond == p_mm_mod->period_bcond.NO_PERIODICITY ) 
			igb = 2;
		else
			igb = 0;
	}
	else if( p_mm_model->electr_method == p_mm_model->electr_method.DIST_DEP_DIEL)
	{
		igb = 3;
	}
	else if( p_mm_model->electr_method == p_mm_model->electr_method.GEN_BORN)
	{
		igb = 1;
	}
	else if( p_mm_model->electr_method == p_mm_model->electr_method.SCREENED_COULOMB)
	{
		igb = 20;
	}

	if(igb != 0)
	{
		if (igb == 1) 
		{
			gb_alpha = 1.0;
			gb_beta = 0.0;
			gb_gamma = 0.0;
		}
		else if (igb == 2)
		{
			// Use our best guesses for Onufriev/Case GB  (GB^OBC I):
			gb_alpha = 0.8;
			gb_beta = 0.0;
			gb_gamma = 2.9091250;
		}
		else if (igb == 5) 
		{
			// Use our second best guesses for Onufriev/Case GB (GB^OBC II):
			gb_alpha = 1.0;
			gb_beta = 0.8;
			gb_gamma = 4.851;
		}
		else if (igb == 7)
		{
			// Use parameters for Mongan et al. CFA GBNECK:
			gb_alpha = 1.09511284;
			gb_beta = 1.90792938;
			gb_gamma = 2.50798245;
			gb_neckscale = 0.361825;
		}
		else if (igb == 20)
		{
			// Screened coulomb interactions - do not use gb_alpha,gb_beta,gb_gamma
			gb_alpha = 0.0;
			gb_beta = 0.0;
			gb_gamma = 0.0;
		}

		// Set gb_kappa as long as saltcon is .ge. 0.d0 (and catch bug below if not).

		if (saltcon >= 0.0)
		{
			// Get Debye-Huckel kappa (A**-1) from salt concentration (M), assuming:
			//  T = 298.15, epsext=78.5,
			gb_kappa = sqrt(0.10806 * saltcon);

			// Scale kappa by 0.73 to account(?) for lack of ion exclusions:

			if( igb != 20) // do not scale kappa for screen electrostatics model
			{       
				gb_kappa = 0.73 * gb_kappa;
			}
		}
	}
	// Here we reset some possibly errant variables without killing the run.

	if(dielc <= 0.0) dielc = 1.0;
	if(scnb == 0.0)  scnb = 2.0;

	if (igb == 0) 
	{
		using_pme_potential = 1;
		using_gb_potential  = 0;
	}
	else if (igb == 1 || igb == 2 || igb == 5 || igb == 7 || igb == 20) 
	{
		using_pme_potential = 0;
		using_gb_potential  = 1;
	}
	else                                
	{
		using_pme_potential  = 0;
		using_gb_potential   = 0;
	}

	if( p_mm_model->ff_type != ForceFieldType::AMOEBA )
	{
		natom = p_mm_model->Atoms.size();
	}

	AtomIntMap& at_idx_map = p_mm_model->GetAtIdxMap(TRUE);

	if( p_mm_model->ff_type == ForceFieldType::AMOEBA )
	{
		SetDistConstrData();
		return ires;
	}

	int i;

	int i1, i2, i3, i4;
	HaAtom *pt1, *pt2, *pt3, *pt4;
	int contain_h;
	HaVec_double bpar(2);
	HaVec_double vang_par(2);
	HaVec_double dang_par;

// Set Valence Bonds arrays
	std::vector<MMBond*>   bonds_h;
	std::vector<MMBond*>   bonds_a;
	loc_bond_params.clear();
	bnd_par_idx_map.clear();

	map<HaVec_double, int, less<HaVec_double> >::iterator mbitr;
	numbnd = 0;

	set<MMBond, less<MMBond> >::iterator bndset_itr = p_mm_model->MBonds.begin();

	for( ; bndset_itr != p_mm_model->MBonds.end(); bndset_itr++)
	{
		MMBond& mb = (MMBond&) (*bndset_itr);
        
		bpar[0] = mb.r0;
		bpar[1] = mb.fc;

		mbitr = bnd_par_idx_map.find(bpar);
		if(mbitr == bnd_par_idx_map.end() )
		{
			loc_bond_params.push_back(bpar);
			numbnd++;
			bnd_par_idx_map[bpar] = numbnd;
		}

		contain_h = FALSE;
		pt1 = mb.pt1;
		pt2 = mb.pt2;

		if( pt1 == NULL || pt2 == NULL)
		{
			ErrorInMod(" AmberMMModel::UpdateAmberData() ",
				       " One of the HaAtom of the Bond is NULL ");
			continue;
		}

		if( pt1->IsHydrogen() || pt2->IsHydrogen() )
		{
			contain_h = TRUE;
		}

		if( contain_h)
		{
			bonds_h.push_back(&mb);
		}
		else
		{
			bonds_a.push_back(&mb);
		}
	}

// Set Valence Angles Arrays

	std::vector<MMValAngle*> val_h;
	std::vector<MMValAngle*> val_a;
	loc_val_angle_params.clear();
	vang_par_idx_map.clear();
	numang = 0;

	map<HaVec_double,int, less<HaVec_double> >::iterator vaitr;
	set<MMValAngle, less<MMValAngle> >::iterator vaset_itr = p_mm_model->ValAngles.begin();

	for( ; vaset_itr != p_mm_model->ValAngles.end(); vaset_itr++)
	{
		MMValAngle& va =  (MMValAngle&) (*vaset_itr);

		vang_par[0] = va.a0;
		vang_par[1] = va.fc;
		
		vaitr = vang_par_idx_map.find(vang_par);
		if( vaitr == vang_par_idx_map.end() )
		{
			loc_val_angle_params.push_back(vang_par);
			numang++;
			vang_par_idx_map[vang_par]= numang;
		}

		contain_h = FALSE;
		pt1 = va.pt1;
		pt2 = va.pt2;
		pt3 = va.pt3;

		if( pt1 == NULL || pt2 == NULL || pt3 == NULL )
		{
			ErrorInMod(" AmberMMModel::UpdateAmberData() ",
				" One of the HaAtom of the Valence Angle is NULL ");
			continue;
		}

		if( pt1->IsHydrogen() || pt2->IsHydrogen() || pt3->IsHydrogen() )
		{
			contain_h = TRUE;
		}

		if( contain_h)
		{
			val_h.push_back(&va);
		}
		else
		{
			val_a.push_back(&va);
		}
	}
// Set Dihedral Angles Arrays

	std::vector<MMDihedral*> dih_h;  
	std::vector<MMDihedral*> dih_a; 
	loc_dih_ang_par.clear(); 
	dang_par_idx_map.clear(); 

	nptra = 0; 
	nphih = 0;
	nphia = 0; 

	map<HaVec_double, int, less<HaVec_double> >::iterator ditr;

	vector<MMDihedral>::iterator lditr;
	vector<MMDihedral> *dih_list;

	for(int im= 0 ; im < 2; im++)
	{
		dih_list = &p_mm_model->Dihedrals;
		if( im == 1)
			dih_list = &p_mm_model->ImprDihedrals;

		for( lditr = dih_list->begin(); lditr != dih_list->end(); lditr++)
		{
			MMDihedral& dang = *lditr;
			int nt = dang.GetNTerms();
			dang_par.newsize(4*nt);
			int it;
			for(it = 0; it < nt; it++)
			{
				dang_par[4*it]   = dang.pn[it];
				dang_par[4*it+1] = dang.phase[it];
				dang_par[4*it+2] = dang.pk[it];
				dang_par[4*it+3] = dang.idivf[it];
			}

			ditr = dang_par_idx_map.find( dang_par );
			if( ditr == dang_par_idx_map.end() )
			{
				loc_dih_ang_par.push_back(dang_par);
				dang_par_idx_map[dang_par] = nptra + 1;
				nptra += nt;
			}
			
			contain_h = FALSE;
			pt1 = dang.pt1;
			pt2 = dang.pt2;
			pt3 = dang.pt3;
			pt4 = dang.pt4;
			
			if( pt1 == NULL || pt2 == NULL || pt3 == NULL || pt4 == NULL)
			{
				ErrorInMod(" AmberMMModel::UpdateAmberData() ",
					" One of the HaAtom of the Dihedral Angle is NULL ");
				continue;
			}
			
			if( pt1->IsHydrogen() || pt2->IsHydrogen() || pt3->IsHydrogen() || pt4->IsHydrogen())
			{
				contain_h = TRUE;
			}
			
			MMDihedral* ptr_dih = &dang;
			if( contain_h)
			{
				dih_h.push_back(ptr_dih);
				nphih += nt;
			}
			else
			{
				dih_a.push_back(ptr_dih);
				nphia += nt;
			}
		}
	}

// Distribute atoms over residues and molecules 
// ( Do not quite correspond to HARLEM classification)

	PrintLog("Before FindResMolPartition() \n");
	FindResMolPartition();
	PrintLog("After FindResMolPartition() \n");

// Set Van-Der-Waals parameters of the Atoms 

	loc_point_params.clear();
	ppar_idx_map.clear();

	HaVec_double ppar(2);
	
	atm_iac.clear();
	atm_iac.newsize(natom);
	
	map<HaVec_double,int, less<HaVec_double> >::iterator pitr;
	ntypes = 0;
	for( i =0; i < natom ; i++) 
	{
		ppar[0] = p_mm_model->Atoms[i]->vdw_rad;
		ppar[1] = p_mm_model->Atoms[i]->ew;

		pitr = ppar_idx_map.find(ppar);
		int idx;
		if( pitr == ppar_idx_map.end())
		{
			loc_point_params.push_back(ppar);
			ntypes++;
			ppar_idx_map[ppar]= ntypes;
			idx = ntypes;
		}
		else
		{
			idx= (*pitr).second;
		}
		atm_iac[i] = idx;
	}

// Set H-Bonds arrays

//	int idx_fst_hb_atom_type = ntypes + 1; // index of the first atom type involved in H-Bond   

//	map<HaVec_double, int, less<HaVec_double> > hbnd_par_idx_map; // map H-Bond Parameters to indexes in loc_hbond_params
//	HaVec_double hbpar(2);

//	map<HaVec_double, int, less<HaVec_double> >::iterator hbitr;

//	vector<AtomContact>::iterator cnt_itr = DistConstraints.begin();

	nbonh  = bonds_h.size();      
	nbona  = bonds_a.size();        
	ntheth = val_h.size();          
	ntheta = val_a.size(); 
	nphb = 0;
	gbl_asol.clear();
	gbl_bsol.clear();             
	nttyp    = ntypes * (ntypes + 1) / 2;
	bonda_idx  = nbonh + 1;
	anglea_idx = ntheth + 1;
	diheda_idx = nphih + 1;  
	gbl_bond_allocsize  = nbonh + nbona;
	if( nbona == 0 ) gbl_bond_allocsize++; 
	gbl_angle_allocsize = ntheth + ntheta;
	if( ntheta == 0 ) gbl_angle_allocsize++; 
	gbl_dihed_allocsize = nphih + nphia;
	if( nphia == 0 ) gbl_dihed_allocsize++;
	
//	for( ; cnt_itr != DistConstraints.end(); cnt_itr++)
//	{
//		AtomContact& hb = (AtomContact&) (*cnt_itr);
//      if( hb.cnt_type != 1) continue;
//      
//		hbpar[0] = hb.cf[0];
//		hbpar[1] = hb.cf[1];

//		loc_hbond_params.push_back(hbpar);
//		iprm.nphb++;
			
//		ppar[0] = hb.pt1->vdw_rad;
//		ppar[1] = hb.pt1->ew;

//		loc_point_params.push_back(ppar);
//		ntypes++;

//		int iat = pt_indx_map[hb.pt1];
//		atm_iac[iat] = ntypes;


//		ppar[0] = hb.pt2->vdw_rad;
//		ppar[1] = hb.pt2->ew;

//		loc_point_params.push_back(ppar);
//		ntypes++;

//      iat = pt_indx_map[hb.pt2];
//		atm_iac[iat] = ntypes;
//	}

// Non-bonded iteractions parameters index 

// Atom Charges

	atm_charge.newsize(natom);

	for(i =0; i < natom ; i++)
	{
		atm_charge[i] = p_mm_model->Atoms[i]->GetCharge() * 18.2223; // to transform atomic charges so for distances 
	                                               // expressed in Ang electrostatic energy 
	                                               // will be in kcal/mol
		if( p_amber_driver->emulate_ext_amber )
		{
			MMDriverAmber::ModifyFormatVal(atm_charge[i],FLOAT_E16_8);
		}
	}
	double dlc = 1.0; 
	if( using_pme_potential && p_mm_model->diel_const > 1.0 )  dlc = p_mm_model->diel_const;
	for( i = 0; i < natom; i++)
	{
		atm_charge[i] = atm_charge[i]/dlc;
	}

// Atom Masses:

	atm_mass.resize(natom);
	atm_mass_inv.resize(natom);
	tmass = 0.0;
	for(i =0; i < natom ; i++)
	{
		 atm_mass[i] = p_mm_model->Atoms[i]->GetMass();
		if( atm_mass[i] < 1e-6)
		{
			PrintLog("Warning in AmberMMModel::UpdateAmberData() \n");
			PrintLog(" Atom Mass of atom %s is invalid = %16.9f   it is set to 1.0 \n",
				       (p_mm_model->Atoms[i]->GetRef()).c_str(),atm_mass[i]);
			atm_mass[i] = 1.0;
		}
		tmass += atm_mass[i];
		atm_mass_inv[i] = 1.0/atm_mass[i];
	}

// Atom Names:

	int len,j,idiff;
	atm_igraph.resize(natom);

	for(i = 0; i < natom; i++)
	{
		HaAtom* aptr = (HaAtom*) p_mm_model->Atoms[i];
		atm_igraph[i] = aptr->GetName();
		len = atm_igraph[i].size();
		idiff = 4 - len;
		if( idiff > 0) 
		{
			for(j = 0; j < idiff; j++)
				atm_igraph[i] += ' '; 
		}
	}

// Atom Force Field symbols

	atm_isymbl.resize(natom);
	for( i=0; i < natom; i++)
	{
		atm_isymbl[i] = p_mm_model->Atoms[i]->GetFFSymbol();
		len = atm_isymbl[i].size();
		idiff = 4 - len;
		if( idiff > 0) 
		{
			for(j = 0; j < idiff; j++)
				atm_isymbl[i] += ' '; 
		}
	}

// ATOM TREE SYMBOLS
	
	atm_itree.resize(natom);
	for( i=0; i < natom; i++)
	{
		atm_itree[i] = " M  ";
	}

// Atom Generalized Born radii 

	atm_gb_radii.newsize(natom);

	for(i =0; i < natom ; i++)
	{
		double rad = MolMechModel::GetGBAtomRad(p_mm_model->Atoms[i],p_mm_model->gb_param_type);
		atm_gb_radii[i] = rad;
	}

// Atom Generalized Born screening Factors 

	atm_gb_fs.newsize(natom);

	for(i =0; i < natom ; i++)
	{
		double screen = MolMechModel::GetGBAtomScreening(p_mm_model->Atoms[i],p_mm_model->gb_param_type);
		atm_gb_fs[i] = screen;
	}

// VdW parameters

	typ_ico.newsize(ntypes*ntypes);
	for( i = 0; i < ntypes; i++)
	{
		for( j = 0; j < ntypes; j++)
		{
			int i_fort = i+1;
			int j_fort = j+1;
			
			if( i < j) 
			{
				i_fort = j+1;
				j_fort = i+1;
			}

			int idx = ntypes* i + j ;
			typ_ico[idx] = (i_fort*(i_fort - 1))/2 + j_fort;
		}
	}

	gbl_cn1.newsize(nttyp);
	gbl_cn2.newsize(nttyp);

	for( i = 0; i < ntypes ; i++)
	{
		for( j =0; j <= i; j++)
		{
			int idx_lin = i*(i+1)/2 + j;
			
			HaVec_double ppar1 = loc_point_params[i];
			HaVec_double ppar2 = loc_point_params[j];

			double e0 = sqrt( fabs( ppar1[1] * ppar2[1]) ); // Minimum VdW energy of 2 atoms in kcal/mol
			double r0 = (ppar1[0] + ppar2[0]); // Minimum energy distance between  
                                                                        // two atoms in Angstroms
			double r6= (r0 * r0);
			r6 = r6 * r6 * r6;   // r0^6 

			gbl_cn1[idx_lin] = e0 * r6 * r6;
			gbl_cn2[idx_lin] = 2.0 * e0 * r6;

			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(gbl_cn1[idx_lin],FLOAT_E16_8);
				MMDriverAmber::ModifyFormatVal(gbl_cn2[idx_lin],FLOAT_E16_8);
			}
		}
	}

// Excluded Atoms Arrays

	atm_numex.newsize(natom);
	for( i =0; i < natom ; i++) 
	{
		atm_numex[i] = p_mm_model->excluded_atom_list[i].size();
		if(atm_numex[i] == 0) atm_numex[i] = 1;
	}

// Bond parameters

	gbl_rk.resize(numbnd);
	gbl_req.resize(numbnd);

	for( i = 0; i < numbnd; i++)
	{
		gbl_rk[i] = loc_bond_params[i][1];
		gbl_req[i] = loc_bond_params[i][0];

		if( p_amber_driver->emulate_ext_amber )
		{
			MMDriverAmber::ModifyFormatVal(gbl_rk[i],FLOAT_E16_8);
			MMDriverAmber::ModifyFormatVal(gbl_req[i],FLOAT_E16_8);
		}
	}

// Valence Angle parameters

	gbl_tk.newsize(numang);
	gbl_teq.newsize(numang);

	for( i = 0; i < numang; i++)
	{
		gbl_tk[i] = loc_val_angle_params[i][1];
		gbl_teq[i] = DEG_TO_RAD*loc_val_angle_params[i][0];

		if( p_amber_driver->emulate_ext_amber )
		{
			MMDriverAmber::ModifyFormatVal(gbl_tk[i],FLOAT_E16_8); 
			gbl_teq[i] *= (3.141594/PI);
			MMDriverAmber::ModifyFormatVal(gbl_teq[i],FLOAT_E16_8);
		}
	}

// Dihedral Angle parameters

	gbl_pk.resize(nptra);
	gbl_phase.resize(nptra);
	gbl_pn.resize(nptra);

	vector<HaVec_double>::iterator ldpar_itr;

	i =0;
	for(ldpar_itr= loc_dih_ang_par.begin(); ldpar_itr != loc_dih_ang_par.end(); ldpar_itr++ )
	{
		HaVec_double& dpar = *ldpar_itr;
		int nt = dpar.size()/4;
		for(int it = 0; it < nt; it++)
		{
			gbl_pn[i] = dpar[4*it];
			gbl_phase[i] = DEG_TO_RAD * dpar[4*it+1];
			gbl_pk[i] = dpar[4*it + 2];
			
			if( p_amber_driver->emulate_ext_amber )
			{
				MMDriverAmber::ModifyFormatVal(gbl_pn[i],FLOAT_E16_8);
				gbl_phase[i] *= (3.141594/PI);
				MMDriverAmber::ModifyFormatVal(gbl_phase[i],FLOAT_E16_8); // tleap is using incorrect PI value
				MMDriverAmber::ModifyFormatVal(gbl_pk[i],FLOAT_E16_8);
			}
			i++;
		}
	}

// Valence Bonds

	gbl_bond.resize(3*gbl_bond_allocsize);
	for( i=0; i < nbonh; i++)
	{
		MMBond* bnd = bonds_h[i];
		i1 = at_idx_map[bnd->pt1] + 1;
		i2 = at_idx_map[bnd->pt2] + 1;
		gbl_bond[3*i]   = i1;
		gbl_bond[3*i+1] = i2;
		bpar[0] = bnd->r0;
		bpar[1] = bnd->fc;
		gbl_bond[3*i+2] = bnd_par_idx_map[bpar];
	}

	for( i=0; i < nbona; i++)
	{
		int idx = i + nbonh;
		MMBond* bnd = bonds_a[i];
		i1 = at_idx_map[bnd->pt1] + 1;
		i2 = at_idx_map[bnd->pt2] + 1;
		gbl_bond[3*idx]   = i1;
		gbl_bond[3*idx+1] = i2;
		bpar[0] = bnd->r0;
		bpar[1] = bnd->fc;
		gbl_bond[3*idx+2] = bnd_par_idx_map[bpar];
	}

	gbl_angle.resize(4*gbl_angle_allocsize);
	for( i=0; i < ntheth; i++)
	{
		MMValAngle* ang = val_h[i];
		i1 = at_idx_map[ang->pt1] + 1;
		i2 = at_idx_map[ang->pt2] + 1;
		i3 = at_idx_map[ang->pt3] + 1;

		gbl_angle[4*i]   = i1;
		gbl_angle[4*i+1] = i2;
		gbl_angle[4*i+2] = i3;
		vang_par[0] = ang->a0;
		vang_par[1] = ang->fc;
		gbl_angle[4*i+3] = vang_par_idx_map[vang_par];
	}

	for( i=0; i < ntheta; i++)
	{
		int idx = ntheth + i;
		MMValAngle* ang = val_a[i];
		i1 = at_idx_map[ang->pt1] + 1;
		i2 = at_idx_map[ang->pt2] + 1;
		i3 = at_idx_map[ang->pt3] + 1;

		gbl_angle[4*idx]   = i1;
		gbl_angle[4*idx+1] = i2;
		gbl_angle[4*idx+2] = i3;

		vang_par[0] = ang->a0;
		vang_par[1] = ang->fc;
		gbl_angle[4*idx+3] = vang_par_idx_map[vang_par];
	}

	gbl_dihed.resize(5*gbl_dihed_allocsize);
	i =0;
	int idx_d;

	int nd_tot = dih_h.size() + dih_a.size();
	int dih_h_size = dih_h.size();
	for( idx_d = 0; idx_d < nd_tot; idx_d++) 
	{
		MMDihedral* dih = NULL;

		if( idx_d < dih_h_size )
		{
			dih = dih_h[idx_d];
		}
		else
		{
			dih = dih_a[idx_d - dih_h_size];
		}		
		i1 = at_idx_map[dih->pt1] + 1;
		i2 = at_idx_map[dih->pt2] + 1;
		i3 = at_idx_map[dih->pt3] + 1;
		i4 = at_idx_map[dih->pt4] + 1;

		if( i3 < i2 )  // pmemd require i3 > i2 otherwise it will not set 1-4 interactions flag correctly
		{
			int i1_save = i1;
			int i2_save = i2;
			i1 = i4;
			i2 = i3;
			i3 = i2_save;
			i4 = i1_save;
		}

		int nt = dih->GetNTerms();

		HaVec_double dang_par(nt*4);
		int it;
		for( it = 0; it < nt; it++)
		{
			dang_par[4*it]     = dih->pn[it];  
			dang_par[4*it + 1] = dih->phase[it]; 
			dang_par[4*it + 2] = dih->pk[it]; 
			dang_par[4*it + 3] = dih->idivf[it]; 
		}

		for( it = 1; it <= nt; it++)
		{
			gbl_dihed[5*i]   = i1;         // IPH : atom involved in dihedral "i", dihedral contains hydrogen
			gbl_dihed[5*i+1] = i2;         // JPH : atom involved in dihedral "i", dihedral contains hydrogen
			gbl_dihed[5*i+2] = i3;         // KPH : atom involved in dihedral "i", dihedral contains hydrogen
			gbl_dihed[5*i+3] = i4;         // LPH : atom involved in dihedral "i", dihedral contains hydrogen
			
			gbl_dihed[5*i+4] = dang_par_idx_map[dang_par]+ it - 1;  // ICPH : index into parameter arrays PK, PN, and PHASE for
	                                                             // dihedral IPH(i)-JPH(i)-KPH(i)-LPH(i)
			if( dih->improper )                                  // take into account multiterm dihedrals
			{
				gbl_dihed[5*i+3] = -gbl_dihed[5*i+3];
			}
			if( (nt > 1 && it < nt) || dih->improper || (!dih->calc_14)) // indicate not to calculate 1-4 interactions
			{
				gbl_dihed[5*i+2] = -gbl_dihed[5*i+2];
			}
	//		PrintLog("Dihedral index %d   third index %d \n", i,iarr[5*i + 3]);
	//		PrintLog("it %d   dih->improper %d dih->calc_14 %d \n", it,dih->improper, dih->calc_14);
			i++;
		}
	}

// Nonbonded interations atom excluded list

	// Compute length of excluded non-bonded list

	next = 0; 
	for(i = 0; i < p_mm_model->excluded_atom_list.size(); i++)
	{
		int ns = p_mm_model->excluded_atom_list[i].size();
		if(ns != 0)
		{
			next +=  p_mm_model->excluded_atom_list[i].size();
		}
		else
		{
			next += 1;
		}
	}
	gbl_natex.resize(next);

	int idx = 0;
	set<HaAtom*, less<HaAtom*> >::iterator si_itr;

	for(i = 0; i < natom; i++)
	{
		if(!p_mm_model->excluded_atom_list[i].empty())
		{	
			int nex = p_mm_model->excluded_atom_list[i].size();
		    
			list<int> axx_list;
			list<int>::iterator litr;

			si_itr = p_mm_model->excluded_atom_list[i].begin();
			for(; si_itr != p_mm_model->excluded_atom_list[i].end(); si_itr++)
			{
				HaAtom* pt1 = *si_itr;
				i1 = at_idx_map[pt1] + 1;
				axx_list.push_back(i1);
			}
			axx_list.sort();

			for(litr = axx_list.begin(); litr != axx_list.end(); litr++)
			{
				gbl_natex[idx] = (*litr);
				idx++;
			}
		}
		else // no excluded atoms for the atom   ( set list = 0 )
		{
//			gbl_natex[idx] = i+1;
			gbl_natex[idx] = 0;
			idx++;
		}
	}

	CalcAddDihParams();
  
	if ( using_gb_potential )
	{
		// If igb .eq. 7 use special S_x screening params; here we overwrite
		// the tinker values read from the prmtop.

		if (igb == 7) 
		{
			PrintLog(" Replacing prmtop screening parameters with GBn (igb=7) values \n");
			for(i = 0; i < natom; i++)
			{
				int elemno = p_mm_model->Atoms[i]->GetElemNo();
				if (elemno == 6) // Carbon
	  			{
					atm_gb_fs[i] = 4.84353823306e-1;
				}
				else if (elemno == 1) // Hydrogen 
				{
					atm_gb_fs[i] = 1.09085413633e0;
				}
				else if (elemno == 7) // Nitrogen
				{
					atm_gb_fs[i] = 7.00147318409e-1;
				}
				else if (elemno == 8)  // Oxygen
				{
					atm_gb_fs[i] = 1.06557401132e0;
				}
				else if (elemno == 16)  // Sulphur
				{
					atm_gb_fs[i] = 6.02256336067e-1;
				}
				else
				{
					atm_gb_fs[i] = 0.5;
				}
	 		}
				// Put fs(i) * (rborn(i) - offset) into the "fs" array:
		}
		gb_fs_max = 0.0;

		for(i = 0; i < natom; i++)
		{
			atm_gb_fs[i] = atm_gb_fs[i] * (atm_gb_radii[i] - offset);
			gb_fs_max = MaxFun(gb_fs_max, atm_gb_fs[i]);
		}
	}
    
	
	SetAtomPosRestrData();
	SetDistConstrData();
	SetMovingAtomsData();
	  
	to_update_amber_data = FALSE;
	return TRUE;
}

int AmberMMModel::SetAtomPosRestrData()
{
	AtomGroup* p_atgrp_restr = p_mm_model->GetRestrAtoms();
	if( p_atgrp_restr == NULL)
	{
		natc = 0;
		return TRUE;
	}

	natc = p_atgrp_restr->GetNAtoms();
	if( p_mm_model->restr_ref_coords.size() != natc )
	{
		PrintLog("Error in MAmberMMModel::SetAtomPosRestrData() \n");
		PrintLog("Dimension of Reference Coordinates array is not equal \n");
		PrintLog("To the size of restrained atoms group \n");
		PrintLog("Atom Position Restraints will be turned off \n");
		natc = 0;
		return FALSE;

	}
			
	int i;
	AtomIntMap& atm_idx_map = p_mm_model->GetAtIdxMap(); 

	atm_xc.resize(3*natom);
	for(i=0; i < natom; i++)
	{
		HaAtom* aptr = p_mm_model->Atoms[i];
		atm_xc[3*i]   = aptr->GetX_Ang();
		atm_xc[3*i+1] = aptr->GetY_Ang();
		atm_xc[3*i+2] = aptr->GetZ_Ang();
	}

	atm_weight.resize(natc);
	atm_jrc.resize(natc);

	p_mm_model->SetAtomRestrForceConst(p_mm_model->atom_restr_const);

	for(i = 0; i < natc; i++)
	{
		HaAtom* aptr = (*p_atgrp_restr)[i];
		int idx_at = atm_idx_map[aptr];
		atm_jrc[i] = idx_at + 1;

		atm_xc[3*idx_at  ] = p_mm_model->restr_ref_coords[i].GetX();
		atm_xc[3*idx_at+1] = p_mm_model->restr_ref_coords[i].GetY();
		atm_xc[3*idx_at+2] = p_mm_model->restr_ref_coords[i].GetZ();
	}
	return TRUE;	
}

void AmberMMModel::Bcast(MPI_Comm& comm)
{
	int ires;
	int rank;
	int nsize;

//	PrintLog(" AmberMMModel::Bcast() pt 1 \n") ;

	ires = MPI_Comm_rank(comm,&rank);

	ires = MPI_Bcast(&ntf,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&es_cutoff,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&vdw_cutoff,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_mm_model->pme_grid_nx,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&p_mm_model->pme_grid_ny,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&p_mm_model->pme_grid_nz,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&p_mm_model->pme_spline_order,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&p_mm_model->vdw_correction_flag,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&p_mm_model->pme_verbose,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&p_mm_model->skin_nb,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_mm_model->pme_dsum_tol,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_mm_model->pme_ew_coeff,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_mm_model->pme_eedtbdns,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_mm_model->fft_grids_per_ang,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&dielc,1,MPI_DOUBLE,0,comm);	
	ires = MPI_Bcast(&scnb,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&scee,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&iamoeba,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&alpb,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&arad,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&intdiel,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&extdiel,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&saltcon,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&rgbmax,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&rbornstat,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&offset,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&gbsa,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&surften,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&cut_inner,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&gb_cutoff,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&bbox_xmin,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&bbox_ymin,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&bbox_zmin,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&bbox_xmax,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&bbox_ymax,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&bbox_zmax,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&ibelly,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&igb,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&gb_alpha,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&gb_beta,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&gb_gamma,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&gb_fs_max,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&gb_kappa,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&gb_neckscale,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&using_pme_potential,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&using_gb_potential, 1,MPI_INT,0,comm);
	ires = MPI_Bcast(&natom,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&ntypes,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&nbonh,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&ntheth,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&nphih,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&next,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&nres,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&nbona,1,MPI_INT,0,comm);	
	ires = MPI_Bcast(&ntheta,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&nphia,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&numbnd,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&numang,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&nptra,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&nphb,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&nspm,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&nttyp,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&bonda_idx,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&anglea_idx,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&diheda_idx,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&gbl_bond_allocsize,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&gbl_angle_allocsize,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&gbl_dihed_allocsize,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&belly_atm_cnt,1,MPI_INT,0,comm);	
	ires = MPI_Bcast(&max_res_size,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&n_solute_res,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&n_solute_mol,1,MPI_INT,0,comm);

	ires = MPI_Bcast(&do_amoeba_valence,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_amoeba_nonbond,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_bond,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_ureyb,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_reg_angle,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_trig_angle,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_opbend,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_torsion,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_pi_torsion,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_strbend,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_torsion_torsion,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_str_torsion,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_recip,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_adjust,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_direct,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_self,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_vdw,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_induced,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_vdw_taper,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&do_vdw_longrange,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&beeman_integrator,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&p_mm_model->dipole_scf_iter_max,1,MPI_INT,0,comm);
	ires = MPI_Bcast(&amoeba_verbose,1,MPI_INT,0,comm); 
	ires = MPI_Bcast(&p_mm_model->dipole_scf_tol,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_mm_model->ee_dsum_cut,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_mm_model->ee_damped_cut,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_mm_model->sor_coefficient,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_mm_model->thole_expon_coeff,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&p_mm_model->vdw_taper,1,MPI_DOUBLE,0,comm);

	if(rank == 0) nsize =  gbl_res_atms.size();
	ires = MPI_Bcast(&nsize,1,MPI_INT,0,comm);
	if(rank != 0) gbl_res_atms.resize(nsize);
	ires = MPI_Bcast(gbl_res_atms.v(),nsize,MPI_INT,0,comm);

//	res_labels - are not broadcasted - may be we should

	if(rank != 0) atm_nsp.resize(nspm);
	ires = MPI_Bcast(atm_nsp.v(),nspm,MPI_INT,0,comm);
	
	if(rank != 0) atm_iac.resize(natom);
	ires = MPI_Bcast(atm_iac.v(),natom,MPI_INT,0,comm);

	if(rank == 0) nsize =  p_amber_driver->title.size();
	ires = MPI_Bcast(&nsize,1,MPI_INT,0,comm);
	if(rank != 0) p_amber_driver->title.resize(nsize);
	ires = MPI_Bcast(&p_amber_driver->title[0],nsize,MPI_CHAR,0,comm);

	if(rank != 0) atm_charge.resize(natom);
	ires = MPI_Bcast(atm_charge.v(),natom,MPI_DOUBLE,0,comm);

	if(rank != 0) atm_mass.resize(natom);
	ires = MPI_Bcast(atm_mass.v(),natom,MPI_DOUBLE,0,comm);

	ires = MPI_Bcast(&tmass,1,MPI_DOUBLE,0,comm);

	if(rank != 0) atm_mass_inv.resize(natom);
	ires = MPI_Bcast(atm_mass_inv.v(),natom,MPI_DOUBLE,0,comm);

	if(rank != 0) atm_igraph.resize(natom);  // Not broadcasted yet
	if(rank != 0) atm_isymbl.resize(natom);  // Not broadcasted yet
	if(rank != 0) atm_itree.resize(natom);  // Not broadcasted yet

	if( rank == 0 ) 
	{
		nsize =  atm_igroup.size();
		if( nsize != natom )
		{
			atm_igroup.resize(natom);
			atm_igroup = 1;
		}
	}
	ires = MPI_Bcast(&nsize,1,MPI_INT,0,comm);
	if( rank != 0 ) atm_igroup.resize(nsize);
	if( nsize > 0 ) ires = MPI_Bcast(atm_igroup.v(),nsize,MPI_INT,0,comm);

	if(rank != 0) atm_gb_radii.resize(natom);
	ires = MPI_Bcast(atm_gb_radii.v(),natom,MPI_DOUBLE,0,comm);

	if(rank != 0) atm_gb_fs.resize(natom);
	ires = MPI_Bcast(atm_gb_fs.v(),natom,MPI_DOUBLE,0,comm);

	if(rank != 0) typ_ico.newsize(ntypes*ntypes);
	ires = MPI_Bcast(typ_ico.v(),ntypes*ntypes,MPI_INT,0,comm);

	if(rank != 0) gbl_cn1.resize(nttyp);
	ires = MPI_Bcast(gbl_cn1.v(),nttyp,MPI_DOUBLE,0,comm);

	if(rank != 0) gbl_cn2.resize(nttyp);
	ires = MPI_Bcast(gbl_cn2.v(),nttyp,MPI_DOUBLE,0,comm);

	if(rank != 0) gbl_asol.resize(nphb);
	if(nphb > 0) ires = MPI_Bcast(gbl_asol.v(),nphb,MPI_DOUBLE,0,comm);

	if(rank != 0) gbl_bsol.resize(nphb);
	if(nphb > 0) ires = MPI_Bcast(gbl_bsol.v(),nphb,MPI_DOUBLE,0,comm);

	if(rank != 0) gbl_rk.resize(numbnd);
	ires = MPI_Bcast(gbl_rk.v(),numbnd,MPI_DOUBLE,0,comm);

	if(rank != 0) gbl_req.resize(numbnd);
	ires = MPI_Bcast(gbl_req.v(),numbnd,MPI_DOUBLE,0,comm);

	if(rank != 0) gbl_tk.resize(numang);
	ires = MPI_Bcast(gbl_tk.v(),numang,MPI_DOUBLE,0,comm);

	if(rank != 0) gbl_teq.resize(numang);
	ires = MPI_Bcast(gbl_teq.v(),numang,MPI_DOUBLE,0,comm);

	if(rank != 0) gbl_pk.resize(nptra);
	ires = MPI_Bcast(gbl_pk.v(),nptra,MPI_DOUBLE,0,comm);

	if(rank != 0) gbl_phase.resize(nptra);
	ires = MPI_Bcast(gbl_phase.v(),nptra,MPI_DOUBLE,0,comm);

	if(rank != 0) gbl_pn.resize(nptra);
	ires = MPI_Bcast(gbl_pn.v(),nptra,MPI_DOUBLE,0,comm);

	if(rank != 0) gbl_gamc.resize(nptra);
	ires = MPI_Bcast(gbl_gamc.v(),nptra,MPI_DOUBLE,0,comm);

	if(rank != 0) gbl_gams.resize(nptra);
	ires = MPI_Bcast(gbl_gams.v(),nptra,MPI_DOUBLE,0,comm);

	if(rank != 0) gbl_ipn.resize(nptra);
	ires = MPI_Bcast(gbl_ipn.v(),nptra,MPI_INT,0,comm);

	if(rank != 0) gbl_fmn.resize(nptra);
	ires = MPI_Bcast(gbl_fmn.v(),nptra,MPI_DOUBLE,0,comm);

	int i;

// Broadcast Valence Bonds 

	if(rank != 0) gbl_bond.resize(3*gbl_bond_allocsize);
	ires = MPI_Bcast(gbl_bond.v(),3*gbl_bond_allocsize,MPI_INT,0,comm);

// Broadcast Valence Angles

	if(rank != 0) gbl_angle.resize(4*gbl_angle_allocsize);
	ires = MPI_Bcast(gbl_angle.v(),4*gbl_angle_allocsize,MPI_INT,0,comm);

// Broadcast dihedrals

	if(rank != 0) gbl_dihed.resize(5*gbl_dihed_allocsize);
	ires = MPI_Bcast(gbl_dihed.v(),5*gbl_dihed_allocsize,MPI_INT,0,comm);

// Nonbonded interations atom excluded list

	if(rank != 0) atm_numex.resize(natom);
	ires = MPI_Bcast(atm_numex.v(),natom,MPI_INT,0,comm);

	if(rank != 0) gbl_natex.resize(next);
	ires = MPI_Bcast(gbl_natex.v(),next,MPI_INT,0,comm);

// Atom Position restraints info 
  
    ires = MPI_Bcast(&natc,1,MPI_INT,0,comm);

	if(rank != 0) atm_jrc.resize(natc);
	ires = MPI_Bcast(atm_jrc.v(),natc,MPI_INT,0,comm);

	if(rank != 0) atm_weight.resize(natc);
	ires = MPI_Bcast(atm_weight.v(),natc,MPI_DOUBLE,0,comm);

	if(rank != 0) atm_xc.resize(3*natom);
	ires = MPI_Bcast(atm_xc.v(),3*natom,MPI_DOUBLE,0,comm);

// Broadcast Atom-Atom Distance Constraints:

	ires =  MPI_Bcast(&num_dist_constr,1,MPI_INT,0,comm);
	if( rank != 0 ) dist_constr_idx.resize(3*num_dist_constr);
	if( rank != 0 ) dist_constr_params.resize(2*num_dist_constr);
	
	ires = MPI_Bcast(dist_constr_idx.v(),   3*num_dist_constr,MPI_INT,   0,comm);
	ires = MPI_Bcast(dist_constr_params.v(),2*num_dist_constr,MPI_DOUBLE,0,comm);

	if( p_mm_model->IsAmoebaFF() )
	{
		ires = MPI_Bcast(&n_bond_amoeba,1,MPI_INT,0,comm);
		nsize = n_bond_amoeba*3;
		if(rank != 0) gbl_bond_amoeba.resize(nsize);
		ires = MPI_Bcast(gbl_bond_amoeba.v(),nsize,MPI_INT,0,comm);

		ires = MPI_Bcast(&n_bond_amoeba_params,1,MPI_INT,0,comm);
		if(rank != 0) bond_amoeba_params.resize(2);
		if(rank != 0) bond_amoeba_params[0].resize(n_bond_amoeba_params);
		if(rank != 0) bond_amoeba_params[1].resize(n_bond_amoeba_params);
		ires = MPI_Bcast(bond_amoeba_params[0].v(),n_bond_amoeba_params,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast(bond_amoeba_params[1].v(),n_bond_amoeba_params,MPI_DOUBLE,0,comm);

		ires = MPI_Bcast(&bond_amoeba_ftab_degree,1,MPI_INT,0,comm);
		nsize = bond_amoeba_ftab_degree + 1;
		if(rank != 0) bond_amoeba_ftab_coef.resize(nsize);
		ires = MPI_Bcast(bond_amoeba_ftab_coef.v(),nsize,MPI_DOUBLE,0,comm);
//
		ires = MPI_Bcast(&n_urey_bond,1,MPI_INT,0,comm);
		nsize = n_urey_bond*3;
		if(rank != 0) gbl_bond_urey.resize(nsize);
		ires = MPI_Bcast(gbl_bond_urey.v(),nsize,MPI_INT,0,comm);

		ires = MPI_Bcast(&n_urey_bond_params,1,MPI_INT,0,comm);
		if(rank != 0) bond_urey_params.resize(2);
		if(rank != 0) bond_urey_params[0].resize(n_urey_bond_params);
		if(rank != 0) bond_urey_params[1].resize(n_urey_bond_params);
		ires = MPI_Bcast(bond_urey_params[0].v(),n_urey_bond_params,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast(bond_urey_params[1].v(),n_urey_bond_params,MPI_DOUBLE,0,comm);

		ires = MPI_Bcast(&bond_urey_ftab_degree,1,MPI_INT,0,comm);
		nsize = bond_urey_ftab_degree + 1;
		if(rank != 0) bond_urey_ftab_coef.resize(nsize);
		ires = MPI_Bcast(bond_urey_ftab_coef.v(),nsize,MPI_DOUBLE,0,comm);

// 
		ires = MPI_Bcast(&n_angle_amoeba,1,MPI_INT,0,comm);
		nsize = n_angle_amoeba*4;
		if(rank != 0) gbl_angle_amoeba_reg.resize(nsize);
		ires = MPI_Bcast( gbl_angle_amoeba_reg.v(),nsize,MPI_INT,0,comm );

		ires = MPI_Bcast(&n_angle_amoeba_params,1,MPI_INT,0,comm);
		if(rank != 0) angle_amoeba_params.resize(2);
		if(rank != 0) angle_amoeba_params[0].resize(n_angle_amoeba_params);
		if(rank != 0) angle_amoeba_params[1].resize(n_angle_amoeba_params);
		ires = MPI_Bcast(angle_amoeba_params[0].v(),n_angle_amoeba_params,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast(angle_amoeba_params[1].v(),n_angle_amoeba_params,MPI_DOUBLE,0,comm);

		ires = MPI_Bcast(&angle_amoeba_ftab_degree,1,MPI_INT,0,comm);
		nsize = angle_amoeba_ftab_degree + 1;
		if(rank != 0) angle_amoeba_ftab_coef.resize(nsize);
		ires = MPI_Bcast( angle_amoeba_ftab_coef.v(),nsize,MPI_DOUBLE,0,comm);
//
		ires = MPI_Bcast(&n_trig_angles,1,MPI_INT,0,comm);
		nsize = n_trig_angles*5;
		if(rank != 0) gbl_angle_amoeba_trig.resize(nsize);
		ires = MPI_Bcast( gbl_angle_amoeba_trig.v(),nsize,MPI_INT,0,comm );
//
		ires = MPI_Bcast(&n_opbend_angles,1,MPI_INT,0,comm);
		nsize = n_opbend_angles*5;
		if(rank != 0) gbl_opbend_angle.resize(nsize);
		ires = MPI_Bcast( gbl_opbend_angle.v(),nsize,MPI_INT,0,comm );

		ires = MPI_Bcast(&n_opbend_angles_params,1,MPI_INT,0,comm);
		if(rank != 0) opbend_angle_params.resize(n_opbend_angles_params);
		ires = MPI_Bcast(opbend_angle_params.v(),n_opbend_angles_params,MPI_DOUBLE,0,comm);
//
		ires = MPI_Bcast(&n_tors_amoeba,1,MPI_INT,0,comm);
		nsize = n_tors_amoeba*5;
		if(rank != 0) gbl_amoeba_tors_angle.resize(nsize);
		ires = MPI_Bcast( gbl_amoeba_tors_angle.v(),nsize,MPI_INT,0,comm );

		ires = MPI_Bcast(&n_tors_amoeba_params,1,MPI_INT,0,comm);
		if(rank != 0) tors_amoeba_params.resize(3);
		if(rank != 0) tors_amoeba_params[0].resize(n_tors_amoeba_params);
		if(rank != 0) tors_amoeba_params[1].resize(n_tors_amoeba_params);
		if(rank != 0) tors_amoeba_params[2].resize(n_tors_amoeba_params);
		ires = MPI_Bcast(tors_amoeba_params[0].v(),n_tors_amoeba_params,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast(tors_amoeba_params[1].v(),n_tors_amoeba_params,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast(tors_amoeba_params[2].v(),n_tors_amoeba_params,MPI_DOUBLE,0,comm);

		ires = MPI_Bcast(&n_pi_torsions,1,MPI_INT,0,comm);
		nsize = n_pi_torsions*7;
		if(rank != 0) gbl_pi_tors_angle.resize(nsize);
		ires = MPI_Bcast( gbl_pi_tors_angle.v(),nsize,MPI_INT,0,comm );

		ires = MPI_Bcast(&n_pi_torsions_params,1,MPI_INT,0,comm);
		if(rank != 0) pi_tors_params.resize(3);
		if(rank != 0) pi_tors_params[0].resize(n_pi_torsions_params);
		if(rank != 0) pi_tors_params[1].resize(n_pi_torsions_params);
		if(rank != 0) pi_tors_params[2].resize(n_pi_torsions_params);
		ires = MPI_Bcast(pi_tors_params[0].v(),n_pi_torsions_params,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast(pi_tors_params[1].v(),n_pi_torsions_params,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast(pi_tors_params[2].v(),n_pi_torsions_params,MPI_DOUBLE,0,comm);
//
		ires = MPI_Bcast(&n_stretch_bend,1,MPI_INT,0,comm);
		nsize = n_stretch_bend*4;
		if(rank != 0) gbl_str_bend_angle.resize(nsize);
		ires = MPI_Bcast( gbl_str_bend_angle.v(),nsize,MPI_INT,0,comm );

		ires = MPI_Bcast(&n_stretch_bend_params,1,MPI_INT,0,comm);
		if(rank != 0) str_bend_params.resize(4);
		if(rank != 0) str_bend_params[0].resize(n_stretch_bend_params);
		if(rank != 0) str_bend_params[1].resize(n_stretch_bend_params);
		if(rank != 0) str_bend_params[2].resize(n_stretch_bend_params);
		if(rank != 0) str_bend_params[3].resize(n_stretch_bend_params);
		ires = MPI_Bcast(str_bend_params[0].v(),n_stretch_bend_params,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast(str_bend_params[1].v(),n_stretch_bend_params,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast(str_bend_params[2].v(),n_stretch_bend_params,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast(str_bend_params[3].v(),n_stretch_bend_params,MPI_DOUBLE,0,comm);
//
        ires = MPI_Bcast(&n_tors_tors,1,MPI_INT,0,comm);
		nsize = n_tors_tors*6;
		if(rank != 0) gbl_tors_tors.resize(nsize);
		ires = MPI_Bcast( gbl_tors_tors.v(),nsize,MPI_INT,0,comm );

		ires = MPI_Bcast(&n_tors_tors_params,1,MPI_INT,0,comm);
		
		if(rank != 0) tors_tors_id_params.resize(n_tors_tors_params);
		ires = MPI_Bcast(tors_tors_id_params.v(),n_tors_tors_params,MPI_INT,0,comm);

		if(rank != 0) tors_tors_params.resize(n_tors_tors_params);

		for(i = 0; i < n_tors_tors_params; i++ )
		{
			int dim1;
			int dim2;
			if( rank == 0) dim1 = tors_tors_params[i][0].size();
			if( rank == 0) dim2 = tors_tors_params[i][1].size();
			ires = MPI_Bcast(&dim1,1,MPI_INT,0,comm);
			ires = MPI_Bcast(&dim2,1,MPI_INT,0,comm);
			if(rank != 0) tors_tors_params[i].resize(6);
			if(rank != 0) tors_tors_params[i][0].resize(dim1);
			if(rank != 0) tors_tors_params[i][1].resize(dim2);
			if(rank != 0) tors_tors_params[i][2].resize(dim1*dim2);
			if(rank != 0) tors_tors_params[i][3].resize(dim1*dim2);
			if(rank != 0) tors_tors_params[i][4].resize(dim1*dim2);
			if(rank != 0) tors_tors_params[i][5].resize(dim1*dim2);

			ires = MPI_Bcast(tors_tors_params[i][0].v(),dim1,MPI_DOUBLE,0,comm);
			ires = MPI_Bcast(tors_tors_params[i][1].v(),dim2,MPI_DOUBLE,0,comm);
			ires = MPI_Bcast(tors_tors_params[i][2].v(),dim1*dim2,MPI_DOUBLE,0,comm);
			ires = MPI_Bcast(tors_tors_params[i][3].v(),dim1*dim2,MPI_DOUBLE,0,comm);
			ires = MPI_Bcast(tors_tors_params[i][4].v(),dim1*dim2,MPI_DOUBLE,0,comm);
			ires = MPI_Bcast(tors_tors_params[i][5].v(),dim1*dim2,MPI_DOUBLE,0,comm);
		}
//
		if(rank != 0) atm_amoeba_vdw_type.resize(natom);	
		if(rank != 0) atm_parent_id.resize(natom);
		if(rank != 0) atm_parent_weight.resize(natom);

		ires = MPI_Bcast( atm_amoeba_vdw_type.v(),natom,MPI_INT,0,comm );
		ires = MPI_Bcast( atm_parent_id.v(),natom,MPI_INT,0,comm );
		ires = MPI_Bcast( atm_parent_weight.v(),natom,MPI_DOUBLE,0,comm);

		ires = MPI_Bcast(&vdw_buffer_delta,1,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast(&vdw_buffer_gamma,1,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast(&n_vdw_params,1,MPI_INT,0,comm);

		if(rank != 0) amoeba_vdw_rstars.resize(n_vdw_params*n_vdw_params);
		if(rank != 0) amoeba_vdw_depths.resize(n_vdw_params*n_vdw_params);

		ires = MPI_Bcast( amoeba_vdw_rstars.v(),n_vdw_params*n_vdw_params,MPI_DOUBLE,0,comm);
		
		ires = MPI_Bcast( amoeba_vdw_depths.v(),n_vdw_params*n_vdw_params,MPI_DOUBLE,0,comm);
//
		ires = MPI_Bcast(&num_local_multipoles,1,MPI_INT,0,comm);
		ires = MPI_Bcast(&num_chiral_frames,1,MPI_INT,0,comm);
		ires = MPI_Bcast(&num_reg_frames,1,MPI_INT,0,comm);

		if(rank != 0) atm_multipoles.resize(10*num_local_multipoles);
		ires = MPI_Bcast(atm_multipoles.v(),10*num_local_multipoles,MPI_DOUBLE,0,comm);
		
		if(rank != 0) atm_chiral_frames.resize(3*num_chiral_frames);
		ires = MPI_Bcast( atm_chiral_frames.v(),3*num_chiral_frames,MPI_INT,0,comm );
		
		if(rank != 0) atm_reg_frames.resize(5*num_reg_frames);
		ires = MPI_Bcast( atm_reg_frames.v(),5*num_reg_frames,MPI_INT,0,comm );
//
		ires = MPI_Bcast(&num_adjust_list,1,MPI_INT,0,comm);
		if(rank != 0) atm_adjust_list.resize( num_adjust_list*3);
		ires = MPI_Bcast( atm_adjust_list.v(),num_adjust_list*3,MPI_INT,0,comm );

		if(rank != 0) adjust_vdw_weights.resize(9);
		if(rank != 0) adjust_mpole_weights.resize(9);
		if(rank != 0) adjust_direct_weights.resize(9);
		if(rank != 0) adjust_polar_weights.resize(9);
		if(rank != 0) adjust_mutual_weights.resize(9);

		ires = MPI_Bcast( adjust_vdw_weights.v(),9,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast( adjust_mpole_weights.v(),9,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast( adjust_direct_weights.v(),9,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast( adjust_polar_weights.v(),9,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast( adjust_mutual_weights.v(),9,MPI_DOUBLE,0,comm);

		if( rank != 0) atm_polar.resize(natom);
		if( rank != 0) atm_hpolar.resize(natom);
		if( rank != 0) atm_screen_polar.resize(natom);
		if( rank != 0) atm_qterm.resize(natom);
		if( rank != 0) is_polarizable.resize(natom);
		if( rank != 0) damp_polar_strength.resize(natom);
		if( rank != 0) damp_polar_sensitivity.resize(natom);
		if( rank != 0) damp_polar_rad.resize(natom);
		
		ires = MPI_Bcast( atm_polar.v(),natom,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast( atm_hpolar.v(),natom,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast( atm_screen_polar.v(),natom,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast( atm_qterm.v(),natom,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast( is_polarizable.v(),natom,MPI_INT,0,comm);
		ires = MPI_Bcast( damp_polar_strength.v(),natom,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast( damp_polar_sensitivity.v(),natom,MPI_DOUBLE,0,comm);
		ires = MPI_Bcast( damp_polar_rad.v(),natom,MPI_DOUBLE,0,comm);

	}
//	PrintLog(" AmberMMModel::Bcast() pt end \n") ;
}

void AmberMMModel::BcastAtmMass(MPI_Comm& comm)
{
	int ires;
	int rank;

	ires = MPI_Comm_rank(comm,&rank);

	if(rank != 0) atm_mass.resize(natom);
	if(rank != 0) atm_mass_inv.resize(natom);
	ires = MPI_Bcast(atm_mass.v(),natom,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(&tmass,1,MPI_DOUBLE,0,comm);
	ires = MPI_Bcast(atm_mass_inv.v(),natom,MPI_DOUBLE,0,comm);
}

void AmberMMModel::Clear()
{
	natom = 0; 
	natc = 0;   // contraints and belly atoms are not set yet
	belly_atm_cnt = 0;
	ibelly = 0;

	iamoeba = 0;

	num_deg         = 0;
	num_deg_solute  = 0; 
	num_deg_solvent = 0;

	num_dist_constr = 0;
	dist_constr_idx.clear();
	dist_constr_params.clear();

	do_amoeba_valence = 1;
	do_amoeba_nonbond = 1;
	do_bond = 1;
	do_ureyb = 1;
	do_reg_angle = 1;
	do_trig_angle = 1;
	do_opbend = 1;
	do_torsion = 1;
	do_pi_torsion = 1;
	do_strbend = 1;
	do_torsion_torsion = 1;
	do_str_torsion = 1;
	do_recip = 1;
	do_adjust = 1;
	do_direct = 1;
	do_self = 1;
	do_vdw = 1;
	do_induced = 1;

	do_vdw_taper = 1;
	do_vdw_longrange = 1;
	beeman_integrator = 0;
	
	amoeba_verbose = 0;

	to_update_amber_data = TRUE;
}

void AmberMMModel::SetUpdateDataFlag(int to_update_flag_new )
{
	to_update_amber_data = to_update_flag_new;
}

int AmberMMModel::SetMovingAtomsData()
{
	atm_igroup.resize(natom); 
	atm_igroup = 1;
	belly_atm_cnt = 0;
	ibelly = 0;

	AtomGroup* p_atgrp_mv =  p_mm_model->GetMovingAtoms();
	if( p_atgrp_mv == NULL) return TRUE;

	AtomIteratorAtomGroup aitr(p_atgrp_mv);

	std::set<HaAtom*, less<HaAtom*> > atset_mv;
	HaAtom* aptr;
	for( aptr = aitr.GetFirstAtom(); aptr; aptr = aitr.GetNextAtom() )
	{
		atset_mv.insert(aptr);
	}

	int i;
	for(i = 0; i < natom; i++)
	{
		aptr = p_mm_model->Atoms[i];
		if( atset_mv.count(aptr) > 0) 
		{
			atm_igroup[i] = 1;
			belly_atm_cnt++;
		}
	}
	if( belly_atm_cnt > 0 ) ibelly = 1;
	
	//int idx = 0;
	//int nbonh_new = nbonh;
	//int nbona_new = nbona;
	//for(i = 0; i < nbonh + nbona; i++)
	//{
	//	int iat1 = gbl_bond[3*i  ] - 1;
	//	int iat2 = gbl_bond[3*i+1] - 1;
	//	if( atm_igroup[iat1] > 0 || atm_igroup[iat2] > 0 )
	//	{
	//		gbl_bond[3*idx]   = gbl_bond[3*i];
	//		gbl_bond[3*idx+1] = gbl_bond[3*i+1];
	//		gbl_bond[3*idx+2] = gbl_bond[3*i+2];
	//		idx++;
	//	}
	//	else
	//	{
	//		gbl_bond_allocsize--;
	//		if( i < nbonh )
	//		{
	//			nbonh_new--;
	//			bonda_idx--;
	//		}
	//		else
	//		{
	//			nbona_new--;
	//		}
	//	}
	//}
	//nbonh = nbonh_new;
	//nbona = nbona_new;

	//idx = 0;
	//int ntheth_new = ntheth;
	//int ntheta_new = ntheta;
	//for(i = 0; i < ntheth + ntheta; i++)
	//{
	//	int iat1 = gbl_angle[4*i  ] - 1;
	//	int iat2 = gbl_angle[4*i+1] - 1;
	//	int iat3 = gbl_angle[4*i+2] - 1;
	//	if( atm_igroup[iat1] > 0 || atm_igroup[iat2] > 0 || atm_igroup[iat3] > 0)
	//	{
	//		gbl_angle[4*idx  ]   = gbl_angle[4*i  ];
	//		gbl_angle[4*idx+1]   = gbl_angle[4*i+1];
	//		gbl_angle[4*idx+2]   = gbl_angle[4*i+2];
	//		gbl_angle[4*idx+3]   = gbl_angle[4*i+3];
	//		idx++;
	//	}
	//	else
	//	{
	//		gbl_angle_allocsize--;
	//		if( i < ntheth )
	//		{
	//			ntheth_new--;
	//			anglea_idx--;
	//		}
	//		else
	//		{
	//			ntheta_new--;
	//		}
	//	}
	//}
	//ntheth = ntheth_new;
	//ntheta = ntheta_new;

	//idx = 0;
	//int nphih_new = nphih;
	//int nphia_new = nphia;
	//for(i = 0; i < nphih + nphia; i++)
	//{
	//	int iat1 = gbl_dihed[5*i  ] - 1;
	//	int iat2 = gbl_dihed[5*i+1] - 1;
	//	int iat3 = gbl_dihed[5*i+2] - 1;
	//	int iat4 = gbl_dihed[5*i+3] - 1;
	//	if( atm_igroup[iat1] > 0 || atm_igroup[iat2] > 0 || atm_igroup[iat3] > 0 || atm_igroup[iat4] )
	//	{
	//		gbl_dihed[5*idx  ]   = gbl_dihed[5*i  ];
	//		gbl_dihed[5*idx+1]   = gbl_dihed[5*i+1];
	//		gbl_dihed[5*idx+2]   = gbl_dihed[5*i+2];
	//		gbl_dihed[5*idx+3]   = gbl_dihed[5*i+3];
	//		gbl_dihed[5*idx+4]   = gbl_dihed[5*i+4];
	//		idx++;
	//	}
	//	else
	//	{
	//		gbl_dihed_allocsize--;
	//		if( i < nphih )
	//		{
	//			nphih_new--;
	//			diheda_idx--;
	//		}
	//		else
	//		{
	//			nphia_new--;
	//		}
	//	}
	//}
	//nphih = nphih_new;
	//nphia = nphia_new;

	return TRUE;
}

int AmberMMModel::SetDistConstrData()
{
	num_dist_constr  = p_mm_model->DistConstraints.size();
	dist_constr_idx.resize(3*num_dist_constr);     
	dist_constr_params.resize(2*num_dist_constr);
	int i;
	int ii = 0;
	int jj = 0;
	
	AtomIntMap& at_idx_map = p_mm_model->GetAtIdxMap();

	for( i = 0; i < num_dist_constr; i++ )
	{
		int i1 = at_idx_map[p_mm_model->DistConstraints[i].pt1] + 1;
		int i2 = at_idx_map[p_mm_model->DistConstraints[i].pt2] + 1;

		dist_constr_idx[ii] = i1; ii++;
		dist_constr_idx[ii] = i2; ii++;
		dist_constr_idx[ii] = p_mm_model->DistConstraints[i].cnt_type; ii++;

		dist_constr_params[jj] = p_mm_model->DistConstraints[i].cf[0]; jj++;
		dist_constr_params[jj] = p_mm_model->DistConstraints[i].cf[1]; jj++;
	}
	return TRUE;
}

int AmberMMModel::AddAtomFrames(HaAtom* aptr, AtomFFParam* p_at_ff, StrAtomMap* p_templ_atname_to_res_map ) 
{
	std::string at_name = aptr->GetName();
	if( !p_at_ff->HasFrameAtomNames() ) 
	{
		PrintLog(" FF template for atom %s  does not have a frame \n", at_name.c_str() );
		return FALSE;
	}
	
	AtomIntMap& atoms_idx = p_mm_model->GetAtIdxMap();

	HaResidue* pres = aptr->GetHostRes();
	std::string res_fname = pres->GetFullName();

	int nat_fr = p_at_ff->frame_atom_names.size();						

	std::vector<int> atoms_idx_fr(nat_fr);
	int j;
	for(j = 0; j < nat_fr; j++)
	{
		std::string at_name_fr = p_at_ff->frame_atom_names[j];
		if( p_templ_atname_to_res_map->count(at_name_fr) == 0) throw std::runtime_error(" No atom of the residue template " + res_fname +  " with name " + at_name_fr + " is mapped to atoms of the system " );
		HaAtom* aptr_fr = (*p_templ_atname_to_res_map)[at_name_fr];
		int at_idx_fr = atoms_idx[aptr_fr];
		at_idx_fr++;
		atoms_idx_fr[j] = at_idx_fr;
	}

	if( nat_fr == 3 )
	{
		int id1 = atoms_idx_fr[0];
		int id2 = atoms_idx_fr[1];
		int id3 = atoms_idx_fr[2];

		if( p_at_ff->IsBisectFrame() )
		{
			atm_reg_frames.push_back( id1 );
			atm_reg_frames.push_back(  1  );
			atm_reg_frames.push_back( id1 );
			atm_reg_frames.push_back( id2 );
			atm_reg_frames.push_back(  2  );

			atm_reg_frames.push_back( id1 );
			atm_reg_frames.push_back(  1  );
			atm_reg_frames.push_back( id1 );
			atm_reg_frames.push_back( id3 );
			atm_reg_frames.push_back(  2  );

			atm_reg_frames.push_back( id1 );
			atm_reg_frames.push_back(  2  );
			atm_reg_frames.push_back( id1 );
			atm_reg_frames.push_back( id3 );
			atm_reg_frames.push_back(  1  );
		}
		else
		{
			atm_reg_frames.push_back( id1 );
			atm_reg_frames.push_back(  1  );
			atm_reg_frames.push_back( id1 );
			atm_reg_frames.push_back( id2 );
			atm_reg_frames.push_back(  1  );

			atm_reg_frames.push_back( id1 );
			atm_reg_frames.push_back(  2  );
			atm_reg_frames.push_back( id1 );
			atm_reg_frames.push_back( id3 );
			atm_reg_frames.push_back(  1  );
		}
	}
	else if( nat_fr == 4 )
	{
		int id1 = atoms_idx_fr[0];
		int id2 = atoms_idx_fr[1];
		int id3 = atoms_idx_fr[2];
		int id4 = atoms_idx_fr[3];

		if( !p_at_ff->IsChiralFrame() ) throw std::runtime_error(" Non Chiral frame with 4 atoms for atom name " + at_name + " In residue " + res_fname );
		atm_reg_frames.push_back( id1 );
		atm_reg_frames.push_back(  1  );
		atm_reg_frames.push_back( id1 );
		atm_reg_frames.push_back( id2 );
		atm_reg_frames.push_back(  1  );

		atm_reg_frames.push_back( id1 );
		atm_reg_frames.push_back(  2  );
		atm_reg_frames.push_back( id1 );
		atm_reg_frames.push_back( id3 );
		atm_reg_frames.push_back(  1  );

		atm_chiral_frames.push_back( id1 );
		atm_chiral_frames.push_back( id4 );
		atm_chiral_frames.push_back( 1 );
	}
	else
	{
		throw std::runtime_error(" The number of atoms specifying frame is not 3 or 4 for atom name " + at_name + " In residue " + res_fname );
	}
	return TRUE;
}

class AmberProcess : public wxProcess
{
public:
	AmberProcess() { p_mm_mod = NULL; }

	HaMolMechMod* p_mm_mod;

	virtual
		void OnTerminate(int pid, int status)
	{
		if (p_mm_mod) p_mm_mod->StopCalc();
		//		::wxMessageBox("AMBER Process Has Stopped \n");

		PrintLog("AMBER Process Has Stopped \n");
	}
};


int MMDriverAmber::RunAmberProg(int sync)
{
	int ires = DeleteOutputFiles();

	int i;
	StrVec sander_args;

	sander_args.push_back("-O");
	sander_args.push_back("-i"); 
	sander_args.push_back( amber_inp_file );
	sander_args.push_back("-o");
	sander_args.push_back( amber_out_file );
	sander_args.push_back("-p");
	sander_args.push_back( amber_top_file );
	sander_args.push_back("-c");
	sander_args.push_back( amber_crd_file );
	sander_args.push_back("-r");
	sander_args.push_back( amber_rst_file );
	sander_args.push_back("-ref");
	sander_args.push_back( amber_constr_crd_file );
	sander_args.push_back("-x");
	sander_args.push_back( amber_trj_coord_file );
	sander_args.push_back("-v");
	sander_args.push_back( amber_trj_vel_file );
	sander_args.push_back("-e");
	sander_args.push_back( amber_trj_ene_file );

	std::string exe_fname;

	/* if( p_mm_mod->ext_mm_prog == p_mm_mod->ext_mm_prog.PMEMD_9 )
	{
		exe_fname = "pmemd";
	}
	else if( p_mm_mod->ext_mm_prog == p_mm_mod->ext_mm_prog.SANDER_9 )
	{
		exe_fname = "sander";
	} 
	if (p_mm_mod->ext_mm_prog == p_mm_mod->ext_mm_prog.PMEMD_10)
	{
		exe_fname = "pmemd_amba";
	} */

	if (p_mm_mod->ext_mm_prog == p_mm_mod->ext_mm_prog.PMEMD_12 )
	{
		exe_fname = "pmemd_amba";
	}

	if (p_mm_mod->ext_mm_prog == p_mm_mod->ext_mm_prog.PMEMD_18 )
	{
		exe_fname = "pmemd";
	}

#if defined(_MSC_VER)
	exe_fname = pApp->harlem_home_dir + boost::filesystem::path::preferred_separator + exe_fname;
#else
	std::string exe_fname_test = pApp->harlem_home_dir + boost::filesystem::path::preferred_separator + "bin" + boost::filesystem::path::preferred_separator + exe_fname;
	PrintLog("exe_fname_test = %s \n", exe_fname_test.c_str() );
	if (boost::filesystem::exists(exe_fname_test)) exe_fname = exe_fname_test;
	PrintLog("exe_fname_test = %s \n", exe_fname.c_str() );
#endif

	std::string cmd_line = exe_fname.c_str();
	for( i =0; i < sander_args.size(); i++)
	{
		cmd_line += " ";
		cmd_line += sander_args[i];
	}

	PrintLog(" MMDriverAmber::RunAmberProg()  cmd_line:\n");
	PrintLog(" %s \n", cmd_line.c_str() );

	AmberProcess* p_sander_proc = new AmberProcess();
	p_sander_proc->p_mm_mod = p_mm_mod;

	// namespace bp = ::boost::process;

	int res;
	if( sync )
	{
//		bp::context ctx;
//		ctx.stdout_behavior = bp::silence_stream();
//		p_mm_mod->ext_proc_id = bp::launch(exec, args, ctx);

		res = wxExecute(cmd_line,wxEXEC_SYNC,p_sander_proc);
	}
	else
	{
		res = wxExecute(cmd_line,wxEXEC_ASYNC,p_sander_proc);
	}
    
	p_mm_mod->ext_proc_id = res;
	return TRUE;
}

int MolMechModel::UpdateConstraints_2()
{
	// Broadcast Atom-Atom Distance Constraints:
//	PrintLog(" MolMechModel::UpdateConstraints_2 () pt 1 \n");

	if(  p_mm_mod->p_amber_driver->numtasks > 1 )
	{
		MPI_Comm comm = p_mm_mod->p_amber_driver->driver_mpi_comm;
	
		int ires;
		int rank;
		ires = MPI_Comm_rank(comm,&rank);

		ires =  MPI_Bcast(&p_amber_model->num_dist_constr,1,MPI_INT,0,comm);

		if( rank != 0 ) p_amber_model->dist_constr_idx.resize(3*p_amber_model->num_dist_constr);
		if( rank != 0 ) p_amber_model->dist_constr_params.resize(2*p_amber_model->num_dist_constr);
	
		ires = MPI_Bcast(p_amber_model->dist_constr_idx.v(),   3*p_amber_model->num_dist_constr,MPI_INT,   0,comm);
		ires = MPI_Bcast(p_amber_model->dist_constr_params.v(),2*p_amber_model->num_dist_constr,MPI_DOUBLE,0,comm);
	}
	FC_FUNC_MODULE(dist_constr_mod,dist_constr_alloc)( &p_amber_model->num_dist_constr );

	if( p_amber_model->num_dist_constr > 0 )               
	{                                                          
		for( int i = 0; i < p_amber_model->num_dist_constr; i++)
		{
			int ip1 = i+1;
			set1_gbl_dist_constr_(&ip1,&p_amber_model->dist_constr_idx[3*i],&p_amber_model->dist_constr_idx[3*i+1],
				&p_amber_model->dist_constr_params[2*i],&p_amber_model->dist_constr_params[2*i+1],&p_amber_model->dist_constr_idx[3*i+2]);
		}			
	}
	HaVec_int use_atm_map_loc(  this->GetNA(), 1 );

	FC_FUNC_MODULE(dist_constr_mod,dist_constr_setup)( use_atm_map_loc.v(), p_mm_mod->p_amber_driver->gbl_atm_owner_map.v() );

	return TRUE;
}


TimerAmber::TimerAmber(MMDriverAmber* p_mm_driver_new)
{
	p_mm_driver = p_mm_driver_new; 
	p_amber_model = p_mm_driver->p_amber_model;

	detailed_timing_flag = FALSE;

	run_start_cputime = 0;       
    run_setup_end_cputime = 0;
    run_end_cputime = 0;
    run_start_walltime = 0;
    run_setup_end_walltime = 0;
    run_end_walltime = 0;
}

TimerAmber::~TimerAmber()
{


}

void TimerAmber::InitTimers()
{
#	if defined(WITH_LIB_PMEMD)
	run_start_cputime  = ((double)clock())/CLOCKS_PER_SEC;
	run_start_walltime = ((double)time(NULL));
	if(detailed_timing_flag)
	{
		FC_FUNC_MODULE(timers_mod,init_test_timers)();   // For detailed mpi performance monitoring
	}
#	endif
}

void TimerAmber::EndSetupTimers()
{
	run_setup_end_cputime  = ((double)clock())/CLOCKS_PER_SEC;
	run_setup_end_walltime = ((double)time(NULL));
}

void TimerAmber::EndRunTimers()
{
	run_end_cputime  = ((double)clock())/CLOCKS_PER_SEC;
	run_end_walltime = ((double)time(NULL));
}

void TimerAmber::PrintTimings()
{
#	if defined(WITH_LIB_PMEMD)
	if(p_mm_driver->master)
	{
		PrintLog("-----------------\n 5.  TIMINGS \n----------------\n");

		FC_FUNC_MODULE(timers_mod,run_end_cputime)      = ((double) run_end_cputime)/CLOCKS_PER_SEC;
		FC_FUNC_MODULE(timers_mod,run_setup_end_cputime) = ((double) run_setup_end_cputime)/CLOCKS_PER_SEC;

		PrintLog("|  Setup CPU time:  %11.2f  seconds \n",  (run_setup_end_cputime - run_start_cputime)); 
        PrintLog("|  NonSetup CPU time: %11.2f seconds \n", (run_end_cputime - run_setup_end_cputime));
        PrintLog("|  Total CPU time: %11.2f seconds  %9.2f hours \n", 
			 (run_end_cputime - run_start_cputime), (run_end_cputime - run_start_cputime)/3600.0);
        PrintLog("|  Setup wall time: %11.2f  seconds \n", (run_setup_end_walltime - run_start_walltime));
		PrintLog("|  NonSetup wall time: %11.2f seconds\n",(run_end_walltime - run_setup_end_walltime));
		PrintLog("|  Total wall time: %11.2f seconds  %9.2f hours \n",
            (run_end_walltime - run_start_walltime),(run_end_walltime - run_start_walltime)/3600.0);
	}
	
	if(detailed_timing_flag)
	{
	   FC_FUNC_MODULE(timers_mod,profile_cpu)(&p_mm_driver->p_mm_mod->run_type.value(), &p_amber_model->igb);
	   if( p_mm_driver->numtasks > 1 )
	   {
		  FC_FUNC_MODULE(timers_mod,print_test_timers_mpi)();    // Debugging output of performance timings
	   }
	   else
	   {
		  FC_FUNC_MODULE(timers_mod,print_test_timers)();
	   }
	}
#	endif
}

void TimerAmber::ZeroTime()
{
#	if defined(WITH_LIB_PMEMD)
	FC_FUNC_MODULE(timers_mod,zero_time)();
#	endif
}

void TimerAmber::UpdateTime(int itimer_type)
{
#	if defined(WITH_LIB_PMEMD)
	int itype = itimer_type;
	FC_FUNC_MODULE(timers_mod,update_time)(&itype);
#	endif
}


void MMDriverAmber::TestFFT1()
{
	HaVec_double data;
	HaVec_double data_save;

	int nsize = 16;
	data.resize(4*nsize);

	data = 0.0;
	int i;
	for(i = 0; i < nsize; i++)
	{
		data[i] = (double)i;
	}

	data_save = data;

	FC_FUNC_MODULE(fft1d_mod,fft1d_create_test)(&nsize);
	
	FC_FUNC_MODULE(fft1d_mod,fft1d_forward_test)(data.v());

	PrintLog("Data after forward transformation : \n");
	for(i = 0; i < 2*nsize; i++)
	{
		PrintLog(" data(%3d) = %16.9f \n",i,data[i]);
	}

	FC_FUNC_MODULE(fft1d_mod,fft1d_back_test)(data.v());

	FC_FUNC_MODULE(fft1d_mod,fft1d_destroy_test)();

	PrintLog("\n\nData after backward transformation : \n");
	for(i = 0; i < 2*nsize; i++)
	{
		PrintLog(" data(%3d) = %16.9f \n",i,data[i]);
	}
}
