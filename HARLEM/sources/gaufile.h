/*! \file gaufile.h
     
	Classes to interact with Gaussian RWF Files

    \author Igor Kurnikov  
    \date 1997-2003

*/
#ifndef GAUFILE_H
#define GAUFILE_H

#include "g94_globals.h"
#include "g94_protos.h" 

#include "hastring.h"

class GauFile
//! Class to represent Gaussian RWF file and its IO operations:
{
public:
  GauFile();
  GauFile(const char* NewFname, const int NewIunit, const char* mode="unknown");
  ~GauFile(){}
  int opendef();
  int open();
  int ifopen() { return flag_open; }
  int close(char* mode="keep");
  void set_gau_iunit(int unit) {iunit = unit ;}
  void set_alloc_unit(int size) { init_alloc=size ; }
  void set_file_name(const char* cfname) { fname=cfname ; }  
  int set_open_mode(const char* mode);
  int fileio(const char* operation, const int sub_file, 
              int length, void* target, int position);
  static void dump_info();
  
private:
	std::string fname;
  int open_mode;
  int flag_open;
  int iunit;
  int init_alloc;

public:
// Gaussian unit numbers:

	enum GAUSS_RWF_NUM { IO_gen = 501,      /* /GEN/  see description in L1.F */
                          IO_basis = 506,    /* /B/ Basis set information */ 
                          IO_label = 502,
                          IO_ibf   = 508,
                          IO_s     = 514,    /* Overlap matrix */
                          IO_hcore = 515,
                          IO_kinet = 516,
                          IO_fermi = 517,
                          IO_multipole = 518,
                          IO_emo       = 522,
                          IO_ca        = 524,
                          IO_mo_real_a = 524,
                          IO_mo_imag_a = 525,
                          IO_cb        = 526,
                          IO_mo_real_b = 526,
                          IO_mo_imag_b = 527,
                          IO_greenf    = 527,
                          IO_scf_den_real_a = 528,
                          IO_scf_den_imag_a = 529,
                          IO_scf_den_real_b = 530,
						  IO_scf_den_imag_b    =  531,
						IO_scf_den_real_tot  =  532,
						IO_scf_den_imag_tot  =  533,
						IO_scf_den_real_spin =  534,
						IO_scf_den_imag_spin =  535,
						IO_fock_real_a       =  536,
						IO_fock_imag_a       =  537,
						IO_fock_real_b       =  538,
						IO_fock_imag_b       =  539,
						IO_int_sym_a         =  563,  /* Integer symmetry assignments (alpha) */
						IO_int_sym_b         =  564,  /* Integer symmetry assignments (beta)  */
						IO_dypole_vel        =  572,
						IO_dens              =  603,  /* Density Matricies at diff. levels of theory */ 
						IO_excit_dens        =  633,  /* Excited-State CI densities */ 
						IO_cis_ampl          =  635,  /* CIS amplitudes */
						IO_S12_sec_deriv     =  656,  /* Non-symmetric S1 and S2 parts of Lagrangian 
                                               for MP2 or CIS second deribvatives */  
						IO_sub               =  991,
						IO_lnksub            =  992,
						IO_info              =  993,
						IO_phycon            =  994,
						IO_munit             =  995,
						IO_top               =  996,
						IO_mol               =  997,
						IO_ilsw              =  998,
						IO_overlay           =  999  };

};


#endif /* !GAUFILE_H */














