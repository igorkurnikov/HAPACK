/*!  electrostmod.h

     Classes to perform Continuum Electrostatics Calculations in HARLEM 

    \author Igor Kurnikov , University of Pittsburgh 
    \date 1999-2002
*/

#ifndef ELECTROSTMOD_H
#define ELECTROSTMOD_H

#include "hastl.h"
#include "vec3d.h"
#include "hastring.h"
#include "hacompmod.h"
#include "hasurface.h"

class HaSurface;
class HaMolSet;
class AtomGroup;
class HaResidue;
class AltChemState;

const int RDX_SHFT_VS_VAC=0x0001;            //!<  Compute redox-potential shift vs redox-group in vaccum 
const int RDX_SHFT_VS_SOLV=0x0002;           //!<  Compute redox-potential shift vs redox-group in solvent

const int FILL_CHARGE_PROD_REV_DIST = 0;     //!< choices for fill_charge_mode: Atom charge set to the 8 nearest grid points
                                             //!< proportionally to the product of reverse distances to the atom position (as in DELPHI) 
const int FILL_CHARGE_SUM_REV_DIST = 1;      //!< Atom charge set to the 8 nearest grid points
                                             //!< proportionally to the product of reverse distances to the atom position (as in DELPHI) 

enum ELPOT_FILE_FORMAT { HLM_F3D_BIN = 0, HLM_F3D_FORM = 1, ELPOT_DELPHI_BIN = 2}; 


//! Class to perform continuum electrostatic calculations on the grid
class ElectrostMod : public HaCompMod
{
public:
	ElectrostMod(HaMolSet* new_phost_mset = NULL,const int new_mtype=COMP_MOD_ELECTROST);
	virtual ~ElectrostMod();

	virtual void clear();    //!< Clear data and set default parameters of the module 

	bool SetStdParam();      //!< Set default parameters of the module
	bool SaveParamFile();    //!< Save input parameter file for DELPHI
	bool SaveChargeFile();   //!< Save file of atomic charges for DELPHI
	bool SaveRadiusFile();   //!< Save atomic radii file for DELPHI
	bool SaveCoordFile();    //!< Save atomic coordinates file for DELPHI
	bool RunDelphi(RunMode rmode); //!< Run Delphi calculations with saved parameters
    virtual bool run(RunMode rmode=RUN_FOREGROUND);  //!< Run continuum electrostatics caclulations with current parameters

	double CalcETReorgEne();  //!< Calculate outer-sphere Electron Transfer reorganization energy ( -1.0 if the error ) 
	bool CalcAltStatePK(AltChemState* alt_res_st, AtomGroup* active_atoms = NULL); //!< Calculate pK value of tht alternative state, if active atoms not equal NULL only they assumed to be present
	bool CalcRedoxPotShft();  //!< Calculate the redox-potential shift of the donor due to solvation
    double CalcAvgPotOn(PointContainer* ptlist); //!< Calculate Average potential on points

	bool ReadTotEne(double &tot_ene, const char* fname= "DELPHI_RES.OUT"); //!< Read Total Energy From DELPHI_RES.OUT file

	int LoadElPotFromFile(int format = HLM_F3D_BIN); //!< Load Electrostatic potential from file 
	
	bool BuildPotIsoSurface(); 
	bool CalcIndCharge();  //!< Calculate charges induced on the dielectric boundary between the molecule and the solvent
	bool PlotIndCharge();  //!< Plot the distribution of the solvent induced charges
	bool FillChargeMap();  //!< distribute atom point charges on 3D grid charge_map
	bool BuildPotVdwDots();   
	bool ColorMolSurfElPot(); //!< Color the molecular surface by the electrostatic potential
    bool ColorDotStruct(DotStruct* dotstruct);//!< Color dot struct according to el_pot  < mikola 27 july 2006
 
// Set boundaries of the molecule 

	bool AddBoundaryAtoms(); //!< Add Ficticious Atoms to the coordinate file to mark specified boundaries
	
	Vec3D min_coord; //!< Point marking minimal x,y,z coordinates for the calculations
	Vec3D max_coord; //!< Point marking maximal x,y,z coordinates for the calculations

	bool SetBoundaryAtoms(const double xmin, const double ymin, const double zmin,  
						  const double xmax, const double ymax, const double zmax); //!< Set boundaries of the molecule in the calculations

	bool ClearBoundaryAtoms(); //!< Clear boundary atoms

// Delphi run parameters

	std::string param_file_title;  //!< Title of the DELPHI job
	int nx, ny, nz;            //!< dimensions of the grid, default = 65
	double perfil;             //!< percent of the grid occupied by the molecule, default = 60.0
	double offset[3];          //!< dimensions of the box, or scale,  default =0.0
	double epsi, epsout;       //!< internal and external dielectric constant, default espi = 4.0, epsout = 81.0
	double rionst;             //!< ionic strength , default = 0.0
	double exrad, radprb;      //!< Ion exclusion radius, and probe radius, default exrad = 2.0, radprb = 1.4 
	int boundary;              //!< Boundary: 0=Zero,1=Apprx Coul,2=Focus,3=Coul , default = 0
    int iper[3];               //!< periodic boundary flags
	int nlit;                  //!< number of linear Poisson-Boltz. iterations, default = -1
	int nnit;                  //!< number of non-linear Poisson-Boltz. iterations, default = 0
	int iconc, ibios;          //!< potential output control 
	int isite;                 //!< flag to print electric field at spec. points
	int iatout;                //!< flag to to write modif. pdb
	std::string toplbl;           //!< title of the job
	int isph;                  //!< flag for spherical charge distrib
	int ipdbwrt;               //!< flag for writing pf pdb file
    int ifrcwrt;               //!< flag for writing frc file
	std::string enc;              //!< energy print control flags 
	int igraph, ipotent, icon1, icon2;
	int imem;                  //!< flag to model  membrane 
	int phiwrt;                //!< flag to write potential map
	int ihs,isen,ish;          //!< flags to write expanded surface data file, 
                               //!< surface energy file and surface charge file
	double elpot_low_val;      //!< low and high value of the electrostatic potential 
	double elpot_high_val;     //!< for coloring of the molecular surface (in kT)
	
	std::string param_file_name;   //!< Delphi input parameters file name
	std::string charge_file_name;  //!< Delphi atomic charges file name
	std::string radius_file_name;  //!< Delphi atomic radii file name
	std::string coord_file_name;   //!< Delphi atomic coordinates file name
	std::string log_file_name;     //!< Delphi log file name 

	std::string elfield_fname;     //!< name of the file with electrostatic field

	float pot_isolevel;  //!< value of potential for isosurface plot
	int dots_number;     //!< number of dots for dot_surface plot

    double tot_ene;  //!< Current total electrostatic energy of the system (in KT) 

	HaField3D el_pot_map; //!< Electrostatic poteintial Map on 3D Grid 
	list<HaSurface*> Surfaces;

	int fill_charge_mode; //!< Mode to spread atomic charges on the grid (FILL_CHARGE_NEAREST_PT = 1, FILL_CHARGE_CUBE=2)
    
	HaField3D charge_map; //!< Charges of atoms on 3D Grid 
	HaField3D ind_charge_map; //!< Induced charges on 3D Grid 
	HaField3D ConcMap0;//!< Maps of Induced Concentrations for +
	HaField3D ConcMap1;//!< Maps of Induced Concentrations for -
	
	int rdx_shft_mode;   //!< reference state to compute redox-potential shift of the donor: RDX_SHFT_VS_VAC, RDX_SHFT_VS_SOLV
    double rdx_shift;    //!< Value of the calculated redox-potential shift
	double ddG;
	double delta_pK;

	HaVec_double axx_ene; //!< Axxiliary arrays of energies to return intermediate energies from some functions

	static int ActiveElectrMod; //!< Choice of Electrostatic Module Engine (DELPHI-based(old) (COMP_MOD_ELECTROST) or Kolya's (COMP_MOD_EL)
protected:

};


extern "C"
{
	void openpm_(int* pnx_f, int* pny_f, int* pnz_f,
		         freal* pxmin_f, freal* pxmax_f, freal* pymin_f,
				 freal* ymax_f,  freal* zmin_f,  freal* zmax_f,
		         char* fname, int fnmlen); 

	void loadpm_(int* pnx_f, int* pny_f, int* pnz_f, 
		         freal* fmap);

#if defined(INT_DELPHI)

const int NGCRG_Q = 200000; //!< maximum number of grid points that may be assigned charge 

#if defined(ELECTROSTMOD_CPP)


void rdprm2_(logical* iautocon,
		         freal* epsin, freal* epsout,
				 freal* rionst, freal* exrad, freal* radprb,
				 int_4* ibctyp, logical* iper, int_4* nlit, int_4* nnit,
				 logical* iconc, logical* imem, 
				 int_4* icon1, int_4* icon2);

void setatq_(freal* xn2, freal* rad3, freal* chrgv4, int_4* natom, freal* ptr_gr_cent, freal* pscale);

void getpotm_(int_4* ngrid, freal* potmap, 
			 freal* pxmin,  freal* pymin, freal* pzmin, 
             freal* pxmax,  freal* pymax, freal* pzmax);

void geteneq_(freal* ptot_ene);
    
#endif

void delphi_(freal* x, int_4* mgrid, freal* phi, freal* phi1, freal* phi2, freal* phi3, freal* achrg,
			 int_4* ieps, int_4* ieps_2, int_4* ideb, 
			 freal*  sf1, freal*   sf2, freal* qmap1, freal*  qmap2, freal* debmap1, freal* debmap2,
			 freal* bndx1,freal* bndx2,freal* bndx3,freal* bndx4,freal* bndx,freal* bndy,freal* bndz,
			 freal* db, int_4* ibgrd,int_4* idpos,int_4* ioff,int_4* iepsv); 

typedef struct{
	freal scale;
	freal oldmid[3];
	int_4 igrid;
	int_4 ibc;
	freal gten;
	freal cgbp[2*NGCRG_Q];
	freal gval[NGCRG_Q];
	freal rmmin;
	freal rmmax;
} scaleq_type; //!< Non-bonded ligand variables

extern scaleq_type scaleq_;

typedef struct{
	int_4 iphi,iphi1,iphi2,iphi3;
} hadel1_type;

extern hadel1_type hadel1_; 
	

#endif

}


#endif // end if !defined(ELECTROSTMOD_H) 
