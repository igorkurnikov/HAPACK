/*!  \file nuclacidmod.h

    Classes to model Nucleic Acids in HARLEM  

    \author Igor Kurnikov 
    \date  2002

*/
#if !defined(NUCLACIDMOD_H) 
#define NUCLACIDMOD_H

#include "hastl.h"
#include "haconst.h"
#include "halinalg.h"
#include "hastring.h"
#include "hacompmod.h"
#include "vec3d.h"

const int FLEX_FF    = 0;
const int AMBER91_FF = 1;
const int AMBER94_FF = 2;


class HaMolecule;

//! Computational module to model Nucleic Acids
class NuclAcidMod : public HaCompMod
{
public:
	 NuclAcidMod(MolSet* new_pmset=NULL);
	~NuclAcidMod();

	int SetStdParams();	

	int BuildNuclAcid();       //!< Build Nucleic Acid from the sequence
	int UpdateXYZ();           //!< Update Cartesian Coordinates for changed internal coordinates
	int MinEne();              //!< Optimize Energy of the system
	int CalcEne();             //!< Calculate Energy for the current conformation
	int CreateMolFromJumna();  //!< Create Molecule from JUMNA structures

	HaMolecule* FindDNAMol();  //!< Try to find DNA molecule in the current MolSet and assign it to the module 

	int SaveAXEfile(const char* axe_file); //!< save internal variables to an axefile
	int ReadAXEfile(const char* axe_file); //!< read internal variable from axefile

	std::string out_prefix; //!< prefix for output file
	std::string pdb_prefix; //!< prefix for pdb file

	HaVec_int nsym_unit; //!< Symmetry unit sizes for the DNA sections  

	int nbreak_1;  //!< break 1 of the DNA chain to specify different symmetry (see param isym_partN)
	int nbreak_2;  //!< break 2 of the DNA chain to specify different symmetry (see param isym_partN)

    HaVec_int homon_symm_flags; //!< Flags to indicate Homonomous(interstrand) symmetry for DNA sections (JUMNA homo,homo2,homo3)
	HaVec_int homon_symm_offs;  //!< residues offset to apply interstrand symmetry (mhomo,mhomo1,mhomo2).
	int ene_per_unit_flag;      //!< Flag to compute energy per unit cell rather then total energy (JUMNA ecen)

	HaVec_int dir_strand;       //!< Direction of the strand (positive means 5'->3' direction, negative 3'->5') 

	std::string seq[4];                 //!< Nucleotide sequence of the strands
	int SetSeq(const char* seq_str); //!< Set Nucleotide sequence of the DNA assuming the second strand is complimentary to the first
    int GenComplStrand();            //!< Set sequence of the second strand complimentary to the first strand

    int max_iter;                    //!< max number of minimization iterations

	int init_flag;      //!< flag to indicate that DNA(RNA) structure was created and data initialized

	double sup_helix_rad;        //!< radius of superhelix (0.0 for no superhelical constraints DNA)
	double sup_helix_pit;        //!< pitch of superhelix 
	double fst_twist;            //!< Twist of the first residue (when sup_helix_rad > 0.0)

	int SetFFtype(const char* new_ff_type); //!< Set Force Field Type
	int SetFFtypeIdx(const int i_ff_type);  //!< Set Force Field Type
    int force_field;             //!< force field type        

	int SetDielSlope( double slope_new);    //!< Set slope of sigmoidal dielectric constant
	double diel_slope;                      //!< Slope of the sigmoidal dielectric constant
	int SetPhosChrg( double phos_chrg_new); //!< Set Phosphate charge (FLEX Force-field)
	double phos_chrg;                       //!< Phosphate charge     (FLEX Force-field)

	HaVec_double hel_crd;       //!< Vector of Helical internal Coordinates
	HaVec_int    lock_hel;      //!< Lock flags of Helical coordinates

	HaVec_double bb_crd;        //!< Vector of Backbone Internal Coordinates
    HaVec_int    lock_bb;       //!< Lock flags of Backbone Internal Coordinates
    HaVec_int    bb_sidx;       //!< A vector of indexes of the starting positions of residue internal coordinates in int_coord 

	Vec3DValArray     u,h,o;        //!< helix axis parameters
	Vec3DValArray     uho,hho,ul;   //!< helix axis parameters
    HaVec_double  dif;          //!< inter nucleotide irregularity function
	
	HaMat_double hel;
	HaMat_double hold;
	HaMat_double vold;
	HaMat_double bend;
	HaMat_double vkin;

	HaMat_double tor;
	HaMat_double sug;
	HaMat_double bs_bs_pars;

	HaVec_int  ncrd_res;        //!< A number of internal coordinates per residue
	HaVec_int  ncrd_ind_res;    //!< A number of independent internal coordinates per residue
	
	int SetHelCoord(int i_strand, int i_res, int i_crd, double crd_val); //!< Set helicoidal internal coordinate  of a nucleotide
	int LockHelCoord(int i_strand, int i_res, int i_crd, int do_lock);   //!< Lock/unlock helicoidal internal coordinate of a nucleotide
    int IsHelCoordLocked(int ir, int i_crd);   //!< Check if the helical coordinate i_crd(0-based) is locked for the residue ir(0-based)

	int SetBBCoord(int i_strand, int i_res, int i_crd, double crd_val);  //!< Set Backbone internal coordinate  of a nucleotide
	int LockBBCoord(int i_strand, int i_res, int i_crd, int do_lock);    //!< Lock/Unlock Backbone internal coordinate of a nucleotide

	int update_var_flag;        //!< Flag to indicate that Var array needed to be updated due to lock changes or other

    int UpdateVarCoord();            //!< Update the array of non-locked internal coordinates
	int SetCoordsFromJumna();  //!< Set Coordinates from Jumna Structures
	int SetIntCoordsToJumna(); //!< Set Internal Coordinates To Jumna Structures 
	int SaveConfig();          //!< Save current configuration (internal and cartesian coordinates)
	int RestoreConfig();       //!< Restore configuration saved by SaveConfig()

	int BBOffset(int ir);        //!< offset of backbone cordinates of residue ir in sap
	int IsRibRes(int ir);        //!< Check is the residue a ribonucleotide
    int IsThymineRes(int ir);    //!< Check if the residue is a thymine
	int NumIndBBVarRes(int ir);  //!< The number of independent backbone variables of the residue
    int NumAllBBVarRes(int ir);  //!< The number of all backbone variables (indep and dep) of the residue
    int Is3endRes(int ir);       //!< Check is the residue at 3' DNA end
    int IsFstResInChain(int ir); //!< Check if the residue a first residue in chain
	int IdxResInChain(int ir);   //!< Index of the residue in the chain (ir - is 0-based , and return value 0-based)
	int ChainIdxOfRes(int ir);   //!< Index of the chain (1-based) for the residue with index ir( 0-based)
	int IsSupHlxConstr();        //!< Check if Supehelical symmetry is assigned to the DNA
	int NumBBCoord();            //!< Get a number of backbone internal coordinates
	int NumFreeBBCoord();        //!< Get a number of non-locked backbone internal coordinates
	int IdxLastHelCoord();       //!< Index (1-based) of the last helical variable in lock() and others arrays
	int IdxLastFreeHelCoord();   //!< Index (1-based) of the last non-locked helical variable in var() array
	int SetShlxTwist(double shlx_tw);             //!< Set twist of the first residue of a superhelix
	double GetAtomCrd(int at_num, int coord_num); //!< Get cartesian coordinates of the atom
	double GetVarNum(int crd_num);                //!< Get Unlocked variable by index(1-based) in var() array 

// Analysis function:

	int CalcPdistToShlxCnt();              //!< Print to the output distance of phosphorus atoms to the center of super-helix when superhelix constraints in place
	int CalcLocCrdSys();                   //!< Calculate local coordinate systems on bases 
	int CalcLocHlxCrd(int for_bp = TRUE);  //!< Calculate local helical coord of bases or base pairs
    int CalcAxis();                        //!< Calculate helical axis 
	int CalcBend();                        //!< Calc Bending
	int CalcGlobHlxCrd();                  //!< Calc Global helical coordinates
	int CalcBBCrd();                       //!< Calculate Backbone internal coordinated


	int CalcAxisPar1(HaMat_double& ax_hlxc, double& sum, HaVec_double& gra, HaVec_double& scp); //!< Calculate parameters characterizing helical axis
	
	int GetNRes() const;     //!< Get the number of residues in the DNA(RNA)
	
	double tot_energy;       //!< Current energy

	HaMolecule* p_dna_mol;   //!< the pointer to the DNA molecule

};

#if defined(HA_DISTRIB)
#undef INT_JUMNA
#endif

#if defined(INT_JUMNA)

extern "C"{
    extern void jumcall_(); //!< main subroutine
    extern void closejm_(int_4* ipl,int_4* iabas);
    extern void ecomp_(); //!< Calc energy FLEX force field
    extern void ecomp91_();
	extern void ecomp94_();
    extern void minim_();
	extern void disth_();
    extern void microb_();  //!< Conformation update during minimisation, 
	                        //!< var() - contain new helical vars and hel() - old ones corresp to current XYZ coords
    extern void penalty_(); //!< Calculation of restraint energy and gradients
	extern void helix_();   //!< List Helical, kink, ligand unternal variables and locking
	extern void backbo_();  //!< List backbone variables and updates C5' hydrogen positions
    extern void setvar_();  //!< Initial setup of variables arrays sap(), var, lar
    extern void equim_(int_4* newq,int_4* nspc,int_4* nhlc,int_4* nbcc,
		    		   int_4* nklc,int_4* nlgc,int_4* istop); //!< Check and list symmetry contractions
    extern void equivjm_(int_4* nspc,int_4* nbcc,int_4* nhlc,int_4* nklc, int_4* nlgc); //!< setup symmetry contractions of variables
    extern void putbac_(int_4* key); //!< Save(key=0)/restore(key=1) configuration
	extern void reset_(char* axe_file, int len); //!< Read axe file and reset internal coordinates
	extern void axeout_(char* axe_file, int len); //!< save axe file 
    extern void setbac_();                       //!< Set backbone conformation using int backbone conf read in reset(0 and equim()

     const int N0_J = 5;     //!< Maximum number of rings/subunits per ligand
     const int N1_J = 2000;  //!< Maximum number of atoms in the DNA oligomer
     const int N2_J = 74;    //!< Maximum number of nucleotides
     const int N3_J = 300;   //!< Maximum number of additional distance constraints
     const int N4_J = 230;   //!< Maximum number of atoms per nucleotide
     const int N5_J = 20;    //!< Maximum number of types of nucleotide
	 const int N6_J = 12000;  //!< Maximum number of physical backbone variables
	 const int N6A_J = 12000; //!< Maximum number of physical variables (lock only)
	 const int N7_J = 1100;  //!< Maximum number of independent, unlocked variables
	 const int N8_J = 200;   //!< Maximum number of physical variables per nucleotide
     const int N9_J = 20;    //!< Maximum number of non-bonded ligands
	   
//      typedef struct{
//         char finp[80];
//	  } jfinp_type;

//	  extern jfinp_type jfinp_;
 
	 typedef struct{
		  int iret; //!< flag return from the jumna because of the fatal error in a subroutine
	 } hacon_type;

	 extern hacon_type hacon_;

	  typedef struct{
		 char libn[80];
	  } jmdirs_type;

      extern jmdirs_type jmdirs_;


	  typedef struct{
		 char mac[32],lmo[32],lib[32],libm[32],out[32];
		 char axe[32],axl[32],noe[32],nol[32], test[32];
		 char pdb[32],ins[32],bar[32],parm[32];
//         char vc[32*N_CHA];
	  } cha_type;
 
	  extern cha_type cha_;
 
	  typedef struct{
		double acc,phos,epsi,epsr,plat,slope,rhbl,vfac,tfac,rfac,xfac,
               scale,rpiv,tpiv,fad,fan,damp,fnoew,fnoes,fnoea,df1,scnb,scee,
               catd,catr,catc,rad,pit,enit;
	    int_4  nshel,maxn,opt,mhomo,mhomo2,mhomo3,
               nick,limit,
			   nop,  //!< the number of open bases
			   nrib, //!< the number of riboses
			   ncat, //!< the number of cations ?
			   lig,
			   naxo;
	    int_2  sup,rcom,homo,homo2,homo3,diep,sum,link,ecen,
		       cyl,lcat,cent,autos,amber;
//		int_2  amber;
	  } datjm_type;

	  extern datjm_type datjm_;
      
	  typedef struct{
		  double xsum,ysum;
		  int_4  itwl[N2_J],
			     nvs[N7_J],
				 ihm[N2_J],
				 ksym[3],
				 nbrk[2],
			     isym;
		  int_4 isur; //!< position in var array of superhelix radius variable( 0 - no super helicity, -1 - fixed suphlx radius)
		  int_4 isup; //!< position in var array of superhelix pitch variable
		  int_2 ihl[N2_J];  //!< indicator of the break of symmetry section on the residue ?
	  } symjm_type;

	  extern symjm_type symjm_;

      typedef struct{
          double corm[N1_J*3]; //!< atom coordinates 
		  double dmon[N1_J];   //!< atomic charge
		  char mnam[N1_J*4];   //!< atom names
		  char munit[N1_J*4];  //!< residue names of the atoms first placed SUC (sugar atoms)
          int_4 imch[N1_J];    //!< element number of the atom
		  int_4	imty[N1_J];    //!< atom  class in the arrays of the force field such as angc,atclas etc..
		  int_4 icm[N1_J];
		  int_4 matm[3*N1_J];  //!< atom the internal coordinate belongs to
		  int_4 matd[N1_J*7];  //!< matd(i,7) - number of bonds of the atom i, matd(i,k) - atoms given atom bonded to 
		  int_4	nunit[N1_J];   //!< residue number
          int_4 nuc[N2_J+1];   //!< last atom number of the residue ( 1- based (fortran 0-based))
		  int_4 ncen[N0_J*N2_J+1];
		  int_4	kam;           //!< Number of atoms in the molecule
		  int_4	khm,lkm,kcen;  //!< 
	  } mrc_type; //!< DNA geometry, charges and linking DATA 

	  extern mrc_type mrc_;

      typedef struct{
          char seq[120]; //!< nucleotide sequence (all strands in sequence (blank is '-')
		  double hel[N2_J*6]; //!< nucleotide helical coords (X-disp,Y-disp,Z-disp,Rise,Inclination,Tip)
		  double vkink[N2_J*4];  //!< kink variables 4 per residue
		  double set[N2_J*N8_J]; //!< nucleotide backbone coords, set(nr,4) - glycosidic angle (O1'-C1'-N9-C8) or O1'-C1'-N1-C6)
		                         //!< set(nr,5),set(nr,6) - sugar torsions (O1'-C1'-C2'-C3') and (C1'-C2'-C3'-C4')
		                         //!< set(nr,7) - 1st backbone torsion epsilon (C4'-C3'-O3'-P)
		                         //!< set(nr,8) - 2nd backbone torsion dzeta   (C3'-O3'-P-O5')
		                         //!< set(nr,11) - Ribose hydroxyl torsion  ( set(nr,7) - for 3' nucleotide)
		                         //!< - all current values of independent variables
	      char code[N2_J*8]; //!< char[8] string per residue form input file to indicate input coordinates and locks
		  char kode[N2_J*8];
	      int_4 irec[N2_J];  //!< residue template number (in moljm) corresponding to the given residue
		  int_4	itr[N2_J];   //!< type of the topology? for STD res = 1 for 5' residues, =2 for 3' res and 2 for in strand, for nonstd > 6
		  int_4	ito[N2_J];   //!< for STD res = 1 for 5' residues, =3 for 3' res and 2 for in strand (4 - 5' RNA, 6 - 3' RNA , 2- in strand)
		  int_4	nst;         //!< the number of strands (see nstrands param)
		  int_4 nto;         //!< The index of the last nucleotide residue (number of nucleotides)
		  int_4	kseq;        //!< length of the sequence (length of the first strand)
		  int_4 ieq[50*4];   //!< ieq(l,i) - abs residue number of the nucleotide with an index (l)in the chain(i), = 0 - if position empty
		  int_4	ilq[N2_J*2]; //!< ilq(k,1) - index of the residue in the chain(relative to sequence of the 1st chain), ilq(k,2) - index of the chain of the residue
		  int_4	idr[4];      //!< The number of the residues in the strand and its direction (len_strand)
          int_2 kink[N2_J];  //!< flag to indicate the kink on the residue
		  int_2 lthy[N2_J];  //!< Thymine flag (if residue is T or S
	  } strjm_type;

	  extern strjm_type strjm_;

//	  typedef struct{
//	  common/strjm/seq,hel(n2,6),vkink(n2,4),set(n2,n8),code(n2),kode(n2),
//   1 irec(n2),itr(n2),ito(n2),nst,nto,kseq,ieq(50,4),ilq(n2,2),idr(4),
//     1 kink(n2),lthy(n2)
//	  } strjm_type;
//	  extern strjm_type strjm_;

	  typedef struct{
		double thy[N2_J];
		double rsr[N0_J*N2_J]; //!< values for a value of a ring closure contstraint
		int_4 iofs[N2_J];  //!< offset of nucleobase atoms in the nucleotide atom sequence? 
		int_4 iofe[N2_J];  //!< Number of atoms in the residue
		int_4 ithy[6];
		int_4 nith;
        int_4 neq[N7_J];   //!< Array of dimension of var to indicate symmetry group given variable belongs to (values should be equal) 
		                   //!< group is determined by abs(neq[]) all except one variable in the group has negative value of neq(i)
		                   //!< values of these varaible are compared to the value of one that has positive value of neq(i)(in equim.f)
		                   //!< neq set in ?
		int_4 isr[N0_J*2*N2_J]; //!< isr(lr,1,is), isr(lr,2,is) - atoms to define constraints for ring closure
		int_4 nsr[N2_J];        //!< Number of ring constraints in the residue
		int_2 ribose[N2_J];     //!< flag to indicate that nucleotide is ribonucleotide 
		int_2 cation[N2_J];     //!< flag to indicate that the nucleotide has a cation
	  } extjm_type;             //!< Number, Symmetry refernce vectors, etc

	  extern extjm_type extjm_;

	  typedef struct{
		   double rlig[N9_J*6]; //!< ligand variables ?
		   double slig[N9_J*2];
		   int_4 ilig[N9_J*6],lopt[N9_J],lpiv[N9_J];
           int_4 ntlg; //!< index of the last lig coord in sap,lock
		   int_4 nlgi; //!< index of the last lig coord in var
		   int_4 nlig; //!< the number of ligand residues
		   int_4 ntl;  //!< index of the last ligand residue (tot num res)( placed after strjm_.nto nucleotides)
		   char  lnam[N9_J*4];
		   int_2 locr[N2_J*N0_J];  //!<  locr(is,lr) - for res(ir)  lock sugar ring number lr(for nucleotide only one ring lr == 1) 
	 } lgd_type; //!< Non-bonded ligand variables

	 extern lgd_type lgd_;

     typedef struct{
		int_2 hst[N2_J*6];    //!<  hst(n2,6) flags to indicate that helical variable have been set in the input file
		int_2 bst[N2_J*N8_J]; //!<  bst(n2,n8) flags to indicate that backbone variables have been set in the input file
		int_2 vst[N2_J*4],
		lgi[N1_J],            
		lgj[N1_J],            //!< flag to compute interactions to the atom 
		lgu[N2_J];
	 } slg_type;

	 extern slg_type slg_;

     typedef struct{
	     int_4 iend[5]; //!< index of the last residue in the strand
		 int_4 ise[N2_J]; //!< for fst res = -5 ,last res 3 (positive strand dir), first res -3 and last res 5(negative strand dir)  
		 int_4 ienb[N2_J*2];
		 int_4 kapt[N2_J+1]; //!< kapt(is+1) - index of the last int coordinate in (sap) of the residue is
		 int_4 nsph; //!< index of the last position variable in  var array
		 } ind_type;

	 extern ind_type ind_;

	 typedef struct{
	    double sor[N4_J*N5_J*3]; //!< sor(i,k,j) - std coordinate j of atom i in the residue template k 
		                         //!< read from db file
		double smon[N4_J*N5_J];  //!< smon(i,k) - partial charge of the atom number i in residue template k
	    char snam[N4_J*N5_J*4];  //!< snam(i,k) - atom symbol for atom number i in residue template k
	    char suni[N4_J*N5_J*4];  //!< suni(i,k) - subresidue name of the atom ( SUCR, PHOS, ADNA etc)
	    int_4 nuni[N4_J*N5_J];   //!< nuni(i,k) - ? some flag 1 or 0 
	    char sub[N5_J*4];        //!< residue template name
	    int_4 isch[N4_J*N5_J];   //!< isch(i,k) - element of the atom i in the template k
		int_4 isty[N4_J*N5_J];   //!< isty(i,k) - force-field type ? of the atom i in the template k
		int_4 ics[N4_J*N5_J];    //!< ics(i,k)  - ? some flag 1 or 0 
        int_4 mats[3*N4_J*N5_J]; //!< mats(i,k) - internal coords? bonds? graph? for the residue template k 
		                         //!< read from db file
		                         //!< if mats(i,k) > 10000  => matm(i) set to 10000 - end atom in the tree?
		                         //!< if mats(i,k) == 0     => matm(i) set to O5' of the previous res in strand
		                         //!< 
		                         
		int_4 kas[N5_J]; //!< Number of atoms in the residue template
		int_4 khs[N5_J]; //!< dimension of mats for a given residue template k (internal coordinates)
		int_4 ksub;      //!< the number of residue templates in the library
	 } moljm_type;

	 extern moljm_type moljm_;

	 typedef struct{
	    double sap[N6_J]; //!< values of internal independent variables (not helical) go over residues kap(1,in) per residue, 
		                  //!< in = itr(is) - topological residue type, correspond to given coordinates
		double refg,refb,refh,refv;
		int_4 iap[N8_J*N4_J*N5_J];  //!< iap(l,j,in) - atoms that move for a int coord l of the residue, and j - atom of the residue
                                    //!< to calculate changes of the coordinates of atom j angle of rotation divided by iap(l,j,in) 
		int_4 nap[N8_J*7   *N5_J];  //!< nap(l,j,in) indexes of atoms of variables (j=1,2,3,4)
		                            //!< (val ang if(nap(l,4,in)== 0)  (ia,ib.ic,id) otherwise torsion 
		                            //!< nap(l,5,in) - indicate ic(=1,-1), ib(=-2,2) - belongs to residues in the next(=1,2) or prev(=-1,-2) res
		                            //!< nap(l,6,in) - if < 0 indicate the end of dependent vars?(see kaps)                            
                                    //!< type of force-field valence angle (set in kaps)

		int_4 kap[3*N5_J]; //!< from (1 to kap(1,in)) - independent backbone variables, torsions in strjm_.set
		                   //!< from (kap(1,in)+1 to kap(2,in)) - dependent backbone variables in strjm_.set, lock etc 
		                   //!< from (kap(1,in)+1 to kap(1,in)+nsr(is)*5 - dependent ring variables
		                   //!< Number of independent variables in the residue of the given topology type in= (itr(is))
	 } flx_type;

	 extern flx_type flx_;

	 typedef struct{
	    double force[N1_J*3]; //!< Array of forces
		double tor[N1_J*3];   //!< 
		double fot[N6_J];
		double ener; //!< Total energy with penalty term 
		double elec; //!< electrical energy  
        double repl; //!< VdW repulsion energy 
		double disp; //!< VdW dispersion energy (attraction)
		double eang; //!< Valence angles energy
		double etog; //!< Torsion energy
		double epen; //!< Penalty energy term
	    int_2  i23[8*N1_J]; //!< excluded atom pair list ?
		int_2  i34[8*N1_J];
		int_2  elim[N1_J];  //!< add excluded atom pair list ?
	 } enf_type;

	 extern enf_type enf_;

	 typedef struct{ 
	    double var[N7_J]; //!< Non-locked independent variable (before symmetry contraction??)
		                  //!< first go position variable from strjm_.set 1 - (kap(1,in) for every residue is in = itr(is)
		                  //!< then go 6  helical variables for each residue
		                   //!< then may be 4 kink variable per residue with kink
		                   //!< then ligand interaction variables
		                   //!< then thymine methyl variable
		double gra[N7_J];  //!< energy gradient on independent variables
		double con[N0_J*N2_J + N3_J];
		double scl[N7_J];   //!< scale factors of indep variables for conjugate gradient minimization routine(not used?)
		int_4 ncon;
		int_4 nvrc;
		int_4 ntba; //!< idx of last variable in sap,lock etc array before helical variables - total number of backbone vars
		int_4 nbac; //!< idx of last variable in var array before helical variables - the number of backbone unlocked ind vars
		int_4 nthe; //!< idx of last helical variable in sap,lock etc (before kink vars)
		int_4 nhel; //!< idx of last helical variable in var (before kink vars)  
		int_4 ntki; //!< idx of last kink var in sap,lock etc
		int_4 nkin; //!< idx of last kink var in var   (6 above set in setvar)
		int_4 ntri;
		int_4 ntot;
		int_4 nvar; //!< the number of independent variables - dimension of var ?
		int_4 nrin;
	    int_2 lar[N7_J];  //!< indx in var = .false. (x,y disp,rise) .true.(inc, tip, twist) ,for backbone vars .true.(for val angles)
	    int_2 lock[N6A_J]; //!< (dim of sap) lock on the coordinate, if 1 not changed in minimization
	 } mnn_type;

	 extern mnn_type mnn_;
};

	 typedef struct{ 
        double rnoe[N3_J]; //!< Low(?)   bound for restraints for jnoe(k) = 3,28,4 or 16
		double bnoe[N3_J]; //!< Upper(?) bound for restraints for jnoe(k)  = 3,28,4 or 16
		double fnoe[N3_J];
		double snoe;
		char   knam[200*4*2];
		int_4 ktyp[200];
		int_4 icon[N2_J*2];
        int_4 inoe[N3_J*9]; //!< inoe(k,1),inoe(k,2),inoe(k,3),inoe(k,4) - atom numbers that define torsional constrain
                            //!< inoe(k,5),inoe(k,6),inoe(k,7),inoe(k,8) - at num of sec torsion and inoe(k,9) - sign of contr oper for 2 tors
                            //!< (inoe(k,1),inoe(k,2)) - atom number for dist constr
		                    //!< for jnoe(k)= 21 or 22 (sugar pucker restr) inoe(k,1) - res number of the nucleotides, inoe(k,2) - pucker code
		int_4 ih68[N2_J];   
		int_4 jnoe[N3_J];   //!< Type of the constraint (size nnoe)
		                    //!< jnoe(k) = 1(D) or 3(d)     - interatomic distances fixed or range
		                    //!< jnoe(k) = 26(D2) or 28(d2) - interat dist sum or diff, fixed or range
		                    //!< jnoe(k) = 2(T) or 4(t) - torsion constr fixed or range
		                    //!< jnoe(k) = 14(T2) or 16(t2) - sum or diff of two torsions fixed or (range)
		                    //!< jnoe(k) = 21(S)  or 22(s)  - sugar pucker restraints 
		int_4 nnoe;         //!< the number of constrained variables
		int_4 ndcs;
		int_4 ndiv[100];
		int_4 nlin;
	 }nmrjm_type;

	 extern nmrjm_type nmrjm_;

//     Amber force field parameters:
//     common/para/
//     ro(57) - VdW radii for computing vdw foreced
//     reps(57) - energy minimum for identical atom VdW interactions
//     angc(191),fangc(191),btor(88),
//     $ gtor(88),ftor(88),bdtor(30),gdtor(30),fdtor(30),
//     amb(57,57), - A and B coef in (A/r12 + B/r6) : amb(i,j) = sqrt(reps(i)*reps(j))*(ro(i)+(ro(j))**12
//     bmb(57,57)                                     bmb(i,j) = 2*sqrt(reps(i)*reps(j))*(ro(i)+(ro(j))**6
//     amh(9,11),bmh(9,11),mtor(88),itor(88,4),iadt(30,4),
//     $ iangc(191,3),mof,nof,nang,ntor,ntors,nad,
//     atclas(57) - character*2 - of atom class
//

#endif

#endif // end !defined(NUCLACIDMOD_H) 

