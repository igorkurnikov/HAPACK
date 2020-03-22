#ifndef G94_PROTOS_H
#define G94_PROTOS_H
//! \file  g94_protos.h
//!  
//!  Prototypes of functions of GAUSSIAN
//!
//!


extern "C" {

extern void fopen_ (int* iunit, int* imode, char* fname_for,
                    int* extend, int* ialloci, int junk2 );

extern void fclose_ (int* fort_iunit, int* idisp );

extern void fdump_(void);

extern void fileio_(int* ioper, int* ifilno, int* len, 
                    void* q, int* ipos);

extern void getnb6_(int* nbas6d);     

extern void setncw_( int* byte_in_f);

extern int itqry_( int* i_file);

extern void ilsw_( int* ioper, int* i_file, int* i_what);

extern void cpyfil_ (int* ifile1, int* ifile2, char* buf, int* len);

extern int filnum_( int* ifile, int* iunit);

extern void square_( double* a,  double* b, int* max, 
                     int* n,int* k);

extern void renorm_(logical* ToRawP); // Renormalize coef in B
extern void bprint_(int* kop, int* nprims); // print basis set info in /B/

extern void psudag_(logical* tr,logical* fl,int* nbzdo, int* naezdo,
					int* navzdo,double* zero, double* ff,
                    double* cmo_a, double* ene_mo_a,
					double* scr1, double* anorm);



extern void zdofnm_(int* natoms, double* core_ch, double* coord,
				   double* force);

extern void formp_ha_(logical* tr,int* nbzdo, int* nbzdo2,
				   int* naezdo, double* cmo_a, double* pa);

extern void ovlp_(int* iout,int* iprint,double* s,int* isdim, int* itrans,
		  int* ipurd1, int* ipurf1, int* ipurd2, int* ipurf2,
		  double* v, int* mdv);

extern void vibfrq_(int* iout, int* iprint, int* natoms, 
					int* multip, int* ian, double* c, int* nat3,
					double* ffx, double* ddip,
                    logical* dopol, double* dpolar, logical* dovcd, double* aat,
					double* atmass, double* orthog, double* vv, 
					double* vecout, double* trialv, double* e2, 
					double* table, double* phycon, logical* ifatm, logical* ifwrt, 
					double* symms, int* nvib, int* nimag, 
					int* ipg, int* maxnz, int* nz, int* iz, int* ianz, 
					int* iproj, double* fx, double* scrmat, 
					double* tuser, double* puser, double* fcscal, 
					double* gen, double* htherm, double* gtherm, 
					double* trot, double* tanirc,
                    logical* dored, logical* rdcrd,
					int* izred, double* value, int* ntred, 
					int* ntbond,int* ntang, int* ntdih, int* ntrrot,
					double* v, int* mdv);


/* see the detailed description in denbas.F */

#if GAUSSVER == 98
extern void denbas_(int* iout,int* icalc,int* idgst,
	int* mindrv,int* maxdrv,int* minmlt,int* maxmlt,
    int* icnbeg,int* icnend, int* idena,double* dena, double* denb,
    int* iscf, int* nbasis, int* nmatd, int* nbt,
	double* c, logical* useatt, int* iattyp, double* wtgrid, int* natoms,  /* Atomic coordinates, weights of charges on the grid  */ 
	double* cgrid, int* ngrid,              /* Coordinates of grid points and number of points in the grid */
	double* val, double* valmat, double* accinp,
	int* isymcn, int* nop, int* neqshl, double* rotop, /*symmetry stuff */
	int* iprint, int* ipflag, logical* allowp, logical* dospar,
	int* i_meg, double* v_meg, int* lenmeg);
#endif

/* see the detailed description in  fofdir.F */

#if GAUSSVER == 98
extern void fofdir_(int* iout, int* iprint, int* ihmeth,
					int* iopcl, int* icntrl,
	                int* iraf, int* ipflag, logical* allowp,  /* PRISM parameters */
					int* icnbeg, int* icnend, 
					logical* addh, logical* initf, logical* dopurf, 
					double* accdes, double* scahfx, 
					int* nmat, int* nmats, int* nmatt, 
					int* nbasis, 
		int* isym2e, int* nsymop, int* nopab, 
		double* neqatm, double* neqshl, double* rotop, double* neqbas,
	    double* ha, double* hb, 
		double* pza, double* pzb, double* pa, double* pb, 
		double* fa, double* fb, double* fda, double* fdb, 
		int* natoms, int* ian, double* c, logical* useatt, int* iattyp,
		logical* frozen, double* fxyz, double* ffxyz,
		int* lensav, int* nsaved, 
		double* r03, double* r1, double* r2, 
		int* mdv, double* core);

#endif

extern void fofmem_(int* iprint, int* iopt,
					logical* initf, logical* dopurf, int* isym2e,
					int* nbasis,int* nmat, int* nmats,
                    int* nsymop,int* neqbas, int* ij,
					double* da, double* db, 
					double* r1, double* r2, double* r3,
					double* fa, double* fb,
					double* v, int* mdv);



extern int intowp_(int* n_ints);

extern void rootmt_(double* a, double* b, double* aa, double* bb,
					int* mdim, int* nbas, int* inv);

extern void regraf_(int* itype, int* nbasis,
		double* r0, double* r1, double* r2, double* r3,
		double* v, int* mdv);

}



#endif /* !G94_PROTOS_H */
    









