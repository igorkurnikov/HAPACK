/* common/constants.h
 * Warren Dukes
 * 2002-10-27
 */

#ifndef PNPCONSTANTS_H
#define PNPCONSTANTS_H

//#include <math.h>

/* used to convert dielectric constant to e/kT */
//#define EPKT	560.738438485
const double EPKT = 561.0;

#ifndef HARLEM_MOD

const double BOHR_TO_ANG       = 0.529177249;        //!< Bohrs to Angstroms
const double HARTREE_TO_KT     = 1059.642;       //!< HARTREE TO kT at 298K

#endif
//#define EPKT	560.46
/* used to convert mV to kT/e */
const double CONFAC = 25.69;
/** used to convert ions concentration [mol/L] to charge density[e/A^3]
 COANGS=Na[e]/1E27[A^3]*/
const double COANGS = 0.000602;
//#define COANGS	0.000602214199474747
/* used to convert M to k'2 (inverse of Debye Length) */
const double  DFACT = 3.047;

#endif
