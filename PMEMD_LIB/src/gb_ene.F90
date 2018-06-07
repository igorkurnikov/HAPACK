#include "copyright.i"

!*******************************************************************************
!
! Module:  gb_ene_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module gb_ene_mod

  use gbl_datatypes_mod

  implicit none

! Global data definitions:

  double precision, allocatable, save           :: gbl_rbmax(:)
  double precision, allocatable, save           :: gbl_rbmin(:)
  double precision, allocatable, save           :: gbl_rbave(:)
  double precision, allocatable, save           :: gbl_rbfluct(:)

! Private data definitions:

  double precision, allocatable, save, private  :: reff(:)
  double precision, allocatable, save, private  :: psi(:)
  double precision, allocatable, save, private  :: rjx(:)
  double precision, allocatable, save, private  :: r2x(:)
  double precision, allocatable, save, private  :: sumdeijda(:)
  double precision, allocatable, save, private  :: vectmp1(:)
  double precision, allocatable, save, private  :: vectmp2(:)
  double precision, allocatable, save, private  :: vectmp3(:)
  double precision, allocatable, save, private  :: vectmp4(:)
  double precision, allocatable, save, private  :: vectmp5(:)

  logical,          allocatable, save, private  :: skipv(:)

  integer,          allocatable, save, private  :: jj(:)
  integer,          allocatable, save, private  :: neck_idx(:)

  ! gb_neckcut: 2.8d0 (diameter of water) is "correct" value but
  ! larger values give smaller discontinuities at the cut:

  double precision, parameter    :: gb_neckcut = 6.8d0

  ! Lookup tables for position (atom separation, r) and value of the maximum
  ! of the neck function for given atomic radii ri and rj. Values of neck
  ! maximum are already divided by 4*Pi to save time. Values are given
  ! for each 0.05 angstrom between 1.0 and 2.0 (inclusive), so map to index
  ! with dnint((r-1.0)*20)).  Values were numerically determined in
  ! Mathematica; note FORTRAN column-major array storage, so the data below
  ! may be transposed from how you might expect it.

   double precision, parameter  :: neckMaxPos(0:20,0:20) = reshape((/ &
     2.26685,2.32548,2.38397,2.44235,2.50057,2.55867,2.61663,2.67444, &
     2.73212,2.78965,2.84705,2.9043,2.96141,3.0184,3.07524,3.13196, &
     3.18854,3.24498,3.30132,3.35752,3.4136, &
     2.31191,2.37017,2.4283,2.48632,2.5442,2.60197,2.65961,2.71711, &
     2.77449,2.83175,2.88887,2.94586,3.00273,3.05948,3.1161,3.1726, &
     3.22897,3.28522,3.34136,3.39738,3.45072, &
     2.35759,2.41549,2.47329,2.53097,2.58854,2.646,2.70333,2.76056, &
     2.81766,2.87465,2.93152,2.98827,3.0449,3.10142,3.15782,3.21411, &
     3.27028,3.32634,3.3823,3.43813,3.49387, &
     2.4038,2.46138,2.51885,2.57623,2.63351,2.69067,2.74773,2.80469, &
     2.86152,2.91826,2.97489,3.0314,3.08781,3.1441,3.20031,3.25638, &
     3.31237,3.36825,3.42402,3.4797,3.53527, &
     2.45045,2.50773,2.56492,2.62201,2.679,2.7359,2.7927,2.8494,2.90599, &
     2.9625,3.0189,3.07518,3.13138,3.18748,3.24347,3.29937,3.35515, &
     3.41085,3.46646,3.52196,3.57738, &
     2.4975,2.5545,2.61143,2.66825,2.72499,2.78163,2.83818,2.89464, &
     2.95101,3.00729,3.06346,3.11954,3.17554,3.23143,3.28723,3.34294, &
     3.39856,3.45409,3.50952,3.56488,3.62014, &
     2.54489,2.60164,2.6583,2.71488,2.77134,2.8278,2.88412,2.94034, &
     2.9965,3.05256,3.10853,3.16442,3.22021,3.27592,3.33154,3.38707, &
     3.44253,3.49789,3.55316,3.60836,3.66348, &
     2.59259,2.6491,2.70553,2.76188,2.81815,2.87434,2.93044,2.98646, &
     3.04241,3.09827,3.15404,3.20974,3.26536,3.32089,3.37633,3.4317, &
     3.48699,3.54219,3.59731,3.65237,3.70734, &
     2.64054,2.69684,2.75305,2.80918,2.86523,2.92122,2.97712,3.03295, &
     3.0887,3.14437,3.19996,3.25548,3.31091,3.36627,3.42156,3.47677, &
     3.5319,3.58695,3.64193,3.69684,3.75167, &
     2.68873,2.74482,2.80083,2.85676,2.91262,2.96841,3.02412,3.07976, &
     3.13533,3.19082,3.24623,3.30157,3.35685,3.41205,3.46718,3.52223, &
     3.57721,3.63213,3.68696,3.74174,3.79644, &
     2.73713,2.79302,2.84884,2.90459,2.96027,3.01587,3.0714,3.12686, &
     3.18225,3.23757,3.29282,3.34801,3.40313,3.45815,3.51315,3.56805, &
     3.6229,3.67767,3.73237,3.78701,3.84159, &
     2.78572,2.84143,2.89707,2.95264,3.00813,3.06356,3.11892,3.17422, &
     3.22946,3.28462,3.33971,3.39474,3.44971,3.5046,3.55944,3.61421, &
     3.66891,3.72356,3.77814,3.83264,3.8871, &
     2.83446,2.89,2.94547,3.00088,3.05621,3.11147,3.16669,3.22183, &
     3.27689,3.33191,3.38685,3.44174,3.49656,3.55132,3.60602,3.66066, &
     3.71523,3.76975,3.82421,3.8786,3.93293, &
     2.88335,2.93873,2.99404,3.04929,3.10447,3.15959,3.21464,3.26963, &
     3.32456,3.37943,3.43424,3.48898,3.54366,3.5983,3.65287,3.70737, &
     3.76183,3.81622,3.87056,3.92484,3.97905, &
     2.93234,2.9876,3.04277,3.09786,3.15291,3.20787,3.26278,3.31764, &
     3.37242,3.42716,3.48184,3.53662,3.591,3.64551,3.69995,3.75435, &
     3.80867,3.86295,3.91718,3.97134,4.02545, &
     2.98151,3.0366,3.09163,3.14659,3.20149,3.25632,3.3111,3.36581, &
     3.42047,3.47507,3.52963,3.58411,3.63855,3.69293,3.74725,3.80153, &
     3.85575,3.90991,3.96403,4.01809,4.07211, &
     3.03074,3.08571,3.14061,3.19543,3.25021,3.30491,3.35956,3.41415, &
     3.46869,3.52317,3.57759,3.63196,3.68628,3.74054,3.79476,3.84893, &
     3.90303,3.95709,4.01111,4.06506,4.11897, &
     3.08008,3.13492,3.1897,3.2444,3.29905,3.35363,3.40815,3.46263, &
     3.51704,3.57141,3.62572,3.67998,3.73418,3.78834,3.84244,3.8965, &
     3.95051,4.00447,4.05837,4.11224,4.16605, &
     3.12949,3.18422,3.23888,3.29347,3.348,3.40247,3.45688,3.51124, &
     3.56554,3.6198,3.674,3.72815,3.78225,3.83629,3.8903,3.94425, &
     3.99816,4.05203,4.10583,4.15961,4.21333, &
     3.17899,3.23361,3.28815,3.34264,3.39706,3.45142,3.50571,3.55997, &
     3.61416,3.66831,3.72241,3.77645,3.83046,3.8844,3.93831,3.99216, &
     4.04598,4.09974,4.15347,4.20715,4.26078, &
     3.22855,3.28307,3.33751,3.39188,3.4462,3.50046,3.55466,3.6088, &
     3.6629,3.71694,3.77095,3.82489,3.8788,3.93265,3.98646,4.04022, &
     4.09395,4.14762,4.20126,4.25485,4.3084 &
     /), (/21,21/))

   double precision, parameter  :: neckMaxVal(0:20,0:20) = reshape((/ &
     0.0381511,0.0338587,0.0301776,0.027003,0.0242506,0.0218529, &
     0.0197547,0.0179109,0.0162844,0.0148442,0.0135647,0.0124243, &
     0.0114047,0.0104906,0.00966876,0.008928,0.0082587,0.00765255, &
     0.00710237,0.00660196,0.00614589, &
     0.0396198,0.0351837,0.0313767,0.0280911,0.0252409,0.0227563, &
     0.0205808,0.0186681,0.0169799,0.0154843,0.014155,0.0129696, &
     0.0119094,0.0109584,0.0101031,0.00933189,0.0086348,0.00800326, &
     0.00742986,0.00690814,0.00643255, &
     0.041048,0.0364738,0.0325456,0.0291532,0.0262084,0.0236399, &
     0.0213897,0.0194102,0.0176622,0.0161129,0.0147351,0.0135059, &
     0.0124061,0.0114192,0.0105312,0.00973027,0.00900602,0.00834965, &
     0.0077535,0.00721091,0.00671609, &
     0.0424365,0.0377295,0.0336846,0.0301893,0.0271533,0.0245038, &
     0.0221813,0.0201371,0.018331,0.0167295,0.0153047,0.014033, &
     0.0128946,0.0118727,0.0109529,0.0101229,0.00937212,0.00869147, &
     0.00807306,0.00751003,0.00699641, &
     0.0437861,0.0389516,0.0347944,0.0311998,0.0280758,0.0253479, &
     0.0229555,0.0208487,0.0189864,0.0173343,0.0158637,0.0145507, &
     0.0133748,0.0123188,0.0113679,0.0105096,0.0097329,0.00902853, &
     0.00838835,0.00780533,0.0072733, &
     0.0450979,0.0401406,0.0358753,0.0321851,0.0289761,0.0261726, &
     0.0237125,0.0215451,0.0196282,0.017927,0.0164121,0.0150588, &
     0.0138465,0.0127573,0.0117761,0.0108902,0.0100882,0.00936068, &
     0.00869923,0.00809665,0.00754661, &
     0.0463729,0.0412976,0.0369281,0.0331456,0.0298547,0.026978, &
     0.0244525,0.0222264,0.0202567,0.0185078,0.0169498,0.0155575, &
     0.0143096,0.0131881,0.0121775,0.0112646,0.010438,0.00968781, &
     0.00900559,0.00838388,0.00781622, &
     0.0476123,0.0424233,0.0379534,0.034082,0.0307118,0.0277645, &
     0.0251757,0.0228927,0.0208718,0.0190767,0.0174768,0.0160466, &
     0.0147642,0.0136112,0.0125719,0.0116328,0.0107821,0.0100099, &
     0.00930735,0.00866695,0.00808206, &
     0.0488171,0.0435186,0.038952,0.0349947,0.0315481,0.0285324, &
     0.0258824,0.0235443,0.0214738,0.0196339,0.0179934,0.0165262, &
     0.0152103,0.0140267,0.0129595,0.0119947,0.0111206,0.0103268, &
     0.00960445,0.00894579,0.00834405, &
     0.0499883,0.0445845,0.0399246,0.0358844,0.032364,0.0292822, &
     0.0265729,0.0241815,0.0220629,0.0201794,0.0184994,0.0169964, &
     0.0156479,0.0144345,0.0133401,0.0123504,0.0114534,0.0106386, &
     0.00989687,0.00922037,0.00860216, &
     0.0511272,0.0456219,0.040872,0.0367518,0.0331599,0.0300142, &
     0.0272475,0.0248045,0.0226392,0.0207135,0.0189952,0.0174574, &
     0.0160771,0.0148348,0.0137138,0.0126998,0.0117805,0.0109452, &
     0.0101846,0.00949067,0.00885636, &
     0.0522348,0.0466315,0.0417948,0.0375973,0.0339365,0.030729, &
     0.0279067,0.0254136,0.023203,0.0212363,0.0194809,0.0179092, &
     0.016498,0.0152275,0.0140807,0.013043,0.012102,0.0112466, &
     0.0104676,0.00975668,0.00910664, &
     0.0533123,0.0476145,0.042694,0.0384218,0.0346942,0.0314268, &
     0.0285507,0.026009,0.0237547,0.0217482,0.0199566,0.018352, &
     0.0169108,0.0156128,0.0144408,0.0133801,0.0124179,0.011543, &
     0.010746,0.0100184,0.00935302, &
     0.0543606,0.0485716,0.04357,0.0392257,0.0354335,0.0321082, &
     0.02918,0.0265913,0.0242943,0.0222492,0.0204225,0.0187859, &
     0.0173155,0.0159908,0.0147943,0.0137111,0.0127282,0.0118343, &
     0.0110197,0.0102759,0.00959549, &
     0.0553807,0.0495037,0.0444239,0.0400097,0.0361551,0.0327736, &
     0.0297949,0.0271605,0.0248222,0.0227396,0.0208788,0.0192111, &
     0.0177122,0.0163615,0.0151413,0.0140361,0.013033,0.0121206, &
     0.0112888,0.0105292,0.00983409, &
     0.0563738,0.0504116,0.0452562,0.0407745,0.0368593,0.0334235, &
     0.0303958,0.0277171,0.0253387,0.0232197,0.0213257,0.0196277, &
     0.0181013,0.0167252,0.0154817,0.0143552,0.0133325,0.0124019, &
     0.0115534,0.0107783,0.0100688, &
     0.0573406,0.0512963,0.0460676,0.0415206,0.0375468,0.0340583, &
     0.030983,0.0282614,0.0258441,0.0236896,0.0217634,0.020036, &
     0.0184826,0.017082,0.0158158,0.0146685,0.0136266,0.0126783, &
     0.0118135,0.0110232,0.0102998, &
     0.0582822,0.0521584,0.0468589,0.0422486,0.038218,0.0346784, &
     0.0315571,0.0287938,0.0263386,0.0241497,0.0221922,0.0204362, &
     0.0188566,0.0174319,0.0161437,0.0149761,0.0139154,0.0129499, &
     0.0120691,0.0112641,0.0105269, &
     0.0591994,0.0529987,0.0476307,0.042959,0.0388734,0.0352843, &
     0.0321182,0.0293144,0.0268225,0.0246002,0.0226121,0.0208283, &
     0.0192232,0.0177751,0.0164654,0.015278,0.0141991,0.0132167, &
     0.0123204,0.0115009,0.0107504, &
     0.0600932,0.053818,0.0483836,0.0436525,0.0395136,0.0358764, &
     0.0326669,0.0298237,0.0272961,0.0250413,0.0230236,0.0212126, &
     0.0195826,0.0181118,0.0167811,0.0155744,0.0144778,0.0134789, &
     0.0125673,0.0117338,0.0109702, &
     0.0609642,0.0546169,0.0491183,0.0443295,0.0401388,0.036455, &
     0.0332033,0.030322,0.0277596,0.0254732,0.0234266,0.0215892, &
     0.0199351,0.018442,0.0170909,0.0158654,0.0147514,0.0137365, &
     0.0128101,0.0119627,0.0111863 &
     /), (/21,21/))

contains

!*******************************************************************************
!
! Subroutine:  final_gb_setup
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine final_gb_setup(atm_cnt)

  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer                       :: atm_cnt

! Local variables:

  integer                       :: alloc_failed
  integer                       :: i

  if( allocated(reff)) deallocate(reff)
  if( allocated(psi)) deallocate(psi)
  if( allocated(rjx)) deallocate(rjx)
  if( allocated(r2x)) deallocate(r2x)
  if( allocated(sumdeijda)) deallocate(sumdeijda)
  if( allocated(vectmp1)) deallocate(vectmp1)
  if( allocated(vectmp2)) deallocate(vectmp2)
  if( allocated(vectmp3)) deallocate(vectmp3)
  if( allocated(vectmp4)) deallocate(vectmp4)
  if( allocated(vectmp5)) deallocate(vectmp5)
  if( allocated(jj)) deallocate(jj)
  if( allocated(skipv)) deallocate(skipv)

  allocate(reff(atm_cnt), &
           psi(atm_cnt), &
           rjx(atm_cnt), &
           r2x(atm_cnt), &
           sumdeijda(atm_cnt), &
           vectmp1(atm_cnt), &
           vectmp2(atm_cnt), &
           vectmp3(atm_cnt), &
           vectmp4(atm_cnt), &
           vectmp5(atm_cnt), &
           jj(atm_cnt), &
           skipv(atm_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  if (rbornstat .ne. 0) then

	if( allocated(gbl_rbmax)) deallocate(gbl_rbmax)
	if( allocated(gbl_rbmin)) deallocate(gbl_rbmin)
	if( allocated(gbl_rbave)) deallocate(gbl_rbave)
	if( allocated(gbl_rbfluct)) deallocate(gbl_rbfluct)

    allocate(gbl_rbmax(atm_cnt), &
             gbl_rbmin(atm_cnt), &
             gbl_rbave(atm_cnt), &
             gbl_rbfluct(atm_cnt), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    gbl_rbmax(:) = 0.d0
    gbl_rbmin(:) = 999.d0
    gbl_rbave(:) = 0.d0
    gbl_rbfluct(:) = 0.d0

  end if

  if (igb .eq. 7) then

	if( allocated(neck_idx)) deallocate(neck_idx)

    allocate(neck_idx(atm_cnt), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    ! Some final error checking before run start for igb 7

    do i = 1, atm_cnt

      neck_idx(i) = dnint((atm_gb_radii(i) - 1.d0) * 20.d0)

      if (neck_idx(i) .lt. 0 .or. neck_idx(i) .gt. 20) then

        if (master) then
          write(error_msg, '(a,a,i6,a,f7.3,a,a,a,a)') error_hdr, 'Atom ', i, &
            ' has a radius (', atm_gb_radii(i), ') outside the allowed range of', &
               extra_line_hdr, &
            '1.0 - 2.0 angstrom for igb=7. ', &
            'Regenerate prmtop with bondi radii.'
        end if

        call mol_mech_error

      end if

    end do

  end if 

  return

end subroutine final_gb_setup

!*******************************************************************************
!
! Subroutine:  gb_ene
!
! Description: Calculate forces, energies based on Generalized Born.
!
!   Compute nonbonded interactions with a generalized Born model,
!   getting the "effective" Born radii via the approximate pairwise method
!   Use Eqs 9-11 of Hawkins, Cramer, Truhlar, J. Phys. Chem. 100:19824
!   (1996).  Aside from the scaling of the radii, this is the same
!   approach developed in Schaefer and Froemmel, JMB 216:1045 (1990).
!   
!   The input coordinates are in the "x" array, and the forces in "f"
!   get updated; energy components are returned in "egb", "eelt" and
!   "evdw".
!   
!   Input parameters for the generalized Born model are "rborn(i)", the
!   intrinsic dielectric radius of atom "i", and "fs(i)", which is
!   set (in init_prmtop_dat()) to (rborn(i) - offset)*si.
!   
!   Input parameters for the "gas-phase" electrostatic energies are
!   the charges, in the "charge()" array.
!   
!   Input parameters for the van der Waals terms are "cn1()" and "cn2()",
!   containing LJ 12-6 parameters, and "asol" and "bsol" containing
!   LJ 12-10 parameters.  (The latter are not used in 1994 and later
!   forcefields.)  The "iac" and "ico" arrays are used to point into
!   these matrices of coefficients.
!   
!   The "numex" and "natex" arrays are used to find "excluded" pairs of
!   atoms, for which gas-phase electrostatics and LJ terms are skipped;
!   note that GB terms are computed for all pairs of atoms.
!   
!   The code also supports a multiple-time-step facility in which:
!   
!   Pairs closer than sqrt(cut_inner) are evaluated every nrespai steps, pairs
!   between sqrt(cut_inner) and sqrt(cut) are evaluated every nrespa steps,
!   and pairs beyond sqrt(cut) are ignored
!   
!   The forces arising from the derivatives of the GB terms with respect
!   to the effective Born radii are evaluated every nrespa steps.
!   
!   The surface-area dependent term is evaluated every nrespa steps.
!   
!   The effective radii are only updated every nrespai steps
!   
!   (Be careful with the above: what seems to work is dt=0.001,
!    nrespai=2, nrespa=4; anything beyond this seems dangerous.)
!   
!   igb = 20  - screen coulomb interactions (added by Igor Kurnikov)
!
!   Written 1999-2000, primarily by D.A. Case, with help from C. Brooks,
!   T. Simonson, R. Sinkovits  and V. Tsui.  The LCPO implementation
!   was written by V. Tsui.
!   
!   Vectorization and optimization 1999-2000, primarily by C. P. Sosa,
!   T. Hewitt, and D. A. Case.  Work presented at CUG Fall of 2000.
!
!   NOTE - in the old sander code, the Generalized Born energy was calc'd
!          and returned as epol; here we rename this egb and will pass it
!          all the way out to the run* routines with this name.
!
!*******************************************************************************

subroutine gb_ene(crd, frc, rborn, fs, charge, iac, ico, numex, &
                  natex, atm_cnt, natbel, egb, eelt, evdw, irespa)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use timers_mod

  implicit none

! Formal arguments:

  double precision      :: crd(*)
  double precision      :: frc(*)
  double precision      :: rborn(*)
  double precision      :: fs(*)
  double precision      :: charge(*)
  integer               :: iac(*)
  integer               :: ico(*)
  integer               :: numex(*)
  integer               :: natex(*)
  integer               :: atm_cnt
  integer               :: natbel
  double precision      :: egb, eelt, evdw
  integer, intent(in)   :: irespa

! Local variables:

  double precision      :: cut2, cut_inner2             !
  double precision      :: extdiel_inv                  !
  double precision      :: intdiel_inv                  !
  double precision      :: ri, rj                       !
  double precision      :: ri1i
  double precision      :: xij, yij, zij                !
  double precision      :: dij1i, dij2i, dij3i          !
  double precision      :: r2                           !
  double precision      :: dij                          !
  double precision      :: sj, sj2                      !
  double precision      :: frespa                       !
  double precision      :: qi, qiqj                     !
  double precision      :: dumx, dumy, dumz             !
  double precision      :: fgbi                         !
  double precision      :: rinv, r2inv, r6inv, r10inv   !
  double precision      :: fgbk                         !
  double precision      :: expmkf                       !
  double precision      :: dl                           !
  double precision      :: de                           !
  double precision      :: e                            !
  double precision      :: temp1                        !
  double precision      :: temp4, temp5, temp6, temp7   !
  double precision      :: eel                          !
  double precision      :: f6, f12, f10                 !
  double precision      :: dedx, dedy, dedz             !
  double precision      :: qi2h, qid2h                  !
  double precision      :: datmp                        !
  double precision      :: thi, thi2                    !
  double precision      :: f_x, f_y, f_z                !
  double precision      :: f_xi, f_yi, f_zi             !
  double precision      :: xi, yi, zi                   !
  double precision      :: dumbo                        !
  double precision      :: tmpsd                        !
   
  ! Variables needed for smooth integration cutoff in Reff:

  double precision      :: rgbmax1i                     !
  double precision      :: rgbmax2i                     !
  double precision      :: rgbmaxpsmax2                 !

  ! Scratch variables used for calculating neck correction:

  double precision      ::  mdist
  double precision      ::  mdist2
  double precision      ::  mdist3
  double precision      ::  mdist5
  double precision      ::  mdist6

  ! Stuff for alpb:

  double precision              :: alpb_beta
  double precision              :: one_arad_beta
  double precision              :: gb_kappa_inv
  ! Alpha prefactor for alpb_alpha:
  double precision, parameter   :: alpb_alpha = 0.571412d0

  ! Boundary box treatment

  double precision              :: cf_box = 2.d0
 

  integer               :: icount
  integer               :: i, j, k
  integer               :: kk1 
  integer               :: max_i
  integer               :: iaci
  integer               :: iexcl, jexcl
  integer               :: jexcl_last
  integer               :: jjv
  integer               :: ic
  integer               :: j3
  logical               :: onstep
   
  ! FGB taylor coefficients follow
  ! from A to H :
  ! 1/3 , 2/5 , 3/7 , 4/9 , 5/11
  ! 4/3 , 12/5 , 24/7 , 40/9 , 60/11

  double precision, parameter  :: te = 4.d0 / 3.d0
  double precision, parameter  :: tf = 12.d0 / 5.d0
  double precision, parameter  :: tg = 24.d0 / 7.d0
  double precision, parameter  :: th = 40.d0 / 9.d0
  double precision, parameter  :: thh = 60.d0 / 11.d0

  egb = 0.d0
  eelt = 0.d0
  evdw = 0.d0
  
  if (mod(irespa, nrespai) .ne. 0) return
   
  cut2 = gb_cutoff * gb_cutoff
  cut_inner2 = cut_inner * cut_inner
  onstep = mod(irespa, nrespa) .eq. 0
   
  if (alpb .eq. 0 .or. igb .eq. 20) then
    ! Standard Still's GB
    extdiel_inv = 1.d0 / extdiel
    intdiel_inv = 1.d0 / intdiel
  else
    ! Sigalov Onufriev ALPB (epsilon-dependent GB):
    alpb_beta = alpb_alpha * (intdiel / extdiel)
    extdiel_inv = 1.d0 / (extdiel * (1.d0 + alpb_beta))
    intdiel_inv = 1.d0 / (intdiel * (1.d0 + alpb_beta))
    one_arad_beta = alpb_beta / arad
    if (gb_kappa .ne. 0.d0) gb_kappa_inv = 1.d0 / gb_kappa
  end if

  max_i = atm_cnt
  if (natbel .gt. 0) max_i = natbel

  ! Smooth "cut-off" in calculating GB effective radii.  
  ! Implemented by Andreas Svrcek-Seiler and Alexey Onufriev. 
  ! The integration over solute is performed up to rgbmax and includes 
  ! parts of spheres; that is an atom is not just "in" or "out", as 
  ! with standard non-bonded cut.  As a result, calculated effective 
  ! radii are less than rgbmax. This saves time, and there is no 
  ! discontinuity in dReff / drij.

  ! Only the case rgbmax > 5*max(sij) = 5*gb_fs_max ~ 9A is handled; this is 
  ! enforced in mdread().  Smaller values would not make much physical
  ! sense anyway.
      
  rgbmax1i = 1.d0 / rgbmax
  rgbmax2i = rgbmax1i * rgbmax1i
  rgbmaxpsmax2 = (rgbmax + gb_fs_max)**2

  !---------------------------------------------------------------------------
  ! Step 1: loop over pairs of atoms to compute the effective Born radii.
  !---------------------------------------------------------------------------
   
  if ( (irespa .lt. 2 .or. mod(irespa, nrespai) .eq. 0) .and. (igb .ne. 20)) &
    call calc_born_radii(atm_cnt, crd, fs, rborn)

  !--------------------------------------------------------------------------
  ! 
  ! Step 2: Loop over all pairs of atoms, computing the gas-phase
  !         electrostatic energies, the LJ terms, and the off-diagonal
  !         GB terms.  Also accumulate the derivatives of these off-
  !         diagonal terms with respect to the inverse effective radii,
  !         sumdeijda(k) will hold  sum over i, j>i (deij / dak),  where
  !         "ak" is the inverse of the effective radius for atom "k".
  ! 
  !         Update the forces with the negative derivatives of the
  !         gas-phase terms, plus the derivatives of the explicit
  !         distance dependence in Fgb, i.e. the derivatives of the
  !         GB energy terms assuming that the effective radii are constant.
  ! 
  !--------------------------------------------------------------------------
   
  sumdeijda(1:atm_cnt) = 0.d0
  
  ! Note: this code assumes that the belly atoms are the first natbel
  !       atoms...this is checked in mdread.
   
  iexcl = 1
  if( numtasks .gt. 1) then
    do i = 1, mytaskid
      iexcl = iexcl + numex(i)
    end do
  endif

  do i = mytaskid + 1, max_i, numtasks
	
    xi = crd(3 * i - 2)
    yi = crd(3 * i - 1)
    zi = crd(3 * i)
    qi = charge(i)
    ri = reff(i)
    iaci = ntypes * (iac(i) - 1)
    jexcl = iexcl
    jexcl_last = iexcl + numex(i) - 1

    dumx = 0.d0
    dumy = 0.d0
    dumz = 0.d0
      
    ! check the exclusion list for eel and vdw:
      
    do k = i + 1, atm_cnt
      skipv(k) = .false.
    end do
    do jjv = jexcl, jexcl_last
      if( natex(jjv) .gt. 0 ) then
        skipv(natex(jjv)) = .true.
      endif
    end do
      
    icount = 0
    do j = i + 1, atm_cnt
         
      xij = xi - crd(3 * j - 2)
      yij = yi - crd(3 * j - 1)
      zij = zi - crd(3 * j)
      r2 = xij * xij + yij * yij + zij * zij
      if (r2 .gt. cut2) cycle
      if (.not. onstep .and. r2 .gt. cut_inner2) cycle
         
      icount = icount + 1
      jj(icount) = j             ! make a short list of atoms j that are within the cutoff radius
      r2x(icount) = r2           ! (r_ij)^2  for j in the short list
      rjx(icount) = reff(j)      ! effective radii of j atoms in the short list
         
    end do
      
    if( igb .ne.20 ) then  ! not screened coulomb interactions
      vectmp1(1:icount) = 4.d0 * ri * rjx(1:icount)
      call vdinv(icount, vectmp1, vectmp1)
      vectmp1(1:icount) = -r2x(1:icount) * vectmp1(1:icount)
      call vdexp(icount, vectmp1, vectmp1)
      ! vectmp1 now contains exp(-rij^2/[4*ai*aj])
      vectmp3(1:icount) = r2x(1:icount) + rjx(1:icount) * ri * vectmp1(1:icount)
      ! vectmp3 now contains fij
    else
      vectmp3(1:icount) = r2x(1:icount) ! fij^2 = rij^2    for screened coulomb interactions
    endif
    
    call vdinvsqrt(icount, vectmp3, vectmp2)
    ! vectmp2 now contains 1/fij   or 1/rij  if igb == 20

    if (gb_kappa .ne. 0.d0) then
      call vdinv(icount, vectmp2, vectmp3)   
      vectmp3(1:icount) = -gb_kappa * vectmp3(1:icount)
      call vdexp(icount, vectmp3, vectmp4)
      ! vectmp4 now contains exp(-kappa*fij)  or exp(-kappa*rij) for (igb == 20)
    end if

    call vdinvsqrt(icount, r2x, vectmp5) ! 1/rij

    ! vectmp1 = exp(-rij^2/[4*ai*aj])
    ! vectmp2 = 1/fij
    ! vectmp3 = -kappa*fij - if kappa .ne. 0.d0, otherwise .eq. fij
    ! vectmp4 = exp(-kappa*fij)  or exp(-kappa*rij) when (igb == 20)
    ! vectmp5 = 1/rij

!    write(*,*) "1/gb_kappa = ", 1.0d0/gb_kappa 

    ! Start first outer loop
    !dir$ ivdep
    do k = 1, icount

      j = jj(k)
      xij = xi - crd(3 * j - 2)
      yij = yi - crd(3 * j - 1)
      zij = zi - crd(3 * j)
      r2 = r2x(k)
      qiqj = qi * charge(j)

	  if( igb .ne. 20) then  ! not simple screen coulomb interactions

         if (gb_kappa .eq. 0.d0) then
           fgbk = 0.d0
           expmkf = extdiel_inv
         else
           expmkf = vectmp4(k) * extdiel_inv
           fgbk = vectmp3(k)*expmkf !-kappa*fij*exp(-kappa*fij)/Eout
           if (alpb .eq. 1) &
             fgbk = fgbk + (fgbk * one_arad_beta * (-vectmp3(k) * gb_kappa_inv))
             ! (-kappa*fij*exp(-kappa*fij)(1 + fij*ab/A)/Eout)*(1/fij+ab/A)
             ! Note: -vectmp2(k)*kappa_inv = fij
         end if

         dl = intdiel_inv - expmkf
         fgbi = vectmp2(k) ! 1.d0/fij

         if (alpb .eq. 0) then
           e = -qiqj * dl * fgbi
         else
           e = -qiqj * dl * (fgbi + one_arad_beta)
         end if

         egb = egb + e

         temp4 = fgbi * fgbi * fgbi ! 1.d0/fij^3

         ! [here, and in the gas-phase part, "de" contains -(1/r)(dE/dr)]
         
         temp6 = -qiqj * temp4 * (dl + fgbk)

         ! -qiqj/fij^3*[1/Ein - e(-Kfij)/Eout) -kappa*fij*
         ! exp(-kappa*fij)(1 + fij*a*b/A ) /Eout]

         temp1 = vectmp1(k) ! exp(-rij^2/[4*ai*aj])

         de = temp6 * (1.d0 - 0.25d0 * temp1)

         rj = rjx(k)

         temp5 = 0.5d0 * temp1 * temp6 * (ri * rj + 0.25d0 * r2)

         sumdeijda(i) = sumdeijda(i) + ri * temp5
         sumdeijda(j) = sumdeijda(j) + rj * temp5
      
      else ! ( igb == 20)
         de = 0.0d0
      endif  ! -  end of  ( igb != 20) condition
         
      ! skip exclusions for remaining terms:
         
      if (.not. skipv(j)) then
            
        ! gas-phase Coulomb energy:
            
        rinv = vectmp5(k) ! 1.d0/rij
        r2inv = rinv * rinv
        
        if(igb.eq.20) then
			eel = extdiel_inv * qiqj * rinv * vectmp4(k)
			eelt = eelt + eel
			de = de + eel * r2inv + gb_kappa*eel*rinv
		else
			eel = intdiel_inv * qiqj * rinv
			eelt = eelt + eel
			de = de + eel * r2inv
		endif
        
 !       write(*,*) "i,j, d_ecoul = ",i,j,eel * r2inv !  TEMPORAL CHECK IGOR
            
        ! van der Waals energy:
            
        ic = ico(iaci + iac(j))
        if (ic .gt. 0) then
          ! 6-12 potential:
          r6inv = r2inv * r2inv * r2inv
          f6 = gbl_cn2(ic) * r6inv
          f12 = gbl_cn1(ic) * (r6inv * r6inv)
      
		  if( igb .ne. 20) then
			evdw = evdw + (f12 - f6)    
			de = de + (12.d0 * f12 - 6.d0 * f6) * r2inv
		  else  ! TEMPORAL FIX (IGOR) TO MODEL DE PABLO ATOM EXCLUSION INTERACTIONS TYPED FIXED PARAMS: 
		    if( rinv*6.86 .gt. 1.0d0 )then
				evdw = evdw + (f12 - f6) 
				evdw = evdw + 0.26                      
				de = de + (12.d0 * f12 - 6.d0 * f6) * r2inv
			else
				f12 = 0.0d0
				f6  = 0.0d0
			endif
		  endif 
        else if( ic .lt. 0 ) then
          ! 10-12 potential:
          r10inv = r2inv * r2inv * r2inv * r2inv * r2inv
          f10 = gbl_bsol(-ic) * r10inv
          f12 = gbl_asol(-ic) * r10inv * r2inv
          evdw = evdw + f12 - f10
          de = de + (12.d0 * f12 - 10.d0 * f10) * r2inv
        end if  ! (ic .gt. 0)
      end if  ! (.not. skipv(j))
         
      ! derivatives:
         
      if (onstep .and. r2 .gt. cut_inner2 ) then
        de = de * nrespa
      else
        de = de * nrespai
      end if
         
      dedx = de * xij
      dedy = de * yij
      dedz = de * zij
      dumx = dumx + dedx
      dumy = dumy + dedy
      dumz = dumz + dedz
      
 !     if( (i .eq. 15 .and. j .eq. 30) .or. (i .eq. 30 .and. j .eq. 15) ) then
 !       write(*,*) " for i,j = 15,30"
 !       write(*,*) "dedx,dedy,dedz, dumx,dumy,dumz =",dedx,dedy,dedz, dumx,dumy,dumz
 !     endif
      
      frc(3 * j - 2) = frc(3 * j - 2) - dedx
      frc(3 * j - 1) = frc(3 * j - 1) - dedy
      frc(3 * j) = frc(3 * j) - dedz
    end do
      
    frc(3 * i - 2) = frc(3 * i - 2) + dumx
    frc(3 * i - 1) = frc(3 * i - 1) + dumy
    frc(3 * i) = frc(3 * i) + dumz
    
    if( numtasks .gt. 1) then
        do k = i, min(i + numtasks - 1, atm_cnt)
            iexcl = iexcl + numex(k)
        end do
    else
        iexcl = iexcl + numex(i)
    end if
      
! boundary box treatment:
	if( (bbox_xmax - bbox_xmin) .gt. 1.0d0) then
		if( xi .lt. bbox_xmin ) then
			frc(3 * i - 2) = frc(3 * i - 2) + cf_box*(bbox_xmin - xi)
			evdw = evdw + 0.5d0*cf_box*(bbox_xmin - xi)*(bbox_xmin - xi)
		endif
		if( xi .gt. bbox_xmax ) then
			frc(3 * i - 2) = frc(3 * i - 2) - cf_box*(xi - bbox_xmax)
			evdw = evdw + 0.5d0*cf_box*(bbox_xmax - xi)*(bbox_xmax - xi)
		endif
	endif
	if( (bbox_ymax - bbox_ymin) .gt. 1.0d0) then
		if( yi .lt. bbox_ymin ) then
			frc(3 * i - 1) = frc(3 * i - 1) + cf_box*(bbox_ymin - yi)
			evdw = evdw + 0.5d0*cf_box*(bbox_ymin - yi)*(bbox_ymin - yi)
		endif
		if( yi .gt. bbox_ymax ) then
			frc(3 * i - 1) = frc(3 * i - 1) - cf_box*(yi - bbox_ymax)
			evdw = evdw + 0.5d0*cf_box*(bbox_ymax - yi)*(bbox_ymax - yi)
		endif	
	endif
	if( (bbox_zmax - bbox_zmin) .gt. 1.0d0) then
		if( zi .lt. bbox_zmin ) then
			frc(3 * i) = frc(3 * i) + cf_box*(bbox_zmin - zi)
			evdw = evdw + 0.5d0*cf_box*(bbox_zmin - zi)*(bbox_zmin - zi)
		endif
		if( zi .gt. bbox_zmax ) then
			frc(3 * i) = frc(3 * i) - cf_box*(zi - bbox_zmax)
			evdw = evdw + 0.5d0*cf_box*(bbox_zmax - zi)*(bbox_zmax - zi)
		endif
	endif	

  end do  !  i = 1, max_i

  call update_gb_time(calc_gb_offdiag_timer)
   
  if( igb .eq. 20) return
   
  !--------------------------------------------------------------------------
  ! 
  ! Step 3:  Finally, do the reduction over the sumdeijda terms:, adding
  !          into the forces those terms that involve derivatives of
  !          the GB terms (including the diagonal or "self" terms) with
  !          respect to the effective radii.  This is done by computing
  !          the vector dai / dxj, and using the chain rule with the
  !          previously-computed sumdeijda vector.
  ! 
  !          Do these terms only at "nrespa" multiple-time step intervals;
  !          (when igb=2 or 5, one may need to do this at every step)
  ! 
  !--------------------------------------------------------------------------
   
  if (onstep) then

    if( numtasks .gt. 1) then
      
    ! first, collect all the sumdeijda terms:
      
      call mpi_allreduce(sumdeijda, vectmp1, atm_cnt, mpi_double_precision, &
                         mpi_sum, lib_mpi_comm, err_code_mpi)

      sumdeijda(1:atm_cnt) = vectmp1(1:atm_cnt)

      call update_gb_time(dist_gb_rad_timer)
    endif
      
    frespa = nrespa

    ! diagonal egb term, plus off-diag derivs wrt alpha .eq. reff^-1:
      
    do i = mytaskid + 1, max_i, numtasks
         
      f_xi = 0.d0
      f_yi = 0.d0
      f_zi = 0.d0
      qi = charge(i)
      expmkf = exp(-gb_kappa * reff(i)) * extdiel_inv
      dl = intdiel_inv - expmkf
      qi2h = 0.5d0 * qi * qi
      qid2h = qi2h * dl

      if (alpb .eq. 0) then
        egb = egb - qid2h / reff(i)
        temp7 = -sumdeijda(i) + qid2h - gb_kappa * qi2h * expmkf * reff(i)
      else
        egb = egb - qid2h * (1.d0/reff(i) + one_arad_beta)
        temp7 = -sumdeijda(i) + qid2h - gb_kappa * qi2h * expmkf * reff(i) * &
                (1.d0 + one_arad_beta * reff(i))
      end if

      xi = crd(3 * i - 2)
      yi = crd(3 * i - 1)
      zi = crd(3 * i)
      ri = rborn(i) - offset
      ri1i = 1.d0 / ri
      iaci = ntypes * (iac(i) - 1)
         
      if (igb .eq. 2 .or. igb .eq. 5 .or. igb .eq. 7) then
            
        ! new onufriev: we have to later scale values by a
        !               alpha,beta,gamma -dependent factor:
            
        ri = rborn(i) - offset
        thi = tanh((gb_alpha + gb_gamma * psi(i) * psi(i) - &
              gb_beta * psi(i)) * psi(i))

        thi2 = (gb_alpha + 3.d0 * gb_gamma * psi(i) * psi(i) - &
               2.d0 * gb_beta * psi(i)) * (1.d0 - thi * thi) * ri / rborn(i)
      end if
         
      icount = 0
      do j = 1, atm_cnt
        if (i .eq. j) cycle
            
        xij = xi - crd(3 * j - 2)
        yij = yi - crd(3 * j - 1)
        zij = zi - crd(3 * j)
        r2 = xij * xij + yij * yij + zij * zij
        if (r2 .gt. rgbmaxpsmax2) cycle

        ! pairlist contains only atoms within rgbmax + safety margin
            
        icount = icount + 1
        jj(icount) = j
        r2x(icount) = r2
            
      end do

      call vdinvsqrt(icount, r2x, vectmp1)
         
      kk1 = 0
      do k = 1, icount
        j = jj(k)
        r2 = r2x(k)
        sj =  fs(j)

        dij1i = vectmp1(k)
        dij = r2 * dij1i
        sj2 = sj * sj
            
        if (dij .gt. 4.d0 * sj) cycle
        kk1 = kk1 + 1
        vectmp3(kk1) = dij + sj
        if (dij .gt. ri + sj) then
          vectmp2(kk1) = r2 - sj2
          vectmp4(kk1) = dij - sj
        else if (dij .gt. abs(ri - sj)) then
          vectmp2(kk1) = dij + sj
          vectmp4(kk1) = ri
        else if (ri .lt. sj) then
          vectmp2(kk1) = r2 - sj2
          vectmp4(kk1) = sj - dij
        else
          vectmp2(kk1) = 1.d0
          vectmp4(kk1) = 1.d0
        end if
      end do
         
      call vdinv(kk1, vectmp2, vectmp2)
      call vdinv(kk1, vectmp3, vectmp3)
      vectmp4(1:kk1) = vectmp4(1:kk1) * vectmp3(1:kk1)
      call vdln(kk1, vectmp4, vectmp4)
         
      kk1 = 0
      do k = 1, icount
        j = jj(k)
        j3 = 3 * j
        r2 = r2x(k)
        xij = xi - crd(j3 - 2)
        yij = yi - crd(j3 - 1)
        zij = zi - crd(j3)
            
        dij1i = vectmp1(k)
        dij = r2 * dij1i
        sj = fs(j)
        if (dij .gt. rgbmax + sj) cycle
        sj2 = sj * sj
            
        ! datmp will hold (1/r)(dai/dr):
            
        dij2i = dij1i * dij1i
        dij3i = dij2i * dij1i

        if (dij .gt. rgbmax - sj) then 

          temp1 = 1.d0 / (dij - sj)
          datmp = 0.125d0 * dij3i * ((r2 + sj2) * &
                  (temp1 * temp1 - rgbmax2i) - 2.d0 * log(rgbmax * temp1))
            
        else if (dij .gt. 4.d0 * sj) then
               
          tmpsd = sj2 * dij2i
          dumbo = te + tmpsd * (tf + tmpsd * (tg + tmpsd * (th + tmpsd * thh)))
          datmp = tmpsd * sj * dij2i * dij2i * dumbo
               
        else if (dij .gt. ri + sj) then
               
          kk1 = kk1 + 1
          datmp = vectmp2(kk1) * sj * (-0.5d0 * dij2i + vectmp2(kk1)) + &
                  0.25d0 * dij3i * vectmp4(kk1)
               
        else if (dij .gt. abs(ri - sj)) then

          kk1 = kk1 + 1
          datmp = -0.25d0 * (-0.5d0 * (r2 - ri * ri + sj2) * &
                  dij3i * ri1i * ri1i + dij1i * vectmp2(kk1) * &
                  (vectmp2(kk1) - dij1i) - dij3i * vectmp4(kk1))
               
        else if (ri .lt. sj) then

          kk1 = kk1 + 1
          datmp = -0.5d0 * (sj * dij2i * vectmp2(kk1) - &
                  2.d0 * sj * vectmp2(kk1) * vectmp2(kk1) - &
                  0.5d0 * dij3i * vectmp4(kk1))
               
        else

          kk1 = kk1 + 1
          datmp = 0.d0

        end if  ! (dij .gt. 4.d0 * sj)
            
        if (igb .eq. 7) then

          if (dij .lt. rborn(i) + rborn(j) + gb_neckcut) then

            ! Derivative of neck with respect to dij is:
            !                     5
            !              9 mdist
            !   (2 mdist + --------) neckMaxVal gb_neckscale
            !                 5
            ! -(------------------------)
            !                        6
            !             2   3 mdist  2
            !   (1 + mdist  + --------)
            !                    10

            mdist = dij - neckMaxPos(neck_idx(i), neck_idx(j))
            mdist2 = mdist * mdist
            mdist3 = mdist2 * mdist
            mdist5 = mdist2 * mdist3
            mdist6 = mdist3 * mdist3

            ! temp1 will be divisor of above fraction * dij
            ! (datmp is deriv * 1/r)

            temp1 = 1.d0 + mdist2 + (0.3d0) * mdist6
            temp1 = temp1 * temp1 * dij

            ! (Note "+" means subtracting derivative, since above 
            !     expression has leading "-")

            datmp = datmp + ((2.d0 * mdist + (9.d0/5.d0) * mdist5) * &
                    neckMaxVal(neck_idx(i), neck_idx(j)) * &
                    gb_neckscale) / temp1

          end if ! if (dij < rborn(i) +rborn(j) + gb_neckcut)

        end if

        datmp = -datmp * frespa * temp7

        if (igb .eq. 2 .or. igb .eq. 5 .or. igb .eq. 7) datmp = datmp * thi2
            
        f_x = xij * datmp
        f_y = yij * datmp
        f_z = zij * datmp
        frc(j3 - 2) = frc(j3 - 2) + f_x
        frc(j3 - 1) = frc(j3 - 1) + f_y
        frc(j3) = frc(j3) + f_z
        f_xi = f_xi - f_x
        f_yi = f_yi - f_y
        f_zi = f_zi - f_z
            
      end do  !  k = 1, icount
         
      frc(3 * i - 2) = frc(3 * i - 2) + f_xi
      frc(3 * i - 1) = frc(3 * i - 1) + f_yi
      frc(3 * i) = frc(3 * i) + f_zi
         
    end do   ! end loop over atom i
      
    call update_gb_time(calc_gb_diag_timer)

  end if  !  end:  if (onstep) 

  return

end subroutine gb_ene

!*******************************************************************************
!
! Subroutine:  calc_born_radii
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine calc_born_radii(atm_cnt, crd, fs, rborn)

  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use prmtop_dat_mod
  use timers_mod

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3*atm_cnt)
  double precision      :: fs(atm_cnt)
  double precision      :: rborn(atm_cnt)

! Local variables:

  double precision      :: ri, rj
  double precision      :: ri1i, rj1i
  double precision      :: xij, yij, zij
  double precision      :: dij1i, dij2i
  double precision      :: r2
  double precision      :: dij
  double precision      :: si, si2
  double precision      :: sj, sj2
  double precision      :: theta
  double precision      :: uij
  double precision      :: xi, yi, zi
  double precision      :: reff_i
  double precision      :: dumbo
  double precision      :: tmpsd

  ! Variables needed for smooth integration cutoff in Reff:

  double precision      :: rgbmax1i
  double precision      :: rgbmax2i
  double precision      :: rgbmaxpsmax2

  ! Scratch variables used for calculating neck correction:

  double precision      ::  mdist
  double precision      ::  mdist2
  double precision      ::  mdist3
  double precision      ::  mdist6
  double precision      ::  neck

  integer               :: icount

  integer               :: i, j, k
  integer               :: kk1, kk2 

  ! FGB taylor coefficients follow
  ! from A to H :
  ! 1/3 , 2/5 , 3/7 , 4/9 , 5/11
  ! 4/3 , 12/5 , 24/7 , 40/9 , 60/11

  double precision, parameter  :: ta = 1.d0 / 3.d0
  double precision, parameter  :: tb = 2.d0 / 5.d0
  double precision, parameter  :: tc = 3.d0 / 7.d0
  double precision, parameter  :: td = 4.d0 / 9.d0
  double precision, parameter  :: tdd = 5.d0 / 11.d0

  ! Smooth "cut-off" in calculating GB effective radii.  
  ! Implemented by Andreas Svrcek-Seiler and Alexey Onufriev. 
  ! The integration over solute is performed up to rgbmax and includes 
  ! parts of spheres; that is an atom is not just "in" or "out", as 
  ! with standard non-bonded cut.  As a result, calculated effective 
  ! radii are less than rgbmax. This saves time, and there is no 
  ! discontinuity in dReff / drij.

  ! Only the case rgbmax > 5*max(sij) = 5*gb_fs_max ~ 9A is handled; this is 
  ! enforced in mdread().  Smaller values would not make much physical
  ! sense anyway.
      
  rgbmax1i = 1.d0 / rgbmax
  rgbmax2i = rgbmax1i * rgbmax1i
  rgbmaxpsmax2 = (rgbmax + gb_fs_max)**2

  reff(1:atm_cnt) = 0.d0
   
  do i = mytaskid + 1, atm_cnt, numtasks

    xi = crd(3 * i - 2)
    yi = crd(3 * i - 1)
    zi = crd(3 * i)
      
    reff_i = reff(i)
    ri = rborn(i) - offset
    ri1i = 1.d0 / ri
    si = fs(i)
    si2 = si * si
      
    ! Here, reff_i will sum the contributions to the inverse effective
    ! radius from all of the atoms surrounding atom "i"; later the
    ! inverse of its own intrinsic radius will be added in
      
    icount = 0

    do j = i + 1, atm_cnt
      xij = xi - crd(3 * j - 2)
      yij = yi - crd(3 * j - 1)
      zij = zi - crd(3 * j)
      r2 = xij * xij + yij * yij + zij * zij
      if (r2 .gt. rgbmaxpsmax2) cycle
      icount = icount + 1
      jj(icount) = j
      r2x(icount) = r2
    end do
      
    call vdinvsqrt(icount, r2x, vectmp1)
      
    kk1 = 0
    kk2 = 0
    !dir$ ivdep
    do k = 1, icount
         
      j = jj(k)
      r2 = r2x(k)
      sj = fs(j)

      ! don't fill the remaining vectmp arrays if atoms don't see each other:

      dij1i = vectmp1(k)
      dij = r2 * dij1i
      if (dij .gt. rgbmax + si .and. dij .gt. rgbmax + sj) cycle
      rj = rborn(j) - offset

      if (dij .le. 4.d0 * sj) then
        kk1 = kk1 + 1
        vectmp2(kk1) = dij + sj
        if (dij .gt. ri + sj) then
          vectmp4(kk1) = dij - sj
        else if (dij .gt. abs(ri - sj)) then
          vectmp4(kk1) = ri
        else if (ri .lt. sj) then
          vectmp4(kk1) = sj - dij
        else
          vectmp4(kk1) = 1.d0
        end if
      end if
         
      if (dij .le. 4.d0 * si) then
        kk2 = kk2 + 1
        vectmp3(kk2) = dij + si
        if (dij .gt. rj + si) then
          vectmp5(kk2) = dij - si
        else if (dij .gt. abs(rj - si)) then
          vectmp5(kk2) = rj
        else if (rj .lt. si) then
          vectmp5(kk2) = si - dij
        else
          vectmp5(kk2) = 1.d0
        end if
      end if
         
    end do  !  k = 1, icount
      
    call vdinv(kk1, vectmp2, vectmp2)
    call vdinv(kk2, vectmp3, vectmp3)
    vectmp4(1:kk1) = vectmp2(1:kk1) * vectmp4(1:kk1)
    vectmp5(1:kk2) = vectmp3(1:kk2) * vectmp5(1:kk2)
    call vdln(kk1, vectmp4, vectmp4)
    call vdln(kk2, vectmp5, vectmp5)
      
    kk1 = 0
    kk2 = 0
    do k = 1, icount
         
      j = jj(k)
      r2 = r2x(k)
         
      rj = rborn(j) - offset
      rj1i = 1.d0 / rj
      sj = fs(j)
         
      sj2 = sj * sj
         
      xij = xi - crd(3 * j - 2)
      yij = yi - crd(3 * j - 1)
      zij = zi - crd(3 * j)
         
      dij1i = vectmp1(k)
      dij = r2 * dij1i
         
         
      if (dij .le. rgbmax + sj) then
         
        if ((dij .gt. rgbmax - sj)) then

          uij = 1.d0 / (dij - sj)
          reff_i = reff_i - 0.125d0 * dij1i * (1.d0 + 2.d0 * dij *uij + &
                   rgbmax2i * (r2 - 4.d0 * rgbmax * dij - sj2) + &
                   2.d0 * log((dij - sj) * rgbmax1i))

        else if (dij .gt. 4.d0 * sj) then
            
          dij2i = dij1i * dij1i
          tmpsd = sj2 * dij2i
          dumbo = ta + tmpsd *  (tb + tmpsd * (tc + tmpsd * (td + tmpsd * tdd)))

          reff_i = reff_i - tmpsd * sj * dij2i * dumbo
            
          !     ---following are from the Appendix of Schaefer and Froemmel,
          !        J. Mol. Biol. 216:1045-1066, 1990, divided by (4*Pi):

        else if (dij .gt. ri + sj) then
            
          kk1 = kk1 + 1
          reff_i = reff_i - 0.5d0 * (sj / (r2 - sj2) + 0.5d0 * dij1i * &
                   vectmp4(kk1))
            
          !-----------------------------------------------------------------
            
        else if (dij .gt. abs(ri - sj)) then
            
          kk1 = kk1 + 1
          theta = 0.5d0 * ri1i * dij1i * (r2 + ri * ri - sj2)
          reff_i = reff_i - 0.25d0 * (ri1i * (2.d0 - theta) - &
                   vectmp2(kk1) + dij1i * vectmp4(kk1))
            
          !-----------------------------------------------------------------
            
        else if (ri .lt. sj) then

          kk1 = kk1 + 1
          reff_i = reff_i - 0.5d0 * (sj / (r2 - sj2) + 2.d0 * ri1i + &
                   0.5d0 * dij1i * vectmp4(kk1))
            
          !-----------------------------------------------------------------
            
        else

          kk1 = kk1 + 1

        end if  ! (dij .gt. 4.d0 * sj)

        if (igb .eq. 7) then

          if (dij .lt. rborn(i) + rborn(j) + gb_neckcut) then
            mdist = dij - neckMaxPos(neck_idx(i), neck_idx(j))
            mdist2 = mdist * mdist
            mdist3 = mdist2 * mdist
            mdist6 = mdist3 * mdist3
            neck = neckMaxVal(neck_idx(i), neck_idx(j)) / &
                   (1.d0 + mdist2 + 0.3d0 * mdist6)
            reff_i = reff_i - gb_neckscale * neck
          end if

        end if
         
      end if

      ! --- Now the same thing, but swap i and j:
         
      if (dij .gt. rgbmax + si) cycle
             
      if (dij .gt. rgbmax - si) then

        uij = 1.d0 / (dij - si)
        reff(j) = reff(j) - 0.125d0 * dij1i * (1.d0 + 2.d0 * dij * uij + &
                  rgbmax2i * (r2 - 4.d0 * rgbmax * dij - si2) + &
                  2.d0 * log((dij - si) * rgbmax1i))

      else if (dij .gt. 4.d0 * si) then
            
        dij2i = dij1i * dij1i
        tmpsd = si2 * dij2i
        dumbo = ta + tmpsd * (tb + tmpsd * (tc + tmpsd * (td + tmpsd * tdd)))
        reff(j) = reff(j) - tmpsd * si * dij2i * dumbo
            
      else if (dij .gt. rj + si) then
            
        kk2 = kk2 + 1
        reff(j) = reff(j) - 0.5d0 * (si / (r2 - si2) + &
                  0.5d0 * dij1i * vectmp5(kk2))
            
        !-----------------------------------------------------------------
            
      else if (dij .gt. abs(rj - si)) then
            
        kk2 = kk2 + 1
        theta = 0.5d0 * rj1i * dij1i * (r2 + rj * rj - si2)
        reff(j) = reff(j) - 0.25d0 * (rj1i * (2.d0 - theta) - &
                  vectmp3(kk2) + dij1i * vectmp5(kk2))
            
        !-----------------------------------------------------------------
            
      else if (rj .lt. si) then
            
        kk2 = kk2 + 1
        reff(j) = reff(j) - 0.5d0 * (si / (r2 - si2) + 2.d0 * rj1i + &
                  0.5d0 * dij1i * vectmp5(kk2))
            
        !-----------------------------------------------------------------
            
      else

        kk2 = kk2 + 1

      end if  ! (dij .gt. 4.d0 * si)

      if (igb == 7) then
        if (dij .lt. rborn(j) + rborn(i) + gb_neckcut) then
          mdist = dij - neckMaxPos(neck_idx(j), neck_idx(i))
          mdist2 = mdist * mdist
          mdist3 = mdist2 * mdist
          mdist6 = mdist3 * mdist3
          neck = neckMaxVal(neck_idx(j), neck_idx(i)) / &
                 (1.d0 + mdist2 + 0.3d0 * mdist6)
          reff(j) = reff(j) - gb_neckscale * neck
        end if
      end if

    end do  !  k = 1, icount
      
    ! we are ending the do-i-loop, reassign the scalar to the original array:
      
    reff(i) = reff_i
      
  end do  !  i = 1, atm_cnt  - end of cycle on atoms 
   
  if( numtasks .gt. 1) then
    call update_gb_time(calc_gb_rad_timer)
   
  ! Collect the (inverse) effective radii from other nodes:
   
    call mpi_allreduce(reff, vectmp1, atm_cnt, mpi_double_precision, &
                       mpi_sum, lib_mpi_comm, err_code_mpi)

    reff(1:atm_cnt) = vectmp1(1:atm_cnt)

    call update_gb_time(dist_gb_rad_timer)
  endif
   
  if (igb .eq. 2 .or. igb .eq. 5 .or. igb .eq. 7) then
      
    ! apply the new Onufriev "gbalpha, gbbeta, gbgamma" correction:
      
    do i = 1, atm_cnt
      ri = rborn(i) - offset
      ri1i = 1.d0 / ri
      psi(i) = -ri * reff(i)
      reff(i) = ri1i - tanh((gb_alpha + gb_gamma * psi(i) * psi(i) - &
                gb_beta * psi(i)) * psi(i)) / rborn(i)

      if (reff(i) .lt. 0.d0) reff(i) = 1.d0/30.d0

      reff(i) = 1.d0 / reff(i)
    end do
      
  else
      
    ! "standard" GB, including the "diagonal" term here:
      
    do i = 1, atm_cnt
      ri = rborn(i) - offset
      ri1i = 1.d0 / ri
      reff(i) = 1.d0 / (reff(i) + ri1i)
    end do
  end if
   
  if (rbornstat .eq. 1) then
    do i = 1, atm_cnt
      gbl_rbave(i) = gbl_rbave(i) + reff(i)
      gbl_rbfluct(i) = gbl_rbfluct(i) + reff(i) * reff(i)
      if (gbl_rbmax(i) .le. reff(i)) gbl_rbmax(i) = reff(i)
      if (gbl_rbmin(i) .ge. reff(i)) gbl_rbmin(i) = reff(i)
    end do
  end if

  call update_gb_time(calc_gb_rad_timer)
   
  return

end subroutine calc_born_radii


subroutine print_born_radii_stat(tspan)

    use file_io_dat_mod
    use prmtop_dat_mod

    implicit none
    
    ! Formal arguments:
    
    double precision ::  tspan
    
    ! Local variables

    integer  ::  m

    write(mdout, 590)
    do m = 1, natom
        gbl_rbave(m) = gbl_rbave(m) / tspan
        gbl_rbfluct(m) = gbl_rbfluct(m) / tspan - &
                         gbl_rbave(m) * gbl_rbave(m)
        gbl_rbfluct(m) = sqrt(gbl_rbfluct(m))
        write(mdout, 600) m, gbl_rbmax(m), gbl_rbmin(m), &
                              gbl_rbave(m), gbl_rbfluct(m)
     end do
     
     return

  590 format('ATOMNUM     MAX RAD     MIN RAD     AVE RAD     FLUCT')
  600 format(i4, 2x, 4f12.4) 
  
end subroutine

end module gb_ene_mod
