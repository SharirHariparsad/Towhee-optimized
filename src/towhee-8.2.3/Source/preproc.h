/* 
* MCCCS - Towhee: A Monte Carlo molecular simulation program           *
* Copyright (C) 2003-2020 Marcus G. Martin                             *
* see the file license.gpl for the full license information            *
*                                                                      *
* This program is free software; you can redistribute it and/or        *
* modify it under the terms of the GNU General Public License          *
* as published by the Free Software Foundation; either version 2       *
* of the License, or (at your option) any later version.               *
*                                                                      *
* This program is distributed in the hope that it will be useful,      *
* but WITHOUT ANY WARRANTY; without even the implied warranty of       *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        *
* GNU General Public License for more details.                         *
*                                                                      *
* You should have received a copy of the GNU General Public            *
* License along with this program; if not, write to the Free           *
* Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,      *
* MA  02111-1307, USA.                                                 *
*/

/* last significant modification 03-24-2020 by M.G. Martin */

/* ARRAYS */

/* maxbox = maximum number of simulation boxes (default: 3) */
#define MAXBOX 3

/* numax = maximum number of units/atoms in a molecule (default 2002) */
#define NUMAX 2202

/* ntmax = maximum number of types of molecules (default: 6) */
#define NTMAX 7

/* ntpairmax = maximum number of pairs of types of molecules */
/* this is determined using NTMAX */
#define NTPAIRMAX NTMAX*(NTMAX-1)/2

/* cbmcmaxbendpeaks = maximum number of peaks in the bending multiple */
/* angle distribution */
#define CBMC_MAXBENDPEAKS 3

/* cbmcmaxtorpeaks = maximum number of peaks in a torsion multiple */
/* gaussian distribution (default: 6) */
#define CBMC_MAXTORPEAKS 6

/* cbmcmaxbond = maximum number of parameters to generate a bond */
/* distribution for use in the CBMC moves */
#define CBMCMAXBOND 4

/* cbmcmaxbend = maximum number of bending parameters for */
/* use generating trial angle distributions */
#define CBMCMAXBEND 5

/* cbmcmaxtor = maximum number of torsion parameters for */
/* use generation trial dihedral distributions */
#define CBMCMAXTOR 4

/* maxtwobond = maximum number of parameters to describe the two-bond */
/* biasing functions */
#define MAXTWOBOND 5

/* maxthreebond = maximum number of parameters to describe the */
/* three-bond biasing functions */
#define MAXTHREEBOND 5

/* vectormax = the max number of reciprocal vectors for Ewald sum */
#define MAXKMAX 12
#define VECTORMAX 3*MAXKMAX*MAXKMAX*MAXKMAX

/* maxblock = the maximum number of blocks for a block average */
#define MAXBLOCK 20

/* fldmax = the maximum number fields */
#define FLDMAX 20

/* nntype = number of types of atoms for regular potentials */
#define NNTYPE 375 

/* crosstypemax = maximum number of cross types depending on nntype */
#define CROSSTYPEMAX (NNTYPE*NNTYPE+NNTYPE)/2

/* maxtabtype = maximum number of types of atoms for tabular potentials
   , must be less than or equal to nntype */
#define MAXTABTYPE 10

/* maxtable = maximum number of points for a tabular potential */
#define MAXTABLE 3001

/* tvibmax = maximum number of types of bond stretches */
#define TVIBMAX 8200

/* tbenmax = maximum number of types of bond angles */
#define TBENMAX 2700

/* ttormax = maximum number of types of torsions */
#define TTORMAX 930

/* timpmax = maximum number of types of impropers */
#define TIMPMAX 120

/* taamax = maximum number of angle-angle types */
#define TAAMAX 705

/* tofmax = maximum number of types of special one-five interactions */
#define TOFMAX 5

/* thbondmax = maximum number of types of special hydrogen bond 
   interactions */
#define THBONDMAX 40

  /* tbimax = maximum number of types of bond-increments */
#define TBIMAX 500

/* smallest: the smallest double precision number we can handle */
#define SMALLEST 1d-300

/* values that are simple consequenses of the values up above */
/* there is no reason to modify these values */

/* nnbond = maximum number of bonds from any bead. */
#define NNBOND 6

/* mmbond = maximum number of bonds for use computing the maximum */
/* number of the multi-atom intramolecular terms.  It is safest */
/* to set this equal to nnbond, but setting this to a value higher */
/* than 4 often results in some arrays that are too large to */
/* compile properly */
#define MMBOND 4

/* maxbend = max number of bending angles on any one atom
   computed from nnbond */
#define MAXBEND MMBOND*(MMBOND-1) 

/* maxtor = max number of torsions on any one atom
   computed from nnbond */
#define MAXTOR MMBOND*(MMBOND-1)*(MMBOND-1)

/* maxaa = maximum number of angle-angle terms on any one atom
   computed from nnbond */
#define MAXAA ((MMBOND*(MMBOND-1))**2)/8

/* maximprop = maximum number of improper torsion on any one atom 
   computed from nnbond */
#define MAXIMPROP (MMBOND-1)*(MMBOND-2)

/* maxof = maximum number of special one-five interactions on any one
   atom.  computed from nnbond */
#define MAXOF MMBOND*(MMBOND-1)**3

/* maxdirlength = maximum number of characters in the directory length
   this is used with some of the parallel features */
#define MAXDIRLENGTH 120

/* MAXIMPLICITYPE = maximum number of types of atoms for the implicit 
   force fields */
#define MAXIMPLICITTYPE 25

/* MAXPIVOTBIN = maximum number of bins for the pivot bookkeeping */
#define MAXPIVOTBIN 6

/* NDUMPHIST = maximum allowable number of histogram computations 
   before an output to file */
#define NDUMPHIST 1000

/* MAX_FOREIGN_LAMBDA = maximum number of (lambda_lj,lambda_c) pairs
   for which energy evaluations will take place.  Relevant only for
   Scaled Lennard-Jones classical potential */
#define MAX_FOREIGN_LAMBDA 50

/* size of the character strings for describing force field names */
#define FFNAMELEN 10
/* size of the character strings for bond patterns */
#define BONPATLEN 5
/* size of the character strings for mixing rule */
#define MIXRULLEN 30
/* size of the character string for classical potential */
#define CLAPOTLEN 30
/* size of the character string for element name */
#define ELENAMLEN 2

/* array dimensions for the nonbonded coefficients (nbcoeff) */
#define MINNBCOEFF 0
#define MAXNBCOEFF 55

/* array dimensions for the vibration coefficients (vibcoeff)*/
#define MINVIBCOEFF 0
#define MAXVIBCOEFF 3

/* maximum number of vibration names in vibnames */
#define MAXSAMEV 15

/* number of names required for a vibration */
#define MAXVNAME 2

/* array dimension for the bending angle coefficients (bencoeff) */
#define MINBENDCOEFF -1
#define MAXBENDCOEFF 10

/* maximum number of bending triplets in bendnames */
#define MAXSAMEB 30

/* bending names currently come in triplets */
#define MAXBNAME 3

/* array dimension for the torsion coefficients (torcoeff) */
#define MINTORCOEFF -1
#define MAXTORCOEFF 31

/* maximum number of torsion quartets in tornames */
#define MAXSAMET 90

/* tornames contains quartets of atoms */
#define MAXTNAME 4

/* array dimensions for the improper coefficients (impcoeff) */
#define MINIMPCOEFF 0
#define MAXIMPCOEFF 3

/* maximum number of improper quartets in impnames */
#define MAXSAMEIMP 21

/* impnames contains quartets of atoms */
#define MAXIMPNAME 4

/* array dimensions for the angle-angle coefficients (aacoeff) */
#define MINAACOEFF 0
#define MAXAACOEFF 2

/* maximum number of quartets for angle-angle terms */
#define MAXSAMEAA 1

/* array dimension for the special one-five interactions (ofcoeff) */
#define MINOFCOEFF 0
#define MAXOFCOEFF 4

/* ofnames contains quintets of atoms */
#define MAXOFNAME 5

/* special hydrogen bond terms come in pairs */
#define MAXHBONDNAME 2

/* array dimensions for the special h-bond coefficients (hbondcoeff) */
#define MINHBONDCOEFF 1
#define MAXHBONDCOEFF 3

/* maximum number of special adjustment terms to the mixrule */
#define MAXMIXADJUST 3

/* CLASSICAL POTENTIAL ALIASES
   Maps classical potential names to integers to speed
   up inner loops but retain readability.
   These ints can be compared against classical_pot_num
   which is set in readclassical.F
   Changing these settings is NOT advised  */

#define CP_INIT_ERROR 0
#define CP_LENNARD_JONES 1
#define CP_12_6_PLUS_SOLVATION 2
#define CP_UFF_12_6 3
#define CP_SCALED_LENNARD_JONES 4
#define CP_12_6_PLUS_12_10_H_BOND 5
#define CP_9_6 6
#define CP_12_9_6 7
#define CP_EXPONENTIAL_6 8
#define CP_HARD_SPHERE 9
#define CP_REPULSIVE_SPHERE 10
#define CP_EMBEDDED_ATOM_METHOD 111
#define CP_EXPONENTIAL_12_6 12
#define CP_GORDON_N_6 13
#define CP_MULTIWELL 114
#define CP_REPULSIVE_MULTIWELL 115
#define CP_REPULSIVE_WELL 16
#define CP_SQUARE_WELL 17
#define CP_STILLINGER_WEBER 18
#define CP_TABULATED_PAIR 119
#define CP_WEEKS_CHANDLER_ANDERSON 20
#define CP_HARD_2580_MULTISTEP 21
#define CP_REPULSIVE_2580_MULTISTEP 22
#define CP_EAM_PAIR_ONLY 123
#define CP_SW_PAIR_ONLY 24
#define CP_BUFFERED_14_7 25
#define CP_DOUBLE_EXPONENTIAL 26
#define CP_DRUKKER_H_BOND 27

/* eam pair integer conversions, arbitrary values */
#define EAM_PAIR_NONE 0
#define EAM_PAIR_ACKLANDTP 1
#define EAM_PAIR_ACKLANDPOW 2
#define EAM_PAIR_EXPON 3
#define EAM_PAIR_MORSE 4
#define EAM_PAIR_TABLE 5
#define EAM_PAIR_BELASHFE 6

/* eam density integer conversions, arbitrary values */
#define EAM_DENS_NONE 0
#define EAM_DENS_ACKLAND 1
#define EAM_DENS_EXPON 2
#define EAM_DENS_TABLE 3

/* eam embed integer conversion, arbitrary values */
#define EAM_EMBED_BELASH 1
#define EAM_EMBED_BELASHTRHO 2
#define EAM_EMBED_LOGARITHMIC 3
#define EAM_EMBED_POWHT 4
#define EAM_EMBED_POWHTF 5
#define EAM_EMBED_SR 6
#define EAM_EMBED_TABLE 7

/* solvation types */
#define SOLV_NONE 0
#define SOLV_SASA 1
#define SOLV_EEF1 2
#define SOLV_TRAMONTO 3

/*  parallel style definitions */
#define DEFAULT -2
#define NONE -1
#define TRAMONTO 0
#define JOBFARM 1
#define REX 2

#ifdef SAFE_COMPARE
/* define relative and absolute errors.  These are intended for */
/* double precision comparisons */
#define CMP_EPS_REL 1e-12
#define CMP_EPS_ABS 1e-12
#endif

/* constants: the parentheses make sure the full value is used and*/
/* not accidentally truncated at the end of a line */
/* Description of constants:                                      */

/* EETOKCAL converts coulombic potential qi*qj/rij to [kcal/mol], */
/*    where charges are multiples of electron charge e and        */
/*    distances are in angstroms.                                 */
/*       eetokcal = e^2 N_a / (4 pi \epsilon_0 10^-10 4184.0)     */
/* KCALTOK converts [kcal/mol] to [K / molecule];                 */
/*       kcaltok = 1 / (k_B * N_a)                                */
/*    where N_a is Avogadro constant, e is electron charge in [C],*/
/*    \epsilon_0 is permittivity in [C / V m], k_B is Boltzmann   */
/*    constant in [kcal/K]                                        */
/* PRESSCONST is k_B in [J/K] * 10^27; converts [K] to            */
/*    [kPa * angstrom^3].                                         */
#define CON_KCALTOK (503.2166d0)
#define CON_KJTOK (120.2717d0)
#define CON_EVTOK (11604.449d0)
#define CON_EETOKCAL (332.0636d0)
#define CON_ATMTOKPA (101.325d0)
#define CON_KPATOSIM (7.2429d-5)
#define CON_BOHRRADIUS (0.52918d0)
#define CON_RYDBERGTOEV (13.602d0)
#define CON_PRESSCONST (1.380662d4)
#define CON_HARTREETOKJ (2625.5d0)

/* global data structure variables */
/* the pound defines are used so the code is still readable, */
/*  but the use of integers improved efficiency */
/* The numerical values are arbitrary */
#define GLB_SET 2001
#define GLB_GET 2002
#define GLB_FIX 2003
#define GLB_CHK 2004
#define GLB_UST 2005
#define GLB_UGT 2006
#define GLB_INIT 2007
#define GLB_INCR 2008
#define GLB_DECR 2009
#define GLB_SCALE 2010
#define GLB_STORE 2011
#define GLB_EXCHANGE 2012
#define GLB_ALLOC 2013
#define GLB_FREE 2014
#define GLB_GET_TRIPLE 2015
#define GLB_SET_TRIPLE 2016
#define GLB_INCR_TRIPLE 2017
#define GLB_ALLOC_INIT 2018
#define GLB_STORE_ALL 2019
#define GLB_SCALE_ALL 2020
#define GLB_OUTPUT 2021
#define GLB_CLOSE 2022

/* coordinate integer flags */
/* values are arbitrary, but must be nonzero and different */
#define CRD_REAL 3001
#define CRD_TEMP 3002
#define CRD_NONE 3003
#define CRD_SINGLE 3004
#define CRD_R_TO_T 3005
#define CRD_T_TO_R 3006
#define CRD_R_TO_N 3007
#define CRD_N_TO_R 3008
#define CRD_R_TO_O 3009
#define CRD_O_TO_R 3010

/* Indices specifying new or old trials.  Used extensively for
   proposed moves and energy calculations.  Values NOT arbitrary */
#define TRL_OLD 1
#define TRL_NEW 2

/* Indices specifying real or trial box dimensions.  Used for volume
   changing moves.  Values NOT arbitrary */
#define BXD_T_TO_R -1
#define BXD_R_TO_T 0
#define BXD_REAL 1
#define BXD_TEMP 2
#define MAXBXD 2

/* Flags for rwconf */
/* First three are for the 'flag' argument to twh_rwconf */
/* Last two are for the 'rwflag' argument to twh_rwconf and other rw*
 * procedures.  Values arbitrary.  */
#define RWCONF_FINAL 0
#define RWCONF_BACKUP 1
#define RWCONF_RESTART 2
#define READ_FLAG 10
#define WRITE_FLAG 11

/* Values of ensemble -- (numerical values arbitrary) */
#define ENS_NVT 156
#define ENS_NPT 157
#define ENS_muVT 158
#define ENS_BUBBLE 159

/* values for pressurestyle.  Values arbitrary.*/
#define PRES_RADIAL 63116
#define PRES_VIRIAL 63110
#define PRES_THERMO 63118

/* Indices for scaling of hamiltonian. Values *not* arbitrary */
#define LAMBDA_LJ 1
#define LAMBDA_C 2
#define ALPHA_LJ 3

/* scaling style parameters.  Numerical values arbitrary. */
#define SCP_NONE 222
#define SCP_SC 223
#define SCP_WCA 224
#define SCP_NONE_STR 'none'
#define SCP_SC_STR 'soft-core'
#define SCP_WCA_STR 'WCA'

/* cbregrow move types.  Numerical values arbitrary */
#define CBR_ORIGINAL 333
#define CBR_BACKBONE 334
#define CBR_SIDECHAIN 335

/* cutstyle options */
#define CUTSTYLE_FULL 0

/* order integer values.  Values not arbitrary */
#define ORD_FORWARD 1
#define ORD_NULL 0
#define ORD_REVERSE -1

/* angle angle compute style.  Values arbitrary */
#define AA_COORD 0
#define AA_ANGLE 1

/* average flags: the numerical values are arbitrary */
#define AVG_INCREMENT 1 
#define AVG_EXTRACT_AVERAGE 2 
#define AVG_EXTRACT_BLOCK_VAL 3
#define AVG_SET_BLOCK_VAL 4
#define AVG_EXTRACT_BLOCK_AVG 5
#define AVG_EXTRACT_BLOCK_STD 6
#define AVG_EXTRACT_COUNT 7
#define AVG_ASSIGN 8
#define AVG_END_OF_BLOCK 9
#define AVG_EXTRACT_NBLOCK 10
#define AVG_INITIALIZE 11

/* maxener = the maximum number of different kinds of energy terms */
/* tracked in the simulation */
#define MAXENER 10
/* average codes: the numerical values are arbitrary, but must be */
/* larger than MAXENER, except for pressure */
/* pressure goes negative for the stress tensor and must be 0 */
#define AC_PRESSURE (0)
/* the energy accumulators are in spots 1 to MAXENER */
/* next are averages that only occur at most 1 time */
#define AC_H_VAP_DIRECT (MAXENER+1)
#define AC_H_VAP_VAPOR_P (MAXENER+2)
#define AC_H_VAP_PV_NRT (MAXENER+3)
/* averages that only depend on the number of boxes */
#define AC_TOT_ENERGY_SQ (MAXENER+4)
#define AC_DU_DLAMBDA_LJ (MAXENER+5)
#define AC_DU_DLAMBDA_C (MAXENER+6)
#define AC_VOLUME (MAXENER+7)
#define AC_VOLUME_SQ (MAXENER+8)
#define AC_ENTHALPY (MAXENER+9)
#define AC_ENTHALPY_SQ (MAXENER+10)
#define AC_PV (MAXENER+11)
#define AC_DU_DV (MAXENER+12)
#define AC_DU_DV_EXPAND (MAXENER+13)
#define AC_DU_DV_CONTRACT (MAXENER+14)
#define AC_THERMO_PRESSURE (MAXENER+15)
/* averages that only depend on the molecule type */
#define AC_HENRY_LAW (MAXENER+16)
/* averages that depend on the box and molecule type */
#define AC_NUM_MOLECULES (MAXENER+17)
#define AC_MOL_FRACTION (MAXENER+18)
#define AC_NUMBER_DENSITY (MAXENER+19)
#define AC_SPECIFIC_DENSITY (MAXENER+20)
#define AC_RADIUS_OF_GYRATION (MAXENER+21)
#define AC_GIBBS_TOT_CHEMPOT (MAXENER+22)
#define AC_NPT_INSERTION_CHEMPOT (MAXENER+23)
#define AC_NVT_INSERTION_CHEMPOT (MAXENER+24)
#define AC_ISOLATION_CHEMPOT (MAXENER+25)
/* MAXBAPROP = maximum number of properties to average considering */
/* the ones that depend on box and molecule type as multiple avaerage */
/* pressure is 8 terms */
#define MAXBAPROP 8+MAXENER*MAXBOX+3+12*MAXBOX+NTMAX+9*MAXBOX*NTMAX

/* field integers, numerical values arbitrary */
#define FLD_HARMONIC_ATTRACTOR 1
#define FLD_HARD_WALL 2
#define FLD_LJ_9_3_WALL 3
#define FLD_HOOPER_UMBRELLA 4
#define FLD_STEELE_WALL 5
#define FLD_HARD_BALL 6

/* center of mass call integers, numerical values arbitrary */
#define CTR_ROTATE 2
#define CTR_SWAP 3
#define CTR_TRANATOM 4
#define CTR_INITIAL 6
#define CTR_CBREGROW 7
#define CTR_ENERGY 9
#define CTR_INITCONF 10
#define CTR_CONROT 11
#define CTR_PIVOT 12

/* field storage maximums */
#define MAXFIELDINT 5
#define MAXFIELDDP 6
#define MAXFIELDNBPARAM 6

/* hard ball field integers, values not important except that they */
/* are less than the dimensions for the storage arrays */
#define FLD_HRBBOX 1

/* hard ball field double precision */
#define FLD_HRBPOSX 1
#define FLD_HRBPOSY 2
#define FLD_HRBPOSZ 3
#define FLD_HRBRAD 4
#define FLD_HRBFAC 5
#define FLD_HRBENERGY 6

/* cbmc generation codes, values arbitrary */
#define CBMC_AUTOFIT_GAUSSIAN 0
#define CBMC_GLOBAL_GAUSSIAN 1
#define CBMC_R_SQ_WITH_BOUNDS 2
#define CBMC_IDEAL_PLUS_AUTOFIT_GAUSSIAN 3
#define CBMC_IDEAL 4

/* cbmc distribution codes, values arbitrary */
#define DIST_NULL 0
#define DIST_DELTA 1
#define DIST_R_SQ 2
#define DIST_GAUSSIAN 3
#define DIST_SINE 4
#define DIST_SINE_GAUSSIAN 5
#define DIST_UNIFORM 6
#define DIST_UNIFORM_GAUSSIAN 7
#define DIST_ENERGY_BIAS 8
#define DIST_BOUNDED_SINE 9
#define DIST_TRIANGLE 10

/* cbmc formulation codes, values arbitrary */
#define CBMC_FORM_MS1999_MT2004 0
#define CBMC_FORM_MF2006 1
/* cbmc schedule codes, values arbitrary */
#define SCH_CODE_COUNT 0
#define SCH_CODE_GROW 1
#define SCH_CODE_SWAP 2
#define SCH_CODE_BACKBONE 4
#define SCH_CODE_LIST 5
#define SCH_CODE_SIDECHAIN 6

/* current forcefield version */
#define CURRENT_FFVERSION 15

/* bond style integer codes, do not change for compatibility */
#define BOND_STYLE_NULL -1
#define BOND_STYLE_FIXED 1
#define BOND_STYLE_HARMONIC 2
#define BOND_STYLE_GROMOS 3
#define BOND_STYLE_NONLINEAR 4
#define BOND_STYLE_MMTWO 5
#define BOND_STYLE_COMPASS 6
#define BOND_STYLE_NONBOND 7
#define BOND_STYLE_NONE 8
#define BOND_STYLE_MORSE 9
#define BOND_STYLE_SQUARE_WELL 10
#define BOND_STYLE_HARMONIC_NONBOND 11
#define BOND_STYLE_FENE 12

/* angle style integer codes, do not change for compatibility */
#define ANG_STYLE_NULL -1
#define ANG_STYLE_RIGID 0
#define ANG_STYLE_HARMONIC 1
#define ANG_STYLE_DREIDING 2
#define ANG_STYLE_HARM_COS 3
#define ANG_STYLE_COMPASS_AUTO 4
#define ANG_STYLE_CHARMM_UB 5
#define ANG_STYLE_NONBOND 6
#define ANG_STYLE_NB_HARM 7
#define ANG_STYLE_COMPASS_EX 8
#define ANG_STYLE_FOURIER 9
#define ANG_STYLE_TTFOURIER 10
#define ANG_STYLE_NONE 11
#define ANG_STYLE_SIXTH_AUTO 12
#define ANG_STYLE_SQUARE_WELL 13
#define ANG_STYLE_MULTI_RIGID 14
#define ANG_STYLE_MMFF 15
#define ANG_STYLE_HARM_COS_OT 16
#define ANG_STYLE_MMFFC_ANG 17
#define ANG_STYLE_MMFFC_BA 18

/* torsion style integer codes, do not change for compatibility */
#define TOR_STYLE_NULL -1
#define TOR_STYLE_SIM_HARM 1
#define TOR_STYLE_OPLS_COS 2
#define TOR_STYLE_GCA_COS 3
#define TOR_STYLE_GCA_CUB 4
#define TOR_STYLE_COMP_X_AD 5
#define TOR_STYLE_COMP_R_AD 6
#define TOR_STYLE_TRAPPE_COS 7
#define TOR_STYLE_NB 8
#define TOR_STYLE_COMP_X_EX 9
#define TOR_STYLE_COS_POW 10
#define TOR_STYLE_O_OPLS_C 11
#define TOR_STYLE_SUM_COS 12
#define TOR_STYLE_O_OPLS_T 13
#define TOR_STYLE_UFF 14
#define TOR_STYLE_DREIDING 15
#define TOR_STYLE_TF_COS 16
#define TOR_STYLE_TRAPPE_PLAN 17
#define TOR_STYLE_SQUARE_WELL 18
#define TOR_STYLE_AMBER 19
#define TOR_STYLE_OPLS_FL 20
#define TOR_STYLE_MRD 21
#define TOR_STYLE_FL 22
#define TOR_STYLE_DRUKKER 23

/* torsion scaling codes, value arbitrary */
#define CODE_ZERO 0
#define CODE_ONE 1
#define CODE_RESCALE 2

/* swapmove codes */
#define SWAP_CB_INTER 0
#define SWAP_CB_INTRA 1
#define SWAP_AVB_ONE 2
#define SWAP_AVB_TWO 3
#define SWAP_AVB_THREE 4
#define SWAP_CB_GC 5
#define SWAP_RB_INTER 6
#define SWAP_CHEMPOT_BOX 7
#define SWAP_CHEMPOT_ISO 8
#define SWAP_COM_INTER 9
/* interpolate codes */
#define INTERP_NORMAL 0
#define INTERP_DERIV 1

/* random number generator codes, arbitary but do not change */
#define RNG_RANLUX 0
#define RNG_DEBUG 1
#define RNG_DX_1597_2_7 5
#define RNG_KISS99 7
#define RNG_MRG32K3A 8

/* logical integers */
#define LOG_TRUE 1
#define LOG_FALSE 0

/* integer pointers */
#define PNT_GROWPREV 0
#define PNT_GROWFROM 1
#define PNT_GROWNUM 2
#define PNT_GROWVALIDCOUNT 3
#define PNT_GLIST 4
#define PNT_GLOBALPOS 5
#define PNT_LOGICAL_EXSCHED 6
#define PNT_LOGICAL_EXIST 7
#define PNT_ACNVOL 8
#define PNT_ACSVOL 9
#define PNT_EWALD_KMAX 10
#define PNT_LOGICAL_MOVEME 11
#define PNT_LOGICAL_PERIODIC 12
#define PNT_MOLTYP 13
#define PNT_NBOXI 14
#define PNT_PARALL 15
#define PNT_CHAINLIST 16
#define PNT_TOROFCODE 17
#define PNT_BAVOL 18
#define PNT_BNVOL 19
#define PNT_PAIRBOX 20
#define PNT_BNCELL 21
#define PNT_BACELL 22
#define PNT_ACNCELL 23
#define PNT_ACSCELL 24
/* max is one more than the final value */
#define MAX_INT_POINTERS 30

/* double pointers */
#define PNT_GYRATION 0
#define PNT_EAM_RHO_REAL 1
#define PNT_EAM_RHO_TEMP 2
#define PNT_RCMU 3
#define PNT_TMMC_WEIGHT 4
#define PNT_V_SEMIGRAND 5
#define PNT_CUBELETWEIGHT 6
#define PNT_COORDSTORAGE 7
#define PNT_COORDTEMP 8
#define PNT_COORDFIELD 9
#define PNT_ACNROT 10
#define PNT_ACNCOMP 11
#define PNT_ACNSWITCH 12
#define PNT_ACSSWITCH 13
#define PNT_ACNTRAA 14
#define PNT_ACSTRAA 15
#define PNT_ACNTRAC 16
#define PNT_ACSTRAC 17
#define PNT_BATRAA 18
#define PNT_BNTRAA 19
#define PNT_BATRAC 20
#define PNT_BNTRAC 21
#define PNT_RMTRAA 22
#define PNT_RMTRAC 23
#define PNT_PM2COMSWBOXPAIR 24
#define PNT_BNCOM2BSWITCH 25
/* open spot */
/* open spot */
#define PNT_ACSROT 28
#define PNT_ACSCOMP 29
#define PNT_RMCOMROT 30
#define PNT_RMCOMTRA 31
#define PNT_RMROT 32
#define PNT_RMVOL 33
#define PNT_ARBCOMFIELD 34
#define PNT_BAROT 35
#define PNT_BNROT 36
#define PNT_CMATRIX 37
#define PNT_COMFIELD 38
#define PNT_COMTEMPFIELD 39
#define PNT_WRAP_FOREIGN_ENERGY 40
#define PNT_WRAP_FOREIGN_LAMBDA_LJ 41
#define PNT_WRAP_FOREIGN_LAMBDA_C 42
#define PNT_BLOCKVALUE 43
#define PNT_RMCELL 44
#define PNT_PM2CBSWPR 45
#define PNT_PM2RBSWPR 46
#define PNT_PMCELLPR 47
#define PNT_PMCELLPT 48
#define PNT_PMVLPR 49
/* max must be at least one more than the final pointer */
#define MAX_DOUBLE_POINTERS 50

/* report integer, values arbitrary, but should be nonzero */
#define REP_BON_REG 1
#define REP_BON_VDW 2
#define REP_BON_CHG 3
#define REP_ANG_REG 4
#define REP_ANG_BA 5
#define REP_ANG_VDW 6
#define REP_ANG_CHG 7
#define REP_TOR_REG 8
#define REP_TOR_AT 9
#define REP_TOR_BT 10
#define REP_TOR_BB 11
#define REP_TOR_AAT 12
#define REP_TOR_VDW 13
#define REP_TOR_CHG 14

#define REP_MAX 14

/* electrostatic form integer code, values arbitrary */
#define EF_CODE_NONE 0
#define EF_CODE_SASA 1
#define EF_CODE_EFF 2
#define EF_CODE_BUF 3
#define EF_CODE_COUL_NOCUT 4
#define EF_CODE_COUL_EWALD 5

/* thermodyanmic pressure update integer code, values arbitrary */
#define PTHERMO_MOVE_DATA 1
#define PTHERMO_RANGE 2

/* fortran file numbers for input and output */
#define FILE_MOVIE 2
#define FILE_CHEMPOTDATA 3

/* distribution limit style code, values arbitrary */
#define LIM_VALUE 0
#define LIM_NEGINF 1
#define LIM_POSINF 2

/* volume move codes, values arbitrary */
#define VMC_GHOST 0
#define VMC_THREED 1
#define VMC_ONED 2
#define VMC_NVT_ISO 100
#define VMC_NVT_ONED 101
#define VMC_NPT_ISO 102
#define VMC_NPT_ONED 103
#define VMC_TPRESS 200

/* recip sum codes, values arbitrary */
#define RCP_ACCEPT_NEW_VOL 1
#define RCP_ACCEPT_NEW 2
#define RCP_STORE 3
#define RCP_RESTORE 4

/* volumemove codes */
#define BCODE_SINGLE 1
#define BCODE_EXCHANGE 2
#define BCODE_GHOST 3

/* the end */
