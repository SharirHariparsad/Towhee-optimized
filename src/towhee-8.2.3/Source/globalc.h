/* 
* MCCCS - Towhee: A Monte Carlo molecular simulation program     *
* Copyright (C) 2003-2018 Marcus G. Martin                       *
* see the file license.gpl for the full license information      *
*                                                                *
* This program is free software; you can redistribute it and/or  *
* modify it under the terms of the GNU General Public License    *
* as published by the Free Software Foundation; either version 2 *
* of the License, or (at your option) any later version.         *
*                                                                *
* This program is distributed in the hope that it will be useful,*
* but WITHOUT ANY WARRANTY; without even the implied warranty of *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
* GNU General Public License for more details.                   *
*                                                                *
* You should have received a copy of the GNU General Public      *
* License along with this program; if not, write to the Free     *
* Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,*
* MA  02111-1307, USA.                                           
*/

/*
  globalc.h
  this file contains all of the information that is at the top 
  of every C routine in Towhee
  last modified 04-03-2018 by M.G. Martin
   generic c includes 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>

/* Include general towhee definitions */
#include "preproc.h"

/* C specific preprocessor derectives */
/* maximum number of timers for keeping track of mpi_wtime */
#define MAXTIMERS 5

#define TRUE 1
#define FALSE 0

/* Structures */
/*
 * user_args contains information about command line arguments passed to
 * towhee.  Populated by parse_args, and passed around as a pointer to be
 * processed where needed.
 */
typedef struct {
  char rex_fn[MAXDIRLENGTH];
  char towhee_input_fn[MAXDIRLENGTH];
  int quiet_mode;
  int parstyle;
  int random_seed;
} user_args;

typedef struct {
  int num_nodes;
  int num_steps;
  int num_rounds;
  int num_swaps;
  float lambda_lj[MAX_FOREIGN_LAMBDA];
  float lambda_c[MAX_FOREIGN_LAMBDA];
} rex_params;

/* 
   template or function declarations 
*/
/* extern void calltest_(); */

#ifdef INTEL_VISUAL_FORTRAN
#define towheemainloop_ TOWHEEMAINLOOP
#define swrap_ SWRAP
#endif

#ifdef MSVC
#define snprintf _snprintf
#endif

/* external function declarations */
extern void towheemainloop_(const int *,const int *, const int *, const int *);
/* internal function declarations */
void ** twh_allocate2dMatrix(size_t, int, int);
void * twh_allocateVector(size_t, int);
int attempt_exchange(float, int, int, double**, int*, rex_params*, FILE*);
void clearline (FILE *);
void do_jobfarm(user_args*, int);
void do_none(user_args*);
void do_rex(user_args*, int);
void do_tramonto(user_args*, int);
int get_parstyle(int, int);
int get_rex_directory(char*, const char*, int, int, int);
void initialize_args( user_args*);
int initialize_MPI(int, char**);
void initialize_wrapper_strings();
void jobfarm(FILE *, int, int);
void open_towhee_parallel(const char*, FILE **);
void parse_args(int, char**, user_args*);
int populate_directory_initial(char*, const char*, int);
int populate_directory_restart(char*, const char*, int, int);
void print_usage(const char*);
int read_rex_params(const char*, rex_params*);
void run_rex_client(user_args*, rex_params*, int);
void run_rex_server(user_args*, rex_params*);
void set_output_dir(const char*);
void set_towhee_input_file(const char*);
void twh_time_array(int, int, double*);
void timekeeper(int, int, int);
void timestamp(char*);
double tramonto_();
void tramonto_control(int, int);
void tramontoloop();
int write_lambda_file(const int*, const rex_params*, const char*, int, const char*);
/* 
   USEMPI only  
*/
#ifdef USEMPI
#include <mpi.h>
#endif
/* 
   USETRAMONTO only 
*/
#ifdef USETRAMONTO
extern dftmain(double *);
#endif

/* templates internal to c */
int twh_1ddouble(const int, int *, const int, double *);
int twh_1dinteger(const int, int *, const int, int *);
int twh_2dregdouble(const int, int *, const int, const int, double *);
int twh_2dreginteger(const int, int *, const int, const int, int *);
int twh_3dreginteger(const int, int *, const int, const int, const int, int *);
int twh_3dregdouble(const int, int *, const int, const int, const int, double *);
/* data storage */

/* templates callable from Fortran */
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acncell_ TWH_ACNCELL
#endif
void twh_acncell_(int *, int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acncomp_ TWH_ACNCOMP
#endif
void twh_acncomp_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acnrot_ TWH_ACNROT
#endif
void twh_acnrot_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acnswitch_ TWH_ACNSWITCH
#endif
void twh_acnswitch_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acntraa_ TWH_ACNTRAA
#endif
void twh_acntraa_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acntrac_ TWH_ACNTRAC
#endif
void twh_acntrac_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acnvol_ TWH_ACNVOL
#endif
void twh_acnvol_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acscell_ TWH_ACSCELL
#endif
void twh_acscell_(int *, int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acscomp_ TWH_ACSCOMP
#endif
void twh_acscomp_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acsrot_ TWH_ACSROT
#endif
void twh_acsrot_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acsswitch_ TWH_ACSSWITCH
#endif
void twh_acsswitch_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acstraa_ TWH_ACSTRAA
#endif
void twh_acstraa_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acstrac_ TWH_ACSTRAC
#endif
void twh_acstrac_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acsvol_ TWH_ACSVOL
#endif
void twh_acsvol_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_arbcmofield_ TWH_ARBCOMFIELD
#endif
void twh_arbcmofield_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_bacell_ TWH_BACELL
#endif
void twh_bacell_(int *, int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_barot_ TWH_BAROT
#endif
void twh_barot_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_batraa_ TWH_BATRAA
#endif
void twh_batraa_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_batrac_ TWH_BATRAC
#endif
void twh_batrac_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_bavol_ TWH_BAVOL
#endif
void twh_bavol_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_bncell_ TWH_BNCELL
#endif
void twh_bncell_(int *, int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_bnrot_ TWH_BNROT
#endif
void twh_bnrot_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_bntraa_ TWH_BNTRAA
#endif
void twh_bntraa_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_bntrac_ TWH_BNTRAC
#endif
void twh_bntrac_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_bnvol_ TWH_BNVOL
#endif
void twh_bnvol_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_c_matrix_ TWH_C_MATRIX
#endif
void twh_c_matrix_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_chainlist_ TWH_CHAINLIST
#endif
void twh_chainlist_(int *, int *, int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_comfield_ TWH_COMFIELD
#endif
void twh_comfield_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_comtempfield_ TWH_COMTEMPFIELD
#endif
void twh_comtempfield_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_coordfield_ TWH_COORDFIELD
#endif
void twh_coordfield_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_coordstorage_ TWH_COORDSTORAGE
#endif
void twh_coordstorage_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_coordtemp_ TWH_COORDTEMP
#endif
void twh_coordtemp_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_cubeletweight_ TWH_CUBELETWEIGHT
#endif
void twh_cubeletweight_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_eam_rho_real_ TWH_EAM_RHO_REAL
#endif
void twh_eam_rho_real_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_eam_rho_temp_ TWH_EAM_RHO_TEMP
#endif
void twh_eam_rho_temp_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_ewald_kmax_ TWH_EWALD_KMAX
#endif
void twh_ewald_kmax_(int *, int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_glist_ TWH_GLIST
#endif
void twh_glist_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_globalpos_ TWH_GLOBALPOS
#endif
void twh_globalpos_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_growfrom_ TWH_GROWFROM 
#endif
void twh_growfrom_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_grownum_ TWH_GROWNUM
#endif
void twh_grownum_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_growprev_ TWH_GROWPREV
#endif
void twh_growprev_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_growvalidcount_ TWH_GROWVALIDCOUNT
#endif
void twh_growvalidcount_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_gyration_ TWH_GYRATION
#endif
void twh_gyration_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_logical_exist_ TWH_LOGICAL_EXIST
#endif
void twh_logical_exist_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_logical_exsched_ TWH_LOGICAL_EXSCHED
#endif
void twh_logical_exsched_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_logical_moveme_ TWH_LOGICAL_MOVEME
#endif
void twh_logical_moveme_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_logical_periodic_ TWH_LOGICAL_PERIODIC
#endif
void twh_logical_periodic_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_moltyp_ TWH_MOLTYP
#endif
void twh_moltyp_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_nboxi_ TWH_NBOXI
#endif
void twh_nboxi_(int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_pairbox_ TWH_PAIRBOX
#endif
void twh_pairbox_(int*, int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_parall_ TWH_PARALL
#endif
void twh_parall_(int *, int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_pm2cbswpr_ TWH_PM2CBSWPR
#endif
void twh_pm2cbswpr_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_pm2rbswpr_ TWH_PM2RBSWPR
#endif
void twh_pm2comswboxpair_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_pm2comswboxpair_ TWH_PM2COMSWBOXPAIR
#endif
void twh_pm2rbswpr_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_pmcellpr_ TWH_PMCELLPR
#endif
void twh_pmcellpr_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_pmcellpt_ TWH_PMCELLPT
#endif
void twh_pmvlpr_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_pmvlpr_ TWH_PMVLPR
#endif
void twh_pmcellpt_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_rcmu_ TWH_RCMU
#endif
void twh_rcmu_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_rmcell_ TWH_RMCELL
#endif
void twh_rmcell(int *, int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_rmcomrot_ TWH_RMCOMROT
#endif
void twh_rmcomrot_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_rmcomtra_ TWH_RMCOMTRA
#endif
void twh_rmcomtra_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_rmrot_ TWH_RMROT
#endif
void twh_rmrot_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#endif
void twh_rmtraa_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_rmtraa_ TWH_RMTRAA
#define twh_rmtrac_ TWH_RMTRAC
#endif
void twh_rmtrac_(int *, int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_rmvol_ TWH_RMVOL
#endif
void twh_rmvol_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_torofcode_ TWH_TOROFCODE
#endif
void twh_torofcode_(int *, int *, int *, int *, int *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_tmmc_weight_ TWH_TMMC_WEIGHT
#endif
void twh_tmmc_weight_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_v_semigrand_ TWH_V_SEMIGRAND
#endif
void twh_v_semigrand_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_wrap_foreign_energy_ TWH_WRAP_FOREIGN_ENERGY
#endif
void twh_wrap_foreign_energy_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_wrap_foreign_lambda_lj_ TWH_WRAP_FOREIGN_LAMBDA_LJ
#endif
void twh_wrap_foreign_lambda_lj_(int *, int *, double *);
#ifdef INTEL_VISUAL_FORTRAN
#define twh_wrap_foreign_lambda_c_ TWH_WRAP_FOREIGN_LAMBDA_C
#endif
void twh_wrap_foreign_lambda_c_(int *, int *, double *);
