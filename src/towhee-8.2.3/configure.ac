#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.59)
# Project name, version, and POC
AC_INIT(towhee,8.2.3,towhee-bugs@lists.sourceforge.net)
AC_CONFIG_SRCDIR(Source)
AC_CONFIG_AUX_DIR(install)
AM_INIT_AUTOMAKE([subdir-objects])

# Do not require the Makefiles to be rebuilt using the autotools
# so the end user will not have to have the same version of autoconf
# and automake.
AM_MAINTAINER_MODE

# Initialize stuff for automake
AM_INIT_AUTOMAKE

# Define the argument to configure script to use internal erfs
AC_ARG_ENABLE(internal,
     [  --enable-interal   Compile version using internal erfs],
     [case "${enableval}" in
       yes) UseInternal=true ;;
       no)  UseInternal=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-internal) ;;
     esac;echo "UseInternal is $UseInternal"],[UseInternal=false])

# Define the argument to configure script to enable tramonto
AC_ARG_ENABLE(tramonto,
     [  --enable-tramonto   Compile version using tramonto libraries],
     [case "${enableval}" in
       yes) UseTramonto=true ;;
       no)  UseTramonto=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-tramonto) ;;
     esac;echo "UseTramonto is $UseTramonto"],[UseTramonto=false])

# Define the argument to enable Safe Compare
 AC_ARG_ENABLE(safe_compare,
      [  --enable-safe-compare    Emphasize cross-platform consistency 
           over performance.],
      [case "${enableval}" in
        yes) UseSafeCompare=true ;;
        no)  UseSafeCompare=false ;;
        *) AC_MSG_ERROR(bad value ${enableval} for --enable-safe-compare) ;;
      esac;echo "UseSafeCompare is $UseSafeCompare"],[UseSafeCompare=false])

# Define the argument to enable fix GNU
 AC_ARG_ENABLE(fix_GNU,
      [  --enable-fix-GNU   Fixes the problem with double underscore.],
      [case "${enableval}" in
        yes) UseFixGNU=true ;;
        no)  UseFixGNU=false ;;
        *) AC_MSG_ERROR(bad value ${enableval} for --enable-fix-GNU) ;;
      esac;echo "UseFixGNU is $UseFixGNU"],[UseFixGNU=false])

# Define the argument to configure script to enable MPI
AC_ARG_ENABLE(mpi,
     [  --enable-mpi    Compile parallel version using MPI],
     [case "${enableval}" in
       yes) UseMPI=true ;;
       no)  UseMPI=false ;;
       *) AC_MSG_ERROR(bad value ${enableval} for --enable-mpi) ;;
     esac;echo "UseMPI is $UseMPI"],[UseMPI=false])

# Turn off Command Line Arguments, if necessary
AC_ARG_ENABLE(command-line-args,
    [  --disable-command-line-args   Disable parsing of command line arguments by towhee executable],
    [case "${enableval}" in
      yes) DisableCLA=false ;;
      no)  DisableCLA=true ;;
      *) AC_MSG_ERROR(bad value ${enableval} for --disable-command-line-args) ;;
    esac;echo "DisableCLA is $DisableCLA"],[DisableCLA=false])

# Tell automake about the conditionals we just defined
AM_CONDITIONAL(UseMPI, test x$UseMPI = xtrue)
AM_CONDITIONAL(UseInternal, test x$UseInternal = xtrue)
AM_CONDITIONAL(UseTramonto, test x$UseTramonto = xtrue)
AM_CONDITIONAL(UseSafeCompare, test x$UseSafeCompare = xtrue)
AM_CONDITIONAL(UseFixGNU, test x$UseFixGNU = xtrue)
AM_CONDITIONAL(DisableCLA, test x$DisableCLA = xtrue)

# Checks for compilers
# It seems really silly to check for cc and when we know
#   where it is, and will set it, but otherwise configure
#   breaks.
# This has to be done oustide a conditional.  Silly, huh?
AC_PROG_CC
# fix for redhat icc
if test x${CC} = xicc ; then	
   CFLAGS="${CFLAGS} -pie -fpic"
fi

# Compiler Tests
# Already did CC test
AC_PROG_F77

# Link normally
AC_SUBST(F77LD,${CC})

#incorporate ADDLIB into our library
LIBS="${LIBS} ${ADDLIB}"

# extra items for UseTramonto
if test x$UseTramonto = xtrue ; then
  # Include the Tramonto libraries
  LIBS="${LIBS} -L../lib -ldft -laztecoo -ltriutils -lifpack \
  -lepetraext -lepetra -llapack -lblas"
fi
# extra items for UseInternal
if test x$UseInternal = xtrue ; then
  # Nothing happens here
  LIBS="${LIBS}"
fi
# extra items for UseSafeCompare
if test x$UseSafeCompare = xtrue ; then
  # Nothing happens here
  LIBS="${LIBS}"
fi
# extra items for UseFixGNU
if test x$UseFixGNU = xtrue ; then
  # Nothing happens here
  LIBS="$(LIBS)"
fi
# If enable-mpi, check for MPI compile scripts
if test x${UseMPI} = xtrue ;then
  # Using default of C language
  ACX_MPI
  # Set default language to Fortran 77, to make ACX_MPI look for mpif77
  AC_LANG_PUSH(Fortran 77)
  ACX_MPI
  AC_LANG_POP(Fortran 77)

  # Use the MPI compilers and libraries
  CC="${MPICC}"
  F77="${MPIF77}"
  LIBS="${MPILIBS} ${LIBS}"
  FFLAGS="${FFLAGS}"
  CFLAGS="${CFLAGS}"
fi

# Checks for libraries.
# Find the F77 libraries needed for linking
ACX_F77_LIBRARY_LDFLAGS
# Then, link using cc, not f77
AC_SUBST(F77LD,${CC})


# Checks for header files.
# Only needed if using MPI (for control.c)
if test x${UseMPI} = xtrue ;then
  AC_HEADER_STDC
  AC_CHECK_HEADERS([stdlib.h string.h])
fi

# Checks for typedefs, structures, and compiler characteristics.
# Checks for library functions.
AC_CHECK_FUNCS([strstr])

# Generate the config
# Have to list all the Makefiles
AC_CONFIG_FILES([Makefile \
Source/Makefile \
ForceFields/Makefile \
Manual/Makefile \
Manual/projects/Makefile \
Manual/projects/towhee/Makefile \
Manual/projects/towhee/algorithm/Makefile \
Manual/projects/towhee/code/Makefile \
Manual/projects/towhee/examples/Makefile \
Manual/projects/towhee/forcefields/Makefile \
Manual/projects/towhee/inpstyle/Makefile \
Manual/projects/towhee/input/Makefile \
Manual/projects/towhee/potentype/Makefile \
Manual/projects/towhee/utils/Makefile \
Utils/Makefile \
Examples/Makefile \
Examples/DFT_Field/Makefile \
Examples/Canonical_Ensemble/Makefile \
Examples/Canonical_Ensemble/AVB1_Methane/Makefile \
Examples/Canonical_Ensemble/Amber_Villin/Makefile \
Examples/Canonical_Ensemble/Benzene_AA_autofit/Makefile \
Examples/Canonical_Ensemble/Charmm19_ubiquitin/Makefile \
Examples/Canonical_Ensemble/Charmm27_Benzene/Makefile \
Examples/Canonical_Ensemble/Charmm27_Heme/Makefile \
Examples/Canonical_Ensemble/Charmm27_Nucleic_Acid/Makefile \
Examples/Canonical_Ensemble/Charmm27_Polyalanine/Makefile \
Examples/Canonical_Ensemble/Compass_Methanol/Makefile \
Examples/Canonical_Ensemble/DREIDING/Makefile \
Examples/Canonical_Ensemble/Formamide_Scaled/Makefile \
Examples/Canonical_Ensemble/Fris_Walls/Makefile \
Examples/Canonical_Ensemble/Hard_Sphere/Makefile \
Examples/Canonical_Ensemble/NaCl_1x1x1/Makefile \
Examples/Canonical_Ensemble/Shukla_Gasses/Makefile \
Examples/Canonical_Ensemble/Small_Peptide/Makefile \
Examples/Canonical_Ensemble/Square_Well_Chain/Makefile \
Examples/Canonical_Ensemble/TraPPE_Isomers/Makefile \
Examples/Canonical_Ensemble/TraPPE_Molecules/Makefile \
Examples/Canonical_Ensemble/Triglycerol/Makefile \
Examples/Canonical_Ensemble/UFF/Makefile \
Examples/Canonical_Ensemble/Wall_Water/Makefile \
Examples/Convert/Makefile \
Examples/Convert/LAMMPS_class2/Makefile \
Examples/Convert/LAMMPS_decane/Makefile \
Examples/Convert/LAMMPS_lc/Makefile \
Examples/Histogram/Makefile \
Examples/Histogram/Phase/Makefile \
Examples/Histogram/PVT/Makefile \
Examples/Histogram/Weights/Makefile \
Examples/Isobaric_Isothermal_Ensemble/Makefile \
Examples/Isobaric_Isothermal_Ensemble/Au_Cu_Switch/Makefile \
Examples/Isobaric_Isothermal_Ensemble/Charmm22_Ethanethiol/Makefile \
Examples/Isobaric_Isothermal_Ensemble/Cu_Pb_EAM/Makefile \
Examples/Isobaric_Isothermal_Ensemble/Dick1994_PETN/Makefile \
Examples/Isobaric_Isothermal_Ensemble/Gromos_Methylpropylsulfide/Makefile \
Examples/Isobaric_Isothermal_Ensemble/Henry_Law/Makefile \
Examples/Isobaric_Isothermal_Ensemble/Ideal_Chain/Makefile \
Examples/Isobaric_Isothermal_Ensemble/MMFF94_Ethyleneoxide/Makefile \
Examples/Isobaric_Isothermal_Ensemble/OPLS_Propanamide/Makefile \
Examples/Isobaric_Isothermal_Ensemble/SMMKmain_2244688nonane/Makefile \
Examples/Isobaric_Isothermal_Ensemble/Solid_LJium/Makefile \
Examples/Isobaric_Isothermal_Ensemble/Wielop1985_Ethyleneoxide/Makefile \
Examples/Gibbs_Ensemble/Makefile \
Examples/Gibbs_Ensemble/Amber_IsoPropanol/Makefile \
Examples/Gibbs_Ensemble/Catlow_Zeolite_4a/Makefile \
Examples/Gibbs_Ensemble/Cu_VLE/Makefile \
Examples/Gibbs_Ensemble/Cui2002_Hexane/Makefile \
Examples/Gibbs_Ensemble/Dubb_Zeolite/Makefile \
Examples/Gibbs_Ensemble/EPM_VLCC/Makefile \
Examples/Gibbs_Ensemble/FENE_Hexamer/Makefile \
Examples/Gibbs_Ensemble/Gordon/Makefile \
Examples/Gibbs_Ensemble/Gromos_Isobutane/Makefile \
Examples/Gibbs_Ensemble/MM2_Ethane/Makefile \
Examples/Gibbs_Ensemble/Potter_CF2H2/Makefile \
Examples/Gibbs_Ensemble/SKS_Pentane/Makefile \
Examples/Gibbs_Ensemble/SMMKnaip_Ethylpentane/Makefile \
Examples/Gibbs_Ensemble/TraPPE_Pentane/Makefile \
Examples/Gibbs_Ensemble/Vink2001_Silicon/Makefile \
Examples/Gibbs_Ensemble/Weiner1984/Makefile \
Examples/Grand_Canonical_Ensemble/Makefile \
Examples/Grand_Canonical_Ensemble/Amber96_Ethane/Makefile \
Examples/Grand_Canonical_Ensemble/Energy_Biasing/Makefile \
Examples/Grand_Canonical_Ensemble/Steele_Wall/Makefile \
Examples/Grand_Canonical_Ensemble/Walt2001_Nanotube/Makefile \
Examples/Parallel_Test/Makefile \
Examples/TMMC/Makefile \
Examples/TMMC/LJ/Makefile \
Examples/TMMC/SPC_E_WATER/Makefile \
Examples/VLCC_Fit/Makefile \
])
AC_OUTPUT
