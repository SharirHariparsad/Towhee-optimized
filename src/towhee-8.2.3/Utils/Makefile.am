#     $Id: Makefile.am,v 1.13 2009/06/17 15:56:58 marcus_martin Exp $
# Compiles various utlity programs
bin_PROGRAMS = \
		analyse_histogram \
		analyse_movie \
		car2towhee \
		charmm2pdb \
		faux2towhee \
		fitcoex \
		jre_to_towhee \
		maftodensity \
		pdb2towhee \
		rdf2pmfpair \
		unitcell \
		xmd2towhee \
		xtl2towhee \
		xyz2towhee

AM_FFLAGS = 

if UseFixGNU
FG_FFLAGS = -fno-second-underscore
endif

AM_CPPFLAGS = $(MPI_CPPFLAGS) $(T_CPPFLAGS) $(LCAO_CPPFLAGS) $(I_CPPFLAGS)
AM_CPPFLAGS += $(SC_CPPFLAGS) $(CLA_CPPFLAGS) $(DM_CPPFLAGS)
AM_FFLAGS += $(SC_FFLAGS) $(FG_FFLAGS) $(DM_FFLAGS)
LIBS = $(MPI_LIBS)

# Add the Fortran libraries, so we can link with C
LDADD = @FLIBS@

analyse_histogram_SOURCES = analyse_histogram.F
analyse_movie_SOURCES = analyse_movie.F
car2towhee_SOURCES = car2towhee.c
charmm2pdb_SOURCES = charmm2pdb.F
faux2towhee_SOURCES = faux2towhee.F
fitcoex_SOURCES = fitvlcc.F leastsquares.F
jre_to_towhee_SOURCES = jre_to_towhee.F
maftodensity_SOURCES = maftodensity.F
pdb2towhee_SOURCES = pdb2towhee.F
rdf2pmfpair_SOURCES = rdf2pmfpair.F
unitcell_SOURCES = unitcell.F
xmd2towhee_SOURCES = xmd2towhee.F
xtl2towhee_SOURCES = xtl2towhee.F
xyz2towhee_SOURCES = xyz2towhee.F

###Fortran dependencies
analyse_histogram.o: analyse_histogram.F
analyse_movie.o: analyse_movie.F
car2towhee.o: car2towhee.c
charmm2pdb.o: charmm2pdb.F
faux2towhee.o: faux2towhee.F
fitvlcc.o: fitvlcc.F leastsquares.F
jre_to_towhee.o: jre_to_towhee.F
leastsquares.o: leastsquares.F
maftodensity.o: maftodensity.F
pdb2towhee.o: pdb2towhee.F
rdf2pmfpair.o: rdf2pmfpair.F
unitcell.o: unitcell.F
xmd2towhee.o: xmd2towhee.F
xtl2towhee.o: xtl2towhee.F
xyz2towhee.o: xyz2towhee.F

EXTRA_DIST = \
parse_vlcc_plots.script
