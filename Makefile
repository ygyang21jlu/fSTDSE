goal: makefile.dep
	  make spherical_tdse.x
	  make build_pes.x

MAKEFLAGS = -r

.SUFFIXES: .f90 .o .x .c .dep

#
#  This is the default; a config file may override it.
#
#ACT = sed -e 's/^!\*qd/    /' # Enable quad-math statements
ACT = sed -e 's/^!\*nq/    /' # Disable quad-math statements

#
# System-specific overrides
#
#  include vanilla.mak
# include configs/shelob-ifort_opt.mak
 include configs/shelob-ifort_noquad_opt.mak
# include configs/shelob-gfortran_opt.mak
# include configs/shelob-gfortran_dbg.mak
# include configs/shelob-gfortran_noquad_opt.mak
# include configs/shelob-gfortran_noquad_dbg.mak
# include configs/mic-ifort_opt.mak

#
# Finish the set-up
#
LIBS = $(LAPACK) $(LAPACK) $(LIBEXTRA)

#
# Compiling and archiving rules
#
.f90.o:
	$(ACT) $< >preprocess/$<
	$(F90) -c preprocess/$<

dgefa.o:        dgefa.f
	$(F90) -c dgefa.f

dgedi.o:        dgedi.f
	$(F90) -c dgedi.f

clean:
	-/bin/rm -f *.{o,mod,x,il,a} *__genmod.f90 checkpoint_{field,main}.* makefile.dep *.optrpt ./preprocess/*.f90

makefile.dep: $(shell echo *.f90)
	./make-depend.sh $^ > $@

#
# Explicit dependencies
#

LIBSPHERICAL += accuracy.o
LIBSPHERICAL += bicg_tools.o
LIBSPHERICAL += cap_tools.o
LIBSPHERICAL += composition_analysis.o
LIBSPHERICAL += constants.o
LIBSPHERICAL += dgefa.o
LIBSPHERICAL += dgedi.o
LIBSPHERICAL += lapack.o
LIBSPHERICAL += math.o
LIBSPHERICAL += potential_tools.o
LIBSPHERICAL += propagator_tools.o
LIBSPHERICAL += rotation_tools.o
LIBSPHERICAL += sort_tools.o
LIBSPHERICAL += spherical_data.o
LIBSPHERICAL += test_tools.o
LIBSPHERICAL += timer.o
LIBSPHERICAL += tridiagonal_tools.o
LIBSPHERICAL += vectorpotential_tools.o
LIBSPHERICAL += wavefunction_tools.o

#
# Building the binaries
#
spherical_tdse.x: spherical_tdse.o $(LIBSPHERICAL)
	$(F90) -o spherical_tdse.x spherical_tdse.o $(LIBSPHERICAL) $(LIBS)
	-hugeedit --text --data spherical_tdse.x

build_pes.x: build_pes.o $(LIBSPHERICAL)
	$(F90) -o build_pes.x build_pes.o $(LIBSPHERICAL) $(LIBS)
#
# Automatically-generated dependencies
#
include makefile.dep
