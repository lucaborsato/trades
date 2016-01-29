# MY MAKEFILE
#
# gfortran compiler
# CC=gfortran
# CMPI=mpif90
# needed because some 'mpi_utils' has `use mpi` w/o #ifdef MPI #endif
CC = mpif90

# compiler options
CFLAGS=-cpp -Og -g -Wall -fcheck=all -ffree-line-length-none -Wsurprising -Waliasing -Wimplicit-interface -fbacktrace -ffpe-trap=invalid,zero,overflow
#CFLAGS+=-fmax-errors=1
#CFLAGS+=-Warray-temporaries
#CFLAGS+=-Wcharacter-truncation
#CFLAGS+=-Wconversion-extra
#CFLAGS+=-Wextra
# CFLAGS+=-pg # profiling 
COPT=-cpp -O3 -ffree-line-length-none -ffpe-trap=invalid,zero,overflow
# COPT+=-pg # profiling

# with openmp
COMP=-fopenmp

CMPI= -DMPI

LIBCHORD = libchord.a
AR = ar r
CCHORD = -L. -lchord

# executable
TARGET_SER=../bin/trades_s
TARGET_OMP=../bin/trades_o
TARGET_MPIOMP=../bin/trades_mo

# this will be used by everything, even PolyChord modules.
CONSTSRC = constants.f90
CONSTOBJ = constants.o

POLYCHORD = PolyChord/

#  OLD POLYCHORD V1.0 src
# PCSRC = $(addprefix $(POLYCHORD), utils.f90 mpi_utils.F90 random_utils.F90 priors.f90 grades.F90 settings.f90 calculate.f90 evidence.f90 feedback.f90 read_write.f90 clustering.f90 chordal_sampling.f90 generate.F90 nested_sampling.F90)
# PCOBJ = utils.o mpi_utils.o random_utils.o priors.o grades.o settings.o calculate.o evidence.o feedback.o read_write.o clustering.o chordal_sampling.o generate.o nested_sampling.o

PCSRC = utils.f90 abort.F90 settings.f90 params.f90 array_utils.f90 priors.f90 mpi_utils.F90 calculate.f90 random_utils.F90 chordal_sampling.f90 run_time_info.f90 clustering.f90 read_write.f90 feedback.f90 generate.F90 ini.f90 nested_sampling.F90
PCOBJ_1 = $(PCSRC:.f90=.o)
PCOBJ = $(PCOBJ_1.F90=.o)
# FILES_OUT = $(patsubst %.xls,%.xlsx,$(patsubst %.doc,%.docx,$(FILES_IN)))
PCOBJ = $(patsubst %.F90,%.o,$(patsubst %.f90,%.o, $(PCSRC)))

SOURCES = parameters.f90 random_trades.f90 convert_type.f90 lin_fit.f90 celestial_mechanics.f90  init_trades.f90 statistics.f90 timing.f90 rotations.f90 sort.f90 eq_motion.f90 output_files.f90 numerical_integrator.f90 radial_velocities.f90 transits.f90 ode_run.f90 derived_parameters_mod.f90 grid_search.f90 lm.f90 pikaia.f90  util_sort.f90 util_qmc.f90 opti_pso.f90 gaussian.f90 bootstrap.f90 PolyChord_driver.f90
OBJECTS = $(SOURCES:.f90=.o)

# .PHONY: print_obj
# print_obj:
# 	echo " "
# 	echo $(PCSRC)
# 	echo " "
# 	echo $(PCOBJ)
# 	echo " "
# 	

serial_debug:
	$(CC) $(CFLAGS) -c $(CONSTSRC)
	$(CC) $(CFLAGS) -c $(addprefix $(POLYCHORD), $(PCSRC))
	$(AR) $(LIBCHORD) $(PCCOBJ)
	$(CC) $(CFLAGS) -I. $(CCHORD) -c $(SOURCES)
	$(CC) $(CFLAGS) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. trades.f90 $(CCHORD) -o $(TARGET_SER)

serial_release:
	$(CC) $(COPT) -c $(CONSTSRC)
	$(CC) $(COPT) -c $(addprefix $(POLYCHORD), $(PCSRC))
	$(AR) $(LIBCHORD) $(PCCOBJ)
	$(CC) $(COPT) -I. $(CCHORD) -c $(SOURCES)
	$(CC) $(COPT) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades.f90 -o $(TARGET_SER)

omp_debug:
	$(CC) $(CFLAGS) -c $(CONSTSRC)
	$(CC) $(CFLAGS) $(COMP) -c $(addprefix $(POLYCHORD), $(PCSRC))
	$(AR) $(LIBCHORD) $(PCCOBJ)
	$(CC) $(CFLAGS) $(COMP) -I. $(CCHORD) -c $(SOURCES)
	$(CC) $(CFLAGS) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades.f90 -o $(TARGET_OMP)

omp_release:
	$(CC) $(COPT) -c $(CONSTSRC)
	$(CC) $(COPT) $(COMP) -c $(addprefix $(POLYCHORD), $(PCSRC))
	$(AR) $(LIBCHORD) $(PCCOBJ)
	$(CC) $(COPT) $(COMP) -I. $(CCHORD) -c $(SOURCES)
	$(CC) $(COPT) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades.f90 -o $(TARGET_OMP)

mpi_omp_debug:
	$(CC) $(CFLAGS) -c $(CONSTSRC)
	$(CC) $(CFLAGS) $(CMPI) $(COMP) -c $(addprefix $(POLYCHORD), $(PCSRC))
	$(AR) $(LIBCHORD) $(PCCOBJ)
	$(CC) $(CFLAGS) $(CMPI) $(COMP) -I. $(CCHORD) -c $(SOURCES)
	$(CC) $(CFLAGS) $(CMPI) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades.f90 -o $(TARGET_MPIOMP)

mpi_omp_release:
	$(CC) $(COPT) -c $(CONSTSRC)
	$(CC) $(COPT) $(CMPI) $(COMP) -c $(addprefix $(POLYCHORD), $(PCSRC))
	$(AR) $(LIBCHORD) $(PCCOBJ)
	$(CC) $(COPT) $(CMPI) $(COMP) $(CONSTOBJ) -I. $(CCHORD) -c $(SOURCES)
	$(CC) $(COPT) $(CMPI) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades.f90 -o $(TARGET_MPIOMP)

cleanall:
	rm -f $(POLYCHORD)*.o $(POLYCHORD)*.mod $(LIBCHORD) *.o *.mod $(TARGET_SER) $(TARGET_OMP)

clean_libchord:
	rm -f $(LIBCHORD) $(POLYCHORD)*.o $(POLYCHORD)*.mod

clean   :
	rm -f $(POLYCHORD)*.o $(POLYCHORD)*.mod *.o *.mod




