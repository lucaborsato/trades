# MY MAKEFILE
#
# gfortran compiler
# needed because some 'mpi_utils' has `use mpi` w/o #ifdef MPI #endif
CC = mpif90
# CC = /usr/bin/mpif90
# CC = /proper/path/to/mpif90/compiler
# f2py
F2PYCC = f2py
# WARNING: file .f2py_f2cmap must be present (it is within the src/ folder)

# compiler options
CFLAGS=-cpp -Og -g  -Wfatal-errors -Wall -fcheck=all -ffree-line-length-none -Wsurprising -Waliasing -Wimplicit-interface -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow
#CFLAGS+=-fmax-errors=1
#CFLAGS+=-Warray-temporaries
#CFLAGS+=-Wcharacter-truncation
#CFLAGS+=-Wconversion-extra
#CFLAGS+=-Wextra
# CFLAGS+=-pg # profiling 
COPT=-cpp -O3 -ffree-line-length-none -ffpe-trap=invalid,zero,overflow,underflow
# COPT+=-pg # profiling

# with openmp
COMP=-fopenmp
# for MPI
CMPI= -DMPI
# for F2PY
F2PYCF = -fPIC
TRADESLIB=pytrades_lib

LIBCHORD = libchord.a
AR = ar r
CCHORD = -L. -lchord

# executable
TARGET_SER_1=../bin/trades_s
TARGET_SER_2=../bin/trades_s_int_lm_bootstrap
TARGET_SER_3=../bin/trades_s_pik_pso
TARGET_SER_4=../bin/trades_s_grid
TARGET_OMP_1=../bin/trades_o
TARGET_OMP_2=../bin/trades_o_int_lm_bootstrap
TARGET_OMP_3=../bin/trades_o_pik_pso
TARGET_OMP_4=../bin/trades_o_grid
TARGET_MPIOMP=../bin/trades_mo
TARGET_PC=../bin/trades_polychord

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

# SOURCES = parameters.f90 random_trades.f90 convert_type.f90 lin_fit.f90 celestial_mechanics.f90  init_trades.f90 statistics.f90 timing.f90 rotations.f90 sort.f90 eq_motion.f90 output_files.f90 numerical_integrator.f90 radial_velocities.f90 transits.f90 ode_run.f90 derived_parameters_mod.f90 grid_search.f90 lm.f90 pikaia.f90 util_sort.f90 util_qmc.f90 opti_pso.f90 gaussian.f90 bootstrap.f90 PolyChord_driver.f90
SOURCES = parameters.f90 random_trades.f90 convert_type.f90 lin_fit.f90 celestial_mechanics.f90 eq_motion.f90 parameters_conversion.f90 transits.f90 init_trades.f90 statistics.f90 timing.f90 rotations.f90 sort.f90 gls_module.f90 output_files.f90 numerical_integrator.f90 radial_velocities.f90 ode_run.f90 derived_parameters_mod.f90 fitness_module.f90 grid_search.f90 lm.f90 pikaia.f90 util_sort.f90 util_qmc.f90 opti_pso.f90 gaussian.f90 bootstrap.f90 driver.f90 PolyChord_driver.f90
OBJECTS = $(SOURCES:.f90=.o)

# .PHONY: print_obj
# print_obj:
# 	echo " "
# 	echo $(PCSRC)
# 	echo " "
# 	echo $(PCOBJ)
# 	echo " "
# 	

constants_debug:
	$(CC) $(CFLAGS) $(F2PYCF) -c $(CONSTSRC)

constants_release:
	$(CC) $(COPT) $(F2PYCF) -c $(CONSTSRC)

polychord_debug: constants_debug
	$(CC) $(CFLAGS) $(F2PYCF) -c $(addprefix $(POLYCHORD), $(PCSRC))
	$(AR) $(LIBCHORD) $(PCCOBJ)
	
polychord_omp_debug: constants_debug
	$(CC) $(CFLAGS) $(COMP) $(F2PYCF) -c $(addprefix $(POLYCHORD), $(PCSRC))
	$(AR) $(LIBCHORD) $(PCCOBJ)
	
polychord_release: constants_release
	$(CC) $(COPT) $(F2PYCF) -c $(addprefix $(POLYCHORD), $(PCSRC))
	$(AR) $(LIBCHORD) $(PCCOBJ)

polychord_omp_release: constants_release
	$(CC) $(COPT) $(COMP) $(F2PYCF) -c $(addprefix $(POLYCHORD), $(PCSRC))
	$(AR) $(LIBCHORD) $(PCCOBJ)

serial_debug: constants_debug polychord_debug
	$(CC) $(CFLAGS) $(F2PYCF) -I. $(CCHORD) -c $(SOURCES)
	$(CC) $(CFLAGS) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades.f90 -o $(TARGET_SER_1)
	$(CC) $(CFLAGS) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_int_lm_bootstrap.f90 -o $(TARGET_SER_2)
	$(CC) $(CFLAGS) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_pik_pso.f90 -o $(TARGET_SER_3)
	$(CC) $(CFLAGS) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_grid.f90 -o $(TARGET_SER_4)
	$(CC) $(CFLAGS) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_polychord.f90 -o $(TARGET_PC)

serial_release: constants_release polychord_release
	$(CC) $(COPT) $(F2PYCF) -I. $(CCHORD) -c $(SOURCES)
	$(CC) $(COPT) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades.f90 -o $(TARGET_SER_1)
	$(CC) $(COPT) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_int_lm_bootstrap.f90 -o $(TARGET_SER_2)
	$(CC) $(COPT) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_pik_pso.f90 -o $(TARGET_SER_3)
	$(CC) $(COPT) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_grid.f90 -o $(TARGET_SER_4)
	$(CC) $(COPT) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_polychord.f90 -o $(TARGET_PC)
	
omp_debug: constants_debug polychord_omp_debug
	$(CC) $(CFLAGS) $(F2PYCF) $(COMP) -I. $(CCHORD) -c $(SOURCES)
	$(CC) $(CFLAGS) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades.f90 -o $(TARGET_OMP_1)
	$(CC) $(CFLAGS) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_int_lm_bootstrap.f90 -o $(TARGET_OMP_2)
	$(CC) $(CFLAGS) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_pik_pso.f90 -o $(TARGET_OMP_3)
	$(CC) $(CFLAGS) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_grid.f90 -o $(TARGET_OMP_4)
	$(CC) $(CFLAGS) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_polychord.f90 -o $(TARGET_PC)

only_pytrades_debug:
	$(F2PYCC) --f90flags='$(CFLAGS) $(F2PYCF) $(COMP)' -lgomp -c $(TRADESLIB).f90 *.o -m $(TRADESLIB)
	
pyomp_debug: omp_debug only_pytrades_debug
	cp $(TRADESLIB).so ../pytrades/$(TRADESLIB).so
	
omp_release: constants_release polychord_omp_release
	$(CC) $(COPT) $(COMP) $(F2PYCF) -I. $(CCHORD) -c $(SOURCES)
	$(CC) $(COPT) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades.f90 -o $(TARGET_OMP_1)
	$(CC) $(COPT) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_int_lm_bootstrap.f90 -o $(TARGET_OMP_2)
	$(CC) $(COPT) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_pik_pso.f90 -o $(TARGET_OMP_3)
	$(CC) $(COPT) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_grid.f90 -o $(TARGET_OMP_4)
	$(CC) $(COPT) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_polychord.f90 -o $(TARGET_PC)

only_pytrades_release: 
	$(F2PYCC) --f90flags='$(COPT) $(F2PYCF) $(COMP)' -lgomp -c $(TRADESLIB).f90 *.o -m $(TRADESLIB)
	
pyomp_release: omp_release only_pytrades_release
	cp $(TRADESLIB).so ../pytrades/$(TRADESLIB).so

mpi_omp_debug: constants_debug
	$(CC) $(CFLAGS) $(CMPI) $(COMP) -c $(addprefix $(POLYCHORD), $(PCSRC))
	$(AR) $(LIBCHORD) $(PCCOBJ)
	$(CC) $(CFLAGS) $(CMPI) $(COMP) -I. $(CCHORD) -c $(SOURCES)
	$(CC) $(CFLAGS) $(CMPI) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades.f90 -o $(TARGET_MPIOMP)
	$(CC) $(CFLAGS) $(CMPI) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_polychord.f90 -o $(TARGET_PC)

mpi_omp_release: constants_release
	$(CC) $(COPT) -c $(CONSTSRC)
	$(CC) $(COPT) $(CMPI) $(COMP) -c $(addprefix $(POLYCHORD), $(PCSRC))
	$(AR) $(LIBCHORD) $(PCCOBJ)
	$(CC) $(COPT) $(CMPI) $(COMP) $(CONSTOBJ) -I. $(CCHORD) -c $(SOURCES)
	$(CC) $(COPT) $(CMPI) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades.f90 -o $(TARGET_MPIOMP)
	$(CC) $(COPT) $(CMPI) $(COMP) $(CONSTOBJ) $(PCOBJ) $(OBJECTS) -I. $(CCHORD) trades_polychord.f90 -o $(TARGET_PC)

# cleanall:
# 	rm -f $(POLYCHORD)*.o $(POLYCHORD)*.mod $(LIBCHORD) *.o *.mod $(TARGET_SER_1) $(TARGET_SER_2) $(TARGET_SER_3) $(TARGET_OMP) $(TARGET_OMP_2) $(TARGET_OMP_3) $(TARGET_PC)

cleanall: clean clean_libchord clean_pytrades

clean_libchord:
	rm -f $(LIBCHORD) $(POLYCHORD)*.o $(POLYCHORD)*.mod

clean   :
	rm -f $(POLYCHORD)*.o $(POLYCHORD)*.mod *.o *.mod

clean_pytrades: clean
	rm -f $(TRADESLIB).so ../pytrades/$(TRADESLIB).so


