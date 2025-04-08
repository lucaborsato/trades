# MY MAKEFILE
#
# UNSET LDFLAGS, DUE TO ANACONDA/ENV ISSUE
undefine LDFLAGS

# gfortran compiler
CC = gfortran
# CC = ifort

F77 = gfortran
F90 = gfortran
# F77 = ifort

# f2py
FVENDOR=$(CC)
# F2PYCC = f2py --fcompiler=$(FVENDOR) --f90exec=$(CC) --f77exec=$(F77)
# F2PYCC = FC=$(CC) f2py
F2PYCC = FC=$(CC) python -m numpy.f2py

# WARNING: file .f2py_f2cmap must be present (it is within the src/ folder)

# compiler options
# CFLAGS=-cpp -ffree-line-length-none -Wfatal-errors -Wall -fcheck=all -Wsurprising -Waliasing -Wimplicit-interface -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -Og -g 
# ====
CFLAGS=-ffree-line-length-none
CFLAGS+=-Wfatal-errors
CFLAGS+=-Wall
CFLAGS+=-fcheck=all
CFLAGS+=-Wsurprising
CFLAGS+=-Waliasing
CFLAGS+=-Wimplicit-interface
CFLAGS+=-fbacktrace
CFLAGS+=-ffpe-trap=invalid,zero,overflow,underflow 
CFLAGS+=-Og
CFLAGS+=-g
CFLAGS+=-cpp
# ====
# CFLAGS+=-pg # profiling 

COPT=-cpp
COPT+=-ffree-line-length-none
COPT+=-Ofast
# COPT+=-O3

# with openmp
GOMP=-lgomp
# GOMP=
COMP=-fopenmp
# for MPI
CMPI= -DMPI
# for F2PY
F2PYCF = -fPIC -shared
TRADESLIB=pytrades_lib
# TRADESLIB0=pytrades_lib.cpython-38-x86_64-linux-gnu.so
TRADESLIB0=pytrades_lib.cpython-*-x86_64-linux-gnu.so

# executable
TARGET_0=../bin/trades_integrates
TARGET_SER_2=../bin/trades_s_int_lm_bootstrap
TARGET_SER_3=../bin/trades_s_pik_pso
TARGET_SER_4=../bin/trades_s_grid
TARGET_SER_5=../bin/trades_s_test
TARGET_SER_6=../bin/trades_s_timing_fitness
TARGET_OMP_2=../bin/trades_o_int_lm_bootstrap
TARGET_OMP_3=../bin/trades_o_pik_pso
TARGET_OMP_4=../bin/trades_o_grid
TARGET_OMP_5=../bin/trades_o_test
TARGET_OMP_6=../bin/trades_o_timing_fitness

FSRCS = constants.f90 sort.f90 custom_type.f90 parameters.f90 random_trades.f90 convert_type.f90 statistics.f90 lin_fit.f90 rotations.f90 celestial_mechanics.f90 eq_motion.f90 parameters_conversion.f90 linear_ephem.f90 utils.f90 init_trades.f90 timing.f90 gls_module.f90 output_files.f90 numerical_integrator.f90 transits.f90 radial_velocities.f90 ode_run.f90 derived_parameters_mod.f90 fitness_module.f90 grid_search.f90 lm.f90 pikaia.f90 util_sort.f90 util_qmc.f90 opti_pso.f90 gaussian.f90 bootstrap.f90 driver.f90
FOBJS = $(FSRCS:.f90=.o)

# .PHONY: print_obj
# print_obj:
# 	echo " "
# 	echo $(PCSRC)
# 	echo " "
# 	echo $(PCOBJ)
# 	echo " "
# 	


# -------
# TESTING
# testing_release:
# 	$(CC) $(COPT) $(F2PYCF) -I. -c $(FSRCS)
# 	$(CC) $(COPT) $(F2PYCF) $(FOBJS) -I. trades_integrates.f90 -o $(TARGET_0)
# 	$(CC) $(COPT) $(F2PYCF) $(FOBJS) -I. trades_int_lm_bootstrap.f90 -o $(TARGET_SER_2)
# -------

# serial debug
FSERDBG= $(CC) $(CFLAGS)
PSERDBG= $(F2PYCF) $(CFLAGS)
# serial release
FSERREL= $(CC) $(COPT)
PSERREL= $(F2PYCF) $(COPT)
# parallel debug
FPARDBG= $(CC) $(CFLAGS) $(COMP)
PPARDBG= $(F2PYCF) $(CFLAGS) $(COMP)
# parallel release
FPARREL=$(CC) $(COPT) $(COMP)
PPARREL=$(F2PYCF) $(COPT) $(COMP)

echoflag:
	@echo "F90  SERIAL DEBUG    : " $(FSERDBG)
	@echo "F2PY SERIAL DEBUG    : " $(PSERDBG)
	@echo "F90  SERIAL RELEASE  : " $(FSERREL)
	@echo "F2PY SERIAL RELEASE  : " $(PSERREL)
	@echo "F90  PARALLEL DEBUG  : " $(FPARDBG)
	@echo "F2PY PARALLEL DEBUG  : " $(PPARDBG)
	@echo "F90  PARALLEL RELEASE: " $(FPARREL)
	@echo "F2PY PARALLEL RELEASE: " $(PPARREL)
	@echo ""
	@echo "FFLAGS='$(PSERDBG)' $(F2PYCC) --verbose -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)"
	@echo ""
	@echo "$(F2PYCC) --verbose --f90exec=$(F90) --f77exec=$(F77) --f90flags='$(F2PYCF) $(CFLAGS)' -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)"
	@echo ""

# compile the sources to objects
obj_compile_test:
	$(FSERDBG) $(F2PYCF) -I. -c constants.f90
	$(FSERDBG) $(F2PYCF) -I. -c sort.f90
	$(FSERDBG) $(F2PYCF) -I. -c custom_type.f90
	$(FSERDBG) $(F2PYCF) -I. -c parameters.f90
	$(FSERDBG) $(F2PYCF) -I. -c random_trades.f90
	$(FSERDBG) $(F2PYCF) -I. -c convert_type.f90
	$(FSERDBG) $(F2PYCF) -I. -c statistics.f90
	$(FSERDBG) $(F2PYCF) -I. -c lin_fit.f90
	$(FSERDBG) $(F2PYCF) -I. -c rotations.f90
	$(FSERDBG) $(F2PYCF) -I. -c celestial_mechanics.f90
	$(FSERDBG) $(F2PYCF) -I. -c eq_motion.f90
	$(FSERDBG) $(F2PYCF) -I. -c parameters_conversion.f90
	$(FSERDBG) $(F2PYCF) -I. -c linear_ephem.f90
	$(FSERDBG) $(F2PYCF) -I. -c utils.f90
	$(FSERDBG) $(F2PYCF) -I. -c init_trades.f90
	$(FSERDBG) $(F2PYCF) -I. -c timing.f90
	$(FSERDBG) $(F2PYCF) -I. -c gls_module.f90
	$(FSERDBG) $(F2PYCF) -I. -c output_files.f90
	$(FSERDBG) $(F2PYCF) -I. -c numerical_integrator.f90
	$(FSERDBG) $(F2PYCF) -I. -c transits.f90
	$(FSERDBG) $(F2PYCF) -I. -c radial_velocities.f90
	$(FSERDBG) $(F2PYCF) -I. -c ode_run.f90
	$(FSERDBG) $(F2PYCF) -I. -c derived_parameters_mod.f90
	$(FSERDBG) $(F2PYCF) -I. -c fitness_module.f90
	$(FSERDBG) $(F2PYCF) -I. -c grid_search.f90
	$(FSERDBG) $(F2PYCF) -I. -c lm.f90
	$(FSERDBG) $(F2PYCF) -I. -c pikaia.f90
	$(FSERDBG) $(F2PYCF) -I. -c util_sort.f90
	$(FSERDBG) $(F2PYCF) -I. -c util_qmc.f90
	$(FSERDBG) $(F2PYCF) -I. -c opti_pso.f90
	$(FSERDBG) $(F2PYCF) -I. -c gaussian.f90
	$(FSERDBG) $(F2PYCF) -I. -c bootstrap.f90
	$(FSERDBG) $(F2PYCF) -I. -c driver.f90
obj_serial_debug:
	$(FSERDBG) $(F2PYCF) -I. -c $(FSRCS)
obj_serial_release:
	$(FSERREL) $(F2PYCF) -I. -c $(FSRCS)
obj_parallel_debug:
	$(FPARDBG) $(F2PYCF) -I. -c $(FSRCS)
obj_parallel_release:
	$(FPARREL) $(F2PYCF) -I. -c $(FSRCS)

# compile the executables FORTRAN only
trades_serial_debug: obj_serial_debug
	@mkdir -p ../bin
	$(FSERDBG) $(FOBJS) -I. trades_grid.f90 -o $(TARGET_SER_4)
	$(FSERDBG) $(FOBJS) -I. trades_int_lm_bootstrap.f90 -o $(TARGET_SER_2)
	$(FSERDBG) $(FOBJS) -I. trades_integrates.f90 -o $(TARGET_0)
	$(FSERDBG) $(FOBJS) -I. trades_pik_pso.f90 -o $(TARGET_SER_3)
	$(FSERDBG) $(FOBJS) -I. trades_timing_fitness.f90 -o $(TARGET_SER_6)

trades_serial_release: obj_serial_release
	@mkdir -p ../bin
	$(FSERREL) $(FOBJS) -I. trades_grid.f90 -o $(TARGET_SER_4)
	$(FSERREL) $(FOBJS) -I. trades_int_lm_bootstrap.f90 -o $(TARGET_SER_2)
	$(FSERREL) $(FOBJS) -I. trades_integrates.f90 -o $(TARGET_0)
	$(FSERREL) $(FOBJS) -I. trades_pik_pso.f90 -o $(TARGET_SER_3)
	$(FSERREL) $(FOBJS) -I. trades_timing_fitness.f90 -o $(TARGET_SER_6)

trades_parallel_debug: obj_parallel_debug
	@mkdir -p ../bin
	$(FPARDBG) $(FOBJS) -I. trades_grid.f90 -o $(TARGET_OMP_4)
	$(FPARDBG) $(FOBJS) -I. trades_int_lm_bootstrap.f90 -o $(TARGET_OMP_2)
	$(FPARDBG) $(FOBJS) -I. trades_integrates.f90 -o $(TARGET_0)
	$(FPARDBG) $(FOBJS) -I. trades_pik_pso.f90 -o $(TARGET_OMP_3)
	$(FPARDBG) $(FOBJS) -I. trades_timing_fitness.f90 -o $(TARGET_OMP_6)

trades_parallel_release: obj_parallel_release
	@mkdir -p ../bin
	$(FPARREL) $(FOBJS) -I. trades_grid.f90 -o $(TARGET_OMP_4)
	$(FPARREL) $(FOBJS) -I. trades_int_lm_bootstrap.f90 -o $(TARGET_OMP_2)
	$(FPARREL) $(FOBJS) -I. trades_integrates.f90 -o $(TARGET_0)
	$(FPARREL) $(FOBJS) -I. trades_pik_pso.f90 -o $(TARGET_OMP_3)
	$(FPARREL) $(FOBJS) -I. trades_timing_fitness.f90 -o $(TARGET_OMP_6)

# compile the python library only

copy_pytrades:
	cp $(TRADESLIB0) ../pytrades/$(TRADESLIB).so

pylib_serial_debug:
	FFLAGS='$(PSERDBG)' $(F2PYCC) --verbose -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)
#	FFLAGS='$(F2PYCF) $(CFLAGS)' $(F2PYCC) --verbose -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)
#	$(F2PYCC) --verbose --f90exec=$(F90) --f77exec=$(F77) --f90flags='$(F2PYCF) $(CFLAGS)' -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)

pytrades_serial_debug: obj_serial_debug pylib_serial_debug

pylib_serial_release:
	FFLAGS='$(PSERREL)' $(F2PYCC) --quiet -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)
#	FFLAGS='$(F2PYCF) $(COPT)' $(F2PYCC) --quiet -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)
#	$(F2PYCC) --quiet --f90flags='$(F2PYCF) $(COPT)' -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)

pytrades_serial_release: obj_serial_release pylib_serial_release

pylib_parallel_debug:
	FFLAGS='$(PPARDBG)' $(F2PYCC) --verbose $(GOMP) -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)
#	$(F2PYCC) --verbose --f90flags='$(F2PYCF) $(CFLAGS) $(COMP)' $(GOMP) -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)

pytrades_parallel_debug: obj_parallel_debug pylib_parallel_debug copy_pytrades

pylib_parallel_release:
	FFLAGS='$(PPARREL)' $(F2PYCC) --quiet $(GOMP) -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)
#	$(F2PYCC) --quiet --f90flags='$(F2PYCF) $(COPT) $(COMP)' $(GOMP) -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)

pytrades_parallel_release: obj_parallel_release pylib_parallel_release copy_pytrades

# compile both FORTRAN and python library
full_serial_debug: trades_serial_debug pylib_serial_debug copy_pytrades

full_serial_release: trades_serial_release pylib_serial_release copy_pytrades

full_parallel_debug: trades_parallel_debug pylib_parallel_debug copy_pytrades

full_parallel_release: trades_parallel_release pylib_parallel_release copy_pytrades

test_debug: pytrades_serial_debug
	$(CC) $(PSERDBG) -I. -c $(TRADESLIB).f90
	@mkdir -p ../bin
	$(FSERDBG) $(FOBJS) $(TRADESLIB).o -I. trades_test.f90 -o $(TARGET_SER_5)

test_release: pytrades_serial_release
	$(CC) $(PSERREL) -I. -c $(TRADESLIB).f90
	@mkdir -p ../bin
	$(FSERREL) $(FOBJS) $(TRADESLIB).o -I. trades_test.f90 -o $(TARGET_SER_5)

#test_meson: obj_parallel_release
#	FFLAGS='$(COPT) $(COMP)' $(F2PYCC) --quiet $(GOMP) -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)
#	FFLAGS='' $(F2PYCC) --quiet $(GOMP) -c $(TRADESLIB).f90 $(FOBJS) -m $(TRADESLIB)
#	copy_pytrades

# clean objs, mod, and lib
clean:
	rm -f *.o *.mod

clean_pytrades: clean
	rm -f $(TRADESLIB0) ../pytrades/$(TRADESLIB).so

clean_exe:
	rm -f $(TARGET_0)
	rm -f $(TARGET_SER_2) $(TARGET_SER_3) $(TARGET_SER_4) $(TARGET_SER_5)
	rm -f $(TARGET_OMP_2) $(TARGET_OMP_3) $(TARGET_OMP_4) $(TARGET_OMP_5)
	
cleanall: clean clean_pytrades clean_exe

