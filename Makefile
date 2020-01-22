# the directories containing the libraries
LIBDIR =

# fortran 90 compiler and compiler flags
#F90 = mpiifort
F90 = mpif90
#F90FLAGS = -qsmp=omp -O3 -qsmallstack
F90FLAGS = -qopenmp -fpp -O3 -mcmodel=large 


#----------- end of user configuration parameters ------------

all: cosmology.mod functions.mod main.f90 vars.mod io.mod tr_ly_utils.mod tr_ly.mod tspin.mod
	$(F90) $(F90FLAGS) main.f90 cosmology.f90 vars.f90 tr_ly_utils.f90 tr_ly.f90 functions.f90 io.f90 tspin.f90 $(LIBDIR) -o licorice 


io.mod: io.f90 vars.mod
	$(F90) $(F90FLAGS)  -c io.f90 $(LIBDIR)

cosmology.mod: cosmology.f90  vars.mod
	$(F90) $(F90FLAGS)  -c cosmology.f90 $(LIBDIR)

tr_ly_utils.mod: tr_ly_utils.f90 vars.mod  
	$(F90) $(F90FLAGS)  -c tr_ly_utils.f90 $(LIBDIR)

tr_ly.mod: tr_ly.f90 tr_ly_utils.mod vars.mod tspin.mod
	$(F90) $(F90FLAGS)  -c tr_ly.f90 $(LIBDIR)

functions.mod: functions.f90 vars.mod
	$(F90) $(F90FLAGS)  -c functions.f90 $(LIBDIR)

vars.mod: vars.f90
	$(F90) $(F90FLAGS)  -c vars.f90 $(LIBDIR)

tspin.mod: tspin.f90 vars.mod
	$(F90) $(F90FLAGS)  -c tspin.f90 $(LIBDIR)


clean:
	rm *.mod *.o licorice
