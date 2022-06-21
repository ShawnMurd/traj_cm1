#################################################################################
#                         MAKEFILE FOR TRAJ_CM1 CODE                            #
#                                                                               #
#-------------------------------------------------------------------------------#
# Shawn Murdzek, 15 June 2022                                                   #
#-------------------------------------------------------------------------------#
# Based on Ryan Hastings' cm1_pp code                                           #
#################################################################################
#################################################################################
# CONFIGURATION VARIABLES                                                       #

F90=gfortran # compiler to use
NCINC=-I/gpfs/group/pmm116/default/sw/netcdf-4.4.1.1-gcc-5.3.1/include # location of include files for netcdf
NCLIB=-L/gpfs/group/pmm116/default/sw/netcdf-4.4.1.1-gcc-5.3.1/lib # location of libraries for netcdf

#################################################################################
# YOU PROBABLY DON'T WANT TO CHANGE ANYTHING BELOW THIS                         #
#################################################################################

ALL_PROGRAMS = traj

CM1PP_code: $(ALL_PROGRAMS)

.PHONY : CM1PP_code

OBJT = globals.o nc_inout.o traj_utils.o traj.o

#################################################################################

traj :	$(OBJT)
	$(F90) -o ../bin/traj $(OBJT) -fcheck=bounds $(NCLIB) -lm -lnetcdff

traj.o :	traj.f90
		$(F90) -c traj.f90 -fcheck=bounds $(OPTS)

nc_inout.o :	nc_inout.f90
		$(F90) -c nc_inout.f90 -fcheck=bounds $(OPTS) $(NCINC) $(NCLIB)

traj_utils.o :	traj_utils.f90
		$(F90) -c traj_utils.f90 -fcheck=bounds $(OPTS)

globals.o :	globals.f90
		$(F90) -c globals.f90 -fcheck=bounds $(OPTS)

clean:
		rm -rf *o *mod
