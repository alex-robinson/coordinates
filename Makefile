.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make ccsm3      : compiles the main program test_ccsm3.x"
	@echo " make clean      : cleans object files"
	@echo ""

objdir = .obj

## Local ##
FC			 = gfortran
FLAGS        = -I$(objdir) -J$(objdir) -I/opt/local/include
DEBUGFLAGS   = -w -pg -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
RELEASEFLAGS = -O3
LFLAGS		 = -L/opt/local/lib -lnetcdff -lnetcdf

## Cluster ##
# FC			 = ifort
# FLAGS        = -module $(objdir) -L$(objdir) -I/home/robinson/apps/netcdf/netcdf/include
# DEBUGFLAGS   = -w -C -traceback -ftrapuv -fpe0 -check all -vec-report0
# RELEASEFLAGS = -vec-report0 -O3
# LFLAGS		 = -L/home/robinson/apps/netcdf/netcdf/lib -lnetcdf


#DIRFLAGS = -I./.obj -J./.obj
#FFLAGS	= $(DIRFLAGS) -O3 
#FFLAGS	= $(DIRFLAGS) -w -pg -fbacktrace -fbounds-check #-g
#FFLAGS	= $(DIRFLAGS) -w -pg -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
#LDFLAGS = $(FFLAGS) -I/opt/local/include
# LIB = -L/opt/local/lib -lnetcdff -lnetcdf

# FC	= ifort
# DIRFLAGS = -module .obj -L./.obj
# FFLAGS	= $(DIRFLAGS) -w -vec-report0 -O3 -xSSE4.1 #-pg -mcmodel medium -shared-intel
# #FFLAGS	= $(DIRFLAGS) -w -C -traceback -vec-report0 -xSSE4.1
# #FFLAGS	= $(DIRFLAGS) -w -C -traceback -ftrapuv -fpe0 -check all -vec-report0 -xSSE4.1
# LDFLAGS = $(FFLAGS) -I/home/robinson/apps/netcdf/netcdf/include
# LIB = -L/home/robinson/apps/netcdf/netcdf/lib -lnetcdf

## Individual libraries or modules ##
$(objdir)/ncio3.o: ../ncio/ncio3.f90
	$(FC) $(FLAGS) -c -o $@ $<

$(objdir)/planet.o: planet.f90
	$(FC) $(FLAGS) -c -o $@ $<

$(objdir)/geodesic.o: geodesic.f90
	$(FC) $(FLAGS) -c -o $@ $<

$(objdir)/projection_oblimap2.o: projection_oblimap2.f90
	$(FC) $(FLAGS) -c -o $@ $<

$(objdir)/coordinates.o: coordinates.f90
	$(FC) $(FLAGS) -c -o $@ $<

## Complete programs

# Program to test interpolations of CCSM3 data
ccsm3: $(objdir)/ncio3.o $(objdir)/geodesic.o $(objdir)/planet.o $(objdir)/projection_oblimap2.o $(objdir)/coordinates.o
	$(FC) $(FLAGS) -o test_CCSM3.x $^ test_CCSM3.f90 $(LFLAGS)
	@echo " "
	@echo "    test_CCSM3.x is ready."
	@echo " "

# Program to test distance calculations using the geographiclib library
geodinverse: $(objdir)/planet.o $(objdir)/geodesic.o
	$(FC) $(FLAGS) -o geodinverse.x $^ geodinverse.f90
	@echo " "
	@echo "    geodinverse.x is ready."
	@echo " "

clean:
	rm -f test_ccsm3.x *.o *.mod

# cleanall: cleansico cleanrembo cleansicoX
