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

ifort ?= 0
debug ?= 0 

ifeq ($(ifort),1)
    FC = ifort 
else
    FC = gfortran
endif 

ifeq ($(ifort),1)
	## IFORT OPTIONS ##
	FLAGS        = -module $(objdir) -L$(objdir) -I/home/robinson/apps/netcdf/netcdf/include
	LFLAGS		 = -L/home/robinson/apps/netcdf/netcdf/lib -lnetcdf

	ifeq ($(debug), 1)
	    DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0
	    # -w 
	else
	    DFLAGS   = -vec-report0 -O3
	endif
else
	## GFORTRAN OPTIONS ##
	FLAGS        = -I$(objdir) -J$(objdir) -I/opt/local/include
	LFLAGS		 = -L/opt/local/lib -lnetcdff -lnetcdf

	ifeq ($(debug), 1)
	    DFLAGS   = -w -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow -fbacktrace -fcheck=all
	else
	    DFLAGS   = -O3
	endif
endif


## Individual libraries or modules ##
$(objdir)/ncio.o: ../ncio/ncio.f90
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
ccsm3: $(objdir)/ncio.o $(objdir)/geodesic.o $(objdir)/planet.o $(objdir)/projection_oblimap2.o $(objdir)/coordinates.o
	$(FC) $(FLAGS) -o test_ccsm3.x $^ test_ccsm3.f90 $(LFLAGS)
	@echo " "
	@echo "    test_ccsm3.x is ready."
	@echo " "

# Program to test distance calculations using the geographiclib library
geodinverse: $(objdir)/planet.o $(objdir)/geodesic.o
	$(FC) $(FLAGS) -o geodinverse.x $^ geodinverse.f90
	@echo " "
	@echo "    geodinverse.x is ready."
	@echo " "

clean:
	rm -f test_ccsm3.x $(objdir)/*.o $(objdir)/*.mod

# cleanall: cleansico cleanrembo cleansicoX
