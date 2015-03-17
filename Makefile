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

# PATH options
objdir = .obj

# netcdf_inc = /usr/include
# netcdf_lib = /usr/lib
netcdf_inc = /opt/local/include
netcdf_lib = /opt/local/lib
netcdf_inc_ifort = /home/robinson/apps/netcdf/netcdf/include
netcdf_lib_ifort = /home/robinson/apps/netcdf/netcdf/lib

# Command-line options at make call
ifort ?= 0
debug ?= 0 

ifeq ($(ifort),1)
    FC = ifort 
else
    FC = gfortran
endif 

ifeq ($(ifort),1)
	## IFORT OPTIONS ##
	FLAGS        = -module $(objdir) -L$(objdir) -I$(netcdf_inc_ifort)
	LFLAGS		 = -L$(netcdf_lib_ifort) -lnetcdf
	SFLAGS       = 

	ifeq ($(debug), 1)
	    DFLAGS   = -C -traceback -ftrapuv -fpe0 -check all -vec-report0
	    # -w 
	else
	    DFLAGS   = -vec-report0 -O3
	endif
else
	## GFORTRAN OPTIONS ##
	FLAGS        = -I$(objdir) -J$(objdir) -I$(netcdf_inc)
	LFLAGS		 = -L$(netcdf_lib) -lnetcdff -lnetcdf
	SFLAGS       = 

	ifeq ($(debug), 1)
	    DFLAGS   = -w -p -ggdb -ffpe-trap=invalid,zero,overflow,underflow \
	               -fbacktrace -fcheck=all -fbackslash
	else
	    DFLAGS   = -O3 -fbackslash
	endif
endif

## Individual libraries or modules ##
$(objdir)/ncio.o: ncio.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/index.o: index.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/interp1D.o: interp1D.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/interp2D.o: interp2D.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/loess.o: loess.f90 $(objdir)/interp1D.o $(objdir)/index.o
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/gaussian_filter.o: gaussian_filter.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/mod_toms526.o: mod_toms526.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/interp_time.o: interp_time.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/polygons.o: polygons.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/planet.o: planet.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/geodesic.o: geodesic.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/projection_oblimap2.o: projection_oblimap2.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/coordinates.o: coordinates.f90 $(objdir)/ncio.o $(objdir)/planet.o $(objdir)/geodesic.o \
						 $(objdir)/projection_oblimap2.o $(objdir)/gaussian_filter.o
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/subset.o: subset.f90 $(objdir)/coordinates.o
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/subset2.o: subset2.f90 $(objdir)/coordinates.o
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

## Complete programs

# coordinates static library - using subset2
coord-static: $(objdir)/ncio.o $(objdir)/index.o $(objdir)/polygons.o \
	$(objdir)/geodesic.o $(objdir)/planet.o $(objdir)/projection_oblimap2.o \
	$(objdir)/interp1D.o $(objdir)/interp2D.o $(objdir)/mod_toms526.o $(objdir)/gaussian_filter.o \
	$(objdir)/loess.o $(objdir)/interp_time.o \
	$(objdir)/subset2.o $(objdir)/coordinates.o
	ar rc libcoordinates.a $^
	@echo " "
	@echo "    libcoordinates.a is ready."
	@echo " "

# coordinates shared library - using subset2
coord-shared: $(objdir)/ncio.o $(objdir)/index.o $(objdir)/polygons.o \
	$(objdir)/geodesic.o $(objdir)/planet.o $(objdir)/projection_oblimap2.o \
	$(objdir)/interp1D.o $(objdir)/interp2D.o $(objdir)/interp_time.o \
	$(objdir)/subset2.o $(objdir)/coordinates.o
	$(FC) $(DFLAGS) $(FLAGS) -shared -fPIC -o libcoordinates.so $^ $(LFLAGS)
	@echo " "
	@echo "    libcoordinates.so is ready."
	@echo " "

# coordinates shared library - using subset
coord0-shared: $(objdir)/ncio.o $(objdir)/index.o $(objdir)/polygons.o \
	$(objdir)/geodesic.o $(objdir)/planet.o $(objdir)/projection_oblimap2.o \
	$(objdir)/interp1D.o $(objdir)/interp2D.o $(objdir)/interp_time.o \
	$(objdir)/subset.o $(objdir)/coordinates.o
	$(FC) $(DFLAGS) $(FLAGS) -shared -fPIC -o libcoordinates0.so $^ $(LFLAGS)
	@echo " "
	@echo "    libcoordinates0.so is ready."
	@echo " "

# Program to test interpolations of CCSM3 data
ccsm3: coord-static
	$(FC) $(DFLAGS) $(FLAGS) -o test_ccsm3.x test_ccsm3.f90 libcoordinates.a -L. $(LFLAGS)
	@echo " "
	@echo "    test_ccsm3.x is ready."
	@echo " "

etopo: coord-static
	$(FC) $(DFLAGS) $(FLAGS) -o test_etopo.x test_etopo.f90 libcoordinates.a -L. $(LFLAGS)
	@echo " "
	@echo "    test_etopo.x is ready."
	@echo " "

test_subset: coord0
	$(FC) $(DFLAGS) $(FLAGS) -o test_subset.x test_subset.f90 -L. -lcoordinates0 $(LFLAGS)
	@echo " "
	@echo "    test_subset.x is ready."
	@echo " "

test_subset2: coord
	$(FC) $(DFLAGS) $(FLAGS) -o test_subset2.x test_subset2.f90 -L. -lcoordinates $(LFLAGS)
	@echo " "
	@echo "    test_subset2.x is ready."
	@echo " "

test_interp: coord-static
	$(FC) $(DFLAGS) $(FLAGS) -o test_interp.x test_interp.f90 libcoordinates.a -L. $(LFLAGS)
	@echo " "
	@echo "    test_interp.x is ready."
	@echo " "

test_climber: coord-static
	$(FC) $(DFLAGS) $(FLAGS) -o test_climber.x test_climber.f90 libcoordinates.a -L. $(LFLAGS)
	@echo " "
	@echo "    test_climber.x is ready."
	@echo " "

ccsm3diff: coord-static
	$(FC) $(DFLAGS) $(FLAGS) -o test_ccsm3_diffusion.x test_ccsm3_diffusion.f90 libcoordinates.a -L. $(LFLAGS)
	@echo " "
	@echo "    test_ccsm3_diffusion.x is ready."
	@echo " "

test_loess: $(objdir)/ncio.o $(objdir)/interp1D.o $(objdir)/index.o $(objdir)/loess.o 
	$(FC) $(DFLAGS) $(FLAGS) -o test_loess.x $^ test_loess.f90 $(LFLAGS)
	@echo " "
	@echo "    test_loess.x is ready."
	@echo " "

test_nat: tile/t4gen.f tile/t4man.f tile/t4int.f tile/t4que.f
	$(FC) $(DFLAGS) $(FLAGS) -o test_nat.x $^ test_nat.f90
	@echo " "
	@echo "    test_nat.x is ready."
	@echo " "

# Program to test distance calculations using the geographiclib library
geodinverse: $(objdir)/planet.o $(objdir)/geodesic.o
	$(FC) $(DFLAGS) $(FLAGS) -o geodinverse.x $^ geodinverse.f90
	@echo " "
	@echo "    geodinverse.x is ready."
	@echo " "

# Program to test polygon calculations
Pointpoly: $(objdir)/polygons.o
	$(FC) $(DFLAGS) $(FLAGS) -o Pointpoly.x $^ Pointpoly.f90 $(LFLAGS)
	@echo " "
	@echo "    Pointpoly.x is ready."
	@echo " "

# Program to test custom polygon calculations
poly: $(objdir)/polygons.o $(objdir)/index.o
	$(FC) $(DFLAGS) $(FLAGS) -o test_polygon.x $^ test_polygon.f90 $(LFLAGS)
	@echo " "
	@echo "    test_polygon.x is ready."
	@echo " "

clean:
	rm -r -f *.x *.dSYM $(objdir)/*.o $(objdir)/*.mod *.so 

# cleanall: cleansico cleanrembo cleansicoX
