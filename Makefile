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
srcdir = src
testdir = tests

#netcdf_inc = /usr/include
#netcdf_lib = /usr/lib
netcdf_inc = /opt/local/include
netcdf_lib = /opt/local/lib
netcdf_inc_ifort = /home/robinson/apps/netcdf/netcdf/include
netcdf_lib_ifort = /home/robinson/apps/netcdf/netcdf/lib

NETCDF_FORTRANROOT = /home/fispalma22/work/librairies/netcdflib
INC_NC  = -I${NETCDF_FORTRANROOT}/include
LIB_NC  = -L${NETCDF_FORTRANROOT}/lib -lnetcdff -L${NETCDF_CROOT}/lib -lnetcdf 

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
	FLAGS        = -heap-arrays -module $(objdir) -L$(objdir) $(INC_NC) 
	LFLAGS	     = $(LIB_NC) 
	SFLAGS       = 

	ifeq ($(debug), 1)
	    DFLAGS   = -C -g -traceback -ftrapuv -fpe0 -check all 
	    # -w 
	else
	    DFLAGS   = -O3
	endif
else
	## GFORTRAN OPTIONS ##
	FLAGS        = -I$(objdir) -J$(objdir) -I$(netcdf_inc)
	LFLAGS		 = -L$(netcdf_lib) -lnetcdff -lnetcdf
	SFLAGS       = 

	ifeq ($(debug), 1)
	    DFLAGS   = -w -pg -ggdb -ffpe-trap=invalid,zero,overflow,underflow \
	               -fbacktrace -fcheck=all -fbackslash
	else
	    DFLAGS   = -O3 -fbackslash
	endif
endif

## Individual libraries or modules ##
$(objdir)/ncio.o: $(srcdir)/ncio.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/coord_constants.o: $(srcdir)/coord_constants.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/index.o: $(srcdir)/index.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/interp1D.o: $(srcdir)/interp1D.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/interp2D.o: $(srcdir)/interp2D.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/loess.o: $(srcdir)/loess.f90 $(objdir)/interp1D.o $(objdir)/index.o
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/gaussian_filter.o: $(srcdir)/gaussian_filter.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/mod_toms526.o: $(srcdir)/mod_toms526.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/interp_time.o: $(srcdir)/interp_time.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/polygons.o: $(srcdir)/polygons.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/planet.o: $(srcdir)/planet.f90 $(objdir)/geodesic.o
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/geodesic.o: $(srcdir)/geodesic.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/projection_oblimap2.o: $(srcdir)/projection_oblimap2.f90
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/coordinates.o: $(srcdir)/coordinates.f90 $(objdir)/ncio.o $(objdir)/planet.o $(objdir)/geodesic.o \
						 $(objdir)/projection_oblimap2.o $(objdir)/gaussian_filter.o
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/coordinates_mapping.o: $(srcdir)/coordinates_mapping.f90 $(objdir)/coordinates.o
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/coordinates_mapping_conservative.o: $(srcdir)/coordinates_mapping_conservative.f90 \
								 $(objdir)/coordinates.o $(objdir)/coordinates_mapping.o
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/interp2D_conservative.o: $(srcdir)/interp2D_conservative.f90 \
								   $(objdir)/coordinates.o $(objdir)/coordinates_mapping.o
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/subset.o: $(srcdir)/subset.f90 $(objdir)/coordinates.o
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/subset2.o: $(srcdir)/subset2.f90 $(objdir)/coordinates.o
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

$(objdir)/grid_gen.o: $(srcdir)/grid_gen.f90 $(objdir)/coordinates.o $(objdir)/interp2D_conservative.o
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

coord_libs = $(objdir)/ncio.o

coord_obj = $(objdir)/coord_constants.o \
			$(objdir)/index.o \
		    $(objdir)/interp1D.o \
		    $(objdir)/interp2D.o \
		    $(objdir)/loess.o \
		    $(objdir)/gaussian_filter.o \
		    $(objdir)/mod_toms526.o \
		    $(objdir)/interp_time.o \
		    $(objdir)/polygons.o \
		    $(objdir)/planet.o \
		    $(objdir)/geodesic.o \
		    $(objdir)/projection_oblimap2.o \
		    $(objdir)/coordinates.o \
		    $(objdir)/coordinates_mapping.o \
		    $(objdir)/coordinates_mapping_conservative.o \
		    $(objdir)/subset2.o \
		    $(objdir)/grid_gen.o \
		    $(objdir)/interp2D_conservative.o

coord0_obj = $(objdir)/ncio.o \
		    $(objdir)/coord_constants.o \
			$(objdir)/index.o \
		    $(objdir)/interp1D.o \
		    $(objdir)/interp2D.o \
		    $(objdir)/loess.o \
		    $(objdir)/gaussian_filter.o \
		    $(objdir)/mod_toms526.o \
		    $(objdir)/interp_time.o \
		    $(objdir)/polygons.o \
		    $(objdir)/planet.o \
		    $(objdir)/geodesic.o \
		    $(objdir)/projection_oblimap2.o \
		    $(objdir)/coordinates.o \
		    $(objdir)/coordinates_mapping.o \
		    $(objdir)/coordinates_mapping_conservative.o \
		    $(objdir)/subset.o \
		    $(objdir)/grid_gen.o 

# The final library wrapper 
$(objdir)/coord.o: $(srcdir)/coord.f90 $(coord_libs) $(coord_obj)
	$(FC) $(DFLAGS) $(FLAGS) $(SFLAGS) -c -o $@ $<

## Complete programs

# coordinates static library - using subset2
coord-static: $(coord_libs) $(objdir)/coord.o $(coord_obj)
	ar rc libcoordinates.a $(objdir)/coord.o $(coord_obj)
	@echo " "
	@echo "    libcoordinates.a is ready."
	@echo " "

# coordinates shared library - using subset2
coord-shared: $(coord_libs) $(objdir)/coord.o $(coord_obj)
	$(FC) $(DFLAGS) $(FLAGS) -shared -fPIC -o libcoordinates.so $(objdir)/coord.o $(coord_obj) $(LFLAGS)
	@echo " "
	@echo "    libcoordinates.so is ready."
	@echo " "

# coordinates shared library - using subset
coord0-shared: $(coord_libs) $(objdir)/coord.o $(coord0_obj)
	$(FC) $(DFLAGS) $(FLAGS) -shared -fPIC -o libcoordinates0.so $(objdir)/coord.o $(coord0_obj) $(LFLAGS)
	@echo " "
	@echo "    libcoordinates0.so is ready."
	@echo " "

# Program to test interpolations of CCSM3 data
test_ccsm3: coord-static
	$(FC) $(DFLAGS) $(FLAGS) -o test_ccsm3.x $(testdir)/test_ccsm3.f90 libcoordinates.a -L. $(LFLAGS) $(coord_libs)
	@echo " "
	@echo "    test_ccsm3.x is ready."
	@echo " "

test_etopo: coord-static
	$(FC) $(DFLAGS) $(FLAGS) -o test_etopo.x $(testdir)/test_etopo.f90 libcoordinates.a -L. $(LFLAGS) $(coord_libs)
	@echo " "
	@echo "    test_etopo.x is ready."
	@echo " "

test_MAR: coord-static
	$(FC) $(DFLAGS) $(FLAGS) -o test_MAR.x $(testdir)/test_MAR.f90 libcoordinates.a -L. $(LFLAGS) $(coord_libs)
	@echo " "
	@echo "    test_MAR.x is ready."
	@echo " "

test_subset: coord0-shared
	$(FC) $(DFLAGS) $(FLAGS) -o test_subset.x $(testdir)/test_subset.f90 -L. -lcoordinates0 $(LFLAGS) $(coord_libs)
	@echo " "
	@echo "    test_subset.x is ready."
	@echo " "

test_subset2: coord-shared
	$(FC) $(DFLAGS) $(FLAGS) -o test_subset2.x $(testdir)/test_subset2.f90 -L. -lcoordinates $(LFLAGS) $(coord_libs)
	@echo " "
	@echo "    test_subset2.x is ready."
	@echo " "

test_multigrid: coord-shared
	$(FC) $(DFLAGS) $(FLAGS) -o test_multigrid.x $(testdir)/test_multigrid.f90 -L. -lcoordinates $(LFLAGS) $(coord_libs)
	@echo " "
	@echo "    test_multigrid.x is ready."
	@echo " "

test_proj: coord-static
	$(FC) $(DFLAGS) $(FLAGS) -o test_proj.x $(testdir)/test_proj.f90 libcoordinates.a -L. $(LFLAGS) $(coord_libs)
	@echo " "
	@echo "    test_proj.x is ready."
	@echo " "

test_proj_etopo1: coord-static
	$(FC) $(DFLAGS) $(FLAGS) -o test_proj_etopo1.x $(testdir)/test_proj_etopo1.f90 libcoordinates.a -L. $(LFLAGS) $(coord_libs)
	@echo " "
	@echo "    test_proj_etopo1.x is ready."
	@echo " "

test_interp: coord-static
	$(FC) $(DFLAGS) $(FLAGS) -o test_interp.x $(testdir)/test_interp.f90 libcoordinates.a -L. $(LFLAGS) $(coord_libs)
	@echo " "
	@echo "    test_interp.x is ready."
	@echo " "

test_climber: coord-static
	$(FC) $(DFLAGS) $(FLAGS) -o test_climber.x $(testdir)/test_climber.f90 libcoordinates.a -L. $(LFLAGS) $(coord_libs)
	@echo " "
	@echo "    test_climber.x is ready."
	@echo " "

test_ccsm3diff: coord-static
	$(FC) $(DFLAGS) $(FLAGS) -o test_ccsm3_diffusion.x $(testdir)/test_ccsm3_diffusion.f90 libcoordinates.a -L. $(LFLAGS) $(coord_libs)
	@echo " "
	@echo "    test_ccsm3_diffusion.x is ready."
	@echo " "

test_loess: $(objdir)/ncio.o $(objdir)/interp1D.o $(objdir)/index.o $(objdir)/loess.o 
	$(FC) $(DFLAGS) $(FLAGS) -o test_loess.x $^ $(testdir)/test_loess.f90 $(LFLAGS)
	@echo " "
	@echo "    test_loess.x is ready."
	@echo " "

test_nat: extra/tile/t4gen.f extra/tile/t4man.f extra/tile/t4int.f \
	      extra/tile/t4que.f extra/tile/t4int.f90
	$(FC) $(DFLAGS) $(FLAGS) -o test_nat.x $^ $(testdir)/test_nat.f90
	@echo " "
	@echo "    test_nat.x is ready."
	@echo " "

test_nat_new: extra/tile/t4gen.f extra/tile/t4man.f extra/tile/t4int.f \
			  extra/tile/t4que.f extra/tile/t4int.f90
	$(FC) $(DFLAGS) $(FLAGS) -o test_nat_new.x $^ $(testdir)/test_nat_new.f90
	@echo " "
	@echo "    test_nat_new.x is ready."
	@echo " "

# Program to test distance calculations using the geographiclib library
test_geodinverse: $(objdir)/planet.o $(objdir)/geodesic.o
	$(FC) $(DFLAGS) $(FLAGS) -o test_geodinverse.x $^ $(testdir)/test_geodinverse.f90
	@echo " "
	@echo "    test_geodinverse.x is ready."
	@echo " "

# Program to test polygon calculations
test_pointpoly: $(objdir)/polygons.o
	$(FC) $(DFLAGS) $(FLAGS) -o test_pointpoly.x $^ $(testdir)/test_pointpoly.f90 $(LFLAGS)
	@echo " "
	@echo "    test_pointpoly.x is ready."
	@echo " "

# Program to test custom polygon calculations
test_polygon: $(objdir)/polygons.o $(objdir)/index.o
	$(FC) $(DFLAGS) $(FLAGS) -o test_polygon.x $^ $(testdir)/test_polygon.f90 $(LFLAGS)
	@echo " "
	@echo "    test_polygon.x is ready."
	@echo " "

clean:
	rm -r -f *.x *.dSYM $(objdir)/*.o $(objdir)/*.mod *.so 

# cleanall: cleansico cleanrembo cleansicoX
