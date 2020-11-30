.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

# PATH options
srcdir  = src
objdir  = libcoordinates/include
bindir  = libcoordinates/bin
testdir = tests

# Command-line options at make call
debug  ?= 0  
openmp ?= 1

## COMPILER CONFIGURATION ##
# (should be loaded from config directory)

<COMPILER_CONFIGURATION>

# Determine whether to use normal flags (default), 
# debugging flags (debug=1) or profiling flags (debug=2)
DFLAGS   = $(DFLAGS_NODEBUG)
ifeq ($(debug), 1)
	DFLAGS   = $(DFLAGS_DEBUG)
endif
ifeq ($(debug), 2)
	DFLAGS   = $(DFLAGS_PROFILE)
endif

# Determine whether to use openmp flags 
FFLAGS = $(FFLAGS_DEFAULT)
ifeq ($(openmp), 1)
	FFLAGS = $(FFLAGS_OPENMP)
endif

# Compilation shortcut
LDFLAGS = $(FFLAGS) $(DFLAGS)

# Specify library related inc/lib commands 
COORDROOT = ${CURDIR}/libcoordinates
INC_COORD = -I${COORDROOT}/include
LIB_COORD = -L${COORDROOT}/include -lcoordinates

###############################################
##							
## List of rules and source files
##
###############################################

include config/Makefile_coordinates.mk

###############################################
##
## Compilation of complete programs
##
###############################################

## Complete programs

# coordinates static library - using subset2
coord-static: $(coord_libs) $(objdir)/coord.o $(coord_obj)
	ar rc $(objdir)/libcoordinates.a $(objdir)/coord.o $(coord_obj)
	@echo " "
	@echo "    $(objdir)/libcoordinates.a is ready."
	@echo " "

# coordinates shared library - using subset2
coord-shared: $(coord_libs) $(objdir)/coord.o $(coord_obj)
	$(FC) $(LDFLAGS) -shared -fPIC -o $(objdir)/libcoordinates.so $(objdir)/coord.o $(coord_obj) $(LFLAGS)
	@echo " "
	@echo "    $(objdir)/libcoordinates.so is ready."
	@echo " "

# coordinates shared library - using subset
coord0-shared: $(coord_libs) $(objdir)/coord.o $(coord0_obj)
	$(FC) $(LDFLAGS) -shared -fPIC -o $(objdir)/libcoordinates0.so $(objdir)/coord.o $(coord0_obj) $(LFLAGS)
	@echo " "
	@echo "    $(objdir)/libcoordinates0.so is ready."
	@echo " "

# Program to test interpolations of CCSM3 data
test_ccsm3: coord-static
	$(FC) $(LDFLAGS) $(INC_COORD) -o $(bindir)/test_ccsm3.x $(testdir)/test_ccsm3.f90 $(LIB_COORD) $(coord_libs) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_ccsm3.x is ready."
	@echo " "

test_etopo: coord-static
	$(FC) $(LDFLAGS) -o $(bindir)/test_etopo.x $(testdir)/test_etopo.f90 $(objdir)/libcoordinates.a -L$(objdir) $(coord_libs) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_etopo.x is ready."
	@echo " "

test_MAR: coord-static
	$(FC) $(LDFLAGS) -o $(bindir)/test_MAR.x $(testdir)/test_MAR.f90 $(objdir)/libcoordinates.a -L$(objdir) $(coord_libs) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_MAR.x is ready."
	@echo " "

test_subset: coord0-shared
	$(FC) $(LDFLAGS) -o $(bindir)/test_subset.x $(testdir)/test_subset.f90 -L. -lcoordinates0 $(coord_libs) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_subset.x is ready."
	@echo " "

test_subset2: coord-static
	$(FC) $(LDFLAGS) -o $(bindir)/test_subset2.x $(testdir)/test_subset2.f90 -L. -lcoordinates $(coord_libs) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_subset2.x is ready."
	@echo " "

test_multigrid: coord-static
	$(FC) $(LDFLAGS) -o $(bindir)/test_multigrid.x $(testdir)/test_multigrid.f90 -L. -lcoordinates $(coord_libs) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_multigrid.x is ready."
	@echo " "

test_proj: coord-static
	$(FC) $(LDFLAGS) -o $(bindir)/test_proj.x $(testdir)/test_proj.f90 $(objdir)/libcoordinates.a -L$(objdir) $(coord_libs) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_proj.x is ready."
	@echo " "

test_proj_etopo1: coord-static
	$(FC) $(LDFLAGS) -o $(bindir)/test_proj_etopo1.x $(testdir)/test_proj_etopo1.f90 $(objdir)/libcoordinates.a -L$(objdir) $(coord_libs) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_proj_etopo1.x is ready."
	@echo " "

test_interp: coord-static
	$(FC) $(LDFLAGS) -o $(bindir)/test_interp.x $(testdir)/test_interp.f90 $(objdir)/libcoordinates.a -L$(objdir) $(coord_libs) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_interp.x is ready."
	@echo " "

test_climber: coord-static
	$(FC) $(LDFLAGS) -o $(bindir)/test_climber.x $(testdir)/test_climber.f90 $(objdir)/libcoordinates.a -L$(objdir) $(coord_libs) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_climber.x is ready."
	@echo " "

test_ccsm3diff: coord-static
	$(FC) $(LDFLAGS) -o $(bindir)/test_ccsm3_diffusion.x $(testdir)/test_ccsm3_diffusion.f90 $(objdir)/libcoordinates.a -L$(objdir) $(coord_libs) $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_ccsm3_diffusion.x is ready."
	@echo " "

test_loess: $(objdir)/ncio.o $(objdir)/interp1D.o $(objdir)/index.o $(objdir)/loess.o 
	$(FC) $(LDFLAGS) -o $(bindir)/test_loess.x $^ $(testdir)/test_loess.f90 $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_loess.x is ready."
	@echo " "

test_nat: extra/tile/t4gen.f extra/tile/t4man.f extra/tile/t4int.f \
	      extra/tile/t4que.f extra/tile/t4int.f90
	$(FC) $(LDFLAGS) -o $(bindir)/test_nat.x $^ $(testdir)/test_nat.f90
	@echo " "
	@echo "    $(bindir)/test_nat.x is ready."
	@echo " "

test_nat_new: extra/tile/t4gen.f extra/tile/t4man.f extra/tile/t4int.f \
			  extra/tile/t4que.f extra/tile/t4int.f90
	$(FC) $(LDFLAGS) -o $(bindir)/test_nat_new.x $^ $(testdir)/test_nat_new.f90
	@echo " "
	@echo "    $(bindir)/test_nat_new.x is ready."
	@echo " "

# Program to test distance calculations using the geographiclib library
test_geodinverse: $(objdir)/planet.o $(objdir)/geodesic.o
	$(FC) $(LDFLAGS) -o $(bindir)/test_geodinverse.x $^ $(testdir)/test_geodinverse.f90
	@echo " "
	@echo "    $(bindir)/test_geodinverse.x is ready."
	@echo " "

# Program to test polygon calculations
test_pointpoly: $(objdir)/polygons.o
	$(FC) $(LDFLAGS) -o $(bindir)/test_pointpoly.x $^ $(testdir)/test_pointpoly.f90 $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_pointpoly.x is ready."
	@echo " "

# Program to test custom polygon calculations
test_polygon: $(objdir)/polygons.o $(objdir)/index.o
	$(FC) $(LDFLAGS) -o $(bindir)/test_polygon.x $^ $(testdir)/test_polygon.f90 $(LFLAGS)
	@echo " "
	@echo "    $(bindir)/test_polygon.x is ready."
	@echo " "


.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make test_ccsm3 : compiles test_ccsm3.x"
	@echo " make clean      : cleans object files"
	@echo ""

clean:
	rm -f $(bindir)/*.x
	rm -f  *.x gmon.out $(objdir)/*.o $(objdir)/*.mod $(objdir)/*.a $(objdir)/*.so
	rm -rf *.x.dSYM
