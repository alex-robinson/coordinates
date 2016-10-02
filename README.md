
# coordinates

The coordinates library provides functions to manage grids
and arbitrary sets of points, including interpolation and mapping
between different coordinate systems. It also contains several helper
modules that are useful generally for the manipulation of sets of
points and grids. 

## coordinates data types

There are 3 main derived data types: `points`, `grid` and `map`. All
calculations are based on the underlying 'points_class':

```fortran
type points_class

    character (len=128) :: name     ! name of this domain (world, GRL, etc)
    character (len=128) :: mtype    ! map type: latlon, cartesian, stereographic, etc
    character (len=128) :: units    ! units of the axes
    type(planet_class)  :: planet   ! which planet are we on?

    ! Projection parameters
    logical :: is_cartesian, is_projection
    logical :: is_lon180
    type(projection_class) :: proj

    ! Points information
    integer :: npts
    real(dp), allocatable, dimension(:)   :: x, y, lon, lat, area
    integer,  allocatable, dimension(:)   :: border
    real(dp) :: xy_conv

end type
```

The points_class stores information about the coordinate system (name, mtype, units),
what type of projection it is (currently only oblique stereographic is possible),
the planet definition (is it Earth?), the locations of the points in the coordinate system
and, if it is a projection, their equivalent longitude and latitude values.

With this information, the coordinates library can achieve any mapping between
one set of locations and another.

A `points` object is initialized with the subroutine `points_init()`.

A `grid` object is equivalent to a `points` object, except the location information
is contained in 2D arrays, and it includes additional grid axis information. It can
be initialized with the subroutine `grid_init()`.

`grid_to_points(grid,pts)` and `points_to_grid(pts,grid)` can be used to convert
between one representation and the other.

The `map` object is used to store information about the nearest neighbors of one
grid or set of points to each point in another grid or set of points. The map_class
looks similar to the points_class and grid_class, except it contains neighbor arrays:

```fortran
type map_class

    ...

    ! Neighbor arrays (will be allocated to size npts,nmax)
    integer :: nmax
    integer,  dimension(:,:), allocatable :: i, quadrant, border
    real(sp), dimension(:,:), allocatable :: dist, weight

end type
```

where `nmax` is the maximum number of neighbors to be stored,
 `i` is the index of a neighbor in the source points/grid, `quadrant` is the location
of that neighbor relative to the target point (I, II, III, IV) and `border` shows whether
the neighbor is located on the border of the source grid. `dist` and `weight` store
the distance and Shepard's inverse distance weighting between the neighbor and the
target point. A `map` object is initialized with the subroutine `map_init()`, which can
take several combinations of arguments:
```fortran
call map_init(pts1,pts2,...)    ! Map from points to points
call map_init(grid1,grid2,...)  ! Map from grid to grid
call map_init(pts1,grid2,...)   ! Map from points to grid
call map_init(grid1,pts2,...)   ! Map from grid to points
```

This mapping information is saved to a NetCDF map file after it is calculated (via `map_write()`).
Initial mapping calculations can take time (typically > 1min for large grids), however it is
only done once. Because the map files are saved with names corresponding to the source
and target grids, the map is then always available. Because the neighbor distances
are saved in the map file, essentially any mapping algorithm can be supported without
the need to recalculate the map.

## A practical example

An underlying theme of the coordinates library is that all necessary information
to represent and manipulate a set of points / grid is contained in the coordinates
objects. Thus to perform operations in a user's program, no temporary variables
are needed. Let's say you are working with 2 gridded domains (a 2° resolution latlon grid
used with CCSM3 and a 20KM resolution regional projection over Greenland) and you need to
map a variable from one to the other. The steps to acheive this are:

1. Initialize each grid.
2. Initialize the map.
3. Map the variable.

```fortran
type(grid_class) :: gCCSM3, gGRL
type(map_class)  :: map1

! Global CCSM3 grid definition
call grid_init(gCCSM3,name="CCSM3-T42",mtype="latlon",units="degrees",
               x0=0.d0,dx=2.d0,nx=180,y0=-90.d0,dx=2.d0,ny=90)

! Oblique stereographic regional grid centered at -40°E and 72°N
call grid_init(gGRL,name="GRL-20KM",mtype="stereographic",units="kilometers", &
                 dx=20.d0,nx=76,dy=20.d0,ny=151, &
                 lambda=-40.d0,phi=72.d0,alpha=7.5d0)

! Generate map
call map_init(map1,gCCSM3,gGRL,max_neighbors=10)

! Map variable using 'quadrant' method
call map_field(map1,"var_name",var_ccsm3,var_grl,mask_map,method="quadrant")

```

In the above example, the 2D arrays `var_ccsm3` and `var_grl` should already be defined in the program.
The only additional variables that need to be defined to complete a mapping are:
`gCCSM3`, `gGRL` and `map1`. Furthermore, defining the mapping between the two grids is trivial,
since all the necessary grid information is accessible in the grid objects.
The user simply specifies how many nearest neighbors to check for. In this way,
the additional code needed to incorporate mapping via the coordinates
library is as minimal as possible. Note: the optional 2D array `mask_map` returned can also be
used to track which points in the target array were overwritten by the mapping.

The mapping is only performed once and saved to a NetCDF file, by default in the `maps/` subdirectory.
If later another mapping between the same grids with the same number of neighbors is performed,
the map is loaded into memory from the file.

For a more complete example based on similar definitions, see the test program **test_ccsm3.f90**,
which includes mapping tests to various regional domains and back to the global grid, as defined
in the paper about oblique stereographic projections by
Reerink et al. (2010, www.geosci-model-dev.net/3/13/2010/).
The example can be compiled by calling `make ccsm3` (see below in **Makefile** for details).

## Algorithms

### Oblique stereographic projections

The coordinates library has the ability to manipulate projected grids - Cartesian
representations of curvilinear coordinates. Currently only oblique stereographic
(including polar stereographic) projections are supported.
The algorithms implemented by Reerink et al. (2010, see www.geosci-model-dev.net/3/13/2010/)
are used to transform coordinates
between latlon and stereographic coordinates. These algorithms were ported from
the OBLIMAP2 package.

To define a given projected coordinate system, three parameters are needed: the
latlon center of the grid (phi, lambda) and the angle alpha, which represents the
level at which the Cartesian plane slices through the spherical surface. Typical
values are (as in Reerink et al, 2010):

- Greenland:  lambda=320°, phi= 72°, alpha= 7.5°
- Antarctica: lambda=  0°, phi=-90°, alpha=19.0°

Typically the smaller the extent of the Cartesian domain, the smaller the value of
alpha that should be used.

### Distance calculations

Unless the sets of points are located on the same Cartesian/projected coordinate system,
all distance calculations are performed using geodesic calculations. The coordinates
library incorporates geographiclib to calculate distances
between points. The principal advantages of these geodesic algorithms
over previous ones (e.g., Vincenty, 1975) are
- accurate to round off for |*f*| < 1/50;
- the solution of the inverse problem is always found;
- differential and integral properties of geodesics are computed.

The algorithms were derived by Charles F. Karney (doi:10.1007/s00190-012-0578-z).
More information can be found here:
http://geographiclib.sourceforge.net/html/Fortran/

geographiclib has been converted to Fortran90 and is included in this library.

### Interpolation

Documentation needed.

## Getting started

### Installation

The coordinates library depends on the Fortran **NetCDF** library (http://www.unidata.ucar.edu/software/netcdf/),
which must already be installed on the user's system. To interface with NetCDF,
coordinates makes use of the library **NCIO** (https://github.com/alex-robinson/ncio),
which has been incorporated into the coordinates library - therefore no additional
installation is necessary.

With **NetCDF** already installed, coordinates should work on any system with
a Fortran compiler. Simply download the latest stable version to begin using it.

### Makefile

The main coordinates module depends on other internal modules
(NCIO, planet, geodesic and projection_oblimap2) that must
be compiled together with the coordinates module itself. This is accomplished
here via a Makefile. Once the coordinates.o object file is compiled along
with these dependencies, it is possible to use it in a Fortran program by
adding the statement `use coordinates`.

A test program is available to make sure the library is working properly.
To be able to compile it, first make sure that the paths to the NetCDF
include and lib directories are properly specified for your system in the Makefile:

```Makefile
netcdf_inc = /opt/local/include
netcdf_lib = /opt/local/lib
```

Note: by default the compiler is `gfortran`. If you will use `ifort`, make sure to
change paths `netcdf_inc_ifort` and `netcdf_lib_ifort`.

To compile the test program call: `make ccsm3` (or `make ccsm3 ifort=1` if using ifort).  
If it compiles without error, run the test program: `./test_ccsm3.x`.

A matrix with the following results for mapping to a regional Antarctic projection
and back to the global CCSM3 grid should be output to the screen at the end of the run:

    Ts     240.2     277.6     259.6      0.0114      0.0360      0.0304
    MB       0.0       0.6       0.2      0.0003      0.0008      0.0484
    Hs    -109.9    3628.4    1337.6      0.9612      2.2264      0.0257

If your results match those above, the program is working correctly.

Fore more details, the individual source code compilation rules along with
dependencies can be found in the Makefile:

```Makefile
## Individual libraries or modules ##
$(objdir)/ncio.o: ncio.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp1D.o: interp1D.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp2D.o: interp2D.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/interp_time.o: interp_time.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/planet.o: planet.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/geodesic.o: geodesic.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/projection_oblimap2.o: projection_oblimap2.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/coordinates.o: coordinates.f90 $(objdir)/ncio.o $(objdir)/planet.o $(objdir)/geodesic.o \
						 $(objdir)/projection_oblimap2.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

$(objdir)/subset.o: subset.f90 $(objdir)/coordinates.o
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<
```
