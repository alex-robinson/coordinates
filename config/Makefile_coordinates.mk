## Individual libraries or modules ##
$(objdir)/ncio.o: $(srcdir)/ncio.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/coord_constants.o: $(srcdir)/coord_constants.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/index.o: $(srcdir)/index.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/interp1D.o: $(srcdir)/interp1D.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/interp2D.o: $(srcdir)/interp2D.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/loess.o: $(srcdir)/loess.f90 $(objdir)/interp1D.o $(objdir)/index.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/gaussian_filter.o: $(srcdir)/gaussian_filter.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/mod_toms526.o: $(srcdir)/mod_toms526.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/interp_time.o: $(srcdir)/interp_time.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/polygons.o: $(srcdir)/polygons.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/planet.o: $(srcdir)/planet.f90 $(objdir)/geodesic.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/geodesic.o: $(srcdir)/geodesic.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/projection_oblimap2.o: $(srcdir)/projection_oblimap2.f90
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/coordinates.o: $(srcdir)/coordinates.f90 $(objdir)/ncio.o $(objdir)/planet.o $(objdir)/geodesic.o \
						 $(objdir)/projection_oblimap2.o $(objdir)/gaussian_filter.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/coordinates_mapping.o: $(srcdir)/coordinates_mapping.f90 $(objdir)/coordinates.o
	$(FC) $(LDFLAGS) -c -o $@ $<


$(objdir)/grid_to_cdo.o: $(srcdir)/grid_to_cdo.f90 $(objdir)/coordinates.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/coordinates_mapping_conservative.o: $(srcdir)/coordinates_mapping_conservative.f90 \
								 $(objdir)/coordinates.o $(objdir)/coordinates_mapping.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/coordinates_mapping_scrip.o: $(srcdir)/coordinates_mapping_scrip.f90 \
								 $(objdir)/coordinates.o $(objdir)/grid_to_cdo.o \
								 $(objdir)/ncio.o $(objdir)/index.o $(objdir)/interp2D.o $(objdir)/gaussian_filter.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/interp2D_conservative.o: $(srcdir)/interp2D_conservative.f90 \
								   $(objdir)/coordinates.o $(objdir)/coordinates_mapping.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/subset.o: $(srcdir)/subset.f90 $(objdir)/coordinates.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/subset2.o: $(srcdir)/subset2.f90 $(objdir)/coordinates.o
	$(FC) $(LDFLAGS) -c -o $@ $<

$(objdir)/grid_gen.o: $(srcdir)/grid_gen.f90 $(objdir)/coordinates.o $(objdir)/interp2D_conservative.o
	$(FC) $(LDFLAGS) -c -o $@ $<

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
		    $(objdir)/coordinates_mapping_scrip.o \
		    $(objdir)/grid_to_cdo.o \
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
	$(FC) $(LDFLAGS) -c -o $@ $<
