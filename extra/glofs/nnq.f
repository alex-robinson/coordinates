c Downloaded from: 
c     http://www.nco.ncep.noaa.gov/pmb/codes/nwprod/glofs.v1.0.0/sorc/nos_interp_sfcmarobs.fd/
c
c------------------------------------------------------------------------
c
c	nn2d_setup - Performs all setup procedures for natural neighbour 
c		     interpolation routine nn2D.
c
c	Input:
c	        np			number of nodes
c	        nt_max			maximum number of triangles
c	        nh_max			maximum number of triangles on convex 
c					hull (max size of array hulltriangles)
c		np_max			maximum number of nodes 
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                        set to ~20 in calling program)
c		nmax			maximum sum of the number of neighbours 
c					per node (should be set to 
c					3*nt_max + np_max in calling program)
c		points(2,np)		array of node co-ordinates
c	        dmode			Delaunay calculation mode (integer)
c	        nmode			NN setup mode
c	        clockwise		logical for the vertices input order  
c		data(np)		data values at each node
c		nnn			Integer work array used by build_nv
c		nnlist			Integer work array used by build_nv
c		ntrilist		Integer work array used by build_nv
c		nohalt_hull		determines error response in routine
c					calculate_hulltriangles
c		eps 			tolerance parameter used by delaun
c					(see delaun for details)
c		vis_tlist		Integer work array used by delaun
c		vis_elist		Integer work array used by delaun
c		add_tlist		Integer work array used by delaun
c		nv_max			size of delaun work arrays vis_tlist
c					vis_elist, & add_tlist (passed to 
c					delaun for error checking)
c
c	Output:
c	        nt			number of triangles
c	        vertices(3,nt)		array of triangle vertices 
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c	        nh			number of triangles with an edge
c					on the convex hull
c		hulltriangles(nh)	array of triangles with an edge 
c					on the convex hull
c		loc			an initial guess triangle for point 
c					location routine `triloc' used by nn2D
c					(set somewhere near the centre)
c
c	Operation modes:
c
c		 The setup routine will perform different tasks depending
c		 on the input parameters dmode and nmode (see table below).
c		 Depending on the modes used some work arrays may be set 
c		 to size 1 to save memory. The "Memory Savings" column in the
c		 table below shows the dimension statement that may
c		 be used in the calling program if the routine is ONLY EVER 
c		 CALLED IN THE CORRESPONDING MODE. 
c
c		 PARAMETERS	ACTION	 		MEMORY SAVINGS
c
c		 nmode = 1	Delaunay only 		real*8 centres(3,1)
c							integer hulltriangles(1)
c		 nmode = 0	Delaunay + nn setup	
c		 nmode = -1	nn setup only		
c
c		 dmode > 0	Delaunay read in from   integer vis_tlist(1)
c				logical unit dmode.	integer vis_elist(1) 
c							integer add_tlist(1) 
c		 dmode = 0      Qhull used		Same as dmode > 0.
c		 dmode = -1	Delaun + X-sort         integer nnn(1)
c							integer nnlist(1)
c							integer ntrilist(1)
c		 dmode = -2	Delaun + no sort        Same as dmode=-1
c
c		 dmode = 0 & nmode=1			integer neighbour(3,1)
c		 
c		 A call with nmode = -1 can only be made after a call 
c		 with nmode = 1.
c
c	Comments:
c
c		 If the arrays are used then they should be dimensioned
c		 in the calling program in the following way:
c
c		 real*8		points(2,np_max)
c		 real*8		centres(3,nt_max)
c		 integer	vertices(3,nt_max)
c		 integer	neighbour(3,nt_max)
c		 integer	hulltriangles(nh_max)
c		 
c		 Note: nh_max can usually be dimensioned much less than nt_max
c		 because hulltriangles, stores in a compact form, all 
c		 triangles with an edge on the convex hull. Except for
c		 very irregular point distributions nh << nt. If nh is 
c		 determined to be > nh_max then an error is reported and the
c		 program is halted (unless nohalt parameter is set to 1). 
c		 The array hulltriangles is only used by nn2Do see routine 
c		 calculate_hulltriangles. If nh_max = 1 then hulltriangles 
c		 is not calculated.
c
c		 The initial guess triangle 'loc' is set in nn_setup but
c		 at each call it will be set to the triangle found during
c		 the previous call to nn2D. The user may modify its value
c		 if the input point (x,y) is known to be in, or near, a
c		 particular triangle.
c
c		 If dmode > 0 the the deluanay tessellation is read in from
c		 logical unit `dmode' instead of being calculated internally 
c		 This can be useful if qhullf fails because
c		 of precision errors. The Deluanay may be determined
c		 externally to this program using a double precision version
c		 or another algorithm, e.g. Fortune's sweepline method.
c
c		 If 50 > dmode > 0 then:
c		 It is assumed that the read in format has one triangle per
c		 line represented as a triplet of nodes numbered from ZERO, 
c		 which is the standard output format of codes qhull 
c		 (quickhull method) and voronoi (sweepline method).
c		 If clockwise = .true. (.false.) then the vertices are assumed 
c		 to be in clockwise (anti-clockwise) order. Note program
c		 qhull outputs vertices in anti-clockwise order while 
c		 voronoi in clockwise order. The internal format is 
c		 anti-clockwise and nodes numbered from ONE.
c
c		 If dmode => 50 then:
c		 It is assumed that the read in format has one triangle per
c		 line represented as a triplet of nodes numbered from ONE,
c		 which is the output format of program del (using delaun). 
c
c		 Three other work arrays are produced as a `by product'
c		 of the routine build_nv which calculates the neighbour
c		 array. These must be dimensioned in the calling program in 
c		 the following way (unless delaun is used for calculating the
c		 Delaunay because it already determines the neighbour array)
c
c		 integer nnn(np_max+1)  : number of neighbours per node
c		 integer nnlist(nmax)   : natural neighbours per node
c		 integer ntrilist(nmax) : triangles attached to each node 
c
c		 The value of nmax should be set to (3*nt_max + np_max)
c		 in the calling program.
c
c		 Each of these are useful lists that describe features of
c		 the Voronoi diagram. Both nnlist and ntrilist are stored in
c		 a compact format to avoid zeros. They are only used 
c		 in the setup routine and the memory may be freed once
c		 initialization is completed.
c		 
c
c		 Calls are made to: qhullf, ccentres, build_nv and 
c				    calculate_hulltriangles, delaun.
c		 
c					M. Sambridge, RSES, April 1994.
c  					        (Last modified 10/4/96)
c
c------------------------------------------------------------------------
c
	subroutine nn2dsetup
     &             (np,nt_max,nh_max,np_max,nnpn_max,nmax,
     &              points,dmode,nmode,clockwise,data,nt,vertices,
     &              centres,neighbour,nh,hulltriangles,nohalt_hull,
     &              loc,nnn,nnlist,ntrilist,
     &              eps,nv_max,vis_tlist,vis_elist,add_tlist,
     &		    lt_work,ln_work)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		eps
	integer		vertices(3,*)
	integer		neighbour(3,*)
        integer		hulltriangles(*)
	integer		nnn(*)
	integer		nnlist(*)
	integer		ntrilist(*)
        integer         vis_tlist(*)
        integer         vis_elist(*)
        integer         add_tlist(*)
        integer         dmode,nmode
	logical*1	lt_work(*)
	logical*1	ln_work(*)
	logical		nnwrite
	logical		clockwise
	logical		ldummy

	common/nnswitches/lud,nnwrite


        if(nmode.eq.1.or.nmode.eq.0)then

           if(dmode.eq.0)then
c                                       calculate Delaunay using qhull 
 
              call qhullf(np,2,2,nt_max,0,points,nt,vertices)

           else if(dmode.eq.-1.or.dmode.eq.-2)then

c					sort the points in ascending x order
	      if(dmode.eq.-1)then

c		  write(*,*)' X sort in progress'

	          call hpsort_d(np,1,points)

c		  write(*,*)' X sort done'

	      end if



c                                       calculate Delaunay using delaun 
 
              call delaun (points,np,neighbour,vertices,nt,nt_max,
     &                     vis_tlist,vis_elist,add_tlist,eps,nv_max,
     &                     0,ldummy,0,0,0)



	   else
c                                       read in Delaunay vertices

              nt = 0
              i1 = 1
              i2 = 2
              if(clockwise)then
                 i1 = 2
                 i2 = 1
              end if
              read(dmode,*)
  1           read(dmode,*,end=3,err=2)
     &        vertices(i1,nt+1),vertices(i2,nt+1),vertices(3,nt+1)
              nt = nt + 1
              if(nt.ge.nt_max)then
                 write(*,*) 'Error in nn_setup: too many triangles'
                 write(*,*) 'Remedy: increase size of parameter nt_max'
                 write(*,*) '        in calling program.'
                 stop 
              end if
              go to 1
  2           write(*,*)
     &        'Error in nn_setup: read error in Delaunay input file'
              stop
  3           continue
     
           end if
 
c					adjust array vertices to
c					range from nodes 1 to np
c
	   if(dmode.ge.0.and.dmode.lt.50)then
	      do 5 i = 1,nt
	         vertices(1,i) = vertices(1,i) + 1
	         vertices(2,i) = vertices(2,i) + 1
	         vertices(3,i) = vertices(3,i) + 1
 5            continue
	   end if

	end if
c
c					Perform set up for nn interpolation
c
        if(nmode.eq.0.or.nmode.eq.-1)then

c					set initial guess for 
c					triangle location procedure
           loc = nt/2

c                                       Calculate Circumcentres




           call ccentres(points,vertices,nt,centres)

c                                       Build neighbour matrix
c					(if not already built)

           if(dmode.ge.0)then
              call build_nv
     &        (np,vertices,nt,np_max,nmax,
     &         neighbour,nnn,nnlist,ntrilist)
           end if

c					calculate hulltriangles

           if(nh_max.gt.1) call calculate_hulltriangles
     &     (neighbour,nt,nh_max,hulltriangles,nh,nohalt_hull)

c					initialize logical work arrays 
           do i=1,nt
              lt_work(i) = .false.
           end do
           do i=1,np
              ln_work(i) = .false.
           end do

	end if

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2D - calculates natural neighbour interpolation at point x,y
c
c	Input:
c		x,y			co-ordinates of input point
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c		hulltriangles(nh)	array of triangles with an edge 
c					on the convex hull
c	        nh			number of triangles with an edge
c					on the convex hull
c		loc			an initial guess triangle for point 
c					location routine `triloc'
c		data			data values at each node
c		nnext			logical if true then extend nn outside
c					of convex hull
c	        int_method		determines the type of interpolation
c					method. 0= Watson's method
c						1= recursive method
c						2= recursive method + df(2)
c						3= Linear interpolation + df(2)
c						4= closed formula method 
c						5= same as 4 + df(2)
c						6= same as 4 + df(2) + ddf(4)
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                        set to ~50 in calling program)
c		work...                 work arrays of the size of 
c					nnpn_max.
c		lt_work			logical work array of size nt_max 
c		ln_work			logical work array of size np_max 
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		v			area of voronoi cell about (x,y)
c		out			logical determining if (x,y) is i
c					inside the convex hull of points
c
c	Comments:
c		 The area of the voronoi cell about (x,y) can be used 
c		 as an inverse measure of the quality of the interpolation.
c		 If (x,y) is outside of the hull the voronoi cell is
c		 unbounded and v is the sum of the natural neighbour
c		 co-ordinates calculated with Watson's method.
c
c		 Note: plot_tc and plot_c are user supplied routines
c		 that can be dummy routines, or used with a plot
c		 program for plotting circum-triangles and circum-circles 
c		 found during the calculation. The file nnplot.f contains
c		 a set of dummy routines to be compiled with nn.f. 
c		 Alternatively all references to plot_tc and plot_c may
c		 commented out.
c
c		 A recursive method is available, as an alternative to
c		 Watson's method, for calculation of natural neighbour
c		 co-ordinates. The action of nn_setup is the same for 
c		 either method although the parameter nnpn_max is also 
c		 used to size arrays used in the recursive mode.
c
c		 Calls are made to: triloc, nn2Dr, nn2Drd, nn2Do or nn2Di
c				    nn2Df, nn2Dfd, nn2Dfdd.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine nn2D
     &  (x,y,points,vertices,neighbour,nnext,int_method,centres,
     &   hulltriangles,nh,loc,data,nnpn_max,work_r1,work_r2,
     &   work_d1,work_i1,work_i2,work_i3,work_d2,work_d3,work_d4,
     &   work_d5,work_d6,work_d7,lt_work,ln_work,out,f,df,ddf,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		xd(2),x,y
	real*8		data(*)
	real*8		f,df(2),ddf(3),v
	integer		vertices(3,*)
	integer		neighbour(3,*)
        integer		hulltriangles(*)
	integer		int_method
	logical		out
	logical		nnext
	logical		nnwrite

        real*8          work_d1(2,nnpn_max)
        real*8          work_d2(2,nnpn_max)
        real*8          work_d3(2,nnpn_max)
        real*8          work_d4(2,nnpn_max)
        real*8          work_d5(2,nnpn_max)
        real*8          work_d6(2,nnpn_max)
        real*8          work_d7(2,nnpn_max)
        real            work_r1(nnpn_max,2)
        real            work_r2(nnpn_max)
	integer		work_i1(nnpn_max)
        integer         work_i2(2,nnpn_max)
        integer         work_i3(nnpn_max)
	logical*1	lt_work(*)
	logical*1	ln_work(*)

	common/nnswitches/lud,nnwrite
 
        xd(1) = x
        xd(2) = y
c					Use triangle walking routine
c					to locate triangle of input
c					point

        call triloc(x,y,points,vertices,neighbour,loc,out)
c
c					If point is outside convex hull
c					call routine nn2Do
	if(out)then
 	   if(nnext.and.int_method.eq.0)then
	        if(nnwrite)write(*,*)' calling nn2Do'
                call nn2Do
     &               (x,y,points,vertices,neighbour,
     &                centres,hulltriangles,nh,data,f,df,v)
           else
               f = 0.d0
               df(1) = 0.d0
               df(2) = 0.d0
               v = 0.d0
	   end if

	else

c					use modified Watson's method 
	   if(int_method.eq.0)then

 	      if(nnwrite)write(*,100)

              call nn2Di
     &             (x,y,points,vertices,neighbour,
     &              centres,loc,data,f,df,v)

           else if(int_method.eq.1)then

 	      if(nnwrite)write(*,200)
c					use recursive method
              call nn2Dr
     &             (x,y,points,vertices,neighbour,centres,
     &              loc,data,nnpn_max,work_r1,work_r2,work_d1,
     &              work_i1,work_i2,work_i3,f,df,v)


           else if(int_method.eq.2)then

 	      if(nnwrite)write(*,300)
c					use recursive method and 
c					calculate derivatives
              call nn2Drd
     &             (x,y,points,vertices,neighbour,centres,
     &              loc,data,nnpn_max,work_r1,work_r2,work_d1,
     &              work_i1,work_i2,work_i3,f,df,v)

           else if(int_method.eq.3)then

 	      if(nnwrite)write(*,400)
c					use linear interpolation in triangles
	      call nn2DL
     &             (x,y,points,vertices,loc,data,f,df,v)

           else if(int_method.eq.4)then

 	      if(nnwrite)write(*,500)
c					use closed formula method
c					no derivatives
c
              call nn2Df
     &             (xd,points,vertices,neighbour,centres,
     &              loc,data,nnpn_max,work_d1,work_d2,work_i1,
     &              work_i2,work_r2,work_i3,lt_work,ln_work,f,v)

              df(1) = 0.d0
              df(2) = 0.d0

           else if(int_method.eq.5)then

 	      if(nnwrite)write(*,600)
c					use closed formula method
c					with 1st derivatives
c
              call nn2Dfd
     &             (xd,points,vertices,neighbour,centres,
     &              loc,data,nnpn_max,work_d1,work_d2,work_i1,
     &              work_i2,work_r2,work_i3,work_d3,work_d4,
     &              lt_work,ln_work,f,df,v)

           else if(int_method.eq.6)then

 	      if(nnwrite)write(*,600)
c					use closed formula method
c					with 1st and 2nd derivatives
c
              call nn2Dfdd
     &             (xd,points,vertices,neighbour,centres,
     &              loc,data,nnpn_max,work_d1,work_d2,work_i1,
     &              work_i2,work_r2,work_i3,work_d3,work_d4,
     &	            work_d5,work_d6,work_d7,
     &              lt_work,ln_work,f,df,ddf,v)

           else

              write(*,*)' '
              write(*,*)
     &        ' Error in nn2D: interpolation method not defined'
              write(*,*)'                int_method set to ',int_method
              write(*,*)' '
              write(*,*)' Valid options are:'
              write(*,*)'   0=NN-Watson'
              write(*,*)'   1=NN-recursive'
              write(*,*)'   2=NN-recursive+derivatives'
              write(*,*)'   3=Linear interpolation'
              write(*,*)'   4=NN-closed formula'
              write(*,*)'   5=NN-closed formula + df(2)'
              write(*,*)'   6=NN-closed formula + df(2) + ddf(3)'
              stop


           end if

	end if

c	if(nnwrite)then
c	   write(*,*)' Overlap voronoi area = ',v
c	   write(*,*)' Normalized value of f   =',f
c	   write(*,*)' Normalized value of dfx =',df(1)
c	   write(*,*)' Normalized value of dfy =',df(2)
c	end if

 100 	format(' calling nn2Di:   Watson method')
 200 	format(' calling nn2Dr:   recursive method')
 300 	format(' calling nn2Drd:  recursive method with derivatives')
 400 	format(' calling nn2DL:   Linear interpolation')
 500 	format(' calling nn2Df:   Closed formulae method')
 600 	format(' calling nn2Dfd:   Closed formulae method',
     &         'with first derivatives')
 700 	format(' calling nn2Dfdd:   Closed formulae method',
     &         'with first and second derivatives')

	return
	end
c
c------------------------------------------------------------------------
c
c       nn2DL - calculates linear triangle interpolation at point x,y
c               (This routine is an alternative to nn2D)
c
c       Input:
c               x,y                     co-ordinates of input points
c               points(2,np)            array of node points
c               vertices(3,nt)          array of triangle vertices
c               loc                     an initial guess triangle for point
c                                       location routine `triloc'
c               data                    data values at each node
c
c       Output:
c               f                       interpolated data value
c               df(1)                   interpolated derivative df/dx
c               df(2)                   interpolated derivative df/dy
c               v                       area of voronoi cell about (x,y)
c
c       Comments:
c
c                Performs simple linear interpolation in triangles
c
c                                       M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
        Subroutine nn2DL
     &  (x,y,points,vertices,loc,data,f,df,v)

        real*8          points(2,*)
        real*8          x,y
        real*8          x1,y1,x2,y2,x3,y3,A,a1,a2,a3
        real*8          dx12,dx13,dy12,dy13,df12,df13
        real*8          data(*)
        real*8          f,df(2),v
        integer         vertices(3,*)


c                                       calculate area of triangle loc
        i = vertices(1,loc)
        j = vertices(2,loc)
        k = vertices(3,loc)
        x1 = points(1,i)
        y1 = points(2,i)
        x2 = points(1,j)
        y2 = points(2,j)
        x3 = points(1,k)
        y3 = points(2,k)
        A = ((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
        if(A.lt.0.0)
     &  write(6,*)' Error in nn2DL: negative area for T=',loc

c                                       calculate f and derivatives
        df12 = data(i) - data(j)
        df13 = data(i) - data(k)
        dx12 = x1 - x2
        dx13 = x1 - x3
        dy12 = y1 - y2
        dy13 = y1 - y3
        a1 = (df12*dy13 - df13*dy12)/A
        a2 = (df13*dx12 - df12*dx13)/A
        a3 = data(i) - a1*x1 - a2*y1
        f = a1*x + a2*y + a3
        df(1) = a1
        df(2) = a2
        v = 0.d0
c       write(*,*)' New value of f =',f
c       write(*,*)' New value of dfx =',df(1)
c       write(*,*)' New value of dfy =',df(2)

        return
        end
c
c
c------------------------------------------------------------------------
c
c	build_nv - Builds neighbour array for Delaunay triangulation in 2-D.
c
c	Input:	
c	        vertices(3,nt)		array of triangle vertices	
c	        nt			number of triangles
c		np_max			maximum number of nodes
c		nmax			maximum total number of neighbours 
c					per node (should be set to 
c					3*nt_max + np_max in calling program)
c
c	Output:
c		neighbour(3,nt)		array of neighbouring triangles
c
c	Comments:
c		 Assumes input list of vertices in anticlockwise sequence
c		 and produces an anticlockwise list of neighbour triangles.
c		 The value of neighbour(i,j) is the index of the neighbouring
c		 triangle opposite node i in triangle j.
c
c		 Three temporary work arrays are used and must be dimensioned
c		 in the calling program in the following way:
c
c		 integer nnn(np_max+1)  : number of neighbours per node
c		 integer nnlist(nmax)   : natural neighbours per node
c		 integer ntrilist(nmax) : triangles attached to node 
c
c		 The value of nmax should be set to (3*nt_max + np_max)
c		 in the calling program.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1994.
c					(using ideas by J.Braun)
c
c------------------------------------------------------------------------
c
	Subroutine build_nv
     &             (np,vertices,nt,np_max,nmax,
     &              neighbour,nnn,nnlist,ntrilist)
c
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		nnn(*)
	integer		nnlist(*)
	integer		ntrilist(*)
	logical		nnwrite

	common/nnswitches/lud,nnwrite

	if(nnwrite)write(*,*)' Building neighbour v ...'
c
c					initialize neighbour list
	do 5 i = 1,3
	   do 5 j = 1,nt
	      neighbour(i,j) = 0
 5      continue
c					initialize work arrays
        do 6 i = 1,nmax
	   nnlist(i) = 0
	   ntrilist(i) = 0
 6      continue

        do 7 i = 1,np
	   nnn(i) = 0
 7      continue

	do 10 it = 1,nt
	   i1 = vertices(1,it)
	   i2 = vertices(2,it)
	   i3 = vertices(3,it)
	   nnn(i1) = nnn(i1) + 1
	   nnn(i2) = nnn(i2) + 1
	   nnn(i3) = nnn(i3) + 1
 10     continue

c					turn nnn into a running sum
	itemp = nnn(1)+1
	nnn(1) = 1
	do 20 j = 2,np+1
	   itemp2  = itemp 
	   itemp   = itemp + nnn(j)+1
	   nnn(j) = itemp2 + 1
 20     continue
c       write(*,*)' size of array =',nnn(np+1)-1
c       write(*,*)' 3nt+np        =',3*nt+np

	if(nnn(np+1).ge.nmax)then
           write(*,*)'Error: array sizes too small in subroutine '
     &               ,'build_nv'
           write(*,*)'       maximum number of neighbours for all nodes'
           write(*,*)'       is too small: current value =',nmax
           write(*,*)'       Increase size of parameter nmax'
           write(*,*)'       to at least',nnn(np+1)
           write(*,*)'       This will be satisfied if nmax is set'
           write(*,*)'       to 3*nt_max+np_max in calling program' 
	   stop
	end if

	do 25 it = 1,nt
	   i1 = vertices(1,it) 
	   i2 = vertices(2,it) 
	   i3 = vertices(3,it) 
c						compare neighbours i1 i2
c						(remove go to ?)
	   j1 = nnn(i1)
	   j2 = nnn(i1+1) - 1 
	   jt = 0
	   do 30 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i2
	         ntrilist(j) = it
c						if we have recorded connection
c						then jump out of loop
	         go to 31
	      else if(nnlist(j).eq.i2.and.ntrilist(j).ne.it)then
                 jt = ntrilist(j)
	         go to 31
	      end if
  30       continue
  31       continue
c						if neighbours are found then
c						skip second loop 
	   if(jt.eq.0)then
	      j1 = nnn(i2)
	      j2 = nnn(i2+1) - 1 
	      do 32 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i1
	            ntrilist(j) = it
c						if we have inserted connection
c						then jump out of loop
	            go to 33
	         end if
  32          continue
	   end if
  33       continue

	   if(jt.ne.0)then
c						found neighbours it,jt with
c						common nodes i1 and i2
	      neighbour(3,it) = jt
	      k1 = vertices(1,jt)
	      k2 = vertices(2,jt)
	      k3 = vertices(3,jt)
	      if(k1.ne.i1.and.k1.ne.i2)then 
	         neighbour(1,jt) = it
	      else if(k2.ne.i1.and.k2.ne.i2)then 
	         neighbour(2,jt) = it
	      else
	         neighbour(3,jt) = it
	      end if
           end if
c						compare neighbours i1 i3
	   jt = 0
	   j1 = nnn(i1)
	   j2 = nnn(i1+1) - 1 
	   do 130 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i3
	         ntrilist(j) = it
c						if we have recorded connection
c						then jump out of loop
	         go to 131
	      else if(nnlist(j).eq.i3.and.ntrilist(j).ne.it)then
                 jt = ntrilist(j)
	         go to 131
	      end if
  130      continue
  131      continue
c						if neighbours are found then
c						skip second loop 
	   if(jt.eq.0)then
	      j1 = nnn(i3)
	      j2 = nnn(i3+1) - 1 
	      do 132 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i1
	            ntrilist(j) = it
c						if we have inserted connection
c						then jump out of loop
	            go to 133
	         end if
  132         continue
	      end if
  133      continue
	   if(jt.ne.0)then
c						found neighbours it,jt with
c						common nodes i1 and i3
	     neighbour(2,it) = jt
	     k1 = vertices(1,jt) 
	     k2 = vertices(2,jt)
	     k3 = vertices(3,jt)
	     if(k1.ne.i1.and.k1.ne.i3)then 
	        neighbour(1,jt) = it
	     else if(k2.ne.i1.and.k2.ne.i3)then 
	        neighbour(2,jt) = it
	     else
	        neighbour(3,jt) = it
	     end if
           end if
c						compare neighbours i2 i3
	   jt = 0
	   j1 = nnn(i2)
	   j2 = nnn(i2+1) - 1 
	   do 230 j = j1,j2
	      if(nnlist(j).eq.0)then
	         nnlist(j) = i3
	         ntrilist(j) = it
c						if we have recorded connection
c						then jump out of loop
	         go to 231
 	      else if(nnlist(j).eq.i3.and.ntrilist(j).ne.it)then
                 jt = ntrilist(j)
 	         go to 231
	      end if
  230      continue
  231      continue
c						if neighbours are found then
c						skip second loop 
	   if(jt.eq.0)then
	      j1 = nnn(i3)
	      j2 = nnn(i3+1) - 1 
	      do 232 j = j1,j2
	         if(nnlist(j).eq.0)then
	            nnlist(j) = i2
	            ntrilist(j) = it
c						if we have inserted connection
c						then jump out of loop
	            go to 233
	         end if
  232         continue
	      end if
  233      continue
	   if(jt.ne.0)then
c						found neighbours it,jt with
c						common nodes i2 and i3
	     neighbour(1,it) = jt
	     k1 = vertices(1,jt) 
	     k2 = vertices(2,jt)
	     k3 = vertices(3,jt)
	     if(k1.ne.i2.and.k1.ne.i3)then 
	        neighbour(1,jt) = it
	     else if(k2.ne.i2.and.k2.ne.i3)then 
	        neighbour(2,jt) = it
	     else
	        neighbour(3,jt) = it
	     end if
           end if

 25     continue

	if(nnwrite)write(*,*)' built neighbour v'

	return
	end
c
c------------------------------------------------------------------------
c
c       Calculate_hulltriangles - finds all triangles with a face on
c                                 the convex hull by searching through
c                                 the entries in the array neighbour.
c
c       Input:
c               neighbour(3,nt)         array of neighbouring tetrahedra
c               nt                      number of tetrahedra
c	        nh_max			maximum number of triangles on convex 
c					hull (max size of array hulltriangles)
c		nohalt			determines error response
c       Output:
c		hulltriangles(nh)	array of triangles with an edge
c					on the convex hull
c               nh                      number of tetrahedra with an edge
c                                       on the convex hull
c       Comments:
c
c                This routine fills up the array hulltriangles which
c                is only used by routine nn2Do, i.e the `pseudo-extension' 
c		 Watson's nn-interpolation method to points outside of the 
c		 convex hull. If nnext is set to false then hulltriangles
c		 is never used and the array can be set to size 1.
c
c		 If nohalt = 0 then the routine will stop with an error
c		 message if nh > nh_max. If nohalt .ne. 0 and nh > nh_max
c		 then it will return nh = -1. 
c
c                No calls to other routines.
c
c                                       M. Sambridge, RSES, May 1995.
c
c------------------------------------------------------------------------
c
	Subroutine calculate_hulltriangles
     &             (neighbour,nt,nh_max,hulltriangles,nh,nohalt)
c
	integer		neighbour(3,*)
	integer		hulltriangles(*)

c                                               store list of triangles
c                                               which have an edge on the
c                                               convex hull.
c                                               (used by routine nn2D)
        nh = 1
        do 100 j = 1,nt
           if(neighbour(1,j).eq.0.or.
     &        neighbour(2,j).eq.0.or.
     &        neighbour(3,j).eq.0)then
              hulltriangles(nh) = j
              nh = nh + 1
              if(nh.gt.nh_max.and.nohalt.eq.0)then
                  write(*,*)' Error array storing outward facing '
                  write(*,*)' triangles on convex hull is too small.'
                  write(*,*)' Increase size of parameter nh_max'
                  stop
              else if(nh.gt.nh_max.and.nohalt.ne.0)then
                  nh = -1
                  return
              end if
           end if
 100    continue
        nh = nh -1
 
	return
	end
c
c------------------------------------------------------------------------
c
c	Triloc - locates the triangle containing point x,y
c
c	Input:
c		x,y			co-ordinates of input points	
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c	        loc			first guess of triangle containing
c					(x, y).
c
c	Output:
c	        loc			index of triangle containing 
c					input point.
c	        out			=true if (x,y) is outside of
c					the convex hull, otherwise = false. 
c
c	Comments:
c		 If (x,y) is outside convex hull loc is a `nearby' triangle
c		 on the hull and out is set to true.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine Triloc(x,y,points,vertices,neighbour,loc,out)
c
	real*8		points(2,*)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		p1,p2
        logical		out
	real*8		x,y,del1,del2
        integer		c1(3)
        data		c1/2,3,1/
c
	out = .false.

 10     continue
c					point is outside convex hull
        if(out)return

        do 20 i=1,3
	   j = c1(i)
	   k = c1(j)
           p1 = vertices(i,loc)
           p2 = vertices(j,loc)
	   del1 = (points(2,p1)-y)*(points(1,p2)-x)
	   del2 = (points(1,p1)-x)*(points(2,p2)-y)
	   if(del1.gt.del2)then
	      if(neighbour(k,loc).eq.0)then
                 out = .true.
	      else
	         loc = neighbour(k,loc)
	      end if
	      go to 10
	   end if
 20     continue
	
	return
	end
c
c------------------------------------------------------------------------
c
c	Ccentres - calculates centres of all Delaunay circumcircles
c
c
c	Input:
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        nt			number of triangles
c
c	Output:
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c
c	Comments:
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine ccentres(points,vertices,nt,centres)
c
	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		x1,x2,x3,y1,y2,y3,x,y
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
	integer		vertices(3,*)
c						Find centres of all
c						Delaunay Circumcircles
	do 5 i= 1,nt

	   x1 = points(1,vertices(1,i))
	   x2 = points(1,vertices(2,i))
	   x3 = points(1,vertices(3,i))
	   y1 = points(2,vertices(1,i))
	   y2 = points(2,vertices(2,i))
	   y3 = points(2,vertices(3,i))

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1
	   x = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1*0.5d0
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0*dy2m1)/
     &         (denom)

	   y = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1)*0.5d0 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1)*0.5d0)/ 
     &         (denom)

	   centres(1,i) = x
	   centres(2,i) = y
           x1 = x - x1
           y1 = y - y1
	   centres(3,i) = x1*x1 + y1*y1

 5	continue

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2Di - calculates natural neighbour interpolation at point x,y
c		when point is inside convex hull 
c
c	Input:
c		x,y			co-ordinates of input point
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c	        loc			index of triangle containing 
c					input point.
c		data			data values at each node
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		v			area of voronoi cell about (x,y)
c
c	Comments:
c		 On input point (x,y) must be in triangle loc. 
c		 Note: The input triangle contains (x,y) and so it's 
c		 circumcircle must also contain (x,y).
c
c		 See comments on plot routine plot_tc above.
c
c		 This routine uses Sloan's LIFO stack technique to 
c		 avoid an expensive global search over all triangles.
c		 (see Sambridge, Braun and McQueen, 1995; 
c		  Geophys. J. Int., vol 122, 837-857).
c
c		 Calls are made to: nn_tri. plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty. 
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Di
     &             (x,y,points,vertices,neighbour,
     &              centres,loc,data,f,df,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		x,y,dist,dx,dy
	real*8		f,df(2),v,dvx,dvy
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tri,pos
	logical		nnwrite

	common/nnswitches/lud,nnwrite

 	if(nnwrite)write(lud,*)' Result of nn2Di search:'

c					Initialize variables
	f = 0.d0
	df(1) = 0.d0
	df(2) = 0.d0
	v = 0.d0
	dvx = 0.d0
	dvy = 0.d0
c
c					Calculate natural neighbour
c					contribution from input triangle.
c
 	call nn_tri(x,y,points,vertices,centres,loc,
     &              data,f,df,v,dvx,dvy)

	if(nnwrite)write(lud,*)loc
c
c					Find all other circumcircles 
c					containing input point and update
c					natural neighbour contribution
c					from these triangles
c
c					plot circum-triangle and circum-circle

 	call plot_tc(loc,points,vertices,centres)

	   
c					initialize stack
	call stackpairinit

c					put input triangle's neighbouring 
c					triangles on LIFO stack
c					together with position of input
c					triangle in their neighbour list

	i = neighbour(1,loc)
	j = neighbour(2,loc)
	k = neighbour(3,loc)
	if(i.ne.0)then
	   if(neighbour(1,i).eq.loc)call pushpair(i,1)
	   if(neighbour(2,i).eq.loc)call pushpair(i,2)
	   if(neighbour(3,i).eq.loc)call pushpair(i,3)
	end if
	if(j.ne.0)then
	   if(neighbour(1,j).eq.loc)call pushpair(j,1)
	   if(neighbour(2,j).eq.loc)call pushpair(j,2)
	   if(neighbour(3,j).eq.loc)call pushpair(j,3)
	end if
	if(k.ne.0)then
	   if(neighbour(1,k).eq.loc)call pushpair(k,1)
	   if(neighbour(2,k).eq.loc)call pushpair(k,2)
	   if(neighbour(3,k).eq.loc)call pushpair(k,3)
	end if

 10	call stackpairempty(k)
c					if stack empty then finish
	if(k.eq.1)go to 100
c					take triangle from stack
	call poppair(tri,pos)

c					test if (x,y) is in circumcircle

         
        dx = x-centres(1,tri)
        dy = y-centres(2,tri)
	dist = dx*dx + dy*dy 

 	if(dist.lt.centres(3,tri))then
           if(nnwrite)write(lud,*)tri
c					Calculate natural neighbour
c					contribution from current triangle.
c
 	   call nn_tri(x,y,points,vertices,centres,tri
     &                 ,data,f,df,v,dvx,dvy)

c					plot triangle and circumcircle

 	   call plot_tc(tri,points,vertices,centres)
	   if(pos.eq.1)then
	      new1 = neighbour(2,tri)
	      new2 = neighbour(3,tri)
	   end if
	   if(pos.eq.2)then
	      new1 = neighbour(1,tri)
	      new2 = neighbour(3,tri)
	   end if
	   if(pos.eq.3)then
	      new1 = neighbour(1,tri)
	      new2 = neighbour(2,tri)
	   end if
	   if(new1.ne.0)then
	      if(neighbour(1,new1).eq.tri)call pushpair(new1,1)
	      if(neighbour(2,new1).eq.tri)call pushpair(new1,2)
	      if(neighbour(3,new1).eq.tri)call pushpair(new1,3)
	   end if
	   if(new2.ne.0)then
	      if(neighbour(1,new2).eq.tri)call pushpair(new2,1)
	      if(neighbour(2,new2).eq.tri)call pushpair(new2,2)
	      if(neighbour(3,new2).eq.tri)call pushpair(new2,3)
	   end if
	end if
	go to 10
 
 100    continue
        call stackpairflush()
	if(v.ne.0.0)then
	   f = f/v
	   df(1) = (df(1) - f*dvx)/v
	   df(2) = (df(2) - f*dvy)/v
c	   if(nnwrite)then
c	   write(*,*)' Overlap voronoi area = ',v
c	   write(*,*)' Normalized value of f   =',f
c	   write(*,*)' Normalized value of dfx =',df(1)
c	   write(*,*)' Normalized value of dfy =',df(2)
c	   end if
c	else
c	   if(nnwrite)then
c	   write(*,*)' Subroutine nn2Di: overlap voronoi area = 0'
c	   write(*,*)' Unnormalized value of f   =',f
c	   write(*,*)' Current value of dfx      =',df(1)
c	   write(*,*)' Current value of dfy      =',df(2)
c	   end if
	end  if

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2Do - calculates natural neighbour interpolation at point x,y
c		when point is outside convex hull 
c
c	Input:
c		x,y			co-ordinates of input point
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c		hulltriangles(nh)	array of triangles with an edge
c					on the convex hull
c	        nh			number of triangles with an edge
c					on the convex hull
c		data			data values at each node
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		v			area of voronoi cell about (x,y)
c
c	Comments:
c		The following procedure works on the principle that if (x,y) is
c		inside the union of circumcircles then it must be in one 
c		of the circumcircles whose Delaunay triangle has an edge 
c		on the convex hull. It then proceeds to find all circumcircles
c		containing the point by repeating the triangle walking
c		algorithm described by Sloan (1987, and attributed to 
c		D. Watson) starting from each triangle on the hull whose
c		circumcircles contain the point. Note the repeated use of
c		the procedure is necessary because when (x,y) is outside of 
c		the hull the required triangles may not form a contiguous set.
c
c		See comments on plot routine plot_tc above.
c
c		Calls are made to: nn_tri. plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty. 
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Do
     &  (x,y,points,vertices,neighbour,centres,
     &   hulltriangles,nh,data,f,df,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		x,y,xc,yc,dist
	real*8		data(*)
	real*8		f,df(2),v,dvx,dvy
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tri,pos
        integer		hulltriangles(*)
	logical		inside_union
	logical		nnwrite

	common/nnswitches/lud,nnwrite

 	if(nnwrite)write(lud,*)' Result of nn2Do search:'
        inside_union = .false.

c					Initialize variables
	f = 0.d0
	df(1) = 0.d0
	df(2) = 0.d0
	v = 0.d0
	dvx = 0.d0
	dvy = 0.d0
c					Find all circumcircles containing
c					input point

c					test all circumcircles whose
c					Delaunay triangles have an edge
c					on the convex hull
	do 5 list = 1,nh
           tri = hulltriangles(list)
	   xc = centres(1,tri)
	   yc = centres(2,tri)
	   dist = (xc-x)*(xc-x) + (yc-y)*(yc-y)
	   if(dist.lt.centres(3,tri))then
	      if(nnwrite)write(lud,*)' circle test found ',tri
	      call plot_tc(tri,points,vertices,centres)
	      inside_union = .True.
	      loc = tri


	      if(nnwrite)write(lud,*)loc
c					Calculate natural neighbour
c					contribution from current triangle.

 	      call nn_tri(x,y,points,vertices,centres,loc,
     &                    data,f,df,v,dvx,dvy)


c					plot triangle and circumcircle

	      call plot_tc(loc,points,vertices,centres)
	   
c					initialize stack
	      call stackpairinit
c					put input triangle's neighbouring 
c					triangles on LIFO stack
c					together with position of input
c					triangle in their neighbour list
	      i = neighbour(1,loc)
	      j = neighbour(2,loc)
	      k = neighbour(3,loc)
	      if(i.ne.0)then
	         if(neighbour(1,i).eq.loc)call pushpair(i,1)
	         if(neighbour(2,i).eq.loc)call pushpair(i,2)
	         if(neighbour(3,i).eq.loc)call pushpair(i,3)
	      end if
	      if(j.ne.0)then
	         if(neighbour(1,j).eq.loc)call pushpair(j,1)
	         if(neighbour(2,j).eq.loc)call pushpair(j,2)
	         if(neighbour(3,j).eq.loc)call pushpair(j,3)
	      end if
	      if(k.ne.0)then
	         if(neighbour(1,k).eq.loc)call pushpair(k,1)
	         if(neighbour(2,k).eq.loc)call pushpair(k,2)
	         if(neighbour(3,k).eq.loc)call pushpair(k,3)
	      end if
         
 10	      call stackpairempty(k)
c					if stack empty then finish
	      if(k.eq.1)go to 100
c					take triangle from stack
	      call poppair(tri,pos)

c					test if (x,y) is in circumcircle

	      dist = (x-centres(1,tri))*(x-centres(1,tri)) + 
     &               (y-centres(2,tri))*(y-centres(2,tri))
	      if(dist.lt.centres(3,tri))then
                 if(nnwrite)write(lud,*)tri
c					Calculate natural neighbour
c					contribution from current triangle.

 	         call nn_tri(x,y,points,vertices,centres,tri,
     &                       data,f,df,v,dvx,dvy)

c					plot triangle and circumcircle
         
	         call plot_tc(tri,points,vertices,centres)
	         if(pos.eq.1)then
	            new1 = neighbour(2,tri)
	            new2 = neighbour(3,tri)
	         end if
	         if(pos.eq.2)then
	            new1 = neighbour(1,tri)
	            new2 = neighbour(3,tri)
	         end if
	         if(pos.eq.3)then
	            new1 = neighbour(1,tri)
	            new2 = neighbour(2,tri)
	         end if
	         if(new1.ne.0)then
	            if(neighbour(1,new1).eq.tri)call pushpair(new1,1)
	            if(neighbour(2,new1).eq.tri)call pushpair(new1,2)
	            if(neighbour(3,new1).eq.tri)call pushpair(new1,3)
	         end if
	         if(new2.ne.0)then
	            if(neighbour(1,new2).eq.tri)call pushpair(new2,1)
	            if(neighbour(2,new2).eq.tri)call pushpair(new2,2)
	            if(neighbour(3,new2).eq.tri)call pushpair(new2,3)
	         end if
	      end if
	      go to 10
 100          continue
              call stackpairflush()
	   end if
 5      continue
	if(.not.inside_union.and.nnwrite)
     &     write(lud,*)' point is outside all circumcircles'
	if(.not.inside_union.and.nnwrite)
     &     write(*,*)' point is outside all circumcircles'

	if(v.ne.0.0)then
	   f = f/v
	   df(1) = (df(1) - f*dvx)/v
	   df(2) = (df(2) - f*dvy)/v
	end  if

	return
	end
c
c------------------------------------------------------------------------
c
c	nn_tri - calculates natural neighbour contribution from 
c		 triangle tri at point x,y
c
c	Input:
c		x,y			co-ordinates of input point
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c		tri			index of current triangle
c	Output:
c		f			partial f interpolation
c		df(2)			partial derivative interpolation
c		dvx			partial sum of derivatives
c		dvy			partial sum of derivatives
c		v			partial area of voronoi cell about (x,y)
c
c	Comments:
c		 Note this routine only adds on the contribution to
c		 f and derivatives to sum from last call. All summation
c		 variables are initialized in the calling routine. 
c		 Nodes are assumed to be in counterclockwise direction.
c
c		 The terms dvx, dvy are needed in the calling routine to 
c		 complete the calculation of f and df. 
c
c		 See comments on plot routine plot_c above.
c
c		 Calls are made to: plot_c.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
	Subroutine nn_tri
     &             (x1,x2,points,vertices,centres,tri,
     &              data,f,df,v,dvx,dvy)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		x1,x2,f,df(2),v
	real*8		dvx,dvy,a,b,dxw,dyw,det
	real*8		p(2,3),c(2,3),f1(2,3),f2(2,3),cx,cy
	real*8		dcx(2,3),dcy(2,3)
	real*8		data(*)
	integer		vertices(3,*)
	integer		tri
	logical		nnwrite
        integer		c1(3)
        data		c1/2,3,1/

	common/nnswitches/lud,nnwrite



c					set position vectors of nodes
	do 5 i=1,2
	   do 6 j=1,3
	      p(i,j) = points(i,vertices(j,tri))
 6         continue
 5      continue
c					set circumcentre nodes
	cx = centres(1,tri)
	cy = centres(2,tri)
	do 10 i = 1,3
           j = c1(i)
           k = c1(j)

c					calculate centres of three 
c					new triangles

           c(1,i) = (((p(1,k)-p(1,j))*(p(1,k)+p(1,j))/2. 
     &               +(p(2,k)-p(2,j))*(p(2,k)+p(2,j))/2.)
     &               *(x2-p(2,j))
     &             -(((x1-p(1,j))*(x1+p(1,j))/2. 
     &               +(x2-p(2,j))*(x2+p(2,j))/2.)
     &               *(p(2,k)-p(2,j))))/
     &               ((p(1,k)-p(1,j))*(x2-p(2,j))
     &               -(x1-p(1,j))*(p(2,k)-p(2,j)))
 
           c(2,i) = ((p(1,k)-p(1,j))*((x1-p(1,j))
     &              *(x1+p(1,j))/2. 
     &              +(x2-p(2,j))*(x2+p(2,j))/2.)
     &             -((x1-p(1,j))*((p(1,k)-p(1,j))
     &              *(p(1,k)+p(1,j))/2. 
     &              +(p(2,k)-p(2,j))*(p(2,k)+p(2,j))/2.)))/
     &              ((p(1,k)-p(1,j))*(x2-p(2,j))
     &              -(x1-p(1,j))*(p(2,k)-p(2,j)))
 
           xs = c(1,i)
           ys = c(2,i)
           call plot_c(xs,ys,x1,x2)

	   f1(1,i) = c(1,i) - x1
	   f1(2,i) = f1(1,i) 
	   f2(1,i) = c(2,i) - x2
	   f2(2,i) = f2(1,i) 

	   dcx(1,i) = ((p(2,k)-x2)*f1(1,i)-(p(2,j)-x2)*f1(2,i))/
     &               ((p(1,j)-x1)*(p(2,k)-x2)-(p(2,j)-x2)*(p(1,k)-x1))
	   dcx(2,i) = ((p(1,j)-x1)*f1(2,i)-(p(1,k)-x1)*f1(1,i))/
     &               ((p(1,j)-x1)*(p(2,k)-x2)-(p(2,j)-x2)*(p(1,k)-x1))
	   dcy(1,i) = ((p(2,k)-x2)*f2(1,i)-(p(2,j)-x2)*f2(2,i))/
     &               ((p(1,j)-x1)*(p(2,k)-x2)-(p(2,j)-x2)*(p(1,k)-x1))
	   dcy(2,i) = ((p(1,j)-x1)*f2(2,i)-(p(1,k)-x1)*f2(1,i))/
     &               ((p(1,j)-x1)*(p(2,k)-x2)-(p(2,j)-x2)*(p(1,k)-x1))
 10	continue

	do 20 i = 1,3
           j = c1(i)
           k = c1(j)
	   node = vertices(i,tri)
	   det  = ((c(1,j)-cx)*(c(2,k)-cy)
     &            -(c(1,k)-cx)*(c(2,j)-cy))/2.d0
	   v = v + det
	   f = f + det*data(node)
	   a = (dcx(1,j)*(c(2,k)-cy)-(c(1,k)-cx)*dcx(2,j))
	   b = ((c(1,j)-cx)*dcx(2,k)-dcx(1,k)*(c(2,j)-cy))
	   dxw = (a + b)/2.d0
	   a = (dcy(1,j)*(c(2,k)-cy)-(c(1,k)-cx)*dcy(2,j))
	   b = ((c(1,j)-cx)*dcy(2,k)-dcy(1,k)*(c(2,j)-cy))
	   dyw = (a + b)/2.d0
	   df(1) = df(1) + dxw*data(node)
	   df(2) = df(2) + dyw*data(node)
	   dvx =  dvx  + dxw
	   dvy =  dvy  + dyw
 20	continue

c	if(nnwrite)write(*,*)'x',x1,' y',x2,' f ',f,' t ',tri
c	if(nnwrite)write(*,*)' df1',df(1),' df2 ',df(2),' v',v

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2Dr - calculates natural neighbour interpolation at point x,y
c		when point is inside convex hull 
c		(using method based on recursive formula)
c
c	Input:
c		x,y			co-ordinates of input point
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c	        loc			index of triangle containing 
c					input point.
c		data			data values at each node
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                        set to ~50 in calling program)
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		v			area of voronoi cell about (x,y)
c
c	Comments:
c		 On input point (x,y) must be in triangle loc. 
c		 Note: The input triangle contains (x,y) and so it's 
c		 circumcircle must also contain (x,y).
c
c		 In order to use the recursive Lasserre formula to
c		 calculate the nn co-ordinate of the input point with
c		 respect to its ith neighbour we must find the set of
c		 nodes that are both neighbours of i and (x,y), where 
c		 the neighbours of i are determined before addition of 
c		 (x,y). This is done using an approach similar to Watson's
c		 algorithm for updating a Delaunay triangulation, although
c		 we avoid the complication of having to remove `double
c		 entires in the list' and also use Sloan's LIFO stack
c		 approach to avoid an expensive global search over
c		 all triangles.
c
c		 This version does not calculate derivatives.
c
c		 Calls are made to: plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty. 
c
c					M. Sambridge, RSES, April 1995.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Dr
     &             (x,y,points,vertices,neighbour,centres,loc,
     &              data,nnpn_max,a,b,p,nodes,pairs,tstore,f,df,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		x,y,dist,dx,dy
	real*8		f,df(2),v
	real*8		xo(2)
        real*8  	xs(2)
        real*8          p(2,nnpn_max)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tri,pos,pcount
	logical		nnwrite
	integer		pairs(2,nnpn_max)
	integer		nodes(nnpn_max)
	integer		tstore(nnpn_max)
        real            a(nnpn_max,2)
        real            b(nnpn_max)

        

	common/nnswitches/lud,nnwrite

 	if(nnwrite)write(lud,*)' Result of nn2Dr search:'

c					Initialize variables
	f = 0.d0
	df(1) = 0.d0
	df(2) = 0.d0
	v = 0.d0

c					Determine list of 
c					neighbours to input point.
	
c
c					Initialize pair list 
c

        pairs(1,1) = vertices(1,loc)
        pairs(2,1) = vertices(2,loc)
        pairs(1,2) = vertices(2,loc)
        pairs(2,2) = vertices(3,loc)
        pairs(1,3) = vertices(3,loc)
        pairs(2,3) = vertices(1,loc)
        pcount = 3

	if(nnwrite)write(lud,*)loc

c					plot input circum-triangle 
c					and circum-circle

 	call plot_tc(loc,points,vertices,centres)

c					record circum-triangle 
        tstore(1) = loc
        tstore(2) = loc
        tstore(3) = loc
c
c					Find all other circumcircles 
c					containing input point and update
c					pair list from these triangles
c
c					initialize stack
	call stackpairinit


c					put input triangle's neighbouring 
c					triangles on LIFO stack
c					together with position of input
c					triangle in their neighbour list

	i = neighbour(1,loc)
	j = neighbour(2,loc)
	k = neighbour(3,loc)
	if(i.ne.0)then
	   if(neighbour(1,i).eq.loc)call pushpair(i,1)
	   if(neighbour(2,i).eq.loc)call pushpair(i,2)
	   if(neighbour(3,i).eq.loc)call pushpair(i,3)
	end if
	if(j.ne.0)then
	   if(neighbour(1,j).eq.loc)call pushpair(j,1)
	   if(neighbour(2,j).eq.loc)call pushpair(j,2)
	   if(neighbour(3,j).eq.loc)call pushpair(j,3)
	end if
	if(k.ne.0)then
	   if(neighbour(1,k).eq.loc)call pushpair(k,1)
	   if(neighbour(2,k).eq.loc)call pushpair(k,2)
	   if(neighbour(3,k).eq.loc)call pushpair(k,3)
	end if

 10	call stackpairempty(k)
c					if stack empty then finish
	if(k.eq.1)go to 100
c					take triangle from stack
	call poppair(tri,pos)

c					test if (x,y) is in circumcircle

        dx = x-centres(1,tri)
        dy = y-centres(2,tri)
	dist = dx*dx + dy*dy 

 	if(dist.lt.centres(3,tri))then
           if(nnwrite)write(lud,*)tri

c					plot circum-triangle and circum-circle

 	   call plot_tc(tri,points,vertices,centres)

c					record each edge of current
c					triangles and triangle
c
	   if(pos.eq.1)then
	      new1 = neighbour(2,tri)
	      new2 = neighbour(3,tri)
              n1 = 3
              n2 = 1
              n3 = 1
              n4 = 2
	   else if(pos.eq.2)then
	      new1 = neighbour(1,tri)
	      new2 = neighbour(3,tri)
              n1 = 2
              n2 = 3
              n3 = 1
              n4 = 2
	   else 
	      new1 = neighbour(1,tri)
	      new2 = neighbour(2,tri)
              n1 = 2
              n2 = 3
              n3 = 3
              n4 = 1
	   end if
           pcount = pcount + 1
           pairs(1,pcount) = vertices(n1,tri)
           pairs(2,pcount) = vertices(n2,tri)
           tstore(pcount) = tri
           pcount = pcount + 1
           pairs(1,pcount) = vertices(n3,tri)
           pairs(2,pcount) = vertices(n4,tri)
           tstore(pcount) = tri

	   if(new1.ne.0)then
	      if(neighbour(1,new1).eq.tri)call pushpair(new1,1)
	      if(neighbour(2,new1).eq.tri)call pushpair(new1,2)
	      if(neighbour(3,new1).eq.tri)call pushpair(new1,3)
	   end if
	   if(new2.ne.0)then
	      if(neighbour(1,new2).eq.tri)call pushpair(new2,1)
	      if(neighbour(2,new2).eq.tri)call pushpair(new2,2)
	      if(neighbour(3,new2).eq.tri)call pushpair(new2,3)
	   end if
	end if
	go to 10
 
 100    continue
        call stackpairflush()

c						check size of pcount
        if(pcount.gt.nnpn_max)then
           write(*,*)' '
           write(*,*)' Error: work arrays in subroutine nn2Dr',
     &               ' not big enough.'
           write(*,*)' Remedy: Increase size of parameter'
           write(*,*)'         nnpn_max in calling program'
           write(*,*)'         current value = ',nnpn_max
           write(*,*)'         required value >= ',pcount
	   stop
        end if
c						write out neighbouring
c						nodes of input points
c
c       do 110 i=1,pcount
c          write(*,*)' neighbouring pair of nodes',pairs(1,i),pairs(2,i)
c110    continue
c
c						find set of common neighbours
c						between input point, and each
c						of its natural neighbours
c						and use this information to
c						set up p matrix of co-ordinates.
c
        j=0
        xs(1) = x
        xs(2) = y
        do 120 i=1,pcount
	   tri = tstore(i)
           do 121 ii=1,2
	      node = pairs(ii,i)
              do 122 kk=1,j
                 if(node.eq.nodes(kk))go to 121
 122          continue
              j = j + 1
	      nodes(j) = node
              p(1,1) = points(1,node)
              p(2,1) = points(2,node)
              jj = 1
              do 130 k=1,pcount
                 if(pairs(1,k).eq.node)then
c                   write(*,*)' Node:',node,' neighbour',pairs(2,k)
                    jj = jj + 1
                    p(1,jj) = points(1,pairs(2,k))
                    p(2,jj) = points(2,pairs(2,k))
                 end if
                 if(pairs(2,k).eq.node)then
c                   write(*,*)' Node:',node,' neighbour',pairs(1,k)
                    jj = jj + 1
                    p(1,jj) = points(1,pairs(1,k))
                    p(2,jj) = points(2,pairs(1,k))
                 end if
 130          continue
              
c						calculate volume of 
c						second-order voronoi cell 
c						for current node using
c						recursive formula and reset
c						origin to circumcentre
c
              xo(1) = centres(1,tri)
              xo(2) = centres(2,tri)

              call second_voronoi (xs,p,2,jj,nnpn_max,2,xo,a,b,vol)
c
              v = v + vol
              f = f + vol*data(node)
 121       continue
 120    continue

        do 140 k=1,j
           nodes(k) = 0
 140    continue

 	if(v.ne.0.0) f = f/v

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2Drd - calculates natural neighbour interpolation at point x,y
c	 	 when point is inside convex hull 
c		 (using method based on recursive formulae)
c
c	Input:
c		x,y			co-ordinates of input point
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c	        loc			index of triangle containing 
c					input point.
c		data			data values at each node
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                        set to ~50 in calling program)
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		v			area of voronoi cell about (x,y)
c
c	Comments:
c		 On input point (x,y) must be in triangle loc. 
c		 Note: The input triangle contains (x,y) and so it's 
c		 circumcircle must also contain (x,y).
c
c		 In order to use the recursive Lasserre formula to
c		 calculate the nn co-ordinate of the input point with
c		 respect to its ith neighbour we must find the set of
c		 nodes that are both neighbours of i and (x,y), where 
c		 the neighbours of i are determined before addition of 
c		 (x,y). This is done using an approach similar to Watson's
c		 algorithm for updating a Delaunay triangulation, although
c		 we avoid the complication of having to remove `double
c		 entires in the list' and also use Sloan's LIFO stack
c		 approach to avoid an expensive global search over
c		 all triangles.
c
c		 This version also calculates derivatives.
c
c		 Calls are made to: plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty. 
c
c					M. Sambridge, RSES, April 1995.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Drd
     &             (x,y,points,vertices,neighbour,centres,loc,
     &              data,nnpn_max,a,b,p,nodes,pairs,tstore,f,df,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		x,y,dist,dx,dy
	real*8		f,df(2),v,dvx,dvy
	real*8		xo(2)
        real*8  	xs(2)
        real*8          p(2,nnpn_max)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tri,pos,pcount
	logical		nnwrite
	integer		pairs(2,nnpn_max)
	integer		nodes(nnpn_max)
	integer		tstore(nnpn_max)
        real  		deriv(2)
        real            a(nnpn_max,2)
        real            b(nnpn_max)

        

	common/nnswitches/lud,nnwrite

 	if(nnwrite)write(lud,*)' Result of nn2Dr search:'

c					Initialize variables
	f = 0.d0
	df(1) = 0.d0
	df(2) = 0.d0
	v = 0.d0
	dvx = 0.d0
	dvy = 0.d0
c					Initialize pair list 
c
        pairs(1,1) = vertices(1,loc)
        pairs(2,1) = vertices(2,loc)
        pairs(1,2) = vertices(2,loc)
        pairs(2,2) = vertices(3,loc)
        pairs(1,3) = vertices(3,loc)
        pairs(2,3) = vertices(1,loc)
        pcount = 3

	if(nnwrite)write(lud,*)loc

c					plot input circum-triangle 
c					and circum-circle

 	call plot_tc(loc,points,vertices,centres)

c					record circum-triangle 
        tstore(1) = loc
        tstore(2) = loc
        tstore(3) = loc
c
c					Find all other circumcircles 
c					containing input point and update
c					pair list from these triangles
c
c					initialize stack
	call stackpairinit


c					put input triangle's neighbouring 
c					triangles on LIFO stack
c					together with position of input
c					triangle in their neighbour list

	i = neighbour(1,loc)
	j = neighbour(2,loc)
	k = neighbour(3,loc)
	if(i.ne.0)then
	   if(neighbour(1,i).eq.loc)call pushpair(i,1)
	   if(neighbour(2,i).eq.loc)call pushpair(i,2)
	   if(neighbour(3,i).eq.loc)call pushpair(i,3)
	end if
	if(j.ne.0)then
	   if(neighbour(1,j).eq.loc)call pushpair(j,1)
	   if(neighbour(2,j).eq.loc)call pushpair(j,2)
	   if(neighbour(3,j).eq.loc)call pushpair(j,3)
	end if
	if(k.ne.0)then
	   if(neighbour(1,k).eq.loc)call pushpair(k,1)
	   if(neighbour(2,k).eq.loc)call pushpair(k,2)
	   if(neighbour(3,k).eq.loc)call pushpair(k,3)
	end if

 10	call stackpairempty(k)
c					if stack empty then finish
	if(k.eq.1)go to 100
c					take triangle from stack
	call poppair(tri,pos)

c					test if (x,y) is in circumcircle

        dx = x-centres(1,tri)
        dy = y-centres(2,tri)
	dist = dx*dx + dy*dy 

 	if(dist.lt.centres(3,tri))then
           if(nnwrite)write(lud,*)tri

c					plot input circum-triangle 
c					and circum-circle

 	   call plot_tc(tri,points,vertices,centres)

c					record each edge 
c					of current triangle.
c
	   if(pos.eq.1)then
	      new1 = neighbour(2,tri)
	      new2 = neighbour(3,tri)
              n1 = 3
              n2 = 1
              n3 = 1
              n4 = 2
	   else if(pos.eq.2)then
	      new1 = neighbour(1,tri)
	      new2 = neighbour(3,tri)
              n1 = 2
              n2 = 3
              n3 = 1
              n4 = 2
	   else 
	      new1 = neighbour(1,tri)
	      new2 = neighbour(2,tri)
              n1 = 2
              n2 = 3
              n3 = 3
              n4 = 1
	   end if
           pcount = pcount + 1
           pairs(1,pcount) = vertices(n1,tri)
           pairs(2,pcount) = vertices(n2,tri)
           tstore(pcount) = tri
           pcount = pcount + 1
           pairs(1,pcount) = vertices(n3,tri)
           pairs(2,pcount) = vertices(n4,tri)
           tstore(pcount) = tri

	   if(new1.ne.0)then
	      if(neighbour(1,new1).eq.tri)call pushpair(new1,1)
	      if(neighbour(2,new1).eq.tri)call pushpair(new1,2)
	      if(neighbour(3,new1).eq.tri)call pushpair(new1,3)
	   end if
	   if(new2.ne.0)then
	      if(neighbour(1,new2).eq.tri)call pushpair(new2,1)
	      if(neighbour(2,new2).eq.tri)call pushpair(new2,2)
	      if(neighbour(3,new2).eq.tri)call pushpair(new2,3)
	   end if
	end if
	go to 10
 
 100    continue
        call stackpairflush()
c						check size of pcount
        if(pcount.gt.nnpn_max)then
           write(*,*)' '
           write(*,*)' Error: work arrays in subroutine nn2Drd',
     &               ' not big enough.'
           write(*,*)' Remedy: Increase size of parameter'
           write(*,*)'         nnpn_max in calling program'
           write(*,*)'         current value = ',nnpn_max
           write(*,*)'         required value >= ',pcount
	   stop
        end if

c						write out neighbouring
c						nodes of input points
c
c       do 110 i=1,pcount
c          write(*,*)' neighbouring pair of nodes',pairs(1,i),pairs(2,i)
c110    continue
c
c						find set of common neighbours
c						between input point, and each
c						of its natural neighbours.
c						and use this information to
c						set up p matrix of co-ordinates.
c
        j=0
        xs(1) = x
        xs(2) = y
        do 120 i=1,pcount
           tri = tstore(i)
           do 121 ii=1,2
	      node = pairs(ii,i)
              do 122 kk=1,j
                 if(node.eq.nodes(kk))go to 121
 122          continue
              j = j + 1
	      nodes(j) = node
              p(1,1) = points(1,node)
              p(2,1) = points(2,node)
              jj = 1
              do 130 k=1,pcount
                 if(pairs(1,k).eq.node)then
c                   write(*,*)' Node:',node,' neighbour',pairs(2,k)
                    jj = jj + 1
                    p(1,jj) = points(1,pairs(2,k))
                    p(2,jj) = points(2,pairs(2,k))
                 end if
                 if(pairs(2,k).eq.node)then
c                   write(*,*)' Node:',node,' neighbour',pairs(1,k)
                    jj = jj + 1
                    p(1,jj) = points(1,pairs(1,k))
                    p(2,jj) = points(2,pairs(1,k))
                 end if
 130          continue
              
c						calculate volume and derivatives
c						of second-order voronoi cell 
c						for current node using
c						recursive formula and reset
c                                               origin to circumcentre
c
              xo(1) = centres(1,tri)
              xo(2) = centres(2,tri)

              call second_voronoi_d 
     &             (xs,p,2,jj,nnpn_max,2,xo,a,b,vol,deriv)
c
              v = v + vol
              f = f + vol*data(node)
              df(1) = df(1) + data(node)*deriv(1)
              df(2) = df(2) + data(node)*deriv(2)
              dvx = dvx + deriv(1)
              dvy = dvy + deriv(2)
 121       continue
 120    continue

        do 140 k=1,j
           nodes(k) = 0
 140    continue

 	if(v.ne.0.0)then
 	   f = f/v
 	   df(1) = (df(1) - f*dvx)/v
 	   df(2) = (df(2) - f*dvy)/v
 	end  if

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2Df - calculates natural neighbour interpolation at point x
c		when point is inside convex hull 
c		(using method based on a closed 2-D formula)
c
c	Input:
c		x(2)			co-ordinates of input point
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c	        loc			index of triangle containing 
c					input point.
c		data			data values at each node
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                        set to ~50 in calling program)
c		tstore			integer work array of size nnpn_max
c		index			integer work array of size nnpn_max
c		vpairs			integer work array of size 2*nnpn_max
c		angles			real work array of size nnpn_max
c		pverts			real*8 work array of size 2*nnpn_max
c		vverts			real*8 work array of size 2*nnpn_max
c		ltri			logical work array size of nt_max
c		lnode			logical work array size of nt_max
c
c	Output:
c		f			interpolated data value
c		v			area of voronoi cell about x
c
c	Comments:
c		 On input point x must be in triangle loc. 
c		 Note: The input triangle contains x and so it's 
c		 circumcircle must also contain x.
c
c		 In order to use the closed formula to
c		 calculate the nn co-ordinate of the input point with
c		 respect to its ith neighbour we must find the set of
c		 nodes that are both neighbours of i and x, where 
c		 the neighbours of i are determined before addition of 
c		 x. This is done using a LIFO stack to avoid an 
c		 expensive global search over all triangles.
c
c		 It is assumed that the logical arrays ltri(nt_max) 
c		 and lnode are initialized to false on input.
c
c		 This version does not calculate derivatives.
c
c		 Calls are made to: plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty
c				    indexx & second_v_area.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Df
     &  (x,points,vertices,neighbour,centres,loc,data,nnpn_max,
     &   pverts,vverts,tstore,vpairs,angles,index,ltri,lnode,f,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		x(2),dist,dx,dy
	real*8		f,v
        real*8          pverts(2,nnpn_max)
        real*8          vverts(2,nnpn_max)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tri,pos,tcount
	integer		tstore(nnpn_max)
	integer		index(nnpn_max)
	integer		vpairs(2,nnpn_max)
	logical		nnwrite
	logical*1	ltri(*)
	logical*1	lnode(*)
        real            angles(nnpn_max)
c
	integer		cyc(3)
        data		cyc/2,3,1/

        

	common/nnswitches/lud,nnwrite

c					Initialize variables
	f = 0.d0
	v = 0.d0
        tcount = 1

	if(nnwrite)write(lud,*)loc
c					plot input circum-triangle 
c					and circum-circle

 	call plot_tc(loc,points,vertices,centres)

c					record circum-triangle 
        tstore(tcount) = loc
        ltri(loc) = .true.
c
c					Find all other circumcircles 
c					containing input point using a 
c					directed walk
c
c					initialize stack
	call stackpairinit
c					put input triangle's neighbouring 
c					triangles on LIFO stack
c					together with position of input
c					triangle in their neighbour list

	i = neighbour(1,loc)
	j = neighbour(2,loc)
	k = neighbour(3,loc)
	if(i.ne.0)then
	   if(neighbour(1,i).eq.loc)then
              call pushpair(i,1)
	   else if(neighbour(2,i).eq.loc)then
              call pushpair(i,2)
	   else if(neighbour(3,i).eq.loc)then
              call pushpair(i,3)
           end if
	end if
	if(j.ne.0)then
	   if(neighbour(1,j).eq.loc)then
              call pushpair(j,1)
	   else if(neighbour(2,j).eq.loc)then
              call pushpair(j,2)
	   else if(neighbour(3,j).eq.loc)then
              call pushpair(j,3)
           end if
	end if
	if(k.ne.0)then
	   if(neighbour(1,k).eq.loc)then
              call pushpair(k,1)
	   else if(neighbour(2,k).eq.loc)then
              call pushpair(k,2)
	   else if(neighbour(3,k).eq.loc)then
              call pushpair(k,3)
           end if
	end if

 10	call stackpairempty(k)
c					if stack empty then finish
	if(k.eq.1)go to 100
c					take triangle from stack
	call poppair(tri,pos)

c					test if x is in circumcircle
        dx = x(1)-centres(1,tri)
        dy = x(2)-centres(2,tri)
	dist = dx*dx + dy*dy 

 	if(dist.lt.centres(3,tri))then
           if(nnwrite)write(lud,*)tri

c					plot circum-triangle and circum-circle

 	   call plot_tc(tri,points,vertices,centres)

	   if(pos.eq.1)then
	      new1 = neighbour(2,tri)
	      new2 = neighbour(3,tri)
	   else if(pos.eq.2)then
	      new1 = neighbour(1,tri)
	      new2 = neighbour(3,tri)
	   else 
	      new1 = neighbour(1,tri)
	      new2 = neighbour(2,tri)
	   end if
c					record current circum-triangle
c					and pseudo angle to x
           tcount = tcount + 1
           tstore(tcount) = tri
           ltri(tri) = .true.

	   if(new1.ne.0)then
	      if(neighbour(1,new1).eq.tri)then
                 call pushpair(new1,1)
	      else if(neighbour(2,new1).eq.tri)then
                 call pushpair(new1,2)
	      else if(neighbour(3,new1).eq.tri)then
                 call pushpair(new1,3)
              end if
	   end if
	   if(new2.ne.0)then
	      if(neighbour(1,new2).eq.tri)then
                 call pushpair(new2,1)
	      else if(neighbour(2,new2).eq.tri)then
                 call pushpair(new2,2)
	      else if(neighbour(3,new2).eq.tri)then
                 call pushpair(new2,3)
              end if
	   end if
	end if
	go to 10
 
 100    continue
        call stackpairflush()

c						check size of tcount
        if(tcount+2.gt.nnpn_max)then
           write(*,*)' '
           write(*,*)' Error: work arrays in subroutine nn2Df',
     &               ' not big enough.'
           write(*,*)' Too many circum-triangles for this node'
           write(*,*)' Remedy: Increase size of parameter'
           write(*,*)'         nnpn_max in calling program'
           write(*,*)'         current value = ',nnpn_max
           write(*,*)'         required value >= ',tcount+2
	   stop
        end if
c
c       call xsave
c
c
c						find list of vertices
c						of new voronoi cell about
c						x and store the pair of
c						nodes associated with each
c						vertex
c       write(*,*)' Circum-triangles of x:'
c
        nv = 0
        do 110 i=1,tcount
           it = tstore(i)
c          write(*,*)it
           do j=1,3
              node = vertices(j,it)
              j1 = cyc(j)
              j2 = cyc(j1)
              na = vertices(j1,it)
              nb = vertices(j2,it)
              iopp = neighbour(j,it)
              if(iopp.eq.0.or..not.ltri(iopp))then
                 nv = nv + 1
c
c						find circum-centre of
c						x, node, and nodeb 
c
                 call circum (x,points(1,na),
     &                points(1,nb),vverts(1,nv))

                 vpairs(1,nv) = na
                 vpairs(2,nv) = nb
                 angles(nv) = pangle(x,vverts(1,nv))
                  
              end if
           end do
           
 110    continue
        
c						sort vertices 
c
        call indexx(nv,angles,index)
c
c						debug I/O
c						
c
c       write(*,*)' Number of vertices = ',nv
c	do i=1,nv
c          write(*,*)' i',i,' index',index(i),' a',angles(i)
c       end do
c
c						find vertices of each 
c						second-order Voronoi cell
c
	ncount = 0
c						loop over neighbours of x
        do 120 it=1,tcount
	   tri = tstore(it)
           do 130 in = 1,3
              node = vertices(in,tri)
              if(lnode(node))go to 130
              lnode(node) = .true.
              nverts=0
c             write(*,*)' Neighbour node =',node
c
c						find pair of `external' vertices
c						of second-order voronoi cell  
c
              do i=1,nv
                 k = index(i)
                 if(vpairs(1,k).eq.node.or.
     &              vpairs(2,k).eq.node)then
                    nverts = nverts + 1
                    pverts(1,nverts) = vverts(1,k) - x(1)
                    pverts(2,nverts) = vverts(2,k) - x(2)
                    ip = i+1
                    im = i-1
                    if(i.eq.nv)ip=1
                    if(i.eq.1)im=nv
                    kp = index(ip)
                    km = index(im)
                    if(vpairs(1,kp).eq.node.or.
     &                 vpairs(2,kp).eq.node)then
                       nverts = nverts + 1
                       pverts(1,nverts) = vverts(1,kp) - x(1)
                       pverts(2,nverts) = vverts(2,kp) - x(2)
                    else if(vpairs(1,km).eq.node.or.
     &                      vpairs(2,km).eq.node)then
                       nverts = nverts + 1
                       pverts(1,1) = vverts(1,km) - x(1)
                       pverts(2,1) = vverts(2,km) - x(2)
                       pverts(1,2) = vverts(1,k) - x(1)
                       pverts(2,2) = vverts(2,k) - x(2)
                    else
                       write(*,*)' Error in nn2Df'
                       write(*,*)
     &                 ' Can not find neighbouring vertices' 
                       write(*,*)' node=',node,' i =',i
                    end if
                    go to 150
                 end if
              end do
 150          continue
c
c						find `internal' vertices of 
c						second-order voronoi cell  
	      do 140 jt = 1,tcount
                 jtri = tstore(jt)
                 do jn = 1,3
                    jnode = vertices(jn,jtri)
                    if(jnode.eq.node)then
                       nverts = nverts + 1
c						record circum-centre
c
                       pverts(1,nverts) = centres(1,jtri) - x(1)
                       pverts(2,nverts) = centres(2,jtri) - x(2)
                       go to 140
                    end if
                 end do
 140          continue
c						Calculate volume of 
c						second-order voronoi cell 
c						for current node using
c						direct 2-D formula 
c
              call second_v_area
     &             (x,nverts,pverts,angles,area)
c
           v = v + area
           f = f + area*data(node)

 130       continue
 120    continue
 	if(v.ne.0.0) f = f/v
        v = v*0.5
c
c						reset logical arrays
        do i=1,tcount
           tri = tstore(i)
           ltri(tri) = .false.
           do j = 1,3
              lnode(vertices(j,tri)) = .false.
           end do
        end do 

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2Dfd - calculates natural neighbour interpolation and 
c		 1st derivatives at point x when point is inside convex hull 
c		 (using method based on a closed 2-D formula)
c
c	Input:
c		x(2)			co-ordinates of input point
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c	        loc			index of triangle containing 
c					input point.
c		data			data values at each node
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                        set to ~50 in calling program)
c		tstore			integer work array of size nnpn_max
c		index			integer work array of size nnpn_max
c		vpairs			integer work array of size 2*nnpn_max
c		angles			real work array of size nnpn_max
c		pverts			real*8 work array of size 2*nnpn_max
c		vverts			real*8 work array of size 2*nnpn_max
c		dx_verts		real*8 work array of size 2*nnpn_max
c		dy_verts		real*8 work array of size 2*nnpn_max
c		ltri			logical work array of size nt_max
c		lnode			logical work array of size nt_max
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		v			area of voronoi cell about x
c
c	Comments:
c		 On input point x must be in triangle loc. 
c		 Note: The input triangle contains x and so it's 
c		 circumcircle must also contain x.
c
c		 In order to use the direct formula to
c		 calculate the nn co-ordinate of the input point with
c		 respect to its ith neighbour we must find the set of
c		 nodes that are both neighbours of i and x, where 
c		 the neighbours of i are determined before addition of x.
c		 This is done using a LIFO stack to avoid an 
c		 expensive global search over all triangles.
c
c		 It is assumed that the logical arrays ltri(nt_max) 
c		 and lnode are initialized to false on input.
c
c		 This version does not calculate derivatives.
c
c		 Calls are made to: plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty 
c				    indexx & second_v_area_d.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Dfd
     &  (x,points,vertices,neighbour,centres,loc,data,nnpn_max,
     &   pverts,vverts,tstore,vpairs,angles,index,dx_verts,dy_verts,
     &   ltri,lnode,f,df,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		x(2),dist,dx,dy
	real*8		f,v,df(2),deriv(2),dvx,dvy
        real*8          pverts(2,nnpn_max)
        real*8          vverts(2,nnpn_max)
        real*8          dx_verts(2,nnpn_max)
        real*8          dy_verts(2,nnpn_max)
        real*8          dpverts(4,2)
        real            angles(nnpn_max)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tri,pos,tcount
	integer		tstore(nnpn_max)
	integer		index(nnpn_max)
	integer		vpairs(2,nnpn_max)
	logical		nnwrite
	logical*1	ltri(*)
	logical*1	lnode(*)
c
	integer		cyc(3)
        data		cyc/2,3,1/

        

	common/nnswitches/lud,nnwrite

c					Initialize variables
	f = 0.d0
	df(1) = 0.d0
	df(2) = 0.d0
	v = 0.d0
	dvx = 0.d0
	dvy = 0.d0
        tcount = 1

	if(nnwrite)write(lud,*)loc
c					plot input circum-triangle 
c					and circum-circle

 	call plot_tc(loc,points,vertices,centres)

c					record circum-triangle 
        tstore(tcount) = loc
        ltri(loc) = .true.
c
c					Find all other circumcircles 
c					containing input point using a 
c					directed walk
c
c					initialize stack
	call stackpairinit


c					put input triangle's neighbouring 
c					triangles on LIFO stack
c					together with position of input
c					triangle in their neighbour list

	i = neighbour(1,loc)
	j = neighbour(2,loc)
	k = neighbour(3,loc)
	if(i.ne.0)then
	   if(neighbour(1,i).eq.loc)then
              call pushpair(i,1)
	   else if(neighbour(2,i).eq.loc)then
              call pushpair(i,2)
	   else if(neighbour(3,i).eq.loc)then
              call pushpair(i,3)
           end if
	end if
	if(j.ne.0)then
	   if(neighbour(1,j).eq.loc)then
              call pushpair(j,1)
	   else if(neighbour(2,j).eq.loc)then
              call pushpair(j,2)
	   else if(neighbour(3,j).eq.loc)then
              call pushpair(j,3)
           end if
	end if
	if(k.ne.0)then
	   if(neighbour(1,k).eq.loc)then
              call pushpair(k,1)
	   else if(neighbour(2,k).eq.loc)then
              call pushpair(k,2)
	   else if(neighbour(3,k).eq.loc)then
              call pushpair(k,3)
           end if
	end if

 10	call stackpairempty(k)
c					if stack empty then finish
	if(k.eq.1)go to 100
c					take triangle from stack
	call poppair(tri,pos)

c					test if x is in circumcircle
        dx = x(1)-centres(1,tri)
        dy = x(2)-centres(2,tri)
	dist = dx*dx + dy*dy 

 	if(dist.lt.centres(3,tri))then
           if(nnwrite)write(lud,*)tri

c					plot circum-triangle and circum-circle

 	   call plot_tc(tri,points,vertices,centres)

	   if(pos.eq.1)then
	      new1 = neighbour(2,tri)
	      new2 = neighbour(3,tri)
	   else if(pos.eq.2)then
	      new1 = neighbour(1,tri)
	      new2 = neighbour(3,tri)
	   else 
	      new1 = neighbour(1,tri)
	      new2 = neighbour(2,tri)
	   end if
c					record current circum-triangle
c					and pseudo angle to x
           tcount = tcount + 1
           tstore(tcount) = tri
           ltri(tri) = .true.

	   if(new1.ne.0)then
	      if(neighbour(1,new1).eq.tri)then
                 call pushpair(new1,1)
	      else if(neighbour(2,new1).eq.tri)then
                 call pushpair(new1,2)
	      else if(neighbour(3,new1).eq.tri)then
                 call pushpair(new1,3)
              end if
	   end if
	   if(new2.ne.0)then
	      if(neighbour(1,new2).eq.tri)then
                 call pushpair(new2,1)
	      else if(neighbour(2,new2).eq.tri)then
                 call pushpair(new2,2)
	      else if(neighbour(3,new2).eq.tri)then
                 call pushpair(new2,3)
              end if
	   end if
	end if
	go to 10
 
 100    continue
        call stackpairflush()

c						check size of tcount
        if(tcount+2.gt.nnpn_max)then
           write(*,*)' '
           write(*,*)' Error: work arrays in subroutine nn2Dfd',
     &               ' not big enough.'
           write(*,*)' Too many circum-triangles for this node'
           write(*,*)' Remedy: Increase size of parameter'
           write(*,*)'         nnpn_max in calling program'
           write(*,*)'         current value = ',nnpn_max
           write(*,*)'         required value >= ',tcount+2
	   stop
        end if
c						find list of vertices
c						of new voronoi cell about
c						x and store the pair of
c						nodes associated with each
c						vertex
c       write(*,*)' Circum-triangles of x:'
        nv = 0
        do 110 i=1,tcount
           it = tstore(i)
c          write(*,*)it
           do j=1,3
              node = vertices(j,it)
              j1 = cyc(j)
              j2 = cyc(j1)
              na = vertices(j1,it)
              nb = vertices(j2,it)
              iopp = neighbour(j,it)
              if(iopp.eq.0.or..not.ltri(iopp))then
                 nv = nv + 1
c
c						find circum-centre of
c						x, node, and nodeb 
c
                 call circum_d (x,points(1,na),points(1,nb),
     &                          vverts(1,nv),
     &                          dx_verts(1,nv),dy_verts(1,nv))

                 vpairs(1,nv) = na
                 vpairs(2,nv) = nb
                 angles(nv) = pangle(x,vverts(1,nv))
c
c                write(*,*)' na =',na,' nb=',nb
c                write(*,*)' dx_verts =',dx_verts(1,nv)
c                write(*,*)' dx_verts =',dx_verts(2,nv)
c                write(*,*)' dy_verts =',dy_verts(1,nv)
c                write(*,*)' dy_verts =',dy_verts(2,nv)
c
              end if
           end do
           
 110    continue
        
c						sort vertices 
c
        call indexx(nv,angles,index)
c						debug I/O
c       write(*,*)' Number of vertices = ',nv
c       write(*,*)' index '
c	do i=1,nv
c          write(*,*)' i',i,' index',index(i),' a',angles(i)
c       end do
c
c						find vertices of each 
c						second-order Voronoi cell
c
	ncount = 0
c						loop over neighbours of x
        do 120 it=1,tcount
	   tri = tstore(it)
           do 130 in = 1,3
              node = vertices(in,tri)
              if(lnode(node))go to 130
              lnode(node) = .true.
c             write(*,*)' Neighbour node =',node
c
c						find pair of `external' vertices
c						of second-order voronoi cell  
c
              do i=1,nv
                 k = index(i)
                 if(vpairs(1,k).eq.node.or.
     &              vpairs(2,k).eq.node)then
                    pverts(1,1) = vverts(1,k) - x(1)
                    pverts(2,1) = vverts(2,k) - x(2)
                    dpverts(1,1) = dx_verts(1,k)
                    dpverts(2,1) = dx_verts(2,k)
                    dpverts(3,1) = dy_verts(1,k)
                    dpverts(4,1) = dy_verts(2,k)
                    ip = i+1
                    im = i-1
                    if(i.eq.nv)ip=1
                    if(i.eq.1)im=nv
                    kp = index(ip)
                    km = index(im)
                    if(vpairs(1,kp).eq.node.or.
     &                 vpairs(2,kp).eq.node)then
                       pverts(1,2) = vverts(1,kp) - x(1)
                       pverts(2,2) = vverts(2,kp) - x(2)
                       dpverts(1,2) = dx_verts(1,kp)
                       dpverts(2,2) = dx_verts(2,kp)
                       dpverts(3,2) = dy_verts(1,kp)
                       dpverts(4,2) = dy_verts(2,kp)
                    else if(vpairs(1,km).eq.node.or.
     &                      vpairs(2,km).eq.node)then
                       pverts(1,1) = vverts(1,km) - x(1)
                       pverts(2,1) = vverts(2,km) - x(2)
                       pverts(1,2) = vverts(1,k) - x(1)
                       pverts(2,2) = vverts(2,k) - x(2)
                       dpverts(1,1) = dx_verts(1,km)
                       dpverts(2,1) = dx_verts(2,km)
                       dpverts(3,1) = dy_verts(1,km)
                       dpverts(4,1) = dy_verts(2,km)
                       dpverts(1,2) = dx_verts(1,k)
                       dpverts(2,2) = dx_verts(2,k)
                       dpverts(3,2) = dy_verts(1,k)
                       dpverts(4,2) = dy_verts(2,k)
                    else
                       write(*,*)' Error in nn2Dfd'
                       write(*,*)
     &                 ' Can not find neighbouring vertices' 
                       write(*,*)' node=',node,' i =',i
                    end if
c
                    go to 150
                 end if
              end do
 150          continue
              nverts = 2
c
c						find `internal' vertices of 
c						second-order voronoi cell  
	      do 140 jt = 1,tcount
                 jtri = tstore(jt)
                 do jn = 1,3
                    jnode = vertices(jn,jtri)
                    if(jnode.eq.node)then
                       nverts = nverts + 1
c
c						record circum-centre
c
                       pverts(1,nverts) = centres(1,jtri) - x(1)
                       pverts(2,nverts) = centres(2,jtri) - x(2)

c                      write(*,*)' found common triangle',jtri

                       go to 140
                    end if
                 end do
 140          continue
c						Calculate volume of 
c						second-order voronoi cell 
c						for current node using
c						direct 2-D formula 
c
              call second_v_area_d
     &        (x,nverts,pverts,dpverts,angles,area,deriv)
c
              v = v + area
              f = f + area*data(node)
              df(1) = df(1) + data(node)*deriv(1)
              df(2) = df(2) + data(node)*deriv(2)
              dvx = dvx + deriv(1)
              dvy = dvy + deriv(2)

 130       continue
 120    continue
 	if(v.ne.0.0)then
           f = f/v
           df(1) = (df(1) - f*dvx)/v
           df(2) = (df(2) - f*dvy)/v
        end if
        v = v*0.5
c
c						reset logical arrays
        do i=1,tcount
           tri = tstore(i)
           ltri(tri) = .false.
           do j = 1,3
              lnode(vertices(j,tri)) = .false.
           end do
        end do 

	return
	end
c
c------------------------------------------------------------------------
c
c	nn2Dfdd - calculates natural neighbour interpolation 
c		  and 1st & 2nd  derivatives at point x when 
c		  point is inside convex hull 
c		  (using method based on a closed 2-D formula)
c
c	Input:
c		x(2)			co-ordinates of input point
c		points(2,np)		array of node points	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,nt)		array of neighbouring triangles.	
c					Neighbour(i,j) is the triangle
c					opposite node i in triangle j,
c					stored counterclockwise about j.
c               centres(3,nt)           centres(j,i) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circumcircle about Delaunay 
c                                       triangle i, (i=1,...,nt),
c                                       j=3 contains squared radius of circle.
c	        loc			index of triangle containing 
c					input point.
c		data			data values at each node
c		nnpn_max		maximum number of neighbours per node
c					(depends on the point distribution,
c                                        set to ~50 in calling program)
c		tstore			integer work array of size nnpn_max
c		index			integer work array of size nnpn_max
c		vpairs			integer work array of size 2*nnpn_max
c		angles			real work array of size nnpn_max
c		pverts			real*8 work array of size 2*nnpn_max
c		vverts			real*8 work array of size 2*nnpn_max
c		dx_verts		real*8 work array of size 2*nnpn_max
c		dy_verts		real*8 work array of size 2*nnpn_max
c		dxx_verts		real*8 work array of size 2*nnpn_max
c		dyy_verts		real*8 work array of size 2*nnpn_max
c		dxy_verts		real*8 work array of size 2*nnpn_max
c		ltri			logical work array of size nt_max
c		lnode			logical work array of size nt_max
c
c	Output:
c		f			interpolated data value
c		df(1)			interpolated derivative df/dx
c		df(2)			interpolated derivative df/dy
c		ddf(1)			interpolated derivative d2f/dxx
c		ddf(2)			interpolated derivative d2f/dyy
c		ddf(3)			interpolated derivative d2f/dxy
c		v			area of voronoi cell about x
c
c	Comments:
c		 On input point x must be in triangle loc. 
c		 Note: The input triangle contains x and so it's 
c		 circumcircle must also contain x.
c
c		 In order to use the direct formula to
c		 calculate the nn co-ordinate of the input point with
c		 respect to its ith neighbour we must find the set of
c		 nodes that are both neighbours of i and x, where 
c		 the neighbours of i are determined before addition of x.
c		 This is done using a LIFO stack to avoid an 
c		 expensive global search over all triangles.
c
c		 It is assumed that the logical arrays ltri(nt_max) 
c		 and lnode are initialized to false on input.
c
c		 This version does not calculate derivatives.
c
c		 Calls are made to: plot_tc, stackpairinit,
c				    poppair,pushpair & stackpairempty 
c				    indexx & second_v_area_dd.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
	Subroutine nn2Dfdd
     &  (x,points,vertices,neighbour,centres,loc,data,nnpn_max,
     &   pverts,vverts,tstore,vpairs,angles,index,dx_verts,dy_verts,
     &   dxx_verts,dyy_verts,dxy_verts,ltri,lnode,f,df,ddf,v)

	real*8		points(2,*)
	real*8		centres(3,*)
	real*8		data(*)
	real*8		x(2),dist,dx,dy
	real*8		f,v,df(2),ddf(3)
        real*8		deriv(2),dderiv(3),dvx,dvy,dvxx,dvyy,dvxy
        real*8          pverts(2,nnpn_max)
        real*8          dpverts(4,2)
        real*8          ddpverts(6,2)
        real*8          vverts(2,nnpn_max)
        real*8          dx_verts(2,nnpn_max)
        real*8          dy_verts(2,nnpn_max)
        real*8          dxx_verts(2,nnpn_max)
        real*8          dyy_verts(2,nnpn_max)
        real*8          dxy_verts(2,nnpn_max)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		tri,pos,tcount
	integer		tstore(nnpn_max)
	integer		index(nnpn_max)
	integer		vpairs(2,nnpn_max)
	logical		nnwrite
	logical*1	ltri(*)
	logical*1	lnode(*)
        real            angles(nnpn_max)
	integer		cyc(3)
        data		cyc/2,3,1/

	common/nnswitches/lud,nnwrite
c					Initialize variables
	f = 0.d0
	df(1) = 0.d0
	df(2) = 0.d0
	ddf(1) = 0.d0
	ddf(2) = 0.d0
	ddf(3) = 0.d0
	v = 0.d0
	dvx = 0.d0
	dvy = 0.d0
	dvxx = 0.d0
	dvyy = 0.d0
	dvxy = 0.d0
        tcount = 1

	if(nnwrite)write(lud,*)loc
c					plot input circum-triangle 
c					and circum-circle

 	call plot_tc(loc,points,vertices,centres)

c					record circum-triangle 
        tstore(tcount) = loc
        ltri(loc) = .true.
c
c					Find all other circumcircles 
c					containing input point using a 
c					directed walk
c
c					initialize stack
	call stackpairinit


c					put input triangle's neighbouring 
c					triangles on LIFO stack
c					together with position of input
c					triangle in their neighbour list

	i = neighbour(1,loc)
	j = neighbour(2,loc)
	k = neighbour(3,loc)
	if(i.ne.0)then
	   if(neighbour(1,i).eq.loc)then
              call pushpair(i,1)
	   else if(neighbour(2,i).eq.loc)then
              call pushpair(i,2)
	   else if(neighbour(3,i).eq.loc)then
              call pushpair(i,3)
           end if
	end if
	if(j.ne.0)then
	   if(neighbour(1,j).eq.loc)then
              call pushpair(j,1)
	   else if(neighbour(2,j).eq.loc)then
              call pushpair(j,2)
	   else if(neighbour(3,j).eq.loc)then
              call pushpair(j,3)
           end if
	end if
	if(k.ne.0)then
	   if(neighbour(1,k).eq.loc)then
              call pushpair(k,1)
	   else if(neighbour(2,k).eq.loc)then
              call pushpair(k,2)
	   else if(neighbour(3,k).eq.loc)then
              call pushpair(k,3)
           end if
	end if

 10	call stackpairempty(k)
c					if stack empty then finish
	if(k.eq.1)go to 100
c					take triangle from stack
	call poppair(tri,pos)

c					test if x is in circumcircle
        dx = x(1)-centres(1,tri)
        dy = x(2)-centres(2,tri)
	dist = dx*dx + dy*dy 

 	if(dist.lt.centres(3,tri))then
           if(nnwrite)write(lud,*)tri

c					plot circum-triangle and circum-circle

 	   call plot_tc(tri,points,vertices,centres)

	   if(pos.eq.1)then
	      new1 = neighbour(2,tri)
	      new2 = neighbour(3,tri)
	   else if(pos.eq.2)then
	      new1 = neighbour(1,tri)
	      new2 = neighbour(3,tri)
	   else 
	      new1 = neighbour(1,tri)
	      new2 = neighbour(2,tri)
	   end if
c					record current circum-triangle
c					and pseudo angle to x
           tcount = tcount + 1
           tstore(tcount) = tri
           ltri(tri) = .true.

	   if(new1.ne.0)then
	      if(neighbour(1,new1).eq.tri)then
                 call pushpair(new1,1)
	      else if(neighbour(2,new1).eq.tri)then
                 call pushpair(new1,2)
	      else if(neighbour(3,new1).eq.tri)then
                 call pushpair(new1,3)
              end if
	   end if
	   if(new2.ne.0)then
	      if(neighbour(1,new2).eq.tri)then
                 call pushpair(new2,1)
	      else if(neighbour(2,new2).eq.tri)then
                 call pushpair(new2,2)
	      else if(neighbour(3,new2).eq.tri)then
                 call pushpair(new2,3)
              end if
	   end if
	end if
	go to 10
 
 100    continue
        call stackpairflush()

c						check size of tcount
        if(tcount+2.gt.nnpn_max)then
           write(*,*)' '
           write(*,*)' Error: work arrays in subroutine nn2Dfdd',
     &               ' not big enough.'
           write(*,*)' Too many circum-triangles for this node'
           write(*,*)' Remedy: Increase size of parameter'
           write(*,*)'         nnpn_max in calling program'
           write(*,*)'         current value = ',nnpn_max
           write(*,*)'         required value >= ',tcount+2
	   stop
        end if
c
c						find list of vertices
c						of new voronoi cell about
c						x and store the pair of
c						nodes associated with each
c						vertex
c       write(*,*)' Circum-triangles of x:'
        nv = 0
        do 110 i=1,tcount
           it = tstore(i)
c          write(*,*)it
           do j=1,3
              node = vertices(j,it)
              j1 = cyc(j)
              j2 = cyc(j1)
              na = vertices(j1,it)
              nb = vertices(j2,it)
              iopp = neighbour(j,it)
              if(iopp.eq.0.or..not.ltri(iopp))then
                 nv = nv + 1
c
c						find circum-centre of
c						x, node, and nodeb 
c
                 call circum_dd (x,points(1,na),points(1,nb),
     &                          vverts(1,nv),
     &                          dx_verts(1,nv),dy_verts(1,nv),
     &                          dxx_verts(1,nv),dyy_verts(1,nv),
     &                          dxy_verts(1,nv))

                 vpairs(1,nv) = na
                 vpairs(2,nv) = nb
                 angles(nv) = pangle(x,vverts(1,nv))
c                write(*,*)' na =',na,' nb=',nb
 
              end if
           end do
           
 110    continue
        
c						sort vertices 
c
        call indexx(nv,angles,index)
c						debug I/O
c       write(*,*)' Number of vertices = ',nv
c	do i=1,nv
c          write(*,*)' i',i,' index',index(i),' a',angles(i)
c       end do
c
c						find vertices of each 
c						second-order Voronoi cell
c
	ncount = 0
c						loop over neighbours of x
        do 120 it=1,tcount
	   tri = tstore(it)
           do 130 in = 1,3
              node = vertices(in,tri)
              if(lnode(node))go to 130
              lnode(node) = .true.
c             write(*,*)' Neighbour node =',node
c
c						find pair of `external' vertices
c						of second-order voronoi cell  
c
              do i=1,nv
                 k = index(i)
                 if(vpairs(1,k).eq.node.or.
     &              vpairs(2,k).eq.node)then
                    pverts(1,1) = vverts(1,k) - x(1)
                    pverts(2,1) = vverts(2,k) - x(2)
                    dpverts(1,1) = dx_verts(1,k)
                    dpverts(2,1) = dx_verts(2,k)
                    dpverts(3,1) = dy_verts(1,k)
                    dpverts(4,1) = dy_verts(2,k)
                    ddpverts(1,1) = dxx_verts(1,k)
                    ddpverts(2,1) = dxx_verts(2,k)
                    ddpverts(3,1) = dyy_verts(1,k)
                    ddpverts(4,1) = dyy_verts(2,k)
                    ddpverts(5,1) = dxy_verts(1,k)
                    ddpverts(6,1) = dxy_verts(2,k)
                    ip = i+1
                    im = i-1
                    if(i.eq.nv)ip=1
                    if(i.eq.1)im=nv
                    kp = index(ip)
                    km = index(im)
                    if(vpairs(1,kp).eq.node.or.
     &                 vpairs(2,kp).eq.node)then
                       pverts(1,2) = vverts(1,kp) - x(1)
                       pverts(2,2) = vverts(2,kp) - x(2)
                       dpverts(1,2) = dx_verts(1,kp)
                       dpverts(2,2) = dx_verts(2,kp)
                       dpverts(3,2) = dy_verts(1,kp)
                       dpverts(4,2) = dy_verts(2,kp)
                       ddpverts(1,2) = dxx_verts(1,kp)
                       ddpverts(2,2) = dxx_verts(2,kp)
                       ddpverts(3,2) = dyy_verts(1,kp)
                       ddpverts(4,2) = dyy_verts(2,kp)
                       ddpverts(5,2) = dxy_verts(1,kp)
                       ddpverts(6,2) = dxy_verts(2,kp)
                    else if(vpairs(1,km).eq.node.or.
     &                      vpairs(2,km).eq.node)then
                       pverts(1,1) = vverts(1,km) - x(1)
                       pverts(2,1) = vverts(2,km) - x(2)
                       pverts(1,2) = vverts(1,k) - x(1)
                       pverts(2,2) = vverts(2,k) - x(2)
                       dpverts(1,1) = dx_verts(1,km)
                       dpverts(2,1) = dx_verts(2,km)
                       dpverts(3,1) = dy_verts(1,km)
                       dpverts(4,1) = dy_verts(2,km)
                       dpverts(1,2) = dx_verts(1,k)
                       dpverts(2,2) = dx_verts(2,k)
                       dpverts(3,2) = dy_verts(1,k)
                       dpverts(4,2) = dy_verts(2,k)
                       ddpverts(1,1) = dxx_verts(1,km)
                       ddpverts(2,1) = dxx_verts(2,km)
                       ddpverts(3,1) = dyy_verts(1,km)
                       ddpverts(4,1) = dyy_verts(2,km)
                       ddpverts(5,1) = dxy_verts(1,km)
                       ddpverts(6,1) = dxy_verts(2,km)
                       ddpverts(1,2) = dxx_verts(1,k)
                       ddpverts(2,2) = dxx_verts(2,k)
                       ddpverts(3,2) = dyy_verts(1,k)
                       ddpverts(4,2) = dyy_verts(2,k)
                       ddpverts(5,2) = dxy_verts(1,k)
                       ddpverts(6,2) = dxy_verts(2,k)
                    else
                       write(*,*)' Error in nn2Dfdd'
                       write(*,*)
     &                 ' Can not find neighbouring vertices' 
                       write(*,*)' node=',node,' i =',i
                    end if
c
                    go to 150
                 end if
              end do
 150          continue
              nverts = 2
c
c						find `internal' vertices of 
c						second-order voronoi cell  
	      do 140 jt = 1,tcount
                 jtri = tstore(jt)
                 do jn = 1,3
                    jnode = vertices(jn,jtri)
                    if(jnode.eq.node)then
                       nverts = nverts + 1
c
c						record circum-centre
c
                       pverts(1,nverts) = centres(1,jtri) - x(1)
                       pverts(2,nverts) = centres(2,jtri) - x(2)

c                      write(*,*)' found common triangle',jtri

                       go to 140
                    end if
                 end do
 140          continue
c						Calculate volume of 
c						second-order voronoi cell 
c						for current node using
c						direct 2-D formula 
c
              call second_v_area_dd
     &             (x,nverts,pverts,dpverts,ddpverts,
     &              angles,area,deriv,dderiv)
c
              v = v + area
              f = f + area*data(node)
              df(1) = df(1) + data(node)*deriv(1)
              df(2) = df(2) + data(node)*deriv(2)
              dvx = dvx + deriv(1)
              dvy = dvy + deriv(2)
              dvxx = dvxx + dderiv(1)
              dvyy = dvyy + dderiv(2)
              dvxy = dvxy + dderiv(3)
              ddf(1) = ddf(1) + data(node)*dderiv(1)
              ddf(2) = ddf(2) + data(node)*dderiv(2)
              ddf(3) = ddf(3) + data(node)*dderiv(3)

c             write(*,*)' Second order Voronoi-cell between x and ',node,
c    &                  ' has ',nverts,' vertices'
 130       continue
 120    continue
 	if(v.ne.0.0)then
           f = f/v
           df(1) = (df(1) - f*dvx)/v
           df(2) = (df(2) - f*dvy)/v
           ddf(1) = (ddf(1) - 2.*dvx*df(1) - dvxx*f)/v
           ddf(2) = (ddf(2) - 2.*dvy*df(2) - dvyy*f)/v
           ddf(3) = (ddf(3) - dvy*df(1) - dvx*df(2) - dvxy*f)/v
        end if
        v = v*0.5
c
c						reset logical arrays
        do i=1,tcount
           tri = tstore(i)
           ltri(tri) = .false.
           do j = 1,3
              lnode(vertices(j,tri)) = .false.
           end do
        end do 

	return
	end
c
c------------------------------------------------------------------------
c
c	first_voronoi - calculates the area or volume of the voronoi	
c			polygon (polyhedron) about point x.
c
c	Input:
c		x(n)			input point	
c		p(n,m)			m neighbouring points about x
c		n			number of dimensions
c		m			number of neighbours
c		nmax			maximum number of dimensions
c		mmax			maximum number of neighbours
c		a(m,n)			work array 
c		b(m)			work array 
c
c	Output:
c		vol			volume (n=3) or area (n=2) of 
c					voronoi cell about x.
c
c	Comments:
c		 The first-order voronoi cell is defined by the perpendicular
c		 bisectors of x and each of its neighbours p(n,i).
c		 Each bisector forms a linear inequality constraint and 
c		 the voronoi cell is the region which satisfies all 
c		 constraints.
c
c		 This routine determines the matrix a(i,j) and vector b(j)   
c		 of the linear system and calls the recursive C routine
c		 `volume' to calculate the volume (area) of the voronoi cell.
c
c
c					J. Braun,
c					& M. Sambridge, RSES, 1995.
c
c------------------------------------------------------------------------
c
      Subroutine first_voronoi (x,p,n,m,mmax,nmax,a,b,vol)

      real  x(nmax),p(nmax,mmax)
      real  a(mmax,nmax),b(mmax)

c					Set up coefficients 
c					of linear system

      do 5 j=1,m
         do 6 i=1,n
            a(j,i)=p(i,j)-x(i)
 6       continue
         b(j)=0.
         do 7 i=1,n
            b(j)=b(j)+p(i,j)*p(i,j)-x(i)*x(i)
 7       continue
         b(j)=b(j)/2.
 5    continue

      call volume (a,b,m,n,mmax,nmax,vol)

      if(vol.eq.-1)then
         write(*,*)' Warning in routine first_voronoi'
         write(*,*)' Voronoi cell is unbounded'
         write(*,*)' volume is set to zero'
         vol = 0.
      else if(vol.eq.0.)then
         write(*,*)' Warning in routine first_voronoi'
         write(*,*)' An inconsistent set of constraints has been'
         write(*,*)' detected in calculating the volume of '
         write(*,*)' a Voronoi cell.'
         write(*,*)' volume is set to zero'
      end if

      return
      end
c
c------------------------------------------------------------------------
c
c	second_voronoi - calculates the area or volume of the second-order
c			 voronoi polygon (polyhedron) between x and its
c			 neighbour p(n,1).
c
c	Input:
c		x(n)			input point	
c		p(n,m)			neighbours shared by x and p(n,1)
c		n			number of dimensions
c		m			number of neighbours
c		nmax			maximum number of dimensions
c		mmax			maximum number of neighbours
c		xo(n)			a local origin point.
c		a(m,n)			work array 
c		b(m)			work array 
c
c	Output:
c		vol			volume (n=3) or area (n=2) of 
c					voronoi cell about x.
c
c	Comments:
c		 The second-order voronoi cell between x and its neighbour p 
c		 is defined by the perpendicular bisector between x and p, 
c		 and the bisectors between p and the nodes which are neighbours 
c		 of both p and x. Each bisector forms a linear inequality 
c		 constraint and the second-order voronoi cell is the region 
c		 which satisfies all constraints.
c
c		 This routine determines the matrix a(i,j) and vector b(j)   
c		 of the linear system and calls the recursive C routine
c		 `volumef' to calculate the volume (area) of the voronoi cell.
c
c		 The point p is stored in the first row of input array p 
c		 (p(i,1),i=1,n). xo is any `nearby' point used as the origin
c		 for the calculation. It's use is only to help reduce
c		 roundoff error. For perfect arithmetic it will not affect  
c		 the calculated voronoi volume.
c
c					J. Braun, 
c					& M. Sambridge, RSES, 1995.
c
c------------------------------------------------------------------------
c
      Subroutine second_voronoi (x,p,n,m,mmax,nmax,xo,a,b,vol)

      real*8	x(nmax)
      real*8	p(nmax,mmax)
      real*8    xo(nmax)
      real  	a(mmax,nmax)
      real	b(mmax)
      real*8	dx1,dx2

c					Set up coefficients 
c					of linear system

      do 5 i=1,n
         a(1,i)=sngl(p(i,1)-x(i))
 5    continue
      b(1)=0.
      do 6 i=1,n
         dx1 = p(i,1)-xo(i)
         dx2 = x(i)-xo(i)
         b(1)=b(1) + sngl(dx1*dx1 - dx2*dx2)
 6    continue
      b(1)=b(1)/2.

      do 7 j=2,m
         do 8 i=1,n
            a(j,i)=sngl(p(i,j)-p(i,1))
 8       continue 
         b(j)=0.
         do 9 i=1,n
            dx1 = p(i,j)-xo(i)
            dx2 = p(i,1)-xo(i)
            b(j)=b(j) + sngl(dx1*dx1 - dx2*dx2)
 9       continue
         b(j)=b(j)/2.
 7    continue


      call volumef (a,b,m,n,mmax,nmax,vol)

      if(vol.eq.-1)then
         write(*,*)' Warning in routine second_voronoi'
         write(*,*)' Second order Voronoi cell is unbounded'
         write(*,*)' volume is set to zero'
         vol = 0.
      else if(vol.eq.0.)then
         write(*,*)' Warning in routine second_voronoi'
         write(*,*)' An inconsistent set of constraints has been'
         write(*,*)' detected in calculating the volume of '
         write(*,*)' the second-order Voronoi cell.'
         write(*,*)' volume is set to zero'
      end if

      return
      end
c
c------------------------------------------------------------------------
c
c	second_voronoi_d - calculates the area or volume of the second-order
c		 	   voronoi polygon (polyhedron) between x and its
c			   neighbour p(n,1). It also calculates the
c			   derivatives of the area (volume) with respect to
c			   the co-ordinates of x. 
c
c	Input:
c		x(n)			input point	
c		p(n,m)			neighbours shared by x and p(n,1)
c		n			number of dimensions
c		m			number of neighbours
c		nmax			maximum number of dimensions
c		mmax			maximum number of neighbours
c		xo(n)			a local origin point.
c		a(m,n)			work array 
c		b(m)			work array 
c
c	Output:
c		vol			volume (n=3) or area (n=2) of 
c					voronoi cell about x.
c		deriv(n)		first derivatives dv/dx_i (i=1,n)
c
c	Comments:
c		 See comments for routine second_voronoi.
c
c		 Calls C routines volumef, volumeb and dvdaf.
c
c					J. Braun, 
c					& M. Sambridge, RSES, 1995.
c
c------------------------------------------------------------------------
c
      subroutine second_voronoi_d (x,p,n,m,mmax,nmax,xo,a,b,vol,deriv)

      real*8	x(nmax)
      real*8	p(nmax,mmax)
      real*8    xo(nmax)
      real  	a(mmax,nmax)
      real	b(mmax)
      real	deriv(nmax)
      real*8    dx1,dx2

c					Set up coefficients 
c					of linear system

      do 5 i=1,n
         a(1,i)=sngl(p(i,1)-x(i))
 5    continue
      b(1)=0.
      do 6 i=1,n
         dx1 = p(i,1)-xo(i)
         dx2 = x(i)-xo(i)
         b(1)=b(1) + sngl(dx1*dx1 - dx2*dx2)
 6    continue
      b(1)=b(1)/2.

      do 7 j=2,m
         do 8 i=1,n
            a(j,i)=sngl(p(i,j)-p(i,1))
 8       continue 
         b(j)=0.
         do 9 i=1,n
            dx1 = p(i,j)-xo(i)
            dx2 = p(i,1)-xo(i)
            b(j)=b(j) + sngl(dx1*dx1 - dx2*dx2)
 9       continue
         b(j)=b(j)/2.
 7    continue

c     write(6,*)' A and b'
c     do 100 j=1,m
c        write(6,*)(a(j,i),i=1,n),' : ',b(j)
c100  continue
c     write(6,*)' '
c						calculate volume and
c						derivative d(vol)/db(1)

      call volumeb (a,b,m,n,mmax,nmax,1,vol,dvdb)

      if(vol.eq.-1)then
         write(*,*)' Warning in routine second_voronoi_d'
         write(*,*)' Second order Voronoi cell is unbounded'
         write(*,*)' volume and derivatives set to zero'
         vol = 0.
         do 10 k=1,n
            deriv(k) = 0.
 10      continue
      else if(vol.eq.0.)then
         write(*,*)' Warning in routine second_voronoi_d'
         write(*,*)' An inconsistent set of constraints has been'
         write(*,*)' detected in calculating the volume of '
         write(*,*)' the second-order Voronoi cell.'
         write(*,*)' volume and derivatives set to zero'
         do 11 k=1,n
            deriv(k) = 0.
 11      continue
      else
c						calculate derivatives
c						d(vol)/da(1,j), (j=1,n),
c						and use formulae to find
c						d(vol)/dx(k).
         do 20 k=1,n
   
           call dvdaf (a,b,m,n,mmax,nmax,k,deriva,icode)

           if(icode.eq.-1)then
              write(*,*)' Warning in routine second_voronoi_d'
              write(*,*)' Second order Voronoi cell is unbounded'
              write(*,*)' derivatives set to zero'
           else if(icode.eq.1)then
              write(*,*)' Warning in routine second_voronoi_d'
              write(*,*)' An inconsistent set of constraints has been'
              write(*,*)' detected in calculating the volume of '
              write(*,*)' the second-order Voronoi cell.'
              write(*,*)' All derivatives set to zero'
              do 30 i=1,n
                 deriv(i) = 0.
 30           continue
           end if

           deriv(k) = -1.*(dvdb*(x(k)-xo(k)) + deriva)
 20      continue 
      end if

      return
      end
c
c------------------------------------------------------------------------
c
c	Heapsort - Sorts an array into ascending order.
c		   Modified from numerical recipes 2, to use a 
c		   two dimensional double precision array.
c
c		   Sorting is done on index ka (ka=1 or 2).
c
c------------------------------------------------------------------------
c
c
      SUBROUTINE hpsort_d(n,ka,ra)
      INTEGER n
c     REAL ra(n)
      REAL*8 ra(2,n)
      INTEGER i,ir,j,l,ka
c     REAL rra
      REAL*8 rra(2)
      kb = 1
      if(ka.eq.1)kb = 2
      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra(ka)=ra(ka,l)
          rra(kb)=ra(kb,l)
        else
          rra(ka)=ra(ka,ir)
          rra(kb)=ra(kb,ir)
          ra(ka,ir)=ra(ka,1)
          ra(kb,ir)=ra(kb,1)
          ir=ir-1
          if(ir.eq.1)then
            ra(ka,1)=rra(ka)
            ra(kb,1)=rra(kb)
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(ka,j).lt.ra(ka,j+1))j=j+1
          endif
          if(rra(ka).lt.ra(ka,j))then
            ra(ka,i)=ra(ka,j)
            ra(kb,i)=ra(kb,j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ra(ka,i)=rra(ka)
        ra(kb,i)=rra(kb)
      goto 10
      END
c
c------------------------------------------------------------------------
c
c	pangle - pseudo angle routine 
c
c		 returns a number between 0 and 360 which is NOT the
c		 angle made by the line from p1 to p2 with the horizontal
c		 but which has the same order properties as that angle,
c		 i.e. has the same order of angles as arctan dy/dx.
c	 	 This function involves only simple products and quotients.
c
c		 From Sedgewick (1990) `Algorithms in C' (Addison Wesley)
c
c						M. Sambridge 1996.
c
c------------------------------------------------------------------------
c
	Function pangle(p1,p2)

	real*8		p1(2)
	real*8		p2(2)

	dx = p2(1) - p1(1)
        ax = abs(dx)
	dy = p2(2) - p1(2)
        ay = abs(dy)
        t = 0.
        a = ax+ay
	if(a.ne.0.)then
           t = dy/a
        end if
        if(dx.lt.0.)then
          t = 2-t
        else if(dy.lt.0.)then 
          t = 4+t
        end if
        pangle = t*90

	return
	end
c
c------------------------------------------------------------------------
c
c	Circum - calculates circum-centre of three points
c
c	Input:
c		pa,pb,pc		array of input points	
c
c	Output:
c               centre(3)               centre(j) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circle passing through the
c					three input points. 
c	Comments:
c
c		 Solves 3x3 linear system of equations.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
	Subroutine circum(pa,pb,pc,centre)
c
	real*8		pa(2),pb(2),pc(2)
	real*8		centre(2)
	real*8		x1,x2,x3,y1,y2,y3
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
c						Find centre of circum-circle
	   x1 = pa(1)
	   x2 = pb(1)
	   x3 = pc(1)
	   y1 = pa(2)
	   y2 = pb(2)
	   y3 = pc(2)

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1

	   centre(1) = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*dy2m1)/
     &         (denom)

	   centre(2) = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1) 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1))/ 
     &         (denom)
c
           centre(1) = centre(1)*0.5d0
           centre(2) = centre(2)*0.5d0

	return
	end
c
c------------------------------------------------------------------------
c
c	Circum_d - calculates circum-centre of three points and 
c	           derivatives of circum-centre with respect to 
c		   co-ordinates of first point.
c
c	Input:
c		pa,pb,pc		array of input points	
c
c	Output:
c               centre(3)               centre(j) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circle passing through the
c					three input points. 
c		vx(2)			derivative of centre with respect
c					to x-component of pa
c		vy(2)			derivative of centre with respect
c					to y-component of pa
c	Comments:
c
c		 Solves 3x3 linear system of equations and calculates
c		 derivatives of circum-centre with respect to co-ordinates
c		 of input vector pa(2).
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
	Subroutine circum_d(pa,pb,pc,centre,vx,vy)
c
	real*8		pa(2),pb(2),pc(2)
	real*8		centre(2)
	real*8		x1,x2,x3,y1,y2,y3
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
	real*8		vx(2),vy(2)
c						Find centre of circum-circle
	   x1 = pa(1)
	   x2 = pb(1)
	   x3 = pc(1)
	   y1 = pa(2)
	   y2 = pb(2)
	   y3 = pc(2)

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1

	   centre(1) = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*dy2m1)/
     &         (denom)

	   centre(2) = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1) 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1))/ 
     &         (denom)
c
           centre(1) = centre(1)*0.5d0
           centre(2) = centre(2)*0.5d0
c						X-derivative
c
           dcxm1 = centre(1) - x1
           dcym1 = centre(2) - y1

           denum1 = (dy3m1 - dy2m1)/denom
           denum2 = (dx2m1 - dx3m1)/denom

	   vx(1) = dcxm1*denum1
	   vx(2) = dcxm1*denum2

c						Y-derivative
c
	   vy(1) = dcym1*denum1
	   vy(2) = dcym1*denum2

	return
	end
c
c------------------------------------------------------------------------
c
c	Circum_dd - calculates circum-centre of three points and 
c		    1st and 2nd derivatives of circum-centre with 
c		    respect to co-ordinates of first point.
c
c	Input:
c		pa,pb,pc		array of input points	
c
c	Output:
c               centre(3)               centre(j) (j=1,2) contains the
c                                       co-ordinates of the centre of 
c                                       circle passing through the
c					three input points. 
c		vx(2)			derivative of centre with respect
c					to x-component of pa
c		vy(2)			derivative of centre with respect
c					to y-component of pa
c	Comments:
c
c		 Solves 3x3 linear system of equations and calculates
c		 derivatives of circum-centre with respect to co-ordinates
c		 of input vector pa(2).
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
	Subroutine circum_dd(pa,pb,pc,centre,vx,vy,vxx,vyy,vxy)
c
	real*8		pa(2),pb(2),pc(2)
	real*8		centre(2)
	real*8		x1,x2,x3,y1,y2,y3
	real*8		dx2m1,dx2p1,dy2m1,dy2p1
	real*8		dx3m1,dx3p1,dy3m1,dy3p1
	real*8		denom
	real*8		vx(2),vy(2),vxx(2),vyy(2),vxy(2)
c
c						Find centre of circum-circle
	   x1 = pa(1)
	   x2 = pb(1)
	   x3 = pc(1)
	   y1 = pa(2)
	   y2 = pb(2)
	   y3 = pc(2)

           dx2m1 = x2-x1
           dx2p1 = x2+x1
           dy2m1 = y2-y1
           dy2p1 = y2+y1
           dx3m1 = x3-x1
           dx3p1 = x3+x1
           dy3m1 = y3-y1
           dy3p1 = y3+y1
           denom = dx2m1*dy3m1-dx3m1*dy2m1

	   centre(1) = ((dx2m1*dx2p1 + dy2m1*dy2p1)*dy3m1
     &         -(dx3m1*dx3p1 + dy3m1*dy3p1)*dy2m1)/
     &         (denom)

	   centre(2) = (dx2m1*(dx3m1*dx3p1 + dy3m1*dy3p1) 
     &         -dx3m1*(dx2m1*dx2p1 + dy2m1*dy2p1))/ 
     &         (denom)
c
           centre(1) = centre(1)*0.5d0
           centre(2) = centre(2)*0.5d0
c						X-derivative
c
           dcxm1 = centre(1) - x1
           dcym1 = centre(2) - y1

           denum1 = (dy3m1 - dy2m1)/denom
           denum2 = (dx2m1 - dx3m1)/denom

	   vx(1) = dcxm1*denum1
	   vx(2) = dcxm1*denum2
c						Y-derivative
c
	   vy(1) = dcym1*denum1
	   vy(2) = dcym1*denum2
c						Second derivatives
	   f11 = 2*vx(1) - 1.
	   f22 = 2*vy(2) - 1.
	   f12 = vy(1) + vx(2)

           vxx(1) = f11*denum1
           vxx(2) = f11*denum2
           vyy(1) = f22*denum1
           vyy(2) = f22*denum2
           vxy(1) = f12*denum1
           vxy(2) = f12*denum2

	return
	end
c
c------------------------------------------------------------------------
c
c						heapsort modified to 
c						sort two arrays  
c
c------------------------------------------------------------------------
c
      SUBROUTINE hpsort_two(n,ra,rb)
      INTEGER n
      REAL ra(n)
      REAL*8 rb(2,n)
      INTEGER i,ir,j,l
      REAL rra
      REAL*8 rrb1,rrb2
      if (n.lt.2) return
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          rra=ra(l)
          rrb1=rb(1,l)
          rrb2=rb(2,l)
        else
          rra=ra(ir)
          rrb1=rb(1,ir)
          rrb2=rb(2,ir)
          ra(ir)=ra(1)
          rb(1,ir)=rb(1,1)
          rb(2,ir)=rb(2,1)
          ir=ir-1
          if(ir.eq.1)then
            ra(1)=rra
            rb(1,1)=rrb1
            rb(2,1)=rrb2
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(ra(j).lt.ra(j+1))j=j+1
          endif
          if(rra.lt.ra(j))then
            ra(i)=ra(j)
            rb(1,i)=rb(1,j)
            rb(2,i)=rb(2,j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        goto 20
        endif
        ra(i)=rra
        rb(1,i)=rrb1
        rb(2,i)=rrb2
      goto 10
      END
c
c------------------------------------------------------------------------
c
c					Numerical Recipes routine index
c
c------------------------------------------------------------------------
c
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
      REAL arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) stop 'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END
C  (C) Copr. 1986-92 Numerical Recipes Software '%1&9p#!.
c
c------------------------------------------------------------------------
c
c	second_v_area - calculates the area of a second-order Voronoi 
c			cell using an un-ordered list of vertices and a
c		        a closed formula. 
c
c	Input:
c		x(2)			input point	
c		p(2,n)			vertices of polygon.
c		n			number of vertices
c		a(n)			pseudo angles of nodes with 
c					respect to input point
c
c	Output:
c		area			area of polygon X 2
c
c	Comments:
c		 The vertices may be input in any order but they are 
c		 re-ordered anti-clockwise upon output.
c		
c		 Calls are made to pangle, hpsort_two.
c
c					M. Sambridge, RSES, May 1996.
c
c------------------------------------------------------------------------
c
      Subroutine second_v_area(x,n,p,a,area)

      real*8	x(2)
      real*8	p(2,*)
      real*4    theta
      real*4    a(*)
c	      					calculate pseudo angles
c						from first node
c
      do i=2,n
        a(i) = pangle(p(1,1),p(1,i))
      end do
      a(1) = -1
c
      theta = a(2)
      do i=2,n
         a(i) = a(i) - theta
        if(a(i).lt.0)a(i) = a(i) + 360.
      end do
      theta = a(2)
c
c     write(*,*)' unsorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						sort these nodes by 
c						pseudo angle from first edge
      call hpsort_two(n,a,p)

      if(a(2).ne.theta)write(*,*)' ERROR: first angle moved ?'

c     write(*,*)' sorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						calculate area of 
c						second-order Voronoi cell.
      area = 0.
      do i=1,n-1
         j = i+1
         area = area + p(1,i)*p(2,j) - p(2,i)*p(1,j)
      end do
      area = area + p(1,n)*p(2,1) - p(2,n)*p(1,1)
      area = abs(area)
      
      return
      end
c
c------------------------------------------------------------------------
c
c	second_v_area_d - calculates the area of a second-order Voronoi 
c			  cell using an un-ordered list of vertices and a
c		          a closed formula. 
c
c	Input:
c		x(2)			input point	
c		p(2,n)			vertices of polygon.
c		dp(4,2)			derivatives of first two vertices 
c					of polygon.
c		n			number of vertices
c		a(n)			pseudo angles of nodes with 
c					respect to input point
c
c	Output:
c		area			area of polygon X 2
c		df(2)			first derivatives of area X 2
c					      df(1) = df/dx,
c					      df(2) = df/dy
c
c	Comments:
c		 The first two vertices are assumed to be the vertices
c		 dependent on x, i.e. the vertices on the voronoi cell about x.
c		 These two vertices must be in clockwise order.
c		 The remaining vertices may be input in any order but they are 
c		 re-ordered anti-clockwise upon output.
c		
c		 Calls are made to pangle, hpsort_two.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
      Subroutine second_v_area_d(x,n,p,dp,a,area,df)

      real*8	x(2)
      real*8	p(2,*)
      real*8	dp(4,2)
      real*8	df(2)
      real*4    theta
      real*4    a(*)
c
c	      					calculate pseudo angles
c						from first node
c
      do i=2,n
        a(i) = pangle(p(1,1),p(1,i))
      end do
      a(1) = -1
c
c     write(*,*)' first pseudo angle =',theta 
c
      theta = a(2)
      do i=2,n
         a(i) = a(i) - theta
        if(a(i).lt.0)a(i) = a(i) + 360.
      end do
      theta = a(2)
c
c     write(*,*)' unsorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						sort these nodes by 
c						pseudo angle from first edge
      call hpsort_two(n,a,p)

      if(a(2).ne.theta)write(*,*)' ERROR: first angle moved ?'

c						plot nodes and v-cell
c     write(*,*)' sorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						calculate area of 
c						second-order Voronoi cell.
      area = 0.
      do i=1,n-1
         j = i+1
         area = area + p(1,i)*p(2,j) - p(2,i)*p(1,j)
      end do
      area = area + p(1,n)*p(2,1) - p(2,n)*p(1,1)
c
c						calculate 1st derivatives
      d222n = p(2,2) - p(2,n) 
      d1n12 = p(1,n) - p(1,2) 
      d2321 = p(2,3) - p(2,1)
      d1113 = p(1,1) - p(1,3)
      df(1) =   dp(1,1)*d222n + dp(2,1)*d1n12 
     &        + dp(1,2)*d2321 + dp(2,2)*d1113
      df(2) =   dp(3,1)*d222n + dp(4,1)*d1n12 
     &        + dp(3,2)*d2321 + dp(4,2)*d1113

      if(area.lt.0.)then
        df(1) = -df(1)
        df(2) = -df(2)
      end if

      area = abs(area)
       
      return
      end
c
c------------------------------------------------------------------------
c
c	second_v_area_dd - calculates the area of a second-order Voronoi 
c			   cell using an un-ordered list of vertices and a
c		           a closed formula. 
c
c	Input:
c		x(2)			input point	
c		p(2,n)			vertices of polygon.
c		dp(4,2)			derivatives of first 
c					two vertices of polygon.
c		n			number of vertices
c		a(n)			pseudo angles of nodes with 
c					respect to input point
c
c	Output:
c		area			area of polygon X 2
c		df(2)			first derivatives of area X 2
c					      df(1) = df/dx
c					      df(2) = df/dy
c		ddf(3)			second derivatives of area X 2
c					      ddf(1) = d2f/dxx
c					      ddf(2) = d2f/dyy
c					      ddf(3) = d2f/dxy
c
c	Comments:
c		 The first two vertices are assumed to be the vertices
c		 dependent on x, i.e. the vertices on the voronoi cell about x.
c		 These two vertices must be in clockwise order.
c		 The remaining vertices may be input in any order but they are 
c		 re-ordered anti-clockwise upon output.
c		
c		 Calls are made to pangle, hpsort_two and xplot routines.
c
c					M. Sambridge, RSES, April 1996.
c
c------------------------------------------------------------------------
c
      Subroutine second_v_area_dd(x,n,p,dp,ddp,a,area,df,ddf)

      real*8	x(2)
      real*8	p(2,*)
      real*8	dp(4,2)
      real*8	ddp(6,2)
      real*8	df(2),ddf(3)
      real*4    theta
      real*4    a(*)
c	      					calculate pseudo angles
c						from first node
c
      do i=2,n
        a(i) = pangle(p(1,1),p(1,i))
      end do
      a(1) = -1
c
c     write(*,*)' first pseudo angle =',theta 
c
      theta = a(2)
      do i=2,n
        a(i) = a(i) - theta
        if(a(i).lt.0)a(i) = a(i) + 360.
      end do
      theta = a(2)
c
c     write(*,*)' unsorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						sort these nodes by 
c						pseudo angle from first edge
      call hpsort_two(n,a,p)

      if(a(2).ne.theta)write(*,*)' ERROR: first angle moved ?'

c						plot nodes and v-cell
c     write(*,*)' sorted angles'
c     do i=1,n
c        write(*,*)i,' :',a(i),p(1,i),p(2,i)
c     end do
c						calculate area of 
c						second-order Voronoi cell.
      area = 0.
      do i=1,n-1
         j = i+1
         area = area + p(1,i)*p(2,j) - p(2,i)*p(1,j)
      end do
      area = area + p(1,n)*p(2,1) - p(2,n)*p(1,1)
c
c						calculate 1st derivatives
      d222n = p(2,2) - p(2,n) 
      d1n12 = p(1,n) - p(1,2) 
      d2321 = p(2,3) - p(2,1)
      d1113 = p(1,1) - p(1,3)
      df(1) =   dp(1,1)*d222n + dp(2,1)*d1n12 
     &        + dp(1,2)*d2321 + dp(2,2)*d1113
      df(2) =   dp(3,1)*d222n + dp(4,1)*d1n12 
     &        + dp(3,2)*d2321 + dp(4,2)*d1113

c						calculate 2nd derivatives
     
      ddf(1) =   ddp(1,1)*d222n + ddp(2,1)*d1n12 
     &         + ddp(1,2)*d2321 + ddp(2,2)*d1113
     &         + 2.*(dp(1,1)*dp(2,2) - dp(2,1)*dp(1,2))
      ddf(2) =   ddp(3,1)*d222n + ddp(4,1)*d1n12 
     &         + ddp(3,2)*d2321 + ddp(4,2)*d1113
     &         + 2.*(dp(3,1)*dp(4,2) - dp(4,1)*dp(3,2))
      ddf(3) =   ddp(5,1)*d222n + ddp(6,1)*d1n12 
     &         + ddp(5,2)*d2321 + ddp(6,2)*d1113
     &         + dp(1,1)*dp(4,2) - dp(2,1)*dp(3,2)
     &         + dp(3,1)*dp(2,2) - dp(4,1)*dp(1,2)

      if(area.lt.0.)then
        df(1) = -df(1)
        df(2) = -df(2)
        ddf(1) = -ddf(1)
        ddf(2) = -ddf(2)
        ddf(3) = -ddf(3)
      end if

      area = abs(area)
       
      return
      end

c
c------------------------------------------------------------------------
c
c       plot_c - dummy routine
c
c
c------------------------------------------------------------------------
c
        Subroutine plot_c(xs,ys,xs2,ys2)

	real*8		xs2,ys2
c                                               do nothing
        return
        end

c
c------------------------------------------------------------------------
c
c       plot_tc - dummy routine
c
c------------------------------------------------------------------------
c
        Subroutine plot_tc(n,points,vertices,centres)

        real*8          points(2,*)
        real*8          centres(3,*)
        integer         vertices(3,*)

c                                               do nothing
        return
        end
 
c
c------------------------------------------------------------------------
c
c	delaun - calculates delaunay triangulation incrementally 
c	 	 for a set of points in 2-D using a variation of
c		 Lawson's algorithm.
c
c	Input:
c		points(2,np)		array of node co-ordinates
c		num			number of nodes to be used
c               vis_tlist(nv_max)       List of triangles visible from 
c					current point.
c               vis_elist(nv_max)       List of edges visible from 
c					current point.
c               add_tlist(nv_max)       work array used by routine addpoint
c               eps                     distance from an interface for a
c                                       a point to be considered on an
c                                       interface (real*8). Prevents zero
c                                       area triangles resulting from rounding
c                                       error when nodes are co-linear.
c		nv_max			size of work arrays
c		mode			(=0,1,2,3) operation mode (see below)
c		inactive(np)		logical array. If mode=1 then the i-th
c					node is ignored if active(i) = .true.
c		nfirst			If mode=3 then nfirst is the first
c					node added to an existing triangulation
c		itstart			If mode=3 then itstart is a first
c					guess triangle containing node first node
c		subset(np)		logical array. If mode=2 then only
c					the nodes (subset(i),i=1,num) are used.
c
c	Output:
c               v(3,*)           	array of triangle vertices
c               numtri                  number of triangles in current
c                                       triangulation.
c               e(3,*)                  adjacency matrix of neighbouring
c                                       triangles. e(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2) 
c					in triangle j, stored counterclockwise 
c					about j.  
c                                       (This is the `opposite' definition)
c
c	Comments:
c
c       This routine calculates the Delaunay triangulation of a set of nodes 
c	using a variation of Lawson's method.  Each node is added sequentially 
c	and the Delaunay triangulation is updated. If the new node is inside 
c	the convex hull of the existing triangulation then the standard Lawson 
c	method is used. If it is outside then the list of triangle edges 
c	which are visible from the new point is calculated using routine 
c	visiblelist and each of these is used as the start of the swapping 
c	routine addpoint.
c
c	Four different operation modes are allowed.
c
c	MODE = 0:
c	The `standard' mode. All nodes from 1 to num are included. The arrays 
c	`subset' and `inactive' are unused and may be set to dummy variables
c       (saving memory). The variables nfirst and itstart are also unused.
c
c	MODE = 1:
c	All nodes from 1 to num are included except those for which
c	inactive(i) is set to true. The array `subset' is unused and may
c	be set to a dummy variable (saving memory). The variables nfirst 
c	and itstart are also unused.
c
c	MODE = 2:
c	Only nodes from subset(1) to subset(num) are included. The array
c	`inactive' is unused and may be set to a dummy variable
c	(saving memory). The variables nfirst and itstart are also unused.
c
c	MODE = 3:
c	Used to add nodes from nfirst to num to an existing triangulation.
c	Nodes for which inactive(i) is set to true are ignored.
c	The array `subset' is unused and may be set to a dummy variable.
c
c       The performance may be sensitive to the order in which the nodes are
c       added so these can be sorted before calling this routine if desired.
c
c	This routine was converted to use the `opposite' definition of
c	the adjacency matrix on 30/1/96.
c
c
c	Calls are made to triloc_del,visiblelist,insert_point,addpoint.
c
c					         M. Sambridge, Dec. 1994.
c					Modified by J. Braun, Sept. 1995.
c					(last change 30/1/96: multiple modes,
c					 and uses opposite definition of
c					 adjacency matrix)
c
c------------------------------------------------------------------------
c
	Subroutine delaun (points,num,e,v,numtri,numtri_max,
     &                     vis_tlist,vis_elist,add_tlist,eps,nv_max,
     &                     mode,inactive,nfirst,itstart,subset)

	real*8		points(2,*)
	real*8		x,y
        real*8          eps,del1,del2,del
 	integer		vis_tlist(*),vis_elist(*),add_tlist(*)
 	integer		v(3,*)
 	integer		e(3,*)
 	integer		subset(*)
	integer		t,p,ccw
        logical         out
	logical		newpoint
	logical		inactive(*)



        if (mode.eq.0.or.mode.eq.1.or.mode.eq.2) then

c					We are calculating Delaunay 
c					of all input points or a
c					subset of all input points

c					find first two active nodes
           if(mode.eq.0)then
              i1=1 
              i2=2 
              nodestart=3
           else if(mode.eq.1)then
              i1=0
              i2=0
              do i=1,num
                 if (i2.ne.0) goto 2222
                 if (i1.ne.0.and..not.inactive(i)) i2=i
                 if (i1.eq.0.and..not.inactive(i)) i1=i
              end do
 2222         continue
              nodestart = i2+1
           else if(mode.eq.2)then
              i1 = subset(1)
              i2 = subset(2)
              nodestart = 3
           end if
c                                       Find three non-colinear points
c                                       to form the first triangle 
           v(1,1) = i1
           v(2,1) = i2
           do 10 j=nodestart,num
              i = j
              if(mode.eq.2)then
                  i = subset(j)
              else if(mode.eq.1.and.inactive(i))then
                  go to 10
              end if
              istart=i
	      del1 = (points(2,i1)-points(2,i))
     &              *(points(1,i2)-points(1,i))
	      del2 = (points(1,i1)-points(1,i))
     &              *(points(2,i2)-points(2,i))
              del = del1-del2
              if(dabs(del).gt.eps) goto 11111
 10        continue
           stop 'all input data are in a line...'
11111      v(3,1) = istart

c					Initialize adjacency matrix
 	   e(1,1) = 0
 	   e(2,1) = 0
 	   e(3,1) = 0
c					Ensure initial triangle 
c					is in ccw order
c					
 	   if(ccw(points(1,v(1,1)),
     &            points(1,v(2,1)),
     &            points(1,v(3,1)),k).eq.-1)then
                  itemp = v(1,1)
                  v(1,1) = v(2,1)
                  v(2,1) = itemp
c	          write(*,*)' initial triangle was cw'
 	    end if
	
c					Initialize variables
 	    numtri = 1
	    t = 1

        else if (mode.eq.3) then
c					We are adding nodes to an 
c					existing triangulation
c					Perform initialization
           nodestart=nfirst
           t = itstart
           istart = 0
           if(t.le.0.or.t.gt.numtri)t=1

        end if 
c					Incrementally update the 
c					Delaunay triangulation

 	do 100 j=nodestart,num


           p = j
           if(mode.eq.1.or.mode.eq.3)then
             if (inactive(j)) goto 100
           else if(mode.eq.2)then
	     p = subset(j)
           end if
             
           if(p.eq.istart)go to 100

	   x = points(1,p)
	   y = points(2,p)

c					locate triangle 
c					containing current node

	   call triloc_del(x,y,points,v,e,t,eps,out,ipos,iface)


	   if(out)then

c					point is outside of convex hull, 
c					so find list of edges that are 
c					visible from current point 

	      call visiblelist(points,e,v,x,y,t,ipos,eps,
     &                         vis_tlist,vis_elist,nvis)

c					for each visible edge 
c	      				start swapping algorithm


              newpoint = .true.

	      if(nvis.gt.nv_max)then
                 write(*,*)' Error in subroutine delaun:'
                 write(*,*)' Too many visible triangles '//
     &                     ' from current point'
                 write(*,*)' Remedy: increase size of parameter nv_max'
                 write(*,*)'         in calling program'
                 write(*,*)'         Number of visible triangles '
                 write(*,*)'         for this point             =',nvis
                 write(*,*)'         Current value of nv_max    =',nv_max
                 stop
              end if


	      do 60 i=1,nvis
                 t = vis_tlist(i)
                 ipos = vis_elist(i)
                 jpos = mod(vis_elist(i),3)+1
c	         write(6,*)' visible t =',t,' node',v(ipos,t),v(jpos,t)


 	         call addpoint
     &           (points,e,v,p,t,ipos,numtri,newpoint,add_tlist) 


                 newpoint = .false.
 60           continue

	   else
 
c	      write(6,*)' point located in triangle',t

c					add node to inside of convex hull
c					using swapping algorithm

  	      call insertpoint(points,e,v,p,t,numtri,iface) 

           end if

           if (numtri.gt.numtri_max) then
              write (*,*) 'Error in subroutine delaun:'
              write (*,*) 'Too many triangles'
              write(*,*)' Remedy: increase size of parameter numtri_max'
              write(*,*)'         in calling program'
              write(*,*)'         Number of triangles '
              write(*,*)'         for this point             =',numtri
              write(*,*)'         Current value of numtri_max    =',
     &                  numtrimax_max
              stop
           endif

 100    continue


	return
	end

c
c-----------------------------------------------------------------------------
c
c	Subroutine qhullf (np,i,j,nt_max,k,points,nt,vertices)
c
c-----------------------------------------------------------------------------
c
	Subroutine qhullf (np,i,j,nt_max,k,points,nt,vertices)
c
	real*8		points(2,*)
	integer		vertices(3,*)

	write(*,*)' '
	write(*,*)' Error in subroutine nn2d_setup'
	write(*,*)' qhull is not installed'
	write(*,*)' Delaunay triangulation must be either'
	write(*,*)' calculated with routine delaun (dmode=-1 or -2)'
	write(*,*)' or read in from a file (dmode>0; logical unit=dmode)'
	write(*,*)' '

	stop
	end
c------------------------------------------------------------------------
c
c	Subroutine visiblelist - calculates all sides of triangles 
c		                 visible from the point (x,y), which is
c			         outside of the convex hull.
c			
c	Input:
c		points(2,*)		array of node co-ordinates	
c	        vertices(3,*)		array of triangle vertices	
c	        neighbour(3,*)		array of neighbouring triangles.	
c                                       neighbour(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c		x,y			Co-ordinates of test point p
c		t			Index of any triangle on hull 
c					that is visible from point p.
c					(Usually given by routine Triloc_del.)
c		tpos			Position of edge in triangle t
c					(using Sloan's adjacency convention)
c		eps			distance from an interface for a
c					a point to be considered on an 
c					interface (real*8). Prevents zero
c					area triangles resulting from rounding
c					error when nodes are co-linear.
c
c	Output:
c		nvis			Number of triangles visible from point p
c		vis_tlist		List of triangles visible from p
c		vis_elist		List of edges visible from p
c
c	Comments:
c		 Assumes point p is outside of the convex hull and vertices
c		 are in ccw order. Uses Sloan's definition of adjacency matrix.
c		
c		 This routine was converted from using Sloan's definition of
c		 the adjacency matrix to the `opposite' definition on 30/1/96.
c
c	Calls routine visible.
c
c					M. Sambridge, RSES, Nov 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c		
	Subroutine visiblelist
     &             (points,neighbour,vertices,x,y,t,tpos,eps,
     &              vis_tlist,vis_elist,nvis)

	real*8		points(2,*)
	real*8		x,y
	real*8		eps
 	integer		vertices(3,*)
 	integer		neighbour(3,*)
 	integer		vis_tlist(*),vis_elist(*)
 	integer		t,tpos,pos,t1,t2,edg,tnew
	logical		visible
	logical		special
	integer		c1(3),c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/

	nvis = 1
	vis_tlist(1) = t
	vis_elist(1) = tpos
    	inode = vertices(tpos,t)
cd      write(6,100)t,inode,vertices(mod(tpos,3)+1,t)
    	pos = c1(tpos)
      	jnode = vertices(pos,t)
    	t1 = neighbour(c2(pos),t)
        special = .false.
	if(t1.eq.0)then
	  t1 = t
          tnew = 0
          special = .true.
	end if

  5     continue
        if(.not.special)then
           pos = edg(t1,jnode,vertices)
           tnew = neighbour(c2(pos),t1)
cd         write(6,*)' tnew =',tnew,' t1',t1,' jnode',jnode
        end if
        special = .false.
        if(tnew.eq.0)then
  6        continue
  	   if(visible(x,y,points,vertices,t1,pos,eps))then
	      nvis = nvis + 1
	      vis_tlist(nvis) = t1
	      vis_elist(nvis) = pos
cd            write(6,100)t1,jnode,vertices(mod(pos,3)+1,t1)
	   else
cd	      write(6,200)t1,jnode,vertices(mod(pos,3)+1,t1)
              go to 10
	   end if
           pos = c1(pos)
	   jnode = vertices(pos,t1)
    	   tnew = neighbour(c2(pos),t1)
	   if(tnew.eq.0) go to 6
           t1 = tnew
	   go to 5 
        else
cd	   write(6,300)t1,jnode,vertices(mod(pos,3)+1,t1)
	   t1 = tnew
	   go to 5 
	end if

  10	jnode = inode
    	pos = c2(tpos)
    	t2 = neighbour(c2(pos),t)
        special = .false.
	if(t2.eq.0)then
	  t2 = t
          tnew = 0
          special = .true.
	end if

  15    continue
        if(.not.special)then
           pos = c2(edg(t2,jnode,vertices))
           tnew = neighbour(c2(pos),t2)
        end if
        special = .false.
        if(tnew.eq.0)then
  16       continue
  	   if(visible(x,y,points,vertices,t2,pos,eps))then
	      nvis = nvis + 1
	      vis_tlist(nvis) = t2
	      vis_elist(nvis) = pos
cd            write(6,100)t2,vertices(pos,t2),vertices(mod(pos,3)+1,t2)
	   else
cd	      write(6,200)t2,vertices(pos,t2),vertices(mod(pos,3)+1,t2)
              go to 20
	   end if
	   jnode = vertices(pos,t2)
    	   pos = c2(pos)
    	   tnew = neighbour(c2(pos),t2)
	   if(tnew.eq.0)go to 16
	   t2 = tnew
	   go to 15
        else
cd	   write(6,300)t2,jnode,vertices(mod(pos,3)+1,t2)
	   t2 = tnew
	   go to 15
	end if

 20	continue
	      
 100    format
     &  (1x,'Triangle',i6,' edge',i6,1x,i6,' is visible')
 200    format
     &  (1x,'Triangle',i6,' edge',i6,1x,i6,' is not visible')
 300    format
     &  (1x,'Triangle',i6,' edge',i6,1x,i6,' is not on convex hull')

	return
	end
c
c------------------------------------------------------------------------
c
c	Function visible - determines whether the triangle t is visible
c		           from the point p on edge tpos.
c	Input:
c		points(2,*)		array of node co-ordinates	
c	        vertices(3,*)		array of triangle vertices	
c	        t			Triangle to be tested
c	        tpos			Edge to be tested in triangle t 
c		eps			distance from an interface for a
c					a point to be considered on an 
c					interface (real*8). Prevents zero
c					area triangles resulting from rounding
c					error when nodes are co-linear.
c
c	Output:
c	        visible			Logical: = true if edge is visible  
c	                                         = false if edge is not visible
c
c	Comments:
c		 Assumes point p is outside of the convex hull and vertices
c		 are in ccw order. Uses Sloan's definition of adjacency matrix.
c		
c	Calls no other routines.
c
c					M. Sambridge, RSES, Nov 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Function visible(x,y,points,vertices,t,tpos,eps)

	real*8		points(2,*)
	real*8		del1,del2
	real*8		x,y
 	integer		vertices(3,*)
 	integer		t,tpos
	logical		visible
	real*8		eps,del
	integer		c1(3)
        save            c1
	data		c1/2,3,1/

        j = c1(tpos)
c						test edge tpos in triangle t
        i1 = vertices(tpos,t)
        i2 = vertices(j,t)
        del1 = (points(2,i1)-y)*(points(1,i2)-x)
        del2 = (points(1,i1)-x)*(points(2,i2)-y)
	del = del1-del2
        if(del.gt.eps)then
           visible = .true.
	else
           visible = .false.
	end if

	return
	end
c
c------------------------------------------------------------------------
c
c	addpoint - inserts a point into an existing delaunay triangulation 
c		   when point is outside of triangulation (but attached to
c		   triangle t) using the stacking procedure of Sloan.
c
c	Input:
c		points(2,np)		array of node co-ordinates	
c	        v(3,*)			array of triangle vertices	
c	        e(3,*)		        array of neighbouring triangles.	
c                                       e(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c		p			index of input point
c		t			triangle on convex hull visible 
c					from input point 
c		numtri			number of triangles in current
c					triangulation.
c		tpos			position of start node in triangle t
c		tri			list of triangles visible from point p
c		newpoint		logical = true if t is the first
c					triangle on the hull visible from p
c
c	Output:
c               v			updated
c               e			updated
c		numtri			updated
c
c	Comments:
c
c	The input nodes are in the form of a subset of an existing set
c	of points. This is so that extra arrays do not need to be used.
c	On input the vertices are assumed to be in ccw order.
c
c	When newpoint = false then there are multiple triangles
c	from the new point to the convex hull, and addpoint must
c	be called once for each attached triangle. In this case
c	the initialization of the adjacency list includes the
c	neighouring triangles already processed by addpoint, i.e.
c	those from point p to the hull.
c
c	This routine was converted from using Sloan's definition of
c	the adjacency matrix to the `opposite' definition on 30/1/96.
c
c					M. Sambridge, RSES, Nov. 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Subroutine addpoint (points,e,v,p,t,tpos,numtri,newpoint,tri) 

	real*8		points(2,*)
 	integer		v(3,*)
 	integer		e(3,*)
	integer		erl,era,erb,edg,v1,v2,v3,a,b,c,l,r
	integer		p,t,tpos,ccw
	logical		swap
	integer		tri(*)
	logical		newpoint
        save 		ip,lp
	integer		c1(3)
	integer		c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/



	if(newpoint)then
	   ip = 0
	   lp = 0
        end if

c			Add new node to existing triangulation

c					create new triangle
        numtri = numtri + 1
        v1 = v(tpos,t)
        v2 = v(c1(tpos),t)
        if(ccw(points(1,v1),
     &         points(1,v2),
     &         points(1,p),k).eq.-1)then
               itemp = v1
               v1 = v2
               v2 = itemp
        end if
        v(1,numtri) = p
        v(2,numtri) = v1
        v(3,numtri) = v2

c					initialize adjacency list including
c					neighbouring triangles attached
c					from the point to the hull.

        e(c2(1),numtri) = 0
        e(c2(2),numtri) = t
        e(c2(3),numtri) = 0

c				
        if(.not.newpoint)then
           do 10 j=1,lp
              k = tri(j)
              if(v(2,k).eq.v1)then
c                write(6,*)' v1 match with node 2'
c                write(6,*)' current triangle',numtri,' new',k
c                write(6,*)' nodes:',v(1,k),v(2,k),v(3,k) 
c                write(6,*)' e mat:',e(c2(1),k),e(c2(2),k),e(c2(3),k) 
                 e(c2(1),numtri) = k
                 e(c2(1),k) = numtri
              else if(v(3,k).eq.v1)then
                 e(c2(1),numtri) = k
                 e(c2(3),k) = numtri
              end if
              if(v(2,k).eq.v2)then
                 e(c2(3),numtri) = k
                 e(c2(1),k) = numtri
              else if(v(3,k).eq.v2)then
                 e(c2(3),numtri) = k
                 e(c2(3),k) = numtri
              end if
 10        continue
        end if



c
c					initialize stack

 	call stackinit



c                                       update adjacency list
c                                       for triangle on old boundary
        e(c2(tpos),t) = numtri


c					add new triangle on stack


	call push(numtri)


c					loop while stack is not empty

 50     continue

	call pop(L)
	r = e(c2(2),l)
c
c					check if new point is in circumcircle
c
	erl=edg(r,l,e)
        erl = c1(erl) 
	era=c1(erl)
	erb=c1(era)
	v1 = v(erl,r)
	v2 = v(era,r)
	v3 = v(erb,r)
	
	if(swap(points(1,v1),points(1,v2),
     &          points(1,v3),points(1,p)))then

c					new point is inside circumcircle
c					for triangle r so swap diagonal
           a=e(c2(era),r)
           b=e(c2(erb),r)
           c=e(c2(3),l)
c					update adjacency list for triangle l
	   v(3,l) = v3
	   e(c2(2),l) = a
	   e(c2(3),l) = r

c					update adjacency list for triangle r
	   v(1,r)=p
	   v(2,r)=v3
	   v(3,r)=v1
	   e(c2(1),r)=l
	   e(c2(2),r)=b
	   e(c2(3),r)=c

c					put edges l-a and r-b on stack
c					update adjacency list for 
c					triangles a and c
	   if(a.ne.0)then
	      e(edg(a,r,e),a)=l
              call push(l)
	   else
c					record triangles 
c					attached to new point
              ip = ip + 1
              tri(ip) = l
	   end if

	   if(b.ne.0)then
              call push(r)
	   else
c					record triangles 
c					attached to new point
              ip = ip + 1
              tri(ip) = r
           end if
	   if(c.ne.0) e(edg(c,l,e),c)=r

	else

c					record triangle attached to p
	   ip = ip + 1
           tri(ip) = l

	end if
	call stackempty(k)
	if(k.ne.1)go to 50
        call stackflush()

	lp = ip

c       write(6,*)' Number of triangles attached to last point',ip
c1	write(6,*)(tri(i),i=1,ip)
c	write(6,*)' triangles attached to last point on hull'
c       do 100 i=1,ip
c          it=tri(i)
c          do 101 k=1,3
c             l=mod(k,3)+1
c             if(e(k,it).eq.0)then
c                write(6,*)' t',it,' edge ',v(k,it),v(l,it)
c             end if
c 101      continue 
c 100   continue 
c

	return
	end
c
c------------------------------------------------------------------------
c
c	insertpoint - inserts a point into an existing delaunay triangulation 
c		      (when new point is inside triangle t) using the stacking 
c		      procedure of Sloan.
c
c	Input:
c		points(2,np)		array of node co-ordinates	
c	        v(3,*)			array of triangle vertices	
c	        e(3,*)		        array of neighbouring triangles.	
c                                       e(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c		p			index of input point
c		t			triangle containing input point
c		numtri			number of triangles in current
c					triangulation.
c	        iface			index of the face containing the
c					input point in triangle loc
c					(if point is on a face)
c
c	Output:
c
c               v			updated
c               e			updated
c		numtri			updated
c
c	Comments:
c
c	The new point is assumed to be inside the convex hull of the
c	existing triangulation.
c
c	The input nodes are in the form of a subset of an existing set
c	of points. This is so that extra arrays do not need to be used.
c	On input the vertices are assumed to be in ccw order.
c
c	This routine was converted from using Sloan's definition of
c	the adjacency matrix to the `opposite' definition on 30/1/96.
c
c					M. Sambridge, RSES, Nov. 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Subroutine insertpoint (points,e,v,p,t,numtri,iface) 

	real*8		points(2,*)
 	integer		v(3,*)
 	integer		e(3,*)
	integer		erl,era,erb,edg,v1,v2,v3,a,b,c,l,r
	integer		p,t
	logical		swap
	integer		c1(3)
	integer		c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/


c					add new node to existing triangulation
        if(iface.eq.0)then
	   a = e(c2(1),t)
	   b = e(c2(2),t)
	   c = e(c2(3),t)
	   v1 = v(1,t)
	   v2 = v(2,t)
	   v3 = v(3,t)
	   v(1,t) = p	
	   v(2,t) = v1	
	   v(3,t) = v2	
	   e(c2(1),t) = numtri+2	
	   e(c2(2),t) = a	
	   e(c2(3),t) = numtri+1	
        

c					create new triangles
c
           numtri = numtri + 1
	   v(1,numtri)=p
	   v(2,numtri)=v2
	   v(3,numtri)=v3
	   e(c2(1),numtri)=t
	   e(c2(2),numtri)=b
	   e(c2(3),numtri)=numtri+1
	   numtri = numtri + 1
	   v(1,numtri)=p
	   v(2,numtri)=v3
	   v(3,numtri)=v1
	   e(c2(1),numtri)=numtri-1
	   e(c2(2),numtri)=c
	   e(c2(3),numtri)=t
        else
           j = iface
           k = c1(j)
           i = c1(k)
	   a = e(c2(i),t)
	   b = e(c2(j),t)
	   c = e(c2(k),t)
	   v1 = v(i,t)
	   v2 = v(j,t)
	   v3 = v(k,t)
	   v(1,t) = p	
	   v(2,t) = v1	
	   v(3,t) = v2	
	   e(c2(1),t) = numtri+1
	   e(c2(2),t) = a	
	   e(c2(3),t) = 0 	

c					create new triangle
c
           numtri = numtri + 1
	   v(1,numtri)=p
	   v(2,numtri)=v3
	   v(3,numtri)=v1
	   e(c2(1),numtri)=0
	   e(c2(2),numtri)=c
	   e(c2(3),numtri)=t
        end if
  
c
c					initialize stack

 	call stackinit

c					add new triangles on stack
c					and update adjacency list

	if(a.ne.0)call push(t)

	if(b.ne.0)then
          e(edg(b,t,e),b)=numtri-1
          call push(numtri-1)
        end if

	if(c.ne.0)then
          e(edg(c,t,e),c)=numtri
          call push(numtri)
        end if
c					loop while stack is not empty

        if(a.eq.0.and.b.eq.0.and.c.eq.0)go to 100

 50     continue

	call pop(L)
	r = e(c2(2),l)
c
c					check if new point is in circumcircle
c
	erl=edg(r,l,e)
        erl = c1(erl)
	era=c1(erl)
	erb=c1(era)
	v1 = v(erl,r)
	v2 = v(era,r)
	v3 = v(erb,r)
	
	if(swap(points(1,v1),points(1,v2),
     &          points(1,v3),points(1,p)))then

c					new point is inside circumcircle
c					for triangle r so swap diagonal
           a=e(c2(era),r)
           b=e(c2(erb),r)
           c=e(c2(3),l)
c					update adjacency list for triangle l
	   v(3,l) = v3
	   e(c2(2),l) = a
	   e(c2(3),l) = r

c					update adjacency list for triangle r
	   v(1,r)=p
	   v(2,r)=v3
	   v(3,r)=v1
	   e(c2(1),r)=l
	   e(c2(2),r)=b
	   e(c2(3),r)=c

c					put edges l-a and r-b on stack
c					update adjacency list for 
c					triangles a and c
	   if(a.ne.0)then
	      e(edg(a,r,e),a)=l
              call push(l)
	   end if
	   if(b.ne.0) call push(r)
	   if(c.ne.0) e(edg(c,l,e),c)=r

	end if
	call stackempty(k)
	if(k.ne.1)go to 50
 100    continue
        call stackflush()

	return
	end
c
c------------------------------------------------------------------------
c
c	Function edg - finds edge in triangle l which is adjacent 
c		       to triangle k.
c
c		       (From Sloan 1987)
c
c------------------------------------------------------------------------
c
	Function edg(l,k,e)
c
	integer		l,k,i,e(3,*),edg
c
	do 10 i=1,3
	   if(e(i,l).eq.k)then
              edg = i
              return
           end if
 10     continue

	write(*,*)' ***Error in function edg***'
	write(*,*)' ***Triangles not adjacent***'
	write(*,*)' triangle = ',l,' looking for triangle',k

	stop
	end
c
c------------------------------------------------------------------------
c
c	logical function swap - checks to see if point p lies 
c			        inside circumcircle about points p1,p2,p3
c				using the algorithm of Cline and Renka
c				(see Sloan 1987).
c
c------------------------------------------------------------------------
c
	Function swap(p1,p2,p3,p)

	logical		swap

	real*8		p(2),p1(2),p2(2),p3(2)
	real*8		x13,y13,x23,y23,x1p,y1p,x2p,y2p
	real*8		cosa,cosb,sina,sinb

	x13=p1(1)-p3(1)
	y13=p1(2)-p3(2)
	x23=p2(1)-p3(1)
	y23=p2(2)-p3(2)
	x1p=p1(1)-p(1)
	y1p=p1(2)-p(2)
	x2p=p2(1)-p(1)
	y2p=p2(2)-p(2)

	cosa = x13*x23 + y13*y23
	cosb = x2p*x1p + y1p*y2p

	if((cosa.ge.0.d0).and.(cosb.ge.0.d0))then
            swap = .false.
        else if((cosa.lt.0.d0).and.(cosb.lt.0.d0))then
            swap = .true.
        else
            sina=x13*y23-x23*y13
            sinb=x2p*y1p-x1p*y2p
	    if((sina*cosb+sinb*cosa).lt.0.d0)then
                swap = .true.
            else
                swap = .false.
            end if
        end if

	return
	end
c
c------------------------------------------------------------------------
c
c	Triloc_del - locates the triangle containing point x,y
c
c	Input:
c		x,y			co-ordinates of input points	
c		points(2,np)		array of node co-ordinates	
c	        vertices(3,nt)		array of triangle vertices	
c	        neighbour(3,*)		array of neighbouring triangles.	
c                                       neighbour(i,j) is the triangle
c                                       which shares the face containing
c                                       node (mod(i,3)+1) and (mod(i,3)+2)
c                                       in triangle j, stored counterclockwise
c                                       about j.
c                                       (This is the `opposite' definition)
c	        loc			first guess of triangle containing
c					(x, y).
c		eps			distance from an interface for a
c					a point to be considered on an 
c					interface (real*8). Prevents zero
c					area triangles resulting from rounding
c					error when nodes are co-linear.
c
c	Output:
c	        loc			index of triangle containing 
c					input point.
c	        out			=true if (x,y) is outside of
c					the convex hull, otherwise = false. 
c	        k			index of face through which the
c					algorithm last passed (used by
c					routine visbilelist if out = .true.)
c	        iface			index of the face containing the
c					input point in triangle loc
c					(if point is on a face)
c
c	Comments:
c		 If (x,y) is outside convex hull loc is a visible triangle
c		 on the hull, out is set to .true., and k is set to the
c		 index of the face of triangle loc visible from the input point
c		 (used as a starting point by the routine visiblelist)
c
c		 This version also returns the parameter iface. 
c		 If iface .ne. 0 then the input point is on the face of 
c		 triangle t between nodes iface and mod(iface,3)+1 and 
c		 it is also on the convex hull.
c
c		 A point is assumed to be on the edge (or its extension)
c		 between two nodes if it is inside the triangle at a 
c		 distance >= eps.
c
c		 Can be extended to higher dimensions using a similar
c		 stepping mechanism but without angular test.
c
c	         This routine was converted from using Sloan's definition of
c	         the adjacency matrix to the `opposite' definition on 30/1/96.
c
c		 No calls to other routines.
c
c					M. Sambridge, RSES, Nov. 1994.
c					(Last updated 30/1/96)
c
c------------------------------------------------------------------------
c
	Subroutine Triloc_del
     &                       (x,y,points,vertices,neighbour,loc,eps,
     &                        out,k,iface)
c
	real*8		points(2,*)
	integer		vertices(3,*)
	integer		neighbour(3,*)
	integer		p1,p2
        logical		out
	real*8		x,y,del1,del2,del
	real*8		eps
	integer		c1(3),c2(3)
        save            c1,c2
	data		c1/2,3,1/
	data		c2/3,1,2/
        logical         new

	out = .false.
        new = .true.
        ic = 0

 10     continue
c					point is outside convex hull
        if(out)return
        iface = 0

        do 20 i=1,3
	   j = c1(i)
c	   k = c1(j)
c					definition of adjacency matrix

c					use Sloan's 
c					definition of adjacency matrix
	   k = i

           p1 = vertices(i,loc)
           p2 = vertices(j,loc)
	   del1 = (points(2,p1)-y)*(points(1,p2)-x)
	   del2 = (points(1,p1)-x)*(points(2,p2)-y)
           del = del1-del2
 	   if(dabs(del).le.eps)then
              iface = i
	   else if(del.gt.0.d0)then
	      if(neighbour(c2(k),loc).eq.0)then
                 out = .true.
	      else
	         loc = neighbour(c2(k),loc)
	      end if
              if(.not.new.and.loc.eq.loc1)then
                 write(*,100) 
                 write(*,*)' Current triangle:',loc, 
     &           ' last three:',loc1,loc2,loc3
                 write(*,*)' New point      x:',x,' y:',y
                 write(*,*)' Triangle ',loc,
     &           ' v:',(vertices(j,loc),j=1,3),
     &           ' n:',(neighbour(c2(j),loc),j=1,3)
c                write(*,*)' del',del,' del1',del1,' del2',del2
                 write(*,101) 
                 stop
              end if
              if(new)then
                ic = ic + 1
                if(ic.eq.3)new = .false.
              end if
              loc1 = loc2
              loc2 = loc3
              loc3 = loc
	      go to 10
	   end if
 20     continue
	
c						check if input point is
c						on the convex hull
c
        if(neighbour(c2(iface),loc).ne.0)iface = 0

c       if(iface.ne.0)then
c          j = mod(iface,3)+1
c          jj = vertices(iface,loc)
c          kk = vertices(j,loc)
c          write(*,*)' point on triangle between nodes ',
c    &               jj,' and',kk
c          write(*,*)' point is on the convex hull'
c       end if

 100    format(/'Error in subroutine Triloc_del:',//
     &  ' Infinite loop detected in walking triangle algorithm',/,
     &  ' Probably due to rounding error creating a flat triangle'/)

 101    format(/1x,'Remedy: '/
     &  ' Either increase size of parameter eps in calling routine '/
     &  ' or re-order input points by running program nn_hull '/)

	return
	end
c
c------------------------------------------------------------------------
c
c	Function ccw - used to test the orientation of three points
c
c		Input : 	points p1,p2 and p3 (vectors 2x1)
c				(e.g. p(1,2) = x co-ordinate of p2)
c
c		Output: 	ccw,I
c
c		ccw    k
c	  	 1     0	:The direction p1,p2,p3 is ccw (+ve)   
c	  	-1     0	:The direction p1,p2,p3 is  cw (-ve)   
c	  	 1     1	:p1,p2,p3 are colinear & p2 in middle  
c	  	-1     1	:p1,p2,p3 are colinear & p1 in middle
c	  	 0     1	:p1,p2,p3 are colinear & p3 in middle 
c
c
c				Calls no other routines.
c
c					M. Sambridge, RSES, April 1994.
c
c------------------------------------------------------------------------
c
      Integer Function ccw(p1,p2,p3,k)
c     
      real*8		p1(2),p2(2),p3(2)
      real*8		dx1,dx2,dy1,dy2,a,b
c
      dx1 = p2(1) - p1(1)
      dx2 = p3(1) - p1(1)
      dy1 = p2(2) - p1(2)
      dy2 = p3(2) - p1(2)
      a = dx1*dy2
      b = dy1*dx2
      if (a.gt.b)then
         k = 0
         ccw = 1
      else if(a.lt.b)then
         k = 0
         ccw = -1
      else if(dx1*dx2.lt.0.0.or.dy1*dy2.lt.0.0)then
         k = 1
         ccw = -1 
      else if((dx1*dx1+dy1*dy1).lt.(dx2*dx2+dy2*dy2))then
         k = 1
         ccw = 1 
      else
         k = 1
         ccw = 0 
      end if
      return
      end
c
c------------------------------------------------------------------------
c
c	function theta - returns a real number between 0 and 360
c		         which has the same ordering as the angle
c			 between the line (a,b) and the horizontal.
c
c------------------------------------------------------------------------
c
	function theta(a,b)

	real*8		a(2),b(2)
	real*4		theta

	dx = b(1) - a(1)  
	ax = abs(dx)
	dy = b(2) - a(2)  
	ay = abs(dy)

	theta = 0.
	d = ax+ay
	if(d.ne.0)theta=dy/d

	if(dx.lt.0.0)then
           theta = 2.-theta
	else if(dy.lt.0.0)then
           theta = 4.+theta
	end if
	theta = theta*90.

	return
	end

c
c------------------------------------------------------------------------
c
c       check_neighbour - checks whether the array neighbour is
c                         consistent.
c
c       Input:
c               neighbour(3,nt)         array of neighbouring tetrahedra.
c                                       Neighbour(i,j) is the tetrahedron
c                                       opposite node i in tetrahedron j,
c                                       stored counterclockwise about j.
c               nd                      number of dimensions
c               vertices(3,nt)          array of tetrahedron vertices
c               nt                      number of tetrahedra
c
c       Output:
c               consistent              logical set to true if the
c                                       test is passed
c
c       Comments:
c                This is a utility routine used to test whether the
c                calculated array neighbour is internally consistent.
c                A 2-D and 3-D version is contained with the one
c                routine. The routine is not call by any nn routines
c                but can be placed after the nn setup routine in the
c                calling program to check the consistency of the
c                neighbour matrix.
c
c                Calls no other routines.
c
c                                       M. Sambridge, RSES, May 1995.
c
c------------------------------------------------------------------------
c
        Subroutine check_neighbour(neighbour,nd,vertices,nt,consistent)

        integer         neighbour(nd+1,*)
        integer         vertices(nd+1,*)
        logical         consistent

        consistent = .true.

c                                               perform 2-D test
        if(nd.eq.2)then

        do 10 i=1,nt
           do 20 j=1,3
              k = neighbour(j,i)
              j1 = mod(j,3)+1
              j2 = mod(j1,3)+1
              i1 = vertices(j1,i)
              i2 = vertices(j2,i)
              if(k.eq.0)go to 20
              do 30 m=1,3
                 m1 = mod(m,3)+1
                 m2 = mod(m1,3)+1
                 k1 = vertices(m1,k)
                 k2 = vertices(m2,k)
                 if(neighbour(m,k).eq.i)then
                    if(k1.eq.i1.and.k2.eq.i2.or.
     &                 k1.eq.i2.and.k2.eq.i1)go to 20
                 else if(m.eq.3)then
                 write(*,*)' neighbour matrix inconsistent'
                 write(*,*)' triangle',i,' has neighbour',k,' index',j
                 write(*,*)' triangle',k,' has no neighbour',i
                 consistent = .false.
                 end if
 30           continue
              consistent = .false.
              write(*,*)
     &        ' neighbour found but vertices inconsistent'
              write(*,*)
     &        ' neighbours triangle',i,' and',k
              write(*,*)
     &        ' vertices of triangle',i,' :',i1,i2,vertices(j,i)
              write(*,*)
     &        ' vertices of triangle',k,' :',k1,k2,vertices(m,k)
 20        continue
 10     continue

c                                               perform 3-D test
        else if(nd.eq.3)then

        do 110 i=1,nt
           do 120 j=1,4
              k = neighbour(j,i)
              j1 = mod(j,4)+1
              j2 = mod(j1,4)+1
              j3 = mod(j2,4)+1
              i1 = vertices(j1,i)
              i2 = vertices(j2,i)
              i3 = vertices(j3,i)
              if(k.eq.0)go to 120
              do 130 m=1,4
                 m1 = mod(m,4)+1
                 m2 = mod(m1,4)+1
                 m3 = mod(m2,4)+1
                 k1 = vertices(m1,k)
                 k2 = vertices(m2,k)
                 k3 = vertices(m3,k)
                 if(neighbour(m,k).eq.i)then
                    if((k1.eq.i1.and.k2.eq.i2.and.k3.eq.i3).or.
     &                 (k1.eq.i1.and.k2.eq.i3.and.k3.eq.i2).or.
     &                 (k1.eq.i2.and.k2.eq.i1.and.k3.eq.i3).or.
     &                 (k1.eq.i2.and.k2.eq.i3.and.k3.eq.i1).or.
     &                 (k1.eq.i3.and.k2.eq.i1.and.k3.eq.i2).or.
     &                 (k1.eq.i3.and.k2.eq.i2.and.k3.eq.i1))go to 120
                 else if(m.eq.4)then
                 write(*,*)' neighbour matrix inconsistent'
                 write(*,*)' triangle',i,' has neighbour',k,' index',j
                 write(*,*)' triangle',k,' has no neighbour',i
                 consistent = .false.
                 end if
 130          continue
              consistent = .false.
              write(*,*)
     &        ' neighbour found but vertices inconsistent'
              write(*,*)
     &        ' neighbours triangle',i,' and',k
              write(*,*)
     &        ' vertices of triangle',i,' :',i1,i2,i3,vertices(j,i)
              write(*,*)
     &        ' vertices of triangle',k,' :',k1,k2,k3,vertices(m,k)
 120       continue
 110    continue

        end if

        return
        end

