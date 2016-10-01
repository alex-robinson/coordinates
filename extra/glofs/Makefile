INCLUDES = -I.
AR = ar rv
RANLIB = ranlib
# add path to c compiler 
CC = gcc 
FC = gfortran
# add c compiler optimization flags
COPTIONS = -D_UNDERSCORE
OPTFLAGS = -O2
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)
LD = $(CC) -L. 

EXECnos=.

CMD =            $(EXECnos)/nos_glofs_interp_sfcmarobs

OBJS =  stack.o stackpair.o volume.o nos_interp_sfcmarobs.o nnq.o juliansubs.o


.c.o:
	$(CC) $(CFLAGS) -c stack.c stackpair.c volume.c
.F.o:
	$FC -c -o nos_interp_sfcmarobs.o nos_interp_sfcmarobs.f
	$FC -c -o nnq.o nnq.f 
	$FC -c -o juliansubs.o juliansubs.f 

all:		$(CMD)

$(CMD):		$(OBJS)
	$(FC) $(LDFLAGS)  -o $(@) $(OBJS) $(LIBS)
	rm -f $(OBJS)

clean:
	rm -f $(OBJS)

clobber:	clean
	rm -f $(CMD) $(CMD).prof

void:	clobber
	rm -f $(SRCS) makefile

