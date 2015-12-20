############################################################
ANIHOME = ~/ani2D-3.1
include $(ANIHOME)/src/Rules.make

LIBS	= -L $(ANILIB) $(LIBAFT) $(LIBFEM) $(LIBLU) $(LIBLMR) $(LIBMBA) \
	              $(LIBVIEW) $(LIBLAPACK) $(LIBBLAS) -lm 
INCLUDE := $(INCLUDE) -I $(ANIHOME)/include

############################################################
all: main

run: main
	./main

runi: run
	make -C png

OBJ  =  main.o solve_u.o mesh.o


############################################################
main: $(OBJ) 
	$(FLINKER) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)

############################################################
$(OBJ): th.fh


############################################################
.f.o:
	$(F77) $(FFLAGS) -c $*.f -o $*.o $(INCLUDE)
.c.o:
	$(CC)  $(CFLAGS) -c $*.c -o $*.o

clean:
	make -C ps clean
	rm -f *.o main