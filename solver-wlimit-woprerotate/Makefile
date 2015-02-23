#include folder of headers
INCLUDE_PATH = include

# Standard compiler
CC = g++

# optimized and correct code
 CFLAGS = -O3

# code needing debugging
#  CFLAGS = -g -Wno-deprecated

#math library, std c++ lib
  LIBS = -lm -lstdc++

OBJ_OUT = driver.o mesh_out.o mesh_read.o mesh_init.o crs.o Lhash.o residual.o analyticjacobian.o splitflux.o gradients.o mat_build.o lin_solve.o Lu4.o fsub.o bsub.o List.o

all: EULER

EULER: $(OBJ_OUT)
	$(CC) $(C_FLAGS) -o $@ $(OBJ_OUT) $(LIBS)

clean:
	\rm -rf *.o EULER

cleanish:
	\rm -rf *.dat

cleanest:
	\rm -rf *.o *.dat EULER
.cpp.o:
	$(CC) $(CFLAGS) -I$(INCLUDE_PATH) -c $<
