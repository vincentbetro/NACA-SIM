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

OBJ_OUT = driver.o opt_smooth.o cost_funct.o xy_pert.o edge_pert.o tri_edge.o subdivision.o refine_tri.o gradients.o area_calc.o analytic_circ.o analytic_curve.o make_nbrs_HO.o trimesh.o winslow.o linear_elastic.o lin_solve.o mat_build.o green_gauss.o young_mod.o crs.o Lhash.o List.o 

all: SMOOTH

SMOOTH: $(OBJ_OUT)
	$(CC) $(C_FLAGS) -o $@ $(OBJ_OUT) $(LIBS)

clean:
	\rm -rf *.o SMOOTH

cleanish:
	\rm -rf *.dat smoothmesh*

cleanest:
	\rm -rf *.o *.dat smoothmesh* SMOOTH

.cpp.o:
	$(CC) $(CFLAGS) -I$(INCLUDE_PATH) -c $<
