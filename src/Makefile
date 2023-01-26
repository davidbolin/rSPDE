toInclude = ${R_LIBRARY_DIR}/INLA/include/

obj = cgeneric_mvnormdens.o cgeneric_aux_nonstat.o cgeneric_aux_nonstat_fixed.o cgeneric_rspde_stat_frac_model.o cgeneric_rspde_nonstat_general.o cgeneric_rspde_stat_general.o cgeneric_rspde_stat_parsim_gen.o cgeneric_rspde_stat_parsim_fixed.o cgeneric_rspde_stat_int.o cgeneric_rspde_nonstat_gen_fixed.o cgeneric_rspde_nonstat_int.o cgeneric_aux_nonstat_int.o

all : rSPDE.so

CC = gcc
CXX = g++

EIGEN_MAC = /usr/local/include/eigen3/
EIGEN_LINUX = /usr/include/eigen3/

flags = -O2 -Wall -Wextra -fpic

%.o: %.c
	$(CC) $(flags) -Iinclude -I$(toInclude)  -c $^ -o $@

%.o: %.cpp
	$(CXX) $(flags)  -I$(toInclude) -I$(EIGEN_MAC) -I$(EIGEN_LINUX) -c $^ -o $@

rSPDE.so: $(obj)
	$(CXX) -shared *.o -o ../inst/shared/rspde_cgeneric_models.so -lblas -llapack

clean :
	rm -f *.o

.PHONY: all clean