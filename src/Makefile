toInclude = ${R_LIBRARY_DIR}/INLA/include/

obj = cgeneric_mvnormdens.o cgeneric_aux_nonstat.o cgeneric_aux_nonstat_fixed.o \
      cgeneric_rspde_stat_frac_model.o cgeneric_rspde_nonstat_general.o \
      cgeneric_rspde_stat_general.o cgeneric_rspde_stat_parsim_gen.o \
      cgeneric_rspde_stat_parsim_fixed.o cgeneric_rspde_stat_int.o \
      cgeneric_rspde_nonstat_gen_fixed.o cgeneric_rspde_nonstat_int.o \
      cgeneric_aux_nonstat_int.o \
	  cgeneric_gpgraph_alpha1.o cgeneric_gpgraph_alpha2.o \
	  cgeneric_aux_gpgraph_alpha2.o cgeneric_gpgraph_alpha1_directional.o cgeneric_aux_gpgraph_alpha1_directional.o \
	  inla_cgeneric_rspde_intrinsic_int_model.o \
      cgeneric_aux_intrinsic.o cgeneric_rspde_intrinsic_eigen_cache.o

all : rSPDE.so

CC = /usr/local/bin/gcc
CXX = /usr/local/bin/g++

EIGEN_MAC = /usr/local
EIGEN_LINUX = /usr
MAC_LIB = /usr/local/brewlib
MAC_INCLUDE = /usr/local/brewinclude

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    flags = -O2 -Wall -Wextra -fpic -g -fopenmp
endif
ifeq ($(UNAME_S),Darwin)
    flags = -Wl,-ld_classic -O2 -Wall -Wextra -fpic -g -I$(MAC_INCLUDE) -fopenmp
endif

%.o: %.c
	test -f ${R_LIBRARY_DIR}/INLA/include/cgeneric.h || test -f cgeneric.h || wget -O cgeneric.h https://raw.githubusercontent.com/hrue/r-inla/devel/inlaprog/src/cgeneric.h
	$(CC) $(flags) -Iinclude -I$(toInclude)  -c $^ -o $@ -L$(MAC_LIB) 

%.o: %.cpp
	$(CXX) $(flags) -I$(toInclude) -I$(EIGEN_MAC)/include/eigen3/ -I$(EIGEN_LINUX)/include/eigen3/ -c $^ -o $@

rSPDE.so: $(obj)
	$(CXX) $(flags) -shared *.o -o ../inst/shared/rspde_cgeneric_models.so -lblas -llapack

clean :
	rm -f *.o
	rm -f cgeneric.h

.PHONY: all clean