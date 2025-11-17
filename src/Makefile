# Include paths and flags
toInclude = $(R_LIBRARY_DIR)/INLA/include/
PKG_CFLAGS = -I$(R_HOME)/include

# Object files
obj = cgeneric_mvnormdens.o cgeneric_aux_nonstat.o cgeneric_aux_nonstat_fixed.o \
      cgeneric_rspde_stat_frac_model.o cgeneric_rspde_nonstat_general.o \
      cgeneric_rspde_stat_general.o cgeneric_rspde_stat_parsim_gen.o \
      cgeneric_rspde_stat_parsim_fixed.o cgeneric_rspde_stat_int.o \
      cgeneric_rspde_nonstat_gen_fixed.o cgeneric_rspde_nonstat_int.o \
      cgeneric_aux_nonstat_int.o \
	  cgeneric_gpgraph_alpha1.o cgeneric_gpgraph_alpha2.o \
	  cgeneric_aux_gpgraph_alpha2.o cgeneric_gpgraph_alpha1_directional.o cgeneric_aux_gpgraph_alpha1_directional.o \
	  inla_cgeneric_rspde_intrinsic_int_model.o \
      cgeneric_aux_intrinsic.o cgeneric_rspde_intrinsic_eigen_cache.o \
      cgeneric_rspde_spacetime.o cgeneric_aux_spacetime.o \
	  cgeneric_rspde_anisotropic2d.o cgeneric_aux_anisotropic.o \
	  cgeneric_rspde_1d_general.o cgeneric_aux_1d.o 

# Compiler settings
CC = /usr/local/bin/gcc
CXX = /usr/local/bin/g++

# Eigen and library paths
EIGEN_MAC = /opt/homebrew/
EIGEN_MAC_VAR = /usr/local
EIGEN_LINUX = /usr
MAC_LIB = /usr/local/brewlib
MAC_INCLUDE = /usr/local/brewinclude

# Architecture detection
ARCH := $(shell uname -m)

# Compiler flags based on system
ifeq ($(shell uname -s), Linux)
    flags = -O2 -Wall -Wextra -fpic -g -fopenmp
    LINK_FLAGS = -L$(R_HOME)/lib -lR -ldl -lm
else ifeq ($(shell uname -s), Darwin)
    ifeq ($(ARCH), arm64)
        flags = -O2 -Wall -Wextra -fpic -g -I$(MAC_INCLUDE) -fopenmp -arch arm64
        LINK_FLAGS = -L$(R_HOME)/lib -lR -ldl -lm -arch arm64
    else ifeq ($(ARCH), x86_64)
        flags = -O2 -Wall -Wextra -fpic -g -I$(MAC_INCLUDE) -fopenmp -arch x86_64
        LINK_FLAGS = -L$(R_HOME)/lib -lR -ldl -lm -arch x86_64
    endif
endif

# Default target
all: rSPDE.so

# Check or download cgeneric.h
cgeneric.h:
	if [ ! -f ${R_LIBRARY_DIR}/INLA/include/cgeneric.h ]; then \
		if [ ! -f cgeneric.h ]; then \
			echo "Downloading cgeneric.h..."; \
			wget -O cgeneric.h https://raw.githubusercontent.com/hrue/r-inla/devel/inlaprog/src/cgeneric.h || \
			curl -o cgeneric.h https://raw.githubusercontent.com/hrue/r-inla/devel/inlaprog/src/cgeneric.h || \
			(echo "Failed to download cgeneric.h" && exit 1); \
		fi; \
	fi

# Object file compilation rules
%.o: %.c cgeneric.h
	$(CC) $(flags) $(PKG_CFLAGS) -Iinclude -I$(toInclude) -L$(R_HOME)/lib -lR -c $< -o $@

%.o: %.cpp cgeneric.h
	$(CXX) $(flags) -I$(toInclude) -I$(EIGEN_MAC_VAR)/include/eigen3/ -I$(EIGEN_MAC)/include/eigen3/ -I$(EIGEN_LINUX)/include/eigen3/ -c $< -o $@

# Link shared library
rSPDE.so: cgeneric.h $(obj)
	$(CXX) $(flags) -shared $(obj) -o ../inst/shared/rspde_cgeneric_models.so -lblas -llapack $(LINK_FLAGS)

# Clean up
clean:
	rm -f *.o
	rm -f cgeneric.h
	rm -f ../inst/shared/rspde_cgeneric_models.so

.PHONY: all clean