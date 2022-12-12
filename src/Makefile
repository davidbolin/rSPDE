toInclude = ${R_LIBRARY_DIR}/INLA/include/

obj = cgeneric_mvnormdens.o cgeneric_aux_nonstat.o cgeneric_aux_nonstat_fixed.o cgeneric_rspde_stat_frac_model.o cgeneric_rspde_nonstat_general.o cgeneric_rspde_stat_general.o cgeneric_rspde_stat_parsim_gen.o cgeneric_rspde_stat_parsim_fixed.o cgeneric_rspde_stat_int.o cgeneric_rspde_nonstat_gen_fixed.o cgeneric_rspde_nonstat_int.o cgeneric_aux_nonstat_int.o

all : rSPDE.so

flags = -O2 -Wall -Wextra -fpic

%.o: %.c
	gcc $(flags) -lstdc++ -Iinclude -I$(toInclude)  -c $^ -o $@

%.o: %.cpp
	g++ $(flags)  -I$(toInclude) -I"/usr/include/eigen3/" -I"/usr/local/include/eigen3/" -c $^ -o $@

rSPDE.so: $(obj)
	g++ -shared *.o -o ../inst/shared/rspde_cgeneric_models.so -lblas -llapack

clean :
	rm -f *.o

.PHONY: all clean