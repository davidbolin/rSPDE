toInclude = ${R_LIBRARY_DIR}/INLA/include/ 

sources = cgeneric_gpgraph_alpha1.c cgeneric_rspde_stat_frac_model.c cgeneric_rspde_stat_general.c cgeneric_rspde_stat_parsim_gen.c cgeneric_rspde_stat_parsim_fixed.c cgeneric_rspde_stat_int.c

all : $(sources:.c=.o) so

%.o : %.c
	@echo gnerating $@
	gcc -c -O3 -Wall -Wextra -fpic -I$(toInclude)  $^ -o $@

so : $(sources:.c=.o)
	@echo generate .so
	gcc -shared *.o -o ../inst/shared/rspde_cgeneric_models.so -lblas -llapack


clean :
	rm -f *.o

.PHONY: all clean