toInclude = "/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/INLA/include/"

obj = cgeneric_mvnormdens.o cgeneric_aux_nonstat.o cgeneric_aux_nonstat_fixed.o cgeneric_rspde_stat_frac_model.o cgeneric_rspde_nonstat_general.o cgeneric_rspde_stat_general.o cgeneric_rspde_stat_parsim_gen.o cgeneric_rspde_stat_parsim_fixed.o cgeneric_rspde_stat_int.o cgeneric_rspde_nonstat_gen_fixed.o

all : rSPDE.so

flags = -O2 -Wall -Wextra -fpic

%.o: %.c
	clang $(flags) -lstdc++ -Iinclude -I$(toInclude)  -c $^ -o $@

%.o: %.cpp
	clang++ $(flags)  -I$(toInclude) -I"/opt/homebrew/opt/eigen/include/eigen3" -c $^ -o $@

rSPDE.so: $(obj)
	clang++ -shared *.o -o rSPDE.so -lblas -llapack