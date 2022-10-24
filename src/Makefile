toInclude = ${R_LIBRARY_DIR}/INLA/include/

sources = $(wildcard *.c)

all : $(sources:.c=.o) so

%.o : %.c
	@echo gnerating $@
	gcc -c -O3 -Wall -Wextra -fpic -lblas -I$(toInclude) $^ -o $@

so : $(sources:.c=.o)
	@echo generate .so
	gcc -shared -lblas *.o -o rspde_cgeneric_models.so


clean :
	@##-rm -f *.o

.PHONY: all clean